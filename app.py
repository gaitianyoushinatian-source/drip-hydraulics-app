# -*- coding: utf-8 -*-
"""
滴头灌溉水力计算 — 交互式可视化
运行：streamlit run app.py
"""
import streamlit as st
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
from io import BytesIO
from collections import Counter

st.set_page_config(page_title='滴头灌溉水力分析', layout='wide')

# 阻止浏览器自动翻译（"等间距"不会被译成"等浓度"，"Cu"不会被译成"铜"）
st.markdown('<meta name="google" content="notranslate">', unsafe_allow_html=True)
st.markdown("""
<div translate="no" id="no-translate">
""", unsafe_allow_html=True)

# ============ 色彩方案 ============
C_INIT = '#6090c1'        # 初始状态-中蓝
C_REPL = '#f2724d'        # 替换后-暖橙
C_EQ_INIT = '#6090c1'     # 等间距初始
C_EQ_REPL = '#115fa4'     # 等间距替换
C_VAR_INIT = '#f2724d'    # 变间距初始
C_VAR_REPL = '#ab2428'    # 变间距替换
C_HIST_INIT = '#acd2e5'   # 直方图初始
C_HIST_REPL = '#f2724d'   # 直方图替换
C_HIGHLIGHT = '#d7312d'   # 高亮标注
ZONE_COLORS = ['#82c9ff', '#acd2e5', '#fef9b7', '#fee395', '#daa81c', '#f79015', '#fb8275', '#ab2428']
HEATMAP_CU_SCALE = [
    [0.0, '#d73027'], [0.17, '#fc8d59'], [0.33, '#fee08b'],
    [0.5, '#ffffbf'], [0.67, '#d9ef8b'], [0.83, '#91cf60'], [1.0, '#1a9850'],
]
HEATMAP_DELTA_SCALE = [
    [0.0, '#fef9b7'], [0.25, '#fdd835'], [0.5, '#f79015'],
    [0.75, '#f2724d'], [1.0, '#d7312d'],
]

# ============ 嵌入滴头统一管道模块 ============
class _ModuleDict(dict):
    """支持属性访问的字典 — dp.G 等价于 dp['G']"""
    __getattr__ = dict.__getitem__
    __setattr__ = dict.__setitem__
    __delattr__ = dict.__delitem__

_DP_SOURCE = r'''
import numpy as np
from collections import Counter

# ---- 水力参数 (运行时会由 apply_config 覆写) ----
G = 9.81
NU = 1.004e-6
XI = 0.45
FLOW_EXPONENT = 0.5
D_PIPE = 0.016
H_INLET = 10.0
Q_INIT_RATED = 2.5
H_SPEC_RATED = 10.0
N_ZONES = 3
MAX_ITER = 200
TOLERANCE = 1e-8
DAMPING = 0.5
A_PIPE = np.pi * D_PIPE**2 / 4
k_init = Q_INIT_RATED / H_SPEC_RATED**FLOW_EXPONENT
SPEC_K = {}
FLOW_EXPONENT_INIT = 0.5
FLOW_EXPONENT_REPL = 0.5
LAST_SOLVE = {}

# ============ 物理函数 ============

def ff(Re):
    """沿程阻力系数：层流泊肃叶，湍流布拉休斯"""
    return 64.0 / Re if Re < 2000 else 0.3164 / Re**0.25

def cu(q):
    """克里斯琴森均匀度 Cu (%)"""
    mean_q = np.mean(q)
    if len(q) == 0 or mean_q <= 0:
        return 0.0
    return (1 - np.sum(np.abs(q - mean_q)) / (len(q) * mean_q)) * 100

# ============ 递减流量求解器 ============

def solve(positions, k_arr):
    """递减流量法水力求解器（阻尼不动点迭代）"""
    n = len(k_arr)
    h = np.full(n, H_INLET)
    damping = DAMPING
    prev_err = float('inf')
    for _ in range(MAX_ITER):
        h_old = h.copy()
        q_m3s = np.array([k_arr[i] * max(h[i], 0.01)**FLOW_EXPONENT
                         for i in range(n)]) / 3.6e6
        q_cum = np.cumsum(q_m3s[::-1])[::-1]
        h_cur = H_INLET
        h_new = np.zeros(n)
        for i in range(n):
            v_app = max(q_cum[i] / A_PIPE, 1e-12)
            hj_i = XI * v_app**2 / (2 * G)
            h_new[i] = max(h_cur - hj_i, 0.01)
            if i < n - 1:
                v_seg = max(q_cum[i+1] / A_PIPE, 1e-12)
                Re = v_seg * D_PIPE / NU
                lam = ff(Re)
                dl = positions[i+1] - positions[i]
                h_cur = h_new[i] - lam * dl * v_seg**2 / (2 * G * D_PIPE)
        err = np.max(np.abs(h_new - h_old))
        if err > prev_err * 1.1:
            damping = max(damping * 0.7, 0.02)
        else:
            damping = min(damping * 1.02, DAMPING)
        prev_err = err
        h = h + damping * (h_new - h)
        if err < TOLERANCE:
            break
    q_lph = np.array([k_arr[i] * max(h[i], 0.01)**FLOW_EXPONENT for i in range(n)])
    global LAST_SOLVE
    LAST_SOLVE = {
        'converged': bool(err < TOLERANCE),
        'iterations': _ + 1,
        'error': float(err),
        'pressure_limited': bool(np.any(h <= 0.010001)),
    }
    return h, q_lph

# ============ 位置生成 ============

def generate_eq_positions(spacing, total_l):
    """生成等间距滴头位置"""
    n = int(total_l / spacing) + 1
    return np.array([i * spacing for i in range(n)])

def generate_var_positions(first_seg, cycle, total_l):
    """生成变间距滴头位置：首段first_seg, 然后按cycle循环累加"""
    pos = [0.0]
    next_p = first_seg
    if next_p > total_l:
        return np.array(pos)
    pos.append(next_p)
    ci = 0
    while True:
        seg = cycle[ci % len(cycle)]
        ci += 1
        next_p = pos[-1] + seg
        if next_p > total_l + 1e-9:
            break
        pos.append(next_p)
    return np.array(pos)

# ============ 末端水头截取 ============

def find_limit_length(positions, k_uniform):
    """二分搜索：找到末端水头衰减到 ~0.1m 的最大滴头数"""
    global MAX_ITER, TOLERANCE
    total = len(positions)
    k_arr_template = np.full(total, k_uniform)
    lo, hi = 5, total
    best_result = None
    LIMIT_THRESH = 0.1
    old_iter, old_tol = MAX_ITER, TOLERANCE
    MAX_ITER, TOLERANCE = 200, 1e-6
    while lo <= hi:
        mid = (lo + hi) // 2
        pos_sub = positions[:mid].copy()
        k_sub = k_arr_template[:mid].copy()
        h, q = solve(pos_sub, k_sub)
        if h[-1] > LIMIT_THRESH:
            best_result = (pos_sub, k_sub, h, q)
            lo = mid + 1
        else:
            hi = mid - 1
    MAX_ITER, TOLERANCE = old_iter, old_tol
    if best_result is None:
        pos_sub = positions[:5].copy()
        k_sub = k_arr_template[:5].copy()
        h, q = solve(pos_sub, k_sub)
        best_result = (pos_sub, k_sub, h, q)
    pos_sub, k_sub, h, q = best_result
    return pos_sub, k_sub, h, q, len(pos_sub)

def find_nearest_tail(positions, target_h, k_uniform):
    """全范围查找最接近目标末端水头的滴头数；目标不可达时明确报错。"""
    total = len(positions)
    k_arr_template = np.full(total, k_uniform)
    best_n = 1
    best_diff = float('inf')
    best_result = None
    for n in range(1, total + 1):
        pos_sub = positions[:n].copy()
        k_sub = k_arr_template[:n].copy()
        h, q = solve(pos_sub, k_sub)
        diff = abs(h[-1] - target_h)
        if diff < best_diff:
            best_diff = diff
            best_n = n
            best_result = (pos_sub, k_sub, h, q)
    if target_h > H_INLET:
        raise ValueError(f'目标末端水头 {target_h:.4g} m 高于入口水头 {H_INLET:.4g} m')
    h_full, _ = solve(positions, k_arr_template)
    if target_h < h_full[-1] - 1e-6:
        raise ValueError(
            f'目标末端水头 {target_h:.4g} m 低于当前长度的可达值 {h_full[-1]:.4g} m；'
            '请增大管道长度或调整水力参数')
    pos_sub, k_sub, h, q = best_result
    return pos_sub, k_sub, h, q, best_n

# ============ 分区 ============

def get_zone_labels(n_zones):
    """返回各区中文标签"""
    if n_zones == 1:
        return ['全区']
    elif n_zones == 2:
        return ['高压区', '低压区']
    elif n_zones == 3:
        return ['高压区', '中压区', '低压区']
    elif n_zones == 4:
        return ['高压区', '中高压区', '中低压区', '低压区']
    else:
        labels = ['高压区']
        for i in range(2, n_zones):
            labels.append(f'区{i}')
        labels.append('低压区')
        return labels

def zones_by_integral(positions, h, n_zones):
    """按压力积分等分分区"""
    n = len(positions)
    integral = np.zeros(n)
    for i in range(1, n):
        dl = positions[i] - positions[i-1]
        integral[i] = integral[i-1] + (h[i] + h[i-1]) / 2 * dl
    total = integral[-1]
    boundaries = [i * total / n_zones for i in range(n_zones + 1)]
    zone_ids = np.full(n, n_zones - 1, dtype=int)
    for zid in range(n_zones - 2, -1, -1):
        zone_ids[integral <= boundaries[zid + 1]] = zid
    return zone_ids, boundaries

def compute_q_target(q, zone_ids, n_zones):
    """取实际非空分区中的中间压力区平均流量作为目标流量。"""
    valid_zones = np.unique(zone_ids)
    if len(valid_zones) == 0:
        raise ValueError('分区结果为空，请检查管道布局')
    mid_pos = (len(valid_zones) - 1) // 2
    mid = valid_zones[mid_pos]
    if len(valid_zones) % 2 == 1:
        return float(np.mean(q[zone_ids == mid]))
    next_zone = valid_zones[mid_pos + 1]
    return float((np.mean(q[zone_ids == mid]) + np.mean(q[zone_ids == next_zone])) / 2)

# ============ 规格选择 ============

def select_zone_specs(zone_ids, k_need, n_zones, spec_k):
    """为每个分区独立选择统一规格（L1误差最小化）"""
    zone_spec = {}
    k_need_stats = []
    for zid in range(n_zones):
        mask = zone_ids == zid
        cnt = np.sum(mask)
        if cnt == 0:
            continue
        choice_cnt = {s: 0 for s in spec_k}
        for i in np.where(mask)[0]:
            best_s = min(spec_k, key=lambda s: abs(k_need[i] - spec_k[s]))
            choice_cnt[best_s] += 1
        best_s = min(spec_k, key=lambda s: np.sum(np.abs(k_need[mask] - spec_k[s])))
        best_err = np.sum(np.abs(k_need[mask] - spec_k[best_s]))
        zone_spec[zid] = (best_s, spec_k[best_s])
        k_need_stats.append({
            'name': get_zone_labels(n_zones)[zid],
            'cnt': cnt,
            'kn_min': np.min(k_need[mask]),
            'kn_max': np.max(k_need[mask]),
            'kn_avg': np.mean(k_need[mask]),
            'choice_dist': ', '.join(
                f'{s}L/h:{c}个' for s, c in sorted(choice_cnt.items())),
            'winner': best_s,
            'l1_err': best_err,
        })
    return zone_spec, k_need_stats
'''

@st.cache_resource
def _init_dp():
    ns = _ModuleDict(__name__='dp')
    exec(_DP_SOURCE, ns)
    return ns

dp = _init_dp()
# ================================================================================

# ============ 辅助函数 ============
def auto_specs(q_init, n_zones):
    """根据初始流量和分区数自动生成可用规格（步长按流量比例）"""
    step = max(0.05, round(q_init * 0.2, 2))
    # 区数少时生成足够候选，覆盖 q_target 的可能范围
    # q_target 通常低于 q_init（沿程压降），所以向下多覆盖
    if n_zones <= 2:
        n_gen = 7
        step = max(0.03, round(q_init * 0.12, 2))
        specs = [round(q_init + (i - (n_gen - 1) / 2) * step, 2) for i in range(n_gen)]
        specs = [max(0.05, s) for s in specs]
        return sorted(set(specs))
    elif n_zones <= 4:
        n_gen = n_zones + 2
    else:
        n_gen = n_zones
    specs = [round(q_init + (i - (n_gen - 1) / 2) * step, 2) for i in range(n_gen)]
    specs = [max(0.05, s) for s in specs]
    return sorted(set(specs))

def apply_config(q_init_rated, h_spec_rated, h_inlet, d_pipe, g, nu, xi,
                 n_zones, max_iter, tolerance, damping, available_specs,
                 init_n=0.5, repl_n=0.5):
    """将参数写入dp模块全局变量"""
    dp.Q_INIT_RATED = q_init_rated
    dp.H_SPEC_RATED = h_spec_rated
    dp.H_INLET = h_inlet
    dp.D_PIPE = d_pipe
    dp.G = g
    dp.NU = nu
    dp.XI = xi
    dp.N_ZONES = n_zones
    dp.MAX_ITER = max_iter
    dp.TOLERANCE = tolerance
    dp.DAMPING = damping
    dp.A_PIPE = np.pi * d_pipe**2 / 4
    dp.FLOW_EXPONENT_INIT = init_n
    dp.FLOW_EXPONENT_REPL = repl_n
    dp.k_init = q_init_rated / h_spec_rated**init_n
    dp.SPEC_K = {s: s / h_spec_rated**repl_n for s in available_specs}

def compute_case(case_type, positions, target_tail_h, available_specs, init_n=0.5, repl_n=0.5, use_limit=False, skip_replacement=False):
    """完整管道计算（不写Excel），返回结果dict"""
    n = len(positions)

    if n == 0:
        raise ValueError('当前布局没有生成滴头，请检查管道长度和间距')
    if dp.N_ZONES > n:
        raise ValueError(f'分区数 {dp.N_ZONES} 不能大于滴头数 {n}，请减少分区数或增加管道长度')

    # ==== 初始求解（使用初始滴头的流量指数）====
    dp.FLOW_EXPONENT = init_n
    if use_limit:
        positions, k_arr, h_init, q_init, best_n = dp.find_limit_length(
            positions, dp.k_init)
    elif target_tail_h is not None:
        positions, k_arr, h_init, q_init, best_n = dp.find_nearest_tail(
            positions, target_tail_h, dp.k_init)
    else:
        k_arr = np.full(n, dp.k_init)
        h_init, q_init = dp.solve(positions, k_arr)

    n = len(positions)
    if dp.N_ZONES > n:
        raise ValueError(f'分区数 {dp.N_ZONES} 不能大于截取后的滴头数 {n}')
    init_solve_status = dp.LAST_SOLVE.copy()
    zone_ids_init, boundaries = dp.zones_by_integral(positions, h_init, dp.N_ZONES)
    q_target = dp.compute_q_target(q_init, zone_ids_init, dp.N_ZONES)
    zone_labels = dp.get_zone_labels(dp.N_ZONES)

    # 极限长度模式：只出初始状态，跳过替换迭代
    if skip_replacement:
        init_cu = dp.cu(q_init)
        return {
            'case_type': case_type,
            'positions': positions,
            'n_drippers': n,
            'pipe_length': positions[-1],
            'init': {'h': h_init, 'q': q_init, 'zone_ids': zone_ids_init},
            'repl': {'h': h_init, 'q': q_init, 'zone_ids': zone_ids_init, 'k_arr': k_arr},
            'boundaries': boundaries,
            'q_target': q_target,
            'zone_labels': zone_labels,
            'sweeps': [],
            'init_cu': init_cu,
            'repl_cu': init_cu,
            'cu_delta': 0.0,
            'init_q_range': (np.min(q_init), np.max(q_init)),
            'repl_q_range': (np.min(q_init), np.max(q_init)),
            'init_tail_h': h_init[-1],
            'repl_tail_h': h_init[-1],
            'final_specs': {},
            'spec_k': dp.SPEC_K.copy(),
            'k_init': dp.k_init,
            'init_n': init_n,
            'repl_n': repl_n,
            'config': {'nu': dp.NU, 'd_pipe': dp.D_PIPE, 'xi': dp.XI,
                       'h_inlet': dp.H_INLET, 'g': dp.G, 'A_pipe': dp.A_PIPE},
            'init_solve_status': init_solve_status,
        }

    # 迭代替换（切换到替换滴头的流量指数）
    dp.FLOW_EXPONENT = repl_n
    k_arr = k_arr.copy()
    h, q = h_init.copy(), q_init.copy()
    zone_ids = zone_ids_init.copy()
    sweeps_log = []

    for sweep in range(1, 31):
        k_need = np.array([q_target / max(h[i], 0.01)**dp.FLOW_EXPONENT for i in range(n)])
        zone_spec, k_need_stats = dp.select_zone_specs(
            zone_ids, k_need, dp.N_ZONES, dp.SPEC_K)

        changes = 0
        for i in range(n):
            if zone_ids[i] not in zone_spec:
                continue
            target_k = zone_spec[zone_ids[i]][1]
            if abs(k_arr[i] - target_k) < 1e-4:
                continue
            k_arr[i] = target_k
            changes += 1

        if changes > 0:
            h, q = dp.solve(positions, k_arr)
            zone_ids, _ = dp.zones_by_integral(positions, h, dp.N_ZONES)

        # 收集本轮详细统计
        after_stats = []
        for zid in range(dp.N_ZONES):
            mask = zone_ids == zid
            cnt = np.sum(mask)
            if cnt == 0:
                continue
            hz, qz = h[mask], q[mask]
            spec_name = f'{zone_spec[zid][0]}' if zid in zone_spec else '—'
            after_stats.append({
                'name': zone_labels[zid], 'cnt': cnt, 'spec': spec_name,
                'h_min': np.min(hz), 'h_max': np.max(hz), 'h_avg': np.mean(hz),
                'q_min': np.min(qz), 'q_max': np.max(qz), 'q_avg': np.mean(qz),
            })

        sweeps_log.append({
            'round': sweep, 'changes': changes, 'cu': float(dp.cu(q)),
            'solve_status': dp.LAST_SOLVE.copy(),
            'k_need_stats': k_need_stats,
            'after_stats': after_stats,
        })
        if changes == 0:
            break

    # 规格统计
    final_specs = Counter()
    for i in range(n):
        for s_name, s_k in dp.SPEC_K.items():
            if abs(k_arr[i] - s_k) < 1e-4:
                final_specs[s_name] += 1
                break

    init_cu = dp.cu(q_init)
    repl_cu = dp.cu(q)

    return {
        'case_type': case_type,
        'init_n': init_n,
        'repl_n': repl_n,
        'config': {'nu': dp.NU, 'd_pipe': dp.D_PIPE, 'xi': dp.XI,
                   'h_inlet': dp.H_INLET, 'g': dp.G, 'A_pipe': dp.A_PIPE},
        'init_solve_status': init_solve_status,
        'positions': positions,
        'n_drippers': n,
        'pipe_length': positions[-1],
        'init': {'h': h_init, 'q': q_init, 'zone_ids': zone_ids_init},
        'repl': {'h': h, 'q': q, 'zone_ids': zone_ids, 'k_arr': k_arr},
        'boundaries': boundaries,
        'q_target': q_target,
        'zone_labels': zone_labels,
        'sweeps': sweeps_log,
        'init_cu': init_cu,
        'repl_cu': repl_cu,
        'cu_delta': repl_cu - init_cu,
        'init_q_range': (np.min(q_init), np.max(q_init)),
        'repl_q_range': (np.min(q), np.max(q)),
        'init_tail_h': h_init[-1],
        'repl_tail_h': h[-1],
        'final_specs': final_specs,
        'spec_k': dp.SPEC_K.copy(),
        'k_init': dp.k_init,
    }

# ============ 导出辅助 ============

def _ff(Re):
    """沿程阻力系数"""
    return 64.0 / Re if Re < 2000 else 0.3164 / Re**0.25

def build_full_export_df(positions, q, k_vals, nu, d_pipe, xi, h_inlet, g, A_pipe, flow_exponent):
    """构建17列完整水力数据表（与原始递减流量计算.py格式一致）"""
    n = len(positions)
    q_m3s = q / 3.6e6
    total_q_m3s = np.sum(q_m3s)
    if not isinstance(k_vals, np.ndarray):
        k_vals = np.full(n, k_vals)

    rows = []
    h_out_prev = None
    for i in range(n):
        v_i = np.sum(q_m3s[i:]) / A_pipe
        hj_i = xi * v_i**2 / (2 * g)

        if i == 0:
            h_in = h_inlet
            hf_i = 0.0
        else:
            v_seg = np.sum(q_m3s[i:]) / A_pipe
            Re_seg = v_seg * d_pipe / nu
            lam_seg = _ff(Re_seg)
            dl = positions[i] - positions[i-1]
            hf_i = lam_seg * dl * v_seg**2 / (2 * g * d_pipe)
            h_in = h_out_prev - hf_i

        h_out_i = max(h_in - hj_i, 0.01)
        h_out_prev = h_out_i

        Re_i = v_i * d_pipe / nu
        lam_i = _ff(Re_i)
        q_rated_i = k_vals[i] * max(h_in, 0.01)**flow_exponent

        rows.append([
            i + 1, positions[i], nu, h_in, total_q_m3s, k_vals[i],
            q_rated_i, d_pipe, v_i, Re_i, lam_i, xi, hf_i, hj_i,
            h_out_i, h_out_i, q[i],
        ])

    headers = [
        '滴头数', '毛管长度L(m)', '运动粘度ν(m²/s)',
        '入口水头值(m)', '流量(m³/s)', '流量系数k',
        '额定流量qd(L/h)', '管径(m)', '流速v(m/s)',
        '雷诺数Re', '沿程阻力系数λ', '局部水头损失系数ξ',
        '沿程水头损失hf(m)', '局部水头损失hj(m)',
        '出口水头值(m)', '末端水头值(m)', '末端滴头流量(L/h)',
    ]
    df = pd.DataFrame(rows, columns=headers)
    return df

def build_sweeps_df(sweeps_log):
    """构建迭代记录表（含各区详细统计）"""
    rows = []
    for sw in sweeps_log:
        base = {'迭代轮次': sw['round'], '本轮替换数': sw['changes']}
        if sw.get('after_stats'):
            for st_item in sw['after_stats']:
                row = base.copy()
                row['分区'] = st_item['name']
                row['规格(L/h)'] = st_item['spec']
                row['滴头数量'] = st_item['cnt']
                row['水头范围(m)'] = f"{st_item['h_min']:.4f}~{st_item['h_max']:.4f}"
                row['水头均值(m)'] = f"{st_item['h_avg']:.4f}"
                row['流量范围(L/h)'] = f"{st_item['q_min']:.4f}~{st_item['q_max']:.4f}"
                row['流量均值(L/h)'] = f"{st_item['q_avg']:.4f}"
                rows.append(row)
        else:
            rows.append(base)
    return pd.DataFrame(rows)

# ============ 图表函数 ============

AXIS_STYLE = dict(showline=True, mirror=True, linecolor='#444', linewidth=1.2,
                  ticks='inside', tickcolor='#444', ticklen=4, tickwidth=1)

def _axis_style(fig):
    """统一坐标轴线样式"""
    fig.update_xaxes(**AXIS_STYLE)
    fig.update_yaxes(**AXIS_STYLE)
    return fig


def make_flow_chart(result, show_init=True, show_repl=True):
    fig = go.Figure()
    if show_init:
        fig.add_trace(go.Scatter(x=result['positions'], y=result['init']['q'],
            mode='lines', name='初始状态', line=dict(color=C_INIT, width=2.5)))
    if show_repl:
        fig.add_trace(go.Scatter(x=result['positions'], y=result['repl']['q'],
            mode='lines', name='替换后', line=dict(color=C_REPL, width=2.5)))
    fig.update_layout(xaxis_title='距离 (m)', yaxis_title='滴头流量 (L/h)',
        template='plotly_white', hovermode='x',
        legend=dict(orientation='h', yanchor='bottom', y=1.02))
    fig.update_xaxes(tickformat='.0f')
    return _axis_style(fig)

def make_head_chart(result):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=result['positions'], y=result['init']['h'],
        mode='lines', name='初始水头', line=dict(color=C_INIT, width=2.5)))
    fig.add_trace(go.Scatter(x=result['positions'], y=result['repl']['h'],
        mode='lines', name='替换后水头', line=dict(color=C_REPL, width=2.5)))
    fig.update_layout(xaxis_title='距离 (m)', yaxis_title='出口水头 (m)',
        template='plotly_white', hovermode='x',
        legend=dict(orientation='h', yanchor='bottom', y=1.02))
    fig.update_xaxes(tickformat='.0f')
    return _axis_style(fig)

def make_zone_chart(result):
    fig = go.Figure()
    pos = result['positions']
    q_repl = result['repl']['q']
    zids = result['repl']['zone_ids']
    n_zones = len(set(zids))
    labels = result['zone_labels']
    for zid in range(n_zones):
        mask = zids == zid
        if np.sum(mask) == 0:
            continue
        color = ZONE_COLORS[zid % len(ZONE_COLORS)]
        fig.add_trace(go.Scatter(x=pos[mask], y=q_repl[mask],
            mode='markers', name=labels[zid],
            marker=dict(color=color, size=5)))
    fig.update_layout(xaxis_title='距离 (m)', yaxis_title='滴头流量 (L/h)',
        template='plotly_white', hovermode='x',
        legend=dict(orientation='h', yanchor='bottom', y=1.02))
    fig.update_xaxes(tickformat='.0f')
    return _axis_style(fig)

def make_histogram(result):
    fig = go.Figure()
    fig.add_trace(go.Histogram(x=result['init']['q'], name='初始状态',
        marker_color=C_HIST_INIT, opacity=0.7, nbinsx=40))
    fig.add_trace(go.Histogram(x=result['repl']['q'], name='替换后',
        marker_color=C_HIST_REPL, opacity=0.7, nbinsx=40))
    fig.update_layout(xaxis_title='滴头流量 (L/h)', yaxis_title='频数',
        template='plotly_white', barmode='overlay',
        legend=dict(orientation='h', yanchor='bottom', y=1.02))
    return _axis_style(fig)

def make_pie_chart(result):
    specs = result['final_specs']
    if not specs:
        fig = go.Figure()
        fig.update_layout(template='plotly_white')
        return fig
    labels = [f'{s} L/h' for s in sorted(specs.keys())]
    values = [specs[s] for s in sorted(specs.keys())]
    colors = ZONE_COLORS[:len(labels)]
    fig = go.Figure(go.Pie(labels=labels, values=values,
        marker=dict(colors=colors), hole=0.4))
    fig.update_layout(template='plotly_white')
    return fig

def make_cu_comparison(eq_result, var_result):
    categories = ['等间距\n初始', '等间距\n替换', '变间距\n初始', '变间距\n替换']
    cu_vals = [eq_result['init_cu'], eq_result['repl_cu'],
               var_result['init_cu'], var_result['repl_cu']]
    colors = [C_EQ_INIT, C_EQ_REPL, C_VAR_INIT, C_VAR_REPL]
    fig = go.Figure(go.Bar(x=categories, y=cu_vals, marker_color=colors,
        text=[f'{v:.2f}%' for v in cu_vals], textposition='outside'))
    fig.update_layout(yaxis_title='Cu (%)', template='plotly_white',
        yaxis=dict(range=[min(cu_vals)-2, max(cu_vals)+2]))
    return _axis_style(fig)

def make_four_curve(eq_result, var_result):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=eq_result['positions'], y=eq_result['init']['q'],
        mode='lines', name='等间距 初始', line=dict(color=C_EQ_INIT, dash='dash', width=2)))
    fig.add_trace(go.Scatter(x=eq_result['positions'], y=eq_result['repl']['q'],
        mode='lines', name='等间距 替换后', line=dict(color=C_EQ_REPL, width=2.5)))
    fig.add_trace(go.Scatter(x=var_result['positions'], y=var_result['init']['q'],
        mode='lines', name='变间距 初始', line=dict(color=C_VAR_INIT, dash='dash', width=2)))
    fig.add_trace(go.Scatter(x=var_result['positions'], y=var_result['repl']['q'],
        mode='lines', name='变间距 替换后', line=dict(color=C_VAR_REPL, width=2.5)))
    fig.update_layout(xaxis_title='距离 (m)', yaxis_title='滴头流量 (L/h)',
        template='plotly_white', hovermode='x',
        legend=dict(orientation='h', yanchor='bottom', y=1.02))
    fig.update_xaxes(tickformat='.0f')
    return _axis_style(fig)

def make_convergence_chart(result):
    cu_history = [result['init_cu']] + [sw['cu'] for sw in result['sweeps']]
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=list(range(len(cu_history))), y=cu_history,
        mode='lines+markers', name='Cu', line=dict(color=C_INIT, width=2)))
    fig.update_layout(xaxis_title='迭代轮次', yaxis_title='Cu (%)',
        template='plotly_white')
    return _axis_style(fig)

def make_length_search_chart(lengths, cu_vals, best_idx):
    fig = go.Figure()
    fig.add_trace(go.Scatter(x=lengths, y=cu_vals,
        mode='lines+markers', name='Cu', line=dict(color=C_INIT, width=2),
        marker=dict(size=8)))
    fig.add_trace(go.Scatter(x=[lengths[best_idx]], y=[cu_vals[best_idx]],
        mode='markers', name=f'最优={lengths[best_idx]:.1f}m',
        marker=dict(color=C_HIGHLIGHT, size=16, symbol='star')))
    fig.update_layout(xaxis_title='铺设长度 (m)', yaxis_title='替换后Cu (%)',
        template='plotly_white')
    return _axis_style(fig)

def make_heatmap(grid_data, metric='repl_cu'):
    q_vals = sorted(set(d['q_init'] for d in grid_data))
    n_vals = sorted(set(d['n_zones'] for d in grid_data))
    z = np.zeros((len(q_vals), len(n_vals)))
    for d in grid_data:
        qi = q_vals.index(d['q_init'])
        ni = n_vals.index(d['n_zones'])
        z[qi, ni] = d.get(metric, 0)

    if metric == 'cu_delta':
        colorscale = HEATMAP_DELTA_SCALE
        cbar_title = 'Cu 增量 (%)'
    else:
        colorscale = HEATMAP_CU_SCALE
        cbar_title = '平均 Cu (%)' if metric == 'repl_cu' else metric

    fig = go.Figure(go.Heatmap(z=z, x=n_vals, y=q_vals,
        colorscale=colorscale,
        text=[[f'{v:.2f}' for v in row] for row in z],
        texttemplate='%{text}',
        textfont=dict(size=13, color='#333'),
        hovertemplate='Q=%{y} L/h, 分区=%{x}<br>值=%{z:.2f}<extra></extra>',
        colorbar=dict(title=dict(text=cbar_title, side='right'), tickformat='.1f')))
    fig.update_layout(
        xaxis_title='分区数', yaxis_title='初始额定流量 (L/h)',
        template='plotly_white',
        xaxis=dict(dtick=1))
    return _axis_style(fig)

def make_cu_delta_bar(grid_data):
    labels = [f"Q={d['q_init']},N={d['n_zones']}" for d in grid_data]
    deltas = [d['cu_delta'] for d in grid_data]
    colors = [C_REPL if d > 0 else C_INIT for d in deltas]
    fig = go.Figure(go.Bar(x=labels, y=deltas, marker_color=colors,
        text=[f'{d:+.2f}%' for d in deltas], textposition='outside'))
    fig.update_layout(yaxis_title='Cu 增量 (%)', template='plotly_white')
    return _axis_style(fig)

# ============ 指标卡片 ============
def show_metrics(result, prefix=''):
    cols = st.columns(4)
    cols[0].metric(f'{prefix}Cu (均匀度)', f"{result['repl_cu']:.2f}%",
                   delta=f"{result['cu_delta']:+.2f}%")
    cols[1].metric(f'{prefix}平均流量', f"{np.mean(result['repl']['q']):.4f} L/h")
    cols[2].metric(f'{prefix}流量范围',
                   f"{result['repl_q_range'][0]:.4f}~{result['repl_q_range'][1]:.4f}")
    cols[3].metric(f'{prefix}末端水头', f"{result['repl_tail_h']:.4f} m")
    # 显示替换后各规格数量
    if result['final_specs']:
        specs_str = ', '.join(f'{s}L/h×{c}' for s, c in sorted(result['final_specs'].items()))
        st.caption(f'{prefix}替换后规格分布: {specs_str}')

# ============ 侧边栏 ============
def build_sidebar():
    st.sidebar.title('参数配置')

    mode = st.sidebar.radio('运行模式', ['单次', '批量扫描'], key='mode',
        help='单次: 一组参数; 批量: 多组 Q×分区 网格扫描')
    analysis_mode = st.sidebar.radio('分析模式',
        ['只看初始状态', '初始+替换对比'], key='analysis_mode',
        help='只看初始: 仅计算统一规格下的水力状态，不执行分区替换')
    do_repl = (analysis_mode == '初始+替换对比')

    # --- 滴头参数 ---
    with st.sidebar.expander('滴头参数', expanded=True):
        st.caption('初始（基准方案）')
        c1, c2 = st.columns(2)
        init_type = c1.radio('类型', ['普通滴头', '压力补偿滴头'],
            horizontal=True, key='init_type', label_visibility='collapsed')
        if init_type == '压力补偿滴头':
            init_n = c2.slider('初始 n', 0.0, 0.5, 0.1, 0.01, key='init_n')
        else:
            c2.caption('n=0.5')
            init_n = 0.5
        q_label = '额定流量' if mode == '单次' else '额定流量 (逗号分隔)'
        q_init_rated = st.text_input(q_label, value='2.5' if mode == '单次' else '2.0, 2.5, 3.0',
            key='q_init_rated')
        h_spec_rated = st.number_input('规格额定水头 (m)', value=10.0, format='%.2f', key='h_spec_rated')

        if do_repl:
            st.markdown('---')
            st.caption('替换（优化方案）')
            c1, c2 = st.columns(2)
            repl_type = c1.radio('类型', ['普通滴头', '压力补偿滴头'],
                horizontal=True, key='repl_type', label_visibility='collapsed')
            if repl_type == '压力补偿滴头':
                repl_n = c2.slider('替换 n', 0.0, 0.5, 0.1, 0.01, key='repl_n')
            else:
                c2.caption('n=0.5')
                repl_n = 0.5
            spec_mode = st.radio('规格', ['自动生成', '手动指定'],
                horizontal=True, key='spec_mode')
            if spec_mode == '手动指定':
                specs_str = st.text_input('可用规格 (逗号分隔)', value='2.0, 2.5, 3.0',
                    key='specs_str')
            else:
                specs_str = ''
                try:
                    qs = [float(x.strip()) for x in q_init_rated.replace('，', ',').split(',') if x.strip()]
                    nzs = [int(x.strip()) for x in ('3' if mode == '单次' else '2, 3, 4').replace('，', ',').split(',') if x.strip()]
                    actual_nzs = [int(x.strip()) for x in str(st.session_state.get('n_zones', nzs[len(nzs)//2])).replace('，', ',').split(',') if x.strip()]
                    preview = auto_specs(qs[0], actual_nzs[len(actual_nzs)//2])
                    st.caption(f'自动: {", ".join(f"{s:.2f}" for s in preview)} L/h')
                except (ValueError, IndexError):
                    pass
        else:
            repl_type = '普通滴头'
            repl_n = 0.5
            spec_mode = 'auto'
            specs_str = ''

    # --- 管道 & 分区 ---
    with st.sidebar.expander('管道 & 分区', expanded=True):
        c1, c2 = st.columns(2)
        d_pipe = c1.number_input('管径 (m)', value=0.016, format='%.4f', key='d_pipe')
        h_inlet = c2.number_input('入口水头 (m)', value=10.0, format='%.4f', key='h_inlet')
        n_label = '分区数' if mode == '单次' else '分区数 (逗号分隔)'
        n_zones = st.text_input(n_label, value='3' if mode == '单次' else '2, 3, 4',
            key='n_zones')

    # --- 管道布局 ---
    with st.sidebar.expander('管道布局', expanded=True):
        # ==== 等间距 ====
        with st.container(border=True):
            eq_enabled = st.checkbox('等间距', value=True, key='eq_enabled')
            eq_disabled = not eq_enabled
            eq_spacing = st.number_input('间距 (m)', value=0.3, format='%.2f',
                key='eq_spacing', disabled=eq_disabled)
            eq_trim_options = ['指定长度', '目标水头']
            if do_repl:
                eq_trim_options.append('Cu达标最长')
            eq_trim_labels = {'指定长度': 'full', '目标水头': 'target', 'Cu达标最长': 'cu_best'}
            eq_trim_mode = eq_trim_labels[
                st.radio('截取模式', eq_trim_options, horizontal=True,
                    key='eq_trim_mode_radio', disabled=eq_disabled)]
            eq_target = None
            eq_length = 120.0
            eq_cu_threshold = 90.0
            eq_scan_min = eq_scan_max = eq_scan_step = None
            if eq_trim_mode == 'full':
                eq_length = st.number_input('铺设长度 (m)', value=120.0, format='%.1f',
                    key='eq_length', disabled=eq_disabled)
            elif eq_trim_mode == 'target':
                eq_target = st.number_input('目标末端水头 (m)', value=5.0, format='%.2f',
                    min_value=0.01, key='eq_target_num', disabled=eq_disabled,
                    help='填极小值(如0.01)等效于搜索极限长度')
            elif eq_trim_mode == 'cu_best':
                eq_cu_threshold = st.number_input('Cu 底线 (%)', value=90.0, format='%.1f',
                    key='eq_cu_threshold', disabled=eq_disabled,
                    help='Cu在此值之上取最长铺设长度')
                c1, c2, c3 = st.columns(3)
                eq_scan_min = c1.number_input('最小 (m)', value=50.0, format='%.0f',
                    key='eq_scan_min_cu', disabled=eq_disabled)
                eq_scan_max = c2.number_input('最大 (m)', value=200.0, format='%.0f',
                    key='eq_scan_max_cu', disabled=eq_disabled)
                eq_scan_step = c3.number_input('步长 (m)', value=10.0, format='%.0f',
                    key='eq_scan_step_cu', disabled=eq_disabled)

        st.markdown('')
        # ==== 变间距 ====
        with st.container(border=True):
            var_enabled = st.checkbox('变间距', value=True, key='var_enabled')
            var_disabled = not var_enabled
            c1, c2 = st.columns(2)
            var_first_seg = c1.number_input('首段间距 (m)', value=0.2, format='%.2f',
                key='var_first_seg', disabled=var_disabled)
            var_cycle_str = c2.text_input('循环间距', value='1.6, 0.4',
                key='var_cycle_str', disabled=var_disabled,
                help='逗号分隔, 如 1.6, 0.4')
            var_trim_map = {'指定长度': 'full', '目标水头': 'target'}
            var_trim_mode = var_trim_map[
                st.radio('截取模式', ['指定长度', '目标水头'], horizontal=True,
                    key='var_trim_mode_radio', disabled=var_disabled)]
            var_length = 500.0
            var_target_mode = 'eq'
            var_target_val = None
            if var_trim_mode == 'full':
                var_length = st.number_input('生成总长度 (m)', value=500.0, format='%.1f',
                    key='var_length', disabled=var_disabled)
            elif var_trim_mode == 'target':
                vm = st.radio('目标来源', ['匹配等间距', '自定义数值'], horizontal=True,
                    key='var_target_mode_radio', disabled=var_disabled)
                if vm == '自定义数值':
                    var_target_mode = 'custom'
                    var_target_val = st.number_input('目标值 (m)', value=5.0, format='%.4f',
                        min_value=0.01, key='var_target_val_num', disabled=var_disabled)
                else:
                    var_target_mode = 'eq'
                    var_target_val = 'eq'

    # --- 高级参数 ---
    with st.sidebar.expander('高级参数', expanded=False):
        c1, c2, c3 = st.columns(3)
        g = c1.number_input('G (m/s²)', value=9.81, format='%.3f', key='g')
        nu = c2.number_input('ν (m/s²)', value=1.004e-6, format='%.4e', key='nu')
        xi = c3.slider('ξ', 0.0, 1.0, 0.45, 0.01, key='xi')
        c1, c2, c3 = st.columns(3)
        max_iter = c1.number_input('最大迭代', value=200, key='max_iter')
        tolerance = c2.number_input('收敛容差', value=1e-8, format='%.1e', key='tolerance')
        damping = c3.slider('阻尼', 0.1, 1.0, 0.5, 0.05, key='damping')

    st.sidebar.markdown('---')
    run_btn = st.sidebar.button('运行计算', type='primary', width='stretch', key='run_btn')

    return {
        'mode': mode,
        'analysis_mode': analysis_mode, 'do_repl': do_repl,
        'g': g, 'nu': nu, 'xi': xi,
        'd_pipe': d_pipe, 'h_inlet': h_inlet,
        'q_init_rated': q_init_rated, 'h_spec_rated': h_spec_rated,
        'init_type': init_type, 'init_n': init_n,
        'repl_type': repl_type, 'repl_n': repl_n,
        'spec_mode': spec_mode, 'specs_str': specs_str,
        'n_zones': n_zones,
        'max_iter': max_iter, 'tolerance': tolerance, 'damping': damping,
        'eq_enabled': eq_enabled,
        'eq_spacing': eq_spacing,
        'eq_trim_mode': eq_trim_mode,
        'eq_length': eq_length,
        'eq_target': eq_target,
        'eq_cu_threshold': eq_cu_threshold,
        'eq_scan_min': eq_scan_min, 'eq_scan_max': eq_scan_max, 'eq_scan_step': eq_scan_step,
        'var_enabled': var_enabled,
        'var_first_seg': var_first_seg, 'var_cycle_str': var_cycle_str,
        'var_trim_mode': var_trim_mode,
        'var_length': var_length,
        'var_target_mode': var_target_mode, 'var_target_val': var_target_val,
        'run_btn': run_btn,
    }


# ============ 单次计算 ============
def run_single_compute(args, specs_list, var_cycle):
    """执行单次模式计算，返回结果存入session_state"""
    q_val = float(args['q_init_rated'])
    n_zones_val = int(args['n_zones'])
    init_n = args['init_n']
    repl_n = args['repl_n']
    do_repl = args['do_repl']
    eq_enabled = args['eq_enabled']
    var_enabled = args['var_enabled']

    apply_config(q_val, args['h_spec_rated'], args['h_inlet'],
                 args['d_pipe'], args['g'], args['nu'], args['xi'],
                 n_zones_val, args['max_iter'], args['tolerance'],
                 args['damping'], specs_list,
                 init_n=init_n, repl_n=repl_n)

    length_search_data = None
    eq_result = None
    var_result = None
    eq_tail_h = None
    config_snapshot = {'nu': dp.NU, 'd_pipe': dp.D_PIPE, 'xi': dp.XI,
                       'h_inlet': dp.H_INLET, 'g': dp.G, 'A_pipe': dp.A_PIPE}

    # 等间距
    if eq_enabled:
        eq_trim = args['eq_trim_mode']
        if eq_trim == 'cu_best':
            lengths = np.arange(args['eq_scan_min'],
                                args['eq_scan_max'] + 1e-9,
                                args['eq_scan_step'])
            cu_vals = []
            for L in lengths:
                pos = dp.generate_eq_positions(args['eq_spacing'], L)
                r = compute_case('eq', pos, None, specs_list,
                               init_n=init_n, repl_n=repl_n, skip_replacement=not do_repl)
                cu_vals.append(r['repl_cu'])
            # Cu 达标下的最长
            threshold = args['eq_cu_threshold']
            best_idx = None
            best_L = None
            for i, L in enumerate(lengths):
                if cu_vals[i] >= threshold:
                    best_idx = i
                    best_L = L
            if best_idx is None:
                st.warning(f'未找到 Cu ≥ {threshold}% 的长度，取 Cu 最高者')
                best_idx = int(np.argmax(cu_vals))
                best_L = lengths[best_idx]
            pos = dp.generate_eq_positions(args['eq_spacing'], best_L)
            eq_result = compute_case('eq', pos, None, specs_list,
                                   init_n=init_n, repl_n=repl_n, skip_replacement=not do_repl)
            length_search_data = {'lengths': list(lengths), 'cu_vals': cu_vals,
                                  'best_idx': best_idx, 'threshold': threshold}
        elif eq_trim == 'target':
            eq_target = args['eq_target']
            if eq_target is not None and eq_target >= dp.H_INLET:
                raise ValueError(f'目标末端水头必须小于入口水头 {dp.H_INLET:g} m')
            if eq_target is not None and eq_target < 0.5:
                # 极小目标 → 自动走极限长度二分搜索
                eq_long_l = max(args['eq_length'] * 10, 3000)
                pos_long = dp.generate_eq_positions(args['eq_spacing'], eq_long_l)
                eq_result = compute_case('eq', pos_long, eq_target, specs_list,
                                       init_n=init_n, repl_n=repl_n, use_limit=True,
                                       skip_replacement=not do_repl)
            elif eq_target is not None:
                # 正常目标水头 → eq_length 作为参考长度，find_nearest_tail 扫描
                pos = dp.generate_eq_positions(args['eq_spacing'], args['eq_length'])
                eq_result = compute_case('eq', pos, eq_target, specs_list,
                                       init_n=init_n, repl_n=repl_n, skip_replacement=not do_repl)
            else:
                # 无目标 → 完整长度
                pos = dp.generate_eq_positions(args['eq_spacing'], args['eq_length'])
                eq_result = compute_case('eq', pos, None, specs_list,
                                       init_n=init_n, repl_n=repl_n, skip_replacement=not do_repl)
        else:  # 'full'
            pos = dp.generate_eq_positions(args['eq_spacing'], args['eq_length'])
            eq_result = compute_case('eq', pos, None, specs_list,
                                   init_n=init_n, repl_n=repl_n, skip_replacement=not do_repl)

        eq_tail_h = eq_result['init_tail_h'] if eq_result else None

    # 变间距
    if var_enabled:
        var_trim = args['var_trim_mode']
        if var_trim == 'target':
            if args['var_target_val'] == 'eq':
                var_target = eq_tail_h
            else:
                var_target = args['var_target_val']
            if var_target is not None and isinstance(var_target, (int, float)) and var_target < 0.5:
                var_long_l = max(args['var_length'] * 5, 3000)
                var_pos = dp.generate_var_positions(args['var_first_seg'], var_cycle, var_long_l)
                var_result = compute_case('var', var_pos, var_target, specs_list,
                                        init_n=init_n, repl_n=repl_n, use_limit=True,
                                        skip_replacement=not do_repl)
            elif var_target is not None:
                var_pos = dp.generate_var_positions(args['var_first_seg'], var_cycle, args['var_length'])
                var_result = compute_case('var', var_pos, var_target, specs_list,
                                        init_n=init_n, repl_n=repl_n, skip_replacement=not do_repl)
            else:
                var_pos = dp.generate_var_positions(args['var_first_seg'], var_cycle, args['var_length'])
                var_result = compute_case('var', var_pos, None, specs_list,
                                        init_n=init_n, repl_n=repl_n, skip_replacement=not do_repl)
        else:  # 'full'
            var_pos = dp.generate_var_positions(args['var_first_seg'], var_cycle, args['var_length'])
            var_result = compute_case('var', var_pos, None, specs_list,
                                    init_n=init_n, repl_n=repl_n, skip_replacement=not do_repl)

    # 预生成Excel（渲染时直接取，避免下载时重算导致卡死）
    excel_bytes = _build_excel_bytes(eq_result, var_result)

    st.session_state['single_result'] = {
        'eq': eq_result, 'var': var_result,
        'length_search': length_search_data,
        'specs_list': specs_list,
        'n_zones': n_zones_val,
        'config': config_snapshot,
        'init_type': args['init_type'], 'init_n': init_n,
        'repl_type': args['repl_type'], 'repl_n': repl_n,
        'repl_enabled': do_repl,
        'eq_enabled': eq_enabled, 'var_enabled': var_enabled,
        '_excel_bytes': excel_bytes,
    }


def _build_excel_bytes(eq_result, var_result):
    """预生成Excel字节流，避免渲染时重复计算"""
    output = BytesIO()
    with pd.ExcelWriter(output, engine='openpyxl') as writer:
        for label, r in [('等间距', eq_result), ('变间距', var_result)]:
            if r is None:
                continue
            df = build_full_export_df(r['positions'], r['init']['q'],
                                      r['k_init'], **r['config'], flow_exponent=r['init_n'])
            df.to_excel(writer, sheet_name=f'{label}-初始', index=False)
            df = build_full_export_df(r['positions'], r['repl']['q'],
                                      r['repl']['k_arr'], **r['config'], flow_exponent=r['repl_n'])
            df.to_excel(writer, sheet_name=f'{label}-替换', index=False)
            build_sweeps_df(r['sweeps']).to_excel(writer, sheet_name=f'{label}-迭代记录', index=False)
    return output.getvalue()


def run_batch_compute(args, specs_list, var_cycle, q_list, n_list):
    """执行批量模式计算，返回结果存入session_state"""
    all_results = []
    total = len(q_list) * len(n_list)
    progress_bar = st.progress(0)
    status_text = st.empty()
    idx = 0
    init_n = args['init_n']
    repl_n = args['repl_n']
    do_repl = args['do_repl']

    for q_val in q_list:
        for n_val in n_list:
            idx += 1
            status_text.text(f'正在计算 Q={q_val} L/h, 分区数={n_val}  ({idx}/{total})')

            cur_specs = auto_specs(q_val, n_val) if args['spec_mode'] == 'auto' else specs_list
            apply_config(q_val, args['h_spec_rated'], args['h_inlet'],
                         args['d_pipe'], args['g'], args['nu'], args['xi'],
                         n_val, args['max_iter'], args['tolerance'],
                         args['damping'], cur_specs,
                         init_n=init_n, repl_n=repl_n)

            eq_r = None
            if args['eq_enabled']:
                eq_pos = dp.generate_eq_positions(args['eq_spacing'], args['eq_length'])
                eq_target = args['eq_target'] if args['eq_trim_mode'] == 'target' else None
                eq_use_limit = eq_target is not None and eq_target < 0.5
                if eq_use_limit:
                    eq_pos = dp.generate_eq_positions(args['eq_spacing'], max(args['eq_length'] * 10, 3000))
                eq_r = compute_case('eq', eq_pos, eq_target, cur_specs,
                                  init_n=init_n, repl_n=repl_n, use_limit=eq_use_limit,
                                  skip_replacement=not do_repl)

            var_r = None
            if args['var_enabled']:
                var_pos = dp.generate_var_positions(args['var_first_seg'], var_cycle,
                                                    args['var_length'])
                var_target = None
                var_use_limit = False
                if args['var_trim_mode'] == 'target':
                    var_target = eq_r['init_tail_h'] if args['var_target_val'] == 'eq' else args['var_target_val']
                    if var_target is None:
                        raise ValueError('变间距匹配等间距时，请同时启用等间距计算')
                    var_use_limit = var_target < 0.5
                    if var_use_limit:
                        var_pos = dp.generate_var_positions(args['var_first_seg'], var_cycle,
                                                            max(args['var_length'] * 5, 3000))
                var_r = compute_case('var', var_pos, var_target, cur_specs,
                                   init_n=init_n, repl_n=repl_n, use_limit=var_use_limit,
                                   skip_replacement=not do_repl)
            active = [r for r in (eq_r, var_r) if r is not None]

            all_results.append({
                'q_init': q_val, 'n_zones': n_val, 'specs': cur_specs,
                'eq': eq_r, 'var': var_r,
                'repl_cu': sum(r['repl_cu'] for r in active) / len(active),
                'cu_delta': sum(r['cu_delta'] for r in active) / len(active),
            })

            progress_bar.progress(idx / total)

    status_text.text('计算完成!')
    st.session_state['batch_result'] = all_results
    progress_bar.empty()
    status_text.empty()


# ============ 渲染函数 ============
def _render_result_tab(result, label, show_repl=True):
    """渲染单个结果的图表"""
    st.subheader(f"{label} {result['pipe_length']:.1f}m, {result['n_drippers']}个滴头")
    solve_status = (result.get('sweeps') or [{}])[-1].get('solve_status') or result.get('init_solve_status', {})
    if not solve_status.get('converged', True):
        st.warning(f'{label}水力求解达到最大迭代次数仍未收敛，结果仅供参考；请调整阻尼或增加迭代次数。')
    if solve_status.get('pressure_limited'):
        st.warning(f'{label}部分节点水头已降至 0.01 m 下限，当前长度下结果可能失真。')
    if show_repl and result['sweeps']:
        show_metrics(result, f'{label} ')
    else:
        c1, c2, c3, c4 = st.columns(4)
        c1.metric(f'{label}Cu (均匀度)', f"{result['init_cu']:.2f}%")
        c2.metric(f'{label}平均流量', f"{np.mean(result['init']['q']):.4f} L/h")
        c3.metric(f'{label}流量范围',
                  f"{result['init_q_range'][0]:.4f}~{result['init_q_range'][1]:.4f}")
        c4.metric(f'{label}末端水头', f"{result['init_tail_h']:.4f} m")

    c1, c2 = st.columns(2)
    with c1:
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=result['positions'], y=result['init']['q'],
            mode='lines', name='初始状态', line=dict(color=C_INIT, width=2.5)))
        if show_repl and result['sweeps']:
            fig.add_trace(go.Scatter(x=result['positions'], y=result['repl']['q'],
                mode='lines', name='替换后', line=dict(color=C_REPL, width=2.5)))
        fig.update_layout(xaxis_title='距离 (m)', yaxis_title='滴头流量 (L/h)',
            template='plotly_white', hovermode='x',
            legend=dict(orientation='h', yanchor='bottom', y=1.02))
        fig.update_xaxes(**AXIS_STYLE, tickformat='.0f')
        fig.update_yaxes(**AXIS_STYLE)
        st.plotly_chart(fig, width='stretch', key=f'{label}_flow')
    with c2:
        fig = go.Figure()
        fig.add_trace(go.Scatter(x=result['positions'], y=result['init']['h'],
            mode='lines', name='初始水头', line=dict(color=C_INIT, width=2.5)))
        if show_repl and result['sweeps']:
            fig.add_trace(go.Scatter(x=result['positions'], y=result['repl']['h'],
                mode='lines', name='替换后水头', line=dict(color=C_REPL, width=2.5)))
        fig.update_layout(xaxis_title='距离 (m)', yaxis_title='出口水头 (m)',
            template='plotly_white', hovermode='x',
            legend=dict(orientation='h', yanchor='bottom', y=1.02))
        fig.update_xaxes(**AXIS_STYLE, tickformat='.0f')
        fig.update_yaxes(**AXIS_STYLE)
        st.plotly_chart(fig, width='stretch', key=f'{label}_head')
    c3, c4 = st.columns(2)
    with c3:
        st.plotly_chart(make_zone_chart(result), width='stretch', key=f'{label}_zone')
    with c4:
        st.plotly_chart(make_histogram(result), width='stretch', key=f'{label}_hist')
    if show_repl and result['sweeps']:
        c5, c6 = st.columns(2)
        with c5:
            st.plotly_chart(make_pie_chart(result), width='stretch', key=f'{label}_pie')
        with c6:
            st.plotly_chart(make_convergence_chart(result), width='stretch', key=f'{label}_conv')
        # 迭代详情表
        with st.expander(f'{label} - 迭代详情', expanded=False):
            for sw in result['sweeps']:
                st.markdown(f"**第 {sw['round']} 轮** — 替换 {sw['changes']} 个滴头")
                if sw.get('after_stats'):
                    tbl = []
                    for z in sw['after_stats']:
                        tbl.append({
                            '分区': z['name'], '规格': z['spec'], '数量': z['cnt'],
                            '水头均值(m)': f"{z['h_avg']:.4f}",
                            '水头范围(m)': f"{z['h_min']:.4f}~{z['h_max']:.4f}",
                            '流量均值(L/h)': f"{z['q_avg']:.4f}",
                            '流量范围(L/h)': f"{z['q_min']:.4f}~{z['q_max']:.4f}",
                        })
                    st.dataframe(tbl, width='stretch', hide_index=True)
                if sw.get('k_need_stats'):
                    st.caption('k_需分析:')
                    ktbl = []
                    for ks in sw['k_need_stats']:
                        ktbl.append({
                            '区': ks['name'], '数量': ks['cnt'],
                            'k_需均值': f"{ks['kn_avg']:.4f}",
                            'k_需范围': f"{ks['kn_min']:.4f}~{ks['kn_max']:.4f}",
                            '选定规格': f"{ks['winner']} L/h",
                            '原选择分布': ks['choice_dist'],
                        })
                    st.dataframe(ktbl, width='stretch', hide_index=True)


def render_single_tabs(current_do_repl=False):
    """渲染单次模式的结果 — current_do_repl 取当前侧边栏设置，非存入结果时的旧值"""
    data = st.session_state['single_result']
    eq_result = data['eq']
    var_result = data['var']
    repl_enabled = data.get('repl_enabled', False)
    # 侧边栏分析与计算结果不一致时，提示重新运行
    need_rerun = current_do_repl != data.get('repl_enabled', False)

    # --- 状态提示 ---
    if need_rerun:
        st.warning('当前分析模式与已保存结果不一致，请点击“运行计算”重新计算。')
    if data['length_search']:
        d = data['length_search']
        st.success(f"最优长度: {d['lengths'][d['best_idx']]:.1f}m, Cu={d['cu_vals'][d['best_idx']]:.2f}%")
        st.plotly_chart(make_length_search_chart(d['lengths'], d['cu_vals'], d['best_idx']),
                        width='stretch', key='length_search')
    elif repl_enabled:
        st.success(f"计算完成! 初始={data['init_type']}(n={data['init_n']}) → 替换={data['repl_type']}(n={data['repl_n']})")
    else:
        st.success(f"初始状态计算完成! 初始={data['init_type']}(n={data['init_n']})，未执行替换")

    # --- 动态Tab ---
    has_eq = eq_result is not None
    has_var = var_result is not None
    both = has_eq and has_var

    tab_labels = []
    if has_eq:
        tab_labels.append('等间距结果')
    if has_var:
        tab_labels.append('变间距结果')
    if both:
        tab_labels.append('汇总对比')

    if not tab_labels:
        st.warning('没有启用的计算模式，请在侧边栏勾选至少一项。')
        return

    tabs = st.tabs(tab_labels)
    tab_idx = 0

    if has_eq:
        with tabs[tab_idx]:
            _render_result_tab(eq_result, '等间距', show_repl=repl_enabled)
        tab_idx += 1

    if has_var:
        with tabs[tab_idx]:
            _render_result_tab(var_result, '变间距', show_repl=repl_enabled)
        tab_idx += 1

    if both:
        with tabs[tab_idx]:
            st.subheader('Cu 对比')
            st.plotly_chart(make_cu_comparison(eq_result, var_result),
                            width='stretch', key='tab3_cu_comp')
            st.subheader('流量分布对比')
            show_overlay = st.checkbox('叠加到一张图', value=False, key='overlay_checkbox')
            if show_overlay:
                st.plotly_chart(make_four_curve(eq_result, var_result),
                                width='stretch', key='tab3_four_curve')
            else:
                c1, c2 = st.columns(2)
                with c1:
                    st.caption('等间距')
                    st.plotly_chart(make_flow_chart(eq_result), width='stretch', key='tab3_eq_flow')
                with c2:
                    st.caption('变间距')
                    st.plotly_chart(make_flow_chart(var_result), width='stretch', key='tab3_var_flow')

    # --- 导出 ---
    st.markdown('---')
    st.subheader('导出数据')
    st.download_button('下载结果 Excel', data=data['_excel_bytes'],
                       file_name='滴头水力计算结果.xlsx',
                       mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                       key='single_download')


def render_batch_tabs():
    """渲染批量模式的四个Tab"""
    all_results = st.session_state['batch_result']

    st.success(f'批量计算完成! 共 {len(all_results)} 组结果')
    # 组合选择器 — 影响各详情页显示哪组
    combo_labels = [f"Q={d['q_init']} L/h, 分区={d['n_zones']} (规格: {', '.join(str(s) for s in d['specs'])})" for d in all_results]
    sel_idx = st.selectbox('选择要查看的组合', range(len(all_results)),
        format_func=lambda i: combo_labels[i], key='batch_combo_selector')

    selected = all_results[sel_idx]
    tab_specs = []
    if selected['eq'] is not None:
        tab_specs.append(('eq', '等间距详情', 'batch_eq'))
    if selected['var'] is not None:
        tab_specs.append(('var', '变间距详情', 'batch_var'))
    both = selected['eq'] is not None and selected['var'] is not None
    if both:
        tab_specs.append(('both', '汇总对比', 'batch_both'))
    tab_specs.append(('grid', '批量总览', 'batch_grid'))
    tabs = st.tabs([label for _, label, _ in tab_specs])
    tab_by_kind = {kind: tab for (kind, _, _), tab in zip(tab_specs, tabs)}

    for kind, label, key_prefix in tab_specs:
        if kind not in ('eq', 'var'):
            continue
        r = selected[kind]
        with tab_by_kind[kind]:
            st.subheader(f"{label[:-2]} Q={selected['q_init']}, N={selected['n_zones']}, "
                         f"规格={', '.join(str(s) for s in selected['specs'])}")
            _render_result_tab(r, label[:-2], show_repl=bool(r['sweeps']))

    if both:
        with tab_by_kind['both']:
            r_eq, r_var = selected['eq'], selected['var']
            st.subheader(f"Cu 对比 — Q={selected['q_init']}, N={selected['n_zones']}")
            st.plotly_chart(make_cu_comparison(r_eq, r_var), width='stretch', key='batch_cu_comp')
            st.subheader('流量分布四线对比')
            st.plotly_chart(make_four_curve(r_eq, r_var), width='stretch', key='batch_four_curve')

    with tab_by_kind['grid']:
        st.subheader('热力图矩阵')
        metric = st.selectbox('指标', ['repl_cu', 'cu_delta', 'repl_q_range_span'],
            format_func=lambda x: {'repl_cu': '平均 Cu (%) — 越蓝越好',
                'cu_delta': 'Cu 增量 (%) — 越红提升越显著',
                'repl_q_range_span': '流量范围宽度 (L/h)'}.get(x, x),
            key='batch_heatmap_metric')
        grid_data = []
        for d in all_results:
            active = [r for r in (d['eq'], d['var']) if r is not None]
            grid_data.append({
                'q_init': d['q_init'], 'n_zones': d['n_zones'],
                'repl_cu': d['repl_cu'],
                'cu_delta': d['cu_delta'],
                'repl_q_range_span': sum(r['repl_q_range'][1]-r['repl_q_range'][0] for r in active)/len(active),
            })
        st.plotly_chart(make_heatmap(grid_data, metric), width='stretch', key='batch_heatmap')
        st.caption('热力图：横轴=分区数，纵轴=初始流量。点击上方选择器切换指标。')

        st.subheader('Cu 提升幅度排序')
        grid_sorted = sorted(grid_data, key=lambda d: d['cu_delta'], reverse=True)
        st.plotly_chart(make_cu_delta_bar(grid_sorted), width='stretch', key='batch_cu_delta_bar')

        st.subheader('结果汇总表')
        table_data = []
        for d in all_results:
            row = {
                'Q (L/h)': d['q_init'], '分区数': d['n_zones'],
                '规格列表': ', '.join(f'{s}' for s in d['specs']),
                '平均 Cu': f"{d['repl_cu']:.2f}%",
                '平均 ΔCu': f"{d['cu_delta']:+.2f}%",
            }
            if d['eq'] is not None:
                row['等间距 Cu'] = f"{d['eq']['repl_cu']:.2f}%"
            if d['var'] is not None:
                row['变间距 Cu'] = f"{d['var']['repl_cu']:.2f}%"
            table_data.append(row)
        st.dataframe(table_data, width='stretch')

        output = BytesIO()
        with pd.ExcelWriter(output, engine='openpyxl') as writer:
            # 汇总表
            pd.DataFrame(table_data).to_excel(writer, sheet_name='汇总', index=False)
            # 当前选中组合的详细17列数据 + 迭代记录
        sel = all_results[sel_idx]
        for label, r in [('等间距', sel['eq']), ('变间距', sel['var'])]:
            if r is None:
                continue
            build_full_export_df(r['positions'], r['init']['q'], r['k_init'],
                                 **r['config'], flow_exponent=r['init_n']).to_excel(
                                     writer, sheet_name=f'选中-{label}-初始', index=False)
            build_full_export_df(r['positions'], r['repl']['q'], r['repl']['k_arr'],
                                 **r['config'], flow_exponent=r['repl_n']).to_excel(
                                     writer, sheet_name=f'选中-{label}-替换', index=False)
            build_sweeps_df(r['sweeps']).to_excel(writer, sheet_name=f'选中-{label}-迭代记录', index=False)
        st.download_button('下载汇总 Excel', data=output.getvalue(),
                           file_name='批量结果汇总.xlsx',
                           mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet',
                           key='batch_download')


# ============ 主函数 ============
def main():
    st.title('滴头灌溉水力分析')
    st.caption('基于递减流量法 — 初始状态生成 → 分区迭代替换 → 均匀度分析')

    args = build_sidebar()

    # ==== 解析并校验参数（在计算前完成，避免界面直接异常） ====
    def _parse_floats(s):
        return [float(x.strip()) for x in s.replace('，', ',').split(',') if x.strip()]

    def _parse_ints(s):
        return [int(x.strip()) for x in s.replace('，', ',').split(',') if x.strip()]

    def _positive(values, label):
        if any(x <= 0 for x in values):
            raise ValueError(f'{label}必须全部大于 0')

    try:
        q_list = _parse_floats(args['q_init_rated'])
        n_list = _parse_ints(args['n_zones'])
        var_cycle = _parse_floats(args['var_cycle_str'])
        if not q_list:
            raise ValueError('初始额定流量不能为空')
        if not n_list:
            raise ValueError('分区数不能为空')
        _positive(q_list, '初始额定流量')
        if any(n < 1 for n in n_list):
            raise ValueError('分区数必须是大于等于 1 的整数')
        for key, label in [('h_spec_rated', '规格额定水头'), ('d_pipe', '管径'),
                           ('h_inlet', '入口水头'), ('g', '重力加速度'), ('nu', '运动粘度')]:
            if args[key] <= 0:
                raise ValueError(f'{label}必须大于 0')
        if args['max_iter'] < 1 or args['tolerance'] <= 0:
            raise ValueError('最大迭代次数必须大于 0，收敛容差必须大于 0')
        if not args['eq_enabled'] and not args['var_enabled']:
            raise ValueError('请至少启用等间距或变间距中的一种布局')
        if args['eq_enabled']:
            _positive([args['eq_spacing']], '等间距间距')
            if args['eq_trim_mode'] != 'target':
                _positive([args['eq_length']], '等间距铺设长度')
            if args['eq_trim_mode'] == 'cu_best':
                _positive([args['eq_scan_min'], args['eq_scan_max'], args['eq_scan_step']], 'Cu扫描参数')
                if args['eq_scan_min'] > args['eq_scan_max']:
                    raise ValueError('Cu扫描最小长度不能大于最大长度')
        if args['var_enabled']:
            if not var_cycle:
                raise ValueError('变间距循环间距不能为空')
            _positive([args['var_first_seg'], *var_cycle], '变间距间距')
            _positive([args['var_length']], '变间距生成长度')
            if args['var_trim_mode'] == 'target' and args['var_target_val'] == 'eq' and not args['eq_enabled']:
                raise ValueError('变间距匹配等间距时，请同时启用等间距计算')
        if args['do_repl']:
            if args['spec_mode'] == '手动指定' and not args['specs_str'].strip():
                raise ValueError('手动指定规格不能为空')
            specs_list = _parse_floats(args['specs_str']) if args['spec_mode'] == '手动指定' else auto_specs(q_list[0], n_list[len(n_list)//2])
            _positive(specs_list, '可用规格')
        else:
            specs_list = auto_specs(q_list[0], n_list[len(n_list)//2])
    except ValueError as e:
        st.error(f'参数错误：{e}')
        st.stop()

    # 按下"运行计算"时执行计算并存入session_state
    if args['run_btn']:
        st.session_state.pop('single_result', None)
        st.session_state.pop('batch_result', None)
        try:
            if args['mode'] == '单次':
                with st.spinner('正在计算…'):
                    ect = args['eq_trim_mode']
                    vct = args['var_trim_mode']
                    et = args.get('eq_target')
                    vt = args.get('var_target_val')
                    if (ect == 'target' and et is not None and et < 0.5) or \
                       (ect == 'cu_best') or \
                       (vct == 'target' and isinstance(vt, (int, float)) and vt < 0.5):
                        st.info('搜索模式可能需要较长时间，请耐心等待…')
                    run_single_compute(args, specs_list, var_cycle)
            else:
                st.info(f'运行 {len(q_list)}×{len(n_list)}={len(q_list)*len(n_list)} 组参数')
                run_batch_compute(args, specs_list, var_cycle, q_list, n_list)
        except Exception as e:
            st.error(f'计算出错: {e}')
            import traceback
            with st.expander('详细错误信息'):
                st.code(traceback.format_exc())

    # 从session_state读取并渲染结果（持久化，控件切换不丢失）
    if 'single_result' in st.session_state:
        render_single_tabs(current_do_repl=args['do_repl'])
    elif 'batch_result' in st.session_state:
        render_batch_tabs()
    else:
        st.info('请在左侧配置参数后点击"运行计算"')

    # 关闭translate div
    st.markdown('</div>', unsafe_allow_html=True)


if __name__ == '__main__':
    main()
