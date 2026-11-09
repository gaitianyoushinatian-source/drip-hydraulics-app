# 怎么把项目传到 GitHub（第一次用，照着点就行）

不用记命令，全部用鼠标点。

---

## 先看这一条：打开 GitHub Desktop 是英文的「Let's get started!」页面

右侧有 4 个按钮，**只点最下面那个**：

- ✅ **Add an Existing Repository from your local drive...**（第 4 个，文件夹图标）—— 点这个
- ❌ Create a tutorial repository...（蓝色那个）—— 教程，别点
- ❌ Clone a repository from the Internet... —— 下载别人的项目，别点
- ❌ Create a New Repository on your local drive... —— 会另外新建一个子文件夹，别点

左边「Your repositories」列出的是已有仓库，也**都不要点**。然后照下面第 2 步继续。

### 接着会弹出「Create a new repository」对话框

只要最下面写着 **The repository will be created at: D:\Project\滴灌优化**，就说明位置对了，**直接点蓝色 Create repository**，别的都不用动。这一步只是给这个文件夹加一个版本记录本，不会删文件、不会新建别的文件夹。

- `Initialize this repository with a README` 保持**不勾**（你已经有一份 README.md）
- `Git ignore` / `License` 保持 `None`（你已经有一份 .gitignore）
- ⚠️ **这里不要把 Name 改成英文**：这个框是靠「Local path + Name」拼出文件夹位置的，一改就会新建 `D:\Project\英文名` 的空文件夹，你的 5 个文件反而没被登记。中文名只影响本机，等后面点 Publish repository 时再换英文名，两边互不影响

---

## 事先说明

你的项目已经在这个文件夹里了：

```
D:\Project\滴灌优化
```

里面有 5 个文件（app.py、README.md 等）。**不用再复制、不用改名、不用新建文件夹**，直接用这个文件夹。

---

## 第 1 步：准备两个账号/软件

1. 去 https://github.com 注册一个账号（网名、邮箱、密码），记住账号密码。
2. 去 https://desktop.github.com 下载 **GitHub Desktop**，一路点"下一步"装好。
3. 打开 GitHub Desktop，它会让你登录，点 **Sign in to GitHub.com**，浏览器会跳出来让你确认，按提示点绿色的确认按钮就行。

> 如果英文看不懂，可以把界面截给我，我告诉你点哪个。

---

## 第 2 步：把你的文件夹交给 GitHub Desktop

1. 打开 GitHub Desktop，看左上角菜单 **File（文件）**，点它 → 选 **Add local repository…（添加本地仓库）**。
2. 弹出的窗口里点 **Choose…（选择）**，找到并选中这个文件夹：

   ```
   D:\Project\滴灌优化
   ```

   选好后点"选择文件夹"。
3. 它大概率会弹出一句红字，说"这个文件夹还不是 Git 仓库"。**别慌，这是正常的**，旁边有个蓝色小字 **create a repository**，点它。
4. 接着出现一个填名字的框：
   - Name（名字）**只能用英文字母、数字、点、横线，不能用中文**，比如 `drip-hydraulics-2026`
   - ⚠️ 如果填了中文（例如「滴灌优化-2026」），GitHub 会把中文吃掉，只剩 `--2026`，你换几个中文名其实都变成同一个，就会一直提示"重复"。看到"你的新仓库将创建为 --2026"就是这个问题
   - 如果提示"名称重复 / 已存在"，换个英文名再点一次（后面加个数字就行）。这个名字只是网上的网址，跟电脑上的文件夹、跟程序都没关系，不用改文件夹名
   - 其他都不用动
   - 点 **Create repository（创建仓库）**

做完这步，你的文件夹就被"登记"了，但**还没传到网上**。

---

## 第 3 步：存一次（相当于拍个快照）

1. 现在看软件左边，会列出 5 个文件，每个前面打着勾。勾表示"这些要保存"，**别取消**。
2. 看左下角有个输入框（写着 Summary），在里面打一行字，比如：

   ```
   第一次上传：滴灌水力计算程序
   ```

3. 打完点下面那个蓝色按钮 **Commit to main**。

这一步只是在你自己电脑上存了一次记录，**还没上网**，放心点。

---

## 第 4 步：传到 GitHub 网站

1. 点软件上方中间或右上角的 **Publish repository（发布仓库）**。
2. 会弹出一个小框：
   - 名字保持刚才那个就行
   - **想只给自己看**：把 **Keep this code private** 那个勾打上
   - **想让别人也能看到**：把那个勾去掉
   - 点 **Publish repository**
3. 等几秒，按钮会变成 **Push origin** 之类，点一下最上面的 **View on GitHub**，浏览器就会打开你的项目页面。

**到这一步就成功了。** 以后打开 https://github.com 登录，就能看到它。

---

## 以后改了代码，怎么留记录

每次改完一个地方，做两步：

1. 回到 GitHub Desktop，左边会列出你改过的文件，左下角 Summary 写一句人话，比如：

   ```
   把末端水头目标改成 9.99
   ```

   点 **Commit to main**。
2. 再点上面的 **Push origin**（或 Fetch/Push 那个按钮），就同步到网上了。

一句话记：**Commit = 存到本机；Push = 传到网上。**

---

## 卡住了怎么办

| 看到的情况 | 怎么办 |
|---|---|
| 说"不是 Git 仓库" | 正常，点 **create a repository** |
| 左边文件列表是空的 | 说明选错文件夹了，重新 File → Add local repository，选 `D:\Project\滴灌优化` |
| 登录一直转圈 | 看浏览器有没有弹出授权页面没点确认 |
| 提示要填名字/邮箱 | 随便填个昵称和注册时的邮箱就行 |
| 已经在网页上建好了一个空仓库 | 那个页面**什么都不用点**，直接回 GitHub Desktop 做 Commit → Publish；Publish 时换一个没被占用的英文名，网页上那个空仓库留着不影响 |
| 提示"仓库名称已存在 / 重复" | 换一个**纯英文**名字，比如 `drip-hydraulics-app`。中文名会被 GitHub 吃掉，多个中文名最后都变成同一个，才显得一直重复 |
| 完全不知道点哪 | **直接截屏发我**，我按图告诉你点哪个 |

---

## 一句话总结

装 GitHub Desktop → 登录 → 把 `D:\Project\滴灌优化` 加进去 → 写句话点 Commit → 点 Publish = 完成。

不用怕点错，这些操作都能重来，也不会删你电脑上的文件。
