# Chapter 0 误差分析 Error Analysis

## 0.1 误差的来源

数值计算几乎都是近似计算，因此与理论值存在误差。数值分析的任务之一是将误差控制在一定的容许范围内或者至少对误差有所估计。

- 模型误差：问题的建立时候，因为忽略了一些次要的因素，而产生的误差。
- 观测误差：受实验观测手段所限，而产生的误差。
- **截断误差**：算法的近似解和模型的理论解之间存在的误差。
- **舍入误差**：由计算机四舍五入而产生的误差。

## 0.2 绝对误差，相对误差，有效数字

**绝对误差**：设 $x^*$ 为准确值，$x$ 是 $x^*$ 的一个近似值，则称 **$e=x^*-x$** 为近似值的绝对误差。若 $|e|\le\epsilon$，则称 $\epsilon$ 为绝对误差限。

**相对误差**：近似值的绝对误差与准确值的比值
$$
\frac{e}{x^*}=\frac{x^*-x}{x^*}\approx\frac{e}{x}
$$
称为相对误差, 记为 $e_r$。若 $|e_r|\le\epsilon_r$，则称 $\epsilon$ 为相对误差限。

**有效数字**：当准确值 $x^*$ 有多位时，常按四舍五入的原则得到 $x^*$ 的前几位近似值 $x$。

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 0.1 有效数字
</blockquote>

设
$$
\begin{gather*}
x=\pm\overline{a_1.a_2\dots a_n}\times10^m\\
a_1\in\{1,\dots,9\},a_i\in\{0,\dots,9\},(i\in\{2,\dots,n\}),m\in\mathbb Z
\end{gather*}
$$
且
$$
\color{red}|x-x^*|\le\frac{1}{2}\times10^{m-n+1}
$$
则称 $x$ 具有 $n$ 位有效数字。

<blockquote style="border-left: 5px solid #bb4545; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(188, 70, 70, 0.1)">
    例题
</blockquote>
> 用 $22/7$ 来近似 $\pi$，则有效数字和误差界为：
> $$\begin{align*}
> 22/7&=3.142857\dots\\
> \pi&=3.141593\dots\\
> |22/7-\pi|&\le0.0015\le0.5\times10^{-2}
> \end{align*}
> $$
> 因此有效数字为 $3$ 位，绝对误差限为 $0.5\times10^{-2}$。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 0.1
</blockquote>

若近似值 $x$ 具有 $n$ 位有效数字，则其相对误差限为：
$$
\color{red}\epsilon_r\le\frac{1}{2a_1}\times10^{-n+1}
$$
反之，若 $x$ 的相对误差限为：
$$
\epsilon_r\le\frac{1}{2(a_1+1)}\times10^{-n+1}
$$
则 $x$ 至少具有 $n$ 位有效数字。

> $$
> \begin{gather*}
> |x|=|\overline{a_1.a_2\dots a_n}\times10^m|\ge10^m\times a_1\\
> |e_r|\le|\frac{x-x^*}{10^m\times a_1}|\le\frac{1}{2a_1}\times10^{-n+1}
> \end{gather*}
> $$
>
> $$
> \begin{gather*}
> |x|\le10^m\times (a_1+1)\\
> |\epsilon|\le|10^m\times (a_1+1)\epsilon_r|\le\frac{1}{2}\times10^{m-n+1}
> \end{gather*}
> $$

近似数 $x$ 的有效位数越多，它的相对误差限越小；反之，$x$ 的相对误差越小，它的有效位数越多。

<blockquote style="border-left: 5px solid #bb4545; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(188, 70, 70, 0.1)">
    例题
</blockquote>

(1) 设某数的近似数具有 4 位有效数字，求其相对误差限。

(2) 某数的近似数的相对误差限为 0.3%，问其至少具有多少位有效数字？

(3) 计算 $\sin1.2$，需要取几位有效数字才能保证相对误差限不大于 0.01%？($\sin1.2=0.932039086$)

> (1)
> $$
> \epsilon_r\le\frac{1}{2a_1}\times10^{-n+1}\le\frac{1}{2}\times10^{-3}
> $$
> (2)
> $$
> \epsilon_r\le0.3\times10^{-2}<\frac{1}{2(9+1)}\times10^{-1}\\
> \Rightarrow n=2
> $$
> (3)
> $$
> \epsilon_r\le\frac{1}{2(9+1)}\times10^{-n+1}\le0.1\times10^{-3}\\
> \Rightarrow n=4
> $$

## 0.3 算数运算的误差

### 四则运算

设两个近似数 $x, y$, 其误差分别为 $e(x), e(y)$​, 则两数经四则运算得到的绝对误差：
$$
\begin{gather*}
e(x\pm y)=(x^*\pm y^*)-(x\pm y)=e(x)\pm e(y)\\
e(xy)=(x^*y^*)-(xy)=\cancel{e(x)e(y)}+xe(y)+ye(x)=xe(y)+ye(x)\\
e\left(\frac{x}{y}\right)=\frac{x^*}{y^*}-\frac{x}{y}=\frac{x^*y-xy^*}{y^*y}=\frac{ye(x)-xe(y)}{y^2}
\end{gather*}
$$
相对误差为：
$$
\begin{gather*}
e_r(x\pm y)=\frac{(x^*\pm y^*)-(x\pm y)}{(x^*\pm y^*)}=\frac{(x^*-x)\pm(y^*-y)}{(x^*\pm y^*)}=\frac{xe_r(x)\pm ye_r(y)}{x\pm y}\\
e_r(xy)=\frac{(x^*y^*)-(xy)}{x^*y^*}=\frac{xe(y)+ye(x)}{x^*y^*}=e_r(x)+e_r(y)\\
e\left(\frac{x}{y}\right)=\frac{x^*y-xy^*}{x^*y}=\frac{ye(x)-xe(y)}{xy}=e_r(x)-e_r(y)
\end{gather*}
$$

### 函数

如果函数 $A^*=f(x_1^*,\dots,x_n^*)$，其中 $f$ 为可微函数，则：
$$
\begin{align*}
e(A)=A^*-A&=f(x_1^*,\dots,x_n^*)-f(x_1,\dots,x_n)\\
&=\sum_{k=1}^n\left(\frac{\partial f(x_1,\dots,x_n)}{\partial x_k}\right)dx_k\\
&\approx\sum_{k=1}^n\left(\frac{\partial f(x_1,\dots,x_n)}{\partial x_k}\right)(x^*_k-x_k)\\
&=\sum_{k=1}^n\left(\frac{\partial f}{\partial x_k}\right)e(x_k)\\
\end{align*}
$$
绝对误差限为：
$$
\epsilon(A)=\sum_{k=1}^n|\frac{\partial f}{\partial x_k}|\epsilon(x_k)\\
$$
相对误差限为：
$$
\epsilon_r(A)=\frac{\epsilon(A)}{|A|}=\sum_{k=1}^n|\frac{\partial f}{\partial x_k}|\frac{\epsilon(x_k)}{|A|}\\
$$

## 0.4 适定性与稳定性

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 0.2 病态问题
</blockquote>

一个数值计算问题如果其输入数据的微小变动会引起输出数据（即计算结果）的很大扰动（误差）则称该数值问题属于病态问题。

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 0.3 条件数
</blockquote>

一个数值计算问题 $W$ 的条件数定义为问题的输出 $y$ 的相对误差 $e_r(y)$ 与输入 $x$ 的相对误差 $e_r(x)$ 之比的绝对值，记为：
$$
\mathrm{Cond}_W(x)=\frac{|e_r(y)|}{|e_r(x)|}
$$
当问题的条件数比较大时候，认为该问题属于病态的。

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 0.4 稳定性
</blockquote>
如果一个数值计算方法（算法）在计算过程中，将输入数据的初始误差放大得较大，则称该方法（算法）为数值不稳定的，反之，称为数值稳定的。

<blockquote style="border-left: 5px solid #bb4545; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(188, 70, 70, 0.1)">
    例题
</blockquote>

数列 $\{y_n\}$ 满足递推关系：
$$
y_n=10y_{n-1}-1
$$
若 $y_0=\sqrt{2}\approx1.41$，计算到 $y_{10}$ 时候误差有多大，这个计算稳定吗？

> $$
> \begin{gather*}
> y_n-a=10(y_{n-1}-a)\Rightarrow a=\frac{1}{9}\\
> y_{10}=10^{10}(y_0-\frac{1}{9})+\frac{1}{9}\\
> y_0=\sqrt{2}\Rightarrow y_{10}=12988888888.88889\\
> y_0=1.41\Rightarrow y_{10}=13031024512.61984\\
> \end{gather*}
> $$
>
> 相对误差和条件数分别为：
> $$\begin{gather*}
> e_r(x)=0.002979\\
> e_r(y)=0.003233\\
> \mathrm{Cond}_W(x)=1.085
> \end{gather*}
> $$
> 比较稳定吧。

## 0.5 避免和减少误差的若干计算原则

避免两相近数做减法。

避免除数过小。

防止大数"吃"小数。

简化计算步骤, 减少重复运算。

<blockquote style="border-left: 5px solid #bb4545; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(188, 70, 70, 0.1)">
    例题
</blockquote>

已知 $x=44.9222,y=44.9110$ 分别为 $x^*=\sqrt{2018},y^*=\sqrt{2017}$ 的近似值。估计近似数 $x-y$ 的有效数字位数。

> $$
> \begin{gather*}
> |e(x)|\le\frac{1}{2}\times10^{-4}\quad |e(y)|\le\frac{1}{2}\times10^{-4}\\
> |e(x-y)|\le|e(x)|+|e(y)|=10^{-4}\le\frac{1}{2}\times10^{-3}
> \end{gather*}
> $$
>
> 因为 $m=-1$，则 $n=3$，即有效数字仅有 $3$ 位。

## 0.6 向量、矩阵和连续函数的范数

### 0.6.1 向量范数 

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 0.5 向量范数
</blockquote>

向量 $\mathbf x$ 的实值函数 $N(\mathbf x)$ 满足：
$$
\begin{align*}
\text{Nonnegativity: }&N(\mathbf x)\ge0,N(\mathbf x)=0\Leftrightarrow\mathbf x=\mathbf 0\\
\text{Homogeneity: }&N(k\mathbf x)=|k|N(\mathbf x)\\
\text{Triangle inequation: }&N(\mathbf x+\mathbf y)\le N(\mathbf x)+N(\mathbf y)
\end{align*}
$$

几个常用的向量范数：
$$
\begin{align*}
\|x\|_1&=\sum_{i=1}^n|x_i|&\|x\|_2&=(\sum_{i=1}^n|x_i|^2)^{1/2}\\
\|x\|_p&=(\sum_{i=1}^n|x_i|^p)^{1/p}&\|x\|_{\infty}&=\max_{1\le i\le n}|x_i|\\
\end{align*}
$$

*向量范数是对向量长度的一种衡量。*

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 0.2
</blockquote>


设 $N(\mathbf x)$ 是 $\mathbb R^n$ 上的任意向量范数，则 $N(\mathbf x)$ 是关于其分量 $x_1,\dots,x_n$ 的连续函数。

> $$
> \lim_{\Delta\mathbf x\to\mathbf 0}N(\mathbf x+\Delta\mathbf x)-N(\mathbf x)\le\lim_{\Delta\mathbf x\to\mathbf 0}N(\Delta\mathbf x)=0
> $$

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 0.6 向量范数等价
</blockquote>

设 $\|\mathbf x\|_s,\|\mathbf x\|_t$ 是 $\mathbb R^n$ 上的任意两种向量范数，则存在常数 $c_1,c_2>0$，使得 $c_1\|\mathbf x\|_s\le\|\mathbf x\|_t\le c_2\|\mathbf x\|_s$，则两种范数等价

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 0.7 向量收敛
</blockquote>

设 $\{\mathbf x^{(k)}=(x_1^{(k)},\dots,x_n^{(k)})^T\}, \mathbf x^{*}=(x_1^{*},\dots,x_n^{*})^T$。如果：
$$
\lim_{n\to\infty}x_k^{(n)}=x_k^*
$$
则称 $\{\mathbf x^{(k)}\}$ 收敛于 $\mathbf x^*$，记为 $\lim_{n\to\infty}\mathbf x^{(n)}=\mathbf x^*$。

<blockquote style="border-left: 5px solid #bb4545; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(188, 70, 70, 0.1)">
    例题
</blockquote>

设 $P\in\mathbb R^{n\times n}$，$N(\mathbf x)$ 是 $\mathbb R^n$ 上的范数，定义范数 $M(\mathbf x)=N(P\mathbf x)$，证明 $M(\mathbf x)$ 也是 $\mathbb R^n$ 上的范数。

> $$
> \begin{align*}
> \text{Nonnegativity: }&M(\mathbf x)=N(P\mathbf x)\ge0,M(\mathbf x)=N(P\mathbf x)=0\Leftrightarrow P\mathbf x=\mathbf x=\mathbf 0\\
> \text{homogeneity: }&M(k\mathbf x)=N(kP\mathbf x)=N(Pk\mathbf x)=|k|M(\mathbf x)\\
> \text{Triangle inequation: }&M(\mathbf x+\mathbf y)=N(P\mathbf x+P\mathbf y)\le N(P\mathbf x)+N(P\mathbf y)=M(\mathbf x)+M(\mathbf y)
> \end{align*}
> $$
>
> 

### 0.6.2 矩阵范数

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 0.8 矩阵范数
</blockquote>

矩阵 $A\in\mathbb R^{n\times n}$ 的实值函数 $N(A)$ 满足：

$$
\begin{align*}
\text{Nonnegativity: }&N(A)\ge0,N(A)=0\Leftrightarrow A=\mathbf 0\\
\text{Homogeneity: }&N(kA)=|k|N(A)\\
\text{Triangle inequation: }&N(A+B)\le N(A)+N(B)\\
\text{Multiplication: }&N(AB)\le N(A)N(B)
\end{align*}
$$

矩阵对向量变换的能力体现在向量变换前后其长度的缩放程度，矩阵范数越大，越可能把一个向量拉得更长。

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 0.9 相容
</blockquote>

如果向量范数 $\|\mathbf x\|_{\alpha}$ 与矩阵范数 $\|A\|_{\beta}$ 满足：
$$
\|A\mathbf x\|_{\alpha}\le\|A\|_{\beta}\|\mathbf x\|_{\alpha}
$$
则这两个相容。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 0.3 算子范数
</blockquote>

定义：
$$
\|A\|_{\alpha}=\max_{\mathbf x\ne\mathbf 0}\frac{\|A\mathbf x\|_{\alpha}}{\|\mathbf x\|_{\alpha}}=\max_{\|\mathbf x\|_{\alpha}=1}\|A\mathbf x\|_{\alpha}
$$
则 $\|A\|_{\alpha}$ 为矩阵范数，且满足相容条件：
$$
\|A\mathbf x\|_{\alpha}\le\|A\|_{\alpha}\|\mathbf x\|_{\alpha}
$$

计算 $\|A^{-1}\|$ 的范数：
$$
\|A^{-1}\|_{\alpha}=\max_{\mathbf y\ne\mathbf 0}\frac{\|A^{-1}\mathbf y\|_{\alpha}}{\|\mathbf y\|_{\alpha}}=1/\left(\min_{\mathbf y\ne\mathbf 0}\frac{\|\mathbf y\|_{\alpha}}{\|A^{-1}\mathbf y\|_{\alpha}}\right)=1/\left(\min_{\mathbf x\ne\mathbf 0}\frac{\|A\mathbf x\|_{\alpha}}{\|\mathbf x\|_{\alpha}}\right)
$$
$A$ 对向量的拉伸能力，等同于 $A^{-1}$ 对向量的压缩能力。

可得不等式：
$$
\|A\mathbf x\|_{\alpha}\le\|A\|_{\alpha}\|\mathbf x\|_{\alpha}\\
\|\mathbf x\|_{\alpha}\le\|A^{-1}\|_{\alpha}\|A\mathbf x\|_{\alpha}
$$

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 0.10 Frobenius范数
</blockquote>

令：
$$
\|A\|_F=\left(\sum_{i=1}^{n}\sum_{j=1}^{n}|a_{ij}|^2\right)^{1/2}
$$
则 $\|A\|_F$ 是矩阵的 Frobenius 范数。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 0.4
</blockquote>

矩阵的行范数：
$$
\|A\|_{\infty}=\max_{1\le i\le n}\sum_{j=1}^{n}|a_{ij}|
$$
矩阵的列范数：
$$
\|A\|_{1}=\max_{1\le j\le n}\sum_{i=1}^{n}|a_{ij}|
$$
矩阵的 2 范数：
$$
\|A\|_{2}=\sqrt{\lambda_{\max}(A^TA)}
$$

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 0.11 矩阵收敛
</blockquote>

在定义了范数的空间 $\mathbf R^{n\times n}$ 中，如果存在 $A\in\mathbf R^{n\times n}$，满足
$$
\lim_{k\to\infty}\|A^{(k)}-A\|=0
$$
则称矩阵列 $\{A^{(k)}\}$ 收敛于矩阵 $A$。

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 0.12 谱半径
</blockquote>

设 $\lambda_i$ 为 $n$ 阶矩阵 $A$ 的特征值，称
$$
\rho(A)=\max_{i}|\lambda_i|
$$
为矩阵A的谱半径。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 0.5 谱半径定理
</blockquote>

设 $A\in\mathbb{R}^{n\times n}$，则 $\rho(A)\le\|A\|$。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 0.6 可逆推论
</blockquote>

设矩阵 $B$ 满足 $\|B\|\le1$，则 $I\pm B$ 为非奇异矩阵，且：
$$
\frac{1}{1+\|B\|}\le\|(I+B)^{-1}\|\le\frac{1}{1-\|B\|}
$$

