# Chapter 3 插值 Interpolation

## 3.1 数据和插值函数

### 3.1.0 插值定义

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 3.1
</blockquote>

对于函数 $y=P(x)$ 和一组数据点 $\{(x_i,y_i)\mid x_i\ne x_j\}$，如果对于每一个 $i$，都有 $P(x_i)=y_i$，则称该函数**插值**数据点 $\{(x_i,y_i)\}$。这样的函数总是存在的。

设已知区间 $[a,b]$ 上的实值函数 $f(x)$ 在 $n+1$ 个点 $x_0,x_1,\dots,x_n$ 上的函数值 $f(x_0),\dots,f(x_n)$。在函数类
$$
\Phi_n=\mathrm{span}(\phi_0,\dots,\phi_n)
$$
中找一个函数 $\phi(x)$ 满足：
$$
\phi(x_i)=f(x_i)
$$
称其为插值条件。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 3.0
</blockquote>

设 $f(x),x\in[a, b],\{x_i\}$ 是 $f(x)$ 在 $[a,b]$ 上 $n+1$ 个互异的节点，$y_i=f(x_i)$，则满足插值条件
$$
P(x_i)=y_i
$$
的插值多项式 $P_n(x)\in\Phi_n$ 存在且唯一。

### 3.1.1 拉格朗日插值 Lagrange's Interpolation

设多项式空间 $\Phi_n=\mathrm{span}(l_0,\dots,l_n)$，其中 $l_k(x)$ 为 $n$ 次多项式，令 $L_n(x)\in\Phi(x)$，则
$$
L_n(x)=y_0l_0(x)+\dots+y_nl_n(x)
$$

如果 $l_i(x_j)=\delta_{ij}$ 则满足插值多项式的条件。因此设
$$
l_i(x)=a_k(x-x_0)\dots(x-x_n)=a_i\prod_{j\ne i}^{n}(x-x_j)
$$
利用 $l_i(x_i)=1$ 可得
$$
a_i=\frac{1}{\prod_{j\ne i}^{n}(x_i-x_j)}
$$
因此 $n+1$ 组点可以定义 $n$ 阶拉格朗日插值多项式：
$$
\color{red}L_n(x)=\sum_{i=0}^{n}y_i\left(\prod_{j=0,j\ne i}^{n}\frac{x-x_j}{x_i-x_j}\right)
$$

记 $\omega_{n+1}(x)=\prod_{j=0}^n(x-x_j)$，则：
$$
\omega_{n+1}'(x_i)=\prod_{j=0,j\ne i}^{n}(x_i-x_j)\quad l_i(x)=\frac{\omega_{n+1}(x)}{(x-x_i)\omega_{n+1}'(x_i)}
$$
则：
$$
\color{red}L_n(x)=\sum_{i=0}^{n}y_i\frac{\omega_{n+1}(x)}{(x-x_i)\omega_{n+1}'(x_i)}
$$

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 3.1
</blockquote>

对于一组数据点 $\{(x_i,y_i)\mid x_i\ne x_j\}$，有且仅有一个 $n$ 次或者更低次的多项式 $P$ 满足 $P(x_i)=y_i$。


> 存在性通过上述构造得到。下面证明唯一性。
>
> 假设存在两个这样的多项式，设为 $P(x)$ 和 $Q(x)$，他们都为至多 $n$ 次并插值过 $n+1$ 个点。定义 $H(x)=P(x)-Q(x)$，则在插值点上 $H(x_i)=0$，这样的零点有 $n+1$ 个。由代数学基本定律可得，$H(x)$ 必为 $n+1$ 次多项式或者零多项式。因为 $n$ 次多项式的和不能构造出 $n+1$ 次多项式，所以 $H(x)\equiv0$，所以 $P(x)=Q(x)$。

<blockquote style="border-left: 5px solid #b94263; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(185, 66, 110, 0.1)">
    例题
</blockquote>

$f(x)\in C^2[a,b]$ 证明
$$
\max_{a\le x\le b}\left|f(x)-f(a)-\frac{f(b)-f(a)}{b-a}(x-a)\right|\le\frac{1}{8}(b-a)^2\max_{a\le x\le b}|f''(x)|
$$

> $$
> \max_{a\le x\le b}|R_1(x)|=\max_{a\le x\le b}\left|\frac{f^{(2)}(\xi)}{2!}(x-a)(x-b)\right|\le\max_{a\le x\le b}\left|f^{(2)}(x)\right|\frac{(b-a)^2}{4\times2!}
> $$

### 3.1.2 牛顿差商 Newton's divided difference

当增减一个插值点，拉格朗日插值需要变动所有的基函数。

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 3.3
</blockquote>

记 $N_i(x)$ 为函数 $f$ 在节点 $x_0,\dots,x_i$ 上的插值多项式，由于：
$$
N_n(x_i)=N_{n-1}(x_i)=y_i,i\in\{1,\dots,n-1\}
$$
从而：
$$
N_n(x)=N_{n-1}(x)+c\omega_n(x)
$$

利用插值条件 $N_n(x_n)=y_n$ 和 Lagrange 插值 $N_{n-1}(x)=\sum_{k=0}^{n-1}y_kl_k(x)$ 可得：
$$
\begin{align*}
c&=\frac{y_n}{\omega_n(x_n)}-\frac{N_{n-1}(x_n)}{\omega_n(x_n)}\\
&=\frac{y_n}{\prod_{i=0}^{n-1}(x_n-x_i)}-\frac{\sum_{k=0}^{n-1}y_kl_k(x_n)}{\prod_{i=0}^{n-1}(x_n-x_i)}\\
&=\frac{y_n}{\prod_{i=0}^{n-1}(x_n-x_i)}+\frac{\sum_{k=0}^{n-1}y_k}{\prod_{i=0,i\ne k}^{n-1}(x_i-x_k)}\\
&=\sum_{k=0}^{n}\frac{y_k}{\prod_{i=0,i\ne k}^{n}(x_k-x_i)}\\
&=\sum_{k=0}^{n}\frac{y_k}{\omega'_{n+1}(x_k)}=\sum_{k=0}^{n}\frac{f(x_k)}{\omega'_{n+1}(x_k)}\\
\end{align*}
$$
用 $f([x_0\dots x_n])$ 表示唯一多项式的 $x^{n}$ 项的系数，这些点插值 $\{(x_i,y_i)\}$。

用上述定义来定义以下的牛顿差商公式：
$$
N_n(x)=f[x_0]+f[x_0x_1](x-x_0)+\dots+f[x_0\dots x_n](x-x_0)\cdots(x-x_{n-1})
$$
如何计算差商（稍后证明）：
$$
f[x_k]=f(x_k)=y_k\\
f[x_kx_{k+1}]=\frac{f[x_{k+1}]-f[x_k]}{x_{k+1}-x_k}\\
\dots\\
f[x_k\dots x_{k+m}]=\frac{f[x_{k+1}\dots x_{k+m}]-f[x_k\dots x_{k+m-1}]}{x_{k+m}-x_k}
$$

以 4 个点为例：
$$
\begin{matrix}
x_0\vline&f[x_0]& &\\
\quad\vline& &f[x_0x_1]&\\
x_1\vline&f[x_1]& &f[x_0x_1x_2]\\
\quad\vline& &f[x_1x_2]& &f[x_0x_1x_2x_3]\\
x_2\vline&f[x_2]& &f[x_1x_2x_3]\\
\quad\vline& &f[x_2x_3]&\\
x_3\vline&f[x_3]& \\
\end{matrix}
$$
<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 3.2
</blockquote>

$k$ 阶差商 $f[x_0\dots x_k]$ 是函数值 $f(x_0),\dots,f(x_k)$ 的线性组合。
$$
f[x_0\dots x_k]=\sum_{j=0}^{k}\frac{f(x_j)}{\omega_{n+1}'(x_j)}
$$
> 归纳法。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 3.3
</blockquote>

由插值多项式的唯一性可得，差商对于定义它的节点而言是对称的，即任意改变节点次序，$f[x_0\dots x_k]$ 的值不变。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 3.4
</blockquote>

如果已知 $f(x)\in C^n[a,b]$，$x_i\in[a,b]$，则：
$$
f[x_0\dots x_n]=\frac{1}{n!}f^{(n)}(\xi)\quad \xi\in(a,b)
$$
> 由插值多项式的唯一性可得：
> $$
> f[x_0\dots x_nx]\omega_{n+1}(x)=\frac{f^{(n+1)}(\xi)}{(n+1)!}\omega_{n+1}(x)
> $$
> 因为 $x$ 是任意的：
> $$
> \begin{align*}
> f[x_0\dots x_nx_{n+1}]&=\frac{f^{(n+1)}(\xi)}{(n+1)!}\\
> f[x_0\dots x_n]&=\frac{f^{n}(\xi)}{n!}\\
> \end{align*}
> $$

> 对于某个固定的，非 $x_i$ 的 $x$，设
> $$
> \phi(t)=f(t)-N_n(t)-f[x_0\dots x_nx](t-x_0)\cdots(t-x_n)
> $$
> 则 $\phi(t)$ 具有 $n+2$ 个零点。对 $t$ 求 $n$ 阶导数后，$\phi^{(n)}(t)$ 具有 $2$ 个零点，此时：
> $$
> \phi^{(n+1)}(\xi)=f^{(n+1)}(\xi)-0-f[x_0\dots x_nx](n+1)!=0
> $$
> 因为 $x$ 是任意的：
> $$
> \begin{align*}
> f[x_0\dots x_nx_{n+1}]&=\frac{f^{(n+1)}(\xi)}{(n+1)!}\\
> f[x_0\dots x_n]&=\frac{f^{n}(\xi)}{n!}\\
> \end{align*}
> $$

<blockquote style="border-left: 5px solid #bb4545; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(188, 70, 70, 0.1)">
    例题
</blockquote>

已知 $f(x)=-6x^9+8x^7+5x-4$，计算 $f[2^02^1\dots2^9]$ 和 $f[2^02^1\dots2^{10}]$

> $$\begin{gather*}
> f[2^02^1\dots2^9]=\frac{1}{9!}f^{(n)}(\xi)=-6\\
> f[2^02^1\dots2^{10}]=0
> \end{gather*}
> $$

### 3.1.3 经过 n 个点的 d 阶多项式有几个

n 阶和 n 阶以上的每一阶都有无数个。

n - 1 阶和以下的可能为 0 个，1 个和无数个。

### 3.1.4 Hermite 插值

既要求插值函数在节点处与已知的函数值相等, 同时要求插值函数的某些导数与已知函数的导数值相等。

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 3.4 两点三次 Hermite 插值
</blockquote>

已知 $f(x_0)=y_0,f(x_1)=y_1,f'(x_0)=m_0,f'(x_1)=m_1$。设
$$
H_3(x)=\alpha_0(x)y_0+\alpha_1(x)y_1+\beta_0(x)m_0+\beta_1(x)m_1
$$
其中：
$$\begin{gather*}
\alpha_i(x_j)=\delta_{ij}\quad\alpha'_i(x_j)=0\\
\beta_i(x_j)=0\quad\beta'_i(x_j)=\delta_{ij}
\end{gather*}
$$
利用 Lagrange 插值基函数 $l_i(x_j)=\delta_{ij}$，则设
$$
\alpha_i(x)=[a(x-x_i)+1]l_i^2(x)
$$
则：
$$
\alpha'_i(x)=al_i^2(x)+[a(x-x_i)+1]2l_i(x)l'_i(x)\Rightarrow a=-2l'_i(x)
$$
$\alpha_i(x)=[1-2l'_i(x)(x-x_i)]l_i^2(x)$。同理得到 $\beta_i(x)=(x-x_i)l_i^2(x)$。

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 3.5 一般 Hermite 插值
</blockquote>

$n+1$ 个点得到 $2n+1$ 阶多项式。

已知在节点 $a\le x_0<x_1<\cdots<x_n\le b$ 上，$y_j=f(x_j),m_j=f'(x_j)$，求插值多项式满足 Hermite 条件。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 3.4
</blockquote>

满足上述插值条件的 Hermite 多项式 $H_{2n+1}(z)$ 存在且唯一，且
$$\begin{gather*}
H_{2n+1}=\sum_{j=0}^n\alpha_j(x)y_j+\sum_{j=0}^n\beta_j(x)m_j\\
\alpha_j(x)=[1-2(x-x_j)l'_j(x_j)]l^2_j(x)\quad \beta_j(x)=(x-x_j)l^2_j(x)\end{gather*}
$$
这里 $\alpha,\beta$ 同上。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 3.5
</blockquote>

设 $f(x)$ 在 $[a,b]$ 上具有 $2n+2$ 阶导数，则：
$$
R_{2n+1}=f(x)-H_{2n+1}(x)=\frac{f^{(2n+2)}(\xi)}{(2n+2)!}\omega^2_{n+1}(x)
$$
### 3.1.5 Newton-Hermite 插值

设 $x_1<x_2<\dots<x_s$，$y^{(h)}_k$ 是给定的实数。其中 $h\in\{0,\dots,a_{k}-1\},k\in\{1, \dots,s\},\{a_s\}\in\mathbb Z^+,\sum a_s=n+1$。

求一个次数不超过 $n$ 次的多项式 $P_n(x)$ 满足
$$
P^{(h)}_n(x_k)=y^{(h)}_k\quad h\in\{0,\dots,a_{k}-1\},k\in\{1, \dots,s\}
$$
1. 将求节点按重复出现次数重排，即求节点为
$$
    \left\{\underbrace{x_1,\dots,x_1}_{a_1},\dots,\underbrace{x_s,\dots,x_s}_{a_s}\right\}
$$
    再将上述点重新按顺序编号，得到 $z_0,\dots,z_n$，称为具有重节点的插值节点组。

2. 根据插值节点组及插值条件，构造差商表。这里，需要利用重节点差商公式
    $$
    f[\underbrace{x_kx_k\dots x_k}_h]=\frac{f^{(h-1)}(x_k)}{(h-1)!}
    $$

3. 写出具有重节点的牛顿插值多项式
    $$
    P_n(x)=f(z_0)+f[z_0z_1](x-z_0)+\cdots+f[z_0,\dots,z_n](x-z_0)\cdots(x-z_{n-1})
    $$

<blockquote style="border-left: 5px solid #bb4545; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(188, 70, 70, 0.1)">
    例题
</blockquote>

求满足 $P(x_j)=f(x_j)(j=0,1,2)$ 和 $P'(x_1)=f'(x_1)$ 的插值多项式 $P(x)$ 及其余项表达式。

> 重排：
> $$
> x_0,x_1,x_1,x_2
> $$
> 系数：
> $$
> \begin{matrix}
> x_0\vline&f[x_0]& &\\
> \quad\vline& &f[x_0x_1]&\\
> x_1\vline&f[x_1]& &f[x_0x_1x_1]\\
> \quad\vline& &f[x_1x_1]& &f[x_0x_1x_1x_2]\\
> x_1\vline&f[x_1]& &f[x_1x_1x_2]\\
> \quad\vline& &f[x_1x_2]&\\
> x_2\vline&f[x_2]& \\
> \end{matrix}
> $$
> 多项式：
> $$
> P_3(x)=f(x_0)+f[x_0x_1](x-x_0)+f[x_0x_1x_1](x-x_0)(x-x_1)\\+f[x_0x_1x_1x_2](x-x_0)(x-x_1)(x-x_1)
> $$
> 余项：
> $$
> R_3(x)=f(x)-P_3(x)=\frac{f^{(4)}(\xi)}{4!}(x-x_0)(x-x_1)^2(x-x_2)
> $$
> 

## 3.2 插值误差

### 3.2.1 Lagrange 多项式的余项

称 $R_n(x)=f(x)-P_n(x)$ 为 Lagrange 多项式的余项。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 3.6
</blockquote>


如果 $f(x)\in C^n[a,b]$，$f^{(n+1)}$ 在 $(a,b)$ 上存在，则对于 $\forall x\in[a,b],\exists\xi\in(a,b)$，使得：
$$
R_n(x)=\frac{f^{(n+1)}(\xi)}{(n+1)!}\omega_{n+1}(x)
$$

> $R_n$ 一定满足
> $$
> R_n(x)=A(x)(x-x_0)\dots(x-x_n)
> $$
> 其中 $A(x)$ 是个系数函数。
>
> 对于某个固定的，非 $x_i$ 的 $x$，设
> $$
> \phi(t)=f(t)-P_n(t)-A(x)(t-x_0)\dots(t-x_n)
> $$
> 则
> $$\begin{gather*}
> \phi(x)=f(x)-P_n(x)-R_n(x)=0\\
> \phi(x_i)=f(x_i)-P_n(x_i)-0=0
> \end{gather*}
> $$
> 因此 $\phi(t)$ 具有 $n+2$ 个零点。对 $t$ 求 $n$ 阶导数后，$\phi^{(n)}(t)$ 具有 $2$ 个零点，此时：
> $$
> \phi^{(n+1)}(t)=f^{(n+1)}(t)-0-A(x)(n+1)!=0\\
> \Longrightarrow A(x)=\frac{f^{(n+1)}(t)}{(n+1)!}
> $$
> 因此 $R_n$ 满足该形式。

如果 $f(x)\in\Phi_n$，则 $L_n\equiv f(x)$，$\sum_{k=0}^{n}l_k(x)\equiv1$ (取 $y_i=1$ 即可得)。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 3.7
</blockquote>



记 $h=\max_{1\le i\le n}(x_i-x_{i-1})$，则：
$$
\|f-L_n\|_{\infty}\le\frac{h^{n+1}}{4(n+1)}\|f^{(n+1)}\|_{\infty}
$$

> 已知：
> $$
> f-L_n=\frac{f^{(n+1)}(\xi)}{(n+1)!}\omega_{n+1}(x)
> $$
> 取最大段的中点：
> $$
> \omega_{n+1}(x)\le\frac{1}{4}h^2n!h^{n-1}
> $$
> 得证。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 3.8
</blockquote>



牛顿插值公式余项为：
$$
R_n(x)=f(x)-N_n(x)=f[x_0\dots x_nx]\omega_{n+1}(x)
$$

> $$
> \begin{align*}
> f(x)&=f(x_0)+f[x_0x_1](x-x_0)+\cdots\\
> &+f[x_0x_1\dots x_{n}](x-x_0)\cdots(x-x_{n-1})\\
> &+f[x_0x_1\dots x_{n}x](x-x_0)\cdots(x-x_{n})\\
> &=N_n(x)+R_n(x)
> \end{align*}
> $$
> 
> 对于任意 $x$ 都成立。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 3.9
</blockquote>


若 $f\in C^4[a,b]$ 则利用 Rolle 定理, 可得 Hermite 插值余项为：
$$
R_3(x)=f(x)-H_3(x)=\frac{1}{4!}f^{(4)}(\xi)(x-x_0)^2(x-x_1)^2
$$

> $$\begin{gather*}
> R_3(x_i)=f(x_i)-H_3(x_i)=0\\
> R'_3(x_i)=f'(x_i)-H'_3(x_i)=0
> \end{gather*}
> $$
>则 $R_3(x)$ 具有如下形式：
> $$
> R_3(x)=A(x)(x-x_0)^2(x-x_1)^2
> $$
> 对于某个固定的，非 $x_i$ 的 $x$，设
> $$
> \phi(t)=f(t)-H_3(t)-A(x)(t-x_0)^2(t-x_1)^2
> $$
> 则其有 3 个零点。
> 
>而 $\phi'(t)$ 具有 4 个零点，则 $\phi^{(4)}(\xi)=0,\xi\in(x,x_0,x_1)$。
> $$
> \phi^{(4)}(\xi)=f^{(4)}(\xi)-4!A(x)=0\\
> \Longrightarrow A(x)=\frac{f^{(4)}(\xi)}{4!}
> $$

