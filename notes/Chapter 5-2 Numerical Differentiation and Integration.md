# Chapter 5 数值微分和积分 Numerical Differentiation and Integration

## 5.2 Newton-Cotes 公式

### 5.2.1 Newton-Cotes 公式

由积分中值定理
$$
\int_b^af(x)\mathrm{d}x=f(\xi)(b-a),\xi\in[a,b]
$$
利用 $f(x)$ 在节点 $a=x_0<x_1<\cdots<x_n=b$ 处函数值的加权平均来近似代替 $f(\xi)$。

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 机械求积公式
</blockquote>

$$
\int_a^bf(x)\mathrm{d}x\approx\sum_{i=0}^{n}A_if(x_i)
$$

为机械求积公式，$a=x_0<x_1<\cdots<x_n=b$ 为求积点，$A_i$ 为求积系数。

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 代数精度
</blockquote>

如果某个求积公式对于次数不超过 $m$ 次的多项式均能精确地成立，而对于次数为 $m+1$ 的多项式不成立，则称该求积公式具有 **$m$ 次代数精度**。

一般地，欲使求积公式具有 $m$ 次代数精度，只要令它对于 $f(x)\in\{1,x,\dots,x^m\}$ 都能准确成立，即
$$
\sum_{i=0}^{n}A_i=b-a\\
\vdots\\
\sum_{i=0}^{n}A_ix_i^m=\frac{1}{m+1}(b^{m+1}-a^{m+1})
$$
这样可以求解线性方程组。

<blockquote style="border-left: 5px solid #b94263; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(185, 66, 110, 0.1)">
    例题
</blockquote>

确定下列求积公式中的常数，使其代数精度尽可能高
$$
\int_{-1}^{1}f(x)\mathrm{d}x\approx A_0f(-1)+A_1f(0)+A_2f(1)
$$

> $$
> \begin{bmatrix}
> 1&1&1\\
> -1&0&1\\
> 1&0&1\\
> -1&0&1\\
> \end{bmatrix}
> \begin{bmatrix}
> A_0\\A_1\\A_2
> \end{bmatrix}=
> \begin{bmatrix}
> 2\\0\\\dfrac{2}{3}\\0
> \end{bmatrix}
> $$
>
> 则精度 $m=3$：
> $$
> \begin{bmatrix}
> A_0\\A_1\\A_2
> \end{bmatrix}=
> \begin{bmatrix}
> \dfrac{1}{3}\\\dfrac{4}{3}\\\dfrac{1}{3}
> \end{bmatrix}
> $$

<blockquote style="border-left: 5px solid #b94263; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(185, 66, 110, 0.1)">
    例题
</blockquote>

推导如下数值积分公式
$$
\int_{0}^{1}f(x)\mathrm{d}x\approx A_0f(x_0)+A_1f'(1)
$$
使它具有尽可能高的代数精度。

> 先算 01 阶：
> $$
> \begin{bmatrix}1&0\\x_0&1\\x^2_0&2\end{bmatrix}
> \begin{bmatrix}A_0\\A_1\end{bmatrix}=
> \begin{bmatrix}1\\1/2\\1/3\end{bmatrix}
> $$
> 则：
> $$
> \int_{0}^{1}f(x)\mathrm{d}x\approx f\left(\frac{3-\sqrt{3}}{3}\right)+\frac{2\sqrt{3}-3}{6}f'(1)
> $$

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 插值型求积公式
</blockquote>

对于给定节点 $a=x_0<x_1<\cdots<x_n=b$，记 $f(x)$ 的 Lagrange 插值多项式为 $L_n(x)$：
$$
f(x)=L_n(x)+R_n(x)
$$
则：
$$
\begin{align*}
I&=\int_{a}^{b}f(x)\mathrm{d}x=\int_{a}^{b}L_n(x)\mathrm{d}x+\int_{a}^{b}R_n(x)\mathrm{d}x\\
&=\int_{a}^{b}\sum_{i=0}^{n}l_i(x)f(x_i)\mathrm{d}x+\int_{a}^{b}\frac{f^{(n+1)}(\xi)}{(n+1)!}\omega_{n+1}(x)\mathrm{d}x\\
&=\sum_{i=0}^{n}\int_{a}^{b}l_i(x)\mathrm{d}x\cdot f(x_i)+\int_{a}^{b}\frac{f^{(n+1)}(\xi)}{(n+1)!}\omega_{n+1}(x)\mathrm{d}x\\
&=\sum_{i=0}^{n}A_if(x_i)+\frac{f^{(n+1)}(\xi)}{(n+1)!}\int_{a}^{b}\omega_{n+1}(x)\mathrm{d}x\\
\end{align*}
$$
*由余项可得，具有 $n+1$ 个节点的插值求积分公式的代数精度至少为 $n$。*

<blockquote style="border-left: 5px solid #b94263; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(185, 66, 110, 0.1)">
    例题
</blockquote>

已知 $f(-1)=-3,f(1)=0,f(2)=4$，求以 $-1,1,2$ 为节点的求积公式, 并判断其代数精度。

> 先求 Lagrange 多项式：
> $$
> l_0(x)=\frac{(x-1)(x-2)}{(-1-1)(-1-2)}=\frac{(x-1)(x-2)}{6}\\
> l_1(x)=\frac{(x+1)(x-2)}{(1+1)(1-2)}=-\frac{(x+1)(x-2)}{2}\\
> l_2(x)=\frac{(x+1)(x-1)}{(2+1)(2-1)}=\frac{(x+1)(x-1)}{3}\\
> $$
> 再求系数：
> $$
> A_0=\int_{-1}^{2}l_0(x)\mathrm{d}x=\frac{3}{4}\\
> A_1=\int_{-1}^{2}l_1(x)\mathrm{d}x=\frac{9}{4}\\
> A_2=\int_{-1}^{2}l_2(x)\mathrm{d}x=0\\
> $$
> 再看代数精度
> $$
> -A_0+A_1+9A_2=\frac{3}{2}\ne\frac{15}{4}\\
> $$
> 为 3。

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 Newton-Cotes 公式
</blockquote>

等距节点下的插值型求积公式称为 Newton-Cotes 公式, 其求积系数称为 Cotes 系数。

令 $x=a+th$，$h=(b-a)/n$
$$
\begin{align*}
A_i=\int_{a}^{b}l_i(x)\mathrm{d}x&=\int_{a}^{b}\prod_{j=0,j\ne i}^{n}\left(\frac{x-x_j}{x_i-x_j}\right)\mathrm{d}x\\
&=\int_{a}^{b}\prod_{j=0,j\ne i}^{n}\left(\frac{x-a-jh}{ih-jh}\right)\mathrm{d}x\\
&=\int_{0}^{n}\prod_{j=0,j\ne i}^{n}\left(\frac{t-j}{i-j}\right)h\mathrm{d}t\quad(t=(x-a)/h)\\
&=\frac{(b-a)}{n}\frac{(-1)^{n-i}}{i!(n-i)!}\int_{0}^{n}\prod_{j=0,j\ne i}^{n}(t-j)\mathrm{d}t\\
\end{align*}
$$
令
$$
\color{red}C^{(n)}_{k}=\frac{(-1)^{n-k}}{k!(n-k)!n}\int_{0}^{n}\prod_{j=0,j\ne k}^{n}(t-j)\mathrm{d}t
$$
为 Cotes 系数。

1. Cotes 系数仅与 $k$ 和 $n$ 有关, 而与被积函数及求积区间无关。
2. $C^{(n)}_{k}=C^{(n)}_{n-k}$
3. $\sum_{k=0}^{n}C^{(n)}_{k}=1$

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    n 阶 Newton-Cotes 公式的代数精度
</blockquote>

对于具有 $n+1$ 个节点的 $n$ 阶 Newton-Cotes 公式, 如果 $n$ 为偶数, 则相应公式的代数精度为 $n+1$；如果 $n$ 为奇数, 则相应求积公式的代数精度为 $n$。

> 取 $f(x)=x^{n+1}$​，则
> $$
> R_n[f]=Ch^{n+2}\int_0^n\prod_{j=0}^n(t-j)\mathrm{d}t
> $$
> 当 $n=2k$ 时，记 $t=u+k$：
> $$
> \begin{align*}
> R_n[f]&=Ch^{n+2}\int_0^n\prod_{j=0}^n(t-j)\mathrm{d}t\\
> &=Ch^{n+2}\int_{-k}^{k}\prod_{j=0}^n(u+k-j)\mathrm{d}u\\
> &=Ch^{n+2}\int_{-k}^{k}\prod_{j=-k}^k(u-j)\mathrm{d}u\\
> &=0
> \end{align*}
> $$
> $\prod_{j=-k}^k(u-j)$ 为奇函数。

### 5.2.2 Newton-Cotes 公式误差分析

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    梯形法则/n = 1 Cotes 公式
</blockquote>

$$
C_0^{(1)}=C_1^{(1)}=\frac{1}{2}
$$

因此：
$$
I=\frac{(b-a)}{2}[f(x_0)+f(x_1)]+\int_{a}^{b}\frac{f^{(2)}(\xi(x))}{2!}(x-a)(x-b)\mathrm{d}x
$$
由推广积分中值定理可得：
$$
\begin{align*}
\int_{a}^{b}\frac{f^{(2)}(\xi(x))}{2!}(x-a)(x-b)\mathrm{d}x&=\frac{f^{(2)}(\xi)}{2!}\int_{a}^{b}(x-a)(x-b)\mathrm{d}x\\
&=-\frac{f^{(2)}(\xi)}{12}(b-a)^3
\end{align*}
$$
因此：
$$
I=\frac{(b-a)}{2}[f(x_0)+f(x_1)]-\frac{f^{(2)}(\xi)}{12}(b-a)^3\quad\xi\in(a,b)
$$

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    Simpson 法则/n = 2 Cotes 公式
</blockquote>

$$
C_1^{(2)}=-\frac{1}{2}\int_{0}^{2}t(t-2)\mathrm{d}t=\frac{2}{3}\\
C_0^{(2)}=C_2^{(2)}=\frac{1}{2}\left(1-\frac{2}{3}\right)=\frac{1}{6}
$$

因此：
$$
I=\frac{(b-a)}{6}[f(x_0)+4f(x_1)+f(x_2)]-\int_{a}^{b}\frac{f^{(3)}(\xi(x))}{3!}(x-a)\left(x-\frac{a+b}{2}\right)(x-b)\mathrm{d}x
$$
推广积分中值定理需要被积函数在区间上不变号，因此我们在 $(a+b)/2$ 处插值两次，注意到 Simpson 公式的代数精度为 3，因此余项变为：
$$
\begin{align*}
R[f]&=\int_{a}^{b}\frac{f^{(4)}(\xi(x))}{4!}(x-a)\left(x-\frac{a+b}{2}\right)^2(x-b)\mathrm{d}x\\
&=\frac{f^{(4)}(\xi)}{4!}\int_{a}^{b}(x-a)\left(x-\frac{a+b}{2}\right)^2(x-b)\mathrm{d}x\\
&=-\frac{f^{(4)}(\xi)}{4!}\frac{(b-a)^5}{120}\\
&=-\frac{f^{(4)}(\xi)}{90}\frac{(b-a)^5}{2^5}\\
\end{align*}
$$

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    一般 Newton-Cotes 公式的余项
</blockquote>

对于 $n$ 为偶数，$f\in C^{n+2}[a,b]$，则：
$$
R_n[f]=C_nh^{n+3}f^{(n+2)}(\eta)\quad\eta\in[a,b]\\
C_n=\frac{1}{(n+2)!}\int_0^nu^2(u-1)\cdots(u-n)\mathrm{d}u
$$
对于 $n$ 为奇数，$f\in C^{n+1}[a,b]$，则：

$$
R_n[f]=C_nh^{n+2}f^{(n+1)}(\eta)\quad\eta\in[a,b]\\
C_n=\frac{1}{(n+1)!}\int_0^nu(u-1)\cdots(u-n)\mathrm{d}u
$$

<blockquote style="border-left: 5px solid #b94263; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(185, 66, 110, 0.1)">
    例题
</blockquote>

求
$$
\int_0^1f(x)\mathrm{d}x=A_0f(0)+A_1f(1)+A_2f'(0)+A_3f'(1)
$$

和余项。

> 带入 $f(x)=1,x,x^2,x^3$：
> $$
> \begin{bmatrix}1&1&0&0\\0&1&1&1\\0&1&0&2\\0&1&0&3\end{bmatrix}\begin{bmatrix}A_0\\A_1\\A_2\\A_3\end{bmatrix}=\begin{bmatrix}1\\1/2\\1/3\\1/4\end{bmatrix}
> $$
> 得到：
> $$
> \begin{bmatrix}A_0\\A_1\\A_2\\A_3\end{bmatrix}=\begin{bmatrix}1/2\\1/2\\1/12\\-1/12\end{bmatrix}
> $$
> 验算得到代数精度为 3。
>
> 余项为：
> $$
> \int_0^1\frac{f^{(4)}(\xi(x))}{4!}x^2(x-1)^2\mathrm{d}x=\frac{f^{(4)}(\xi)}{720}
> $$

### 5.2.3 Newton-Cotes 公式的计算稳定性问题

假定初始数据 $f(x_k)$ 有舍入误差, 设 $f(x_k)\approx\tilde{f}(x_k)$，$\epsilon_k=\max_{0\le k\le n}|f(x_k)-\tilde{f}(x_k)|$
$$
\left|(b-a)\sum_{k=0}^{n}c_{k}^{(n)}f(x_k)-(b-a)\sum_{k=0}^{n}c_{k}^{(n)}\tilde{f}(x_k)\right|\le(b-a)\epsilon\sum_{k=0}^{n}|c_{k}^{(n)}|
$$
当 Cotes 系数全为正时候：
$$
\left|(b-a)\sum_{k=0}^{n}c_{k}^{(n)}f(x_k)-(b-a)\sum_{k=0}^{n}c_{k}^{(n)}\tilde{f}(x_k)\right|\le(b-a)\epsilon\sum_{k=0}^{n}c_{k}^{(n)}=(b-a)\epsilon
$$
当 Cotes 系数出现负数时候：
$$
\sum_{k=0}^{n}|c_{k}^{(n)}|>\sum_{k=0}^{n}c_{k}^{(n)}=1
$$
误差会增大，数值计算可能不稳定。

### 5.2.4 Composite Newton-Cotes Formulas 复合求积公式

设将积分区间 $[a, b]$ 划分为 $n$ 等分，小区间用低阶的 Newton-Cotes 公式计算。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    梯形法则
</blockquote>

小区间上的积分用梯形公式计算
$$
\begin{align*}
I&=\int_a^bf(x)\mathrm{d}x=\sum_{k=0}^{n-1}\int_{x_k}^{x_{k+1}}f(x)\mathrm{d}x\\
&=T_n\approx\sum_{k=0}^{n-1}\frac{h}{2}[f(x_k)+f(x_{k+1})]\\
&=\frac{h}{2}[f(a)+f(b)+2\sum_{k=1}^{n-1}f(x_k)]\\
&=\frac{b-a}{n}[f(a)+f(b)+2\sum_{k=1}^{n-1}f(x_k)]
\end{align*}
$$
余项为：
$$
\begin{align*}
R_T[f]&=I-T_n=\sum_{k=0}^{n-1}\left(\int_{x_k}^{x_{k+1}}f(x)\mathrm{d}x-\frac{h}{2}[f(x_k)+f(x_{k+1})]\right)\\
&=\sum_{k=0}^{n-1}-\frac{h^3f''(\eta_k)}{12}\quad(x_k\le\eta_k\le x_{k+1})\\
&=-\frac{h^3}{12}nf''(\xi)=-\frac{(b-a)^3}{12n^2}f''(\xi)\quad(x_0\le\xi\le x_n)
\end{align*}
$$
最后一个等号为推广中值定理。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    复合中矩形公式
</blockquote>

当精度不够时，可以增加节点，减少步长。
$$
\begin{align*}
T_{2n}&=\frac{(b-a)}{2n}[f(a)+f(b)+2\sum_{k=1}^{2n-1}f(x_{\frac{1}{2}k})]\\
&=\frac{h}{4}[f(a)+f(b)+2\sum_{k=1}^{n-1}f(x_{\frac{1}{2}k})]+\frac{h}{2}\sum_{k=0}^{n-1}f(x_{\frac{1}{2}+k})\\
&=\frac{1}{2}T_n+\frac{1}{2}M_n
\end{align*}
$$
其中
$$
M_n=h\sum_{k=0}^{n-1}f(x_{\frac{1}{2}+k})
$$
为复合中矩形公式。
$$
R_T[f]=-\frac{(b-a)^3}{12n^2}f''(\xi_1)\\
R_{2T}[f]=-\frac{(b-a)^3}{48n^2}f''(\xi_2)
$$
因此：
$$
\frac{I-T_{2n}}{I-T_{n}}\approx\frac{1}{4}
$$
因此：
$$
I\approx\frac{4}{3}T_{2n}-\frac{1}{3}T_n
$$

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    Simpson 法则
</blockquote>

$$
\begin{align*}
I&=\int_a^bf(x)\mathrm{d}x=\sum_{k=0}^{n-1}\int_{x_k}^{x_{k+1}}f(x)\mathrm{d}x\\
&\approx \sum_{k=0}^{n-1}\frac{h}{6}[f(x_k)+4f(x_{k+\frac{1}{2}})+f(x_{k+1})]-\sum_{k=0}^{n-1}\frac{1}{90}\left(\frac{h}{2}\right)^5f^{(4)}(\eta_k)\quad(x_k\le\eta_k\le x_{k+1})\\
&=\sum_{k=0}^{n-1}\frac{h}{6}[f(x_k)+4f(x_{k+\frac{1}{2}})+f(x_{k+1})]-\frac{b-a}{2880}h^4f^{(4)}(\eta)\quad(x_0\le\eta\le x_{n})\\
&=\sum_{k=0}^{n-1}\frac{b-a}{6n}[f(x_k)+4f(x_{k+\frac{1}{2}})+f(x_{k+1})]-\frac{(b-a)^5}{2880n^4}f^{(4)}(\eta)\quad(x_0\le\eta\le x_{n})
\end{align*}
$$

而：
$$
\begin{align*}
S_n&=\sum_{k=0}^{n-1}\frac{h}{6}[f(x_k)+4f(x_{k+\frac{1}{2}})+f(x_{k+1})]\\
&=\sum_{k=0}^{n-1}\frac{h}{6}[f(x_k)+f(x_{k+1})]+\sum_{k=0}^{n-1}\frac{2h}{3}f(x_{k+\frac{1}{2}})\\
&=\frac{1}{3}T_n+\frac{2}{3}M_n
\end{align*}
$$
由 $T_{2n}=\dfrac{1}{2}T_n+\dfrac{1}{2}M_n$ 得到 $I\approx\dfrac{4}{3}T_{2n}-\dfrac{1}{3}T_n$。

