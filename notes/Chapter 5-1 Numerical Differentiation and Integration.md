# Chapter 5 数值微分和积分 Numerical Differentiation and Integration

## 5.1 Numerical Differentiation 数值微分

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定理 推广中值定理
</blockquote>

设 $f\in C[a,b]$，$x_1,\dots,x_n\in[a,b]$，$a_1,\dots,a_n>0$，则 $\exist c\in[a,b]$，满足：
$$
(a_1+\dots+a_n)f(c)=a_1f(x_1)+\dots+a_nf(x_n)
$$

> 设 $f(x_i),f(x_j)$ 是 $f(x_1),\dots,f(x_n)$ 之间的最小和最大值，则必然有：
> $$
> \sum_{k=1}^{n}a_kf(x_i)\le\sum_{k=1}^{n}a_kf(x_k)\le\sum_{k=1}^{n}a_kf(x_j)\\
> f(x_i)\le\frac{\sum_{k=1}^{n}a_kf(x_k)}{\sum_{k=1}^{n}a_k}\le f(x_j)
> $$
> 根据介值定理，$f(x)$ 能在 $x_i$ 与 $x_j$ 之间取到 $\xi$，使得
> $$
> \color{red}f(\xi)=\frac{\sum_{k=1}^{n}a_kf(x_k)}{\sum_{k=1}^{n}a_k}
> $$

### 5.1.1 Finite Difference Formulas

设 $f\in C^{n+1}$，对其在 $x+h$ 处 Taylor 展开有：
$$
f(x+h)=f(x)+hf'(x)+\frac{h^2}{2!}f''(x)+\dots+\frac{h^n}{n!}f^{(n)}(x)+\frac{h^{n+1}}{(n+1)!}f^{(n+1)}(\xi)
$$
$\xi$ 位于 $x$ 与 $x+h$ 之间。

**$n$ 阶近似：误差是 $O(h^n)$。**

> Two-point forward-difference formula 二点前向差分公式
>
> 设 $f\in C^{2}$，
> $$
> \color{red}f'(x)=\frac{f(x+h)-f(x)}{h}-\frac{h}{2}f''(\xi)
> $$
> $\xi$ 位于 $x$ 与 $x+h$ 之间，且依赖于 $h$。该公式是一阶近似。

> Three-point centered-difference formula 三点中心差分公式
>
> 设 $f\in C^{3}$，由
> $$
> f(x+h)=f(x)+hf'(x)+\frac{h^2}{2!}f''(x)+\frac{h^3}{3!}f'''(\xi_1)\\
> f(x-h)=f(x)-hf'(x)+\frac{h^2}{2!}f''(x)-\frac{h^3}{3!}f'''(\xi_2)
> $$
> 相减可得：
> $$
> f'(x)=\frac{f(x+h)-f(x-h)}{2h}-\frac{h^2}{12}f'''(\xi_1)-\frac{h^2}{12}f'''(\xi_2)
> $$
> $x-h<\xi_2<x<\xi_1<x+h$。
>
> 由推广中值定理可得到：
> $$
> \color{red}f'(x)=\frac{f(x+h)-f(x-h)}{2h}-\frac{h^2}{6}f'''(\xi)
> $$
> 其中 $x-h<\xi<x+h$。该公式是二阶近似。

> Three-point centered-difference formula for second derivative 二阶导数的三点中心差分公式
>
> 设 $f\in C^{4}$，由
> $$
> f(x+h)=f(x)+hf'(x)+\frac{h^2}{2!}f''(x)+\frac{h^3}{3!}f'''(x)+\frac{h^4}{4!}f^{(4)}(\xi_1)\\
> f(x-h)=f(x)-hf'(x)+\frac{h^2}{2!}f''(x)-\frac{h^3}{3!}f'''(x)+\frac{h^4}{4!}f^{(4)}(\xi_2)
> $$
> 相加，并应用推广中值定理可得：
> $$
> \color{red}f''(x)=\frac{f(x+h)-2f(x)+f(x-h)}{h^2}-\frac{h^2}{12}f^{(4)}(\xi)
> $$
> 其中 $x-h<\xi<x+h$。该公式是二阶近似。

### 5.1.2 Rouding Error 舍入误差

分析三点中心差分公式的舍入误差：

设 $f(x+h)$ 在机器中的存储的形式为 $\hat{f}(x+h)$，且误差是与机器精度同阶的，即 $\hat{f}(x+h)=f(x+h)+\epsilon_1$ 因此：
$$
\begin{align*}
f'(x)_{\text{correct}}-f'(x)_{\text{machine}}&=f'(x)_{\text{correct}}-\frac{\hat{f}(x+h)-\hat{f}(x-h)}{2h}\\
&=f'(x)_{\text{correct}}-\frac{f(x+h)+\epsilon_1-f(x-h)-\epsilon_2}{2h}\\
&=\left[f'(x)_{\text{correct}}-\frac{f(x+h)-f(x-h)}{2h}\right]+\frac{\epsilon_2-\epsilon_1}{2h}\\
&=\left[f'(x)_{\text{correct}}-f'(x)_{\text{formula}}\right]+\text{error}_{\text{rounding}}
\end{align*}
$$
因此误差的绝对值的上界就是：
$$
\begin{align*}
|f'(x)_{\text{correct}}-f'(x)_{\text{machine}}|&=|\left(f'(x)_{\text{correct}}-f'(x)_{\text{formula}}\right)+\text{error}_{\text{rounding}}|\\
&\le|f'(x)_{\text{correct}}-f'(x)_{\text{formula}}|+\left|\frac{\epsilon_2-\epsilon_1}{2h}\right|\\
&\le\frac{h^2}{6}f'''(\xi)+\frac{2\epsilon}{2h}\\
&=\frac{h^2}{6}f'''(\xi)+\frac{\epsilon}{h}
\end{align*}
$$
其中 $\epsilon$ 为机器精度，$x-h<\xi<x+h$。

### 5.1.3 Richardson Extrapolation 外推法

假设 $F(h)$ 逼近一个给定的量 $F^*$ ($F^*$ 与 $h$ 无关) ，其余项为：
$$
F^*-F(h)=\sum_{k=1}^{\infty}\alpha_kh^{p_k}
$$
其中 $p_k,\alpha_k$ 为与 $h$ 无关的常数，$\alpha_k\ne0,k>0$。定义外推公式：
$$
\begin{align*}
F_1(h)&=F(h)\\
F_{m+1}(h)&=\frac{F_m(qh)-q^{p_m}F_m(h)}{1-q^{p_m}},q\in(0,1)
\end{align*}
$$

> 以 $q=1/2$ 为例，由于
> $$
> F^*-F(h)=\sum_{k=1}^{\infty}\alpha_kh^{p_k}
> $$
>
> 步长减半：
> $$
> F^*-F\left(\frac{h}{2}\right)=\sum_{k=1}^{\infty}\alpha_k\left(\frac{h}{2}\right)^{p_k}
> $$
> 提取右端第一项并整理形式：
> $$
> 2^{p_1}F^*-2^{p_1}F\left(\frac{h}{2}\right)=\alpha_1h^{p_1}+\sum_{k=2}^{\infty}\alpha_k\frac{h^{p_k}}{2^{p_k-p_1}}=\alpha_1h^{p_1}+\sum_{k=2}^{\infty}\frac{\alpha_k}{2^{p_k-p_1}}h^{p_k}
> $$
> 减去最上面的式子，消去右端第一项：
> $$
> \begin{gather*}
> (2^{p_1}-1)F^*-\left[2^{p_1}F\left(\frac{h}{2}\right)-F(h)\right]=\sum_{k=2}^{\infty}\frac{\alpha_k}{2^{p_k-p_1}}h^{p_k}\\
> F^*-\frac{2^{p_1}F\left(\frac{h}{2}\right)-F(h)}{(2^{p_1}-1)}=\sum_{k=1}^{\infty}\beta_kh^{p_k}\\
> F^*-\frac{F\left(\frac{h}{2}\right)-\left(\frac{1}{2}\right)^{p_1}F(h)}{(1-\left(\frac{1}{2}\right)^{p_1})}=\sum_{k=1}^{\infty}\beta_kh^{p_k}\\
> \end{gather*}
> $$
> 可以发现误差的形式保持。

<blockquote style="border-left: 5px solid #b94263; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(185, 66, 110, 0.1)">
    例题
</blockquote>
设 $I(h)$ 逼近 $I$ 的余项为:
$$
I-I(h)=\alpha_1h^2+\alpha_2h^4+\cdots
$$
其中 $\alpha_k$ 为与 $h$ 无关的常数，且 $I(0.1)=0.99,I(0.04)=1.20$，在此基础上用外推法得到 $I$ 的更高精度的近似值为多少？

> $$
> \begin{gather*}
> I-I(0.1)=\alpha_10.1^2+\alpha_20.1^4+\cdots\\
> I-I(0.04)=\alpha_10.04^2+\alpha_20.04^4+\cdots\\
> \end{gather*}
> $$
>
> 消去 $\alpha_1$:
> $$
> \begin{gather*}
> 5.25I-6.25I(0.04)+I(0.1)=0.1^2(0.1^2-0.04^2)\alpha_2+\cdots\\
> \frac{6.25I(0.04)-I(0.1)}{5.25}=1.19238
> \end{gather*}
> $$

<blockquote style="border-left: 5px solid #b94263; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(185, 66, 110, 0.1)">
    例题
</blockquote>

外推二阶精度的三点中心差分公式。

> 已知：
> $$
> F(h)=f'(x)=\frac{f(x+h)-f(x-h)}{2h}-\frac{h^2}{6}f'''(\xi)
> $$
> 步长减半：
> $$
> F(h/2)=f'(x)=\frac{f(x+h/2)-f(x-h/2)}{h}-\frac{h^2}{24}f'''(\xi)
> $$
> 消去三阶导数：
> $$
> \begin{align*}
> f'(x)&=\frac{2^2F(h/2)-F(h)}{2^2-1}\\
> &=\frac{-f(x+h)+8f(x+h/2)-8f(x-h/2)+f(x-h)}{6h}
> \end{align*}
> $$
> 由外推法可知该公式的精度至少是 3 阶。

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
\begin{gather*}
\sum_{i=0}^{n}A_i=b-a\\
\vdots\\
\sum_{i=0}^{n}A_ix_i^m=\frac{1}{m+1}(b^{m+1}-a^{m+1})
\end{gather*}
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
> 则精度 $m=3$：yt
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
\begin{gather*}
R_T[f]=-\frac{(b-a)^3}{12n^2}f''(\xi_1)\\
R_{2T}[f]=-\frac{(b-a)^3}{48n^2}f''(\xi_2)\\
\end{gather*}
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

