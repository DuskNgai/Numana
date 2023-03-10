# Chapter 5 数值微分和积分 Numerical Differentiation and Integration

## 5.3 Romberg 积分

将 Richardson 外推和复合梯形公式结合

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    Romberg 积分
</blockquote>

设 $f(x)\in C^{2m+2}[a,b]$，则复合梯形公式：
$$
T(h)=I(f)+\sum_{k=1}^{m}d_kh^{2k}+d_{md+1}(h)h^{2m+2}
$$
其中：$T(h)=T_n$，$I(f)=\int_a^bf(x)\mathrm{d}x,d_j(j\in\{1,\dots,n\})$ 为与 $h$ 无关的常数，$d_{n+1}(h)$ 为 $h$ 的函数。

> 证明略。

记：
$$
T_i^0=T\left(\frac{h}{2^i}\right)
$$
则：
$$
T_i^0=I(f)+\sum_{k=1}^{n}d_k\left(\frac{h}{2^i}\right)^{2k}+d_{n+1}\left(\frac{h}{2^i}\right)\left(\frac{h}{2^i}\right)^{2n+2}
$$
此时就可以外推了：
$$
T_i^1=\frac{4T_i^0-T_{i-1}^0}{4-1}
$$
Romberg 算法的递推表：
$$
\begin{matrix}
T_0^0\\
T_1^0&T_1^1\\
T_2^0&T_2^1&T_2^2\\
\vdots&\ddots&\ddots&\ddots\\
T_n^0&T_n^1&\cdots&\cdots&T_n^n
\end{matrix}
$$
其中：
$$
\begin{align*}
T_0^0&=[f(a)+f(b)]\frac{h}{2}\quad(h=b-a)\\
T_i^0&=\frac{1}{2}T_{i-1}^0+\frac{h}{2^i}\sum_{k=0}^{i-1}f\left(x_{k+\frac{1}{2}}\right)\\
T_i^j&=\frac{4^jT_{i}^{j-1}-T_{i-1}^{j-1}}{4^j-1}
\end{align*}
$$

第二个式子的作用是采样，第三个式子的作用是外推

## 5.5 Gauss 积分

形如
$$
\int_a^bf(x)\mathrm{d}x=\sum_{k=0}^nA_kf(x_k)
$$
的求积公式，有 $2n+2$ 个待定参数 $A_k,x_k$。选取合适的参数可以使得求积公式具有 $2n+1$ 次代数精度。

### 5.5.1 Hermite 积分

已知在节点 $a\le x_0<x_1<\cdots<x_n\le b$ 上，节点处的函数值和导数值 $f(x_i)$ 和 $f'(x_i)$ 已知。满足上述插值条件的 Hermite 多项式 $H_{2n+1}(x)$ 存在且唯一，且
$$
H_{2n+1}=\sum_{i=0}^n\alpha_i(x)f(x_i)+\sum_{i=0}^n\beta_i(x)f'(x_i)\\
\alpha_i(x)=[1-2(x-x_i)l'_i(x_i)]l^2_i(x)\quad\beta_i(x)=(x-x_i)l^2_i(x)
$$
则对应的带权 $\rho(x)$ Hermite 积分为：
$$
\begin{align*}
\int_a^b\rho(x)f(x)\mathrm{d}x&\approx\sum_{i=0}^nA_if(x_i)+\sum_{i=0}^nB_if'(x_i)\\
A_i&=\int_a^b\rho(x)[1-2(x-x_i)l'_i(x_i)]l^2_i(x)\mathrm{d}x\\
B_i&=\int_a^b\rho(x)(x-x_i)l^2_i(x)\mathrm{d}x
\end{align*}
$$
由于 Hermite 多项式的误差为：
$$
R_{2n+1}(x)=f(x)-H_{2n+1}(x)=\frac{f^{(2n+2)}(\xi)}{(2n+2)!}\omega^2_{n+1}(x)
$$
那么积分的误差也应该为：
$$
E_{2n+1}(x)=\int_a^b\rho(x)R_{2n+1}(x)\mathrm{d}x=\int_a^b\frac{f^{(2n+2)}(\xi)}{(2n+2)!}\rho(x)\omega^2_{n+1}(x)\mathrm{d}x
$$
这里权函数 $\rho(x)$ 在该区间上是非负的，由第二积分中值定理可得
$$
E_{2n+1}(x)=\frac{f^{(2n+2)}(\xi)}{(2n+2)!}\int_a^b\rho(x)\omega^2_{n+1}(x)\mathrm{d}x
$$
因此其代数精度为 $2n+1$，但不是 Gauss 积分，因为需要导数值。不过只要选取合适的 $x_0,\dots,x_n$，就可以使得 $B_i(x)\equiv0$，从而提高精度的同时减少计算量。

### 5.5.2 Gauss 积分

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    Gauss 点
</blockquote>

含有 $n+1$ 个求积点 $x_0,\dots,x_n$，而代数精度为 $2n+1$ 的求积公式称为 Gauss 型求积公式，$x_0,\dots,x_n$ 为 Gauss 点。由于上式为非线性方程组，实际计算非常困难。

<blockquote style="border-left: 5px solid #b94263; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(185, 66, 110, 0.1)">
    例题
</blockquote>

$$
\int_{-1}^1f(x)\mathrm{d}x=Af\left(-\frac{\sqrt{3}}{3}\right)+Bf\left(\frac{\sqrt{3}}{3}\right)
$$

> $$
> \begin{bmatrix}
> 1&1\\
> -\frac{\sqrt{3}}{3}&\frac{\sqrt{3}}{3}\\
> \end{bmatrix}
> \begin{bmatrix}
> A\\B
> \end{bmatrix}=
> \begin{bmatrix}
> 2\\0\\
> \end{bmatrix}
> $$
>
> 得到 $A=B=1$。且代数精度为 3。则是 Gauss 型求积公式。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    从 Hermite 积分发展得来
</blockquote>

设以 $x_0,\dots,x_n$ 为节点的 Hermite 积分能够约化成以下形式
$$
\int_a^b\rho(x)f(x)\mathrm{d}x\approx\sum_{i=0}^nA_if(x_i)+E(x)
$$
的充要条件是 $n+1$ 次的多项式 $\omega_{n+1}(x)$ 与任意不超过 $n$ 次的多项式在区间 $[a,b]$ 上带权 $\rho(x)$ 正交，其中
$$
E_{2n+1}(x)=\frac{f^{(2n+2)}(\xi)}{(2n+2)!}\int_a^b\rho(x)\omega^2_{n+1}(x)\mathrm{d}x
$$
此时 Hermite 积分转变为 Gauss 型求积公式，代数精度保持 $2n+1$。

> **必要性**
>
> 由于 $B_i(x)$​ 可以表示为：
> $$
> B_i=\int_a^b\rho(x)(x-x_i)l^2_i(x)\mathrm{d}x=\frac{1}{\omega'_{n+1}(x_i)}\int_a^b\rho(x)l_i(x)\omega_{n+1}(x)
> \mathrm{d}x
> $$
> $\omega_{n+1}(x)$ 与任意不超过 $n$ 次的多项式在区间 $[a,b]$ 上带权 $\rho(x)$ 正交，且 $l_i(x)$ 为 $n$ 次多项式，因此 $B_i\equiv0$，可以约化成上述形式。
>
> **充分性**
>
> 将 $f(x)$ 分解为 $u_n(x)\omega_{n+1}(x)+v_n(x)$，其中 $u_n(x)$ 和 $v_n(x)$ 为不超过 $n$ 次的多项式。由于正交，我们有：
> $$
> \begin{align*}
> \int_a^b\rho(x)f(x)\mathrm{d}x&=\int_a^b\rho(x)u_n(x)\omega_{n+1}(x)\mathrm{d}x+\int_a^b\rho(x)v_n(x)\mathrm{d}x\\
> &=0+\int_a^b\rho(x)v_n(x)\mathrm{d}x\\
> &=\sum_{i=0}^{n}A_iv_n(x_i)
> \end{align*}
> $$
> 由于 $\omega_{n+1}(x_i)=0$，所以 $f(x_i)=u_n(x_i)\omega_{n+1}(x_i)+v_n(x_i)=v_n(x_i)$。
>
> 因此
> $$
> \int_a^b\rho(x)f(x)\mathrm{d}x=\sum_{i=0}^{n}A_iv_n(x_i)=\sum_{i=0}^{n}A_if(x_i)
> $$

由于 $B_i\equiv0$，则
$$
\begin{align*}
A_i&=\int_a^b\rho(x)[1-2(x-x_i)l'_i(x_i)]l^2_i(x)\mathrm{d}x\\
&=\int_a^b\rho(x)l^2_i(x)\mathrm{d}x-2l'_i(x_i)\int_a^b\rho(x)\frac{\omega_{n+1}(x)}{\omega'_{n+1}(x_i)}l_i(x)\mathrm{d}x\\
&=\int_a^b\rho(x)l^2_i(x)\mathrm{d}x-2l'_i(x_i)B_i\\
&=\int_a^b\rho(x)l^2_i(x)\mathrm{d}x
\end{align*}
$$

<blockquote style="border-left: 5px solid #b94263; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(185, 66, 110, 0.1)">
    例题
</blockquote>

构造 Gauss 型求积公式
$$
\int_0^1xf(x)\mathrm{d}x=Af(x_0)+Bf(x_1)
$$

> 设 $[0,1]$ 上带权 $\rho(x)=x$ 正交，且首系数为 $1$ 的多项式系为 $\{P_n(x)\}_0^{\infty}$。设 $P_2(x)=x^2+ax+b$，则 $P_2(x)$ 由下式可得：
> $$
> \int_0^1x\cdot 1\cdot P_2(x)\mathrm{d}x=0\\
> \int_0^1x\cdot x\cdot P_2(x)\mathrm{d}x=0
> $$
> 得到 $a=-6/5,b=3/10$。则两根为：
> $$
> \frac{6\pm\sqrt{6}}{10}
> $$

### 5.5.2 Gauss 求积法的稳定性和收敛性

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    Gauss 求积法的稳定性和收敛性
</blockquote>

Gauss 求积公式中的求积系数全大于零，即 $A_k>0$。

> 设 $l_j(x)$ 是以 Gauss 点为插值点的拉格朗日插值基函数, 由于 Gauss 求积公式的代数精度为 $2n+1$，从而：
> $$
> 0<\int_a^b\rho(x)l^2_j(x)\mathrm{d}x=\sum_{k=0}^{n}A_kl_j^2(x_k)=A_j
> $$

这表明 Gauss 求积方法的数值计算是稳定。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    Gauss 型求积公式余项
</blockquote>

若 $f\in C^{2n+2}[a,b]$，则 Gauss 求积公式具有余项：
$$
\int_a^b\rho(x)f(x)\mathrm{d}x-\sum_{k=0}^nA_kf(x_k)=\frac{f^{(2n+2)}(\xi)}{(2n+2)!}\int_a^b\rho(x)\omega^2_{n+1}(x)\mathrm{d}x
$$

> 构造 Hermite 插值多项式 $H_{2n+1}(x)$：
> $$
> f(x)-H_{2n+1}(x)=\frac{f^{(2n+2)}(\xi)}{(2n+2)!}\omega^2_{n+1}(x)
> $$
> 则得到余项。

### 5.5.3 Gauss-Legendre 求积公式

若 $\rho(x)=1$，$[a,b]=[-1,1]$，则求积点可取为 Legendre 多项式的零点。

对于一般区间 $[a,b]\ne[-1,1]$ 的情形，可通过坐标变换
$$
x=\frac{b+a}{2}+\frac{b-a}{2}t
$$
将积分 $\int_a^bf(x)\mathrm{d}x$ 变换为：
$$
\frac{b-a}{2}\int_{-1}^1f\left(\frac{b+a}{2}+\frac{b-a}{2}t\right)\mathrm{d}t
$$

### 5.5.4 Gauss-Chebyshev 求积公式

若 $\rho(x)=1/\sqrt{1-x^2}$，$[a,b]=[-1,1]$，则求积点可取为 Chebyshev 多项式的零点：
$$
x_k=\cos\frac{2k+1}{2(n+1)}\pi
$$
求积系数为：
$$
A_k=\frac{\pi}{n+1}
$$
