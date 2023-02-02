## Question 1

> 确定下列求积公式中的待定参数，使其代数精度尽量高，并指明所构造出的求积公式所具有的代数精度。
>
> (1)
> $$
> \int_{-2h}^{2h}f(x)\mathrm{d}x\approx A_{-1}f(-h)+A_0f(0)+A_1f(h)
> $$
> (2)
> $$
> \int_{0}^{1}xf(x)\mathrm{d}x\approx Af(0)+Bf(1)+Cf'(0)+Df'(1)
> $$

(1)

先算 012 阶：
$$
\begin{bmatrix}1&1&1\\-h&0&h\\h^2&0&h^2\\\end{bmatrix}
\begin{bmatrix}A_{-1}\\A_0\\A_1\end{bmatrix}=
\begin{bmatrix}4h\\0\\\frac{16}{3}h^3\end{bmatrix}
$$
解得：
$$
A_{-1}=\frac{8}{3}h,A_0=-\frac{4}{3}h,A_1=\frac{8}{3}h
$$
求积公式为：
$$
\int_{-2h}^{2h}f(x)\mathrm{d}x\approx \frac{8}{3}hf(-h)-\frac{4}{3}hf(0)+\frac{8}{3}hf(h)
$$
验算可得：
$$
-h^3A_{-1}+h^3A_1=0\\
h^4A_{-1}+h^4A_1=\frac{16}{3}h^5\ne\frac{64}{5}h^5
$$
因此代数精度为 3。

(2)

先算 0123 阶：
$$
\begin{bmatrix}1&1&0&0\\0&1&1&1\\0&1&0&2\\0&1&0&3\end{bmatrix}
\begin{bmatrix}A\\B\\C\\D\end{bmatrix}=
\begin{bmatrix}\frac{1}{2}\\\frac{1}{3}\\\frac{1}{4}\\\frac{1}{5}\end{bmatrix}
$$
解得：
$$
A=\frac{3}{20},B=\frac{7}{20},C=\frac{1}{30},D=-\frac{1}{20}
$$
求积公式为：
$$
\int_{0}^{1}f(x)\mathrm{d}x\approx \frac{1}{2}f(0)+\frac{1}{2}f(1)+\frac{1}{12}f'(0)-\frac{1}{12}f'(1)
$$
验算可得：
$$
B+4D=\frac{3}{20}\ne\frac{1}{6}
$$
因此代数精度为 3。

## Question 2

> 推导下列三种矩形求积公式：
>
> (1)
> $$
> \int_{a}^{b}f(x)\mathrm{d}x=f(a)(b-a)+\frac{f'(\eta)}{2}(b-a)^2
> $$
> (2)
> $$
> \int_{a}^{b}f(x)\mathrm{d}x=f(b)(b-a)-\frac{f'(\eta)}{2}(b-a)^2
> $$
> (3)
> $$
> \int_{a}^{b}f(x)\mathrm{d}x=f\left(\frac{a+b}{2}\right)(b-a)+\frac{f''(\eta)}{24}(b-a)^3
> $$

这里需要假设 $f(x)$ 在 $[a,b]$ 上连续，在 $(a,b)$ 上可导。

(1)

在 $a$ 处 Taylor 展开到 $x$ 点：
$$
f(x)=f(a)+f'(\xi)(x-a)\quad\xi\in(a,x)
$$
由于 $x-a$ 不变号，因此：
$$
\begin{align*}
&\quad\ \int_{a}^{b}[f(x)-f(a)]\mathrm{d}x\\
&=\int_{a}^{b}f'(\xi)(x-a)\mathrm{d}x\\
&=\frac{f'(\eta)}{2}(b-a)^2\quad(\eta\in(a,b))\\
\end{align*}
$$
因此：
$$
\int_{a}^{b}f(x)\mathrm{d}x=f(a)(b-a)+\frac{f'(\eta)}{2}(b-a)^2
$$
(2)

在 $b$ 处 Taylor 展开到 $x$ 点：
$$
f(x)=f(b)-f'(\xi)(x-b)\quad\xi\in(x,b)
$$
由于 $b-x$ 不变号，因此：
$$
\begin{align*}
&\quad\ \int_{a}^{b}[f(x)-f(b)]\mathrm{d}x\\
&=\int_{a}^{b}f'(\xi)(x-b)\mathrm{d}x\\
&=-\frac{f'(\eta)}{2}(b-a)^2\quad(\eta\in(a,b))\\
\end{align*}
$$
因此：
$$
\int_{a}^{b}f(x)\mathrm{d}x=f(b)(b-a)-\frac{f'(\eta)}{2}(b-a)^2
$$
(3)

在 $\left(\frac{a+b}{2}\right)$ 处 Taylor 展开到 $x$ 点：
$$
f(x)=f\left(\frac{a+b}{2}\right)+f'\left(\frac{a+b}{2}\right)\left(x-\frac{a+b}{2}\right)+\frac{f''(\eta)}{2}\left(x-\frac{a+b}{2}\right)^2
$$
因此：
$$
\begin{align*}
&\quad\ \int_{a}^{b}f(x)\mathrm{d}x-f\left(\frac{a+b}{2}\right)(b-a)\\
&=\int_{a}^{b}f(x)-f\left(\frac{a+b}{2}\right)\mathrm{d}x\\
&=\int_{a}^{b}f'\left(\frac{a+b}{2}\right)\left(x-\frac{a+b}{2}\right)+\frac{f''(\eta)}{2}\left(x-\frac{a+b}{2}\right)^2\mathrm{d}x\\
&=\frac{f''(\eta)}{6}\left(x-\frac{a+b}{2}\right)^3\biggl|_a^b\\
&=\frac{f''(\eta)}{24}(b-a)^3
\end{align*}
$$
因此：
$$
\int_{a}^{b}f(x)\mathrm{d}x=f\left(\frac{a+b}{2}\right)(b-a)+\frac{f''(\eta)}{24}(b-a)^3
$$

## Question 3

> 已知 $x_0=1/4,x_1=1/2,x_2=3/4$。
>
> (1) 推导以上述 3 点作为求积节点在 $[0,1]$ 上的插值型求积公式。
>
> (2) 指明求积公式具有的代数精度。
>
> (3) 用所求公式计算 $\int_0^1x^2\mathrm{d}x$。

(1)

设积分公式为：
$$
\int_0^1f(x)\mathrm{d}x=A_0f\left(\frac{1}{4}\right)+A_1f\left(\frac{1}{2}\right)+A_2f\left(\frac{3}{4}\right)
$$
先算 012 阶：
$$
\begin{bmatrix}1&1&1\\\frac{1}{4}&\frac{1}{2}&\frac{3}{4}\\\frac{1}{16}&\frac{1}{4}&\frac{9}{16}\\\end{bmatrix}\begin{bmatrix}A_0\\A_1\\A_2\end{bmatrix}=\begin{bmatrix}1\\\frac{1}{2}\\\frac{1}{3}\end{bmatrix}
$$
解得：
$$
A_0=\frac{2}{3},A_1=-\frac{1}{3},A_2=\frac{2}{3}
$$
则积分公式为：
$$
\int_0^1f(x)\mathrm{d}x=\frac{2}{3}f\left(\frac{1}{4}\right)-\frac{1}{3}f\left(\frac{1}{2}\right)+\frac{2}{3}f\left(\frac{3}{4}\right)
$$
(2)

验算可得：
$$
\frac{1}{64}A_0+\frac{1}{8}A_1+\frac{27}{64}A_2=\frac{1}{4}\\
\frac{1}{256}A_0+\frac{1}{16}A_1+\frac{81}{256}A_2\ne\frac{1}{5}
$$
因此代数精度为 3。

(3)
$$
\begin{align*}
\int_0^1f(x)\mathrm{d}x&=\frac{2}{3}f\left(\frac{1}{4}\right)-\frac{1}{3}f\left(\frac{1}{2}\right)+\frac{2}{3}f\left(\frac{3}{4}\right)\\
\int_0^1x^2\mathrm{d}x&=\frac{2}{3}\frac{1}{16}-\frac{1}{3}\frac{1}{4}+\frac{2}{3}\frac{9}{16}\\
&=\frac{1}{3}
\end{align*}
$$

## Question 4

> 设 $f(x)\in C^4[0,2]$，$I(f)=\int_0^2f(x)\mathrm{d}x$。给定求积公式：
> $$
> I(f)\approx\tilde{I}(f)=Af(x_0)+f(x_1)
> $$
> (1) 求常数 $A,x_0,x_1$，使得求积公式的代数精度尽可能提高。
>
> (2) 导出上述求积公式的误差 $I(f)-\tilde{I}(f)$ 的表达式。

(1)

分别带入 $f(x)=1,x,x^2$：
$$
\begin{bmatrix}
1&1\\
x_0&x_1\\
x_0^2&x_1^2\\
\end{bmatrix}
\begin{bmatrix}A\\1\end{bmatrix}=
\begin{bmatrix}
2\\2\\8/3
\end{bmatrix}
$$
可以解得：
$$
A=1,x_0=1-\frac{\sqrt{3}}{3},x_1=1+\frac{\sqrt{3}}{3}
$$
带入 $f(x)=x^3,x^4$：
$$
\left(1-\frac{\sqrt{3}}{3}\right)^3+\left(1+\frac{\sqrt{3}}{3}\right)^3=4=\int_0^2x^3\mathrm{d}x\\
\left(1-\frac{\sqrt{3}}{3}\right)^4+\left(1+\frac{\sqrt{3}}{3}\right)^4=\frac{56}{9}\ne\int_0^2x^4\mathrm{d}x
$$
因此代数精度为 $3$。

求积公式为：
$$
\tilde{I}(f)=f\left(1-\frac{\sqrt{3}}{3}\right)+f\left(1+\frac{\sqrt{3}}{3}\right)
$$
(2)

余项为：
$$
\begin{align*}
R[f]&=\int_0^2f(x)\mathrm{d}x-f\left(1-\frac{\sqrt{3}}{3}\right)-f\left(1+\frac{\sqrt{3}}{3}\right)\\
&=\frac{f^{(4)}(\xi)}{4!}\int_0^2\left(x-1+\frac{\sqrt{3}}{3}\right)^2\left(x-1-\frac{\sqrt{3}}{3}\right)^2\mathrm{d}x\\
&=\frac{f^{(4)}(\xi)}{135}\quad\xi\in(0,2)
\end{align*}
$$

## Question 5

> 设 $f\in C^3[0,2]$，给定如下的数值积分公式：
> $$
> I=\int_0^2f(x)\mathrm{d}x\approx I_n=A_1f(1)+A_2f(2)+A_3f'(1)
> $$
> (1) 求常数 $A_1,A_2,A_3$ 使上述求积公式的代数精度最高。
>
> (2) 导出上述求积公式的余项（或截断误差）$R[f]=I-I_n$。

先算 012 阶：
$$
\begin{bmatrix}1&1&0\\1&2&1\\1&4&2\\\end{bmatrix}
\begin{bmatrix}A_1\\A_2\\A_3\end{bmatrix}=
\begin{bmatrix}2\\2\\\frac{8}{3}\end{bmatrix}
$$
解得：
$$
A_1=\frac{4}{3},A_2=\frac{2}{3},A_3=-\frac{2}{3}
$$
积分公式为：
$$
I_n=\frac{4}{3}f(1)+\frac{2}{3}f(2)-\frac{2}{3}f'(1)
$$
验证可得：
$$
A_1+8A_2+3A_3=14/3\ne4
$$
因此代数精度为 2。

(2)

构造二次多项式 $H(x)$，满足：
$$
H(1)=f(1),H'(1)=f'(1),H(2)=f(2)
$$
其余项满足：
$$
R(x)=f(x)-H(x)=\frac{f^{(3)}(\xi)}{3!}(x-1)^2(x-2)\quad(\xi\in(0,2))
$$
其求积公式余项为：
$$
\begin{align*}
R[f]&=\int_0^2\frac{f^{(3)}(\xi)}{3!}(x-1)^2(x-2)\\
&=-\frac{f^{(3)}(\eta)}{9}\quad\eta\in(0,2)
\end{align*}
$$

## Question 6

> 记 $T_n$ 表示将区间 $[a,b]$ $n$ 等分后，按复化梯形公式计算积分所得的近似值，已知 $T_4=0.94451$, $T_8=0.94568$。
>
> (1) 计算 $I-T_8$ 的近似值。
>
> (2) 计算经过一次 Richardson 外推所得 $I$ 的近似值。

(1)
$$
\begin{align*}
I&\approx\frac{1}{3}(4T_8-T_4)\\
I-T_8&\approx\frac{1}{3}(T_8-T_4)=0.00039
\end{align*}
$$
(2)

$T_4$ 是第 3 次近似，$T_8$ 是第 4 次近似
$$
I\approx\frac{4T_8-T_4}{4-1}=0.94607
$$


## Question 7

> 设 $f(x)\in C^2[a,b]$, $I[f]=\int_a^bf(x)\mathrm{d}x$, $h=(b-a)/n$, $x_i=a+ih$, $i\in\{0,n\}$. $T_n(f)$ 为计算 $I(f)$ 的复化梯形公式，证明：
> $$
> \lim_{h\to0}\frac{I(f)-T_n(f)}{h^2}=\frac{1}{12}[f'(a)-f'(b)]
> $$

已知：

$$
I(f)-T_n(f)=-\sum_{k=0}^{n-1}\frac{h^3}{12}f''(\eta_k)\quad\eta_k\in[x_k,x_{k+1}]
$$

那么：
$$
\begin{align*}
\lim_{h\to0}\frac{I(f)-T_n(f)}{h^2}&=-\lim_{h\to0}\sum_{k=0}^{n-1}\frac{h}{12}f''(\eta_k)\\
&=\lim_{n\to\infty}\sum_{k=0}^{n-1}\frac{(a-b)}{12n}f''(\eta_k)\\
&=\frac{1}{12}\int_a^bf''(x)\mathrm{d}x\\
&=\frac{1}{12}[f'(a)-f'(b)]\\
\end{align*}
$$

## Question 8

> 考虑求积公式 $\int_0^1f(x)\mathrm{d}x\approx af(0)+bf(1)+c[f'(0)-f'(1)]$.
>
> (1) 确定 $a,b,c$ 的值，使得上述求积公式的代数精度尽可能高，并求它的最高的代数精度。
>
> (2) 将积分区间 $[0,1]$ 二等分，用上述公式的复化形式求 $\int_0^1\ln(1+x)\mathrm{d}x$ 的近似值。

(1)

带入 $1,x,x^2$：
$$
\begin{bmatrix}1&1&0\\0&1&0\\0&1&-2\end{bmatrix}\begin{bmatrix}a\\b\\c\end{bmatrix}=\begin{bmatrix}1\\1/2\\1/3\end{bmatrix}
$$
解得：
$$
a=\frac{1}{2},b=\frac{1}{2},c=\frac{1}{12}
$$
求积公式为：
$$
\frac{1}{2}[f(0)+f(1)]+\frac{1}{12}[f'(0)-f'(1)]
$$
带入 $x^3,x^4$：
$$
\int_0^1x^3\mathrm{d}x=\frac{1}{4}=\frac{1}{2}-\frac{3}{12}\\
\int_0^1x^4\mathrm{d}x\ne\frac{1}{5}=\frac{1}{2}-\frac{4}{12}\\
$$
因此代数精度为 3。

(2)

一般区间的求积公式为：
$$
\int_a^bf(x)\mathrm{d}x\approx \frac{b-a}{2}[f(0)+f(1)]+\left[\frac{1}{12}(a^2+b^2)-\frac{1}{6}ab\right][f'(0)-f'(1)]
$$
直接带入相加可得：
$$
I=\frac{1}{4}[f(0)+2f(1/2)+f(1)]+\frac{1}{48}[f'(0)-f'(1)]
$$
则带入值可得：

$$
\begin{align*}
I&=\frac{1}{4}[0+2\ln(3/2)+\ln(2)]+\frac{1}{48}\left(\frac{1}{1}-\frac{1}{2}\right)\\
&=\frac{1}{2}\ln3-\frac{1}{4}\ln2+\frac{1}{96}\\
&=0.386436
\end{align*}
$$

## Question 9

> 根据“割圆法”的思想，当圆的外切正 $n$ 边形的边数 $n\to\infty$ 时，其周长与圆的直径之比趋于圆周率 $\pi$。设圆的直径为 1，则圆的外切正 $n$​ 边形的周长 
> $$
> L(n)=n\tan\frac{\pi}{n}
> $$
> 由 Taylor 展式
> $$
> \tan x=x+\frac{1}{3}x^3+\frac{2}{15}x^5+\frac{17}{315}x^7+\cdots\quad|x|\le\frac{\pi}{2}
> $$
> 记 $h=1/n$，$A(h)=L(n)$。
>
> (1) 验证用 $A(h),A(h/2)$ 逼近 $\pi$ 的误差阶数均为 $O(h^2)$。
>
> (2) 导出计算 $\pi$ 的高阶计算公式。

(1)
$$
\begin{align*}
A(h)=L(n)&=n\left(\frac{\pi}{n}+\frac{1}{3}\left(\frac{\pi}{n}\right)^3+\frac{2}{15}\left(\frac{\pi}{n}\right)^5+\cdots\right)\\
&=\pi+\frac{1}{3}\frac{\pi^3}{n^2}+\frac{2}{15}\frac{\pi^5}{n^4}+\cdots\\
&=\pi+\mathcal{O}(h^2)
\end{align*}
$$

$$
\begin{align*}
A\left(\frac{h}{2}\right)=L(2n)&=2n\left(\frac{\pi}{2n}+\frac{1}{3}\left(\frac{\pi}{2n}\right)^3+\frac{2}{15}\left(\frac{\pi}{2n}\right)^5+\cdots\right)\\
&=\pi+\frac{1}{3}\frac{\pi^3}{(2n)^2}+\frac{2}{15}\frac{\pi^5}{(2n)^4}+\cdots\\
&=\pi+\mathcal{O}(h^2)
\end{align*}
$$

(2)

已知：
$$
A(h)-\pi=\frac{1}{3}\frac{\pi^3}{n^2}+\frac{2}{15}\frac{\pi^5}{n^4}+\cdots\\
A\left(\frac{h}{2}\right)-\pi=\frac{1}{3}\frac{\pi^3}{(2n)^2}+\frac{2}{15}\frac{\pi^5}{(2n)^4}+\cdots\\
$$
下式乘 4 减去上式：
$$
4A\left(\frac{h}{2}\right)-A(h)-3\pi=\mathcal{O}(h^4)\\
\pi\approx\frac{4L(2n)-L(n)}{3}=\frac{1}{3}\left(8n\tan\frac{\pi}{4n}-n\tan\frac{\pi}{n}\right)
$$

## Question 10

> 构造带权的高斯积分公式 $\int_{-1}^{1}(x+1)f(x)\mathrm{d}x\approx A_0f(x_0)+A_1f(x_1)$，可分如下两步完成：
>
> (1) 先构造区间 $[-1,1]$ 上带权函数 $\rho(x)=x+1$ 正交的 $n$ 次多项式 $P_n(x)(n=0,1,2)$，当它们的最高次项系数均为 $1$ 时，求 $P_n(x)(n=0,1,2)$。
>
> (2) 在此基础上求高斯积分公式的求积节点和系数。

(1)

设 $P_2(x)=(x-x_0)(x-x_1)=x^2+bx+c$，且与 $1,x$ 带权函数 $\rho(x)=x+1$ 正交。因此：
$$
\int_{-1}^{1}(x+1)P_2(x)\mathrm{d}x=\int_{-1}^{1}(x+1)(x^2+bx+c)\mathrm{d}x=\frac{2}{3}(b+1)+2c=0\\
\int_{-1}^{1}x(x+1)P_2(x)\mathrm{d}x=\int_{-1}^{1}x(x+1)(x^2+bx+c)\mathrm{d}x=\frac{2}{5}+\frac{2}{3}(b+c)=0\\
$$
可以得到：
$$
b=-\frac{2}{5},c=-\frac{1}{5}
$$
(2)

先求 $P_2(x)$ 的零点：
$$
x_0,x_1=\frac{\frac{2}{5}\pm\sqrt{\frac{24}{25}}}{2}=\frac{1\pm\sqrt{6}}{5}
$$
带入 $1,x$ 到公式中：
$$
\int_{-1}^{1}(x+1)\mathrm{d}x=2=A_0+A_1\\
\int_{-1}^{1}(x+1)x\mathrm{d}x=\frac{2}{3}=A_0\frac{1-\sqrt{6}}{5}+A_1\frac{1+\sqrt{6}}{5}\\
$$
得到：
$$
A_0=\frac{9-\sqrt{6}}{9},A_1=\frac{9+\sqrt{6}}{9}
$$
因此求积公式为：
$$
\int_{-1}^{1}(x+1)f(x)\mathrm{d}x\approx \frac{9-\sqrt{6}}{9}f\left(\frac{1-\sqrt{6}}{5}\right)+\frac{9+\sqrt{6}}{9}f\left(\frac{1+\sqrt{6}}{5}\right)
$$

## Question 11

> 设 $f(x)\in C^2[a,b]$，$I(f)=\int_a^bf(x)\mathrm{d}x$。记 $h=(b-a)/n$，$x_i=a+ih$，$i\in\{0,\dots,n\}$，$T_n(f)$ 为计算 $I(f)$ 的复化梯形公式。
>
> (1) 导出计算 $I(f)$ 的一点的 Gauss 型求积公式 $G_0(f)\equiv A_0f(x_0)$ 及其复化公式 $G_n(f)$。
>
> (2) 求参数 $\alpha$，使得 $T_{2n}(f)=\frac{1}{2}T_n(f)+\alpha G_n(f)$。

(1)

带入 $1,x$ 到公式中：
$$
\int_a^b1\mathrm{d}x=(b-a)=A_0\\
\int_a^bx\mathrm{d}x=(b+a)/2=A_0x_0\\
$$
得到：
$$
A_0=b-a,x_0=\frac{b+a}{2}
$$
因此公式为：
$$
G(x)=(b-a)f\left(\frac{b+a}{2}\right)
$$
复化可得：
$$
\begin{align*}
G_n(f)&=\sum_{k=0}^{n-1}\int_{x_k}^{x_{k+1}}f(x)\mathrm{d}x\\
&\approx\sum_{k=0}^{n-1}(x_{k+1}-x_k)f\left(\frac{x_{k+1}+x_k}{2}\right)\\
&=h\sum_{k=0}^{n-1}f(x_{k+\frac{1}{2}})\\
\end{align*}
$$
其中 $h=(b-a)/n$，$x_{k+\frac{1}{2}}=(x_{k+1}+x_k)/2$​。

(2)

带入：
$$
T_{2n}(f)=\frac{1}{2}T_n(f)+\alpha h\sum_{k=0}^{n-1}f(x_{k+\frac{1}{2}})
$$
由 Romberg 公式可得：$\alpha=\frac{1}{2}$。

