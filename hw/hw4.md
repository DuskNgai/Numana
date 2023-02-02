## Question 1

> 当 $x=1,-1,2$ 时，$f(x)=0,-3,4$，求 $f(x)$ 的二次插值多项式。

$$
\begin{align*}
L_2(x)&=0-3\frac{(x-1)(x-2)}{(-1-1)(-1-2)}+4\frac{(x-1)(x+1)}{(2-1)(2+1)}\\
&=-\frac{(x-1)(x-2)}{2}+4\frac{(x-1)(x+1)}{3}\\
&=\frac{5}{6}x^2+\frac{3}{2}x-\frac{7}{3}
\end{align*}
$$

## Question 2

> 取节点 $x_0=0,x_1=0.5,x_2=1$，估计函数 $f(x)=e^{-x}$ 在区间 $[0,1]$ 上的二次插值多项式的误差。

在区间 $[0,1]$ 上存在 $\xi$ 使得
$$
\begin{align*}
R_2(x)&=L_2(x)-f(x)\\
&=\frac{f^{3}(\xi)}{3!}x(x-0.5)(x-1)\\
&=\frac{-e^{-\xi}}{6}x(x-0.5)(x-1)
\end{align*}
$$

$$
\omega(x)=x(x-0.5)(x-1)\\
\omega'(x)=(x-0.5)(x-1)+x(x-1)+x(x-0.5)=3x^2-3x+0.5\\
\omega'(x)=0\Longrightarrow x=\frac{3\pm\sqrt{3}}{6}
$$

对应的 $\omega(x)$ 为 $\pm\sqrt{3}/36$。因此：
$$
-\frac{\sqrt{3}}{216}e^{-\xi}<R_2(x)<\frac{\sqrt{3}}{216}e^{-\xi}
$$

## Question 3

> 根据下列数表，求 Newton 插值多项式。
>
> |  $x_i$   |  1   | 1.5  |  2   | 2.5  |  3   |  3.5  |
> | :------: | :--: | :--: | :--: | :--: | :--: | :--: |
> | $f(x_i)$ | -1 | 0.5 | 2.5 | 5 | 8 | 11.5 |
>

$$
\begin{matrix}
1  &\vline&-1.0&\\
   &\vline&    &3.0\\
1.5&\vline&0.5 &   &1.0\\
   &\vline&    &4.0&   &0.0\\
2  &\vline&2.5 &   &1.0&\\
   &\vline&    &5.0&   &0.0\\
2.5&\vline&5.0 &   &1.0&\\
   &\vline&    &6.0&   &0.0\\
3  &\vline&8.0 &   &1.0\\
   &\vline&    &7.0\\
3.5&\vline&11.5\\
\end{matrix}
$$

因此这个多项式为：
$$
N_5(x)=-1+3(x-1)+(x-1)(x-1.5)
$$

## Question 4

> 若 $f(x)=\sum_{i=0}^na_ix^i$ 有 $n$ 个不同实根 $x_1,\dots,x_n$，证明：
> $$
> \sum_{j=1}^n\frac{x_j^k}{f'(x_j)}=\begin{cases}0&0\le k<n-2\\a_n^{-1}&k=n-1\end{cases}
> $$

由于 $f(x)$ 具有 $n$ 个不同实根，因此 $f(x)$ 可被因式分解为：
$$
f(x)=a_n\prod_{i=1}^n(x-x_i)
$$
求导可得：
$$
f'(x)=a_n\sum_{j=1}^n\prod_{i=1,i\ne j}^n(x-x_i)
$$
带入 $x_j\in\{x_1,\dots,x_n\}$ 可得：
$$
f'(x_j)=a_n\prod_{i=1,i\ne j}^n(x_j-x_i)\ne0
$$
因此：
$$
\sum_{j=1}^n\frac{x_j^k}{f'(x_j)}=\sum_{j=1}^n\frac{x_j^k}{a_n\prod_{i=1,i\ne j}^n(x_j-x_i)}=\frac{1}{a_n}\sum_{j=1}^n\frac{x_j^k}{\prod_{i=1,i\ne j}^n(x_j-x_i)}
$$
设被插值函数为 $g(x)=x^k,k<n-1$，插值点为 $x_1,\dots,x_n$，则对应的 Lagrange 多项式为：
$$
P_{n-1}(x)=\sum_{i=1}^nx_i^k\left(\prod_{j=1,j\ne i}^{n}\frac{x-x_j}{x_i-x_j}\right)
$$
当 $k<n-1$ 时，$P_{n-1}(x)$ 的最高项 ($n-1$ 次项) 的系数必然为 0，即：
$$
\sum_{i=1}^n\left(\frac{x_i^k}{\prod_{j=1,j\ne i}^{n}(x_i-x_j)}\right)=0
$$
当 $k=n-1$​ 时，$P_{n-1}(x)$​ 的最高项 ($n-1$​ 次项) 的系数必然为 1。因此：
$$
\sum_{j=1}^n\frac{x_j^k}{f'(x_j)}=\begin{cases}0&0\le k<n-2\\a_n^{-1}&k=n-1\end{cases}
$$

## Question 5

> 已知
> $$
> f\left(-\frac{1}{2}\right)=f\left(\frac{1}{2}\right)=-\frac{3}{4},f(2)=3,f'\left(\frac{1}{2}\right)=1
> $$
> (1) 求满足插值条件的插值多项式 $P(x)$。
>
> (2) 若 $f(x)$ 不是二次函数，且在 $[-1/2,2]$ 内满足 $|f^{(4)}(x)|<1$，证明 $|f(1)|<1/64$。

(1)
$$
\begin{matrix}
-0.5&\vline&-0.75&\\
    &\vline&     &0  \\
0.5 &\vline&-0.75&   &1\\
    &\vline&     &1  & &0\\
0.5 &\vline&-0.75&   &1\\
    &\vline&     &2.5\\
2   &\vline&3\\
\end{matrix}
$$
因此多项式为：
$$
P(x)=-\frac{3}{4}+\left(x+\frac{1}{2}\right)\left(x-\frac{1}{2}\right)
$$
(2)

余项为：
$$
R(x)=f(x)-P(x)=\frac{f^{(4)}(\xi)}{4!}\left(x+\frac{1}{2}\right)\left(x-\frac{1}{2}\right)^2(x-2)
$$
其中 $\xi\in[-0.5,2]$​。
$$
R(1)=f(1)-P(1)=\frac{f^{(4)}(\xi)}{4!}\left(\frac{3}{2}\right)\left(\frac{1}{2}\right)^2(-1)=-\frac{f^{(4)}(\xi)}{64}\\
P(1)=-\frac{3}{4}+\left(\frac{3}{2}\right)\left(\frac{1}{2}\right)=0
$$
因此：
$$
|f(1)|<\left|\frac{f^{(4)}(\xi)}{64}\right|<\frac{1}{64}
$$

## Question 6

> 设 $f\in C^4[-1,1]$，求一个次数不超过 3 次的插值多项式 $P_3(x)$，满足插值条件：
>
> |   $x$   |  -1  |  0   |  1   |
> | :-----: | :--: | :--: | :--: |
> | $f(x)$  |  10  |  0   |  10  |
> | $f'(x)$ |      |  0   |      |
>
> 并写出插值余项 $f(x)-P_3(x)$。
>
> (2) 设 $|f^{(k)}(x)|\le k(k=2,3,4),x\in[-1,1]$，若用 $P_3(0.75)$ 作为 $f(0.75)$ 的近似值，试分析 $P_3(0.75)$ 具有几位有效数字。

(1)
$$
\begin{matrix}
-1&\vline&10&\\
  &\vline&  &-10\\
0 &\vline&0 &   &10\\
  &\vline&  &0  & &0\\
0 &\vline&0 &   &10\\
  &\vline&  &10 \\
1 &\vline&10\\
\end{matrix}
$$
插值多项式：
$$
P_3(x)=10-10(x+1)+10(x+1)x
$$
余项：
$$
R_3(x)=f(x)-P_3(x)=\frac{f^{(4)}(\xi)}{4!}(x+1)x^2(x-1)\quad\xi\in(-1,1)
$$
(2)
$$
P_3(0.75)=5.625\\
R_3(0.75)=f(0.75)-P_3(0.75)=-\frac{f^{(4)}(\xi)}{4!}=-\frac{21f^{(4)}(\xi)}{2048}\\
|R_3(0.75)|=\left|\frac{f^{(4)}(\xi)}{4!}\right|\le\frac{21}{512}=0.041<\frac{1}{2}\times10^{-1}
$$
因此具有 2 位有效数字。

## Question 7

> 设 $x_1,\dots,x_n$ 为插值节点，$P_{n-1}(x)$ 为满足 $P_{n-1}(x)=f(x_i)(i\in\{1,\dots,n\})$ 的插值多项式。
>
> (1) 证明：$P_{n-1}(x)$ 的最高次项的系数为：
> $$
> \sum_{i=1}^n\frac{f(x_i)}{\prod_{j=1,j\ne i}^{n}(x_i-x_j)}
> $$
> (2) 证明：对 $k=0,\dots,n-2$ 有：
> $$
> \sum_{i=1}^n\frac{i^k}{\prod_{j=1,j\ne i}^{n}(i-j)}=0
> $$

(1)

由 Lagrange 多项式：
$$
P_{n-1}(x)=\sum_{i=1}^nf(x_i)\left(\prod_{j=1,j\ne i}^{n}\frac{x-x_j}{x_i-x_j}\right)
$$
可直接得到。

---

由 Newton 多项式：
$$
P_{n-1}(x)=\sum_{i=1}^{n}f[x_1,\dots,x_i]\prod_{j=1}^{i}(x-x_j)
$$
此时
$$
f[x_1,\dots,x_n]=\sum_{i=1}^n\frac{f(x_i)}{\prod_{j=1,j\ne i}^{n}(x_i-x_j)}
$$
(2)

设被插值函数为 $g(x)=x^k,k<n-1$，插值点为 $1,\dots,n$，则对应的 Lagrange 多项式为：
$$
P_{n-1}(x)=\sum_{i=1}^ni^k\left(\prod_{j=1,j\ne i}^{n}\frac{x-j}{i-j}\right)
$$
当 $k<n-1$ 时，$P_{n-1}(x)$ 的最高项 ($n-1$ 次项) 的系数必然为 0。

## Question 8

> 设 $H_4(x)$ 是满足下列插值条件的 4 次多项式：
> $$
> H_4(0)=f(0),H'_4(0)=f'(0)\\
> H_4(1)=f(1),H'_4(1)=f'(1),H''_4(1)=f''(1)\\
> $$
> 若 $f\in C^5[0,1]$，试证明：
> $$
> f(x)=H_4(x)+\frac{1}{5!}f^{(5)}(\xi)x^2(x-1)^3
> $$
> $\xi$ 介于 0, 1 与 $x$ 之间。

设
$$
R_4(x)=f(x)-H_4(x)
$$
由于：
$$
R_4(0)=H_4(0)-f(0)=0,R'_4(0)=H'_4(0)-f'(0)=0\\
R_4(1)=H_4(1)-f(1)=0,R'_4(x)=H'_4(1)-f'(1)=0\\
R''_4(x)=H''_4(1)-f''(1)=0
$$
则 $R_4(x)$ 具有如下形式：
$$
R_4(x)=A(x)x^2(x-1)^3\quad x\in[0,1]
$$
对于某个固定的，非 0, 1 的 $x$，设：
$$
\phi(t)=f(t)-H_4(t)-A(x)t^2(t-1)^3
$$
$\phi(t)$ 有 $0,x,1$ 三个零点。$\phi'(t)$ 有 $0,1$ 和介于 $a\in(0,x),b\in(x,1)$ 之间的 $4$ 个零点，$\phi''(t)$ 有 $c\in(0,a),d\in(a,b),e\in(b,1)$ 和 1 共 4 个零点。反复运用 Rolle 定理，$\phi^{(4)}(t)$ 就有 2 个零点，介于 $0,x,1$ 之间。
$$
\phi^{(5)}(\xi)=f^{(5)}(\xi)-H_4^{(5)}(\xi)-5!A(x)=f^{(5)}(\xi)-5!A(x)=0
$$
因此：
$$
A(x)=\frac{f^{(5)}(\xi)}{5!}
$$
得到：
$$
f(x)=H_4(x)+\frac{1}{5!}f^{(5)}(\xi)x^2(x-1)^3
$$

## Question 9

> (1) 根据所学知识，试给出四次样条函数的定义
>
> (2) 设函数
> $$
> S(x)=\begin{cases}
> x^4+2x+1&0\le x<1\\
> (x-1)^4+a(x-1)^3+b(x-1)^2+c(x-1)+d&1\le x\le3
> \end{cases}
> $$
> 是节点 $x_0=0,x_1=1,x_2=3$ 上的四次样条函数，求常数 $a,b,c,d$ 的值。

(1)

设区间 $[a,b]$​ 上给定一个剖分:
$$
\Delta:a\le x_0<x_1<\dots<x_n\le b
$$
若函数 $S(x)$ 在区间 $[a,b]$ 上满足

1. $S(x)\in C^3[a,b]$
2. $S(x)$ 在每个小区间 $[x_J,x_{j+1}]$ 上是 4 次多项式。

则称 $S$ 是关于剖分 $\Delta$ 的一个 4 次样条曲线。

(2)

0 阶导数信息
$$
S(1)=4\Longrightarrow d=4
$$
1 阶导数信息
$$
S'(x)=4x^3+2\\
S'(1)=6\Longrightarrow c=6
$$
2 阶导数信息
$$
S''(x)=12x^2\\
S''(1)=12=2b\Longrightarrow b=6
$$
3 阶导数信息
$$
S'''(x)=24x\\
S'''(1)=24=6a\Longrightarrow a=4
$$

## Question 10

> 设 $f(x)=x^{1/3}, x\in[-1,1]$ 记 $\Phi_3=\text{span}\{x,x^3\}$。求 $f$ 在 $\Phi_3$ 上的最佳平方逼近 $P_3(x)$。

$$
\langle\phi_1,\phi_1\rangle=\int_{-1}^1x^2\mathrm{d}x=\frac{2}{3}\\
\langle\phi_2,\phi_2\rangle=\int_{-1}^1x^6\mathrm{d}x=\frac{2}{7}\\
\langle\phi_1,\phi_2\rangle=\int_{-1}^1x^4\mathrm{d}x=\frac{2}{5}\\
\langle\phi_1,f\rangle=\int_{-1}^1x^{4/3}\mathrm{d}x=\frac{6}{7}\\
\langle\phi_2,f\rangle=\int_{-1}^1x^{10/3}\mathrm{d}x=\frac{6}{13}\\
$$

因此矩阵为：
$$
\begin{bmatrix}
\dfrac{2}{3}&\dfrac{2}{5}\\
\dfrac{2}{5}&\dfrac{2}{7}
\end{bmatrix}
\begin{bmatrix}a\\b\end{bmatrix}=
\begin{bmatrix}\dfrac{6}{7}\\\dfrac{6}{13}\end{bmatrix}
$$
得到
$$
a=\frac{180}{91},b=-\frac{15}{13}\\
f(x)=\frac{180}{91}x-\frac{15}{13}x^3
$$

## Question 11

> 以$-1, 0, 1$ 为节点的一次样条函数类 $S_1$​ 定义如下：
> $$
> S_1=\{S(x)|S(x)\in C[-1,1]\}
> $$
> 设 $f(x)=e^x$。
>
> (1) 求 $S(x)\in S_1$ 满足插值条件。
>
> (2) 求 $f$ 在 $S_1$ 中的最佳平方逼近。

(1)
$$
f(-1)=e^{-1},f(0)=1,f(1)=e
$$
因此
$$
S(x)=\begin{cases}
-(x-0)\times e^{-1}+(x+1)\times1&x\in[-1,0]\\
-(x-1)\times1+(x-0)\times e&x\in[0,1]\\
\end{cases}\\
=\begin{cases}
(1-e^{-1})x+1&x\in[-1,0]\\
(e-1)x+1&x\in[0,1]\\
\end{cases}
$$
(2)

设
$$
S(x)=\begin{cases}
mx+b&x\in[-1,0]\\
nx+b&x\in[0,1]
\end{cases}
$$
则
$$
\int_{-1}^{1}f^2(x)\mathrm{d}x=\int_{-1}^{1}e^{2x}\mathrm{d}x=\frac{1}{2}(e^2-e^{-2})
$$

$$
\begin{align*}
\int_{-1}^{1}f(x)S(x)\mathrm{d}x&=\int_{-1}^{0}e^x(mx+b)\mathrm{d}x+\int_{0}^{1}e^x(nx+b)\mathrm{d}x\\
&=m(x-1)e^x\mid_{-1}^0+be^x\mid_{-1}^0+n(x-1)e^x\mid_0^1+be^x\mid_0^1\\
&=m(-1+2e^{-1})+n+b(e-e^{-1})\\
\end{align*}
$$

$$
\begin{align*}
\int_{-1}^{1}S^2(x)\mathrm{d}x&=\int_{-1}^{0}(mx+b)^2\mathrm{d}x+\int_{0}^{1}(nx+b)^2\mathrm{d}x\\
&=\frac{1}{3}m^2-mb+b^2+\frac{1}{3}n^2+nb+b^2\\
&=\frac{1}{3}m^2-mb+\frac{1}{3}n^2+nb+2b^2\\
\end{align*}
$$

则
$$
\begin{align*}
I(m,n,b)&=\|f(x)-S(x)\|^2_2\\
&=\int_{-1}^{1}f^2(x)-2f(x)S(x)+S^2(x)\mathrm{d}x\\
&=C-2m(-1+2e^{-1})-2n-2b(e-e^{-1})+\frac{1}{3}m^2-mb+\frac{1}{3}n^2+nb+2b^2
\end{align*}
$$

$$
\frac{\partial I}{\partial m}=2-4e^{-1}+\frac{2}{3}m-b=0\\
\frac{\partial I}{\partial n}=-2+\frac{2}{3}n+b=0\\
\frac{\partial I}{\partial b}=-2(e-e^{-1})-m+n+4b=0\\
$$

因此
$$
\begin{cases}
m=3e+12e^{-1}-12\\
n=-3e-6e^{-1}+12\\
b=2e+4e^{-1}-6\\
\end{cases}
$$
则：
$$
S(x)=\begin{cases}
(3e+12e^{-1}-12)x+(2e+4e^{-1}-6)&x\in[-1,0]\\
(-3e-6e^{-1}+12)x+(2e+4e^{-1}-6)&x\in[0,1]
\end{cases}
$$

## Question 12

> 给定函数 $y=f(x)$ 在节点 $x_i(i=1\sim5)$ 处的值如下：
>
> | $x_i$ |  -2  |  -1  |  0   |  1   |  2   |
> | :---: | :--: | :--: | :--: | :--: | :--: |
> | $y_i$ | 0.5  | 0.5  | 1.5  | 3.5  | 6.5  |
>
> 用最小二乘法求上述数据的直线拟合。

基函数为：$\phi_0(x)=1,\phi_1(x)=x$。

因此：
$$
\langle\phi_0(x),\phi_0(x)\rangle=\sum_{i=1}^{5}1=5\\
\langle\phi_0(x),\phi_1(x)\rangle=\sum_{i=1}^{5}x_i=0\\
\langle\phi_1(x),\phi_1(x)\rangle=\sum_{i=1}^{5}x_i^2=10\\
\langle\phi_0(x),f(x)\rangle=\sum_{i=1}^{5}y_i=12.5\\
\langle\phi_1(x),f(x)\rangle=\sum_{i=1}^{5}x_iy_i=15\\
$$
法方程为：
$$
\begin{bmatrix}5&0\\0&10\end{bmatrix}\begin{bmatrix}a\\b\end{bmatrix}=\begin{bmatrix}12.5\\15\end{bmatrix}
$$
因此 $a=2.5,b=1.5$，方程为：
$$
y=2.5+1.5x
$$

## Question 13

> 设
> $$
> A=\begin{bmatrix}1&2\\0&-1\\-1&0\end{bmatrix},b=\begin{bmatrix}-1\\3\\2\end{bmatrix}
> $$
> 求最小二乘解。

$$
\begin{align*}
x&=(A^TA)^{-1}A^Tb\\
&=\begin{bmatrix}2&2\\2&5\end{bmatrix}^{-1}\begin{bmatrix}1&0&-1\\2&-1&0\end{bmatrix}\begin{bmatrix}-1\\3\\2\end{bmatrix}\\
&=\frac{1}{6}\begin{bmatrix}5&-2\\-2&2\end{bmatrix}\begin{bmatrix}1&0&-1\\2&-1&0\end{bmatrix}\begin{bmatrix}-1\\3\\2\end{bmatrix}\\
&=\begin{bmatrix}-5/6\\-2/3\end{bmatrix}
\end{align*}
$$

