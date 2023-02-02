# Chapter 3 插值 Interpolation

## 3.4 三次样条 Cubic Spline

样条使用多个公式，每个都是低阶多项式，来连接所有数据点。

一次样条（直线连接）连续但不一阶连续。

二次样条（抛物线线连接）一阶连续但不二阶连续，即斜率光滑。

三次样条二阶连续，即曲率光滑，拟合较好。

### 3.4.1 分段线性插值

给定节点 $a\le x_0<\dots<x_n\le b$ 的函数值 $y_j$，记 $h_j=x_{j+1}-x_j,h=\max_j{h_j}$。

如果函数 $I_h(x)$ 满足在 $[a,b]$ 上连续，插值每个函数点，且在每个小区间上为线性函数，则称 $I(x)$ 为 $f$ 的分段线性插值函数。$I(x)$ 在小区间 $[x_j,x_{j+1}]$​ 可用 Lagrange 插值表示为：
$$
I_h(x)=\frac{x-x_{j+1}}{x_j-x_{j+1}}y_j+\frac{x-x_j}{x_{j+1}-x_j}y_{j+1}\quad x\in[x_j,x_{j+1}]
$$
利用基函数可写为 $I_h(x)=\sum_{j=0}^{n}y_jI_j(x)$，其中：
$$\begin{gather*}
l_0(x)=\frac{x-x_1}{x_0-x_1}\quad x\in[x_0,x_{1}]\\
l_j(x)=\begin{cases}\dfrac{x-x_{j-1}}{x_j-x_{j-1}}\quad x\in[x_{j-1},x_{j}]\\
\dfrac{x-x_{j+1}}{x_j-x_{j+1}}\quad x\in[x_{j},x_{j+1}]\\
\end{cases}\\
l_n(x)=\frac{x-x_{n-1}}{x_n-x_{n-1}}\quad x\in[x_{n-1},x_{n}]
\end{gather*}
$$
分段线性插值有如下的误差估计：
$$
\|R\|_{\infty}=\|f-I_h\|_{\infty}\le\frac{1}{8}\|f''\|_{\infty}h^2
$$

### 3.4.2 分段三次 Hermite 插值

给定节点 $a\le x_0<\dots<x_n\le b$ 的函数值和导数值 $y_j$ 和 $m_j$。

记 $h_i=x_{i+1}-x_i$，$h=\max_j\{h_j\}$ 如果函数 $H_h(x)$ 满足：

1. $H_h(x)\in C^1[a,b]$
2. 满足 Hermite 插值条件
3. $H_h(x)$ 在每个小区间 $[x_j,x_{j+1}]$ 上是三次多项式

利用基函数 $H_h(x)$ 可写为 $H_h(x)=\sum_{j=0}^n[y_j\alpha_j(x)+m_j\beta_j(x)]$。

分段三次 Hermite 插值有如下的余项估计：
$$
\|R\|_{\infty}=\|f-H_h\|_{\infty}\le\frac{35}{27}\|f'\|_{\infty}h
$$

### 3.4.3 样条的性质

给定节点 $a\le x_0<\dots<x_n\le b$，若函数 $S(x)$ 在区间 $[a,b]$ 上满足：

1. $S(x)\in C^{k-1}[a,b]$
2. $S(x)$ 在每个小区间上是 $k$ 次多项式

则称 $S$ 关于区间的一个 $k$ 次样条函数。若插值节点的函数值，则为 $k$ 次样条插值函数。

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义  三次样条
</blockquote>

$$
\begin{align*}
S_0(x)&=y_0+b_0(x-x_0)+c_0(x-x_0)^2+d_0(x-x_0)^3&x\in[x_0,x_1]\\
S_1(x)&=y_1+b_1(x-x_1)+c_1(x-x_1)^2+d_1(x-x_1)^3&x\in[x_1,x_2]\\
&\dots\\
S_{n-1}(x)&=y_{n-1}+b_{n-1}(x-x_{n-1})+c_{n-1}(x-x_{n-1})^2\\
&+d_{n-1}(x-x_{n-1})^3&x\in[x_{n-1},x_n]\\
\end{align*}
$$

且满足 0,1,2 阶导数值相等
$$\begin{gather*}
S_i(x_i)=y_i\\
S_i(x_{i+1})=y_{i+1},i\in[0, n-1]\\
S'_{i-1}(x_i)=S'_i(x_i),i\in[1, n-1]\\
S''_{i-1}(x_i)=S''_i(x_i),i\in[1, n-1]
\end{gather*}
$$

这样总共 $3n$ 个未知数，$3n-2$ 个方程，满足上述条件的三次样条有无数多个。

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    条件 a: 自然三次样条
</blockquote>

因此如果要求两个端点处二阶导为 0，则这样的三次样条被称为*自然三次样条*。

这些方程分别为
$$
y_{i+1}=y_i+b_i(x_{i+1}-x_i)+c_i(x_{i+1}-x_i)^2+d_i(x_{i+1}-x_i)^3,i\in[0, n-1]
$$

$$
\begin{align*}
&\quad S'_{i-1}(x_i)-S'_i(x_i)\\
&=b_{i-1}-b_i+2c_{i-1}(x_i-x_{i-1})+3d_{i-1}(x_i-x_{i-1})^2\\
&=0,i\in[1, n-1]
\end{align*}
$$

$$
\begin{align*}
&\quad S''_{i-1}(x_i)-S''_i(x_i)\\
&=2c_{i-1}-2c_i+6d_{i-1}(x_i-x_{i-1})\\
&=0,i\in[1, n-1]
\end{align*}
$$

$$
S''_{1}(x_1)=2c_1=0\\
S''_{n-1}(x_n)=2c_{n-1}+6d_{n-1}(x_n-x_{n-1})=0\\
$$

记 $\Delta x_i=x_{i+1}-x_i,\Delta y_i=y_{i+1}-y_i$，引入 $c_n=S''_{n-1}(x_n)/2$，则：

带入二阶导：
$$
d_i=\frac{c_{i+1}-c_i}{3\Delta x_i},i\in[1, n-1]
$$
带入零阶导：
$$
\begin{align*}
b_i&=\frac{\Delta y_i}{\Delta x_i}-c_i\Delta x_i-d_i(\Delta x_i)^2\\
&=\frac{\Delta y_i}{\Delta x_i}-c_i\Delta x_i-\frac{(c_{i+1}-c_i)\Delta x_i}{3}\\
&=\frac{\Delta y_i}{\Delta x_i}-\frac{\Delta x_i}{3}(c_{i+1}+2c_{i}),i\in[1, n-1]\\
\end{align*}
$$
带入一阶导：
$$\begin{gather*}
\frac{\Delta y_{i}}{\Delta x_{i}}-\frac{\Delta x_{i}}{3}(c_{i+1}+2c_{i})-\frac{\Delta y_{i+1}}{\Delta x_{i+1}}+\frac{\Delta x_{i+1}}{3}(c_{i+2}+2c_{i+1})+2c_i\Delta x_{i}+(c_{i+1}-c_i)(\Delta x_{i})^2=0\\
\Delta x_{i+1}(c_{i+2}+2c_{i+1})-\Delta x_{i}(c_{i+1}+2c_{i})+6c_i\Delta x_{i}+3(c_{i+1}-c_i)\Delta x_{i}=3(\frac{\Delta y_{i+1}}{\Delta x_{i+1}}-\frac{\Delta y_{i}}{\Delta x_{i}})\\
\Delta x_{i+1}(c_{i+2}+2c_{i+1})+c_i\Delta x_{i}+2c_{i+1}\Delta x_{i}=3(\frac{\Delta y_{i+1}}{\Delta x_{i+1}}-\frac{\Delta y_{i}}{\Delta x_{i}})\\
c_i\Delta x_{i}+2c_{i+1}(\Delta x_{i}+\Delta x_{i+1})+c_{i+2}\Delta x_{i+1}=3(\frac{\Delta y_{i+1}}{\Delta x_{i+1}}-\frac{\Delta y_{i}}{\Delta x_{i}}),i\in[1, n-2]
\end{gather*}
$$
写成矩阵形式就是：
$$
\underbrace{\begin{bmatrix}
1&0&0& & &\\
\Delta x_1&2(\Delta x_1+\Delta x_2)&\Delta x_2&0&\ddots&\\
0&\Delta x_2&2(\Delta x_2+\Delta x_3)&\Delta x_3&\ddots&\\
 &\ddots&\ddots&\ddots&\ddots&\\
 & & &\Delta x_{n-2}&2(\Delta x_{n-2}+\Delta x_{n-1})&\Delta x_{n-1}\\
 & & &0&0&1
\end{bmatrix}}_{A}
\underbrace{\begin{bmatrix}
c_1\\
\vdots\\
c_n
\end{bmatrix}}_{\mathbf c}=
\underbrace{\begin{bmatrix}
0\\
3(\frac{\Delta y_2}{\Delta x_2}-\frac{\Delta y_1}{\Delta x_1})\\
\vdots\\
3(\frac{\Delta y_{n-1}}{\Delta x_{n-1}}-\frac{\Delta y_{n-2}}{\Delta x_{n-2}})\\
0
\end{bmatrix}}_{\mathbf b}
$$

$A$ 显然是可逆矩阵，故一定可以解出 $\mathbf c$。

### 3.4.4 端点条件

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    条件 b: 曲率调整三次样条
</blockquote>

端点处的二阶导数值是输入值，即：
$$
2c_1=u\\
2c_n=v
$$
对应 $A,\mathbf b$ 的首尾两行需修改为
$$
\begin{bmatrix}
2&0&0&\cdots&0\\
0&\cdots&0&0&2
\end{bmatrix}\mathbf c=
\begin{bmatrix}
u\\
v
\end{bmatrix}
$$

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    条件 c: 嵌制三次样条
</blockquote>

端点处的一阶导数值是输入值，即：
$$\begin{gather*}
S_1'(x)=b_1=u\\
2\Delta x_1c_{1}+\Delta x_1c_{2}=3(\frac{\Delta y_1}{\Delta x_1}-u)\\
\frac{\Delta y_{n-1}}{\Delta x_{n-1}}-\frac{\Delta x_{n-1}}{3}(c_{n}+2c_{n-1})+2c_{n-1}\Delta x_{n-1}+3\frac{c_{n}-c_{n-1}}{3\Delta x_{n-1}}(\Delta x_{n-1})^2=v\\
2c_{n}\Delta x_{n-1}+c_{n-1}\Delta x_{n-1}=3(v-\frac{\Delta y_{n-1}}{\Delta x_{n-1}})
\end{gather*}
$$
对应 $A,\mathbf b$ 的首尾两行需修改为
$$
\begin{bmatrix}
2\Delta x_{1}&\Delta x_{1}&0&\cdots&0\\
0&\cdots&0&\Delta x_{n-1}&2\Delta x_{n-1}
\end{bmatrix}\mathbf c=
\begin{bmatrix}
3(\frac{\Delta y_1}{\Delta x_1}-u)\\
3(v-\frac{\Delta y_{n-1}}{\Delta x_{n-1}})
\end{bmatrix}
$$

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    条件 d: 抛物线端点的三次样条
</blockquote>

端点所连的两条曲线 $S_1(x),S_{n-1}(x)$ 为二阶曲线，即 $d_1=d_{n-1}=0$，即 $c_1=c_2, c_{n-1}=c_n$，这里需假设 $n\ge3$。

对应 $A,\mathbf b$ 的首尾两行需修改为
$$
\begin{bmatrix}
1&-1&0&\cdots&0\\
0&\cdots&0&1&-1
\end{bmatrix}\mathbf c=
\begin{bmatrix}
0\\
0
\end{bmatrix}
$$

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    条件 e: 非纽结的三次样条
</blockquote>

要求 $S_1'''(x)=S_2'''(x),S_{n-2}'''(x)=S_{n-1}'''(x)$，即 $d_1=d_2,d_{n-2}=d_{n-1}$。由于 $S_1(x)$ 和 $S_2(x)$ 均为 3 阶或更低的多项式且对应的 0, 1, 2, 3 阶都相同，因此二者其实是同一个表达式。同理 $S_{n-2}(x)$ 和 $S_{n-1}(x)$ 也是同一个表达式。
$$
d_1=d_2\\
\frac{c_2-c_1}{3\Delta x_1}=\frac{c_3-c_2}{3\Delta x_2}\\
\Delta x_2c_1-(\Delta x_1+\Delta x_2)c_2+\Delta x_1c_3=0\\
$$
对应 $A,\mathbf b$​ 的首尾两行需修改为
$$
\begin{bmatrix}
\Delta x_2&-(\Delta x_1+\Delta x_2)&\Delta x_1&\cdots&0\\
0&\cdots&\Delta x_{n-1}&-(\Delta x_{n-1}+\Delta x_{n-2})&\Delta x_{n-2}
\end{bmatrix}\mathbf c=
\begin{bmatrix}
0\\
0
\end{bmatrix}
$$

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定理 3.1
</blockquote>

条件 a, b, c 对于 $n\ge2$ 有唯一解，条件 d 对于 $n\ge3$ 有唯一解，条件 e 对于$n\ge4$  有唯一解。

### 3.4.5 函数的最佳平方逼近

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义
</blockquote>

设 $f\in C[a,b]$，$\Phi=\mathrm{span}\{\phi_0,\dots,\phi_n\}\subset C[a,b]$，若存在函数 $\phi^*\in\Phi$ 使得
$$
\|f-\phi^*\|_2=\inf_{\phi\in\Phi}\|f-\phi\|_2
$$
则 $\phi^*$ 为 $f$ 在 $\Phi$ 中的最佳平方逼近。即求：
$$
\begin{align*}
I(a_0,\dots,a_n)&=\int_a^bw(x)\left(\phi^*(x)-f(x)\right)^2dx\\
&=\int_a^bw(x)\left(\sum_{i=0}^na_i\phi_i(x)-f(x)\right)^2dx\\
\end{align*}
$$
的最小值。其中 $w(x)$ 是权函数。

一阶导必要条件可得：
$$
\begin{align*}
\frac{\partial I}{\partial a_k}&=2\int_a^bw(x)\left(\sum_{i=0}^na_i\phi_i(x)-f(x)\right)\phi_k(x)\mathrm{d}x\\
&=2\left(\sum_{i=0}^na_i\int_a^bw(x)\phi_i(x)\phi_k(x)\mathrm{d}x-\int_a^bw(x)f(x)\phi_k(x)\mathrm{d}x\right)\\
&=0
\end{align*}
$$
即：
$$
\sum_{i=0}^na_i\langle \phi_k,\phi_i\rangle=\langle \phi_k,f\rangle
$$
这样的方程对于每一个导数都有，因此得到一个法方程和对应的系数矩阵：
$$
\begin{bmatrix}
\langle \phi_0,\phi_0\rangle&\langle \phi_0,\phi_1\rangle&\cdots&\langle \phi_0,\phi_n\rangle\\
\langle \phi_1,\phi_0\rangle&\langle \phi_1,\phi_1\rangle&\cdots&\langle \phi_1,\phi_n\rangle\\
\vdots&\vdots&\ddots&\vdots\\
\langle \phi_n,\phi_0\rangle&\langle \phi_n,\phi_1\rangle&\cdots&\langle \phi_n,\phi_n\rangle
\end{bmatrix}
\begin{bmatrix}
a_0\\a_1\\\vdots\\a_n
\end{bmatrix}=
\begin{bmatrix}
\langle \phi_0,f\rangle\\\langle \phi_1,f\rangle\\\vdots\\\langle \phi_n,f\rangle
\end{bmatrix}
$$
由于其可逆，则解存在且唯一。这样解出来的 $a_k$​ 和基函数的线性组合即为 $\phi^*(x)$​。

由法方程得到
$$
\left(\sum_{i=0}^na_i\phi_i(x)-f(x)\right)\phi_k(x)=0\\
\Longrightarrow\langle \phi^*-f,\phi\rangle=0
$$
误差和 $\Phi$ 中任意函数的距离为 0。

以下可得最佳平方逼近和其他函数的距离：
$$
\begin{align*}
\|f-\phi\|^2_2&=\|f-\phi^*+\phi^*-\phi\|^2_2\\
&=\|f-\phi^*\|^2_2+2\langle f-\phi^*,\phi^*-\phi\rangle+\|\phi^*-\phi\|^2_2\\
&=\|f-\phi^*\|^2_2+\|\phi^*-\phi\|^2_2\\
\end{align*}
$$
记 $\delta(x)=f(x)-\phi^*(x)$ ，其二范数为最佳平方逼近的误差。
$$
\|\delta\|^2_2=\langle f-\phi^*,f-\phi^*\rangle=\|f\|^2_2-\langle f,\phi^*\rangle
$$

<blockquote style="border-left: 5px solid #bb4545; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(188, 70, 70, 0.1)">
    例题
</blockquote>

设 $f\in C[0,1]$，$\Phi=\mathrm{span}(1,x,x^2,\dots,x^n)$，求 $f$ 在 $\Phi$ 中的最佳平方逼近多项式的法方程中的系数矩阵。
$$
H=\begin{bmatrix}
1&1/2&\cdots&1/(n+1)\\
1/2&1/3&\cdots&1/(n+2)\\
\vdots&\vdots&\ddots&\vdots\\
1/(n+1)&1/(n+2)&\cdots&1/(2n+1)
\end{bmatrix}
$$
是 Hilbert 矩阵。

### 3.4.6 正交多项式

如果 $[a, b]$ 上的连续函数系 $\{\phi_j\}^n_0$ 满足：
$$
\int^b_aw(x)\phi_j(x)\phi_k(x)\mathrm{d}x=\begin{cases}0&j\ne k\\A^2&j=k\end{cases}
$$
则称 $\{\phi_j\}^n_0$ 是 $[a,b]$ 上带权 $w(x)$ 的正交函数系。

利用 Gram-Schmidt 方法可以构造正交多项式系。

<blockquote style="border-left: 5px solid #bb4545; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(188, 70, 70, 0.1)">
    例题
</blockquote>

$\{\phi_j\}^n_0$ 是 $[0,1]$ 上带权 $w(x)=x$ 且最高项系数为 1 的正交多项式系，$\phi_0(x)=1$，求 $\phi_1(x),\phi_2(x)$。

> 令 $\phi_1(x)=x+a$，则：
> $$
> \langle\phi_0(x),\phi_1(x)\rangle=\int_0^1x(x+a)\ \mathrm{d}x=\frac{1}{3}+\frac{1}{2}a=0\Rightarrow a=-\frac{2}{3}
> $$
> 令 $\phi_2(x)=x^2+a\phi_1(x)+b\phi_1(x)$，则：
> $$\begin{gather*}
> \langle\phi_0(x),\phi_2(x)\rangle=\int_0^1x^3+bx\ \mathrm{d}x=\frac{1}{4}+\frac{1}{2}b\Rightarrow b=-\frac{1}{2}\\
> \langle\phi_1(x),\phi_2(x)\rangle=\int_0^1x^3\left(x-\frac{2}{3}\right)+ax\left(x-\frac{2}{3}\right)^2\ \mathrm{d}x=\\\frac{1}{4}+\frac{1}{2}b\Rightarrow a=-\frac{1}{2}\end{gather*}
> $$
> 

#### Legendre 多项式

$[-1, 1]$ 区间上带权 $w(x)=1$ 的正交多项式
$$
L_0(x)=1\\
L_n(x)=\frac{1}{2^nn!}\frac{\mathrm d^n}{\mathrm dx^n}[(x^2-1)^n]
$$
化最高次项的系数为 1：
$$
P_0(x)=1\\
P_n(x)=\frac{n!}{(2n!)^2}\frac{\mathrm d^n}{\mathrm dx^n}[(x^2-1)^n]
$$
具有正交性。

奇偶性：
$$
P_n(-x)=(-1)^nP_n(x)
$$
递推关系：
$$
(n+1)P_{n+1}(x)=(2n+1)xP_n(x)-nP_{n-1}(x)
$$

#### Chebyshev 多项式

$[-1, 1]$ 区间上带权 $w(x)=1/\sqrt{1-x^2}$ 的正交多项式
$$
T_n(x)=\cos(n\arccos(x))
$$
具有正交性。

奇偶性：
$$
T_n(-x)=(-1)^nT_n(x)
$$
递推关系：
$$
T_{n+1}(x)=2xT_n(x)-P_{n-1}(x)
$$

