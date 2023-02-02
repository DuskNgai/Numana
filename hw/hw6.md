## Question 1

> 用二分法求方程 $x\sin x-1=0$ 在 $[0,2]$ 内的根的近似值，求：
>
> (1) 为使近似根误差不超过 $10^{-4}$ 所需要的等分次数。
>
> (2) 经过 3 次等分后的近似根。

(1)
$$
n>\frac{\ln2-\ln(10^{-4})}{\ln2}=14.28771237954945
$$
因此 15 次。

(2)

第一次等分：$x=1$。其中 $\sin1-1<0$，选取 $[1,2]$ 区间。

第二次等分：$x=3/2$。其中 $(3/2)\sin(3/2)-1>0$，选取 $[1,3/2]$ 区间。

第二次等分：$x=5/4$。即近似根。

## Question 2

> 设方程 $x=\phi(x)$ 的根是 $x^*$，证明：若 $|\phi'(x)|\ge1$，则由迭代格式 $x_{k+1}=\phi(x_k),x_0\ne x^*$ 所产生的序列 $x_k$ 不收敛到 $x^*$。

根据 Lagrange 中值定理，在 $x_k$ 和 $x^*$ 之间存在 $c_i$ 满足：
$$
\frac{\phi(x_k)-\phi(x^*)}{x_k-x^*}=\phi'(c_k)\\
x_{k+1}-x^*=\phi'(c_k)(x_k-x^*)
$$
带入 $e_k = |x_k - x^*|$，则：
$$
e_{k+1}=|\phi'(c_k)|e_k
$$
如果 $S=|\phi'(x^*)|\ge1$，则在 $x^*$​ 的邻域上满足
$$
1<(S+1)/2<|\phi'(x^*)|
$$
则：
$$
e_{k+1}=|\phi'(c_k)|e_k>e_k
$$
因此产生的序列 $x_k$ 不收敛到 $x^*$。

## Question 3

> 记：
> $$
> \phi(x)=\frac{2-e^x+x^2}{3}
> $$
> 用 $x_{k+1}=\phi(x_k)$ 求解方程 $2-e^x+x^2-3x=0$ 在 $[0,0.5]$ 的根，试讨论上述迭代的局部收敛性。

$$
2-e^x+x^2-3x=0\\
x=\frac{2-e^x+x^2}{3}=\phi(x)
$$

而在 $[0,0.5]$ 上：
$$
\phi'(x)=\frac{-e^x+2x}{3}<0,\phi''(x)=\frac{-e^x+2}{3}>0
$$
因此，$\phi'(x)$ 在 $[0,0.5]$ 上单调递增，此时：
$$
-1<\phi'(0)=-\frac{1}{3}<0,-1<\phi'(0.5)=\frac{1-e^{0.5}}{3}<0\\
\phi(0)=\frac{1}{3}<0.5,\phi(0.5)=\frac{2.25-e^{0.5}}{3}>0
$$
因此 $\phi(x)\subset[0, 0.5]$，且 $|\phi'(x)|<1$。上述迭代的局部收敛。

迭代 5 次可得一解：
$$
x_0=0,x_1=0.333,x_2=0.238\\
x_3=0.262,x_4=0.256,x_5=0.256
$$

## Question 4

> 给定方程 $\sin x-x^2+2=0$。
>
> (1) 验证方程在区间 $[1,2]$ 内具有唯一的根
>
> (2) 若用迭代公式
> $$
> x_{k+1}=\sqrt{2+\sin x_k},k\in\{0,1,\dots\}
> $$
> 求解上述方程，并取 $x_0\in[1,2]$，试分析迭代是否具有整体收敛性。

(1)

令
$$
f(x)=\sin x-x^2+2
$$
此时：
$$
f(1)=\sin1+1>0,f(2)=\sin2-2<0,f'(x)=\cos x-2x
$$
在区间 $[1,2]$ 内，$f'(x)<0$。因此区间上单调递减且区间两端异号，则有唯一的根。

(2)

令
$$
\phi(x)=\sqrt{2+\sin x}
$$
则：
$$
\phi'(x)=\frac{\cos x}{2\sqrt{2+\sin x}}
$$
显然 $|\phi'(x)|<1$。则收敛。

## Question 5

> 已知方程 $f(x)=e^{-x}-x=0$ 在区间 $[0.5, 0.6]$ 内有唯一的实根。
>
> (1) 试判断以下两种求上述方程根的迭代格式的局部收敛性，并说明理由。
> $$
> (1)\ x_{n+1}=-\ln{x_n},x_0>0\\
> (2)\ x_{n+1}=e^{-x_n},x_0>0
> $$
> (2) 方程 $f(x)=0$ 的根就是 $y=f(x)$ 的反函数 $x=g(y)$ 在 $y=0$ 时候 $x$ 的值。已知下列数据表是 $y=f(x)=e^{-x}-x$ 的一组数
>
> | $x$  |  0.50   |  0.55   |   0.60   |
> | :--: | :-----: | :-----: | :------: |
> | $y$  | 0.10653 | 0.02695 | -0.05119 |
>
> 求出 $g(y)$ 的插值多项式，在此基础上求方程 $f(x)=e^{-x}-x=0$ 的近似根。

(1)

(1.1)

令 $\phi(x)=-\ln x$，则
$$
\phi’(x)=-\frac{1}{x}
$$
显然在区间 $[0.5, 0.6]$ 内，$|\phi'(x)|>1$，因此不收敛。

(1.2)

令 $\phi(x)=e^{-x}$​，则
$$
\phi'(x)=-e^{-x}
$$
显然在区间 $[0.5, 0.6]$ 内，$|\phi'(x)|<1$。且 $\phi(x)$ 单调递减。$\phi(0.5)=0.607, \phi(0.6)=0.549$，$\phi(0.549)=0.578,\phi(0.607)=0.545$。因此局部收敛。

(2)

用 Lagrange 插值得到多项式为：
$$
g(y)=0.073412y^2-0.638098y+0.567143
$$
因此近似根为：$x=0.567143$。

## Question 6

> 给定 $f(x)=e^{2x}-1-2x-2x^2=0$。
>
> (1) 确定 $x=0$ 为方程的几重根。
>
> (2) 构造二阶收敛的 Newton 迭代公式。

(1)
$$
f'(x)=2e^{2x}-2-4x,f'(0)=0\\
f''(x)=4e^{2x}-4,f''(0)=0\\
f'''(x)=8e^{2x},f'''(0)=8\\
$$
因此为 3 重根。

(2)

$$
x_{k+1}=x_{k}-3\frac{e^{2x}-1-2x-2x^2}{2e^{2x}-2-4x}
$$
为二阶收敛的 Newton 迭代公式。

## Question 7

> 给定 $c>0$，且 Newton 迭代法解方程：
> $$
> f(x)=\frac{1}{x}-c=0
> $$
> (1) 写出 Newton 迭代格式。
>
> (2) 记 $r_k=c(x^*-x_k)$，$x^*=1/c$ 为方程的解，$x_k$ 为 $k$ 次迭代值，证明：$r_{k+1}=r^2_k$。
>
> (3) 给出初始迭代值 $x_0$ 的范围，使得迭代收敛。

(1)
$$
f'(x)=-\frac{1}{x^2}
$$
因此迭代格式为：
$$
x_{k+1}=x_{k}-\frac{\frac{1}{x_k}-c}{-\frac{1}{x_k^2}}=2x_k-cx_k^2
$$
(2)
$$
r_{k+1}=c(x^*-x_{k+1})=c(x^*-2x_k+cx_k^2)=c^2(x^*-x_k)^2=r_k^2
$$
(3)

当 $|2x-cx^2|<1$ 时候，迭代收敛
$$
|2x-cx^2|<1
$$

因此 $1/2c<x_0<3/2c$。

## Question 8

> 应用牛顿法于方程 $f(x)=x^n-a=0$ 和 $f(x)=1-a/x^n=0$，分别导出求 $\sqrt[n]{a}$ 的迭代公式，并求 $\lim_{k\to\infty}(\sqrt[n]{a}-x_{k+1})/(\sqrt[n]{a}-x_k)^2$。

第一种方法：
$$
x_{k+1}=x_k-\frac{f(x_k)}{f'(x_k)}=x_k-\frac{x_k^n-a}{nx_k^{n-1}}=\frac{(n-1)x_k^n+a}{nx_k^{n-1}}
$$
则
$$
\lim_{k\to\infty}\frac{(\sqrt[n]{a}-x_{k+1})}{(\sqrt[n]{a}-x_k)^2}=\frac{f''(\sqrt[n]{a})}{2f'(\sqrt[n]{a})}=\frac{1-n}{2\sqrt[n]{a}}
$$
第二种方法：
$$
x_{k+1}=x_k-\frac{f(x_k)}{f'(x_k)}=x_k-\frac{1-\frac{a}{x_k^n}}{\frac{na}{x_k^{n-1}}}=\frac{(na+a)x_k-x_k^{n+1}}{na}
$$
则
$$
\lim_{k\to\infty}\frac{(\sqrt[n]{a}-x_{k+1})}{(\sqrt[n]{a}-x_k)^2}=\frac{f''(\sqrt[n]{a})}{2f'(\sqrt[n]{a})}=\frac{n+1}{2\sqrt[n]{a}}
$$

## Question 9

> 给定方程 $f(x)=x^{1+a}-x=0$，其中 $0<a<1$，记 $p_1,p_2$ 是用牛顿迭代法求解上述方程两个根 $x_1=0,x_2=1$ 的收敛阶，求 $p_1,p_2$ 的值。

$$
x_{k+1}=x_k-\frac{x_k^{1+a}-x_k}{(1+a)x_k^a-1}=\frac{ax_k^{1+a}}{(1+a)x_k^a-1}
$$

则：
$$
\lim_{k\to\infty}\frac{|x_{k+1}-0|}{|x_k-0|^{p_1}}=\lim_{k\to\infty}\frac{|\frac{ax_k^{1+a}}{(1+a)x_k^a-1}|}{|x_k|^{p_1}}=\lim_{k\to\infty}\frac{|a||x_k|^{1+a-p_1}}{|(1+a)x_k^a-1|}
$$
由于下面趋向 $-1$，上面趋向于 $0$，为了使得比值不为 $0$，则 $p_1=1+a$。
$$
\phi(x)=\frac{ax^{1+a}}{(1+a)x^a-1}\\
\phi'(x)=0,\phi''(x)\ne0
$$
由于下面趋向 $a$，上面趋向于 $0$，为了使得比值不为 $0$，则 $p_2=2$。

## Question 10

> 设 $a>0,b>0,c>0,I_0>0$，
>
> (1) 讨论迭代格式 $I_k=c/(b+aI_{k-1})$ 的敛散性。
>
> (2) 若该格式是线性收敛，再构造一个平方收敛的迭代格式。

(1)
$$
\phi(x)=\frac{c}{b+ax},\phi'(x)=\frac{-ac}{(b+ax)^2}
$$
当：
$$
\left|\frac{-ac}{(b+ax)^2}\right|<1
$$
即：
$$
x>\max\left(\frac{-b+\sqrt{ac}}{a},0\right)
$$
时候收敛。

(2)

若该格式是线性收敛，则：
$$
\lim_{k\to\infty}\frac{|x_{k+1}-x^*|}{|x_k-x^*|}=\lim_{k\to\infty}\frac{|\frac{c}{b+ax_k}-x^*|}{|x_k-x^*|}
$$
令
$$
u(x)=(x-x^*)(x-\frac{c}{b+ax})
$$
因此：
$$
u'(x)=(1+\frac{2ac}{(b+ax)^2})(x-x^*)+(x-\frac{c}{b+ax})\\
u'(x^*)=0\\
u''(x)=-\frac{6a^2c}{(b+ax)^3}(x-x^*)+2(1+\frac{2ac}{(b+ax)^2})\\
u''(x^*)\ne0
$$
因此 $u(x)$ 二阶收敛。

## Question 11

> 写出方程组
> $$
> \begin{cases}
> x^2+y^2=4\\
> x^2-y^2=1
> \end{cases}
> $$
> 的 Newton 迭代格式。

改写得到：
$$
F_1(x,y)=x^2+y^2-4=0\\F_2(x,y)=x^2-y^2-1=0\\
$$
因此 Jacobi 矩阵为：
$$
\begin{bmatrix}
\dfrac{\partial{F_1}}{\partial{x}}&\dfrac{\partial{F_1}}{\partial{y}}\\
\dfrac{\partial{F_2}}{\partial{x}}&\dfrac{\partial{F_2}}{\partial{y}}
\end{bmatrix}=
\begin{bmatrix}
2x&2y\\
2x&-2y
\end{bmatrix}
$$
逆矩阵为：
$$
\begin{bmatrix}
\dfrac{1}{4x}&\dfrac{1}{4y}\\
\dfrac{1}{4x}&-\dfrac{1}{4y}
\end{bmatrix}
$$
因此迭代格式为：
$$
\begin{bmatrix}
{x^{(k+1)}}\\{y^{(k+1)}}
\end{bmatrix}=\begin{bmatrix}
{x^{(k)}}\\{y^{(k)}}
\end{bmatrix}-\begin{bmatrix}
\dfrac{1}{4x^{(k)}}&\dfrac{1}{4y^{(k)}}\\
\dfrac{1}{4x^{(k)}}&-\dfrac{1}{4y^{(k)}}
\end{bmatrix}
\begin{bmatrix}
{x^{(k)}}^2+{y^{(k)}}^2-4\\{x^{(k)}}^2-{y^{(k)}}^2-1
\end{bmatrix}
$$

## Question 12

> 非线性方程组：
> $$
> \begin{cases}
> 4x_1+\cos{x_1}-x_2=0\\
> \frac{1}{8}x_1^2-x_1+4x_2=0\\
> \end{cases}
> $$
> 在 $D=\{(x_1,x_2)|[-1,1]\times[0,2]\}$ 有解 $x^*$。
>
> (1) 写出求解上述非线性方程组的 Newton 迭代公式。
>
> (2) 取 $x^{(0)}=(0.5,0.5)$，用 Newton 迭代计算 $x^{(1)}$。

(1) 改写得到：
$$
F_1(x,y)=4x_1+\cos{x_1}-x_2=0\\F_2(x,y)=\frac{1}{8}x_1^2-x_1+4x_2=0\\
$$
因此 Jacobi 矩阵为：
$$
\begin{bmatrix}
\dfrac{\partial{F_1}}{\partial{x_1}}&\dfrac{\partial{F_1}}{\partial{x_2}}\\
\dfrac{\partial{F_2}}{\partial{x_1}}&\dfrac{\partial{F_2}}{\partial{x_2}}
\end{bmatrix}=
\begin{bmatrix}
4-\sin{x_1}&-1\\
\frac{1}{4}x_1-1&4
\end{bmatrix}
$$
因此迭代格式为：
$$
\begin{bmatrix}
4-\sin{x_1^{(k)}}&-1\\
\frac{1}{4}x_1^{(k)}-1&4
\end{bmatrix}
\begin{bmatrix}
\Delta{x_1^{(k)}}\\
\Delta{x_2^{(k)}}\\
\end{bmatrix}=-
\begin{bmatrix}
4x_1+\cos{x_1}-x_2\\
\frac{1}{8}x_1^2-x_1+4x_2\\
\end{bmatrix}
$$

$$
\begin{bmatrix}
x_1^{(k+1)}\\
x_2^{(k+1)}\\
\end{bmatrix}=\begin{bmatrix}
\Delta{x_1^{(k)}}\\
\Delta{x_2^{(k)}}\\
\end{bmatrix}+\begin{bmatrix}
x_1^{(k)}\\
x_2^{(k)}\\
\end{bmatrix}
$$

(2)

带入 $x^{(0)}=(0.5,0.5)$
$$
\begin{bmatrix}
4-\sin{0.5}&-1\\
-\frac{7}{8}&4
\end{bmatrix}
\begin{bmatrix}
\Delta{x_1^{(0)}}\\
\Delta{x_2^{(0)}}\\
\end{bmatrix}=-
\begin{bmatrix}
1.5+\cos{0.5}\\
\frac{49}{32}\\
\end{bmatrix}
$$
得到：
$$
\begin{bmatrix}
\Delta{x_1^{(0)}}\\
\Delta{x_2^{(0)}}\\
\end{bmatrix}=\begin{bmatrix}
-0.83602114\\
-0.56569213\\
\end{bmatrix}
$$
因此：
$$
\begin{bmatrix}
x_1^{(1)}\\
x_2^{(1)}\\
\end{bmatrix}=\begin{bmatrix}
\Delta{x_1^{(0)}}\\
\Delta{x_2^{(0)}}\\
\end{bmatrix}+\begin{bmatrix}
x_1^{(0)}\\
x_2^{(0)}\\
\end{bmatrix}=\begin{bmatrix}
-0.33602114\\
0.06569213\\
\end{bmatrix}
$$

