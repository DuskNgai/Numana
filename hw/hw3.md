# Homework 3

## Question 1

> 给定线性代数方程组 $Ax=b$，其中，
> $$
> A=\begin{bmatrix}
> 1&2&-2\\
> 1&1&1\\
> 2&2&1
> \end{bmatrix}\quad
> b=\begin{bmatrix}
> 1\\3\\5
> \end{bmatrix}
> $$
> 设求解 $Ax=b$ 的一个迭代格式为：
> $$
> \begin{cases}
> x^{(k+1)}_3=5-2x^{(k)}_1-2x^{(k)}_2\\
> x^{(k+1)}_2=3-x^{(k+1)}_3-x^{(k)}_1\\
> x^{(k+1)}_1=1-2x^{(k+1)}_2+2x^{(k+1)}_3
> \end{cases}
> $$
> 写出上述迭代格式的迭代矩阵，并判断该迭代格式的敛散性？

设 $A=D-L-U$，其中 $D$ 是 $A$ 对角元素组成的矩阵，$L$ 是 $A$ 下三角元素的相反数组成的矩阵，$U$ 是 $A$ 上三角元素的相反数组成的矩阵。

则
$$
\begin{align*}
x^{(k+1)}_i&=\frac{1}{a_{ii}}\left(b_i-\sum_{j=i+1}^{n}a_{ij}x^{(k+1)}_j-\sum_{j=1}^{i-1}a_{ij}x^{(k)}_j\right)\\
Dx^{(k+1)}&=\left(b+Ux^{(k+1)}+Lx^{(k)}\right)\\
(D-U)x^{(k+1)}&=Lx^{(k)}+b\\
x^{(k+1)}&=(D-U)^{-1}Lx^{(k)}+(D-U)^{-1}b\\
\end{align*}
$$

$$
G=(D-U)^{-1}L=\begin{bmatrix}
-6&-8&0\\
1&2&0\\
-2&-2&0
\end{bmatrix}
$$

收敛性：
$$
\begin{align*}
\det(G)&=\begin{vmatrix}
\lambda+6&8&0\\
-1&\lambda-2&0\\
2&2&\lambda
\end{vmatrix}\\
&=\lambda^3+4\lambda^2-4\lambda\\
&=\lambda(\lambda^2+4\lambda+4-8)\\
&=\lambda(\lambda+2+2\sqrt{2})(\lambda+2-2\sqrt{2})
\end{align*}
$$
$\rho(G)=2+2\sqrt{2}>1$ 因此不收敛。

## Question 2

> 设 
> $$
> A=\begin{bmatrix}
> 3&2\\1&2
> \end{bmatrix}\quad
> b=\begin{bmatrix}
> 3\\-1
> \end{bmatrix}
> $$
> 用迭代公式 $x^{(k+1)}=x^{(k)}+\alpha(Ax^{(x)}-b)$ 求解方程组 $Ax=b$，问：当 $\alpha$ 为何值时，迭代收敛？

$$
\begin{align*}
x^{(k+1)}&=x^{(k)}+\alpha(Ax^{(k)}-b)\\
x^{(k+1)}&=(I+\alpha A)x^{(k)}-\alpha b\\
\end{align*}
$$

当 $\rho(I+\alpha A)<1$ 时候，迭代收敛，则：
$$
\begin{align*}
\det(\lambda I-(I+\alpha A))&=\det((\lambda-1)I-\alpha A)\\
&=\begin{bmatrix}
(\lambda-1)-3\alpha&-2\alpha\\-\alpha&(\lambda-1)-2\alpha
\end{bmatrix}\\
&=[(\lambda-1)-3\alpha][(\lambda-1)-2\alpha]-2\alpha^2\\
&=(\lambda-1)^2-5\alpha(\lambda-1)+4\alpha^2\\
&=(\lambda-1-\alpha)(\lambda-1-4\alpha)
\end{align*}
$$
特征值为：
$$
\lambda=1+\alpha\in(-1, 1)\Longrightarrow\alpha\in(-2, 0)\\
\lambda=1+4\alpha\in(-1, 1)\Longrightarrow\alpha\in(-0.5, 0)\\
$$
因此 $\alpha\in(-0.5, 0)$。

## Question 3

> 设有方程组 $Ax=b$，其中 $A$ 为对称正定阵. 用迭代公式 $x^{(k+1)}=x^{(k)}+\omega(b-Ax^{(k)})$ 求解上述方程组，证明当 $0<\omega<2/\beta$ 时迭代收敛。（其中 $0<\alpha\le\lambda(A)\le\beta$ ）

$$
\begin{align*}
x^{(k+1)}&=x^{(k)}+\omega(b-Ax^{(k)})\\
x^{(k+1)}&=(I-\omega A)x^{(k)}+\omega b\\
\end{align*}
$$

当 $\rho(I-\omega A)<1$ 时候，迭代收敛，则：
$$
\det(\lambda I-(I-\omega A))=\det((\lambda-1)I+\omega A)
$$
其中 $-1<\lambda<1$

由于 $\omega Ax=\omega\mu x$，则矩阵数乘 $\omega$ 倍可得到对应的特征值也放缩 $\omega$ 倍。

因此：
$$
0>\lambda-1\ge-\omega\beta>-2\Longrightarrow0<\omega<\frac{2}{\beta}\\
0>\lambda-1\le-\omega\alpha>-2\Longrightarrow0<\omega<\frac{2}{\alpha}\\
$$
由于 $\beta$ 是 $A$ 最大的特征值，则 $2/\beta$ 最小，因此当 $0<\omega<2/\beta$ 时迭代收敛。

## Question 4

> 设用迭代格式：$x^{(k+1)}=Bx^{(k)}+f$ 求解线性代数方程组 $Ax=b$ 收敛，证明：当 $0<\omega<1$ 时，迭代格式 $x^{(k+1)}=[(1-\omega)I+\omega B]x^{(k)}+\omega f$ 也收敛。

迭代格式：$x^{(k+1)}=Bx^{(k)}+f$ 求解线性代数方程组 $Ax=b$ 收敛，则 $\rho(B)<1$。
$$
\det(\lambda I-[(1-\omega)I+\omega B])=\det((\lambda-1+\omega)I-\omega B)
$$
因此：
$$
-\omega\rho(B)<(\lambda-1+\omega)<\omega\rho(B)\\
-\omega\rho(B)+1-\omega<\lambda<\omega\rho(B)+1-\omega\\
-\omega\rho(B)+1-\omega=-\omega(\rho(B)+1)+1>-2\omega+1>-1\\
\omega\rho(B)+1-\omega=\omega(\rho(B)-1)+1<1
$$
因此 $\rho((1-\omega)I+\omega B)<1$，迭代格式 $x^{(k+1)}=[(1-\omega)I+\omega B]x^{(k)}+\omega f$ 也收敛。

## Question 5

> 设 $A$ 为对称正定矩阵，$\omega>0$ 为常数，用以下的迭代公式：
> $$
> x^{(k+1)}=x^{(k)}-\omega\left[A\left(\frac{x^{(k+1)}+x^{(k)}}{2}\right)-b\right]
> $$
> 求解方程 $Ax=b$。
>
> (1) 证明：对任意初始向量，上述迭代格式均收敛。
>
> (2) 取
> $$
> A=\begin{bmatrix}
> 2&1\\1&2
> \end{bmatrix}\quad
> \omega=2
> $$
> 计算上述迭代格式的收敛速度。

(1)
$$
\begin{align*}
x^{(k+1)}&=x^{(k)}-\omega\left[A\left(\frac{x^{(k+1)}+x^{(k)}}{2}\right)-b\right]\\
x^{(k+1)}&=x^{(k)}-\frac{\omega A}{2}\left(x^{(k+1)}+x^{(k)}\right)+\omega b\\
\left(I+\frac{\omega A}{2}\right)x^{(k+1)}&=\left(I-\frac{\omega A}{2}\right)x^{(k)}+\omega b\\
x^{(k+1)}&=\left(I+\frac{\omega A}{2}\right)^{-1}\left(I-\frac{\omega A}{2}\right)x^{(k)}+\left(I+\frac{\omega A}{2}\right)^{-1}\omega b\\
\end{align*}
$$
由于 $A$ 是对称正定矩阵，则 $I+\frac{\omega A}{2}$ 可逆；同时 $A$ 可被对角化为 $Q^TDQ$，其中 $Q$ 是正交矩阵，$D$ 是对角矩阵。迭代矩阵 $B$ 为：
$$
\begin{align*}
B&=\left(I+\frac{\omega A}{2}\right)^{-1}\left(I-\frac{\omega A}{2}\right)\\
&=\left[Q^T\left(I+\frac{\omega}{2}D\right)Q\right]^{-1}\left[Q^T\left(I-\frac{\omega}{2}D\right)Q\right]\\
&=\left[Q^{T}\left(I+\frac{\omega}{2}D\right)^{-1}Q\right]\left[Q^T\left(I-\frac{\omega}{2}D\right)Q\right]\\
&=Q^{T}\left(I+\frac{\omega}{2}D\right)^{-1}\left(I-\frac{\omega}{2}D\right)Q
\end{align*}
$$
这里 $B$ 的每个特征值，即 $(I+\frac{\omega}{2}D)^{-1}(I-\frac{\omega}{2}D)$ 的每个元素都是
$$
\frac{(1-\frac{\omega}{2}\lambda_i)}{(1+\frac{\omega}{2}\lambda_i)}=\frac{(-1-\frac{\omega}{2}\lambda_i)+2}{(1+\frac{\omega}{2}\lambda_i)}=\frac{2}{(1+\frac{\omega}{2}\lambda_i)}-1>-1\\
\frac{(1-\frac{\omega}{2}\lambda_i)}{(1+\frac{\omega}{2}\lambda_i)}=\frac{(1+\frac{\omega}{2}\lambda_i)-\omega\lambda_i}{(1+\frac{\omega}{2}\lambda_i)}=1-\frac{\omega\lambda_i}{(1+\frac{\omega}{2}\lambda_i)}<1
$$
其中 $\lambda_i>0$ 是 $A$ 的特征值。因此 $\rho(B)<1$。

因此迭代公式永远收敛。

(2)

$A$ 的特征值为 $\lambda=1,3$，则 $B$ 为
$$
\left(I+\frac{\omega}{2}D\right)^{-1}\left(I-\frac{\omega}{2}D\right)=\begin{bmatrix}0&0\\0&-\frac{1}{2}\end{bmatrix}
$$
因此 $\rho(B)=\frac{1}{2}$，收敛速度为
$$
-\ln\rho(B)=\ln2
$$

## Question 6

> 设 $A$ 为 $n$ 阶对称正定矩阵，其最大和最小特征值为 $\lambda_1,\lambda_n$，用迭代格式
> $$
> x^{(k+1)}=(I-\omega A)x^{(k)}+\omega b
> $$
> 求解方程 $Ax=b$。
>
> (1) 求 $\omega$ 的取值范围，使得迭代收敛
>
> (2) 求 $\omega$ 的值，使得迭代格式具有最大的收敛速度。

(1)

由于 $A$ 是对称正定矩阵，$A$ 可被对角化为 $Q^TDQ$，其中 $Q$ 是正交矩阵，$D$ 是对角矩阵。
$$
B=I-\omega A=Q^T(I-\omega D)Q
$$
因此 $B$ 最大的特征值和最小的特征值分别为 $1-\omega\lambda_n,1-\omega\lambda_1$，可以得到：
$$
0<\omega<\frac{2}{\lambda_1}
$$
(2)

收敛速度最大则矩阵的谱半径要足够的小，而谱半径为：
$$
\rho(B)=\max(|1-\omega\lambda_1|,|1-\omega\lambda_n|)
$$
则：
$$
1-\omega\lambda_n=-1+\omega\lambda_1\\
\omega=\frac{2}{\lambda_1+\lambda_n}
$$
迭代格式具有最大的收敛速度。

## Question 7

> 给定线性代数方程组 $Ax=b$，其中 
> $$
> A=\begin{bmatrix}
> \alpha&2&1\\
> 2&\alpha&-1\\
> 1&1&0.5
> \end{bmatrix}\quad
> b=\begin{bmatrix}
> 1\\1\\1
> \end{bmatrix}\quad
> \alpha\ne0
> $$
> (1) 写出求解上述方程组的 Jacobi 迭代格式
>
> (2) 试确定 $\alpha$ 的取值范围以确保 Jacobi 迭代格式收敛

(1)
$$
x^{(k+1)}=(I-D^{-1}A)x^{(k)}+D^{-1}b\\
x^{(k+1)}=\begin{bmatrix}
0&-\frac{2}{\alpha}&-\frac{1}{\alpha}\\
-\frac{2}{\alpha}&0&\frac{1}{\alpha}\\
-2&-2&0
\end{bmatrix}x^{(k)}+\begin{bmatrix}
\frac{1}{\alpha}\\
\frac{1}{\alpha}\\
2
\end{bmatrix}\\
$$
(2)
$$
\det(\lambda I-J)=\det(D^{-1})\det(\lambda D-(L+U))=0
$$

$$
\begin{align*}
\det(\lambda D-(L+U))&=\begin{vmatrix}
\lambda\alpha&2&1\\
2&\lambda\alpha&-1\\
1&1&0.5\lambda
\end{vmatrix}\\
&=0.5\alpha^2\lambda^3-2+2+\lambda\alpha-2\lambda-\lambda\alpha\\
&=0.5\alpha^2\lambda^3-2\lambda\\
&=0.5\lambda(\alpha^2\lambda^2-4)
\end{align*}
$$

因此 $\lambda=0,\pm2/\alpha$，因此 $\alpha\in(-\infty,-2)\cup(2,\infty)$。

## Question 8

> 写出方程组
> $$
> \begin{cases}
> x_1+2x_2-2x_3=1\\
> x_1+x_2+x_3=1\\
> 2x_1+2x_2+x_3=1
> \end{cases}
> $$
> 的 Jacobi 迭代法和 Gauss-Seidel 迭代法的计算公式，并判断收敛性？

$$
A=\begin{bmatrix}
1&2&-2\\
1&1&1\\
2&2&1
\end{bmatrix}\quad
b=\begin{bmatrix}
1\\1\\1
\end{bmatrix}
$$

Jacobi
$$
x^{(k+1)}=(I-D^{-1}A)x^{(k)}+D^{-1}b\\
x^{(k+1)}=\begin{bmatrix}
0&-2&2\\
-1&0&-1\\
-2&-2&0
\end{bmatrix}x^{(k)}+\begin{bmatrix}
1\\1\\1
\end{bmatrix}
$$
Gauss-Seidel
$$
x^{(k+1)}=(D-L)^{-1}Ux^{(k)}+(D-L)^{-1}b\\
x^{(k+1)}=\begin{bmatrix}
0&-2&2\\
0&2&-3\\
0&0&2
\end{bmatrix}x^{(k)}+\begin{bmatrix}
1\\0\\-1
\end{bmatrix}
$$
Jacobi
$$
\begin{align*}
\det(\lambda D-(L+U))&=\begin{vmatrix}
\lambda&2&-2\\
1&\lambda&1\\
2&2&\lambda
\end{vmatrix}\\
&=\lambda^3+4-4-2\lambda-2\lambda+4\lambda\\
&=\lambda^3
\end{align*}
$$
$\rho(J)=0$，谱半径小于 1，收敛。

Gauss-Seidel
$$
\begin{align*}
\det(\lambda(D-L)-U)&=\begin{vmatrix}
\lambda&2&-2\\
\lambda&\lambda&1\\
2\lambda&2\lambda&\lambda
\end{vmatrix}\\
&=\lambda^3-4\lambda^2+4\lambda-2\lambda^2-2\lambda^2+4\lambda^2\\
&=\lambda(\lambda-2)^2
\end{align*}
$$
$\rho(G)=2$，谱半径大于 1，不收敛。

## Question 9

> 给定线性代数方程组 $Ax=b$​，其中
> $$
> A=\begin{bmatrix}
> 1&\omega&\omega\\
> 3\omega&1&0\\
> \omega&0&1
> \end{bmatrix}
> $$
> $\omega$ 为实数。证明：用 Jacobi 迭代和 Gauss-Seidel 迭代解方程 $Ax=b$ 时，同时收敛，同时发散。

Jacobi
$$
\begin{align*}
\det(\lambda D-(L+U))&=\begin{vmatrix}
\lambda&\omega&\omega\\
3\omega&\lambda&0\\
\omega&0&\lambda
\end{vmatrix}\\
&=\lambda^3-4\lambda\omega^2\\
&=\lambda(\lambda+2\omega)(\lambda-2\omega)
\end{align*}
$$
$\rho(J)=|2\omega|$

Gauss-Seidel
$$
\begin{align*}
\det(\lambda(D-L)-U)&=\begin{vmatrix}
\lambda&\omega&\omega\\
3\lambda\omega&\lambda&0\\
\lambda\omega&0&\lambda
\end{vmatrix}\\
&=\lambda^3-4\lambda^2\omega^2\\
&=\lambda^2(\lambda-4\omega^2)
\end{align*}
$$
$\rho(G)=4\omega^2=|2\omega|^2$

此时 $\rho(J)<1$ 等价于 $\rho(G)<1$，因此同时收敛，同时发散。

## Question 10

> 记 $h=1/N$ ($N$ 为正整数)，$x_i=ih,y_i\approx y(x_i),f_i\approx f(x_i)$ 。用差分格式
> $$
> \begin{cases}
> \dfrac{y_{i+1}-2y_{i}+y_{i-1}}{h^2}+\dfrac{y_{i+1}-y_{i-1}}{2h}-2y_i=f_i&i\in\{1,\dots,N-1\}\\
> y_0=0,y_N=0
> \end{cases}
> $$
> 求解二阶常微分方程边值问题
> $$
> \begin{cases}
> y''(x)+y'(x)-2y(x)=f(x)\\
> y(0)=0,y(1)=0
> \end{cases}
> $$
> (1) 写出近似解 $y=[y_1,\dots,y_{N-1}]^T$ 所满足的线性代数方程组
>
> (2) 若用 Gauss-Seidel 迭代解（1）中方程组，问迭代是否收敛？为什么？

(1)
$$
\frac{y_{i+1}-2y_{i}+y_{i-1}}{h^2}+\frac{y_{i+1}-y_{i-1}}{2h}-2y_i=f_i\\
\left(\frac{1}{h^2}-\frac{1}{2h}\right)y_{i-1}-\left(\frac{2}{h^2}+2\right)y_{i}+\left(\frac{1}{h^2}+\frac{1}{2h}\right)y_{i+1}=f_i\\
(2-h)y_{i-1}-(4+4h^2)y_{i}+(2+h)y_{i+1}=2h^2f_i
$$
因此线性方程组满足：
$$
Ay=f
$$
其中
$$
A=\begin{bmatrix}
-4-4h^2&2+h&0&0&0&\cdots\\
2-h&-4-4h^2&2+h&0&0&\cdots\\
0&2-h&-4-4h^2&2+h&0&\cdots\\
&\ddots&\ddots&\ddots&\ddots\\
&&0&2-h&-4-4h^2&2+h\\
&&&0&2-h&-4-4h^2\\
\end{bmatrix}
$$

$$
y=\begin{bmatrix}
y_1\\y_2\\\vdots\\y_{N-1}
\end{bmatrix}\qquad
f=2h^2
\begin{bmatrix}
f_1\\f_2\\\vdots\\f_{N-1}
\end{bmatrix}
$$

(2)

由于 $h$ 为很小的正数，因此
$$
\left|2-h\right|+\left|2+h\right|=4<\left|-4-4h^2\right|
$$
因此 $A$ 是严格对角占有阵，因此其 Gauss-Seidel 迭代收敛。

## Question 11

> 给定以下线性代数方程组 $Ax=b$​，其中
> $$
> A=\begin{bmatrix}
> 1&-2&2\\
> -1&1&-1\\
> -2&-2&1
> \end{bmatrix}\quad b=\begin{bmatrix}
> -1\\0\\-5
> \end{bmatrix}
> $$
>
> (1) 讨论用 Jacobi 迭代求解上述方程的敛散性
>
> (2) 取初始迭代值 $x^{(0)}=\mathbf 0^T$，用 Jacobi迭代计算三次，即计算 $x^{(3)}$
>
> (3) 证明：经过有限次 Jacobi 迭代后能得到方程精确解

(1)
$$
J=I-D^{-1}A=\begin{bmatrix}
0&2&-2\\
1&0&1\\
2&2&0
\end{bmatrix}\quad
f=\begin{bmatrix}
-1\\0\\-5
\end{bmatrix}
$$

$$
\det(\lambda I-J)=\begin{vmatrix}
\lambda&-2&2\\
-1&\lambda&-1\\
-2&-2&\lambda
\end{vmatrix}=\lambda^3-4+4-2\lambda-2\lambda+4\lambda=\lambda^3
$$

因此 $\rho(J)=0$，Jacobi 迭代收敛。

(2)
$$
x^{(1)}=\begin{bmatrix}
0&2&-2\\
1&0&1\\
2&2&0
\end{bmatrix}\begin{bmatrix}0\\0\\0\end{bmatrix}+\begin{bmatrix}-1\\0\\-5\end{bmatrix}
=\begin{bmatrix}-1\\0\\-5\end{bmatrix}\\
x^{(2)}=\begin{bmatrix}
0&2&-2\\
1&0&1\\
2&2&0
\end{bmatrix}\begin{bmatrix}-1\\0\\-5\end{bmatrix}+\begin{bmatrix}-1\\0\\-5\end{bmatrix}
=\begin{bmatrix}9\\-6\\-7\end{bmatrix}\\
x^{(3)}=\begin{bmatrix}
0&2&-2\\
1&0&1\\
2&2&0
\end{bmatrix}\begin{bmatrix}9\\-6\\-7\end{bmatrix}+\begin{bmatrix}-1\\0\\-5\end{bmatrix}
=\begin{bmatrix}1\\2\\1\end{bmatrix}\\
$$
(3)

谱半径为 0，收敛到精确值。

验算可得 $x^{(3)}$ 即是 $Ax=b$ 的精确解。因此可以。

## Question 12

> 已知 
> $$
> A=\begin{bmatrix}
> 2&\alpha&1\\
> 2&2&\alpha\\
> 0&-1&2
> \end{bmatrix}
> $$
> (1) 找出参数 $\alpha$ 的最大范围，使得求解以 $A$ 为系数矩阵的线性代数方程组的 Gauss-Seidel 迭代法收敛
>
> (2) 当取 $\alpha$ 何值时，Gauss-Seidel 迭代法经有限次迭代后得到方程的精确解

(1)
$$
\det(\lambda(D-L)-U)=\begin{vmatrix}
2\lambda&\alpha&1\\
2\lambda&2\lambda&\alpha\\
0&-\lambda&2\lambda
\end{vmatrix}=8\lambda^3-2\lambda^2+2\lambda^2\alpha-4\lambda^2\alpha=2\lambda^2(4\lambda-1-\alpha)
$$
因此 $\lambda=0,(1+\alpha)/4$。

收敛则
$$
-1<\frac{1+\alpha}{4}<1\Longrightarrow-5<\alpha<3
$$
(2)

谱半径为 $0$ 即 $\alpha=-1$。