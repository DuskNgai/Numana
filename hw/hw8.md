## Question 1

> 用幂法求矩阵
> $$
> A=\begin{bmatrix}6&2&1\\2&3&1\\1&1&1\end{bmatrix}
> $$
> 按模最大的特征值和对应的特征向量。取初始值 $x^{(0)}=(1,1,1)^T$，迭代两次。

迭代一次
$$
v_1=Ax^{(0)}=\begin{bmatrix}6&2&1\\2&3&1\\1&1&1\end{bmatrix}\begin{bmatrix}1\\1\\1\end{bmatrix}=\begin{bmatrix}9\\6\\3\end{bmatrix}\\
x^{(1)}=\frac{v_1}{\max(v_1)}=\begin{bmatrix}1\\2/3\\1/3\end{bmatrix}
$$
迭代两次
$$
v_2=Ax^{(1)}=\begin{bmatrix}6&2&1\\2&3&1\\1&1&1\end{bmatrix}\begin{bmatrix}1\\2/3\\1/3\end{bmatrix}=\begin{bmatrix}23/3\\13/3\\6/3\end{bmatrix}\\
x^{(2)}=\frac{v_1}{\max(v_1)}=\begin{bmatrix}1\\13/23\\6/23\end{bmatrix}
$$

## Question 2

> 设矩阵
> $$
> A=\begin{bmatrix}3&0\\2&1\end{bmatrix}
> $$
> 取 $x_0=\begin{bmatrix}1\\2\end{bmatrix}$，构造如下迭代：
> $$
> \begin{cases}
> y_k=Ax_k\\
> x_{k+1}=y_k/3
> \end{cases}
> $$
> 并设 $\tau_k=\langle Ax_k,x_k\rangle/\langle x_k,x_k\rangle$。计算 $\lim_{k\to\infty}x_k,\lim_{k\to\infty}\tau_k$，并解释极限值与矩阵 $A$ 的关系。

设 $x_0=(a,b)^T$，则：
$$
x_1=\begin{bmatrix}a\\(2a+b)/3\end{bmatrix}\\
x_2=\begin{bmatrix}a\\(8a+b)/9\end{bmatrix}\\
\vdots\\
x_k=\begin{bmatrix}a\\(1-3^{-k})a+3^{-k}b\end{bmatrix}\\
$$
因此带入值：
$$
\lim_{k\to\infty}x_k=\lim_{k\to\infty}x_k\begin{bmatrix}a\\(1-3^{-k})a+3^{-k}b\end{bmatrix}=\begin{bmatrix}a\\a\end{bmatrix}=\begin{bmatrix}1\\1\end{bmatrix}
$$

$$
\tau_k=\frac{\langle Ax_k,x_k\rangle}{\langle x_k,x_k\rangle}\\
=\frac{3a^2+2a[(1-3^{-k})a+3^{-k}b]+(1-3^{-k})^2a^2+3^{-2k}b^2+2(1-3^{-k})3^{-k}ab}{a^2+(1-3^{-k})^2a^2+3^{-2k}b^2+2(1-3^{-k})3^{-k}ab}
$$

因此：
$$
\lim_{k\to\infty}\tau_k=3
$$

## Question 3

> 已知
> $$
> A=\begin{bmatrix}2&0&2\\1&1&2\\0&-1&-3\end{bmatrix}
> $$
> 的三个特征值分别为：
> $$
> \lambda_1=-\frac{1+\sqrt{17}}{2},\lambda_1=-\frac{1-\sqrt{17}}{2},\lambda_3=1
> $$
> (1) LU 分解。
>
> (2) 取初始向量 $x^{(0)}=(1,1,0)^T$，用反幂法计算 $\lambda_3$ 与对应特征向量的近似值（只迭代一次）。

(1)
$$
\begin{bmatrix}1&0&0\\1/2&1&0\\0&-1&1\end{bmatrix}\begin{bmatrix}2&0&2\\0&1&1\\0&0&-2\end{bmatrix}
$$
(2)
$$
z^{(1)}=A^{-1}x^{(0)}=\begin{bmatrix}1/4&1/2&1/2\\-3/4&3/2&1/2\\1/4&-1/2&-1/2\end{bmatrix}\begin{bmatrix}1\\1\\0\end{bmatrix}\\
=\begin{bmatrix}3/4\\3/4\\-1/4\end{bmatrix}\\
$$

$$
x^{(1)}=\frac{z^{(1)}}{\max(z^{(1)})}=\begin{bmatrix}1\\1\\-1/3\end{bmatrix}
$$

## Question 4

> 已知 $a=[3,4]^T$，$b=[0,1]^T$。求 Householder 阵 $H$，使得 $Ha=kb$，其中 $k>0$。

## Question 5

> 求矩阵
> $$
> A=\begin{bmatrix}3&0&0\\0&3&2\\0&2&3\end{bmatrix}
> $$
> 的 QR 分解。
