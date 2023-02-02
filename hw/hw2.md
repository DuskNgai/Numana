# Homework 2

## Question 1

> 设 $A$ 是对称正定矩阵，经过 Gauss 消元法一步以后，$A$ 约化为 
> $$
> \begin{bmatrix}
> a_{11}     &*\\
> \mathbf 0^T&A_2
> \end{bmatrix}
> $$
> 证明：
>
> (1) $A$ 的对角元素 $a_{ii}>0$。
>
> (2) $A_2$ 是对称正定矩阵。
>
> (3) $a^{(2)}_{ii}\le a_{ii}$。
>
> (4) $A$ 的绝对值最大的元素必在对角线上。
>
> (5) $\max_{2\le i,j\le n}|a^{(2)}_{ij}|\le\max_{2\le i,j\le n}|a_{ij}|$

(1)

取向量 $\mathbf x=(0, \dots, 1,\dots, 0)^T$，其中 $1$ 在第 $i$ 个元素。
$$
\mathbf x^TA\mathbf x=a_{ii}
$$
由正定矩阵定义可得，$a_{ii}>0,i\in\{1, \dots,n\}$。

(2)

对于矩阵的 $ij$ 元素，$i,j\in\{2, \dots,n\}$：
$$
\begin{align*}
a^{(2)}_{ij}&=a_{ij}-\frac{a_{i1}}{a_{11}}a_{1j}=\frac{a_{11}a_{ij}-a_{i1}a_{1j}}{a_{11}}\\
a^{(2)}_{ji}&=a_{ji}-\frac{a_{j1}}{a_{11}}a_{1i}=\frac{a_{11}a_{ji}-a_{j1}a_{1i}}{a_{11}}
\end{align*}
$$
由于 $A$ 是对称矩阵，则 $a^{(2)}_{ij}=a^{(2)}_{ji},i,j\in\{2, \dots,n\}$，因此 $A_2$ 是对称矩阵。

由于 $A$ 是正定矩阵，则 $A$ 的顺序主子式 $\Delta_k>0$，且约化后仍成立。取 $A_2$​ 的顺序余子式，则 $A_2$ 的顺序余子式 $\Delta^{(2)}_k$：
$$
\Delta^{(2)}_k=\frac{\Delta_{k+1}}{a_{11}}>0\quad k\in\{1, \dots,n-1\}
$$
因此 $A_2$ 也是正定矩阵。

(3)

此时 $A_2$ 的对角线的元素为：
$$
a^{(2)}_{ii}=a_{ii}-\frac{a_{i1}}{a_{11}}a_{1i}=a_{ii}-\frac{a_{i1}^2}{a_{11}}\le a_{ii}\quad i\in\{2, \dots,n\}
$$
(4)

假设 $a_{ij}$ 是 $A$ 中绝对值最大的元素，由于 $A$ 是对称正定矩阵，则其包含 $a_{ij}$ 的主子式必须大于 0。而：
$$
\begin{vmatrix}a_{ii}&a_{ij}\\a_{ji}&a_{jj}\end{vmatrix}=a_{ii}a_{jj}-a^2_{ij}>0\Rightarrow a^2_{ij}<a_{ii}a_{jj}
$$
这与 $a_{ij}$ 是 $A$ 中绝对值最大的元素的假设矛盾，因此最大的元素必在对角线上。

(5)

由于对称正定矩阵的绝对值最大的元素必在对角线上，则要证
$$
\max_{2\le i,j\le n}|a^{(2)}_{ij}|\le\max_{2\le i,j\le n}|a_{ij}|
$$
由 (1), (2), (4) 可得即证：
$$
\begin{align*}
\max_{2\le i\le n}a^{(2)}_{ii}&\le\max_{2\le i\le n}a_{ii}\\
\max_{2\le i\le n}(a_{ii}-\frac{a_{i1}^2}{a_{11}})&\le\max_{2\le i\le n}a_{ii}\\
\end{align*}
$$
显然这是恒成立的，得证。

## Question 2

> 设 $A\in\mathbb R^{n\times n}$，$A$ 经过 Gauss 消元法一步以后：
> $$
> A^{(1)}=\begin{bmatrix}
> a_{11}     &*\\
> \mathbf 0^T&A_2
> \end{bmatrix}\quad A_2\in\mathbb R^{(n-1)\times(n-1)}\quad a^{(1)}_{ij}&=a_{ij}-\frac{a_{i1}}{a_{11}}a_{1i}
> $$
> (1) 证明若 $A$ 严格主对角占优，则 $A^{(1)}$ 也严格主对角占优
>
> (2) 若 $n=4$，
> $$
> A=\begin{bmatrix}
> 4&-1&&\\
> -1&4&-1&\\
> &-1&4&-1\\
> &&-1&4
> \end{bmatrix}
> $$
> 求 $A^{(1)}$。

(1)

要证 $A^{(1)}$ 严格主对角占优，即证：
$$
\left|a_{ii}-\frac{a_{i1}}{a_{11}}a_{1i}\right|>\sum_{j=2,j\ne i}^{n}\left|a_{ij}-\frac{a_{i1}}{a_{11}}a_{1j}\right|\quad i,j\in\{2,\dots,n\}\\
$$
而：
$$
\begin{align*}
\left|a_{ii}-\frac{a_{i1}}{a_{11}}a_{1i}\right|&\ge|a_{ii}|-\frac{|a_{i1}|}{|a_{11}|}|a_{1i}|\\
\left|a_{ij}-\frac{a_{i1}}{a_{11}}a_{1j}\right|&\le|a_{ij}|+\frac{|a_{i1}|}{|a_{11}|}|a_{1j}|\quad i,j\in\{2,\dots,n\}\\
\end{align*}
$$
带入，ci'shi即证：
$$
\begin{align*}
|a_{ii}|&>\sum_{j=2,j\ne i}^{n}|a_{ij}|+\sum_{j=2}^{n}\frac{|a_{i1}|}{|a_{11}|}|a_{1j}|\\
|a_{ii}|&>\sum_{j=2,j\ne i}^{n}|a_{ij}|+|a_{i1}|\sum_{j=2}^{n}\frac{|a_{1j}|}{|a_{11}|}\quad i,j\in\{2,\dots,n\}\\
\end{align*}
$$
又因为 $A$ 是严格主对角占优，所以：
$$
\sum_{j=2}^{n}\frac{|a_{1j}|}{|a_{11}|}<1
$$
因此：
$$
\begin{align*}
\sum_{j=2,j\ne i}^{n}|a_{ij}|+|a_{i1}|\sum_{j=2}^{n}\frac{|a_{1j}|}{|a_{11}|}&<\sum_{j=2,j\ne i}^{n}|a_{ij}|+|a_{i1}|\\
&=\sum_{j=1,j\ne i}^{n}|a_{ij}|\\
&<|a_{ii}|
\end{align*}
$$
恒成立，因此 $A^{(1)}$ 也严格主对角占优。

(2)
$$
A^{(1)}=\begin{bmatrix}
4&-1&&\\
0&\frac{15}{4}&-1&\\
&-1&4&-1\\
&&-1&4
\end{bmatrix}
$$

## Question 3

> 用列主元法解线性代数方程：
> $$
> \begin{cases}
> x_1+2x_2+3x_3=1\\
> 2x_1+4x_2+5x_3=1\\
> 3x_1+5x_2+6x_3=1\\
> \end{cases}
> $$
> 若不选主元，Gauss 消元法能否解此方程组？说明选主元的原因。

将线性方程组转化为增广矩阵并进行列主元法解线性方程组：
$$
\begin{bmatrix}
1&2&3&1\\
2&4&5&1\\
3&5&6&1\\
\end{bmatrix}
\to
\begin{bmatrix}
1&2&3&1\\
0&0&-1&-1\\
0&-1&-3&-2\\
\end{bmatrix}
\to
\begin{bmatrix}
1&2&3&1\\
0&-1&-3&-2\\
0&0&-1&-1\\
\end{bmatrix}
\to
\begin{cases}
x_1=0\\
x_2=-1\\
x_3=1\\
\end{cases}
$$
不选主元，Gauss 消元法不能解此方程组，故选主元可以拓展 Gauss 消元法的应用，即选取的主元为 0 或很小的数的时候，通过交换两行，使得 Gauss 消元法得以进行下去。

## Question 4

> 对矩阵 
> $$
> A=\begin{bmatrix}
> 2&-1&1\\
> -1&-2&3\\
> 1&3&1\\
> \end{bmatrix}
> $$
> 用 LU 分解，并解方程 $Ax=b$，其中
> $$
> x=\begin{bmatrix}
> x_1\\x_2\\x_3
> \end{bmatrix}
> \quad
> b=\begin{bmatrix}
> 4\\5\\6
> \end{bmatrix}
> $$

$$
\begin{bmatrix}
2&-1&1\\
-1&-2&3\\
1&3&1\\
\end{bmatrix}
\to
\begin{bmatrix}
2&-1&1\\
-0.5&-2&3\\
0.5&3&1\\
\end{bmatrix}
\to
\begin{bmatrix}
2&-1&1\\
-0.5&-2.5&3.5\\
0.5&-1.4&1\\
\end{bmatrix}
\to
\begin{bmatrix}
2&-1&1\\
-0.5&-2.5&3.5\\
0.5&-1.4&5.4\\
\end{bmatrix}
$$

因此：
$$
L=
\begin{bmatrix}
1&0&0\\
-0.5&1&0\\
0.5&-1.4&1\\
\end{bmatrix}
\quad
U=
\begin{bmatrix}
2&-1&1\\
0&-2.5&3.5\\
0&0&5.4\\
\end{bmatrix}
$$

$$
Ax=b\Rightarrow LUx=b\Rightarrow Lc=b,Ux=c
$$

先求出 $c$
$$
c=\begin{bmatrix}
4\\7\\13.8
\end{bmatrix}
$$
再求出 $x$
$$
x=\begin{bmatrix}
10/9\\7/9\\23/9
\end{bmatrix}
$$

## Question 5

> 对矩阵 
> $$
> A=\begin{bmatrix}
> 1&-1&2\\
> -1&5&2\\
> 2&2&17\\
> \end{bmatrix}
> $$
> 作平方根分解，规定 $G$ 的对角元素全为正。

$$
\begin{bmatrix}
1&-1&2\\
-1&5&2\\
2&2&17\\
\end{bmatrix}
\to
\begin{bmatrix}
1&0&0\\
-1&5&0\\
2&2&17\\
\end{bmatrix}
\to
\begin{bmatrix}
1&0&0\\
-1&2&0\\
2&2&3\\
\end{bmatrix}
$$

因此：
$$
G=\begin{bmatrix}
1&0&0\\
-1&2&0\\
2&2&3\\
\end{bmatrix}
$$
