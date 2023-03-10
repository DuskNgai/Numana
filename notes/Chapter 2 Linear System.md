# Chapter 2 方程组 Linear System

## 2.1 Gauss 消元法 Gauss Elimination

### 2.1.1 上三角矩阵回代

设 $a_{ii}\ne0$：
$$
\begin{cases}
x_n=b_n/a_{nn}\\
x_i=(b_i-\sum_{j=i+1}^{n}a_{ij}x_{j})/a_{ii}
\end{cases}
$$

### 2.1.3 Gauss 消元可行性

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 2.1 顺序主子式
</blockquote>

$$
\Delta_1=a_{11},\Delta_2=\begin{vmatrix}a_{11}&a_{12}\\a_{21}&a_{22}\end{vmatrix},\dots,\Delta_n=\det(A)
$$


<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 2.1
</blockquote>

如果矩阵的顺序主子式都不为 0，则 Gauss 消元可以进行到底。如果某个顺序主子式 $\Delta_k$ 不为 0，则
$$
a^{(1)}_{11}=\Delta_{1}\\
a^{(k)}_{kk}=\Delta_{k}/\Delta_{k-1}
$$

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 2.2 严格对角占优
</blockquote>

如果矩阵 $A$ 满足：
$$
|a_{ii}|>\sum_{j\ne i}^{n}|a_{ij}|\quad i\in\{1, 2,\dots, n\}
$$
则称矩阵 $A$ 为严格对角占优。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 2.2
</blockquote>

如果线性方程组的系数矩阵 $A$ 为严格对角占优阵，经过一步的 Gauss 消元后具有：
$$
\begin{bmatrix}
a_{11}     &*\\
\mathbf 0^T&A_2
\end{bmatrix}
$$
的形式，则 $A_2$ 仍为严格对角占优阵，且 Gauss 消元法能进行到底。

### 2.1.4 Gauss 消元法的矩阵分析

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 2.3 LU 分解
</blockquote>


由于对矩阵A施行初等行变换相当于左乘初等矩阵，消元可用左乘矩阵表示。
$$
L_kA^{(k)}=A^{(k+1)}\quad L_kb^{(k)}=b^{(k+1)}\\
L_k=\begin{bmatrix}
1&      &          & &      & \\
 &\ddots&          & &      & \\
 &      &1         & &      & \\
 &      &-m_{k+1,k}&1&      & \\
 &      &\vdots    & &\ddots& \\
 &      &-m_{n,k}  & &      &1\\
\end{bmatrix}
$$
显然 $L_k$ 是下三角矩阵。而 Gauss 消元法最终会得到一个上三角矩阵 $U$，则：
$$
\begin{align*}
L_{n-1}\dots L_{1}A&=U\\
A&=L^{-1}_{1}\dots L^{-1}_{n-1}U=LU
\end{align*}
$$
其中：
$$
L=\begin{bmatrix}
1     &      &      &           & \\
m_{21}&1     &      &          & \\
m_{31}&m_{32}&1     &          & \\
\vdots&\vdots&\ddots&\ddots    & \\
m_{n1}&m_{n2}&\cdots&m_{n(n-1)}&1\\
\end{bmatrix}
$$

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 2.3
</blockquote>

如果矩阵 $A$ 的各阶顺序主子式非零，则 $A$ 一定可作 $LU$ 分解，其中 $L$ 为单位下三角阵，$U$ 为上三角阵，且**分解唯一**。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    推论
</blockquote>

$$
A=LU=LDU'
$$

$U'$ 是单位上三角阵，$D$ 是对角矩阵。

### 2.1.5 列主元 Gauss 消元法

把一个特别需要的元素放到对角线的位置上选取主元。把列绝对值最大的行交换一下。

设已经进行了 $(k-1)$ 步 Gauss 消元，记 $|a^{(k)}_{i_0k}|=\max_{k\le i\le n}|a^{(k)}_{ik}|$，交换 $i_0$ 和 $k$ 行，再消元。

将若干个初等排列阵的乘积称为排列阵

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 2.4
</blockquote>

设 $A$ 为非奇异阵，则存在排列阵 $P$，单位下三角阵 $L$ 和上三角阵 $U$，使得 $PA=LU$。

## 2.2 三角分解法

### 2.2.1 直接三角分解法

$$
L=\begin{bmatrix}
1     &      &      &      & \\
l_{21}&1     &      &      & \\
l_{31}&l_{32}&1     &      & \\
\vdots&\vdots&\ddots&\ddots& \\
l_{n1}&l_{n2}&l_{n3}&\ddots&1\\
\end{bmatrix}
\quad
U=\begin{bmatrix}
u_{11}&u_{12}&u_{13}&\cdots&u_{1n}\\
      &u_{22}&u_{23}&\cdots&u_{2n}\\
      &      &u_{33}&\cdots&u_{3n}\\
      &      &     &\ddots&\vdots\\
      &      &      &     &u_{nn}\\
\end{bmatrix}
$$

> $i\le j$，即上三角部分
> $$
> a_{ij}=\sum_{k=1}^{i}l_{ik}u_{kj}=\sum_{k=1}^{i-1}l_{ik}u_{kj}+u_{ij}
> $$
> $i>j$，即下三角部分
> $$
> a_{ij}=\sum_{k=1}^{j}l_{ik}u_{kj}=\sum_{k=1}^{j-1}l_{ik}u_{kj}+l_{ij}u_{jj}
> $$
> 因此先得出 $U$ 的一行：
> $$
> u_{ij}=a_{ij}-\sum_{k=1}^{i-1}l_{ik}u_{kj}\quad i\le j
> $$
> 再得出 $L$ 的一列：
> $$
> l_{ij}=(a_{ij}-\sum_{k=1}^{j-1}l_{ik}u_{kj})/u_{jj}\quad i>j
> $$

### 2.2.2 正定方程的 Cholesky 分解

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 2.5
</blockquote>

设 $A$ 为**对称正定矩阵**，则存在三角分解 $A=RR^T$，其中 $R$ 是非奇异下三角阵，且当限定 $G$ 的对角线元素为正时，这种分解惟一。

> $$
> A=LDL^T=(LD^{1/2})(D^{1/2}L^T)=(LD^{1/2})(LD^{1/2})^T=RR^T\\
> $$

设
$$
R=\begin{bmatrix}
r_{11}&      &      &      \\
r_{21}&r_{22}&      &      \\
\vdots&\vdots&\ddots&      \\
r_{n1}&r_{n2}&\cdots&r_{nn}\\
\end{bmatrix}
\quad
R^T=\begin{bmatrix}
r_{11}&r_{21}&\cdots&r_{n1}\\
      &r_{22}&\cdots&r_{n2}\\
      &      &\ddots&\vdots\\
      &      &      &r_{nn}\\
\end{bmatrix}
$$
> $i=j$，即对角部分
> $$
> a_{ii}=\sum_{k=1}^{i}r^2_{ik}=\sum_{k=1}^{i-1}r^2_{ik}+r^2_{ii}
> $$
> $i<j$，即上三角部分
> $$
> a_{ij}=\sum_{k=1}^{i}r_{ik}r_{kj}=\sum_{k=1}^{i-1}r_{ik}r_{kj}+r_{ii}r_{ij}
> $$
> 因此先得出 $R$ 的一个对角线元素：
> $$
> r^2_{ii}=\sqrt{a_{ii}-\sum_{k=1}^{i-1}r^2_{ik}}
> $$
> 再得出 $R$ 的其他元素：
> $$
> r_{ij}=(a_{ij}-\sum_{k=1}^{i-1}r_{ik}r_{kj})/r_{ii}
> $$

### 2.2.3 三对角方程的追赶法

三对角方程 $Ax=f$：
$$
A=\begin{bmatrix}
b_{1}&c_{1} &      &      &       \\
a_{2}&b_{2} &c_{2} &      &       \\
     &\ddots&\ddots&\ddots&       \\
     &      &\ddots&\ddots&c_{n-1}\\
     &      &      &a_{n} &b_{n}  \\
\end{bmatrix}
$$
且 (1) $|b_1|>|c_1|$；(2) $|b_i|\ge|a_i|+|c_i|$；(3) $|b_n|>|a_n|$ 且均非零。

可对其做 Court 分解。

## 2.3 误差分析

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 2.4 条件数
</blockquote>

设 $A$ 为非奇异阵，称 $\mathrm{cond}(A)=\|A^{-1}\|\|A\|$ 为矩阵的条件数。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 2.6
</blockquote>

设 $A$ 为非奇异阵，方程 $Ax=b\ne0$，且
$$
\begin{align*}
A(x+\delta x)&=b+\delta b\\
Ax+A\delta x&=b+\delta b\\
A\delta x&=\delta b\\
\delta x&=A^{-1}\delta b\\
\end{align*}
$$
绝对误差为：
$$
\|\delta x\|=\|A^{-1}\delta b\|\le\|A^{-1}\|\|\delta b\|
$$
相对误差为：
$$
\frac{\|\delta x\|}{\|x\|}\le\frac{\|A^{-1}\|\|\delta b\|}{\|Ax\|/\|A\|}=\|A\|\|A^{-1}\|\frac{\|\delta b\|}{\|b\|}
$$

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 2.7
</blockquote>

设 $A$ 为非奇异阵，方程 $Ax=b\ne0$，且
$$
\begin{align*}
(A+\delta A)(x+\delta x)&=b\\
Ax+\delta A x+A\delta x+\delta A\delta x&=b\\
\delta A x+A\delta x+\delta A\delta x&=0\\
\delta x&=-A^{-1}\delta Ax-A^{-1}\delta A\delta x\\
\end{align*}
$$
相对误差为：
$$
\begin{align*}
\|\delta x\|&=\|A^{-1}\delta Ax+A^{-1}\delta A\delta x\|\\
&\le\|A^{-1}\delta Ax\|+\|A^{-1}\delta A\delta x\|\\
&\le\|A^{-1}\|\|\delta A\|\|x\|+\|A^{-1}\|\|\delta A\|\|\delta x\|\\
&=\frac{\|A^{-1}\|\|\delta A\|\|x\|}{1-\|A^{-1}\|\|\delta A\|}
\end{align*}
$$
即：
$$
\frac{\|\delta x\|}{\|x\|}\le\frac{\|A^{-1}\|\|\delta A\|}{1-\|A^{-1}\|\|\delta A\|}=\frac{\|A^{-1}\|\|A\|(\|\delta A\|/\|A\|)}{1-\|A^{-1}\|\|A\|(\|\delta A\|/\|A\|)}
$$

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 2.8
</blockquote>

设 $A$ 为非奇异阵，方程 $Ax=b\ne0$，且
$$
\begin{align*}
(A+\delta A)(x+\delta x)&=b+\delta b
\end{align*}
$$
误差为：
$$
\frac{\|\delta x\|}{\|x\|}\le\frac{\mathrm{cond}(A)}{1-\mathrm{cond}(A)(\|\delta A\|/\|A\|)}\left(\frac{\|\delta A\|}{\|A\|}+\frac{\|\delta b\|}{\|b\|}\right)
$$

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 2.9
</blockquote>

设 $A$ 为非奇异阵，$x^*$ 是方程 $Ax=b\ne0$ 的精确解，$\bar x$ 是近似解：
$$
\frac{\|\bar x-x^*\|}{\|x\|}\le\mathrm{cond}(A)\frac{\|b-A\bar x\|}{\|b\|}
$$

## 2.4 迭代法

### 2.4.1 迭代法的构造

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 2.5 迭代法
</blockquote>

设 $A$ 非奇异，如果方程组 $Ax=b$ 与 $x=Bx+f$ 具有相同的解，则二者等价。
$$
\begin{align*}
Ax&=b\\
(M-N)x&=b\\
Mx&=Nx+b\\
x&=M^{-1}Nx+M^{-1}b
\end{align*}
$$
或
$$
\begin{align*}
Ax&=b\\
x&=x+\lambda(b-Ax)\\
x&=(I-\lambda A)x+\lambda b
\end{align*}
$$

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 2.10
</blockquote>

若矩阵 $B$ 满足 $\|B\|<1$，则 $x=Bx+f$ 存在唯一解。

### 2.4.2 迭代格式的收敛性

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 2.6 迭代格式收敛
</blockquote>

若迭代序列 $x^{(k)}$ 满足 $\lim_{k\to\infty}x^{(k)}=x^*$，称迭代格式收敛。$B$ 为迭代矩阵，$x^{(0)}$ 为初始迭代向量，$x^{(k)}$ 为第 $k$ 步近似值，$e^{(k)}=x^*-x^{(k)}$ 为第 $k$ 步误差向量。
$$
e^{(k)}=x^*-x^{(k)}=Bx^*-f-Bx^{(k)}+f=B^{k}(x^*-x^{(0)})=B^{k}e^{(0)}
$$
范数可得：
$$
\|e^{(k)}\|\le\|B^k\|\|e^{(0)}\|\le\|B\|^k\|e^{(0)}\|
$$

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 2.11
</blockquote>

- 迭代格式收敛。
- $\rho(B)<1$ ($\|B\|\le\rho(B)+\epsilon$)
- 至少存在一种算子范数使得 $\|B\|<1$ ($\rho(B)\le\|B\|<1$)

<blockquote style="border-left: 5px solid #bb4545; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(188, 70, 70, 0.1)">
    例题
</blockquote>

$$
B=\begin{bmatrix}0&2\\3&0\end{bmatrix}
$$

> $\lambda=\pm\sqrt{6}>1$
>
> 因此发散。

<blockquote style="border-left: 5px solid #bb4545; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(188, 70, 70, 0.1)">
    例题
</blockquote>

$$
A=\begin{bmatrix}2&1\\1&2\end{bmatrix}\quad b=\begin{bmatrix}1\\2\end{bmatrix}\\
x^{(k+1)}=x^{(k)}+\omega(Ax^{(k)}+b)
$$

> $$
> x^{(k+1)}=(\omega A+I)x^{(k)}+\omega b
> $$
>
> 因此：
> $$
> \det((\lambda-1)I-\omega A)=(\lambda-1-2\omega)(\lambda-1-2\omega)-\omega^2\\
> \lambda=1+\omega,1+3\omega\\
> -1\le1+\omega,1+3\omega\le1\\
> -2/3\le\omega\le0
> $$

<blockquote style="border-left: 5px solid #bb4545; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(188, 70, 70, 0.1)">
    例题
</blockquote>

$$
x^{(k+1)}_i=x^{(k)}_i+\frac{w}{a_{ii}}(b_{i}-\sum_{j=1}^na_{ij}x^{(k)}_j)
$$

(1) 试写出其矩阵形式迭代格式及迭代矩阵

(2) 若 $A$ 是严格对角占优阵，且 $w=1/2$，则迭代收敛.

> (1)
> $$
> x^{(k+1)}=x^{(k)}+w\ \mathrm{diag}(A)^{-1}(b-Ax^{(k)})\\
> x^{(k+1)}=(I-w\ \mathrm{diag}(A)^{-1}A)x^{(k)}+w\ \mathrm{diag}(A)^{-1}b
> $$
> (2)
> $$
> I-\frac{1}{2}\mathrm{diag}(A)^{-1}A=\begin{bmatrix}
> \frac{1}{2}&-\frac{1}{2}a_{12}/a_{11}&\dots\\
> -\frac{1}{2}a_{21}/a_{11}&\frac{1}{2}&\dots\\
> \vdots&\vdots&\ddots\\
> \end{bmatrix}
> $$
> 取范数：
> $$
> \|B\|_{\infty}=\max_{i}\sum_{j=1}^{n}|\frac{a_{ij}}{2a_{ii}}|<1
> $$

### 2.4.3 迭代格式的误差估计与收敛速度

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 2.12
</blockquote>

$$
\|x^{(k)}-x^*\|\le\frac{\|B\|}{1-\|B\|}\|x^{(k)}-x^{(k-1)}\|\\
\|x^{(k)}-x^*\|\le\frac{\|B\|^k}{1-\|B\|}\|x^{(1)}-x^{(0)}\|\\
$$

第一项为事后误差估计，可以用很小的数 $\epsilon$ 来控制迭代终止，其前提条件为 $\|B\|<1$。

第二项为先验估计。令 $\|B\|=q$，则：
$$
\|x^{(k)}-x^{(0)}\|<\frac{q^k}{1-q}\|x^{(1)}-x^{(0)}\|<\epsilon
$$
则迭代次数为：
$$
k\ge\ln(\frac{\epsilon(1-\|B\|)}{\|x^{(1)}-x^{(0)}\|})/\ln(\|B\|)
$$

> $$
> \begin{align*}
> \|x^{(k+p)}-x^{(k)}\|&\le\|x^{(k+p)}-x^{(k+p-1)}\|+\dots+\|x^{(k+1)}-x^{(k)}\|\\
> &\le(\|B\|^{p-1}+\dots+\|B\|^0)\|x^{(k+1)}-x^{(k)}\|\\
> &=\frac{\|B\|^{p}}{1-\|B\|}\|x^{(k+1)}-x^{(k)}\|\\
> \lim_{p\to\infty}\|x^{(k+p)}-x^{(k)}\|&=\lim_{p\to\infty}\frac{\|B\|^{p}}{1-\|B\|}\|x^{(k+1)}-x^{(k)}\|\\
> \|x^{(k)}-x^*\|&=\frac{1}{1-\|B\|}\|x^{(k+1)}-x^{(k)}\|\\
> &\le\frac{\|B\|}{1-\|B\|}\|x^{(k)}-x^{(k-1)}\|\\
> \end{align*}
> $$

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 2.7 迭代平均收敛速度
</blockquote>

$$
R_k(B)=k\ln\|B^k\|=-\ln\|B^k\|^{1/k}
$$

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 2.8 渐近收敛速度
</blockquote>

$$
R(B)=-\ln\rho(B)
$$

<blockquote style="border-left: 5px solid #bb4545; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(188, 70, 70, 0.1)">
    例题
</blockquote>

$$
B=-\frac{1}{2}\begin{bmatrix}0&1\\1&0\end{bmatrix}\quad f=\frac{1}{2}\begin{bmatrix}1\\2\end{bmatrix}\\
x^{(k+1)}=Bx^{(k)}+f
$$

(1) $x^{(0)}=(1,1)^T$，求 $x^{(k)}$

(2) 使误差 $\|x^{(k)}-x^*\|_{\infty}<1/2^{10}$，估计迭代次数

(3) 求上述迭代格式的收敛速度

(4) 写出上述迭代格式对应的方程组

(5) 写出误差 $\|x^{(k)}-x^*\|_{\infty}$ 的精确表达式

> (1)
> $$
> \begin{bmatrix}1\\1\end{bmatrix}\to\begin{bmatrix}0\\1/2\end{bmatrix}\to\begin{bmatrix}1/4\\1\end{bmatrix}\to\begin{bmatrix}0\\7/8\end{bmatrix}\to\cdots\begin{bmatrix}0\\1\end{bmatrix}
> $$
> $$
> \begin{bmatrix}0\\1-1/2^k\end{bmatrix}(k\text{ is odd})\quad
> \begin{bmatrix}1/2^k\\1\end{bmatrix}(k\text{ is even})
> $$
>
> (2)
> $$
> \begin{align*}
> k&\ge\ln(\frac{\epsilon(1-\|B\|)}{\|x^{(1)}-x^{(0)}\|})/\ln(\|B\|)\\
> &=-\ln(2^{-11})/\ln(2)\\
> &=11
> \end{align*}
> $$
> (3)
> $$
> \rho(B)=\max|\lambda(B)|=\frac{1}{2}\\
> R(B)=-\ln\rho(B)=\ln2
> $$
> (4)
> $$
> \begin{cases}
> x^{(k+1)}_1=-\frac{1}{2}x^{(k)}_2+\frac{1}{2}\\
> x^{(k+1)}_2=-\frac{1}{2}x^{(k)}_1+1\\
> \end{cases}
> $$
> (5)
>
> $$
> \|x^{(k)}-x^*\|_{\infty}=1/2^k
> $$

## 2.5 Jacobi 迭代和 Gauss-Seidel 迭代

### 2.5.1 Jacobi 迭代

$$
A=D-L-U\\
Ax=b\Rightarrow (D-L-U)x=b\Rightarrow Dx=(L+U)x+b\\
x=D^{-1}(L+U)x+D^{-1}b\\
x=D^{-1}(D-A)x+D^{-1}b\\
x=(I-D^{-1}A)x+D^{-1}b
$$
记 $J=(I-D^{-1}A)$。

其分量形式为：
$$
\begin{align*}
x^{(k+1)}_i&=x^{(k)}_i-\sum_{j=1}^{n}\frac{a_{ij}}{a_{ii}}x^{(k)}_j+\frac{b_i}{a_{ii}}\\
&=\frac{1}{a_{ii}}\left(b_i-\sum_{j=1,j\ne i}^{n}a_{ij}x^{(k)}_j\right)
\end{align*}
$$
### 2.5.2 Gauss-Seidel 迭代

如果迭代收敛，则后一步的近似值比前一步更精确，因此，可将 Jacobi 迭代的分量形式改写为：
$$
\begin{align*}
x^{(k+1)}_i&=\frac{1}{a_{ii}}\left(b_i-\underbrace{\sum_{j=1}^{i-1}a_{ij}x^{(k+1)}_j}_{\text{lower triangle}}-\underbrace{\sum_{j=i+1}^{n}a_{ij}x^{(k)}_j}_{\text{upper triangle}}\right)
\end{align*}
$$
矩阵形式为：
$$
\begin{align*}
x^{(k+1)}_i&=\frac{1}{a_{ii}}\left(b_i-\sum_{j=1}^{i-1}a_{ij}x^{(k+1)}_j-\sum_{j=i+1}^{n}a_{ij}x^{(k)}_j\right)\\
Dx^{(k+1)}&=\left(b+Lx^{(k+1)}+Ux^{(k)}\right)\\
(D-L)x^{(k+1)}&=Ux^{(k)}+b\\
x^{(k+1)}&=(D-L)^{-1}Ux^{(k)}+(D-L)^{-1}b\\
\end{align*}
$$
<blockquote style="border-left: 5px solid #bb4545; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(188, 70, 70, 0.1)">
    例题
</blockquote>

设 
$$
A=\begin{bmatrix}
2&1\\
1&2
\end{bmatrix}\quad
b=\begin{bmatrix}
1\\2
\end{bmatrix}
$$
(1) 写出求解方程 $Ax=b$ 的 Jacobi 迭代和 Gauss-Seidel 迭代格式

(2) 判断 Jacobi 迭代和 Gauss-Seidel 迭代是否收敛

> (1)
> $$
> x^{(k+1)}=\begin{bmatrix}0&-0.5\\-0.5&0\end{bmatrix}x^{(k)}+\begin{bmatrix}1/2\\1\end{bmatrix}\\
> x^{(k+1)}=\begin{bmatrix}0&-1/2\\0&1/4\end{bmatrix}x^{(k)}+\begin{bmatrix}1/2\\3/4\end{bmatrix}
> $$
>
> (2)
> $$
> \rho(J)=0.5\quad\rho(G)=0.5
> $$

### 2.5.3 Jacobi 迭代和 Gauss-Seidel 迭代的收敛性

Jacobi 迭代矩阵满足：
$$
\det(\lambda I-(I-D^{-1}A))=\frac{1}{\det(D)}\det(\lambda D-L-U)
$$
Gauss-Seidel 迭代矩阵满足：
$$\begin{align*}
\det(\lambda I-(D-L)^{-1}U)&=\frac{1}{\det(D-L)}\det(\lambda(D-L)-U)\\
&=\frac{1}{\det(D)}\det(\lambda(D-L)-U)\end{align*}
$$
所以可以很容易得出谱半径。

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 2.9 
</blockquote>

若 $A\in\mathbb R^{n\times n}$ 不能找到排列矩阵 $P$，使得
$$
P^TAP=\begin{bmatrix}A_{11}&A_{12}\\0&A_{22}\end{bmatrix}
$$
称 $A$ 为不可约的。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 2.13
</blockquote>

设 $A$ 为严格对角占优阵或不可约弱对角占优阵，则解方程组 $Ax=b$ 的 Jacobi 迭代和 Gauss-Seidel 迭代均收敛。

> 仅对 $A$ 为严格对角占优的 Jacobi 迭代证明：
> $$
> \|J\|_{\infty}=\max\left(\frac{1}{|a_{ii}|}\sum_{j\ne i}|a_{ij}|\right)<1
> $$

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 2.14
</blockquote>

若 $A$ 的对角元素全为正数，而非对角元素非正，则以下四个命题有且仅有一个成立：
$$
\rho(G)=\rho(J)=0\qquad0<\rho(G)<\rho(J)<1\\
\rho(G)=\rho(J)=1\qquad\rho(G)>\rho(J)>1
$$
## 2.6 超松弛（SOR）迭代法

### 2.6.1

原方程等价于 $\omega Ax=\omega b$，可等价于化为
$$
\begin{align*}
\omega(D-L-U)x&=\omega b\\
\omega(D-L)x&=\omega Ux+\omega b\\
(D-\omega L)x&=[(1-\omega)D+\omega U]x+\omega b\\
\end{align*}
$$
若 $D$ 可逆，则 $D-\omega L$ 也可逆：
$$
x=(D-\omega L)^{-1}[(1-\omega)D+\omega U]x+\omega(D-\omega L)^{-1}b\\
$$
$L_{\omega}=(D-\omega L)^{-1}[(1-\omega)D+\omega U]$ 为 SOR 迭代矩阵。$\omega=1$ 即为 Gauss-Seidel 迭代，$\omega>1$ 为超松弛迭代。

### 2.6.2 SOR 迭代的收敛性

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 2.15
</blockquote>
$$
\rho(L_{\omega})\ge|\omega-1|
$$

> 记 $L_{\omega}$ 的特征值为 $\lambda_i$
> $$
> \begin{align*}
> \prod_{i=1}^n\lambda_i&=\det(L_{\omega})=\det(D-\omega L)^{-1}\det((1-\omega)D+\omega U)\\
> &=\prod_{i=1}^n\frac{1}{a_{ii}}\prod_{i=1}^n(1-\omega){a_{ii}}\\
> &=(1-\omega)^n
> \end{align*}
> $$
> 若收敛则 $|1-\omega|<1$，即 $0<\omega<2$。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 2.16
</blockquote>

设 $A$ 为对称正定，且 $0<\omega<2$，则 SOR 迭代收敛。

> 设 $\lambda,x$ 为 $L_{\omega}$ 的特征值和特征向量：
> $$
> [(1-\omega)D+\omega U]x=\lambda(D-\omega L)x
> $$
> 记 $p=x^TDx>0,a+ib=x^TLx$，则 $a-ib=x^TUx,p-2a=x^TAx$，两边与 $x$ 做内积：
> $$
> \lambda=\frac{(1-\omega)p+\omega a-i\omega b}{p-\omega a-i\omega b}\\
> |\lambda|^2=\frac{((1-\omega)p+\omega a)^2+(\omega b)^2}{(p-\omega a)^2+(\omega b)^2}\\
> ((1-\omega)p+\omega a)^2-(p-\omega a)^2=p\omega(2-\omega)(2a-p)<1
> $$

