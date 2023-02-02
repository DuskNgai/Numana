# Homework 1

## Question 1

> 指出下列近似数的绝对误差限、相对误差限和有效数字位数：
> $$
> 23000, 0.00230, 2300.00, 2.30\times10^{4}
> $$

(1)
$$
\begin{gather*}
23000=2.3000\times10^{4}\Rightarrow m=4,n=5\\
\epsilon\le0.5\times10^{4-5+1}=0.5\\
\epsilon_r\le\frac{1}{2\times2}\times10^{-5+1}=2.5\times10^{-5}
\end{gather*}
$$
$23000$ 具有 $5$ 位有效数字，绝对误差限 $\epsilon$ 为 $0.5$，相对误差限 $\epsilon_r$ 为 $2.5\times10^{-5}$。

(2)
$$
0.00230=2.30\times10^{-3}\Rightarrow m=-3,n=3\\
\epsilon\le0.5\times10^{-3-3+1}=0.5\times10^{-5}\\
\epsilon_r\le\frac{1}{2\times2}\times10^{-3+1}=2.5\times10^{-3}
$$
$0.00230$ 具有 $3$ 位有效数字，绝对误差限 $\epsilon$ 为 $0.5\times10^{-5}$，相对误差限 $\epsilon_r$ 为 $2.5\times10^{-3}$。

(3)
$$
2300.00=2.30000\times10^{3}\Rightarrow m=3,n=6\\
\epsilon\le0.5\times10^{3-6+1}=0.5\times10^{-2}\\
\epsilon_r\le\frac{1}{2\times2}\times10^{-6+1}=2.5\times10^{-6}
$$
$2300.00$ 具有 $6$ 位有效数字，绝对误差限 $\epsilon$ 为 $0.5\times10^{-2}$，相对误差限 $\epsilon_r$ 为 $2.5\times10^{-6}$。

(4)
$$
2.30\times10^{4}\Rightarrow m=4,n=3\\
\epsilon\le0.5\times10^{4-3+1}=0.5\times10^{2}\\
\epsilon_r\le\frac{1}{2\times2}\times10^{-3+1}=2.5\times10^{-3}
$$
$2.30\times10^{4}$ 具有 $3$ 位有效数字，绝对误差限 $\epsilon$ 为 $0.5\times10^{2}$，相对误差限 $\epsilon_r$ 为 $2.5\times10^{-3}$。

## Question 2

> 用某数值方法计算 $I^*=\cos(\pi/3)$ 所得值 $I=0.500012$ 为，求 $I$ 的有效数字位数。

$$
\begin{align*}
I^*=\cos\frac{\pi}{3}&=0.500000\\
I&=0.500012=5.00012\times10^{-1}\Rightarrow m=-1\\
|I^*-I|&=0.000012\le0.5\times10^{-4}\Rightarrow n=4
\end{align*}
$$

$I$ 的有效数字位数为 $4$。

## Question 3

> 计算球体积要使相对误差限为 $1\%$，问度量半径 $R$ 允许的相对误差是多少？

$$
V=\frac{4}{3}\pi R^3
$$

$$
\begin{align*}
\epsilon_r(V)&=\left|\frac{dV}{dR}\right|\frac{\epsilon(R)}{|V|}\\
&=3\epsilon(R)/R\\
&=3\epsilon_r(R)\\
&\le0.01\\
\epsilon_r(R)&\le\frac{1}{3}\times10^{-2}
\end{align*}
$$

度量半径 $R$ 允许的相对误差是 $\frac{1}{3}\times10^{-2}$。

## Question 4

> 已知 $a=3.201,b=0.57$ 是经过四舍五入后得到的近似值，求 $a+b$ 的有效数字位数

$$
\frac{1}{2}\times10^{-3}+\frac{1}{2}\times10^{-2}=0.55\times10^{-2}<0.5\times10^{-1}
$$

因此是 2 位有效数字。

## Question 5

> 求方程 $x^2-56x+1=0$ 的两个根，使之至少具有四位有效数字。已知 $\sqrt{783}\approx27.982$

$$
x^2-56x+1=0\\
(x-28)^2=783\\
\begin{align*}
x_1=28+\sqrt{783}&=55.98\\
x_2=28-\sqrt{783}&=\frac{1}{28+\sqrt{783}}=0.01786\\
\end{align*}
$$

## Question 6

> 设 $Y_0=28$，按递推公式 $Y_n=Y_{n-1}-\frac{1}{100}\sqrt{783}$ 计算到 $Y_{100}$，若取 $\sqrt{783}\approx27.982$，试问计算 $Y_{100}$ 将有多大误差。

$$
x^*=\sqrt{783}/100\\
x=0.27982\\
|x^*-x|\le\frac{1}{2}\times10^{-5}\\
|Y^*_{100}-Y_{100}|=|(Y^*_{99}-Y_{99})-(x^*-x)|\le100|x^*-x|=\frac{1}{2}\times10^{-3}
$$

## Question 7

> 计算 $f=(\sqrt{2}-1)^6$，取 $\sqrt{2}=1.4$，利用下列公式计算，哪一个得到的结果最好？
> $$
> \frac{1}{(\sqrt{2}+1)^6}\quad(3-2\sqrt{2})^3\quad\frac{1}{(3+2\sqrt{2})^3}\quad99-70\sqrt{2}
> $$

利用偏导数：

$$
\begin{align*}
\epsilon(x)&=|x^*-x|\le\frac{1}{2}\times10^{-1}\\
\epsilon\left(\frac{1}{(x+1)^6}\right)&=\frac{6}{(x+1)^7}\epsilon(x)\Rightarrow\frac{6}{2.4^7}\epsilon(x)=0.013\epsilon(x)\\
\epsilon\left((3-2x)^3\right)&=6(3-2x)^2\epsilon(x)\Rightarrow6\times0.2^2\epsilon(x)=0.24\epsilon(x)\\
\epsilon\left(\frac{1}{(3+2x)^3}\right)&=\frac{6}{(3+2\sqrt{2})^4}\epsilon(x)\Rightarrow\frac{6}{5.8^4}\epsilon(x)=0.0053\epsilon(x)\\
\epsilon\left(99-70x\right)&=70\epsilon(x)\\
\end{align*}
$$

因此最好的是第 3 个。

## Question 8

> 已知描述某个实际问题的数学模型为：$u(x,y)=3xy^2+x^2/y$，其中 $x,y$ 由统计方法得到，分别为 $x=4,y=2$。统计方法的绝对误差限均为 $0.01$。求 $u$ 的绝对误差限和相对误差限。

$$
\begin{align*}
\epsilon(u)&=\left|\frac{\partial u}{\partial x}\right|\epsilon(x)+\left|\frac{\partial u}{\partial y}\right|\epsilon(y)\\
&=\left|3y^2+\frac{2x}{y}\right|\epsilon(x)+\left|6xy-\frac{x^2}{y^2}\right|\epsilon(y)\\
&=0.16+0.44\\
&=0.60\\
\epsilon_r(u)&=\frac{\epsilon(u)}{|u|}=\frac{0.60}{56}=0.0107
\end{align*}
$$

## Question 9

> 设
> $$
> A=\begin{bmatrix}
> 0.6&0.5\\
> 0.1&0.3
> \end{bmatrix}
> $$
> 计算 $A$ 的行范数，列范数，$2$-范数及 F-范数。

$$
\|A\|_{\infty}=\max_{1\le i\le n}\sum_{j=1}^{n}|a_{ij}|=1.1\\
\|A\|_{1}=\max_{1\le j\le n}\sum_{i=1}^{n}|a_{ij}|=0.8\\
$$

$$
A^TA=\begin{bmatrix}
0.37&0.33\\
0.33&0.34
\end{bmatrix}\\
\det(\lambda I-A^TA)=(\lambda-0.37)(\lambda-0.34)-0.1089=\lambda^2-0.71\lambda+0.0169\\
\lambda_{\max}=\frac{0.71+\sqrt{0.71^2-4*0.0169}}{2}= 0.6853\\
\|A\|_{2}=\sqrt{\lambda_{\max}}=0.8279
$$

$$
\|A\|_{F}=\left(\sum_{i=1}^{n}\sum_{j=1}^{n}|a_{ij}|^2\right)^{1/2}=\sqrt{0.71}=0.8426
$$

## Question 10

> 设 $A\in\mathbb R^{n\times n}$ 为对称正定矩阵，$\|\mathbf x\|_A=\sqrt{\mathbf x^TA\mathbf x}$，证明 $\|\mathbf x\|_A$ 为 $\mathbb R^n$ 上向量的一种范数。

非负性：
$$
\mathbf x^TA\mathbf x\ge0\Rightarrow\|\mathbf x\|_A=\sqrt{\mathbf x^TA\mathbf x}\ge0\\
\|\mathbf x\|_A=0\Leftrightarrow\mathbf x=\mathbf0
$$
齐次性：
$$
\|k\mathbf x\|_A=\sqrt{k\mathbf x^TAk\mathbf x}=|k|\sqrt{\mathbf x^TA\mathbf x}=|k|\|\mathbf x\|_A
$$
三角不等式：
$$
\begin{align*}
\|\mathbf x+\mathbf y\|_A&=\sqrt{(\mathbf x+\mathbf y)^TA(\mathbf x+\mathbf y)}\\
&=\sqrt{\sum_{i}\sum_{j}a_{ij}(x_j+y_j)(x_i+y_i)}\\
&=\sqrt{\sum_{i}\sum_{j}a_{ij}(x_ix_j+y_ix_j+x_iy_j+y_iy_j)}\\
&\le\sqrt{\sum_{i}\sum_{j}a_{ij}x_ix_j}+\sqrt{\sum_{i}\sum_{j}a_{ij}y_iy_j}\quad(\text{Cauchy inequality})\\
&=\|\mathbf x\|_A+\|\mathbf y\|_A
\end{align*}
$$

因此 $\|\mathbf x\|_A$ 为 $\mathbb R^n$ 上向量的一种范数。

## Question 11

> 设 $A$ 为非奇异矩阵，且 $\|A^{-1}\|\|\delta A\|<1$，求证 $(A+\delta A)^{-1}$ 存在且有估计：
> $$
> \frac{\|A^{-1}-(A+\delta A)^{-1}\|}{\|A^{-1}\|}\le\frac{\mathrm{cond}(A)\frac{\|\delta A\|}{\|A\|}}{1-\mathrm{cond}(A)\frac{\|\delta A\|}{\|A\|}}
> $$

先证明可逆：
$$
A+\delta A=A(I+A^{-1}\delta A)\\
\|A^{-1}\delta A\|\le\|A^{-1}\|\|\delta A\|<1
$$
因为 $\|A^{-1}\delta A\|<1$，所以 $I+A^{-1}\delta A$ 非奇异，所以 $(A+\delta A)^{-1}$ 存在。

设方程 $Ax=b$，则：
$$
\begin{align*}
(A+\delta A)(x+\delta x)&=b=Ax\\
\delta x&=(A+\delta A)^{-1}Ax-x\\
\delta x&=(A+\delta A)^{-1}Ax-A^{-1}Ax\\
\delta x&=((A+\delta A)^{-1}-A^{-1})Ax\\
\end{align*}
$$

$$
\begin{align*}
(A+\delta A)(x+\delta x)&=b=Ax\\
\delta Ax+A\delta x+\delta A\delta x&=0\\
\delta x&=-(A+\delta A)^{-1}\delta Ax\\
\delta x&=-(I+A^{-1}\delta A)^{-1}A^{-1}\delta Ax\\
\end{align*}
$$

两边取等：
$$
\begin{align*}
((A+\delta A)^{-1}-A^{-1})Ax&=-(I+A^{-1}\delta A)^{-1}A^{-1}\delta Ax\\
(A^{-1}-(A+\delta A)^{-1})Ax&=(I+A^{-1}\delta A)^{-1}A^{-1}\delta Ax\\
(A^{-1}-(A+\delta A)^{-1})&=(I+A^{-1}\delta A)^{-1}A^{-1}\delta AA^{-1}\\
\|A^{-1}-(A+\delta A)^{-1}\|&=\|(I+A^{-1}\delta A)^{-1}A^{-1}\delta AA^{-1}\|\\
&\le\|(I+A^{-1}\delta A)^{-1}\|\|A^{-1}\|\|\delta A\|\|A^{-1}\|\\
&\le\frac{\|\delta A\|\|A^{-1}\|}{1-\|A^{-1}\delta A\|}\|A^{-1}\|\\
&\le\frac{\|\delta A\|\|A^{-1}\|}{1-\|A^{-1}\|\|\delta A\|}\|A^{-1}\|\\
\frac{\|A^{-1}-(A+\delta A)^{-1}\|}{\|A^{-1}\|}&\le\frac{\|\delta A\|\|A^{-1}\|}{1-\|A^{-1}\|\|\delta A\|}\\
\frac{\|A^{-1}-(A+\delta A)^{-1}\|}{\|A^{-1}\|}&\le\frac{1}{1-\mathrm{cond}(A)\frac{\|\delta A\|}{\|A\|}}\mathrm{cond}(A)\frac{\|\delta A\|}{\|A\|}\\
\end{align*}
$$

## Question 12

> 设矩阵 
> $$
> A=\begin{bmatrix}
> 1&2&0\\
> 1&-1&-1\\
> 1&-1&1\\
> \end{bmatrix}
> $$
> 求 $\mathrm{cond}(A)_2$。

先求 $A$ 的范数
$$
A^TA=\begin{bmatrix}
3&0&0\\
0&6&0\\
0&0&2\\
\end{bmatrix}
\Longrightarrow
\|A\|_2=\sqrt{6}
$$
再求 $A^{-1}$ 的范数
$$
A^{-1}=\begin{bmatrix}
1/3&1/3&1/3\\
1/3&-1/6&-1/6\\
0&1/2&1/2\\
\end{bmatrix}\\
(A^{-1})^TA^{-1}=\begin{bmatrix}
2/9&1/18&1/18\\
1/18&7/18&-1/9\\
1/18&-1/9&7/18\\
\end{bmatrix}\Longrightarrow
\|A^{-1}\|_2=\sqrt{2}/2
$$
条件数：
$$
\mathrm{cond}(A)_2=\sqrt{3}
$$

