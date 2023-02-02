# Chapter 4 最小二乘法 Least Square Method

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 3.7 最小二乘法
</blockquote>

自变量和因变量的函数关系 $y=f(x)$，要求在给定点上的误差 $\delta_i=f(x_i)-y$ 按某种度量最小。

设 $f$ 为在 $m+1$ 个节点给定的离散函数，$\Phi=\mathrm{span}(\phi_0,\dots,\phi_n)$，最小二乘法就是在 $\Phi$ 中找一个 $\phi^*(x)$ ，使得：
$$
\sum_{i=0}^{m}w(x_i)[y_i-\phi^*(x_i)]^2=\min_{\phi\in\Phi}\sum_{i=0}^{m}w(x_i)[y_i-\phi(x_i)]^2
$$
由于 $\phi(x)$ 是基函数的线性组合，则：
$$
\min_{\phi\in\Phi}\sum_{i=0}^{m}w(x_i)[y_i-\phi(x_i)]^2=\min_{\phi\in\Phi}\sum_{i=0}^{m}w(x_i)\left[y_i-\sum_{j=0}^{n}a_j\phi_j(x_i)\right]^2
$$
相当于求多元函数
$$
I(a_0,\dots,a_n)=\sum_{i=0}^{m}w(x_i)\left[y_i-\sum_{j=0}^{n}a_j\phi_j(x_i)\right]^2
$$
的极小值。由极值的必要性，令
$$
\frac{\partial{I}}{\partial{a_k}}=2\phi_k(x_i)\sum_{i=0}^{m}w(x_i)\left[y_i-\sum_{j=0}^{n}a_j\phi_j(x_i)\right]=0
$$
引入离散内积：
$$
\sum_{j=0}^{n}\langle\phi_{k},\phi_{j}\rangle=\sum_{i=0}^{m}w(x_i)\phi_{k}(x_i)\phi_{j}(x_i)\\
\sum_{j=0}^{n}\langle f,\phi_{k}\rangle=\sum_{i=0}^{m}w(x_i)y_{i}\phi_{k}(x_i)\\
$$
与最佳平方逼近的形式一致。则可得到求解系数的法方程。法方程中系数矩阵可能为奇异阵。需增加 Haar 条件才能保证系数矩阵非奇异，从而法方程存在唯一解。

<blockquote style="border-left: 5px solid #bb4545; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(188, 70, 70, 0.1)">
    例题
</blockquote>

设有方程组 $Ax=b$，其中 $A\in\mathbb{R}^{m\times n},b\in\mathbb{R}^m,m>n$，求 $x^*\in\mathbb{R}^{n}$ 使得
$$
\|Ax^*-b\|_2=\min_{x\in\mathbb{R}^{n}}\|Ax-b\|
$$

> 可得到：
> $$
> A^TAx=A^Tb\\
> x=(A^TA)^{-1}A^Tb
> $$

