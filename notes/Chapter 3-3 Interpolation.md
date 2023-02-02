# Chapter 3 插值 Interpolation

## 3.5 贝塞尔曲线 Bezier Curve

### 3.5.1 Bernstein 多项式

$$
B_{i,n}(t)=\binom{n}{i}t^i(1-t)^{n-i}
$$

这里，$n$ 为多项式阶数，$i\in[0,n]$，$t\in[0,1]$。

Bernstein 多项式具有非负性，归一性，对称性。

#### 端点性质


$$
B_{i,n}(0)=\begin{cases}1&i=0\\0&i\ne0\end{cases}\\
B_{i,n}(1)=\begin{cases}1&i=n\\0&i\ne0\end{cases}\\
$$

#### 递归性


$$
B_{i,n}(t)=(1-t)B_{i,n-1}(t)+tB_{i-1,n-1}(t)
$$

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    证明
</blockquote>



$$
\begin{align*}
B_{i,n}(t)&=\binom{n}{i}t^i(1-t)^{n-i}\\
&=\left(\binom{n-1}{i}+\binom{n-1}{i-1}\right)t^i(1-t)^{n-i}\\
&=\binom{n-1}{i}t^i(1-t)^{n-i}+\binom{n-1}{i-1}t^i(1-t)^{n-i}\\
&=(1-t)\binom{n-1}{i}t^i(1-t)^{(n-1)-i}+t\binom{n-1}{i-1}t^{(i-1)}(1-t)^{n-i}\\
&=(1-t)B_{i,n-1}(t)+tB_{i-1,n-1}(t)
\end{align*}
$$

#### 一阶导数


$$
\frac{dB_{i,n}(t)}{dt}=n(B_{i-1,n-1}(t)-B_{i,n-1}(t))
$$

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    证明
</blockquote>



$$
\begin{align*}
\frac{dB_{i,n}(t)}{dt}&=\binom{n}{i}\frac{dt^i(1-t)^{n-i}}{dt}\\
&=i\binom{n}{i}t^{i-1}(1-t)^{n-i}-(n-i)\binom{n}{i}t^i(1-t)^{n-i-1}\\
&=n\binom{n-1}{i-1}t^{i-1}(1-t)^{n-i}-n\binom{n-1}{i}t^i(1-t)^{n-i-1}\\
&=n[B_{i-1,n-1}(t)-B_{i,n-1}(t)]
\end{align*}
$$

#### 二阶导数


$$
\frac{d^2B_{i}^{n}(t)}{dt^2}=n(n-1)[B_{i-2,n-2}(t)-2B_{i-1,n-2}(t)+B_{i,n-2}(t)]
$$

### 3.5.2 贝塞尔曲线

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义
</blockquote>



给定控制点 $P_0,P_1,\dots,P_n$，贝塞尔曲线上的任意一点 $P_{0,n}(t):[0,1]\to\mathbb R^3$ 定义为：
$$
P_{0,n}(t)=\sum_{i=0}^{n}P_iB_{i,n}(t)\quad t\in[0,1]
$$

Berstein 函数作为权重。控制点连起来形成折线，折线每一段叫做控制脚 leg。

Bezier 曲线阶数 $n$ = 节点数 - 1；曲线经过头尾两个节点；在空间中，任意平面穿过控制点折线的次数比穿过 Bezier 曲线的次数多；具有凸包性；仿射 Bezier 曲线等价于访射控制点然后再构建曲线。

### 3.5.3 De Casteljat 算法

用到了 Bernstein 多项式的递归性质，即
$$
\begin{align*}
P_{0,n}(t)&=\sum_{i=0}^{n}P_iB_{i,n}(t)\\
&=(1-t)\sum_{i=0}^{n-1}P_iB_{i,n-1}(t)+t\sum_{i=1}^{n}P_iB_{i-1,n-1}(t)\\
&=(1-t)\sum_{i=0}^{n-1}P_iB_{i,n-1}(t)+t\sum_{i=0}^{n-1}P_{i+1}B_{i,n-1}(t)\\
&={\color{red}(1-t)P_{0, n-1}(t)+tP_{1,n}(t)}
\end{align*}
$$

### 3.5.4 移动控制点

移动一个控制点，假设给 $P_k$ 一个位移 $\mathbf v$，由定义可得：
$$
\begin{align*}
Q_{0,n}(t)&=\sum_{i=0}^{n}P_iB_{i,n}(t)+\mathbf vB_{k,n}(t)\\
&=P_{0,n}(t)+\mathbf vB_{k,n}(t)
\end{align*}
$$
因此曲线会向 $\mathbf v$ 移动 $B_{k,n}(t)$ 长度，而且由于 $B_{k,n}(t)$ 是一个全域函数，因此作用也是全域的。

### 3.5.5 切向量和曲率

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    一阶导数
</blockquote>

对 $P_{0,n}(t)$ 求一阶导：
$$
\begin{align*}
\frac{dP_{0,n}(t)}{dt}&=\frac{d\sum_{i=0}^{n}P_iB_{i,n}(t)}{dt}\\
&=\sum_{i=0}^{n}P_i\frac{dB_{i,n}(t)}{dt}\\
&=\sum_{i=0}^{n}P_i\left[n(B_{i-1,n-1}(t)-B_{i,n-1}(t))\right]\\
&=n\left[\sum_{i=0}^{n}P_{i}B_{i-1,n-1}(t)-\sum_{i=0}^{n}P_iB_{i,n-1}(t)\right]\\
&=n\left[\sum_{i=0}^{n-1}(P_{i+1}-P_i)B_{i,n-1}(t)\right]\\
&={\color{red}n[P_{1,n}(t)-P_{0,n-1}(t)]}
\end{align*}
$$
特别的，端点切向量为：
$$
\dfrac{dP_{0,n}(t)}{dt}\mid_{t=0}=n(P_1-P_0)\\
\dfrac{dP_{0,n}(t)}{dt}\mid_{t=1}=n(P_n-P_{n-1})
$$
因此，Bezier 曲线和首尾两条折线段相切。

对于两条 Bezier 曲线，只要 $n_1(P_1-P_0)=n_2(P_{n_2}-P_{n_2-1})$，则两条曲线一阶连续。

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    二阶导数
</blockquote>

对 $P_{0,n}(t)$ 求二阶导：
$$
\begin{align*}
\frac{d^2P_{0,n}(t)}{dt^2}&=n\left[\sum_{i=0}^{n-1}(P_{i+1}-P_i)\frac{dB_{i}^{n-1}(t)}{dt}\right]\\
&=n\left[\sum_{i=0}^{n-1}(P_{i+1}-P_i)(n-1)(B_{i-1}^{n-2}(t)-B_{i}^{n-2}(t))\right]\\
&=n(n-1)\left[\sum_{i=0}^{n-1}P_{i+1}B_{i-1}^{n-2}(t)-P_iB_{i-1}^{n-2}(t)-P_{i+1}B_{i}^{n-2}(t)+P_iB_{i}^{n-2}(t)\right]\\
&=n(n-1)\left[\sum_{i=0}^{n-2}P_{i+2}B_{i}^{n-2}(t)-P_{i+1}B_{i}^{n-2}(t)-P_{i+1}B_{i}^{n-2}(t)+P_iB_{i}^{n-2}(t)\right]\\
&=n(n-1)\left[\sum_{i=0}^{n-2}(P_{i+2}-2P_{i+1}+P_i)B_{i}^{n-2}(t)\right]\\
\end{align*}
$$
特别的
$$
\frac{d^2P_{0,n}(t)}{dt^2}\mid_{t=0}=n(n-1)(P_2-2P_1+P_0)\\
\frac{d^2P_{0,n}(t)}{dt^2}\mid_{t=1}=n(n-1)(P_n-2P_{n-1}+P_{n-2})
$$
因此端点曲率为：
$$
\begin{align*}
\kappa(0)&=\frac{|n(P_1-P_0)\times n(n-1)[(P_2-P_1)-(P_1-P_0)]|}{n^3|P_1-P_0|^3}\\
&=\frac{n-1}{n}\frac{|(P_1-P_0)\times (P_2-P_1)|}{|P_1-P_0|^3}\\
\kappa(1)&=\frac{n-1}{n}\frac{|(P_{n-1}-P_{n-2})\times (P_n-P_{n-1})|}{|P_n-P_{n-1}|^3}\\
\end{align*}
$$

<blockquote style="border-left: 5px solid #b94263; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(185, 66, 110, 0.1)">
    注意
</blockquote>

接点处曲率相同，不一定二阶连续；二阶连续一定曲率相同。

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    高阶导数
</blockquote>

$$
\begin{align*}
\frac{d^kP_{0,n}(t)}{dt^k}&=n(n-1)\cdots(n-k+1)\left[\sum_{i=0}^{n-k}(\Delta P_i)^k B_{i}^{n-k}(t)\right]\\
\end{align*}
$$

其中 $(\Delta P_i)^k$ 为 $k$ 阶差分。

### 3.5.6 Bezier 曲线的拆分

给定控制点 $P_0,P_1,\dots,P_n$，将 bezier 曲线拆分为由控制点 $Q_0,Q_1,\dots,Q_n$ 和由控制点 $R_0,R_1,\dots,R_n$ 形成的 bezier 曲线。

这里，拆分点为 $P_{0,n}(t)$。因此，$Q_0,Q_1,\dots,Q_n$ 为 $P_{0,i}(t),i\in[0,n]$，而 $R_0,R_1,\dots,R_n$ 为 $P_{i,n}(t),i\in[0,n]$。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    证明
</blockquote>

证明前一段的拆分正确性即可。

$Q_0,Q_1,\dots,Q_n$ 形成的 bezier 曲线方程为：
$$
\begin{align*}
Q_{0,n}(t)&=\sum_{i=0}^{n}Q_iB_{i,n}(t)\\
&=\sum_{i=0}^{n}P_{0,i}(u)B_{i,n}(t)\\
&=\sum_{i=0}^{n}\left[\sum_{j=0}^{i}P_{j}B_{j,i}(u)\right]B_{i,n}(t)\\
&=\sum_{i=0}^{n}P_{i}\left[\sum_{j=i}^{n}B_{i,j}(u)B_{j,n}(t)\right]\\
&=\sum_{i=0}^{n}P_{i}\left[\sum_{j=i}^{n}\binom{j}{i}u^{i}(1-u)^{j-i}\binom{n}{j}t^{j}(1-t)^{n-j}\right]\\
&=\sum_{i=0}^{n}P_{i}(ut)^i\left[\sum_{j=i}^{n}\binom{j}{i}\binom{n}{j}[(1-u)t]^{j-i}(1-t)^{n-j}\right]\\
&=\sum_{i=0}^{n}P_{i}\binom{n}{i}(ut)^i\left[\sum_{j=i}^{n}\binom{n-i}{j-i}(t-ut)^{j-i}(1-t)^{n-j}\right]\\
&=\sum_{i=0}^{n}P_{i}\binom{n}{i}(ut)^i(1-ut)^{n-i}\\
&=P_{0,n}(ut)\\
\end{align*}
$$
因此拆分正确。

## 3.6 B-样条 B-Spline

### 3.6.1 基函数

用局域非零函数作为曲线的基。

假设 $0=t_0<t_1<\dots<t_m=1$ 是严格单调递增的数列，每个点都叫做节点 knot。

其中基函数定义为：
$$
\begin{align*}
B_{i,0}(t)&=
\begin{cases}
1&t_i\le t<t_{i+1}\\
0&\text{otherwise}
\end{cases}\\
B_{i,n}(t)&=\frac{t-t_i}{t_{i+n}-t_i}B_{i,n-1}(t)+\frac{t_{i+n+1}-t}{t_{i+n+1}-t_{i+1}}B_{i+1,n-1}(t)
\end{align*}
$$
其中，$n$ 是多项式次数。

$B_{i,0}(t)$ 表示了每一段的最低阶函数是闭区间上的方波函数。

$B_{i,n}(t)$ 表示了每一段的函数是 2 个 $n+1$ 个节点、$n$ 个区间上的加权和。

高阶基函数是低阶基函数的卷积。$B_{i,n}(t)$ 会在区间 $[t_i,t_{i+n+1}]$ 上非零。

#### 归一性

在区间 $[t_i,t_{i+1}]$ 上的所有 $n$ 阶基函数和为 $1$。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    证明
</blockquote>



$n=0$ 时：满足

$n=1$ 时：
$$
\begin{align*}
B_{i,1}(t)&=\frac{t-t_i}{t_{i+1}-t_i}B_{i,0}(t)+\frac{t_{i+2}-t}{t_{i+2}-t_{i+1}}B_{i+1,0}(t)\\
&=\frac{t-t_i}{t_{i+1}-t_i}\quad t\in[t_i,t_{i+1}]\\
B_{i-1,1}(t)&=\frac{t-t_{i-1}}{t_{i}-t_{i-1}}B_{{i-1},0}(t)+\frac{t_{i+1}-t}{t_{i+1}-t_{i}}B_{i,0}(t)\\
&=\frac{t_{i+1}-t}{t_{i+1}-t_{i}}\quad t\in[t_i,t_{i+1}]\\
B_{i-1,1}(t)+B_{i,1}(t)&=1
\end{align*}
$$
假设 $n=k$ 时满足：

易证
$$
B_{i-k,k+1}(t)+\dots+B_{i,k+1}(t)=1
$$

### 3.6.2 B-样条

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义
</blockquote>



给定控制顶点 $P_0,P_2,\dots,P_n$ 和节点 $0=t_0<t_1<\dots<t_m=1$，B-Spline 上的任意一点 $P_{0,n}(t):\mathbb R\to\mathbb R^3$ 定义为：
$$
P_{0,n}(t)=\sum_{i=0}^{n}P_iB_{i,p}(t)
$$
具有局部凸包性，但不具有全局凸包性。

假设有 $n+1$ 个控制点，$m+1$ 个节点，$p$ 是 B-Spline 的阶数：
$$
m=n+p+1
$$
可以通过设置节点的值形成：开放的，嵌制的，闭合的 B-Spline。

### 3.6.3 B-样条种类

如果第一个和最后一个节点不重复 $p+1$ 次，则不能形成像 Bezier 曲线一样的，曲线和首尾折线段相切的情况。

#### 开放的

开放的 B-Spline 值域为 $[t_p,t_{m-p}]$。

### 3.6.5 De Boor 算法

De Boor 算法和 De Casteljat 算法类似，都利用了递归的性质：
$$
\begin{align*}
P_{0,n}(t)&=\sum_{i=0}^{n}P_iB_{i,p}(t)\quad t\in[t_i,t_{i+p+1})\\
&=\sum_{i=0}^{n}P_i\left[\frac{t-t_i}{t_{i+p}-t_i}B_{i,p-1}(t)+\frac{t_{i+p+1}-t}{t_{i+n+1}-t_{i+1}}B_{i+1,p-1}(t)\right]\quad t\in[t_{i},t_{i+p+1})\\
&=\sum_{i=0}^{n}\frac{t-t_i}{t_{i+p}-t_i}P_iB_{i,p-1}(t)\quad t\in[t_{i},t_{i+p})\\
&+\sum_{i=0}^{n}\frac{t_{i+p+1}-t}{t_{i+p+1}-t_{i+1}}P_iB_{i+1,p-1}(t)\quad t\in[t_{i+1},t_{i+p+1})\\
\end{align*}
$$







