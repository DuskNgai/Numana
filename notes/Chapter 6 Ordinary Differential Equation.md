# Chapter 6 常微分方程 Ordinary Differential Equation

## 6.1 初值问题 Initial Value Problems

$$
\begin{cases}
y'=f(x,y)&(1)\\
y(x_0)=y_0&(2)
\end{cases}
$$

$f$ 如果在 $D=\{(x,y)|a\le x\le b,|y|<\infty\}$ 关于 $y$ 满足 Lipschitz 条件，则初值问题的解存在、唯一, 且连续依赖于初始条件。

### 6.1.1 数值方法的基本思想

设相邻两个节点的间距采用等距, 步长 **$h=x_{k+1}-x_{k},x_n=x_0+nh$**。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    化导数为差商
</blockquote>

$$
y'(x_i)\approx\frac{y_{i+1}-y_i}{h}\\
$$

得到：
$$
\begin{cases}
\dfrac{y_{i+1}-y_i}{h}=f(x_i,y_i)\\
y(x_0)=y_0
\end{cases}\Longrightarrow y_{n+1}=y_n+hf(x_n,y_n)
$$

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    泰勒展开法
</blockquote>

$$
y(x_{i+1})\approx y(x_i)+y'(x_i)h+\frac{y''(x_i)}{2!}h^2+\cdots+\frac{y^{(p)}}{p!}h^p
$$

各项导数带入 $(1)$。
$$
\begin{cases}
y'(x_i)&=&f(x_i,y(x_i))\\
y''(x_i)&=&f'_x(x_i,y(x_i))+f'_y(x_i,y(x_i))y'(x_i)\\
y'''(x_i)&=&f''_x(x_i,y(x_i))+2f''_{xy}(x_i,y(x_i))y'(x_i)+f''_{yy}(x_i,y(x_i))(y'(x_i))^2\\
&+&f'_y(x_i,y(x_i))y''(x_i)\\
&\vdots&
\end{cases}
$$

带入：
$$
y(x_{i+1})\approx y(x_i)+f(x_i,y(x_i))h+\frac{f'_x(x_i,y(x_i))+f'_y(x_i,y(x_i))y'(x_i)}{2!}h^2+\cdots
$$

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    数值积分法
</blockquote>

对方程在区间 $[x_n,x_{n+1}]$​ 上积分得到：
$$
y(x_{n+1})=y(x_n)+\int_{x_n}^{x_{n+1}}f(x,y(x))\mathrm{d}x
$$
并对右端积分项给出一种数值积分方法。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    微分中值定理
</blockquote>

$$
y(x_{n+1})-y(x_n)=hy'(\xi)=hf(\xi,y(\xi))
$$

且对 $f(\xi,y(\xi))$ 给出适当的数值方法。

### 6.1.2 数值方法的几个基本概念

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    单步多步
</blockquote>

单步：$y_{n+1}=y_n+hf(x_n,y(x_n))$

多步：$y_{n+1}=y_{n-1}+2hf(x_n,y(x_n))$

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    隐式显示
</blockquote>

隐式：$y_{n+1}=y_n+hf(x_{n+1},y(x_{n+1}))$

显示：$y_{n+1}=y_n+hf(x_n,y(x_n))$

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    截断误差和方法的阶数
</blockquote>

称误差 $T_{n+1}=y(x_{n+1})-y_{n+1}$ 为 $x_{n+1}$ 处的整体截断误差。(精确值 - 数值)

如果 $y(x_n)=y_n$ 即精确的，在此条件下，称 $T_{n+1}$ 为局部误差。如果局部截断误差阶为 $O(h^{p+1})$，则称该方法为 $p$ 阶方法。

## 6.2 欧拉方法

数值积分法：$y(x_{n+1})=y(x_n)+\int_{x_n}^{x_{n+1}}f(x,y(x))\mathrm{d}x$ 。

用左矩形公式：$y_{n+1}=y(x_n)+hf(x_n,y(x_n))=y(x_n)+hf(x_n,y_n)$，是显式欧拉法。

用右矩形公式：$y_{n+1}=y(x_n)+hf(x_{n+1},y(x_{n+1}))=y(x_n)+hf(x_{n+1},y_{n+1})$，是隐式欧拉法。

用中矩形公式：$y_{n+1}=y(x_n)+(h/2)[f(x_{n},y_{n})+f(x_{n+1},y_{n+1})]$。

在 $[x_{n-1},x_{n+1}]$ 上求积分，用中矩形公式：$y_{n+1}=y_{n-1}+2hf(x_n,y_n)$，是两步欧拉法。

### 6.2.2 局部截断误差

$$
T_{n+1}=y(x_{n+1})-[y(x_n)+h\phi(x_n,y_n,h)]
$$

#### 显式欧拉法

$$
T_{n+1}=y(x_{n+1})-[y(x_n)+hf(x_n,y_n)]
$$

由 $y(x_n)=y_n$ 可得 $f(x_n,y_n)=f(x_n,y(x_n))=y'(x_n)$，因此 Taylor 展开 $y(x_{n+1})$：
$$
T_{n+1}=y(x_{n+1})-[y(x_n)+hy'(x_n)]=\frac{y''(x_n)}{2}h^2+O(h^3)=O(h^2)
$$
因此显式欧拉方法为一阶方法。

#### 隐式欧拉法

$$
T_{n+1}=y(x_{n+1})-[y(x_n)+hf(x_{n+1},y_{n+1})]
$$

记：$\tilde{y}_{n+1}=y_n+hf(x_{n+1},y(x_{n+1}))$，则：
$$
T_{n+1}=y(x_{n+1})-\tilde{y}_{n+1}+\tilde{y}_{n+1}-y_{n+1}
$$
Taylor 展开：
$$
\tilde{y}_{n+1}=y(x_n)+hy'(x_{n+1})=y(x_n)+hy'(x_n)+h^2y''(x_n)+O(h^3)
$$
前一项：
$$
y(x_{n+1})-\tilde{y}_{n+1}=-h^2\frac{y''(x_n)}{2}+O(h^3)
$$
后一项：
$$
\begin{align*}
\tilde{y}_{n+1}-y_{n+1}&=y_n+hf(x_{n+1},y(x_{n+1}))-y(x_n)-hf(x_{n+1},y_{n+1})\\
&=h[f(x_{n+1},y(x_{n+1}))-f(x_{n+1},y_{n+1})]\\
&=hf'_y(x_{n+1},\xi)(y(x_{n+1})-y_{n+1})\\
&=hf'_y(x_{n+1},\xi)T_{n+1}\\
\end{align*}
$$
因此：
$$
T_{n+1}=\frac{1}{1-hf'_y(x_{n+1},\xi)}\left[-h^2\frac{y''(x_n)}{2}+O(h^3)\right]=O(h^2)
$$
隐式欧拉方法为一阶方法。

#### 梯形方法

$$
T_{n+1}=y(x_{n+1})-[y(x_n)+(h/2)[f(x_{n},y_{n})+f(x_{n+1},y_{n+1})]]
$$

Taylor 展开 $y(x_{n+1})$，得到：
$$
\begin{align*}
y(x_{n+1})&=y(x_n)+f(x_n,y(x_n))h+\frac{h^2}{2!}[f'_x(x_n,y(x_n))+f'_y(x_n,y(x_n))y'(x_n)]+O(h^3)\\
&=y(x_n)+f(x_n,y_n)h+\frac{h^2}{2!}[f'_x(x_n,y_n)+f'_y(x_n,y_n)y'(x_n)]+O(h^3)\\
\end{align*}
$$
二维 Taylor 展开 $f(x_{n+1},y_{n+1})$：
$$
\begin{align*}
f(x_{n+1},y_{n+1})&=f(x_{n}+h,y_{n}+hf(x_{n},y_{n}))\\
&=f(x_{n},y_{n})+hf'_x(x_{n},y_{n})+hf(x_{n},y_{n})f'_y(x_{n},y_{n})+O(h^2)\\
\end{align*}
$$
合并得到：
$$
\begin{align*}
&\quad\ \ y(x_n)+(h/2)[f(x_{n},y_{n})+f(x_{n+1},y_{n+1})]\\
&=y(x_n)+\frac{h}{2}[2f(x_{n},y_{n})+hf'_x(x_{n},y_{n})+hf(x_{n},y_{n})f'_y(x_{n},y_{n})+O(h^2)]\\
&=y(x_n)+hf(x_{n},y_{n})+\frac{h^2}{2}[f'_x(x_{n},y_{n})+f(x_{n},y_{n})f'_y(x_{n},y_{n})]+O(h^3)\\
\end{align*}
$$
因此 
$$
T_{n+1}=O(h^3)
$$

#### 两步欧拉方法

$$
T_{n+1}=y(x_{n+1})-[y_{n-1}+2hf(x_n,y_n)]
$$

这里需要假设 $y_{n-1}=y(x_{n-1}),y_{n}=y(x_n)$。

Taylor 展开 $y(x_{n-1})$
$$
\begin{align*}
y(x_{n-1})&=y(x_{n}-h)\\
&=y(x_n)-hf'(x_n)+\frac{h^2}{2}f''(x_n)-\frac{h^3}{6}f'''(x_n)+O(h^4)
\end{align*}
$$
带入：
$$
\begin{align*}
T_{n+1}&=y(x_{n+1})-[y_{n-1}+2hf(x_n,y_n)]\\
&=y(x_{n+1})-[y(x_n)+hf'(x_n)+\frac{h^2}{2}f''(x_n)-\frac{h^3}{6}f'''(x_n)+O(h^4)]\\
&=\frac{h^3}{3}f'''(x_n)+O(h^4)
\end{align*}
$$
因此 $T_{n+1}=O(h^3)$。

### 6.2.3 隐式方法的迭代计算

隐式欧拉方法和梯形方法均为一种隐式的单步方法。

梯形方法：
$$
y_{n+1}=\frac{h}{2}f(x_{n+1},y_{n+1})+y_n+\frac{h}{2}f(x_n,y_n)
$$
构造迭代格式：
$$
y^{[k+1]}_{n+1}=\frac{h}{2}f(x_{n+1},y^{[k]}_{n+1})+[y_n+\frac{h}{2}f(x_n,y_n)]
$$
当 $f$ 关于变量 $y$ 满足 Lipschitz 条件时候：
$$
\begin{align*}
|y^{[k+1]}_{n+1}-y^{[k]}_{n+1}|&=\frac{h}{2}|f(x_{n+1},y^{[k]}_{n+1})-f(x_{n+1},y^{[k-1]}_{n+1})|\\
&\le\frac{hL}{2}|y^{[k]}_{n+1}-y^{[k-1]}_{n+1}|\\
&\le(\frac{hL}{2})^{k}|y^{[1]}_{n+1}-y^{[0]}_{n+1}|\\
\end{align*}
$$
因此只需要 $h<2/L$ 就可以收敛。

### 6.2.4 改进欧拉方法

对梯形方法作如下修正

先用显式欧拉方法：
$$
\bar{y}_{n+1}=y_n+hf(x_n,y_n)
$$
再带入梯形法：
$$
y_{n+1}=y(x_n)+\frac{h}{2}[f(x_{n},y_{n})+f(x_{n+1},\bar{y}_{n+1})]
$$

## 6.3 Runge-Kutta 方法

微分中值定理的方法，对于 $K^*=y'(\xi)$ 的估计。显式欧拉，隐式欧拉是一阶方法。梯形法是二阶方法。

### 6.3.1 二阶 Runge-Kutta 方法

在 $[x_n,x_{n+1}]$ 内找一点 $x_{n+p}=x_n+ph,p\in(0,1]$。仿照改进欧拉方法，两点斜率的加权平均作为平均斜率, 即
$$
K^*=\lambda_1f(x_n,y_n)+\lambda_2f(x_{n+p},y_n+phf(x_n,y_n))
$$
或：
$$
\begin{cases}
y_{n+1}=y_n+hK^*\\
K^*=\lambda_1K_1+\lambda_2K_2\\
K_1=f(x_n,y_n)\\
K_2=f(x_{n}+ph,y_n+phf(x_n,y_n))
\end{cases}
$$
类似梯形方法可以得到 $T_{n+1}=O(h^3)$。

---

二维 Taylor 展开 $K_2$：
$$
\begin{align*}
K_2&=f(x_{n}+ph,y_n+phf(x_n,y_n))\\
&=f(x_{n},y_{n})+phf'_x(x_{n},y_{n})+phf(x_{n},y_{n})f'_y(x_{n},y_{n})+O(h^2)\\
&=K_1+phy''(x_n)+O(h^2)
\end{align*}
$$
合并得到：
$$
\begin{align*}
K^*&=(\lambda_1+\lambda_2)K_1+\lambda_2phy''(x_n)+O(h^2)\\
\end{align*}
$$
因此 
$$
\begin{align*}
T_{n+1}&=(1-\lambda_1-\lambda_2)y'(x_n)+(\frac{1}{2}-\lambda_2p)h^2y''(x_n)+O(h^3)\\
\end{align*}
$$
因此 $\lambda_1+\lambda_2=1,2\lambda_2 p=1$ 时候，误差为 2 阶。

### 6.3.2 $N$ 阶 Runge-Kutta 方法

把 $[x_n,x_{n+1}]$ 均匀地分成 $n$ 分。$N$ 阶 Runge-Kutta 方法的一般形式为：
$$
y_{n+1}=y_n+h\sum_{i=1}^{N}c_iK_i\\
\begin{cases}
{\displaystyle K_{1}=f(x_n,y_n)}\\
{\displaystyle K_{2}=f(x_n+a_2h,y_n+h(b_{21}K_1))}\\
{\displaystyle K_{3}=f(x_n+a_3h,y_n+h(b_{31}K_1+b_{32}K_2))}\\
\qquad\vdots\\
{\displaystyle K_{N}=f(x_n+a_Nh,y_n+h\sum_{j=1}^{N-1}b_{ij}K_j))}\\
{\displaystyle a_i=\sum_{j=1}^{i-1}b_{ij}, \sum_{i=1}^{N}c_i=1}
\end{cases}
$$
经典 4 阶 Runge-Kutta 格式：
$$
y_{n+1}=y_n+h\left(\frac{1}{6}K_1+\frac{2}{6}K_2+\frac{2}{6}K_3+\frac{1}{6}K_4\right)\\
\begin{cases}
{\displaystyle K_{1}=f\left(x_n,y_n\right)}\\
{\displaystyle K_{2}=f\left(x_n+\frac{1}{2}h,y_n+h\frac{1}{2}K_1\right)}\\
{\displaystyle K_{3}=f\left(x_n+\frac{1}{2}h,y_n+h\frac{1}{2}K_2\right)}\\
{\displaystyle K_{4}=f\left(x_n+h,y_n+hK_3\right)}
\end{cases}
$$

## 6.4 单步方法的收敛性和稳定性

### 6.4.1 收敛性和相容性

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    收敛性
</blockquote>

对于所有初值问题的一种单步法
$$
y_{n+1}=y(x_n)+h\phi(x_n,y_n,h)
$$
所产生的近似解，如果对任一固定的点 $x_n$，均有 $\lim_{h\to\infty}y_n=y(x_n)$，则称该方法为收敛的。

<blockquote style="border-left: 5px solid #b94263; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(185, 66, 110, 0.1)">
    例题
</blockquote>

用显式欧拉方法求解下列问题
$$
y'=ax+b,y(0)=0
$$
证明其收敛性。

> $$
> y(x)=\frac{a}{2}x^2+bx\Rightarrow y(x_n)=\frac{a}{2}x_n^2+bx_n=\frac{an^2h^2}{2}+bnh
> $$
>
> 显示欧拉方法：
> $$
> y_{n+1}=y_n+ahx_n+bh=y_n+anh^2+bh
> $$
>
> $$
> \begin{align*}
> e_{i+1}&=y(x_{n+1})-y_{n+1}\\
> &=\frac{a(n+1)^2h^2}{2}+b(n+1)h-y_n-anh^2-bh\\
> &=\frac{a(n^2+1)h^2}{2}+bnh-y_n\\
> &=y(x_n)-y_n+\frac{a}{2}h^2=e_n+\frac{a}{2}h^2\\
> &=e_0+\frac{a}{2}nh^2=\frac{a}{2}nh^2
> \end{align*}
> $$
>
> 因此 $h\to0$ 时误差为零。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    单步方法收敛
</blockquote>

如果单步方法 $y_{n+1}=y_n+h\phi(x_n,y_n,h)$ 是 $p$​ 阶的，且
$$
|\phi(x,y_1,h)-\phi(x,y_2,h)|\le L|y_1-y_2|
$$
则该方法收敛，且 $y(x_n)-y_n=O(h^p)$。

> 设 $e_n=y(x_n)-y_n$，由局部截断误差定义得
> $$
> \begin{align*}
> y(x_{n+1})&=y(x_{n})+h\phi(x_n,y(x_{n}),h)+T_{n+1}\\
> y(x_{n+1})-y_{n+1}&=y(x_{n})+h\phi(x_n,y(x_{n}),h)-y_n-h\phi(x_n,y_n,h)+T_{n+1}\\
> e_{n+1}&=e_n+h[\phi(x_n,y(x_{n}),h)-\phi(x_n,y_n,h)]+T_{n+1}\\
> |e_{n+1}|&=|e_n|+h|\phi(x_n,y(x_{n}),h)-\phi(x_n,y_n,h)|+Ch^{p+1}\\
> |e_{n+1}|&\le|e_n|+hL|y(x_{n})-y_n|+Ch^{p+1}\\
> &\le(1+hL)|e_n|+Ch^{p+1}\\
> &\le\alpha|e_n|+\beta\\
> &\le\alpha(\alpha|e_{n-1}|+\beta)+\beta=\alpha^2|e_{n-1}|+\alpha\beta+\beta\\
> &\le\alpha^{n+1}|e_0|+\beta\sum_{i=0}^{n}\alpha^{i}
> \end{align*}
> $$
> $\alpha=(1+hL),\beta=Ch^{p+1}$。因此：
> $$
> \begin{align*}
> |e_n|&\le\alpha^{n}|e_0|+\beta\sum_{i=0}^{n-1}\alpha^{i}\\
> &=(1+hL)^{\frac{1}{hL}nhL}|e_0|+\beta\frac{(1+hL)^n-1}{(1+hL)-1}\\
> &\le e^{nhL}|e_0|+\beta\frac{e^{nhL}-1}{hL}\\
> &=e^{nhL}|e_0|+\frac{e^{nhL}-1}{L}Ch^p\\
> \end{align*}
> $$

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    相容性
</blockquote>

方法 $y_{n+1}=y_n+h\phi(x_n,y_n,h)$ 满足 $\phi(x,y,0)=f(x,y)$，则称上述方法与初值问题相容。

### 6.4.2 稳定性

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    稳定性
</blockquote>

如果方法在节点值 $y_n$ 上产生大小为 $\delta$ 的误差，而以后各节点值 $y_{n+k}(k\ge1)$ 上产生的误差均不超过 $\delta$，则称该方法为稳定的。

<blockquote style="border-left: 5px solid #b94263; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(185, 66, 110, 0.1)">
    例题
</blockquote>

讨论显式欧拉方法的稳定性。

> 设计算 $y_n$ 产生的误差为 $\epsilon_n$，则实际得到为 $\bar{y}_n=y_n+\epsilon_n$​，则再计算一步得到的是
> $$
> \bar{y}_{n+1}=\bar{y}_n+hf(x_n,\bar{y}_n)\\
> \begin{align*}
> \epsilon_{n+1}&=\bar{y}_{n+1}-y_{n+1}\\
> &=\bar{y}_n+hf(x_n,\bar{y}_n)-y_n+hf(x_n,y_n)\\
> &=\epsilon_n+h[f(x_n,\bar{y}_n)-f(x_n,y_n)]\\
> &=\epsilon_n+hf'_y(x_n,\eta)(\bar{y}_n-y_n)\\
> &=\epsilon_n[1+hf'_y(x_n,\eta)]\\
> \end{align*}
> $$
> 当 $|1+hf'_y(x_n,\eta)|<1$ 时候，误差不增长, 从而是稳定的。

<blockquote style="border-left: 5px solid #b94263; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(185, 66, 110, 0.1)">
    例题
</blockquote>

用 $y'=\lambda y,\lambda<0$ 讨论欧拉方法、梯形方法的稳定性。

> 显式欧拉公式：$\epsilon_{n+1}=(1+\lambda h)\epsilon_n$ => $0<h<-2/\lambda$。有条件稳定。
>
> 隐式欧拉公式：$\epsilon_{n+1}=\epsilon_n+\lambda h\epsilon_{n+1}$ => $|1/(1-\lambda h)|<1$。无条件稳定。
>
> 梯形公式：$\epsilon_{n+1}=\epsilon_n+(\lambda h)(\epsilon_n+\epsilon_{n+1})/2$ => $|(2+\lambda h)/(2-\lambda h)|<1$。无条件稳定。

## 6.5 线性多步方法

### 6.5.1 一般形式的线性多步法

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    线性多步法
</blockquote>

称求解常微分方程初值问题的如下格式为线性多步方法
$$
y_{n+1}=\sum_{k=0}^ra_ky_{n-k}+h\sum_{k=-1}^rb_ky'_{n-k}
$$
其中 $y'_{n-k}=f(x_{n-k},y_{n-k})$。$b_{-1}=0$ 为显式，否则为隐式。

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    线性多步法的局部截断误差
</blockquote>

$$
T_{n+1}=y(x_{n+1})-\left(\sum_{k=0}^ra_ky_{n-k}+h\sum_{k=-1}^rb_ky'_{n-k}\right)
$$

这里需要规定 $y_{n-k}=y(x_{n-k})$。直接 Taylor 展开就行。

### 6.5.2 基于数值积分的线性多步法

将微分方程在区间 $[x_{n+1-l},x_{n+1}]$ 积分。

### 6.5.3 待定系数法

<blockquote style="border-left: 5px solid #b94263; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(185, 66, 110, 0.1)">
    例题
</blockquote>

确定线性多步方法
$$
y_{n+1}=a(y_n+y_{n-1})+h(b_0f_n+b_1f_{n-1})
$$

> $$
> a=1/2,b_0=7/4,b_1=-1/4
> $$

## 6.6 方程组与高阶方程

### 6.6.1 一阶方程组

$$
\begin{cases}
\mathbf{y}=f(x,\mathbf{y})\\
\mathbf{y}(x_0)=\mathbf{y}_0
\end{cases}
$$

分开构造迭代格式即可。

<blockquote style="border-left: 5px solid #b94263; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(185, 66, 110, 0.1)">
    例题
</blockquote>

$$
\begin{cases}
y'=f(x,y,z)\\
z'=g(x,y,z)\\
y(x_0)=y_0,z(x_0)=z_0
\end{cases}
$$

构造求解上述方程组的欧拉方法、梯形方法、改进欧拉方法。

> 显式欧拉方法
> $$
> \begin{cases}
> y_{n+1}=y_n+hf(x_n,y_n,z_n)\\
> z_{n+1}=z_n+hf(x_n,y_n,z_n)
> \end{cases}
> $$
> 隐式欧拉方法
> $$
> \begin{cases}
> y_{n+1}=y_n+hf(x_n,y_{n+1},z_{n+1})\\
> z_{n+1}=z_n+hf(x_n,y_{n+1},z_{n+1})
> \end{cases}
> $$
> Jacobi型：
> $$
> \begin{cases}
> y_{n+1}^{[k+1]}=y_n+hf(x_n,y_{n+1}^{[k]},z_{n+1}^{[k]})\\
> z_{n+1}^{[k+1]}=z_n+hf(x_n,y_{n+1}^{[k]},z_{n+1}^{[k]})
> \end{cases}
> $$
> Seidel型：
> $$
> \begin{cases}
> y_{n+1}^{[k+1]}=y_n+hf(x_n,y_{n+1}^{[k]},z_{n+1}^{[k]})\\
> z_{n+1}^{[k+1]}=z_n+hf(x_n,y_{n+1}^{[k+1]},z_{n+1}^{[k]})
> \end{cases}
> $$
> 梯形方法
> $$
> \begin{cases}
> y_{n+1}=y_n+\dfrac{h}{2}[f(x_n,y_n,z_n)+f(x_n,y_{n+1},z_{n+1})]\\
> z_{n+1}=z_n+\dfrac{h}{2}[f(x_n,y_n,z_n)+f(x_n,y_{n+1},z_{n+1})]
> \end{cases}
> $$
> 改进欧拉方法
> $$
> \begin{cases}
> \bar{y}_{n+1}=y_n+hf(x_n,y_n,z_n)\\
> \bar{z}_{n+1}=y_n+hf(x_n,y_n,z_n)\\
> y_{n+1}=y_n+\dfrac{h}{2}[f(x_n,y_n,z_n)+f(x_n,\bar{y}_{n+1},\bar{z}_{n+1})]\\
> z_{n+1}=z_n+\dfrac{h}{2}[f(x_n,y_n,z_n)+f(x_n,\bar{y}_{n+1},\bar{z}_{n+1})]
> \end{cases}
> $$

<blockquote style="border-left: 5px solid #b94263; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(185, 66, 110, 0.1)">
    例题
</blockquote>

用梯形公式解
$$
\begin{cases}
y'=-x\\
x'=y\\
x(0)=x_0,y(0)=y_0
\end{cases}
$$

> $$
> yy'=-xy,xx'=xy\\
> xx'+yy'=0\\
> x^2+y^2=x_0^2+y_0^2=C
> $$
>
> 梯形公式可以保证解仍然在相平面上。





