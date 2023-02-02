# Chapter 1 求解方程 Solve Equation

## 1.1 区间搜索法

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    零点定理
</blockquote>

$$
\forall f\in C[a,b],f(a)f(b)<0,\exists\xi\in(a,b),f(\xi)=0
$$

如果 $f(x)=(x-x^*)^mg(x)$，其中 $m$ 为正整数，$g(x^*)\ne0$，称 $x^*$ 为方程 $f(x)=0$ 的 $m$ 重根。

### 1.1.1 二分法

假设方程在 $(a, b)$ 内仅有唯一实根。

```python
fa, fb = f(a), f(b)
tolerance = 2 * np.finfo(float).eps * max(1, abs(a))
while abs(b - a) > tolerance:
    c = (a + b) / 2
    fc = f(c)
    if fc == 0.0:
        return c
    if fa * fc < 0:
        b, fb = c, fc
    else:
        a, fa = c, fc
return (a + b) / 2
```

其中每次迭代的 $x_n$ 满足：
$$
|x_n-x^*|<\frac{b-a}{2^n}
$$
可得二分所需的次数：
$$
n>\frac{\ln(b-a)-\ln\epsilon}{\ln2}
$$
只要求函数连续。不能求偶数重根, 也不能求复根。

弦截法就是用连接区间端点的直线与 $x$ 轴的交点代替二分法中区间的中点。

<blockquote style="border-left: 5px solid #bb4545; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(188, 70, 70, 0.1)">
    例题
</blockquote>

用二分法求方程 $x^3+4x^2-10=0$ 在区间 $[1,2]$ 内的根。求两步等分后的有根区间，误差不超过 $0.5\times10^{-3}$ 的等分次数。

> $f(1)=-5,f(2)=14,f'(x)=3x^2+8x>0$，因此唯一实根。
>
> (1) $[1.25,1.5]$.
>
> (2) 
> $$
> n>\frac{\ln1-\ln0.5\times10^{-3}}{\ln2}=\frac{\ln2000}{\ln2}=\log_{2}2000
> $$
> $n=10$

## 1.2 不动点迭代 Fixed Point Iteration

闭区间上连续可导。

```python
x[i + 1] = g(x[i])
```

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 误差估计
</blockquote>

设迭代函数 $\phi(x)$ 满足唯一不动点条件，则迭代格式 $x_{k+1}=\phi(x_k)$ 对于 $\forall x_0\in[a,b]$ 收敛，且有误差估计：
$$
\begin{gather*}
|x_k-x^*|\le\frac{L}{1-L}|x_k-x_{k-1}|\\
|x_k-x^*|\le\frac{L^k}{1-L}|x_1-x_0|
\end{gather*}
$$
上式为事后估计，下式为先验估计（证明同线性方程组迭代解）。

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义  p 阶收敛
</blockquote>


令 $e_i = |x_i - r|$ 表示迭代过程中第 $i$ 步时的误差，如果：
$$
\lim_{i\to\infty}\frac{e_{i+1}}{e_i^p}=S<1
$$
则称该迭代公式具有 $p$ 阶收敛。$p=1$ 称该方法满足**线性收敛**，收敛速度为 $S$。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 收敛到不动点
</blockquote>

设函数 $g$ 是连续可微函数，$g(r)=r$，$S=|g'(r)|<1$，则不动点迭代对于一个足够接近 $r$ 的初始估计，以速度 $S$ 线性收敛到不动点 $r$。

> 令 $x_i$ 表示第 $i$ 步迭代，根据 Lagrange 中值定理，在 $x_i$ 和 $r$ 之间存在 $c_i$​ 满足：
> $$
> \begin{gather*}
> \frac{g(x_i)-g(r)}{x_i-r}=g'(c_i)\\
> x_{i+1}-r=g'(c_i)(x_i-r)
> \end{gather*}
> $$
> 带入 $e_i = |x_i - r|$，则：
> $$
> e_{i+1}=|g'(c_i)|e_i
> $$
> 如果 $S=|g'(r)|<1$，则在 $r$ 的邻域上满足 $S<|g'(x)|<(S+1)/2<1$。选取恰当的 $x_i$ 使得其在该邻域内，则
> $$
> e_{i+1}\le\frac{S+1}{2}e_i
> $$
> 因此误差逐步下降，迭代收敛。

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定义 局部收敛
</blockquote>

如果迭代方法对于一个足够接近 $r$ 的初始值能够收敛到 $r$，则该迭代方法被称为**局部收敛**到 $r$。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 不动点收敛加速
</blockquote>

设对于二阶连续可微函数 $g$ 使用不动点迭代，对于一个不动点 $r$，有 $g'(r)=0,\dots,g^{(p-1)}(r)=0$，$g^{(p)}(r)\ne0$，若收敛到 $r$，则误差遵从 $\lim_{i\to\infty}(e_{i+1}/e^p_i)=|g''(r)|/2$。

> 对 $g(x)$ 在 $r$ 处进行 Taylor 展开
> $$
> \begin{align*}
> g(x_i)&=r+g'(r)(x_i-r)+\dots+\frac{g^{p}(c_i)}{p!}(x_i-r)^{p}\\
> g(x_i)-r&=\frac{g^{p}(c_i)}{p!}(x_i-r)^{p}\\
> x_{i+1}-r&=\frac{g^{p}(c_i)}{p!}(x_i-r)^{p}\\
> e_{i+1}&=\left|\frac{g^{(r)}(c_i)}{p!}\right|e_i^p
> \end{align*}
> $$
> $c_i$ 为在 $x_i$ 和 $r$ 之间的某个数。
>
> 由于 $g'(r)=0,\dots,g^{(p-1)}(r)=00$，则 $e_{i}$ 收敛到 $0$，因此
> $$
> \lim_{i\to\infty}\frac{e_{i+1}}{e_i^p}=\lim_{i\to\infty}\left|\frac{g^{(p)}(c_i)}{p!}\right|=\left|\frac{g^{(p)}(r)}{p!}\right|
> $$
> 

<blockquote style="border-left: 5px solid #bb4545; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(188, 70, 70, 0.1)">
    例题
</blockquote>

$x^3-x^2-x-1=0$ 在 $[1,2]$ 内的根。

> $\phi_1(x)=x^3-x^2-1$:
> $$
> |\phi_1'(x)|=|3x^2-2x|>1
> $$
> $\phi_2(x)=\sqrt[3]{x^2+x+1}$:
> $$
> |\phi_2'(x)|=\left|\frac{1}{3}(x^2+x+1)^{-\frac{2}{3}}(2x+1)\right|<1
> $$
> $\phi_3(x)=1+1/x+1/x^2$:
> $$
> |\phi_3'(x)|=|-1/x^2-2/x^3|>1
> $$

<blockquote style="border-left: 5px solid #bb4545; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(188, 70, 70, 0.1)">
    例题
</blockquote>

$3-3x-2\sin x=0$ 在 $[0,1]$ 内。证明对于任意初始值，$x_{k+1}=1-(2\sin x_k)/3$ 均收敛。$x_0=0$ 时候经过多少次迭代才能保证误差小于 $10^{-3}$。

> $$
> \begin{gather*}
> \frac{1}{3}<1-\frac{2}{3}\sin x_0<\frac{5}{3}\\
> 0<1-\frac{2}{3}\sin x_1<1\\
> \end{gather*}
> $$
>
> $$
> |\phi'(x)|=\left|-\frac{2}{3}\cos x\right|=\frac{2}{3}=L<1
> $$
>
> $$
> |x_n-x^*|\le\frac{L^n}{1-L}|x_1-x_0|
> $$
>
> 可得。

### 1.2.3 Aitken 加速技术

设迭代格式 $x_{k+1}=\phi(x_k)$ 线性收敛，则：
$$
x_{k+1}-x^*\approx c(x_k-x^*)\\
x_{k+2}-x^*\approx c(x_{k+1}-x^*)\\
$$
消去 $c$：
$$
x_{k+1}-x^*=\frac{x_{k+2}-x^*}{x_{k+1}-x^*}(x_k-x^*)\\
x^*=\frac{x_{k+2}x_k-x_{k+1}^2}{x_{k+2}-2x_{k+1}+x_k}
$$

可以使得原来发散的迭代仍然收敛。

## 1.3 精度的极限 limit of Precision

### 敏感性问题

假设找到了 $f(x)=0$ 的根 $r$，但是对输入 $f$ 做了一个很小的变化 $\epsilon g(x)$。令 $\Delta r$ 是对应根的变化，即：
$$
f(r+\Delta r)+\epsilon g(r+\Delta r)=0
$$
一阶 Taylor 展开：
$$
f(r)+(\Delta r)f'(r)+\epsilon g(r)+\epsilon(\Delta r)g'(r)+O((\Delta r)^2)=0
$$
因此：
$$
\begin{align*}
(\Delta r)(f'(r)+\epsilon g'(r))&\approx-(f(r)+\epsilon g(r))=-\epsilon g(r)\\
\Delta r&\approx\frac{-\epsilon g(r)}{f'(r)+\epsilon g'(r)}
\end{align*}
$$
当 $\epsilon<<f'(r)$ 时，得到根的敏感公式
$$
\Delta r\approx-\epsilon \frac{g(r)}{f'(r)}
$$

### 误差放大因子

设找到了 $f(x)=0$ 的根 $r$ 的一个近似值 $r'$，对应 $f(r')+\epsilon g(r')=0$，则定义**相对后向误差**为 $|f(r)/g(r)|=|\epsilon g(r)/g(r)|$，**相对前向误差**为 $|(r'-r)/r|$。

定义误差放大因子为
$$
|\frac{(r'-r)/r}{\epsilon g(r)/g(r)}|=|\frac{g(r)}{rf'(r)}|
$$

## 1.4 牛顿法 Newton's Method

闭区间上连续二阶可导。

```python
x[i + 1] = x[i] - (f(x[i]) / df(x[i]))
```

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 牛顿法
</blockquote>

令 $f$ 为二阶连续可微函数，$f(r)=0$，如果 $f'(r)\ne0$，则牛顿法局部二次收敛到 $r$。 误差 $e_i$ 遵从
$$
\lim_{i\to\infty}\frac{e_{i+1}}{e^2_i}=M=\frac{f''(r)}{2f'(r)}
$$

> 先证明<u>局部收敛</u>：
>
> 注意到牛顿法是不动点迭代法的一个特殊形式：
> $$
> g(x)=x-\frac{f(x)}{f'(x)}
> $$
> 因此：
> $$
> g'(x)=1-\frac{(f'(x))^2-f(x)f''(x)}{(f'(x))^2}=\frac{f(x)f''(x)}{(f'(x))^2}
> $$
> 由于 $f(r)=0$，则 $g'(r)=0$，该方法**局部收敛**。
>
> 再证明<u>二次收敛</u>：
>
> 对 $f(x)$ 在 $x_i$ 处进行 Taylor 展开：
> $$
> f(r)=f(x_i)+(r-x_i)f'(x_i)+\frac{1}{2}(r-x_i)^2f''(c_i)
> $$
> $c_i$ 在 $x_i$ 和 $r$ 之间。
>
> 由于 $f(r)=0,f'(r)\ne0$，则：
> $$
> \begin{align*}
> 0&=f(x_i)+(r-x_i)f'(x_i)+\frac{1}{2}(r-x_i)^2f''(c_i)\\
> x_i-\frac{f(x_i)}{f'(x_i)}&=r+\frac{1}{2}(r-x_i)^2\frac{f''(c_i)}{f'(x_i)}
> \end{align*}
> $$
> 这里假设 $f'(x_i)\ne0$。因此误差满足：
> $$
> \begin{align*}
> x_i-\frac{f(x_i)}{f'(x_i)}&=r+\frac{1}{2}(r-x_i)^2\frac{f''(c_i)}{f'(x_i)}\\
> x_i-\frac{f(x_i)}{f'(x_i)}-r&=\frac{1}{2}(r-x_i)^2\frac{f''(c_i)}{f'(x_i)}\\
> e_{i+1}&=\left|\frac{f''(c_i)}{2f'(x_i)}\right|e^2_i\\
> \lim_{i\to\infty}\frac{e_{i+1}}{e^2_i}&=\lim_{i\to\infty}\left|\frac{f''(c_i)}{2f'(x_i)}\right|=\frac{f''(r)}{2f'(r)}
> \end{align*}
> $$
> 由于收敛，$x_i$ 趋向于 $r$。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 不动点收敛加速
</blockquote>

在多重根处，牛顿法收敛速度为线性。设 $f\in C^{m+1}[a,b]$ 在 $r$ 处有个 $m$ 阶的多重根，则牛顿法局部收敛到 $r$，第 $i$ 步误差 $e_i$ 满足：
$$
\lim_{i\to\infty}\frac{e_{i+1}}{e_i}=\frac{m-1}{m}
$$
> 设此时的 $f(x)=(x-x^*)^mg(x)$，则
> $$
> \begin{gather*}
> x_{k+1}=x_k-\frac{(x_k-x^*)^mg(x_k)}{m(x_k-x^*)^{m-1}g(x_k)+(x_k-x^*)^mg'(x_k)}\\
> \frac{x_{k+1}-x^*}{x_k-x^*}=1-\frac{g(x_k)}{mg(x_k)+(x_k-x^*)g'(x_k)}\\
> \lim_{k\to\infty}\frac{x_{k+1}-x^*}{x_k-x^*}=1-\frac{1}{m}\\
> \end{gather*}
> $$

因此，我们可以改进在多重根处的牛顿法使其收敛速度更快：

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 不动点收敛加速
</blockquote>

(1) 设 $f\in C^{m+1}[a,b]$ 在 $r$ 处有个 $m$ 阶的多重根，则牛顿法局部收敛到 $r$，加速如下：
$$
x_{k+1}=x_k-m\frac{f(x_k)}{f'(x_k)}
$$
得到二阶收敛效果。

(2) $x^*$ 是 $f(x)$ 的 $m$ 重根，则是 $f'(x)$ 的 $m-1$ 重根：
$$
\begin{align*}
f'(x)&=m(x-x^*)^{m-1}g(x)+(x-x^*)^mg'(x)\\
&=(x-x^*)^{m-1}[mg(x)+(x-x^*)g'(x)]\\
&=(x-x^*)^{m-1}h(x)
\end{align*}
$$
则：
$$
\frac{f(x)}{f'(x)}=(x-x^*)\frac{g(x)}{h(x)}
$$
$x^*$ 是 $u(x)=f(x)/f'(x)$​ 的单根，然后迭代格式为：
$$
x_{k+1}=x_k-\frac{u(x_k)}{u'(x_k)}
$$

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 牛顿法可以不收敛
</blockquote>

设 $g(x)=x-f(x)/f'(x)$ 表示函数 $f$ 的牛顿法迭代，定义 $h(x)=g(g(x))$ 是牛顿法相邻两步的结果，$h'(x)=g'(g(x))g'(x)$。假设 $c$ 是 $h$ 的不动点，但 $c$ 不是 $g$ 的不动点，如果 $c$ 是 $f$ 的拐点，则不动点迭代 $h$ 局部收敛到 $c$。且对于初始估计 $c$，牛顿法不收敛到 $f$ 的根，而是震荡于 $c$ 和 $g(c)$ 两个值之间。

> 已知：
> $$
> h(c)=g(g(c))=c, g(c)\ne c,f''(c)=0
> $$
> 首先证明 $h$ 不动点迭代到 $c$：
> $$
> \begin{align*}
> h'(x)&=g'(g(x))g'(x)\\
> &=(\frac{f(g(x))f''(g(x))}{(f'(g(x)))^2})(\frac{f(x)f''(x)}{(f'(x))^2})\\
> \end{align*}
> $$
> 由于 $f''(c)=0$，因此 $h'(c)=0$，因此不动点迭代 $h$ 局部收敛到 $c$。
>
> 再证明以 $c$ 为初始估计的牛顿法迭代 $f$ 不收敛到根：
> $$
> x_1=c\\
> x_2=g(x_1)=g(c)\\
> x_3=g(x_2)=g(g(c))=h(c)=c
> $$
> 因为 $g(c)\ne c$，所以 $c$ 不是 $f$ 的根。
>
> 因此对于这样的函数 $f$，其牛顿法迭代不收敛，震荡在 $c$ 和 $g(c)$ 之间。

## 1.5 不需要导数的根求解

### 1.5.1 割线法 Secant Method

**`Condition:`** 闭区间上连续二阶可导。

```python
x[i + 1] = x[i] - (f(x[i]) * (x[i] - x[i - 1])) / (f(x[i]) - f(x[i - 1]))
```

该方法是牛顿法的变种，当迭代次数多的时候，$x_i-x_{i-1}$ 会足够接近，和导数值相近。
$$
\begin{gather*}
x_{i+1}=x_i-\frac{f(x_i)}{f'(x_i)}\\
f'(x_i)\approx\frac{f(x_i)-f(x_{i-1})}{x_i-x_{i-1}}\\
x_{i+1}=x_i-\frac{f(x_i)(x_i-x_{i-1})}{f(x_i)-f(x_{i-1})}=\frac{x_if(x_{i-1})-x_{i-1}f(x_i)}{f(x_i)-f(x_{i-1})}\\
\end{gather*}
$$

**`Theorem: `**令 $f$ 为二阶连续可微函数，如果割线法收敛到 $r$，即 $f(r)=0$，如果 $f'(r)\ne0,f''(r)\ne0$，则误差 $e_i$ 遵从：
$$
e_{i+1}\approx|\frac{f''(r)}{2f'(r)}|e_ie_{i-1}
$$
且如果：
$$
\lim_{i\to\infty}\frac{e_{i+1}}{e_i^{\alpha}}
$$
存在且不为 $0$，则 $\alpha=\dfrac{(1+\sqrt{5})}{2},e_{i+1}\approx|\dfrac{f''(r)}{2f'(r)}|^{\alpha-1}e_i^{\alpha}$。

### 1.5.2 试位方法 Regula Falsi

结合了割线法和二分法的一种方法。

```python
for _ in range(step):
    c = a - f(a) * (a - b) / (f(a) - f(b))
    if f(c) == 0:
        break
    if f(a) * f(c) < 0:
        b = c
    else:
        a = c
c if abs(fc) < abs(f((a + b) / 2)) else (a + b) / 2
```

试位方法一开始表现的比二分法和割线法都好，但是二分法可以保证消除 $1/2$ 的不确定性，试位方法却没有能力做出这样的保证，可能收敛更慢。

### 1.5.3 Muller 方法 Muller's Method

割线法的推广，用三个点得到一个形如 $y=f(x)$ 的抛物线，通过抛物线和坐标轴的交点确定新的点。

设三个估计点为 $x_0\ne x_1\ne x_2$，设抛物线为 $y=a(x-x_2)^2+b(x-x_2)+c$，则可得方程组
$$
\begin{cases}
a(x_0-x_2)^2+b(x_0-x_2)+c=f(x_0)\\
a(x_1-x_2)^2+b(x_1-x_2)+c=f(x_1)\\
c=f(x_2)\\
\end{cases}
$$
即：
$$
\begin{align*}
&\begin{cases}
a(x_0-x_2)^2+b(x_0-x_2)=f(x_0)-f(x_2)\\
a(x_1-x_2)^2+b(x_1-x_2)=f(x_1)-f(x_2)\\
\end{cases}\\
\end{align*}
$$

令：
$$
\begin{align*}
h_0&=x_1-x_0&h_1&=x_2-x_1\\
\delta_0&=\frac{f(x_1)-f(x_0)}{h_0}&\delta_1&=\frac{f(x_2)-f(x_1)}{h_1}
\end{align*}
$$
则：
$$
\begin{align*}
&\begin{cases}
a(x_0-x_2)^2+b(x_0-x_2)=f(x_0)-f(x_2)\\
a(x_1-x_2)^2+b(x_1-x_2)=f(x_1)-f(x_2)\\
\end{cases}\\
&\begin{cases}
a(h_0+h_1)^2&-b(h_0+h_1)&=h_0\delta_0+h_1\delta_1\\
a(h_1)^2&-b(h_1)&=h_1\delta_1\\
\end{cases}\\
\end{align*}
$$
因此解得：
$$
\begin{align*}
a&=\frac{\delta_1-\delta_0}{h_1+h_0}\\
b&=ah_1+\delta_1\\
c&=f(x_2)
\end{align*}
$$
考虑到该二次函数的零点和 $x_2$ 很接近，因此使用误差小的求根公式：
$$
x_3-x_2=\frac{-2c}{b\pm\sqrt{b^2-4ac}}\\
$$

上面的求根公式可能得到复数，这个复数也可以是原方程的解。如果得到两个实数，取和 $x_2$ 更接近的作为下一个估计点。

### 1.5.4 逆二次插值 Inverse Quadratic Interpolation

 逆二次插值类似于 Muller 和切线方法，但是用形如 $x=p(y)$ 形式的抛物线：
$$
\begin{gather*}
p(y)=\frac{a(y-B)(y-C)}{(A-B)(A-C)}+\frac{b(y-A)(y-C)}{(B-A)(B-C)}+\frac{c(y-A)(y-B)}{(C-A)(C-B)}\\
A=f(a),B=f(b),C=f(c)
\end{gather*}
$$
带入 $y=0$：
$$
\begin{align*}
x=p(0)&=\frac{aBC}{(A-B)(A-C)}+\frac{bAC}{(B-A)(B-C)}+\frac{cAB}{(C-A)(C-B)}\\
&=\frac{aBC(B-C)+bAC(C-A)+cAB(A-B)}{(A-B)(A-C)(B-C)}\\
&=-\frac{a(\frac{B}{A}-\frac{C}{A})+b(\frac{C}{B}-\frac{A}{B})+c(\frac{A}{C}-\frac{B}{C})}{(\frac{A}{B}-1)(\frac{C}{A}-1)(\frac{B}{C}-1)}\\
\end{align*}
$$

## 1.6 非线性方程组求解方法

非线性方程组：
$$
\begin{bmatrix}
f_1(\mathbf{x})\\
\vdots\\
f_n(\mathbf{x})\\
\end{bmatrix}=\mathbf{0}\quad\mathbf{x}=\begin{bmatrix}
x_1\\
\vdots\\
x_n\\
\end{bmatrix}
$$

### 1.6.1 不动点迭代

改写为 $\mathbf{x}=\Phi(\mathbf{x})$。构造迭代格式：
$$
\mathbf{x}_{k+1}=\Phi(\mathbf{x}_{k})
$$
如果向量序列收敛，则 $\mathbf{x}^*$ 是 $\Phi(\mathbf{x})$ 的一个不动点，也是原方程的解。

<blockquote style="border-left: 5px solid #4545aa; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(70, 70, 188, 0.1)">
    定理 不动点迭代
</blockquote>

设 $\Phi:D\subset\mathbb{R}^n\mapsto\mathbb{R}^n$ 在闭区域 $D_0\subset D$ 内满足：
$$
\|\Phi(\mathbf{x})-\Phi(\mathbf{y})\|\le L\|\mathbf{x}-\mathbf{y}\|\\
\forall{\mathbf{x}}\in D_0\to\Phi(\mathbf{x})\in D_0
$$
则：
$$
\forall{\mathbf{x}_0}\in D_0,\lim_{k\to\infty}\mathbf{x}^{(k)}=\mathbf{x}^*\\
\|\mathbf{x}^{(k)}-\mathbf{x}^*\|\le\frac{L}{1-L}\|\mathbf{x}^{(k)}-\mathbf{x}^{(k-1)}\|
$$
Jacobi $\nabla_{\Phi}(\mathbf{x})$ 代替 $L$。

### 1.6.2 Newton 迭代法

利用多元函数 Taylor 展开：
$$
F(\mathbf{x})=F(\mathbf{x}^{(k)})+\nabla{F}(\mathbf{x}^{(k)})(\mathbf{x}-\mathbf{x}^{(k)})
$$


迭代格式为：
$$
\mathbf{x}^{(k+1)}=\mathbf{x}^{(k)}-[\nabla{F}(\mathbf{x}^{(k)})]^{-1}F(\mathbf{x}^{(k)})
$$
