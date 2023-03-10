# Chapter 5 数值微分和积分 Numerical Differentiation and Integration

## 5.1 Numerical Differentiation 数值微分

<blockquote style="border-left: 5px solid #42b983; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(66, 185, 131, 0.1)">
    定理 推广中值定理
</blockquote>

设 $f\in C[a,b]$，$x_1,\dots,x_n\in[a,b]$，$a_1,\dots,a_n>0$，则 $\exist c\in[a,b]$，满足：
$$
(a_1+\dots+a_n)f(c)=a_1f(x_1)+\dots+a_nf(x_n)
$$

> 设 $f(x_i),f(x_j)$ 是 $f(x_1),\dots,f(x_n)$ 之间的最小和最大值，则必然有：
> $$
> \sum_{k=1}^{n}a_kf(x_i)\le\sum_{k=1}^{n}a_kf(x_k)\le\sum_{k=1}^{n}a_kf(x_j)\\
> f(x_i)\le\frac{\sum_{k=1}^{n}a_kf(x_k)}{\sum_{k=1}^{n}a_k}\le f(x_j)
> $$
> 根据介值定理，$f(x)$ 能在 $x_i$ 与 $x_j$ 之间取到 $\xi$，使得
> $$
> \color{red}f(\xi)=\frac{\sum_{k=1}^{n}a_kf(x_k)}{\sum_{k=1}^{n}a_k}
> $$

### 5.1.1 Finite Difference Formulas

设 $f\in C^{n+1}$，对其在 $x+h$ 处 Taylor 展开有：
$$
f(x+h)=f(x)+hf'(x)+\frac{h^2}{2!}f''(x)+\dots+\frac{h^n}{n!}f^{(n)}(x)+\frac{h^{n+1}}{(n+1)!}f^{(n+1)}(\xi)
$$
$\xi$ 位于 $x$ 与 $x+h$ 之间。

**$n$ 阶近似：误差是 $O(h^n)$。**

> Two-point forward-difference formula 二点前向差分公式
>
> 设 $f\in C^{2}$，
> $$
> \color{red}f'(x)=\frac{f(x+h)-f(x)}{h}-\frac{h}{2}f''(\xi)
> $$
> $\xi$ 位于 $x$ 与 $x+h$ 之间，且依赖于 $h$。该公式是一阶近似。

> Three-point centered-difference formula 三点中心差分公式
>
> 设 $f\in C^{3}$，由
> $$
> f(x+h)=f(x)+hf'(x)+\frac{h^2}{2!}f''(x)+\frac{h^3}{3!}f'''(\xi_1)\\
> f(x-h)=f(x)-hf'(x)+\frac{h^2}{2!}f''(x)-\frac{h^3}{3!}f'''(\xi_2)
> $$
> 相减可得：
> $$
> f'(x)=\frac{f(x+h)-f(x-h)}{2h}-\frac{h^2}{12}f'''(\xi_1)-\frac{h^2}{12}f'''(\xi_2)
> $$
> $x-h<\xi_2<x<\xi_1<x+h$。
>
> 由推广中值定理可得到：
> $$
> \color{red}f'(x)=\frac{f(x+h)-f(x-h)}{2h}-\frac{h^2}{6}f'''(\xi)
> $$
> 其中 $x-h<\xi<x+h$。该公式是二阶近似。

> Three-point centered-difference formula for second derivative 二阶导数的三点中心差分公式
>
> 设 $f\in C^{4}$，由
> $$
> f(x+h)=f(x)+hf'(x)+\frac{h^2}{2!}f''(x)+\frac{h^3}{3!}f'''(x)+\frac{h^4}{4!}f^{(4)}(\xi_1)\\
> f(x-h)=f(x)-hf'(x)+\frac{h^2}{2!}f''(x)-\frac{h^3}{3!}f'''(x)+\frac{h^4}{4!}f^{(4)}(\xi_2)
> $$
> 相加，并应用推广中值定理可得：
> $$
> \color{red}f''(x)=\frac{f(x+h)-2f(x)+f(x-h)}{h^2}-\frac{h^2}{12}f^{(4)}(\xi)
> $$
> 其中 $x-h<\xi<x+h$。该公式是二阶近似。

### 5.1.2 Rouding Error 舍入误差

分析三点中心差分公式的舍入误差：

设 $f(x+h)$ 在机器中的存储的形式为 $\hat{f}(x+h)$，且误差是与机器精度同阶的，即 $\hat{f}(x+h)=f(x+h)+\epsilon_1$ 因此：
$$
\begin{align*}
f'(x)_{\text{correct}}-f'(x)_{\text{machine}}&=f'(x)_{\text{correct}}-\frac{\hat{f}(x+h)-\hat{f}(x-h)}{2h}\\
&=f'(x)_{\text{correct}}-\frac{f(x+h)+\epsilon_1-f(x-h)-\epsilon_2}{2h}\\
&=\left[f'(x)_{\text{correct}}-\frac{f(x+h)-f(x-h)}{2h}\right]+\frac{\epsilon_2-\epsilon_1}{2h}\\
&=\left[f'(x)_{\text{correct}}-f'(x)_{\text{formula}}\right]+\text{error}_{\text{rounding}}
\end{align*}
$$
因此误差的绝对值的上界就是：
$$
\begin{align*}
|f'(x)_{\text{correct}}-f'(x)_{\text{machine}}|&=|\left(f'(x)_{\text{correct}}-f'(x)_{\text{formula}}\right)+\text{error}_{\text{rounding}}|\\
&\le|f'(x)_{\text{correct}}-f'(x)_{\text{formula}}|+\left|\frac{\epsilon_2-\epsilon_1}{2h}\right|\\
&\le\frac{h^2}{6}f'''(\xi)+\frac{2\epsilon}{2h}\\
&=\frac{h^2}{6}f'''(\xi)+\frac{\epsilon}{h}
\end{align*}
$$
其中 $\epsilon$ 为机器精度，$x-h<\xi<x+h$。

### 5.1.3 Richardson Extrapolation 外推法

假设 $F(h)$ 逼近一个给定的量 $F^*$ ($F^*$ 与 $h$ 无关) ，其余项为：
$$
F^*-F(h)=\sum_{k=1}^{\infty}\alpha_kh^{p_k}
$$
其中 $p_k,\alpha_k$ 为与 $h$ 无关的常数，$\alpha_k\ne0,k>0$。定义外推公式：
$$
\begin{align*}
F_1(h)&=F(h)\\
F_{m+1}(h)&=\frac{F_m(qh)-q^{p_m}F_m(h)}{1-q^{p_m}},q\in(0,1)
\end{align*}
$$

> 以 $q=1/2$ 为例，由于
> $$
> F^*-F(h)=\sum_{k=1}^{\infty}\alpha_kh^{p_k}
> $$
>
> 步长减半：
> $$
> F^*-F\left(\frac{h}{2}\right)=\sum_{k=1}^{\infty}\alpha_k\left(\frac{h}{2}\right)^{p_k}
> $$
> 提取右端第一项并整理形式：
> $$
> 2^{p_1}F^*-2^{p_1}F\left(\frac{h}{2}\right)=\alpha_1h^{p_1}+\sum_{k=2}^{\infty}\alpha_k\frac{h^{p_k}}{2^{p_k-p_1}}=\alpha_1h^{p_1}+\sum_{k=2}^{\infty}\frac{\alpha_k}{2^{p_k-p_1}}h^{p_k}
> $$
> 减去最上面的式子，消去右端第一项：
> $$
> (2^{p_1}-1)F^*-\left[2^{p_1}F\left(\frac{h}{2}\right)-F(h)\right]=\sum_{k=2}^{\infty}\frac{\alpha_k}{2^{p_k-p_1}}h^{p_k}\\
> F^*-\frac{2^{p_1}F\left(\frac{h}{2}\right)-F(h)}{(2^{p_1}-1)}=\sum_{k=1}^{\infty}\beta_kh^{p_k}\\
> F^*-\frac{F\left(\frac{h}{2}\right)-\left(\frac{1}{2}\right)^{p_1}F(h)}{(1-\left(\frac{1}{2}\right)^{p_1})}=\sum_{k=1}^{\infty}\beta_kh^{p_k}
> $$
> 可以发现误差的形式保持。

<blockquote style="border-left: 5px solid #b94263; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(185, 66, 110, 0.1)">
    例题
</blockquote>

设 $I(h)$ 逼近 $I$ 的余项为:
$$
I-I(h)=\alpha_1h^2+\alpha_2h^4+\cdots
$$
其中 $\alpha_k$ 为与 $h$ 无关的常数，且 $I(0.1)=0.99,I(0.04)=1.20$，在此基础上用外推法得到 $I$ 的更高精度的近似值为多少？

> $$
> I-I(0.1)=\alpha_10.1^2+\alpha_20.1^4+\cdots\\
> I-I(0.04)=\alpha_10.04^2+\alpha_20.04^4+\cdots
> $$
>
> 消去 $\alpha_1$:
> $$
> 5.25I-6.25I(0.04)+I(0.1)=0.1^2(0.1^2-0.04^2)\alpha_2+\cdots\\
> \frac{6.25I(0.04)-I(0.1)}{5.25}=1.19238
> $$

<blockquote style="border-left: 5px solid #b94263; border-radius: 3px 0 0 3px; padding: 10px 15px; background-color: rgba(185, 66, 110, 0.1)">
    例题
</blockquote>

外推二阶精度的三点中心差分公式。

> 已知：
> $$
> F(h)=f'(x)=\frac{f(x+h)-f(x-h)}{2h}-\frac{h^2}{6}f'''(\xi)
> $$
> 步长减半：
> $$
> F(h/2)=f'(x)=\frac{f(x+h/2)-f(x-h/2)}{h}-\frac{h^2}{24}f'''(\xi)
> $$
> 消去三阶导数：
> $$
> \begin{align*}
> f'(x)&=\frac{2^2F(h/2)-F(h)}{2^2-1}\\
> &=\frac{-f(x+h)+8f(x+h/2)-8f(x-h/2)+f(x-h)}{6h}
> \end{align*}
> $$
> 由外推法可知该公式的精度至少是 3 阶。

### 5.1 Exercises

#### 8

证明：一阶导数的二阶公式：
$$
f'(x)=\frac{-f(x+2h)+4f(x+h)-3f(x)}{2h}+O(h^2)
$$

> $$
> f(x+h)=f(x)+hf'(x)+\frac{h^2}{2}f''(x)+\frac{h^3}{6}f'''(\xi_1)\\
> f(x+2h)=f(x)+2hf'(x)+2h^2f''(x)+\frac{4h^3}{3}f'''(\xi_2)\\
> $$
>
> 因此：
> $$
> \frac{-f(x+2h)+4f(x+h)-3f(x)}{2h}=f'(x)-\frac{1}{3}h^2\left[2f''(\xi_2)-f''(\xi_1)\right]
> $$



