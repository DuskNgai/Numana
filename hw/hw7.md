## Question 1

> 已知常微分方程初值问题
> $$
> \begin{cases}
> y'(x)=x+1&(x>0)\\
> y(0)=0
> \end{cases}
> $$
> 的解析解为 $y(x)=0.5x^2+x$。
>
> (1) 分别导出求解上述问题的显式欧拉方法和隐式欧拉方法的近似解的表达式。
>
> (2) 求局部截断误差表达式。

(1) 显式欧拉方法：
$$
\begin{align*}
y_{n+1}&=y(x_n)+hf(x_n,y_n)\\
&=y_n+h(x_n+1)\\
&=y_n+h(nh+1)\\
&=y_n+nh^2+h\\
&=\frac{n(n-1)}{2}h^2+nh\\
&=0.5x_n^2+x_n-0.5x_nh
\end{align*}
$$
隐式欧拉方法：
$$
\begin{align*}
y_{n+1}&=y(x_n)+hf(x_{n+1},y_{n+1})\\
&=y_n+h(x_{n+1}+1)\\
&=y_n+h(nh+h+1)\\
&=y_n+nh^2+h^2+h\\
&=\frac{(n+1)n}{2}h^2+nh\\
&=0.5x_n^2+x_n+0.5x_nh
\end{align*}
$$
(2) 显式欧拉方法：
$$
\begin{align*}
T_{n+1}&=y(x_{n+1})-y_{n+1}=y(x_{n+1})-y(x_n)-hf(x_n,y_n)\\
&=y(x_{n+1})-y(x_n)-hf(x_n,y(x_n))=y(x_{n+1})-y(x_n)-hy'(x_n)\\
&=\frac{y''(\xi)}{2}h^2\\
&=\frac{1}{2}h^2
\end{align*}
$$
隐式欧拉方法：过程同显式：$T_{n+1}=-\frac{1}{2}h^2$。

## Question 2

> 考虑初值问题 
> $$
> \begin{cases}
> y'+y=0&(x>0)\\
> y(0)=1
> \end{cases}
> $$
> (1) 取 $h=0.1$，分别用显式欧拉公式、隐式欧拉公式和梯形公式计算 $y(0.2)$ 的近似值。
>
> (2) 证明：用梯形公式求上述问题的近似解为 $y_n=((2-h)/(2+h))^n$，并证明当 $h\to0$ 时，它收敛于原初值问题的准确解 $y=e^{-x}$。

(1) 显式欧拉公式：
$$
y_{n+1}=y(x_n)+hf(x_n,y_n)=(1-h)y_n
$$

$$
y_0=1\\y_1=0.9y_0=0.9\\y_2=0.9y_1=0.81\approx y(0.2)
$$

隐式欧拉公式：
$$
y_{n+1}=y(x_n)+hf(x_{n+1},y_{n+1})=y_n-hy_{n+1}\\
y_{n+1}=\frac{1}{1+h}y_n
$$

$$
y_0=1\\y_1=\frac{10}{11}y_0=\frac{10}{11}\\y_2=\frac{10}{11}y_1=\frac{100}{121}=0.8264\approx y(0.2)
$$

梯形公式：
$$
y_{n+1}=y(x_n)+(h/2)[f(x_{n},y_{n})+f(x_{n+1},y_{n+1})]=y_n+(h/2)[-y_n-y_{n+1}]\\
\frac{h+2}{2}y_{n+1}=\frac{2-h}{2}y_{n}\\
y_{n+1}=\frac{2-h}{2+h}y_{n}\\
$$

$$
y_0=1\\y_1=\frac{19}{21}y_0=\frac{19}{21}\\y_2=\frac{19}{21}y_1=\frac{361}{441}=0.8186\approx y(0.2)
$$

(2) 
$$
y_{n+1}=\frac{2-h}{2+h}y_{n}=\left(\frac{2-h}{2+h}\right)^{n+1}y_{0}=\left(\frac{2-h}{2+h}\right)^{n+1}\\
y_{n}=\left(\frac{2-h}{2+h}\right)^n
$$
由于 $x_n=nh$，带入有：
$$
\begin{align*}
y_n=\lim_{h\to0}\left(\frac{2-h}{2+h}\right)^n&=\lim_{n\to\infty}\left(\frac{2n-x_n}{2n+x_n}\right)^n\\
&=\lim_{n\to\infty}\left(1-\frac{x_n}{n+x_n/2}\right)^n\\
&=\lim_{n\to\infty}\left(1-\frac{x_n}{n+x_n/2}\right)^{\frac{n+x_n/2}{x_n}\frac{x_n}{n+x_n/2}n}\\
&=e^{-x_n}\\
\end{align*}
$$
即收敛到 $e^{-x}$。

## Question 3

> 给定求解常微分方程初值问题
> $$
> \begin{cases}
> y'=f(x,y)\\
> y(x_0)=y_0
> \end{cases}
> $$
> 的数值方法：
> $$
> \begin{cases}
> y_{n+1}=y_n+\dfrac{h}{4}(K_1+3K_2)\\
> K_1=f(x_n,y_n)\\
> K_2=f\left(x_{n}+\dfrac{2}{3}h,y_n+\dfrac{2}{3}hK_1\right)
> \end{cases}
> $$
> (1) 求局部截断误差的阶数。
>
> (2) 讨论用上述方法求解模型问题 $f(x,y)=\lambda y$ 稳定的步长 $h$ 的取值范围。

(1) Taylor 展开 $y(x_{n+1})$：
$$
\begin{align*}
y(x_{n+1})&=y(x_n+h)\\
&=y(x_{n})+hy'(x_n)+\frac{h^2}{2}y''(x_n)+O(h^3)
\end{align*}
$$
Taylor 展开 $K_2$：
$$
\begin{align*}
K_2&=f\left(x_{n}+\frac{2}{3}h,y_n+\frac{2}{3}hK_1\right)\\
&=f(x_n,y_n)+\frac{2}{3}h[f'_x(x_n,y_n)+f'_y(x_n,y_n)y'(x_n)]+O(h^2)\\
&=f(x_n,y(x_n))+\frac{2}{3}h[f'_x(x_n,y(x_n))+f'_y(x_n,y(x_n))y'(x_n)]+O(h^2)\\
&=y'(x_n)+\frac{2}{3}hy''(x_n)+O(h^2)\\
\end{align*}
$$
因此：
$$
\begin{align*}
y_{n+1}&=y_n+\frac{h}{4}(K_1+3K_2)\\
&=y_n+\frac{h}{4}(4y'(x_n)+2hy''(x_n)+O(h^2))\\
&=y_n+hy'(x_n)+\frac{h^2}{2}y''(x_n)+O(h^3)\\
\end{align*}
$$
因此：
$$
T_{n+1}=y(x_{n+1})-y_{n+1}=O(h^3)
$$
局部截断误差的阶数 2 阶。

(2) 设 $\epsilon_{n}=y(x_n)-y_n$：
$$
y_{n+1}=y_n+\frac{h}{4}(4\lambda y_n+2h\lambda^2y_n)\\
\epsilon_{n+1}=\epsilon_n(1+h\lambda+\frac{1}{2}h^2\lambda^2)
$$
即：
$$
|1+h\lambda+\frac{1}{2}h^2\lambda^2|\le1\\
-2\le h\lambda\le0\\
0\le h\le-\frac{2}{\lambda}
$$

## Question 4

> 对于求解常微分方程初值问题的单步方法：
> $$
> y_{n+1}=y_{n}+\frac{h}{8}\left[3f(x_n,y_n)+5f\left(x_n+\frac{4}{5}h,y_n+\frac{4}{5}hf(x_n,y_n)\right)\right]
> $$
> (1) 求局部截断误差 $T_{n+1}$ 的主部。
>
> (2) 判断此方法是无条件稳定还是有条件稳定？若是有条件稳定，求出步长 $h$ 的取值范围。

(1) Taylor 展开 $y(x_{n+1})$：
$$
\begin{align*}
y(x_{n+1})&=y(x_n+h)\\
&=y(x_{n})+hy'(x_n)+\frac{h^2}{2}y''(x_n)+\frac{h^3}{6}y'''(x_n)+O(h^4)
\end{align*}
$$
Taylor 展开 $K_2$：
$$
\begin{align*}
K_2&=f\left(x_n+\frac{4}{5}h,y_n+\frac{4}{5}hf(x_n,y_n)\right)\\
&=f(x_n,y_n)+\left(\frac{4}{5}h\right)[f'_x(x_n,y_n)+f'_y(x_n,y_n)y'(x_n)]\\
&+\left(\frac{4}{5}h\right)^2[f''_{xx}(x_n,y_n)+2f''_{xy}(x_n,y_n)y'(x_n)+f''_{yy}(x_n,y_n)(y'(x_n))^2]+O(h^3)\\
&=y'(x_n)+\left(\frac{4}{5}h\right)y''(x_n)+\left(\frac{4}{5}h\right)^2[y'''(x_n)-f'_y(x_n,y(x_n))y''(x_n)]+O(h^3)\\
\end{align*}
$$
因此：
$$
\begin{align*}
y_{n+1}&=y_n+\frac{h}{8}(3K_1+5K_2)\\
&=y_n+\frac{h}{8}\left(8y'(x_n)+4hy''(x_n)+\frac{8h^2}{5}[y'''(x_n)-f'_y(x_n,y(x_n))y''(x_n)]+O(h^3)\right)\\
&=y_n+hy'(x_n)+\frac{h^2}{2}y''(x_n)+\frac{h^3}{5}[y'''(x_n)-f'_y(x_n,y(x_n))y''(x_n)]+O(h^4)\\
\end{align*}
$$
局部截断误差 $T_{n+1}$ 的为：
$$
\begin{align*}
T_{n+1}&=\frac{h^3}{6}y'''(x_n)-\frac{h^3}{5}[y'''(x_n)-f'_y(x_n,y(x_n))y''(x_n)]+O(h^4)\\
&=\left[\frac{1}{5}f'_y(x_n,y(x_n))y''(x_n)-\frac{1}{30}y'''(x_n)\right]h^3+O(h^4)
\end{align*}
$$
主部为 $\left[\frac{1}{5}f'_y(x_n,y(x_n))y''(x_n)-\frac{1}{30}y'''(x_n)\right]h^3$。

(2) 设 $\epsilon_{n}=y(x_n)-y_n$：
$$
y_{n+1}=y_n+\frac{h}{8}(8\lambda y_n+4h\lambda^2y_n)\\
\epsilon_{n+1}=\epsilon_n(1+h\lambda+\frac{1}{2}h^2\lambda^2)
$$
因此为有条件稳定，步长 $0\le h\le-2/\lambda,(\lambda<0)$。

## Question 5

> 证明解 $y'=f(x,y)$ 的差分格式：
> $$
> y_{n+1}=\frac{1}{2}(y_n+y_{n-1})+\frac{h}{4}(4y'_{n+1}-y'_n+3y'_{n-1})
> $$
> 是二阶的，并求出截断误差的首项。

Taylor 展开 $y(x_{n+1})$：
$$
\begin{align*}
y(x_{n+1})&=y(x_n+h)\\
&=y(x_{n})+hy'(x_n)+\frac{h^2}{2}y''(x_n)+\frac{h^3}{6}y'''(x_n)+O(h^4)
\end{align*}
$$
Taylor 展开 $y_{n-1}$​：
$$
\begin{align*}
y_{n-1}&=y(x_{n-1})=y(x_{n}-h)\\
&=y(x_n)-hy'(x_n)+\frac{h^2}{2}y''(x_n)-\frac{h^3}{6}y'''(x_n)+O(h^4)
\end{align*}
$$
Taylor 展开 $y'_{n+1}$：
$$
\begin{align*}
y'_{n+1}&=y'(x_{n+1})=y'(x_{n}+h)\\
&=y'(x_n)+hy''(x_n)+\frac{h^2}{2}y'''(x_n)+O(h^3)
\end{align*}
$$
Taylor 展开 $y'_{n-1}$：
$$
\begin{align*}
y'_{n-1}&=y'(x_{n-1})=y'(x_{n}-h)\\
&=y'(x_n)-hy''(x_n)+\frac{h^2}{2}y'''(x_n)+O(h^3)
\end{align*}
$$
合并得到：
$$
\begin{align*}
y_{n+1}&=\frac{1}{2}(y_n+y_{n-1})+\frac{h}{4}(4y'_{n+1}-y'_n+3y'_{n-1})\\
&=\frac{1}{2}(2y_n-hy'(x_n)+\frac{h^2}{2}y''(x_n)-\frac{h^3}{6}y'''(x_n)+O(h^4))\\
&+\frac{h}{4}(6y'(x_n)+hy''(x_n)+\frac{7h^2}{2}y'''(x_n)+O(h^3))\\
&=y(x_n)+hy'(x_n)+\frac{h^2}{2}y''(x_n)+\frac{19h^3}{24}y'''(x_n)+O(h^4)
\end{align*}
$$
局部截断误差 $T_{n+1}$ 的为：
$$
\begin{align*}
T_{n+1}&=\frac{h^3}{6}y'''(x_n)-\frac{19h^3}{24}y'''(x_n)+O(h^4)\\
&=\frac{5}{8}y'''(x_n)h^3+O(h^4)
\end{align*}
$$
因此是二阶的，首项为 $\frac{5}{8}y'''(x_n)h^3$。

## Question 6

> 试确定常数 $a,b,c$ 的值，使得求解常微分方程初值问题的单步方法：
> $$
> y_{n+1}=y_n+\frac{h}{8}\left[2f(x_n,y_n)+af(x_n+bh,y_n+chf(x_n,y_n))\right]
> $$
> 的局部截断误差为 $T_{n+1}=O(h^3)$。

Taylor 展开 $y(x_{n+1})$：
$$
\begin{align*}
y(x_{n+1})&=y(x_n+h)\\
&=y(x_{n})+hy'(x_n)+\frac{h^2}{2}y''(x_n)+\frac{h^3}{6}y'''(x_n)+O(h^4)
\end{align*}
$$
Taylor 展开 $K_2$：
$$
\begin{align*}
K_2&=f(x_n+bh,y_n+chf(x_n,y_n))\\
&=f(x_n,y_n)+h[bf'_x(x_n,y_n)+cf'_y(x_n,y_n)y'(x_n)]+O(h^2)\\
\end{align*}
$$
因此：
$$
\begin{align*}
y_{n+1}&=y_n+\frac{h}{8}(2K_1+aK_2)\\
&=y_n+\frac{h}{8}\left((2+a)y'(x_n)+ah[bf'_x(x_n,y_n)+cf'_y(x_n,y_n)y'(x_n)]+O(h^2)\right)\\
\end{align*}
$$
为使得局部截断误差 $T_{n+1}=O(h^3)$，可得到：
$$
a=6,b=c=\frac{2}{3}
$$

## Question 7

> 给定求解常微分方程初值问题的线性多步法
> $$
> y_{n+2}=y_{n+1}+h[\alpha_1f_n+\alpha_2f_{n+1}]
> $$
> 确定常数 $\alpha_1,\alpha_2$ 使得上述方法的局部截断误差阶数最高，并求局部截断误差的主部。

化为
$$
y_{n+1}=y_{n}+h[\alpha_1f_{n-1}+\alpha_2f_{n}]
$$
Taylor 展开 $y(x_{n+1})$：
$$
\begin{align*}
y(x_{n+1})&=y(x_n+h)\\
&=y(x_{n})+hy'(x_n)+\frac{h^2}{2}y''(x_n)+\frac{h^3}{6}y'''(x_n)+O(h^4)
\end{align*}
$$
Taylor 展开 $y'(x_{n-1})$​：
$$
\begin{align*}
y'(x_{n-1})&=y'(x_{n}-h)\\
&=y'(x_n)-hy''(x_n)+\frac{h^2}{2}y'''(x_n)+O(h^3)
\end{align*}
$$
带入：
$$
\begin{align*}
&\quad\ \ y(x_n)+h[\alpha_1[y'(x_n)-hy''(x_n)+\frac{h^2}{2}y'''(x_n)]+\alpha_2y'(x_n)+O(h^3)]\\
&=y(x_n)+h[(\alpha_1+\alpha_2)y'(x_n)-\alpha_1hy''(x_n)+\alpha_1\frac{h^2}{2}y'''(x_n)+O(h^3)]\\
\end{align*}
$$
即：
$$
\alpha_1+\alpha_2=1,-\alpha_1=\frac{1}{2}
$$
即：
$$
\alpha_1=-\frac{1}{2},\alpha_2=\frac{3}{2}
$$
此时：
$$
T_{n+1}=\frac{5h^3}{12}y'''(x_n)+O(h^4)
$$

## Question 8

> 给定求解常微分方程组初值问题 $\mathbf{y}=f(x,\mathbf{y}),\mathbf{y}(x_0)=\mathbf{y}_0$​ 的单步方法：
> $$
> \mathbf{y}_{n+1}=\mathbf{y}_n+h\phi(x_n,\mathbf{y}_n,\mathbf{y}_{n+1})
> $$
> 其中，$\mathbf{y}(x)=[y_1(x)\quad y_2(x)\quad\dots\quad y_m(x)]^T$。
>
> 如果 $\|\mathbf{y}_{n+1}(x)\|^2=\dots=\|\mathbf{y}_{1}(x)\|^2=\|\mathbf{y}_{0}(x)\|^2$，则称单步方法具有平方守恒律。
>
> (1) 证明：用梯形公式解二阶方程初始值问题
> $$
> \begin{cases}
> y''(x)+y(x)=0\\
> y(x_0)=y_0,y'(x_0)=y'_0
> \end{cases}
> $$
> 具有平方守恒律。
>
> (2) 若用梯形公式解一阶方程组
> $$
> \begin{cases}
> \mathbf{y}'=A\mathbf{y}\\
> \mathbf{y}(x_0)=\mathbf{y}_0
> \end{cases}
> $$
> 所导出的公式具有平方守恒律，试给出矩阵 $A$ 满足的一个充分条件，并证明你的结论。

(1) 改写得到
$$
\begin{cases}
y_1=y(x)\\
y_2=y'(x)\\
y_3=y''(x)
\end{cases}\Longrightarrow
\begin{cases}
y_1'=y_2\\
y_2'=y_3\\
y_3=-y_1\\
y_1(x_0)=y_0,y_2(x_0)=y'_0
\end{cases}
$$
应用梯形公式：
$$
\begin{cases}
y_{1,n+1}=y_{1,n}+\dfrac{h}{2}[y_{2,n}+y_{2,n+1}]\\
y_{2,n+1}=y_{2,n}+\dfrac{h}{2}[y_{3,n}+y_{3,n+1}]=y_{2,n}-\dfrac{h}{2}[y_{1,n}+y_{1,n+1}]\\
\end{cases}
$$
即：
$$
\begin{cases}
y_{1,n+1}-y_{1,n}=\dfrac{h}{2}[y_{2,n}+y_{2,n+1}]\\
y_{2,n+1}-y_{2,n}=-\dfrac{h}{2}[y_{1,n}+y_{1,n+1}]\\
\end{cases}
$$
乘以系数：
$$
\begin{cases}
(y_{1,n+1}-y_{1,n})(y_{1,n}+y_{1,n+1})=\dfrac{h}{2}(y_{2,n}+y_{2,n+1})(y_{1,n}+y_{1,n+1})\\
(y_{2,n+1}-y_{2,n})(y_{2,n}+y_{2,n+1})=-\dfrac{h}{2}(y_{1,n}+y_{1,n+1})(y_{2,n}+y_{2,n+1})\\
\end{cases}
$$
相加得到：
$$
y_{1,n+1}^2+y_{2,n+1}^2-y_{2,n}^2-y_{1,n}^2=0
$$
即：
$$
y_{1,n+1}^2+y_{2,n+1}^2=y_{1,n}^2+y_{2,n}^2
$$
因此具有平方守恒律。

(2)

应用梯形公式：
$$
\begin{align*}
\mathbf{y}_{n+1}&=\mathbf{y}_n+\frac{h}{2}A(\mathbf{y}_n+\mathbf{y}_{n+1})\\
\mathbf{y}_{n+1}&=\left(I-\frac{h}{2}A\right)^{-1}\left(I+\frac{h}{2}A\right)\mathbf{y}_n
\end{align*}
$$
要使 $\|\|^2=\|\mathbf{y}_{n}(x)\|^2$，则
$$
\mathbf{y}_{n+1}(x)^T\mathbf{y}_{n+1}(x)=\mathbf{y}_{n}(x)^T\mathbf{y}_{n}(x)
$$
即
$$
\left(I+\frac{h}{2}A\right)^T\left(I-\frac{h}{2}A\right)^{-T}\left(I-\frac{h}{2}A\right)^{-1}\left(I+\frac{h}{2}A\right)=I\\
(I+tA)^T(I-tA)^{-T}(I-tA)^{-1}(I+tA)=I\\
$$
因此 $A$ 为反对称矩阵。
