# Chapter 9 Random Numbers and Applications

## 9.4 Stochastic Differential Equations

### 9.4.1 Differential Equations with Noise

ODE 的解是函数，而 SDE 的解是随机过程。特别的，函数可以看作是方差为 0 的随机过程。考虑一维的 ODE：
$$
\mathrm{d} y = f(y, t) \mathrm{d} t
$$
对应的 SDE 可以看作是在每一个点都加入了噪声：
$$
\mathrm{d} y = f(y, t) \mathrm{d} t + g(y, t) \mathrm{d} W(t)
$$
其中 $W$ 是 Wiener 过程。等价的积分形式为：
$$
y(t) = y(0) + \int_{0}^{t} f(y(s), s) \mathrm{d} s + \int_{0}^{t} g(y(s), s) \mathrm{d} W(s)
$$
其中第二个积分是 Ito 积分，定义为：
$$
\int_{a}^{b} h(s) \mathrm{d} W(s) = \lim_{\Delta t \to 0} \sum_{i=1}^{n} h(t_{i-1}) (W(t_i) - W(t_{i-1}))
$$
其中 $a = t_0 < t_1 < \cdots < t_n = b$ 是一个分割，$\Delta t = \max_{i} (t_i - t_{i-1})$ 是分割的最大间隔。需要指出的是，相比于 Riemann 积分取的是分割内的任意点，Ito 积分取的是分割的左端点。此外，定义 Wiener 过程的微分 $\mathrm{d} W$ 为白噪声。

为了得到 SDE 的解析解，对于 $y = y(x, t)$，我们可以使用 Ito 公式：
$$
\mathrm{d} y = \frac{\partial y}{\partial x} \mathrm{d} x + \frac{\partial y}{\partial t} \mathrm{d} t + \frac{1}{2} \frac{\partial^2 y}{\partial x^2} (\mathrm{d} x)^2
$$
其中 $x = x(t)$ 是含有 $W(t)$ 的项，$(\mathrm{d} x)^2$ 可以由 $(\mathrm{d} t)^2 = 0$, $\mathrm{d} t \mathrm{d} W = \mathrm{d} W \mathrm{d} t = 0$, $\mathrm{d} W^2 = \mathrm{d} t$ 做进一步展开：
$$
\begin{align*}
(\mathrm{d} x)^2 &= \left(\frac{\partial x}{\partial t} \mathrm{d} t + \frac{\partial x}{\partial W} \mathrm{d} W\right)^2 \\
&= \left(\frac{\partial x}{\partial t}\right)^2 (\mathrm{d} t)^2 + 2 \frac{\partial^2 x}{\partial t \partial W} \mathrm{d} t \mathrm{d} W + \left(\frac{\partial x}{\partial W}\right)^2 (\mathrm{d} W)^2 \\
&= \left(\frac{\partial x}{\partial W}\right)^2 \mathrm{d} t
\end{align*}
$$

> 验证 $y(t) = y_0 \exp\left((r - \sigma^2 / 2)t + \sigma W\right)$ 是 SDE $\mathrm{d} y = r y \mathrm{d} t + \sigma y \mathrm{d} W(t)$ 的解。

令 $y(x, t) = y_0 \exp(x)$, $x = (r - \sigma^2 / 2) t + \sigma W$，则：
$$
\begin{align*}
\mathrm{d} x &= (r - \sigma^2 / 2) \mathrm{d} t + \sigma \mathrm{d} W \\
\mathrm{d} y &= y_0 \exp(x) \mathrm{d} x + 0 \mathrm{d} t + \frac{1}{2} y_0 \exp(x) (\mathrm{d} x)^2 \\
&= y (r - \sigma^2 / 2) \mathrm{d} t + y \sigma \mathrm{d} W + \frac{1}{2} y \sigma^2 \mathrm{d} t \\
&= r y \mathrm{d} t + \sigma y \mathrm{d} W
\end{align*}
$$

#### Exercises

> 1. 验证 $y(t) = t W(t) + c$ 是 SDE $\mathrm{d} y = W(t) \mathrm{d} t + t \mathrm{d} W$, $y(0) = c$ 的解。

令 $y(x, t) = tx + c$, $x = W$，则：
$$
\mathrm{d} y = x \mathrm{d} t + t \mathrm{d} x = W \mathrm{d} t + t \mathrm{d} W
$$

> 2. 验证 $y(t) = W(t)^2 - t + c$ 是 SDE $\mathrm{d} y = 2 W(t) \mathrm{d} W(t)$, $y(0) = c$ 的解。

令 $y(x, t) = x^2 - t + c$, $x = W$，则：
$$
\begin{align*}
\mathrm{d} y &= 2 x \mathrm{d} x - \mathrm{d} t + (\mathrm{d} x)^2 \\
&= 2 W \mathrm{d} W - \mathrm{d} t + (\mathrm{d} W)^2 \\
&= 2 W \mathrm{d} W
\end{align*}
$$

> 3. 验证 $y(t) = \ln (1 + W^2(t))$ 是 SDE $\mathrm{d} y = (1 - W^2(t)) \exp(-2y) \mathrm{d} t + 2 W(t) \exp(-y) \mathrm{d} W(t)$, $y(0) = 0$ 的解。

令 $y(x, t) = \ln (1 + x^2)$，$x = W$，则：$\exp(-y) = 1 / (1 + x^2)$，$\exp(-2y) = 1 / (1 + x^2)^2$，
$$
\begin{align*}
\mathrm{d} y &= \frac{2 x}{1 + x^2} \mathrm{d} x + 0 \mathrm{d} t + \frac{1}{2} \frac{2 - 4 x^2}{(1 + x^2)^2} (\mathrm{d} x)^2 \\
&= \frac{2 W}{1 + W^2} \mathrm{d} W + \frac{1 - 2 W^2}{(1 + W^2)^2} \mathrm{d} t \\
&= (1 - W^2) \exp(-2y) \mathrm{d} t + 2 W \exp(-y) \mathrm{d} W
\end{align*}
$$

> 4. 验证 $y(t) = W^3(t) / 3$ 是 SDE $\mathrm{d} y = W \mathrm{d} t + \sqrt[3]{9 y^2} \mathrm{d} W(t)$, $y(0) = 0$ 的解。

令 $y(x, t) = x^3 / 3$，$x = W$，则：$\sqrt[3]{9 y^2} = x^2$，
$$
\begin{align*}
\mathrm{d} y &= x^2 \mathrm{d} x + 0 \mathrm{d} t + \frac{1}{2} (2 x) (\mathrm{d} x)^2 \\
&= W \mathrm{d} t + \sqrt[3]{9 y^2} \mathrm{d} W
\end{align*}
$$

> 5. 验证 $y(t) = (1 + W) \exp(t^2 / 2)$ 是 SDE $\mathrm{d} y = ty \mathrm{d} t + \exp(t^2 / 2) \mathrm{d} W(t)$, $y(0) = 1$ 的解。

令 $y(x, t) = (1 + x) \exp(t^2 / 2)$，$x = W$，则：
$$
\begin{align*}
\mathrm{d} y &= \exp(t^2 / 2) \mathrm{d} x + t (1 + x) \exp(t^2 / 2) \mathrm{d} t + 0 (\mathrm{d} x)^2 \\
&= t y \mathrm{d} t + \exp(t^2 / 2) \mathrm{d} W
\end{align*}
$$

> 6. 验证 $y(t) = W^3 - 3 t W$ 是 SDE $\mathrm{d} y = 3(W^2 - t) \mathrm{d} W(t)$, $y(0) = 0$ 的解。

令 $y(x, t) = x^3 - 3 t x$，$x = W$，则：
$$
\begin{align*}
\mathrm{d} y &= (3 x^2 - 3 t) \mathrm{d} x - 3 x \mathrm{d} t + 3 x (\mathrm{d} x)^2 \\
&= 3 (W^2 - t) \mathrm{d} W
\end{align*}
$$

> 7. 验证 $y(t) = \sin(W(t))$ 是 SDE $\mathrm{d} y = - (y / 2) \mathrm{d} t + \sqrt{1 - y^2} \mathrm{d} W(t)$, $y(0) = 0$ 的解。

令 $y(x, t) = \sin(x)$，$x = W$，则：$\sqrt{1 - y^2} = \cos(x)$，
$$
\begin{align*}
\mathrm{d} y &= \cos(x) \mathrm{d} x + 0 \mathrm{d} t + \frac{1}{2} (-\sin(x)) (\mathrm{d} x)^2 \\
&= - \frac{y}{2} \mathrm{d} t + \sqrt{1 - y^2} \mathrm{d} W
\end{align*}
$$

> 8. 验证 $y(t) = \exp(W^2(t))$ 是 SDE $\mathrm{d} y = y(1 + 2\ln y) \mathrm{d}t + 2y W(t) \mathrm{d}W(t)$, $y(0) = 1$ 的解。

令 $y(x, t) = \exp(x^2)$，$x = W$，则：
$$
\begin{align*}
\mathrm{d} y &= 2 x \exp(x^2) \mathrm{d} x + 0 \mathrm{d} t + (1 + 2 x^2) \exp(x^2) (\mathrm{d} x)^2 \\
&= y(1 + 2 \ln y) \mathrm{d} t + 2 y W \mathrm{d} W
\end{align*}
$$

> 9. 验证 $y(t) = \ln(2W + \exp(y_0))$ 是 SDE $\mathrm{y} = -2 \exp(-2y) \mathrm{d}t + 2 \exp(-y) \mathrm{d}W(t)$, $y(0) = y_0$ 的解。

令 $y(x, t) = \ln(x)$，$x = 2W + \exp(y_0)$，则：$\exp(-y) = 1 / x$，$\exp(-2y) = 1 / x^2$，
$$
\begin{align*}
\mathrm{d} x &= 2 \mathrm{d} W \\
\mathrm{d} y &= \frac{1}{x} \mathrm{d} x + 0 \mathrm{d} t + \frac{1}{2} \left(-\frac{1}{x^2}\right) (\mathrm{d} x)^2 \\
&= \exp(-y) (2 \mathrm{d} W) - 2 \exp(-2y) \mathrm{d} t \\
&= -2 \exp(-2y) \mathrm{d} t + 2 \exp(-y) \mathrm{d} W
\end{align*}
$$

### 9.4.2 Numerical Methods for SDEs

用 Ito 积分的离散形式来得到 SDE 的数值解。对于在 $[a, b]$ 上的 SDE：
$$
\begin{cases}
\mathrm{d} y = f(y, t) \mathrm{d} t + g(y, t) \mathrm{d} W(t) \\
y(a) = y_a
\end{cases}
$$
类似 Euler 方法，我们可以使用 Euler-Maruyama 方法：
$$
\Delta y_{i} = f(y_{i}, t_{i}) \Delta t_i + g(y_{i}, t_{i}) \Delta W_{i}
$$

```python
def euler_maruyama(f, g, y0, t0, t1, n):
    dt = (t1 - t0) / n
    t = t0
    y = y0
    for i in range(n):
        dW = np.random.normal(0, 1) * np.sqrt(dt)
        y += f(y, t) * dt + g(y, t) * dW
        t += dt
    return y
```

数值方法的关键在于生成 $\Delta W$。由于 $W(t)$ 是 Wiener 过程，其增量 $W(t + \Delta t) - W(t)$ 是正态分布 $N(0, \Delta t)$ 的随机变量。因此，$\sqrt{\Delta t}$ 就是 $\Delta W$ 的标准差。

定义 SDE 求解器的**阶数**为 $m$，如果误差的绝对值的期望是步长 $\Delta t$ 的 $m$ 次方。即对于任意时刻 $t$，当 $\Delta t \to 0$ 时，有：
$$
\mathbb{E} \left| \hat{y}(t) - y(t) \right| = O((\Delta t)^m)
$$
其中 $\hat{y}(t)$ 是数值解，$y(t)$ 是解析解。Euler-Maruyama 方法是 1/2 阶的方法。

Milstein 方法是 Euler-Maruyama 方法的改进，其误差是 $\Delta t$ 的 1 阶。Milstein 方法把 $g(y, t)$ 的导数考虑进来：
$$
\Delta y_{i} = f(y_{i}, t_{i}) \Delta t + g(y_{i}, t_{i}) \Delta W_{i} + \frac{1}{2} g(y_{i}, t_{i}) \left(\frac{\partial g}{\partial y}(y_{i}, t_{i})\right) (\Delta W_{i}^2 - \Delta t_i)
$$

如果把 $g(y, t)$ 的导数用差分代替，就得到了一阶随机 Runge-Kutta 方法：
$$
\Delta y_{i} = f(y_{i}, t_{i}) \Delta t + g(y_{i}, t_{i}) \Delta W_{i} + \frac{1}{2} \cancel{g(y_{i}, t_{i})} \left(\frac{g(y_{i} + g(y_{i}, t_{i}) \sqrt{\Delta t}, t_{i}) - g(y_{i}, t_{i})}{\cancel{g(y_{i}, t_{i})}\sqrt{\Delta t}}\right) (\Delta W_{i}^2 - \Delta t_i)
$$

高阶 SDE 求解器相当复杂，一般不被采用。
