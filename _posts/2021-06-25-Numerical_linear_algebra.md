---
layout: post
title:  Numerical linear algebra
date: 2021-06-25
description: "数值代数考试复习"
tag: learning
---

### 1. Gram-Schmidt正交化的修正策略

(标准 Gram-Schmidt 正交化) 给定 $n \times m(m \leqslant$$n)$ 阶列满秩矩阵 $\boldsymbol{X}=\left[\boldsymbol{x}_ {1}, \boldsymbol{x}_ {2}, \cdots, \boldsymbol{x}_ {m}\right]$, 本算法产 生 $n \times m$ 阶正交矩阵 $\boldsymbol{Q}=\left[\boldsymbol{q}_ {1}, \boldsymbol{q}_ {2}, \cdots, \boldsymbol{q}_ {m}\right]$和 $m \times m$ 阶非奇异上三角矩阵 $\boldsymbol{R}=\left(r_ {i k}\right)$.

$r_ {11}=\left\|\boldsymbol{x}_ {1}\right\|_ {2} ; \boldsymbol{q}_ {1}=\boldsymbol{x}_ {1} / r_ {11}$

for $k=2: m$  
$\begin{aligned}
    &\text{ for }(i=1: k-1), r_ {i k}=\left(\boldsymbol{x}_ {k}, \boldsymbol{q}_ {i}\right) ; \text { end } \\
    &\widetilde{\boldsymbol{q}}=\boldsymbol{x}_ {k}-\sum_ {i=1}^{k-1} r_ {i k} \boldsymbol{q}_ {i} \\
    &r_ {k k}=\|\widetilde{\boldsymbol{q}}\|_ {2} ; \boldsymbol{q}_ {k}=\widetilde{\boldsymbol{q}} / r_ {k k}
\end{aligned}$  
end

可以看出，当 $\boldsymbol{x}_ {k}$ 与前面的 $\boldsymbol{x}_ {i}$ 接近线性相关时，$r_ {k k}$将接近于 0，用它作分母会导致很大的舍入误差 ，$\boldsymbol{q}_ {i}$ 之间很快失去了正交性。一种改进就是修正的 Gram-Schmidt 交化。它把算法的第 3 行和第 4 行修改为下面的循环:

$$
\widetilde{\boldsymbol{q}}=\boldsymbol{x}_ {k}, \quad r_ {i k}=\left(\widetilde{\boldsymbol{q}}, \boldsymbol{q}_ {i}\right), \quad \widetilde{\boldsymbol{q}}:=\widetilde{\boldsymbol{q}}-r_ {i k} \boldsymbol{q}_ {i}, \quad i=1,2, \cdots, k-1
$$

在不考虑舍入误差时，上述计算结果 $\widetilde{\boldsymbol{q}}$ 与算法 中第 4 行的 $\widetilde{\boldsymbol{q}}$ 是相同的。修正的 Gram-Schmidt 正交化利用了最新的 $\widetilde{\boldsymbol{q}}$，有助于减少舍入误差影响。

### 2. Arnoldi分解和Lanczos分解的区别

给定矩阵$\boldsymbol{A} \in \mathbb{R}^{n\times n}$ 和向量 $\boldsymbol{v} \in \mathbb{R}^n$ , 设 $\boldsymbol{K}_ {k+1}(\boldsymbol{A}, \boldsymbol{v})=\left[\boldsymbol{v}, \boldsymbol{A} \boldsymbol{v}, \cdots, \boldsymbol{A}^{k} \boldsymbol{v}\right] 
$是列满秩的，并假定它的
QR 分解为

$$
\boldsymbol{K}_ {k+1}(\boldsymbol{A}, \boldsymbol{v})=\boldsymbol{V}_ {k+1} \boldsymbol{R}_ {k+1}
$$

式中: $\boldsymbol{V}_ {k+1}=\left[\boldsymbol{v}_ {1}, \boldsymbol{v}_ {2}, \cdots, \boldsymbol{v}_ {k+1}\right] \in \mathbb{R}^{n \times(k+1)}$ 满足 $\boldsymbol{V}_ {k+1}^{\mathrm{T}} \boldsymbol{V}_ {k+1}=$$\boldsymbol{I}_ {k+1}$，$\boldsymbol{R}_ {k+1}$ 为非奇异的上三角矩阵。

 Arnoldi 分解式为

$$
\boldsymbol{A} \boldsymbol{V}_ {k}=\boldsymbol{V}_ {k} \boldsymbol{H}_ {k}+\beta_ {k} \boldsymbol{v}_ {k+1} \boldsymbol{e}_ {k}^{\mathrm{T}}
$$

其中 $\boldsymbol{H}_k$ 为 $(k+1)\times k$ 阶的上 Hessenberg 矩阵。

当 $A \in \mathbb{R}^{n \times n}$ 是对称矩阵时,  对应的 Arnoldi 分解就变成

$$
\boldsymbol{A} \boldsymbol{V}_ {k}=\boldsymbol{V}_ {k} \boldsymbol{T}_ {k}+\beta_ {k} \boldsymbol{v}_ {k+1} \boldsymbol{e}_ {k}^{\mathrm{T}}
$$

式中:

$$
\boldsymbol{T}_ {k}=\left[\begin{array}{cccc}
\alpha_ {1} & \beta_ {1} & & \\
\beta_ {1} & \alpha_ {2} & \ddots & \\
& \ddots & \ddots & \beta_ {k-1} \\
& & \beta_ {k-1} & \alpha_ {k}
\end{array}\right] \text { . }
$$

称上式为一个长度为 $k$ 的 Lanczos 分解。因此，Lanczos 分解是 Arnoldi 的特殊情形。

### 3. 投影方法的基本框架

投影方法的基本框架为：
step 1, 选取一对子空间 $\mathcal{K}$ 和 $\mathcal{L} .$  
step 2, 选取 $\mathcal{K}$ 和 $\mathcal{L}$ 的基 $\boldsymbol{V}=\left[\boldsymbol{v}_ {1}, \cdots, \boldsymbol{v}_ {m}\right]$ 和 $\boldsymbol{W}=\left[\boldsymbol{w}_ {1}, \cdots, \boldsymbol{w}_ {m}\right] .$  
step 3, 计算残差 $\boldsymbol{r}=\boldsymbol{b}-\boldsymbol{A} \boldsymbol{x} .$  
step 4, 解方程组 $\boldsymbol{W}^{\mathrm{T}} \boldsymbol{A} \boldsymbol{V} \boldsymbol{y}=\boldsymbol{W}^{\mathrm{T}} \boldsymbol{r}$ 得到 $\boldsymbol{y}$.  
step 5, 计算近似解 $\boldsymbol{x}:=\boldsymbol{x}+\boldsymbol{V} \boldsymbol{y}$.

### 4. 3阶矩阵的Householder变换和Givens变换

习题1. 设

$$
\boldsymbol{A}=\left[\begin{array}{ccc}
0 & 4 & 1 \\
1 & 1 & 1 \\
0 & 3 & 2
\end{array}\right]
$$

利用 Householder 变换求 $\boldsymbol{A}$ 的 $\mathrm{QR}$ 分解。

解：$\boldsymbol{s}_ 1 = [0,1,0]^ {\top}$,  $e_1 = [1, 0, 0]^ {\top}$,  $c_1 = -\mathrm{sgn}(a_ {11})\Vert s_1 \Vert = 1$, $u_1 = s_1 - c_1 e_1 = [-1,1,0]^ {\top}$,

$$
\boldsymbol{H}_1 = I - 2\frac{u_1 u_1^{\top}}{u_1^{\top}u_1} =
\left[\begin{array}{ccc}
0 & 1 & 0 \\
1 & 0 & 0 \\
0 & 0 & 1
\end{array}\right]
$$

则有

$$
A_1 = H_1A =\left[\begin{array}{ccc}
0 & 1 & 0 \\
1 & 0 & 0 \\
0 & 0 & 1
\end{array}\right]
\left[\begin{array}{ccc}
0 & 4 & 1 \\
1 & 1 & 1 \\
0 & 3 & 2
\end{array}\right]
= \left[\begin{array}{ccc}
1 & 1 & 1 \\
0 & 4 & 1 \\
0 & 3 & 2
\end{array}\right]
$$

$\boldsymbol{s}_ 2 = [0,4,3]^ {\top}$，($H_2A_1$ 不改变 $A_1$ 的第一行)，$e_2 = [0, 1, 0]^ {\top}$,  $c_2 = -\mathrm{sgn}(a_ {22})\Vert s_2 \Vert = -5$,  $u_2 = s_2-c_2e_2 = [0,9,3]^ {\top}$, 

$$
H_2 = I - 2\frac{u_2 u_2^{\top}}{u_2^{\top}u_2} =
\left[\begin{array}{ccc}
1 & 0 & 0 \\
0 & -\frac{4}{5} & -\frac{3}{5} \\
0 & -\frac{3}{5} & \frac{4}{5}
\end{array}\right]
$$

于是

$$
A_2 = H_2A_1 =\left[\begin{array}{ccc}
1 & 0 & 0 \\
0 & -\frac{4}{5} & -\frac{3}{5} \\
0 & -\frac{3}{5} & \frac{4}{5}
\end{array}\right]
\left[\begin{array}{ccc}
1 & 1 & 1 \\
0 & 4 & 1 \\
0 & 3 & 2
\end{array}\right]
= \left[\begin{array}{ccc}
1 & 1 & 1 \\
0 & -5 & -2 \\
0 & 0 & 1
\end{array}\right]
$$

从而

$$
R = A_2 
= \left[\begin{array}{ccc}
1 & 1 & 1 \\
0 & -5 & -2 \\
0 & 0 & 1
\end{array}\right]
$$

故$H_2H_1A=R$，令 $Q^{\top}=H_2H_1$，由于$H_1, H_2$ 为对称矩阵，故

$$
Q = (H_2H_1)^{\top} = H_1^{\top}H_2^{\top}
=H_1H_2 
= \left[\begin{array}{ccc}
0 & -\frac{4}{5} & -\frac{3}{5}\\
1 & 0 & 0 \\
0 & -\frac{3}{5} & \frac{4}{5}
\end{array}\right]
$$

满足 $A=QR$.


习题2. 设

$$
\boldsymbol{A}=\left[\begin{array}{ccc}
2 & 2 & 1 \\
0 & 2 & 2 \\
2 & 1 & 2
\end{array}\right]
$$

利用 Givens 变换求矩阵 $\boldsymbol{A}$ 的 QR 分解。

$s_1 = [2,0,2]^{\top}$，取 $G_ {12}=I_ {3}$ 即可，则有 $A_1=G_ {12}A=A$。

对于 $G_ {13}$，$c=\frac{x_1}{\sqrt{x_1^2+x_3^2}}=\frac{\sqrt{2}}{2}$，$c=\frac{x_3}{\sqrt{x_1^2+x_3^2}}=\frac{\sqrt{2}}{2}$，令

$$
G_ {13}=\left[\begin{array}{ccc}
\frac{\sqrt{2}}{2} & 0 &  \frac{\sqrt{2}}{2}\\
0 & 1 & 0 \\
-\frac{\sqrt{2}}{2} & 0 & \frac{\sqrt{2}}{2}
\end{array}\right]
$$

于是

$$
A_2 = G_ {13}A_1 =\left[\begin{array}{ccc}
2\sqrt{2} & \frac{3\sqrt{2}}{2} &  \frac{3\sqrt{2}}{2}\\
0 & 2 & 2 \\
0 & -\frac{\sqrt{2}}{2} & \frac{\sqrt{2}}{2}
\end{array}\right]
$$

$s_2 = [ \frac{3\sqrt{2}}{2},2, -\frac{\sqrt{2}}{2}]^{\top}$，对于 $G_ {23}$，$c=\frac{x_2}{\sqrt{x_2^2+x_3^2}}=\frac{2\sqrt{2}}{3}$，$s=\frac{x_3}{\sqrt{x_2^2+x_3^2}}=-\frac{1}{3}$，令

$$
G_ {23}=\left[\begin{array}{ccc}
1 & 0 & 0 \\
0 & \frac{2\sqrt{2}}{3} &-\frac{1}{3}\\
0 & \frac{1}{3} & \frac{2\sqrt{2}}{3}\\
\end{array}\right]
$$

于是

$$
R = A_3 = G_ {23}A_2 =\left[\begin{array}{ccc}
2\sqrt{2} & \frac{3\sqrt{2}}{2} &  \frac{3\sqrt{2}}{2}\\
0 & \frac{3\sqrt{2}}{2} & \frac{7\sqrt{2}}{6}\\
0 & 0 & \frac{4}{3}
\end{array}\right]
$$

故$G_ {23}G_ {13}G_ {12}A=R$，令 $Q^{\top}=G_ {23}G_ {13}G_ {12}$，于是

$$
Q = (G_ {23}G_ {13}G_ {12})^{\top} =\left[\begin{array}{ccc}
\frac{\sqrt{2}}{2} & 0 & \frac{\sqrt{2}}{2}  \\
\frac{\sqrt{2}}{6} & \frac{2\sqrt{2}}{3} & -\frac{\sqrt{2}}{6} \\
-\frac{2}{3}  & \frac{1}{3} & \frac{2}{3}\\
\end{array}\right]^{\top}=
\left[\begin{array}{ccc}
\frac{\sqrt{2}}{2} & \frac{\sqrt{2}}{6} & -\frac{2}{3}  \\
0 & \frac{2\sqrt{2}}{3} & \frac{1}{3} \\
\frac{\sqrt{2}}{2}  & -\frac{\sqrt{2}}{6} & \frac{2}{3}\\
\end{array}\right]
$$

满足 $A=QR$。

### 5. HSS迭代中参数$\alpha$与特征值分布

**命题1.** 设 $A \in \mathbb{C}^{n \times n}$ 为正定矩阵, $\boldsymbol{H}=\frac{1}{2}\left(\boldsymbol{A}+\boldsymbol{A}^{\mathrm{H}}\right)$ 和$S=\frac{1}{2}\left(\boldsymbol{A}-\boldsymbol{A}^{\mathrm{H}}\right)$ 为其 Hermite 和反 Hermite 部分, 且 $\alpha$ 是一个正常数. 则 $\mathrm{HSS}$ 迭代矩阵 $\boldsymbol{M}(\alpha)$ 为

$$
\boldsymbol{M}(\alpha)=(\alpha \boldsymbol{I}+\boldsymbol{S})^{-1}(\alpha \boldsymbol{I}-\boldsymbol{H})(\alpha \boldsymbol{I}+\boldsymbol{H})^{-1}(\alpha \boldsymbol{I}-\boldsymbol{S})
$$

且其谱半径 $\rho(\boldsymbol{M}(\alpha))$ 有上界

$$
\sigma(\alpha) \equiv \max _ {\lambda_ {i} \in \lambda(\boldsymbol{H})}\left|\frac{\alpha-\lambda_ {i}}{\alpha+\lambda_ {i}}\right|
$$

式中: $\lambda(\boldsymbol{H})$ 为矩阵 $\boldsymbol{H}$ 的谱集，因此成立

$$
\rho(\boldsymbol{M}(\alpha)) \leqslant \sigma(\alpha)<1, \quad \forall \alpha>0
$$

即 $\mathrm{HSS}$ 迭代收敘到线性方程组 $\boldsymbol{Ax}=\boldsymbol{b}$ 的唯一解 $\boldsymbol{x}^{*} \in \mathbb{C}^{n} .$

**证明**

记

$$
\boldsymbol{M}_ {1}=\alpha \boldsymbol{I}+\boldsymbol{H}, \boldsymbol{N}_ {1}=\alpha \boldsymbol{I}-\boldsymbol{S}, \boldsymbol{M}_ {2}=\alpha \boldsymbol{I}+\boldsymbol{S}, \boldsymbol{N}_ {2}=\alpha \boldsymbol{I}-\boldsymbol{H}
$$

注意到对任意的 $\alpha>0$，$\alpha \boldsymbol{I}+\boldsymbol{H}$ 和 $\alpha \boldsymbol{I}+\boldsymbol{S}$ 是非奇异的，利用矩阵谱的相似不变性，得

$$
\begin{aligned}
\rho(\boldsymbol{M}(\alpha)) &=\rho\left((\alpha \boldsymbol{I}+\boldsymbol{S})^{-1}(\alpha \boldsymbol{I}-\boldsymbol{H})(\alpha \boldsymbol{I}+\boldsymbol{H})^{-1}(\alpha \boldsymbol{I}-\boldsymbol{S})\right) \\
&=\rho\left((\alpha \boldsymbol{I}-\boldsymbol{H})(\alpha \boldsymbol{I}+\boldsymbol{H})^{-1}(\alpha \boldsymbol{I}-\boldsymbol{S})(\alpha \boldsymbol{I}+\boldsymbol{S})^{-1}\right) \\
& \leqslant\left\|(\alpha \boldsymbol{I}-\boldsymbol{H})(\alpha \boldsymbol{I}+\boldsymbol{H})^{-1}\right\|_ {2} \cdot\left\|(\alpha \boldsymbol{I}-\boldsymbol{S})(\alpha \boldsymbol{I}+\boldsymbol{S})^{-1}\right\|_ {2}
\end{aligned}
$$

记 $\boldsymbol{Q}(\alpha)=(\alpha \boldsymbol{I}-\boldsymbol{S})(\alpha \boldsymbol{I}+\boldsymbol{S})^{-1}$， 注意到 $\boldsymbol{S}^{\mathrm{H}}=-\boldsymbol{S}$，有

$$
\begin{aligned}
\boldsymbol{Q}(\alpha)^{\mathrm{H}} \boldsymbol{Q}(\alpha) &=(\alpha \boldsymbol{I}-\boldsymbol{S})^{-1}(\alpha \boldsymbol{I}+\boldsymbol{S})(\alpha \boldsymbol{I}-\boldsymbol{S})(\alpha \boldsymbol{I}+\boldsymbol{S})^{-1} \\
&=(\alpha \boldsymbol{I}-\boldsymbol{S})^{-1}(\alpha \boldsymbol{I}-\boldsymbol{S})(\alpha \boldsymbol{I}+\boldsymbol{S})(\alpha \boldsymbol{I}+\boldsymbol{S})^{-1}=\boldsymbol{I}
\end{aligned}
$$

即 $\boldsymbol{Q}(\alpha)$ 是酉矩阵，故有 $\|Q(\alpha)\|_ {2}=1$。从而有

$$
\rho(\boldsymbol{M}(\alpha)) \leqslant\left\|(\alpha \boldsymbol{I}-\boldsymbol{H})(\alpha \boldsymbol{I}+\boldsymbol{H})^{-1}\right\|_ {2}=\max _ {\lambda_ {i} \in \lambda(\boldsymbol{H})}\left|\frac{\alpha-\lambda_ {i}}{\alpha+\lambda_ {i}}\right| \equiv \sigma(\alpha)
$$

由于 $\lambda_ {i}>0(i=1,2, \cdots, n)$ 及 $\alpha>0$，容易推得 $\rho(\boldsymbol{M}(\alpha)) \leqslant$$\sigma(\alpha)<1$ 。



**命题2.** 设 $A \in \mathbb{C}^{n \times n}$ 为正定矩阵，$\boldsymbol{H}=\frac{1}{2}\left(\boldsymbol{A}+\boldsymbol{A}^{\mathrm{H}}\right)$ 和$\boldsymbol{S}=\frac{1}{2}\left(\boldsymbol{A}-\boldsymbol{A}^{\mathrm{H}}\right)$ 为其 Hermite 和反Hermite 部分，且 $\lambda_ {\min }$ 和 $\lambda_ {\max }$分别为矩阵 $\boldsymbol{H}$ 的最小和最大特征值，$\alpha$ 是一个正常数。则

$$
\alpha^{*} \equiv \arg \min _ {\alpha}\left\{\max _ {\lambda_ {\min } \leqslant \lambda \leqslant \lambda_ {\max }}\left|\frac{\alpha-\lambda}{\alpha+\lambda}\right|\right\}=\sqrt{\lambda_ {\min } \lambda_ {\max }}
$$

且

$$
\sigma\left(\alpha^{*}\right)=\frac{\sqrt{\lambda_ {\max }}-\sqrt{\lambda_ {\min }}}{\sqrt{\lambda_ {\max }}+\sqrt{\lambda_ {\min }}}=\frac{\sqrt{\kappa(\boldsymbol{H})}-1}{\sqrt{\kappa(\boldsymbol{H})}+1},
$$

式中: $\kappa(\boldsymbol{H})=\lambda_ {\max } / \lambda_ {\min }$ 为 $\boldsymbol{H}$ 的谱条件数。

**证明**

对任意的 $\alpha>0$, 函数 $f(\lambda)=\frac{\alpha-\lambda}{\alpha+\lambda}$ 关于 $\lambda$ 是单调递减的 $\left(f^{\prime}(\lambda)=-2 \alpha /(\alpha+\lambda)^2<0\right)$, 故有

$$
\sigma(\alpha)=\max \left\{\left|\frac{\alpha-\lambda_ {\min }}{\alpha+\lambda_ {\min }}\right|,\left|\frac{\alpha-\lambda_ {\max }}{\alpha+\lambda_ {\max }}\right|\right\}.
$$

若 $\alpha^ {*}$ 是 $\sigma(\alpha)$ 的极小点，则必有 $\alpha^ {*}=\lambda_ {\min }>0$， $\alpha^ {*}-\lambda_ {\max }<0$,

$$
\frac{\alpha^ {*}-\lambda_ {\min }}{\alpha^ {*} + \lambda_ {\min }}=-\frac{\alpha^{*}-\lambda_ {\max }}{\alpha {*}+\lambda_ {\max }}
$$

从上式解得

$$
\alpha^ {*}=\sqrt{\lambda_ {\min } \lambda_ {\max }}
$$

从而推论的结论成立。

### 6. 迭代法的外推思想

使用迭代格式求解方程组 $\boldsymbol{Ax}=\boldsymbol{b}$ 时，迭代过程可能收敛，也可能不收敛。我们希望找到一种改进方法，使得不收敛的格式变得收敛，收敛慢的格式变得收敛快。

给定第 $k$ 次迭代后的估计解 $\boldsymbol{x}^{(k)}$，由迭代格式 $\boldsymbol{x}^{k+1}=B\boldsymbol{x}^{k}+\boldsymbol{f},\ k=1,2,\cdots$ 得到了一 个新的估计解 $\boldsymbol{B} \boldsymbol{x}^{(k)}+\boldsymbol{f}$，再将它与 $\boldsymbol{x}^{(k)}$，$\boldsymbol{x}^{(k-1)}$ 等进行适当的组合，就可能得到更快的迭代收敛速度。这就是迭代法加速的基本思想。

在迭代格式 $\boldsymbol{x}^{k+1}=B\boldsymbol{x}^{k}+\boldsymbol{f},\ k=1,2,\cdots$ 中引入一个参数 $\gamma \neq 0$，构造新的迭代格式

$$
\boldsymbol{x}^{(k+1)}=(1-\gamma) \boldsymbol{x}^{(k)}+\gamma\left(\boldsymbol{B} \boldsymbol{x}^{(k)}+\boldsymbol{f}\right):=\boldsymbol{B}_ {\gamma} \boldsymbol{x}^{(k)}+\gamma \boldsymbol{f}, \quad k=0,1, \cdots \label{eq1} \tag{1}
$$

式中:

$$
\boldsymbol{B}_ {\gamma}=(1-\gamma) \boldsymbol{I}+\gamma \boldsymbol{B}
$$

若上述迭代格式收敛到某个 $x^{*}$，则有

$$
\boldsymbol{x}^{*}=(1-\gamma) \boldsymbol{x}^{*}+\gamma\left(\boldsymbol{B} \boldsymbol{x}^{*}+\boldsymbol{f}\right) \Longrightarrow \boldsymbol{x}^{*}=\boldsymbol{B} \boldsymbol{x}^{*}+\boldsymbol{f}
$$

所以，不管 $\gamma \neq 0$ 如何选取，当迭代格式 $\eqref{eq1}$ 收敛时，它必定收敛到原来方程组 $\boldsymbol{x}=\boldsymbol{B} \boldsymbol{x}+\boldsymbol{f}$ 的解。因此希望选取一个比较好的参数 $\gamma$，使得迭代格式 $\eqref{eq1}$ 收敛得尽可能快。这就是外推法。

### 7. 广义极小参量法（GMRES）

广义极小残量法 $(\mathrm{GMRES})$ ：
(1) 令 $\boldsymbol{v}_ {1}=\boldsymbol{r}_ {0} /\left\|\boldsymbol{r}_ {0}\right\|_ {2}$，产生一个长度为 $k$ 的 Arnoldi 分解 $\boldsymbol{A} \boldsymbol{V}_ {k}=\boldsymbol{V}_ {k} \boldsymbol{H}_ {k}+\beta_ {k} \boldsymbol{v}_ {k+1} \boldsymbol{e}_ {k}^{\mathrm{T}}=\boldsymbol{V}_ {k+1} \widetilde{\boldsymbol{H}}_ {k+1, k}$  
(2) 利用 Givens 变换求 $\widetilde{\boldsymbol{H}}_ {k+1, k}$ 的 QR 分解 $\left(\boldsymbol{G}_ {k} \boldsymbol{G}_ {k-1} \cdots \boldsymbol{G}_ {2} \boldsymbol{G}_ {1}\right) \widetilde{\boldsymbol{H}}_ {k+1, k}=\left[\begin{array}{c}\boldsymbol{R}_ {k} \\ \mathbf{0}\end{array}\right]$ 并按式 

$$\left \{ \begin{array}{l} \tau_ {1}=\beta c_ {1} \\ \tau_ {i}=(-1)^ {i-1} \beta s_ {1} s_ {2} \cdots s_ {i-1} c_ {i}, \quad i=2,3, \cdots, k \\ \rho_ {k}=(-1)^{k} \beta s_ {1} s_ {2} \cdots s_ {k} \end{array} \right.
$$
求得向量 $\boldsymbol{t}_ {k}$ 和数 $\rho_ {k} .$  
(3) 用回代法求解上三角方程组 $\boldsymbol{R}_ {k} \boldsymbol{z}_ {k}=\boldsymbol{t}_ {k}$, 得 $\boldsymbol{z}_ {k} .$  
(4) 计算 $\boldsymbol{x}_ {k}=\boldsymbol{x}_ {0}+\boldsymbol{V}_ {k} \boldsymbol{z}_ {k} .$  
(5) 若 $\left|\rho_ {k}\right| / \beta<\varepsilon($ 事先给定的误差界 $)$，则终止；否则增加 $k$ 的值，重复上面的过程。

> 重启技术是为了节省计算机内存，预处理技术是为了使系数矩阵的谱相对集中，从而加快速度。

### 8. BCGSTAB方法避免残差震荡，保持下降的修正策略

$\mathrm{BCGSTAB}$ 方法是为了改进 $\mathrm{CGS}$ 方法之残量的范数剧烈抖动而提出的。这一方法的基本思想是 $\mathrm{CGS}$ 方法的残量 $\boldsymbol{r}_ {k}^{\mathrm{CGS}}$ 满足

$$
\boldsymbol{r}_ {k}^{\mathrm{CGS}}=\varphi_ {k}(\boldsymbol{A}) \boldsymbol{r}_ {k}^{\mathrm{BCG}}=\varphi_ {k}^{2}(\boldsymbol{A}) \boldsymbol{r}_ {0}
$$

可以考虑不用多项式 $\varphi_ {k}$，而是选择一个其他的 $k$ 次多项式 $\widetilde{\varphi}_ {k}$，使

$$
\boldsymbol{r}_ {k}=\widetilde{\varphi}_ {k}(\boldsymbol{A}) \boldsymbol{r}_ {k}^{\mathrm{BCG}}=\widetilde{\varphi}_ {k}(\boldsymbol{A}) \varphi_ {k}(\boldsymbol{A}) \boldsymbol{r}_ {0}
$$

以期待这样选取的相对残差范数 $\left\|\boldsymbol{r}_ {k}\right\|_ {2} /\left\|\boldsymbol{r}_ {0}\right\|_ {2}$ 的振荡性有所改进。

一种选择 $\widetilde{\varphi}_ {k}$ 的方法是根据递推公式

$$
\widetilde{\varphi}_ {0}(t)=1, \quad \widetilde{\varphi}_ {k+1}(t)=\left(1-\omega_ {k+1} t\right) \widetilde{\varphi}_ {k}(t)
$$

式中: $\omega_ {k+1}$ 为待定参数。BCGSTAB 方法正是利用这一参数的可选择性来改进相对残差范数的振荡性。

### 9. 线性最小二乘问题的等价形式

设 $A \in \mathbb{C}^{m \times n}, \boldsymbol{b} \in \mathbb{C}^{m}$, 确定 $\boldsymbol{x} \in \mathbb{C}^{n}$ 使得

$$
\|\boldsymbol{A} \boldsymbol{x}-\boldsymbol{b}\|_ {2}=\min _ {\boldsymbol{z} \in \mathbb{C}^{n}}\|\boldsymbol{A} \boldsymbol{z}-\boldsymbol{b}\|_ {2} \label{eq2}
\tag{2}
$$

问题 $\eqref{eq2}$ 称为线性最小二乘问题 $($ Least Squares, $\mathrm{LS}$ 问题 $)$，而 $\boldsymbol{x}$ 则称为最小二乘解或极小解。称 $\boldsymbol{r}(\boldsymbol{x})=\boldsymbol{b}-\boldsymbol{A} \boldsymbol{x}$ 为残差向量 ，简称残量。

> $Ax=b$ 有解 $\iff$ $AA^{+}b=b$，通解 $x=A^{+}b + (I-A^{+}A)z$，唯一极小范数解 $x_ {0}=A^{+}b$，
>
> $Ax=b$ 无解 $\Longrightarrow$ 唯一极小范数最小二乘解 $x_ {0}=A^{+}b$。

***法方程***

$\boldsymbol{x}$ 是最小二乘问题 $\eqref{eq2}$ 的极小解，即 $x \in S_ {\mathrm{LS}}$ 的充分必要条件是 $\boldsymbol{x}$ 为方程

$$
\boldsymbol{A}^{\mathrm{H}} \boldsymbol{A} \boldsymbol{x}=\boldsymbol{A}^{\mathrm{H}} \boldsymbol{b} \label{eq3} \tag{3}
$$

的解，其中式 $\eqref{eq3}$ 称为最小二乘问题的法方程。  
**证明**：最小二乘问题 $\eqref{eq2}$ 等价于极小化问题

$$
\min \varphi(\boldsymbol{x})=\frac{1}{2}\|\boldsymbol{A} \boldsymbol{x}-\boldsymbol{b}\|_ {2}^{2}=\frac{1}{2} \boldsymbol{x}^{\mathrm{H}}\left(\boldsymbol{A}^{\mathrm{H}} \boldsymbol{A}\right) \boldsymbol{x}-\left(\boldsymbol{b}^{\mathrm{H}} \boldsymbol{A}\right) \boldsymbol{x}+\frac{1}{2}\|\boldsymbol{b}\|_ {2}^{2}
$$

由于 $\boldsymbol{A}^{\mathrm{H}} \boldsymbol{A} \in \mathbb{C}^{n \times n}$ 是半正定矩阵，因此 $n$ 元实函数 $\varphi(\boldsymbol{x})$ 是凸函数，故 $x$ 是最小二乘问题 $\eqref{eq2}$ 的极小解等价于

$$
\nabla \varphi(\boldsymbol{x})=\boldsymbol{A}^{\mathrm{H}}(\boldsymbol{A} \boldsymbol{x}-\boldsymbol{b})=\mathbf{0}
$$

***KKT方程***

设 $A \in \mathbb{C}^{m \times n}$，$\boldsymbol{b} \in \mathbb{R}^{n}$。则 $\boldsymbol{x}$ 和 $\boldsymbol{r}=\boldsymbol{b}-\boldsymbol{A} \boldsymbol{x}$ 分别为最小二乘问题 $\eqref{eq2}$ 的极小解和残量的充分必要条件是 $x$ 和$\boldsymbol{r}$
为鞍点系统

$$
\left[\begin{array}{ll}
\boldsymbol{I} & \boldsymbol{A} \\
\boldsymbol{A}^{\mathrm{H}} & \boldsymbol{O}
\end{array}\right]\left[\begin{array}{l}
\boldsymbol{r} \\
\boldsymbol{x}
\end{array}\right]=\left[\begin{array}{l}
\boldsymbol{b} \\
\mathbf{0}
\end{array}\right] \label{eq4} \tag{4}
$$

的解。上述线性系统称为最小二乘问题的KKT 方程。  
**证明** ：若 $x$ 为最小二乘问题 $\eqref{eq2}$ 的极小解，而 $\boldsymbol{r}=\boldsymbol{b}-\boldsymbol{A} \boldsymbol{x}$ 为 其残量，则

$$
x=A^{\dagger} b+\left(I-A^{\dagger} A\right) z, \quad r=b-A x=b-A A^{\dagger} b=\left(I-A A^{\dagger}\right) b
$$

由广义逆 $A^{\dagger}$ 的性质及等式 $r+A x=b$，得

$$
\boldsymbol{A}^{\mathrm{H}} \boldsymbol{r}=\boldsymbol{A}^{\mathrm{H}}\left(\boldsymbol{I}-\boldsymbol{A} \boldsymbol{A}^{\dagger}\right) \boldsymbol{b}=\left[\left(\boldsymbol{I}-\boldsymbol{A} \boldsymbol{A}^{\dagger}\right) \boldsymbol{A}\right]^{\mathrm{H}} \boldsymbol{b}=\mathbf{0} .
$$

故式 $\eqref{eq4}$ 是相容的线性系统，且 $x$ 和 $r$ 满足式 $\eqref{eq4}$。
反之，通过验证广义逆的四个条件，可以验证

$$
\boldsymbol{B}^{\dagger} \equiv\left[\begin{array}{cc}
\boldsymbol{I} & \boldsymbol{A} \\
A^{\mathrm{H}} & \boldsymbol{O}
\end{array}\right]^{\dagger}=\left[\begin{array}{cc}
\boldsymbol{I}-\boldsymbol{A} \boldsymbol{A}^{\dagger} & \left(\boldsymbol{A}^{\dagger}\right)^{\mathrm{H}} \\
\boldsymbol{A}^{\dagger} & -\boldsymbol{A}^{\dagger}\left(\boldsymbol{A}^{\dagger}\right)^{\mathrm{H}}
\end{array}\right]
$$

故式 $\eqref{eq4}$ 的任一解向量 $x, r$ 有如下形式

$$
\left[\begin{array}{l}r \\ x\end{array}\right]=B^{\dagger}\left[\begin{array}{l}b \\ 0\end{array}\right]+\left(I-B^{\dagger} B\right)\left[\begin{array}{l}y \\ z\end{array}\right]=\left[\begin{array}{c}\left(I-A A^{\dagger}\right) b \\ A^{\dagger} b-\left(I-A^{\dagger} A\right) z\end{array}\right]
$$

式中: $\boldsymbol{y} \in \mathbb{C}^{m}, \boldsymbol{z} \in \mathbb{C}^{n}$ 为任意向量。故满足式 $\eqref{eq4}$ 的任一组向量 $\boldsymbol{x}, \boldsymbol{r}$ 分别为式 $\eqref{eq4}$ 的极小解和残量。证毕。

### 10. 求解LS问题的QR方法

本节考虑最小二乘问题 $\eqref{eq2}$ 中的矩阵 $A \in \mathbb{R}^{m \times n}(m \geqslant n)$ 的情形。根据正交矩阵保持向量 2 范数不变的性质，对于任意的正交矩阵 $Q \in \mathbb{R}^{m \times m}$，最小二乘问题 $\eqref{eq2}$ 等价于

$$
\left\|\boldsymbol{Q}^{\mathrm{T}}(\boldsymbol{A} \boldsymbol{x}-\boldsymbol{b})\right\|_ {2}=\min \left\{\left\|\boldsymbol{Q}^{\mathrm{T}}(\boldsymbol{A} \boldsymbol{z}-\boldsymbol{b})\right\|_ {2}: \boldsymbol{z} \in \mathbb{R}^{n}\right\}
\label{eq5} \tag{5}
$$

这样, 就可望通过适当选取正交矩阵 $Q$，使原问题 $\eqref{eq2}$ 转化为较为容易求解的最小二乘问题 $\eqref{eq5}$ ，这就是正交化方法一 $\mathrm{QR}$ 分解方法的基本思想。

### 11. 求矩阵特征值的基本QR方法下三角矩阵趋于0

**定理1** 设对称矩阵 $A \in \mathbb{R}^{n \times n}$ 满足

$$
\left|\lambda_ {1}\right|>\left|\lambda_ {2}\right| \geqslant \cdots \geqslant\left|\lambda_ {n}\right|>0
$$

对应的规范正交化特征向量为 $\boldsymbol{x}_ {1}, \boldsymbol{x}_ {2}, \cdots, \boldsymbol{x}_ {n}$。如果单位坐标向量

$$
\boldsymbol{e}_ {1}=(1,0, \cdots, 0)^{\mathrm{T}}=\sum_ {i=1}^{n} \alpha_ {i} \boldsymbol{x}_ {i}
$$

中的 $\alpha_ {1} \neq 0$，那么由算法 $7.5$ 生成的矩阵序列 ${\boldsymbol{A}_ {k}}$ 具有收敛性质

$$
\lim _ {k \rightarrow \infty} \boldsymbol{A}_ {k} \boldsymbol{e}_ {1}=\lambda_ {1} \boldsymbol{e}_ {1}
$$

**证明**：记 $\widetilde{Q}_ {k}$ 的第 1 列为 $\widetilde{\boldsymbol{q}}_ {1}^{(k)}=\widetilde{\boldsymbol{Q}}_ {k} \boldsymbol{e}_ {1}$，$\widetilde{\boldsymbol{R}}_ {k}$ 的第 1 个对角元为$\widetilde{r}_ {11}^{(k)}$，则由式 $(7.52)$ 可知

$$
\boldsymbol{A}^{k} \boldsymbol{e}_ {1}=\widetilde{\boldsymbol{Q}}_ {k} \widetilde{\boldsymbol{R}}_ {k} \boldsymbol{e}_ {1}=\widetilde{\boldsymbol{Q}}_ {k}\left(\widetilde{r}_ {11}^{(k)} \boldsymbol{e}_ {1}\right)=\widetilde{r}_ {11}^{(k)} \widetilde{\boldsymbol{q}}_ {1}^{(k)}
$$

注意到 $\left\| \widetilde{\boldsymbol{q}}_ {1}^ {(k)}\right\|_ {2}=1, \widetilde{r}_ {11}^ {(k)} \neq 0$ (因矩阵 $A$ 非奇异 )。从而 $\left|\widetilde{r}_ {11}^ {(k)}\right|=\left\| A^ {k} e_ {1} \right\|_ {2}$。于是，根据幂法的收敛性，有

$$
\lim _ {k \rightarrow \infty} \widetilde{\boldsymbol{q}}_ {1}^{(k)}=\lim _ {k \rightarrow \infty} \frac{\boldsymbol{A}^{k} \boldsymbol{e}_ {1}}{\widetilde{r}_ {11}^{(k)}}=\boldsymbol{z}_ {1}
$$

式中：$\boldsymbol{z}_ {1}$ 为矩阵 $\boldsymbol{A}$ 的对应于 $\lambda_ {1}$ 的规范化特征向量 $($ 可以相差一个常数因子 $\pm 1) .$ 进而, 由式 $(7.51)$，可以把 $\boldsymbol{A}_ {k+1}$ 写成

$$
\begin{aligned}
\boldsymbol{A}_ {k+1} &=\widetilde{\boldsymbol{Q}}_ {k}^{\mathrm{T}} \boldsymbol{A} \widetilde{\boldsymbol{Q}}_ {k}=\left[\widetilde{\boldsymbol{q}}_ {1}^{(k)}, \widetilde{\boldsymbol{Q}}_ {k-1}\right]^{\mathrm{T}} \boldsymbol{A}\left[\widetilde{\boldsymbol{q}}_ {1}^{(k)}, \widetilde{\boldsymbol{Q}}_ {k-1}\right] \\
&=\left[\begin{array}{cc}
a_ {11}^{(k+1)} & * \\
\boldsymbol{\alpha}^{(k+1)} & *
\end{array}\right]
\end{aligned}
$$

式中：$a_ {11}^{(k+1)}=\left(\widetilde{\boldsymbol{q}}_ {1}^{(k)}\right)^{\mathrm{T}} \boldsymbol{A} \widetilde{\boldsymbol{q}}_ {1}^{(k)}, \boldsymbol{\alpha}^{(k+1)}=\left(\widetilde{\boldsymbol{Q}}_ {k-1}^{(k)}\right)^{\mathrm{T}} \boldsymbol{A} \widetilde{\boldsymbol{q}}_ {1}^{(k)}$。 根据式 $(7.55)$，得

$$
\lim _ {k \rightarrow \infty} a_ {11}^{(k+1)}=\lim _ {k \rightarrow \infty}\left(\widetilde{\boldsymbol{q}}_ {1}^{(k)}\right)^{\mathrm{T}} \boldsymbol{A} \widetilde{\boldsymbol{q}}_ {1}^{(k)}=\boldsymbol{z}_ {1}^{\mathrm{T}} \boldsymbol{A} \boldsymbol{z}_ {1}=\lambda_ {1}
$$

故由引理 7.1，得

$$
\begin{aligned}
\lim _ {k \rightarrow \infty}\left\|\boldsymbol{\alpha}^{(k+1)}\right\|_ {2} &=\lim _ {k \rightarrow \infty}\left\|\boldsymbol{A} \widetilde{\boldsymbol{q}}_ {1}^{(k)}-\left[\left(\widetilde{\boldsymbol{q}}_ {1}^{(k)}\right)^{\mathrm{T}} \boldsymbol{A} \widetilde{\boldsymbol{q}}_ {1}^{(k)}\right] \widetilde{\boldsymbol{q}}_ {1}^{(k)}\right\|_ {2} \\
&=\left\|\boldsymbol{A} \boldsymbol{z}_ {1}-\left[\boldsymbol{z}_ {1}^{\mathrm{T}} \boldsymbol{A} \boldsymbol{z}_ {1}\right] \boldsymbol{z}_ {1}\right\|_ {2}=0 .
\end{aligned}
$$

因此, $\lim _ {k \rightarrow \infty} \boldsymbol{\alpha}^{(k+1)}=\mathbf{0}$。于是有

$$
\lim _ {k \rightarrow \infty} \boldsymbol{A}_ {k+1}=\left[\begin{array}{cc}
\lambda_ {1} & * \\
\mathbf{0} & *
\end{array}\right],
$$

即式 $(7.54)$ 成立。证毕。

**定理2** 设 $\boldsymbol{A}=\boldsymbol{X} \boldsymbol{\Lambda} \boldsymbol{X}^{-1}$，其中 $\boldsymbol{\Lambda}=\operatorname{diag}\left(\lambda_ {1}, \lambda_ {2}, \cdots, \lambda_ {n}\right)$。如果 (1) $\left|\lambda_ {1}\right|>\left|\lambda_ {2}\right|>\cdots>\left|\lambda_ {n}\right|>0$；(2) $\boldsymbol{X}^{-1}$ 具有 LU 分解：$\boldsymbol{X}^{-1}=\boldsymbol{L} \boldsymbol{U}$， 则 $\mathrm{QR}$ 方法基本收敛于上三角矩阵。
**证明**： 由式 $(7.51)$ 可知，只需分析 $\widetilde{Q}_ {k}$ 的极限情况。而由式 $(7.52)$，只需分析 $A^{k}$ 的极限情况。注意到

$$
\boldsymbol{A}^{k}=\boldsymbol{X} \boldsymbol{\Lambda}^{k} \boldsymbol{X}^{-1}=\boldsymbol{X} \boldsymbol{\Lambda}^{k} \boldsymbol{L} \boldsymbol{U}=\boldsymbol{X}\left(\boldsymbol{\Lambda}^{k} \boldsymbol{L} \boldsymbol{\Lambda}^{-k}\right) \boldsymbol{\Lambda}^{k} \boldsymbol{U}
$$

令

$$
\boldsymbol{\Lambda}^{k} \boldsymbol{L} \boldsymbol{\Lambda}^{-k}=\boldsymbol{I}+\boldsymbol{E}_ {k},
$$

则

$$
\boldsymbol{A}^{k}=\boldsymbol{X}\left(\boldsymbol{I}+\boldsymbol{E}_ {k}\right) \boldsymbol{\Lambda}^{k} \boldsymbol{U}
$$

由于 $L$ 是单位下三角形,，故

$$
\left(\boldsymbol{E}_ {k}\right)_ {i j}=\left\{\begin{array}{ll}
0, & i \leqslant j \\
l_ {i j}\left(\lambda_ {i} / \lambda_ {j}\right)^{k}, & i>j
\end{array}\right.
$$

由假设条件 $\mathbb{1}$ 知，$\boldsymbol{E}_ {k} \rightarrow \boldsymbol{O}$ 且 $\left(\boldsymbol{E}_ {k}\right)_ {i j}$ 的收敛速度是 $\left|\lambda_ {i} / \lambda_ {j}\right|$.
设 $\boldsymbol{X}=Q R$ 且 $R$ 的对角元均为正数，则有

$$
\begin{aligned}
\boldsymbol{A}^{k} &=\boldsymbol{Q} \boldsymbol{R}\left(\boldsymbol{I}+\boldsymbol{E}_ {k}\right) \boldsymbol{\Lambda}^{k} \boldsymbol{U} \\
&=\boldsymbol{Q}\left(\boldsymbol{I}+\boldsymbol{R} \boldsymbol{E}_ {k} \boldsymbol{R}^{-1}\right) \boldsymbol{R} \boldsymbol{\Lambda}^{k} \boldsymbol{U}
\end{aligned}
$$

因为 $\boldsymbol{E}_ {k} \rightarrow \boldsymbol{O}(k \rightarrow \infty)$, 故当 $k$ 充分大时， $\boldsymbol{I}+\boldsymbol{R} \boldsymbol{E}_ {k} \boldsymbol{R}^{-1}$ 非奇异，所以有唯一的 QR 分解$\boldsymbol{I}+\boldsymbol{R} \boldsymbol{E}_ {k} \boldsymbol{R}^{-1}=\widehat{\boldsymbol{Q}}_ {k} \widehat{\boldsymbol{R}}_ {k}\left(\widehat{\boldsymbol{R}}_ {k}\right.$ 的对角元为正 $)$，而且当 $k \rightarrow \infty$ 时，$\widehat{\boldsymbol{Q}}_ {k} \rightarrow \boldsymbol{I}, \widehat{\boldsymbol{R}}_ {k} \rightarrow \boldsymbol{I}$。 此时，$\boldsymbol{A}^{k}$ 有如下分解

$$
\boldsymbol{A}^{k}=\left(\boldsymbol{Q} \widehat{\boldsymbol{Q}}_ {k}\right)\left(\widehat{\boldsymbol{R}}_ {k} \boldsymbol{R} \boldsymbol{\Lambda}^{k} \boldsymbol{U}\right)
$$

妨碍上式成为 $A^{k}$ 的 $\mathrm{QR}$ 分解的仅仅是上式右端第 2 个因子(上三角矩阵 )​ 的对角元可能非正。为补救这一点，可引入两个对角正交矩阵

$$
\boldsymbol{D}_ {1}=\operatorname{diag}\left(\frac{\lambda_ {1}}{\left|\lambda_ {1}\right|}, \cdots, \frac{\lambda_ {n}}{\left|\lambda_ {n}\right|}\right), \quad \boldsymbol{D}_ {2}=\operatorname{diag}\left(\frac{\boldsymbol{U}_ {11}}{\left|\boldsymbol{U}_ {11}\right|}, \cdots, \frac{\boldsymbol{U}_ {n n}}{\left|\boldsymbol{U}_ {n n}\right|}\right)
$$

式中: $\boldsymbol{U}_ {i i}(i=1,2, \cdots, n)$ 为矩阵 $\boldsymbol{U}$ 的对角元。于是

$$
\boldsymbol{A}^{k}=\left(\left(\boldsymbol{Q} \widehat{\boldsymbol{Q}}_ {k}\right)\left(\boldsymbol{D}_ {2} \boldsymbol{D}_ {1}^{k}\right)\right)\left(\boldsymbol{D}_ {1}^{-k} \boldsymbol{D}_ {2}^{-1} \widehat{\boldsymbol{R}}_ {k} \boldsymbol{R} \boldsymbol{\Lambda}^{k} \boldsymbol{U}\right)
$$

是 $A^{k}$ 的唯一 QR 分解。从而由式 $(7.51)$，有

$$
\boldsymbol{A}_ {k+1}=\left(\boldsymbol{Q} \widehat{\boldsymbol{Q}}_ {k} \boldsymbol{D}_ {2} \boldsymbol{D}_ {1}^{k}\right)^{\mathrm{T}} \boldsymbol{A}\left(\boldsymbol{Q} \widehat{\boldsymbol{Q}}_ {k} \boldsymbol{D}_ {2} \boldsymbol{D}_ {1}^{k}\right)
$$

将 $\boldsymbol{A}=\boldsymbol{X} \boldsymbol{\Lambda} \boldsymbol{X}^{-1}=Q \boldsymbol{R} \boldsymbol{\Lambda} \boldsymbol{R}^{-1} \boldsymbol{Q}^{-1}$ 代入上式，得

$$
\boldsymbol{A}_ {k+1}=\left(\boldsymbol{D}_ {2} \boldsymbol{D}_ {1}^{k}\right)^{\mathrm{T}}\left(\widehat{\boldsymbol{Q}}_ {k}^{-1} \boldsymbol{R} \boldsymbol{\Lambda} \boldsymbol{R}^{-1} \widehat{\boldsymbol{Q}}_ {k}\right)\left(\boldsymbol{D}_ {2} \boldsymbol{D}_ {1}^{k}\right)
$$

因为 $\widehat{\boldsymbol{Q}}_ {k}^{-1} \boldsymbol{R} \boldsymbol{\Lambda} \boldsymbol{R}^{-1} \widehat{\boldsymbol{Q}}_ {k} \rightarrow \boldsymbol{R} \boldsymbol{\Lambda} \boldsymbol{R}^{-1} \equiv \overline{\boldsymbol{R}}($ 上三角矩阵 $)$，所以 $\boldsymbol{A}_ {k}$ 的对角线以下元素收敛于 $0$。因为 $D_ {1}^{k}$ 可能不收敛，故 $\boldsymbol{A}_ {k}$ 基本收敛 于 $\overline{\boldsymbol{R}}$。证毕。

### 12. Krylov子空间求大规模特征值问题的思想

设 $(\lambda, \boldsymbol{x})$ 是矩阵 $\boldsymbol{A}$ 的一个特征对，即

$$
\boldsymbol{A} \boldsymbol{x}=\lambda \boldsymbol{x} \label{eq6} \tag{6}
$$

再假定

$$
\boldsymbol{x}=\boldsymbol{V}_ {k} \boldsymbol{y}+\boldsymbol{v} \approx \boldsymbol{V}_ {k} \boldsymbol{y}
$$

于是，由式 $\eqref{eq6}$有

$$
\boldsymbol{A} \boldsymbol{V}_ {k} \boldsymbol{y} \approx \lambda \boldsymbol{V}_ {k} \boldsymbol{y}
$$

从而有

$$
\boldsymbol{H}_ {k} \boldsymbol{y}=\boldsymbol{V}_ {k}^{\mathrm{T}} \boldsymbol{A} \boldsymbol{V}_ {k} \boldsymbol{y} \approx \lambda \boldsymbol{y}
$$

因此，若 $(\lambda, \boldsymbol{y})$ 是 $\boldsymbol{H}_ {k}=\boldsymbol{V}_ {k}^{\mathrm{T}} \boldsymbol{A} \boldsymbol{V}_ {k}$ 的特征对，则可用 $\left(\lambda, \boldsymbol{V}_ {k} \boldsymbol{y}\right)$近似 $A$ 的特征对。这样就将大型特征值问题 $\eqref{eq6}$ 转化为一个小型特征值问题:

$$
\boldsymbol{H}_ {k} \boldsymbol{y}=\lambda \boldsymbol{y}
$$

这就是用 Krylov 子空间方法求解大规模特征值问题的基本思想。

### 13. 数值代数理论基础

**1.正定矩阵定义**

设 $A \in \mathbb{R}^{n \times n}$， 若对任意的非零向量 $x \in \mathbb{R}^{n}$ 有$(\boldsymbol{A} \boldsymbol{x}, \boldsymbol{x})>0(\geqslant 0)$，则称 $\boldsymbol{A}$ 为 $($ 实 $)$ 正定 矩阵 $($ 半正定矩阵 $)$， 若 $\boldsymbol{A}$还是对称的，则称为 $($ 实 $)$ 对称正定矩阵 $($ 对称半正定矩阵 $)$。

**2.Schur分解**

设 $A \in \mathbb{C}^{n \times n}$，则存在酉矩阵 $\boldsymbol{P} \in \mathbb{C}^{n \times n}$ 使得

$$
\boldsymbol{P}^{\mathrm{H}} \boldsymbol{A} \boldsymbol{P}=\boldsymbol{T}
$$

式中：$\boldsymbol{T}$ 为上三角矩阵，其对角元素是 $\boldsymbol{A}$ 的特征值，而且可以选取 $\boldsymbol{P}$ 使得 $\boldsymbol{T}$ 的对角元可以任意排列。

**3.奇异值定义**

设 $A \in \mathbb{C}_ {r}^{m \times n}(r \geqslant 1)$，记 Hermite 矩阵 $A^{\mathrm{H}} \boldsymbol{A}$ 的 $n$ 个特征值为 $\lambda_ {1} \geqslant \lambda_ {2} \geqslant \cdots \lambda_ {r}>\lambda_ {r+1}=\cdots=\lambda_ {n}=0$，则称

$$
\sigma_ {i}=\sqrt{\lambda_ {i}}, \quad i=1,2, \cdots, n
$$

为 $A$ 的奇异值。

**4.奇异值分解定理**

设 $A \in \mathbb{C}_ {r}^{m \times n}$ 的正奇异值为 $\sigma_ {1}, \sigma_ {2}, \cdots, \sigma_ {r}$。则存在 $m$ 阶酉矩阵 $\boldsymbol{U}$ 和 $n$ 阶酉矩阵 $\boldsymbol{V}$，使得

$$
\boldsymbol{A}=\boldsymbol{U} \boldsymbol{\Sigma} \boldsymbol{V}^{\mathrm{H}}
$$

式中：$\boldsymbol{\Sigma}$ 为 $m \times n$ 阶对角矩阵，且这里 $\boldsymbol{\Sigma}_ {r}=\operatorname{diag}\left(\sigma_ {1}, \sigma_ {2}, \cdots, \sigma_ {r}\right), \sigma_ {i}>0,1 \leqslant i \leqslant r, r \leqslant \min \{m, n\} $。

**5.矩阵范数**

$$
\|\boldsymbol{A}\|_ {\mathrm{F}}=\left(\sum_ {i=1}^{m} \sum_ {j=1}^{n}\left|a_ {i j}\right|^{2}\right)^{1 / 2}=
\left(\sum_ {i=1}^{n} \sigma_ {i}^{2}\right)^{1 / 2}
$$

$$
\|\boldsymbol{A}\|_ {2}=\sigma_ {1}=\sqrt{\lambda_ {\max }\left(\boldsymbol{A}^{\mathrm{H}} \boldsymbol{A}\right)}
$$

**6.谱半径**

设 $A \in \mathbb{C}^{n \times n}$，其特征值为 $\lambda_ {1}, \lambda_ {2}, \cdots, \lambda_ {n}$，则称

$$
\rho(\boldsymbol{A})=\max _ {1 \leqslant i \leqslant n}\left|\lambda_ {i}\right|
$$

为矩阵 $A$ 的谱半径。由上述定义，$\|\boldsymbol{A}\|_ {2}$ 可定义为

$$
\|\boldsymbol{A}\|_ {2}=\sqrt{\rho\left(\boldsymbol{A}^{\mathrm{H}} \boldsymbol{A}\right)}
$$

特别地，当 $A$ 为 Hermite 矩阵时，有

$$
\|\boldsymbol{A}\|_ {2}=\rho(\boldsymbol{A}) .
$$

$\boldsymbol{A}$ 的谱半径不超过 $\boldsymbol{A}$ 的任何范数，即

$$
\rho(\boldsymbol{A}) \leqslant\|\boldsymbol{A}\|.
$$

**7.满秩分解**

设 $A \in \mathbb{C}_ {r}^{m \times n}(r>0)$, 若存在列满秩矩阵 $\boldsymbol{F} \in \mathbb{C}_ {r}^{m \times r}$ 和行满秩矩阵 $G \in \mathbb{C}_ {r}^{r \times n}$，使得 $\boldsymbol{A}=\boldsymbol{F} G$，则称 $\boldsymbol{F} \boldsymbol{G}$ 为 $\boldsymbol{A}$ 的
一个满秩分解。

$$
\boldsymbol{A}=\left[\begin{array}{cccc}-1 & 0 & 1 & 2 \\ 1 & 2 & -1 & 1 \\ 2 & 2 & -2 & -1\end{array}\right]
$$

由

$$
[\boldsymbol{A}, \boldsymbol{I}]=\left[\begin{array}{ccccccc}
-1 & 0 & 1 & 2 & 1 & 0 & 0 \\
1 & 2 & -1 & 1 & 0 & 1 & 0 \\
2 & 2 & -2 & -1 & 0 & 0 & 1
\end{array}\right] \stackrel{\text { 行 }}{\longrightarrow}\left[\begin{array}{cccccc}
-1 & 0 & 1 & 2 & 1 & 0 & 0 \\
0 & 2 & 0 & 3 & 1 & 1 & 0 \\
0 & 0 & 0 & 0 & 1 & -1 & 1
\end{array}\right]
$$

可知

$$
\boldsymbol{P}=\left[\begin{array}{lcc}
1 & 0 & 0 \\
1 & 1 & 0 \\
1 & -1 & 1
\end{array}\right], \quad \boldsymbol{P}^{-1}=\left[\begin{array}{ccc}
1 & 0 & 0 \\
-1 & 1 & 0 \\
-2 & 1 & 1
\end{array}\right]
$$

于是有

$$
\boldsymbol{F}=\left[\begin{array}{cc}
1 & 0 \\
-1 & 1 \\
-2 & 1
\end{array}\right], \quad \boldsymbol{G}=\left[\begin{array}{cccc}
-1 & 0 & 1 & 2 \\
0 & 2 & 0 & 3
\end{array}\right], \quad \boldsymbol{A}=\boldsymbol{F} \boldsymbol{G}
$$

**8.矩阵广义逆的计算**

> 设$A$ 的最大秩分解 $A=FG$ ,  则 $\boldsymbol{A}^{\dagger}=\boldsymbol{G}^{\mathrm{H}}\left(\boldsymbol{G}\boldsymbol{G}^{\mathrm{H}}\right)^{-1} \left(\boldsymbol{F}^{\mathrm{H}}\boldsymbol{F}\right)^{-1}\boldsymbol{F}^{\mathrm{H}}$。

已知

$$
\boldsymbol{A}=\left[\begin{array}{llll}1 & 2 & 1 & 2 \\ 0 & 1 & 0 & 1 \\ 1 & 0 & 1 & 0 \\ 2 & 1 & 2 & 1\end{array}\right]
$$

求 $\boldsymbol{A}$ 的满秩分解和 $A^{\dagger}$。

解

$$
\quad \boldsymbol{A}^{\text {初等行变换 }}\left[\begin{array}{llll}1 & 0 & 1 & 0 \\ 0 & 1 & 0 & 1 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0\end{array}\right]: c_ {1}=1, c_ {2}=2
$$

$$
\boldsymbol{A}=\boldsymbol{F} \boldsymbol{G}: \quad \boldsymbol{F}=\left[\begin{array}{ll}1 & 2 \\ 0 & 1 \\ 1 & 0 \\ 2 & 1\end{array}\right], \quad \boldsymbol{G}=\left[\begin{array}{llll}1 & 0 & 1 & 0 \\ 0 & 1 & 0 & 1\end{array}\right]
$$

因此，有

$$
\boldsymbol{F}^{\dagger}=\left(\boldsymbol{F}^{\mathrm{T}} \boldsymbol{F}\right)^{-1} \boldsymbol{F}^{\mathrm{T}}=\left[\begin{array}{ll}6 & 4 \\ 4 & 6\end{array}\right]^{-1} \boldsymbol{F}^{\mathrm{T}}=\frac{1}{10}\left[\begin{array}{cccc}-1 & -2 & 3 & 4 \\ 4 & 3 & -2 & -1\end{array}\right]
$$

$$
\boldsymbol{G}^{\dagger}=\boldsymbol{G}^{\mathrm{T}}\left(\boldsymbol{G} \boldsymbol{G}^{\mathrm{T}}\right)^{-1}=\boldsymbol{G}^{\mathrm{T}}\left[\begin{array}{ll}2 & 0 \\ 0 & 2\end{array}\right]^{-1}=\frac{1}{2}\left[\begin{array}{ll}1 & 0 \\ 0 & 1 \\ 1 & 0 \\ 0 & 1\end{array}\right]
$$

$$
\boldsymbol{A}^{\dagger}=\boldsymbol{G}^{\dagger} \boldsymbol{F}^{\dagger}=\frac{1}{20}\left[\begin{array}{cccc}-1 & -2 & 3 & 4 \\ 4 & 3 & -2 & -1 \\ -1 & -2 & 3 & 4 \\ 4 & 3 & -2 & -1\end{array}\right]
$$