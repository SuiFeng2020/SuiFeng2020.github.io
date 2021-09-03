---
layout: post
title:  Eckart-Young theorem
date:   2021-09-03 
description: "矩阵的SVD低秩近似"
tag: learning
---

For a matrix $A \in \mathbb{R}^{m \times n}$, it must have an SVD decomposition $A=U \Sigma V^{T}$ with singular values $\sigma_ {1} \geq \sigma_ {2} \geq \ldots \geq \sigma_ {p} \geq 0$ and $p=\min \{m, n\}$ . Using $\left\|U^{T} A V\right\|=\|\Sigma\|$ and $U, V$ are orthogonal, we can express both 2 -norm and Frobenius norm in terms of singular values

$$
\|A\|_ {2}=\sigma_ {1}, \quad\|A\|_ {F}=\sqrt{\sigma_ {1}^{2}+\cdots+\sigma_ {p}^{2}}
$$

since by definition

$$
\begin{gathered}
\|A\|_ {2}=\|\Sigma\|_ {2}=\max _ {\|\mathbf{x}\|=1}\|\Sigma \mathbf{x}\|_ {2}=\max \left\{\left|\sigma_ {i}\right|\right\}=\sigma_ {\max }(A)=\sigma_ {1} \\
\|A\|_ {F}=\|\Sigma\|_ {F}=\sqrt{\sigma_ {1}^{2}+\cdots+\sigma_ {p}^{2}}
\end{gathered}
$$

The Frobenius norm can also be shown using the equivalent definition

$$
\|A\|_ {F}^{2}=\operatorname{Tr}\left(A^{T} A\right)=\operatorname{Tr}\left(\left(U \Sigma V^{T}\right)^{T}\left(U \Sigma V^{T}\right)\right)=\operatorname{Tr}\left(\Sigma^{T} \Sigma\right)=\sigma_ {1}^{2}+\cdots+\sigma_ {p}^{2}
$$

In reality, we would like to use fewer numbers to represent the large matrix $A$ as in data compression, the best low-rank approximation can be obtained from SVD, which is shown below.

## The Eckart-Young Theorem

Suppose a matrix $A \in \mathbb{R}^{m \times n}$ has an SVD-decomposition $A=U \Sigma V^{T}$. Let $k<r=\operatorname{rank}(A)$ and truncated matrix

$$
A_ {k}=\sum_ {i=1}^{k} \sigma_ {i} \mathbf{u}_ {i} \mathbf{v}_ {i}^{T}
$$

then, for any matrix $B$ of rank $k$, the minimal error is achieved with $A_ {k}$ : $\min _ {\operatorname{rank}(B)=k}\|A-B\|_ {2}=\left\|A-A_ {k}\right\|_ {2}=\sigma_ {k+1}$
The same holds for Frobenius norm as well

$$
\min _ {\operatorname{rank}(B)=k}\|A-B\|_ {F}=\left\|A-A_ {k}\right\|_ {F}=\sqrt{\sigma_ {k+1}^{2}+\cdots+\sigma_ {p}^{2}}
$$

- ### Proof. (2-norm case)

Since $U A_ {k} V=\operatorname{diag}\left(\sigma_ {1}, \ldots, \sigma_ {k}, 0, \ldots, 0\right)$ it means that $A_ {k}$ is rank $k$. Moreover, $U^{T}\left(A-A_ {k}\right) V=\operatorname{diag}\left(0, \ldots, 0, \sigma_ {k+1}, \ldots, \sigma_ {p}\right)$ with the largest singular value is $\sigma_ {k+1}$ and thus $\left\|A-A_ {k}\right\|_ {2}=\sigma_ {k+1}$.
Now suppose some matrix $B \in \mathbb{R}^{m \times n}$ has rank $k$, then the dimension of the nullspace of $B$ is $n-k$, which means we can find orthonormal vectors $\mathbf{x}_ {1}, \ldots, \mathbf{x}_ {n-k}$ that span the nullspace, i.e. $null(B) = \operatorname{span} \left\{\mathbf{x}_ {1}, \ldots, \mathbf{x}_ {n-k}\right\}$ Since the sum of dimension of $null(B)$ and $\operatorname{span}\left\{\mathbf{v}_ {1}, \ldots, \mathbf{v}_ {k+1}\right\}$ (dimension $k+1)$ is $n+1>n$, so 

$$
\operatorname{span}\left\{\mathbf{x}_ {1}, \ldots, \mathbf{x}_ {n-k}\right\} \cap \operatorname{span}\left\{\mathbf{v}_ {1}, \ldots, \mathbf{v}_ {k+1}\right\} \neq\{\mathbf{0}\}
$$

There exists a nontrivial unit vector $\mathbf{z}(\|\mathbf{z}\|=1)$ in the intersection set. Using $B \mathbf{z}=\mathbf{0}$ and

$$
A \mathbf{z}=\sum_ {i=1}^{n} \sigma_ {i}\left(\mathbf{v}_ {i}^{T} \mathbf{z}\right) \mathbf{u}_ {i}=\sum_ {i=1}^{k+1} \sigma_ {i}\left(\mathbf{v}_ {i}^{T} \mathbf{z}\right) \mathbf{u}_ {i}
$$

we have

$$
\begin{gathered}
\|A-B\|_ {2}^{2}=\|A-B\|_ {2}^{2} \cdot\|\mathbf{z}\|_ {2}^{2} \geq\|(A-B) \mathbf{z}\|_ {2}^{2} \\
\quad=\|A \mathbf{z}\|_ {2}^{2}=\left\|U \Sigma V^{T} \mathbf{z}\right\|_ {2}^{2}=\left\|\Sigma V^{T} \mathbf{z}\right\|_ {2}^{2} \\
=\sum_ {i=1}^{k+1} \sigma_ {i}^{2}\left(\mathbf{v}_ {i}^{T} \mathbf{z}\right)^{2} \geq \sigma_ {k+1}^{2} \sum_ {i=1}^{k+1}\left(\mathbf{v}_ {i}^{T} \mathbf{z}\right)^{2}=\sigma_ {k+1}^{2}
\end{gathered}
$$

since $\sigma_ {1} \geq \cdots \geq \sigma_ {k+1}$ and $\sum_ {i=1}^{k+1}\left(\mathbf{v}_ {i}^{T} \mathbf{z}\right)^{2}=\|V \mathbf{z}\|_ {2}^{2}=\|\mathbf{z}\|_ {2}^{2}=1$

- ### Proof. (Frobenius norm case)

> Lemma: If $A, B \in \mathbb{R}^{m \times n}$, with $B$ having rank $k$, then $\sigma_ {k+i}(A) \leq \sigma_ {i}(A-B) \quad$ for all $i$

To prove the lemma, first consider the case $i=1$, we have proved that $\sigma_ {k+1}(A) \leq \sigma_ {1}(A-B)=\|A-B\|_ {2}$ in the 2-norm case. Then we do the general case:

$$
\begin{aligned}
\sigma_ {i}(A-B) &=\sigma_ {i}(A-B)+\sigma_ {1}\left(B-B_ {k}\right) \quad \text { since } B=B_ {k} \\
&=\sigma_ {1}\left(A-B-(A-B)_ {i-1}\right)+\sigma_ {1}\left(B-B_ {k}\right) \\
& \geq \sigma_ {1}\left(A-B-(A-B)_ {i-1}+B-B_ {k}\right) \\
&=\sigma_ {1}\left(A-(A-B)_ {i-1}-B_ {k}\right) \\
& \geq \sigma_ {1}\left(A-A_ {k+i-1}\right) \\
&=\sigma_ {k+i}(A)
\end{aligned}
$$

Here, on the second equal sign, we used that removing $\sigma_ {1}, \cdots, \sigma_ {i-1}$ of $(A-B)$ leaves $\sigma_ {i}(A-B)$ as the largest singular value of the remaining matrix, i.e. $\sigma_ {i}(X)=\sigma_ {1}\left(X-X_ {i-1}\right) ;$
on the first $\geq$ sign, we used the triangle inequality for the 2 -norm,

$$
\sigma_ {1}(X+Y) \leq \sigma_ {1}(X)+\sigma_ {1}(Y)
$$

and on the last $\geq$ sign, since $\operatorname{rank}\left((A-B)_ {i-1}+B_ {k}\right) \leq \operatorname{rank}\left((A-B+B)_ {i-1+k}\right)=\operatorname{rank}\left(A_ {k+i-1}\right)$ therefore the residue $A-A_ {k+i-1}$ has smaller $\sigma_ {1}$. Now we have proved the lemma.
Finally, using the lemma, we have for Frobenius norm,

$$
\left\|A-A_ {k}\right\|_ {F}^{2}=\sum_ {i=k+1}^{r} \sigma_ {i}(A)^{2} \leq \sum_ {i=1}^{r-k} \sigma_ {i}(A-B)^{2} \leq\|A-B\|_ {F}^{2}
$$