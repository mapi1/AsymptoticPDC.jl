# Multichannel Autoregression Model

**Definition:** A multichannel or multivariate autoregressive model (MCAR) of order $p$ is defined as

$X(t) = \sum_{i=1}^{p} A(i)X(t-i) + E(t)$

where $X(t) = [X_1(t), X_2(t), ...,X_k(t)]$ is a k-channel set of signals and $E(t)$ a vector of k withe noises at time $t$. $A$ is a $k  \times k \times p$ matrix holding the model parameters.


```@docs
mcar
```

```@docs
MCAR_Model
```