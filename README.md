# binomial_cis
This package computes confidence intervals for the probability of success parameter, $p$, of a binomial distribution.
The confidence intervals computed by this package cover $p$ with exactly the user-specified probability and have minimal excess length.

## Documentation
The documentation for the package can be found at: [https://tri-ml.github.io/binomial_cis](https://tri-ml.github.io/binomial_cis).


## Installation
Install the package with pip:
```
pip install binomial_cis
```


## **What does this package do?**
This package constructs optimal confidence intervals for the probability of success parameter, $p$, of a binomial distribution.

#### **Lower Bounds**
Given user specified miscoverage rate ($\alpha$) and maximum expected shortage ($\text{MES}$), return a lower bound on $p$ that satisfies the following requirements:
1. achieves exact desired coverage: $\mathbb{P}[\underline{p} \le p] = 1-\alpha$,
2. $[\underline{p}, 1]$ is uniformly most accurate,
3. achieves exact desired maximum expected shortage: $\max_p \ \mathbb{E}_p[\max (p - \underline{p}, 0)] = \text{MES}$,
4. uses the minimum number of samples $n$ to achieve requirements 1,2,3.


#### **Upper Bounds**
Given user specified miscoverage rate ($\alpha$) and maximum expected excess ($\text{MEE}$), return an upper bound on $p$ that satisfies the following requirements:
1. achieves exact desired coverage: $\mathbb{P}[p \le \overline{p}] = 1-\alpha$,
2. $[0, \overline{p}]$ is uniformly most accurate,
3. achieves exact desired maximum expected excess: $\max_p \ \mathbb{E}_p[\max (\overline{p} - p, 0)] = \text{MEE}$,
4. uses the minimum number of samples $n$ to achieve requirements 1,2,3.




#### **2-Sided Bounds**
Given the user specified miscoverage rate ($\alpha$) and maximum expected width ($\text{MEW}$), return simultaneous lower and upper bounds on $p$ that satisfy the following requirements:
1. achieves exact desired coverage: $\mathbb{P}[\underline{p} \le p \le \overline{p}] = 1-\alpha$,
2. $[\underline{p}, \overline{p}]$ is uniformly most accurate unbiased,
3. achieves exact desired maximum expected width: $\max_p \ \mathbb{E}_p[\overline{p} - \underline{p}] = \text{MEW}$,
4. uses the minimum number of samples $n$ to achieve requirements 1,2,3.


## **How do I use this package?**
### Lower Bounds

Find a lower bound on $p$:
```
from binomial_cis import binom_ci

k = 5 # number of successes
n = 10 # number of trials
alpha = 0.05 # miscoverage probability

lb = binom_ci(k, n, alpha, 'lb')

```

Find maximum expected shortage given miscoverage rate and number of samples:
```
mes_ub, mes_lb, p_lb, num_iters = max_expected_shortage(alpha, n, tol=1e-3)
```


### Upper Bounds

Find an upper bound on $p$:
```
from binomial_cis import binom_ci

k = 5 # number of successes
n = 10 # number of trials
alpha = 0.05 # miscoverage probability

ub = binom_ci(k, n, alpha, 'ub')

```

Find maximum expected excess given miscoverage rate and number of samples:
```
mee_ub, mee_lb, p_lb, num_iters = max_expected_excess(alpha, n, tol=1e-3)
```


### 2-Sided Bounds

Find simultaneous lower and upper bounds on $p$:
```
from binomial_cis import binom_ci

k = 5 # number of successes
n = 10 # number of trials
alpha = 0.05 # miscoverage probability

lb, ub = binom_ci(k, n, alpha, 'lb,ub')

```

Find maximum expected width given miscoverage rate and number of samples:
```
mew_ub, mew_lb, p_lb, num_iters = max_expected_width(alpha, n, tol=1e-3)
```




