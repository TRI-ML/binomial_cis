# binomial_cis
This package computes confidence intervals for the probability of success parameter, $p$, of a binomial distribution.
The confidence intervals computed by this package cover $p$ with exactly the user-specified probability and have minimal excess length.


## Installation
Install the package with pip:
```
pip install git+https://github.com/TRI-ML/binomial_cis.git
```



## What is a binomial confidence interval?
The binomial distribution represents the likelihood of observing $k$ successes in $n$ trials where the probability of success for each trial is $p$. One often does not know the true value of $p$ and wishes to estimate this value. After observing $k$ successes in $n$ trials with unknown probability of success $p$, a confidence interval (CI) is constructed in such a way that it contains the true value of $p$ with some high probability.

In constructing confidence intervals one has to tradeoff between three quantities:
1. Confidence: The probability that the CI contains the true parameter. Often written as $1-\alpha$ where $\alpha$ is small.
2. Volume: The length of the CI. If the CI is constructed as a lower bound on $p$ (i.e. $[\underline{p}, 1]$), then we care about the length of the CI which is below $p$. This is known as the *shortage*. 
3. Number of samples: In general, with more samples one can construct CIs with higher confidence and smaller volume.




## Why does this package exist?
Existing implementations of binomial CIs fail to optimally control the tradeoffs between coverage, volume, and number of samples. What this means in practice is that if a user specifies $k$ and $n$, existing implementations return CIs with more/less coverage than desired and/or CIs with higher volume than necessary.

Existing software implementations for binomial CIs include:
- [statsmodels.stats.proportion.proportion_confint](https://www.statsmodels.org/devel/generated/statsmodels.stats.proportion.proportion_confint.html) (Python)
- [scipy.stats._result_classes.BinomTestResult.proportion_ci](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats._result_classes.BinomTestResult.proportion_ci.html) (Python)
- [astropy.stats.binom_conf_interval](https://docs.astropy.org/en/stable/api/astropy.stats.binom_conf_interval.html) (Python)
- [scipy.stats.binom](https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.binom.html) (Python)
- [EBCIC](https://github.com/KazKobara/ebcic) (Python)
- [binom.test](https://stat.ethz.ch/R-manual/R-devel/library/stats/html/binom.test.html) (R)
- [binom.confint](https://cran.r-project.org/web/packages/binom/binom.pdf) (R)
- [BinomCI](https://search.r-project.org/CRAN/refmans/DescTools/html/BinomCI.html) (R)
- [HypothesisTests.jl](https://juliastats.org/HypothesisTests.jl/stable/methods/#Confidence-interval) (Julia)
- [RobustStats.jl](https://github.com/mrxiaohe/RobustStats.jl/tree/67cec0ea9ac9c4c8a3fe8f2a4e93cadbe68b1d77) (Julia)
- [ClinicalTrialUtilities.jl](https://docs.juliahub.com/ClinicalTrialUtilities/lEJUs/0.3.1/) (Julia)
- [binofit](https://www.mathworks.com/help/stats/binofit.html) (Matlab)

The methods these packages implement include:
- Wald/Normal [1]
- Agresti-Coull [1]
- Clopper-Pearson (*) [1]
- Wilson [1]
- Modified Wilson [1]
- Wilson with continuity correction (*) [1]
- Jeffreys [1]
- Bayesian uniform prior
- Inverting the binomial test (*) [3]
- Arcsine [1]
- Logit [1]
- Probit [8]
- Complementary log [8]
- Likelihood (Profile) [1]
- Witting (*) [4]
- Pratt [5]
- Mid-p [7]
- Blaker (*) [6]
- Second-order corrected [2]

The reference of each method points to a survey paper (if possible) rather than the original derivation of the method.

Only the methods marked with (*) are guaranteed to provide at least as much coverage as desired. However, none of these methods provide exactly the coverage desired (Witting might, but the listed reference is not freely available online and is only published in German).

It is also worth noting the [ump](http://CRAN.R-project.org/package=ump) R package. This package implements UMP and UMPU hypothesis tests for the binomial distribution. Such tests could be leveraged to construct UMA and UMAU confidence intervals, but it this doesn't appear to be implemented based on the documented functions. In addition, in the documentation the authors of this package note that their implemetation has issues with numerical stability.

Existing software implementations for computing binomial CIs with exact coverage are then either non-existent or unsatisfactory. Unsurprisingly then, there is also no open-source implementation of computing the expected shortage (or expected excess or expected width) of UMA and UMAU CIs, and no implementation for computing the worst-case values of these quantities. This package exists to fill this gap.









## **What exactly does this package do?**
This package constructs optimal confidence intervals for the probability of success parameter of a binomial distribution.

#### **Lower Bounds**
Given user specified miscoverage rate ($\alpha$) and maximum expected shortage ($\text{MES}$), return a lower bound on $p$ that satisfies the following requirements:
1. achieves exact desired coverage: $\mathbb{P}[\underline{p} \le p] = 1-\alpha$
2. $[\underline{p}, 1]$ is uniformly most accurate
3. achieves exact desired maximum expected shortage: $\max_p \ \mathbb{E}_p[\max (p - \underline{p}, 0)] = \text{MES}$
4. uses the minimum number of samples $n$ to achieve requirements 1,2,3.


#### **Upper Bounds**
Given user specified miscoverage rate ($\alpha$) and maximum expected excess ($\text{MEE}$), return an upper bound on $p$ that satisfies the following requirements:
1. achieves exact desired coverage: $\mathbb{P}[p \le \overline{p}] = 1-\alpha$
2. $[0, \overline{p}]$ is uniformly most accurate
3. achieves exact desired maximum expected excess: $\max_p \ \mathbb{E}_p[\max (\overline{p} - p, 0)] = \text{MEE}$
4. uses the minimum number of samples $n$ to achieve requirements 1,2,3.




#### **Simultaneous Lower and Upper Bounds**
Given the user specified miscoverage rate ($\alpha$) and maximum expected width ($\text{MEW}$), return simultaneous lower and upper bounds on $p$ that satisfy the following requirements:
1. achieves exact desired coverage: $\mathbb{P}[\underline{p} \le p \le \overline{p}] = 1-\alpha$
2. $[\underline{p}, \overline{p}]$ is uniformly most accurate unbiased
3. achieves exact desired maximum expected width: $\max_p \ \mathbb{E}_p[\overline{p} - \underline{p}] = \text{MEW}$
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






## **Notebooks**
The `notebooks/` directory has the following notebooks:
- `tradeoff_table.ipynb`: Computes MES vs miscoverage rate $\alpha$ and number of samples $n$. Precomputed values have been stored in `MES_table.csv` which is visualized in a plot from the last cell of the notebook.
- `visualizations.ipynb`: Visualizes the mixed-monotonic forms of expected shortage and expected width. Also visualizes how these functions vary with $p$ and their maxima.

The `tests/` directory has the following notebooks:
- `binom_helper_validation.ipynb`: Tests our implementation of the binomial coefficient, binomial pmf, and binomial cdf against their SciPy counterparts.
- `conf_set_validation.ipynb`: Tests the theoretical guarantees of the CIs (coverage, shortage, excess, width, unbiasedness) using Monte-Carlo simulation.


## Automated Tests
Not implemented yet, but future automated tests may include
- validating binomial_helper.py
- comparing the bounds with those in the appendix of [Table of Neyman-shortest unbiased confidence intervals for the binomial parameter](https://www.jstor.org/stable/2333308) by Blyth and Hutchinson


## **Caution!**
#### **Randomized CIs**
The methods used in this package to construct CIs are based on the inversion of randomized hypothesis tests. This means that calling `binom_ci()` with the same `k,n,alpha` will result in different CIs. For the guarantees of the CI to hold it is critical that the user only construct one CI for the experiment they have. Constructing multiple CIs and choosing the best one invalidates the guarantees of the CI. 

For the 1-sided bounds there is the option to get less efficient but non-randomized CIs:
```
lb = binom_ci(k, n, alpha, 'lb', randomized=False)
ub = binom_ci(k, n, alpha, 'ub', randomized=False)
```

These non-randomized 1-sided bounds are equivalent to 1-sided Clopper-Pearson bounds. We currently don't have an implementation of non-randomized 2-sided bounds.

Randomization allows the CIs to be UMA. Although randomization has been a point of debate amongst statisticians, we take the view (first given by Mark Eudey) that insofar as construction of confidence intervals can be treated as a (von Neumann) game, randomization merely allows the statistician to employ a mixed strategy.



#### **Multiple Tests**
As with all CIs one must take special care when interpreting the results of multiple CIs constructed from independent tests. If one constructs $m$ CIs where the probability of each CI containing the true parameter is $1-\alpha$, then the probability that **all** $m$ CIs contain their respective parameters is less than $1-\alpha$. For more explanation, see the [Wikipedia article](https://en.wikipedia.org/wiki/Multiple_comparisons_problem) on the multiple comparisons problem.



## Building the Package
Activate the virtual environment and run
```
python -m build
```
You will need the `build` package for this.





## Relevant Literature
Below are some of the papers that we found most useful for understanding binomial confidence intervals.
- [Testing Statistical Hypotheses](https://link.springer.com/book/10.1007/978-3-030-70578-7) by Lehmann and Romano
- [On the treatment of discontinuous random variables](https://search.library.berkeley.edu/permalink/01UCS_BER/s4lks2/cdi_proquest_journals_301822201) by Eudey
- [Table of Neyman-shortest unbiased confidence intervals for the binomial parameter](https://www.jstor.org/stable/2333308) by Blyth and Hutchinson
- [Length of Confidence Intervals](https://www.tandfonline.com/doi/abs/10.1080/01621459.1961.10480644) by Pratt
- [More on length of confidence intervals](https://www.tandfonline.com/doi/abs/10.1080/01621459.1962.10500547) by Madansky
- [Binomial confidence intervals](https://www.tandfonline.com/doi/abs/10.1080/01621459.1983.10477938) by Blyth and Still
- [Smallest confidence intervals for one binomial proportion](https://www.sciencedirect.com/science/article/pii/S0378375805002430) by Wang
- [Fuzzy and randomized confidence intervals and p-values](https://www.jstor.org/stable/20061193) by Geyer and Meeden
- [Nonoptimality of Randomized Confidence Sets](https://www.stat.purdue.edu/docs/research/tech-reports/1988/tr88-09.pdf) by Casella and Robert









## **References**
1. [*Interval Estimation for a Binomial Proportion*](https://www.jstor.org/stable/2676784) by Brown, Cai, and DasGupta
2. [*One-sided confidence intervals in discrete distributions*](https://www.sciencedirect.com/science/article/pii/S0378375804000679) by Cai
3. [*Some Remarks on Confidence or Fiducial Limits*](https://www.jstor.org/stable/2333026) by Sterne
4. [*Mathematische Statistik I.*](https://link.springer.com/book/10.1007/978-3-322-90150-7) by Witting
5. [*Binomial Confidence Intervals*](https://www.jstor.org/stable/2287116) by Blyth and Still
6. [*Confidence Curves and Improved Exact Confidence Intervals for Discrete Distributions*](https://www.jstor.org/stable/3315916) by Blaker
7. [*Comment: Randomized Confidence Intervals and the Mid-P Approach*](https://www.semanticscholar.org/paper/Comment%3A-Randomized-Confidence-Intervals-and-the-Agresti-Gottard/596fd2a98486f65599b56ac81aef32e5c1a0d586) by Agresti and Gottard
8. [*binom*](https://cran.r-project.org/web/packages/binom/binom.pdf) by Sundar Dorai-Raj
9. [*Fuzzy and randomized confidence intervals and p-values*](https://www.jstor.org/stable/20061193) by Geyer and Meeden.
