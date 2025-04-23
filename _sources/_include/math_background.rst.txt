Mathematical Background
=======================


Uniformly Most Accurate (UMA)
*****************************
The standard notion of optimality for 1-sided confidence bounds is that of being uniformly most accurate (UMA).
Here we present the UMA condition for a lower confidence bound (Eq. 3.22 of [Lehmann]_), and the condition for upper confidence bounds is analagous.
For test statistic :math:`T`, a lower confidence bound :math:`\underline{p}(T)` satisfying

.. math::

    \mathbb{P}_p [\underline{p}(T) \le p] \ge 1-\alpha \quad \forall p

and, amongst all possible lower bounds :math:`\underline{p}^\prime(T)`, minimizes

.. math::

    \mathbb{P}_p [\underline{p}^\prime(T) \le p_0] \quad \forall p_0 < p

is a UMA lower confidence bound at confidence level :math:`1-\alpha`.




Uniformly Most Accurate Unbiased (UMAU)
***************************************
A UMAU 2-sided confidence interval satisfies the property of unbiasedness (Eq. 5.39 of [Lehmann]_):

.. math::

    \mathbb{P}_p [\underline{p}(T) \le p \le \overline{p}(T)] \ge  \mathbb{P}_p [\underline{p}(T) \le p^\prime \le \overline{p}(T)] \quad \forall p^\prime \neq p

and amongst all unbiased confidence intervals is UMA.
Being unbiased ensures that true unknown parameter has the highest probability of being included in the interval.


UMA Lower Bound
***************
The background for the UMA lower bound is from [Vincent]_.
Let :math:`\textup{bin}(k,n,p)` and :math:`\textup{Bin}(k,n,p)` denote the binomial probability mass function (PMF) and CDF, with :math:`k` successes, :math:`n` trials, and success probability :math:`p`.
The floor function :math:`\lfloor x \rfloor` rounds the argument :math:`x` down to the nearest integer.
Now, define a test statistic :math:`T = U + \sum_{i=1}^n X_i` where :math:`U \sim \mathcal{U}[0,1]` and :math:`X_i` are the observed successes and failures. 
Then :math:`T` is random with probability density

.. math::

    f_p(t) = {n \choose \lfloor t \rfloor} p^{\lfloor t \rfloor} (1-p)^{n - \lfloor t \rfloor}.


Then, we can compute the CDF as follows,

.. math::

    F_p(t) = \int_0^{t} {n \choose \lfloor x \rfloor} p^{\lfloor x \rfloor} (1-p)^{n - \lfloor x \rfloor} dx \\
    = \textup{Bin}(\lfloor t \rfloor - 1, n, p) + (t - \lfloor t \rfloor) \cdot \textup{bin}(\lfloor t \rfloor, n, p). 


:math:`F_p(t)` is continuous in both :math:`t` and :math:`p`.

Then, by Corollary 3.5.1 of [Lehmann]_, we have the UMA lower confidence bound,

.. math::

    \underline{p}(t) = \begin{cases}
                        0 \quad &\text{if } t < 1-\alpha \\
                        1 \quad &\text{if } t > n + 1-\alpha \\
                        p^* \quad &\text{s.t. } F_{p^*}(t) = 1-\alpha \text{ otherwise.}
    \end{cases}

That satisfies :math:`\mathbb{P}[\underline{p} \le p] = 1-\alpha`

The bisection method can be used to solve :math:`F_{p^*}(t) = 1-\alpha`, as :math:`F_p(t)` is monotonically decreasing in :math:`p`. 
Taking limits :math:`p \rightarrow 0^+` and :math:`p \rightarrow 1^-`, one can show that :math:`F_p(t) = 1-\alpha` has no solution when :math:`t<1-\alpha` and :math:`t>n+1-\alpha`, we thus set :math:`\underline{p}` to either zero or one in these cases. 
This bound is well-established (see Example 3.5.2 of [Lehmann]_).
The Clopper-Pearson bound can be obtained by setting :math:`U = 0` in the test statistic.

Next, (employing the Ghosh-Pratt identity) we compute expected shortage. A useful identity is

.. math::

    \mathbb{P}[\underline{p} \le p_0] = \mathbb{P}[t \le t^*(p_0)]

where :math:`t^*(p_0)` is the unique value that satisfies

.. math::

    F_{p_0}(t^*) = 1-\alpha.

Then, we can compute expected shortage as

.. math::

    \textup{ES}(p) = \int_0^p \mathbb{P}_p[\underline{p} \le p_0] dp_0 = \int_0^p \mathbb{P}_p[t \le t^*(p_0)] dp_0 = \int_0^p F_p(t^*(p_0)) dp_0.


We evaluate this integral standard numerical software.

Next, we find MES via global optimization. Specifically, consider the following reformulation of expected shortage,

.. math::

    \textup{ES}(p_1, p_2) = \int_0^{p_1} F_{p_2}(t^*(p_0)) dp_0.

In this form, expected shortage is increasing in :math:`p_1` and decreasing in :math:`p_2`; this property is known as mixed monotonicity.
Note that for :math:`p_1 = p_2` we recover the original function.
Mixed monotonic functions can be globally optimized using the branch and bound procedure given in [Matthiesen]_.



UMA Upper Bound
***************
We construct a UMA upper bound on the probability of success by first computing a UMA lower bound on the probability of failure and then applying the identity

.. math::

    \overline{p} = 1 - \underline{q} \\
    q = 1 - p

To construct a UMA lower bound on the probability of failure, :math:`\underline{q}`, we simply use the procedure mentioned above but with the test statistic counting failures instead of successes.


UMAU 2-Sided Bound
******************
From [Blyth]_ we know that a parameter value is in the 2-sided confidence interval (:math:`p_0 \in [\underline{p}, \overline{p}]`) if and only if 

.. math::
    n_l + \gamma_l \le K + U \le n_u + \gamma_u

where :math:`K` is the number of observed successes in :math:`n` trials, :math:`U` is a uniform random variable (thus the test statistic is :math:`T = K + U`), and :math:`n_l, n_u` are integers and :math:`\gamma_l, \gamma_u \in [0,1]` that uniquely satisfy the relationships:


.. math::
    \begin{gather}
    0 \le n_l \le n_u \le n+1 \\
    \gamma_l = \frac{(n_u - np_0)(\mathbb{P}_{p_0}[n_l \le K \le n_u -1] - (1-\alpha)) - n_l(1-p_0)\mathbb{P}_{p_0}[K=n_l] + n_u(1-p_0)\mathbb{P}_{p_0}[K=n_u]}
    {(n_u - n_l) \mathbb{P}_{p_0}(K = n_l)} \\
    \gamma_u = \frac{(n_l - np_0)(\mathbb{P}_{p_0}[n_l \le K \le n_u -1] - (1-\alpha)) - n_l(1-p_0)\mathbb{P}_{p_0}[K=n_l] + n_u(1-p_0)\mathbb{P}_{p_0}[K=n_u]}
    {(n_u - n_l) \mathbb{P}_{p_0}(K = n_u)}
    \end{gather}

Thus, given a test statistic :math:`T`, we find :math:`\underline{p}` by finding the smallest value of :math:`p_0` for which the above conditions have a solution.
Similarly, we find :math:`\overline{p}` by finding the largest value of :math:`p_0` for which the above conditions have a solution.

Then, (employing the Ghosh-Pratt identity), we can compute expected width as

.. math::
    \begin{gather}
    \textup{EW}(p) = \int_0^1 \mathbb{P}_p[\underline{p} \le p_0 \le \overline{p}] dp_0 = \int_0^1 \mathbb{P}_p[n_l + \gamma_l \le t \le n_u + \gamma_u] dp_0 \\
     = \int_0^1 F_p(n_u + \gamma_u) - F_p(n_l + \gamma_l) dp_0
    \end{gather}

Similar to expected shortage (and expected excess), expected width also has a mixed monotonic form.

.. math::
    \begin{gather}
    \textup{EW}(p_1, p_2) = \int_0^1 F_{p_2}(n_u + \gamma_u) - F_{p_1}(n_l + \gamma_l) dp_0.
    \end{gather}

Thus, we can compute the maximum expected width (MEW) using mixed monotonic programming.





Clopper-Pearson Lower Bound
***************************
In the ``tradeoff_table.ipynb`` notebook we compare the MES of the UMA lower bound with that of the Clopper-Pearson lower bound.
Computation of MES for the Clopper-Pearson bound can be accomplished in a similar manner to that of the UMA bound.
The Clopper-Pearson lower confidence bound is computed as

.. math::
    \underline{p}_{cp} = \min \{p \mid \textup{Bin}(k-1, n, p) \ge 1-\alpha \}

where :math:`k` is the number of successes and :math:`n` is the number of trials.
It follows that

.. math::
    \begin{gather}
    \underline{p}_{cp} \le p_0 \iff k \le t_{cp}^*(p_0) + 1 \\
    t_{cp}^*(p_0) = \max \{t \in [0,1,\ldots,n] \mid \textup{Bin}(t, n, p_0) \le 1-\alpha \}.
    \end{gather}

Thus, (employing the Ghosh-Pratt identity) the expected shortage for the Clopper-Pearson lower bound is

.. math::

    \textup{ES}_{cp}(p) = \int_0^p \mathbb{P}_p[\underline{p}_{cp} \le p_0] dp_0 = \int_0^p \mathbb{P}_p[k \le t_{cp}^*(p_0)] dp_0 = \int_0^p \textup{Bin}(t_{cp}^*(p_0), n, p) dp_0.

It is best to evaluate this integral without integration software because the integrand is piecewise-constant.
In addition, the expected shortage for the Clopper-Pearson lower confidence bound also has a mixed monotonic form.

.. math::

    \textup{ES}_{cp}(p_1, p_2) = \int_0^{p_1} \textup{Bin}(t_{cp}^*(p_0), n, p_2) dp_0.

We can use this form to find the MES for the Clopper-Pearson lower bound, as done in the ``tradeoff_table.ipynb`` notebook.

References
**********
.. [Vincent] Vincent, Joseph A., et al. "How Generalizable Is My Behavior Cloning Policy? A Statistical Approach to Trustworthy Performance Evaluation." arXiv preprint arXiv:2405.05439 (2024).
.. [Lehmann] E. L. Lehmann and J. P. Romano, Testing Statistical Hypotheses. Springer, 2022, vol. 4.
.. [Matthiesen] B. Matthiesen, C. Hellings, E. A. Jorswieck, and W. Utschick, “Mixed Monotonic Programming for Fast Global Optimization,” IEEE Transactions on Signal Processing, vol. 68, pp. 2529–2544, 2020.
.. [Blyth] Blyth, Colin R., and David W. Hutchinson. "Table of Neyman-shortest unbiased confidence intervals for the binomial parameter." Biometrika 47.3/4 (1960): 381-391.

