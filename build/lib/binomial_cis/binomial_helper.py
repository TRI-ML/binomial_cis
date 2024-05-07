import numpy as np
from numba import cfunc

@cfunc("float64(float64, float64)")
def binom_coeff(n, k):
    """
    Inputs
    n: total number of items
    k: number being chosen
    
    Returns
    c: n choose k, the binomial coefficient
    """
    # Use the multiplicative formula for computational efficiency:
    # https://en.wikipedia.org/wiki/Binomial_coefficient#Multiplicative_formula
    c = 1
    upper_lim = int(min(k, n-k))
    for i in range(1, upper_lim + 1):
        c *= (n + 1 - i) / i
    return c


@cfunc("float64(float64, float64, float64)")
def binom_pmf(k, n, p):
    """
    Inputs
    k: number of successes
    n: number of trials
    p: probability of success
    
    Returns
    pmf: binomial pmf evaluated at k, n, p
    """
    pmf = binom_coeff(n, k) * p**k * (1-p)**(n-k)
    
    # safeguard against floating point errors
    pmf = min(1, max(0, pmf)) # unable to use np.clip()
    return pmf


@cfunc("float64(float64, float64, float64)")
def binom_cdf(k, n, p):
    """
    Inputs
    k: number of successes
    n: number of trials
    p: probability of success
    
    Returns
    cdf: binomial cdf evaluated at k, n, p
    """
    # code from https://stackoverflow.com/a/45869209
    if k == n or p == 0:
        cdf = 1
    elif p == 1:
        cdf = 0
    else:
        cdf = 0
        b = 0
        for i in range(int(k) + 1):
            if i > 0:
                b += np.log(n + 1 - i) - np.log(i) 
            log_pmf_i = b + i * np.log(p) + (n-i) * np.log(1-p)
            cdf += np.exp(log_pmf_i)
    
    # safeguard against floating point errors
    cdf = min(1, max(0, cdf)) # unable to use np.clip()
    return cdf

