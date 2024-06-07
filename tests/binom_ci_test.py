import numpy as np
import sys; sys.path.insert(0, '..')
from binomial_cis import binom_ci
from scipy.stats import binomtest
from math import isclose


# define range of test conditions
ns = np.array(range(1,50+1))
alphas = np.array([0.01, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99])








#########################################
############ Helper Functions ###########
#########################################

def scipy_cp(k, n, alpha, side):
    """
    k: number of successes
    n: number of samples
    alpha: miscoverage rate, P(p_ <= p) >= 1-alpha
    side: 'lb' for lower bound, 'ub' for upper bound

    Returns
    p_cp: Clopper-Pearson lower confidence bound
    """
    assert side == 'lb' or side == 'ub'

    alternative = 'greater' if side=='lb' else 'less'
    stat_model = binomtest(k=k, n=n, alternative=alternative)
    ci = stat_model.proportion_ci(confidence_level=1-alpha, method='exact')
    p_cp = ci.low if side=='lb' else ci.high

    return p_cp






#########################################
################# Tests #################
#########################################

def test_lb():
    """
    Tests 
    1. Our non-randomized lower confidence bound is equivalent to Clopper-Pearson
    2. Our randomized lower confidence bound is always greater than the non-randomized version
    """
    for n in ns:
        print("\nn:", n)
        for alpha in alphas:
            print("   alpha:", alpha)
            p_lb_prev = 0.0
            for k in range(n):
                print("      k:", k)
                p_lb = binom_ci(k, n, alpha, side='lb', randomized=True, verbose=False)
                p_lb_cp = binom_ci(k, n, alpha, side='lb', randomized=False, verbose=False)
                p_lb_cp_scipy = scipy_cp(k, n, alpha, side='lb')
                
                # check agreement up to 4 decimals
                assert isclose(p_lb_cp, p_lb_cp_scipy, abs_tol=10**-4)

                # check p_lb is between clopper pearson bounds
                assert round(p_lb, 4) >= round(p_lb_cp, 4)
                assert round(p_lb_prev, 4) <= round(p_lb_cp, 4)
                p_lb_prev = p_lb



def test_ub():
    """
    Tests 
    1. Our non-randomized upper confidence bound is equivalent to Clopper-Pearson
    2. Our randomized upper confidence bound is always less than the non-randomized version
    """
    for n in ns:
        print("\nn:", n)
        for alpha in alphas:
            print("   alpha:", alpha)
            p_ub_prev = 1.0
            for k in range(n, -1, -1):
                print("      k:", k)
                p_ub = binom_ci(k, n, alpha, side='ub', randomized=True, verbose=False)
                p_ub_cp = binom_ci(k, n, alpha, side='ub', randomized=False, verbose=False)
                p_ub_cp_scipy = scipy_cp(k, n, alpha, side='ub')
                
                # check up to 4 decimals
                assert isclose(p_ub_cp, p_ub_cp_scipy, abs_tol=10**-4)
                
                # check p_lb is between clopper pearson bounds
                assert round(p_ub, 4) <= round(p_ub_cp, 4)
                assert round(p_ub_prev, 4) >= round(p_ub_cp, 4)
                p_ub_prev = p_ub



def test_paper_data():
    """
    Test that functions don't error on the data used in the paper, 
    'How Generalizable Is My Behavior Cloning Policy? A Statistical Approach to Trustworthy Performance Evaluation'
    """
    
    # data from hardware experiment
    n = 50
    k1 = 38
    k2 = 4
    lb_pour_good = binom_ci(k1, 50, 0.05, 'lb', verbose=False)
    lb_pour_bad = binom_ci(k2, 50, 0.05, 'lb', verbose=False)


    # data from policy comparison experiment
    n_rollouts = np.array([50, 110, 30])
    vc1_successes = np.array([6, 12, 4])
    rt2_successes = np.array([40, 53, 16])

    alpha = 1-np.sqrt(0.95)

    vc1_ubs = np.array([binom_ci(vc1_successes[i], n_rollouts[i], alpha, 'ub', randomized=True, verbose=False) for i in range(3)])
    rt2_lbs = np.array([binom_ci(rt2_successes[i], n_rollouts[i], alpha, 'lb', randomized=True, verbose=False) for i in range(3)])

    return None
