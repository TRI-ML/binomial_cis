import numpy as np
from numba import cfunc
import numba as nb
from scipy import LowLevelCallable
from binomial_cis.binomial_helper import binom_pmf, binom_cdf

# CDF() acts differently at p=0 and p=1, so we use epsilon to avoid this.
epsilon = 1e-10 # tolerance


def test_stat(num_successes, randomized=True):
    """
    Inputs
    num_successes: number of observed successes in n trials
    
    Returns
    T: number of successes + uniform random variable
    """
    T = num_successes + np.random.rand()
    return T


@cfunc("float64(float64, float64, float64)")
def CDF(t, p, n):
    """
    Inputs
    t: cutoff value of test stat distribution
    p: probability of success
    n: number of samples
    
    Returns
    cdf: integral of pdf from t=0 to t=t
    """
    t_ = np.floor(t)

    # since pdf is piecewise constant where each piece is width 1,
    # we can compute the CDF using the CDF of the binomial plus
    # the leftover area we need to integrate the pdf
    if t_ == 0:
        cdf = (t - t_) * binom_pmf(t_, n, p)
    else:
        cdf = binom_cdf(t_ - 1, n, p) + (t - t_) * binom_pmf(t_, n, p)
    return cdf


@cfunc("float64(intc, CPointer(float64))")
def accept_prob(num_args, args):
    """
    Inputs
    num_args: required for code compilation. https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html
    args: list containing the following
        p_0: hypothesis probability of success (or x_1 for mixed_monotonic version)
        alpha: miscoverage rate, P(p_ <= p) >= 1-alpha
        n: number of samples
        p: true probability of success (or x_2 for mixed_monotonic version)

    Returns
    accept_prob: probability of accepting p_0 given the true p
    """
    # unpack args
    p_0, alpha, n, p = args[0], args[1], args[2], args[3]

    # (can't use helper method if we want to compile bc can't compile functions with function inputs)
    tol = 1e-6
    num_iters = int(np.ceil( np.log2(1 / tol) )) + 10
    def F(t): return CDF(t, p_0, n)
    RHS = 1-alpha

    # initial endpoints
    lb = 0
    ub = n + 1
    lb_val = F(lb)
    ub_val = F(ub)

    # perform bisection iterations
    for _ in range(num_iters):
        mid = lb + (ub - lb)/2
        mid_val = F(mid)

        # refine endpoints
        if mid_val >= RHS:
            ub = mid
            ub_val = mid_val
        else:
            lb = mid
            lb_val = mid_val

        # sanity checks
        assert(lb < ub)
        assert(lb_val <= RHS)
        assert(ub_val >= RHS)

        # return t*
        if ub_val - RHS <= tol:
            t_star = ub
            break

    accept_prob = CDF(t_star, p, n)
    return accept_prob

llc_accept_prob = LowLevelCallable(accept_prob.ctypes) # SciPy LowLevelCallable version of accept_prob


def binom_ci(k, n, alpha, side, verbose=True, randomized=True):
    """
    Inputs
    k: number of observed successes in n trials
    n: number of trials
    alpha: miscoverage rate, P(p in CI) = 1-alpha
    verbose: whether to print intermediate updates
    randomized: if true solves for the UMA bounds, if false then solves for (less efficient) non-randomized bounds
    
    Returns
    CI: either a lower bound, upper bound, or simultaneous lower & upper bounds
    """
    if side == "lb":
        print("Comuting lower confidence bound") if verbose else None
        CI = get_lb(k, n, alpha, randomized=randomized)
    elif side == "ub":
        print("Comuting upper confidence bound") if verbose else None
        # lb on failure prob is 1 - ub on success prob
        q_lb = get_lb(n-k, n, alpha, randomized=randomized)
        CI = 1 - q_lb
    elif side == "lb,ub":
        print("Comuting 2-sided confidence bound") if verbose else None
        CI = get_lb_ub(k, n, alpha)
    else:
       raise ValueError("Invalid argument for 'side' given!")
    
    return CI
    

def get_lb(num_successes, n, alpha, randomized=True):
    """
    Inputs
    num_successes: number of observed successes in n trials
    n: number of trials
    alpha: miscoverage rate, P(p_ <= p) >= 1-alpha
    
    Returns
    p_: a lower bound on p such that P(p_ <= p) >= 1-alpha
    """

    if not randomized:
        # compute Clopper-Pearson bound
        p_ = binom_bisection(num_successes-1, n, alpha)

        # could also find p_ by setting U=0:
        # def cdf(p): return CDF(num_successes, p, n)
        # p_ = bisection(cdf, alpha, tol=1e-8)
        return p_

    t_o = test_stat(num_successes)

    # check edge cases
    if t_o < 1-alpha:
        p_ = 0
    elif t_o > n + 1-alpha:
        p_ = 1
    else:
        # typical case
        def cdf(p): return CDF(t_o, p, n)
        p_ = bisection(cdf, alpha, tol=1e-8)
    
    return p_


def bisection(CDF, alpha, tol=1e-6):
    """
    Inputs
    CDF: a function mapping from p -> F(t | p)
    alpha: miscoverage rate, P(p_ <= p) >= 1-alpha
    tol: upper bound on how far below the returned value of p is to the true solution
    
    Returns
    p_: a lower bound on p such that P(p_ <= p) >= 1-alpha
    """
    # solving -CDF(p) = -(1-alpha) for p where CDF is monotonically decreasing in p
    def F(p): return -CDF(p) # now have monotonically increasing
    RHS = -(1-alpha)

    # minimum number of iterations to guarantee the tolerance
    # see https://en.wikipedia.org/wiki/Bisection_method#Analysis
    # add buffer because floating point arithmetic sometimes requiring more iters
    n = int(np.ceil( np.log2(1 / tol) )) + 10

    lb = 0 + epsilon
    ub = 1 - epsilon
    lb_val = F(lb)
    ub_val = F(ub)

    # perform bisection iterations
    for _ in range(n):
        mid = lb + (ub - lb)/2
        mid_val = F(mid)

        # refine endpoints
        if mid_val >= RHS:
            ub = mid
            ub_val = mid_val
        else:
            lb = mid
            lb_val = mid_val

        # sanity checks
        # print("\nlb: ", lb)
        # print("ub: ", ub)
        # print("lb_val: ", lb_val)
        # print("ub_val: ", ub_val)
        # print("RHS: ", RHS)
        assert(lb < ub)
        assert(lb_val <= RHS)
        assert(ub_val >= RHS)

        # return tight lower bound on true p_
        if ub_val - RHS <= tol:
            p_ = ub
            return p_

    raise ValueError("Exited bisection search with no solution. Error in code.")




















#########################################
##### Functions for Clopper-Pearson #####
#########################################

@cfunc("float64(intc, CPointer(float64))")
def accept_prob_cp(num_args, args):
    """
    Inputs
    num_args: required for code compilation. https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html
    args: list containing the following
        p_0: hypothesis probability of success (or x_1 for mixed_monotonic version)
        alpha: miscoverage rate, P(p_ <= p) >= 1-alpha
        n: number of samples
        p: true probability of success (or x_2 for mixed_monotonic version)

    Returns
    accept_prob: probability of accepting p_0 given the true p
    """
    # unpack args
    # probability that cp lower bound is less than p_0, given our observed data is generated from Bernoulli(p)
    p_0, alpha, n, p = args[0], args[1], args[2], args[3]

    # t_star(p_0) = max{t \in [0, 1, ..., n] | binom_cdf(t,n,p_0) <= 1-alpha}
    t_star = 0
    for i in range(n+1):
        if binom_cdf(i,n,p_0) > 1-alpha:
            t_star = i-1
            break
    
    accept_prob = binom_cdf(min(t_star+1, n), n, p)
    
    return accept_prob

llc_accept_prob_cp = LowLevelCallable(accept_prob_cp.ctypes) # SciPy LowLevelCallable version of accept_prob



def binom_bisection(k, n, alpha):
    """
    Inputs
    k: number of successes
    n: number of samples
    alpha: miscoverage rate, P(p_ <= p) >= 1-alpha

    Returns
    p: smallest p such that binom_cdf(k,n,p) >= 1-alpha
    """
    tol = 1e-5
    num_iters = 100 # int(np.ceil( np.log2(1 / tol) )) + 10

    # binom_cdf is decreasing as p increases
    def F(p): return -binom_cdf(k, n, p) # negative so now monotonically increasing in p
    RHS = -(1-alpha)

    # initial endpoints
    lb = epsilon
    ub = 1 - epsilon
    lb_val = round(F(lb), 6)
    ub_val = round(F(ub), 6)

    if lb_val > RHS:
        return 0
    elif ub_val < RHS:
        return 1

    # perform bisection iterations
    for i in range(num_iters):
        mid = (lb + ub) / 2.
        mid_val = round(F(mid), 6)
        # print("\nlb:", lb, "  lb_val:", lb_val)
        # print("mid:", mid, "  mid_val:", mid_val)
        # print("ub:", ub, "  ub_val:", ub_val)
        assert(lb_val <= mid_val and mid_val <= ub_val)

        # refine endpoints
        if mid_val >= RHS:
            ub = mid
            ub_val = mid_val
        else:
            lb = mid
            lb_val = mid_val

        # sanity checks
        assert(lb < ub)
        assert(lb_val <= RHS)
        assert(ub_val >= RHS)

        # return p
        if ub_val - RHS <= tol:
            p = mid
            return p

        assert(i < num_iters-1) # should always return before hitting end of for loop




def get_ps_cp(p, n, alpha):
    """
    Inputs
    p: true success probability
    n: number of samples
    alpha: miscoverage rate, P(p_ <= p) >= 1-alpha

    Returns
    ps: array of possible lower bounds from CP that are cutoff points in the expected shortage integral
    """
    ps = np.array([0.])

    for i in range(1,n+1): # i = 1,2,...,n
        p_i = binom_bisection(i-1, n, alpha)
        ps = np.append(ps, [p_i])

        if p_i >= p:
            return ps
        
    # make len(ps) = n+2 if ps[n] < p
    ps = np.append(ps, [1.])

    return ps




















#########################################
###### Functions for 2-sided Bound ######
#########################################

def get_lb_ub(num_successes, n, alpha):
    """
    Inputs
    num_successes: number of observed successes in n trials
    n: number of trials
    alpha: miscoverage rate, P(p_lb <= p <= p_ub) >= 1-alpha
    
    Returns
    p_lb, p_ub: A UMAU confidence interval for p

    We know that the bounds of the test statistic acceptance region
    are monotonically increasing in p (see Lehmann textbook).
    Strategy: bisection search to find bounds on p
    
    """
    t_o = test_stat(num_successes)

    # check edge cases
    # see "Nonoptimality of Randomized Confidence Sets" by Casella
    if t_o < alpha:
        p_lb, p_ub = 0, 0
    elif t_o > n + 1-alpha:
        p_lb, p_ub = 1, 1
    else:
        # typical case
        p_lb = UMAU_lb(t_o, n, alpha)
        p_ub = UMAU_ub(t_o, n, alpha)
    
    return p_lb, p_ub



def UMAU_lb(t_o, n, alpha, tol=1e-6):
    """
    Inputs
    t_o: observed test statistic
    n: number of trials
    alpha: miscoverage rate
    
    Returns
    p_lb: lower bound of UMAU confidence interval for p
    """
    # find smallest value of p, such that t_o in [t_lo, t_up]

    # minimum number of iterations to guarantee the tolerance
    # see https://en.wikipedia.org/wiki/Bisection_method#Analysis
    # add buffer because floating point arithmetic sometimes requiring more iters
    n_iter = int(np.ceil( np.log2(1 / tol) )) + 10

    # first find some value in the CI
    p_feas = get_feasible_p(t_o, n, alpha)
    # print("Found feasible p!")

    # set endpoints
    p_lb_up = p_feas
    p_lb_lo = epsilon
    
    # sanity checks
    assert in_UMPU_accept_region(t_o, p_lb_up, n, alpha)

    # sometimes the lower bound is zero
    if p_feas <= p_lb_lo or in_UMPU_accept_region(t_o, p_lb_lo, n, alpha):
        return 0.0
    
    
    # return if within tolerance
    if p_lb_up - p_lb_lo <= tol:
        return p_lb_lo

    # bisection search
    for _ in range(n_iter):       
        p_lb_mid = (p_lb_up - p_lb_lo)/2 + p_lb_lo
        if in_UMPU_accept_region(t_o, p_lb_mid, n, alpha):
            p_lb_up = p_lb_mid
        else:
            p_lb_lo = p_lb_mid
        
        assert p_lb_up > p_lb_lo

        # return if within tolerance
        if p_lb_up - p_lb_lo <= tol:
            return p_lb_lo
    
    raise ValueError("Bisection search ended with no solution!")
        


def UMAU_ub(t_o, n, alpha, tol=1e-6):
    """
    Inputs
    t_o: observed test statistic
    n: number of trials
    alpha: miscoverage rate
    
    Returns
    p_ub: upper bound of UMAU confidence interval for p
    """
    # find largest value of p, such that t_o in [t_lo, t_up]

    # minimum number of iterations to guarantee the tolerance
    # see https://en.wikipedia.org/wiki/Bisection_method#Analysis
    # add buffer because floating point arithmetic sometimes requiring more iters
    n_iter = int(np.ceil( np.log2(1 / tol) )) + 10

    # first find some value in the CI
    p_feas = get_feasible_p(t_o, n, alpha)

    # set endpoints
    p_lb_up = 1 - epsilon
    p_lb_lo = p_feas
    
    # sanity checks
    assert in_UMPU_accept_region(t_o, p_lb_lo, n, alpha)

    # sometimes the upper bound is 1
    if p_feas >= p_lb_up or in_UMPU_accept_region(t_o, p_lb_up, n, alpha):
        return 1.0
    
    # return if within tolerance
    if p_lb_up - p_lb_lo <= tol:
        return p_lb_up

    # bisection search
    for _ in range(n_iter):       
        p_lb_mid = (p_lb_up - p_lb_lo)/2 + p_lb_lo
        if in_UMPU_accept_region(t_o, p_lb_mid, n, alpha):
            p_lb_lo = p_lb_mid
        else:
            p_lb_up = p_lb_mid
        
        assert p_lb_up > p_lb_lo

        # return if within tolerance
        if p_lb_up - p_lb_lo <= tol:
            return p_lb_up
    
    raise ValueError("Bisection search ended with no solution!")



def get_feasible_p(t_o, n, alpha):
    """
    Inputs
    t_o: observed test statistic
    n: number of trials
    alpha: miscoverage rate
    
    Returns
    p_o: a value of p in the UMAU interval
    """
    # first try the empirical p
    p_o = np.floor(t_o) / n
    if in_UMPU_accept_region(t_o, p_o, n, alpha):
        return p_o
    
    # otherwise do grid search: 5 grid points
    grid_5 = np.linspace(0,1,num=5, endpoint=True)
    for p_o in grid_5:
        if in_UMPU_accept_region(t_o, p_o, n, alpha):
            return p_o
    
    # 10 grid points
    grid_10 = np.linspace(0,1,num=10, endpoint=True)
    grid_10 = list(set(grid_10) - set(grid_5))
    for p_o in grid_10:
        if in_UMPU_accept_region(t_o, p_o, n, alpha):
            return p_o
    
    # 100 grid points
    grid_100 = np.linspace(0,1,num=100, endpoint=True)
    grid_100 = list(set(grid_100) - set(grid_10) - set(grid_5))
    for p_o in grid_100:
        if in_UMPU_accept_region(t_o, p_o, n, alpha):
            return p_o
     
    raise ValueError("Could not find feasible point in 2-sided CI!")



def in_UMPU_accept_region(t_o, p_o, n, alpha):
    """
    Inputs
    t_o: observed test statistic
    p_o: null hypothesis prob of success
    n: number of trials
    alpha: miscoverage (Type 1 error) rate, P(accept | p = p_o) >= 1-alpha
    
    Returns
    boolean for whether t_o is in the UMAU acceptance region for the test statistic
    """
    # print(get_UMPU_accept_region(p_o, n, alpha))
    t_lo, t_up = get_UMPU_accept_region(p_o, n, alpha)
    return t_o >= t_lo and t_o <= t_up



def get_UMPU_accept_region(p_o, n, alpha):
    """
    Inputs
    p_o: null hypothesis prob of success
    n: number of trials
    alpha: miscoverage (Type 1 error) rate, P(accept | p = p_o) >= 1-alpha
    
    Returns
    t_lo, t_up: UMAU acceptance region for the test statistic
    """
    
    # see "Nonoptimality of Randomized Confidence Sets" by Casella
    # or "Table of Neyman-Shortest Unbiased Confidence Intervals for the Binomial Parameter"
    
    # first check edge cases
    if p_o == 0:
        return alpha, 2-alpha
    elif p_o == 1:
        return n-1+alpha, n+1-alpha
    
    n_los = range(n) # {0,1,...,n-1}
    for n_lo in n_los:
        n_ups = range(n_lo+1, n+1)
        for n_up in n_ups:
            gamma_lo, gamma_up = get_gammas(n_lo, n_up, p_o, n, alpha)
            if in_unit(gamma_lo) and in_unit(gamma_up):
                t_lo = n_lo + gamma_lo
                t_up = n_up + gamma_up
                assert t_lo <= t_up
                return t_lo, t_up
    
    raise ValueError("No UMPU acceptance region found!")


# @cfunc("(float64, float64, float64, float64, float64)")
@cfunc(nb.types.UniTuple(nb.float64,2)(nb.float64, nb.float64, nb.float64, nb.float64, nb.float64))
def get_gammas(n_l, n_u, p_o, n, alpha):
    """
    Inputs:
    n_l: lower integer cutoff
    n_u: upper integer cutoff
    p_o: null hyp prob of success
    n: number of trials
    alpha: miscoverage rate

    Returns:
    gamma_l: lower continuous variable cutoff
    gamma_u: upper continuous variable cutoff
    """
    k1 = binom_cdf(n_u-1,n,p_o) - binom_cdf(n_l-1,n,p_o) - (1-alpha)
    k2 = -n_l*(1-p_o)*binom_pmf(n_l,n,p_o) + n_u*(1-p_o)*binom_pmf(n_u,n,p_o)

    num_l = (n_u - n*p_o)*k1 + k2
    den_l = (n_u - n_l)*binom_pmf(n_l, n, p_o)

    num_u = (n_l - n*p_o)*k1 + k2
    den_u = (n_u - n_l)*binom_pmf(n_u, n, p_o)

    # check conditions before doing potentially unsafe division
    gl_cond = num_l <= den_l and num_l >= 0 and den_l > 0
    gu_cond = num_u <= den_u and num_u >= 0 and den_u > 0
    
    if gl_cond and gu_cond:
        gamma_l = num_l / den_l
        gamma_u = num_u / den_u
    else:
        # something outside the unit interval so that we don't accept these gammas
        gamma_l = -10.0
        gamma_u = 10.0

    return gamma_l, gamma_u


@cfunc("boolean(float64)")
def in_unit(x):
    return x >= 0 and x <= 1


@cfunc("float64(intc, CPointer(float64))")
def accept_prob_2_sided(num_args, args):
    """
    Inputs
    num_args: required for code compilation. https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.quad.html
    args: list containing the following
        p_0: hypothesis probability of success (or x_1 for mixed_monotonic version)
        alpha: miscoverage rate, P(p_ <= p) >= 1-alpha
        n: number of samples
        p1: true probability of success for CDF at t_u
        p2: true probability of success for CDF at t_l (only allowed to be differenct for mixed-monotonic form)

    Returns
    accept_prob: probability of accepting p_0 given the true p
    """
    # unpack args
    p_0, alpha, n, p1, p2 = args[0], args[1], args[2], args[3], args[4]

    # first check edge cases 
    if p_0 == 0:
        accept_prob = CDF(2-alpha, p1, n) - CDF(alpha, p2, n)
        return accept_prob
    elif p_0 == 1:
        accept_prob = CDF(n+1-alpha, p1, n) - CDF(n-1+alpha, p2, n)
        return accept_prob
   
    # execute search
    n_los = range(n) # {0,1,...,n-1}
    for n_lo in n_los:
        n_ups = range(n_lo+1, n+1) # {n_lo,...,n}
        for n_up in n_ups:
            gamma_lo, gamma_up = get_gammas(n_lo, n_up, p_0, n, alpha)
            if in_unit(gamma_lo) and in_unit(gamma_up):
                t_lo = n_lo + gamma_lo
                t_up = n_up + gamma_up
                assert t_lo <= t_up
                accept_prob = CDF(t_up, p1, n) - CDF(t_lo, p2, n)
                return accept_prob
    
    assert False
    return None

llc_accept_prob_2_sided = LowLevelCallable(accept_prob_2_sided.ctypes) # SciPy LowLevelCallable version of accept_prob


def max_accept_prob_2_sided(alpha, n, p, n_grid):
    """
    Inputs
    alpha: miscoverage rate, P(p_ <= p) >= 1-alpha
    n: sample size
    p: true probability of success
    n_grid: size of discrete grid to search over

    Returns
    max_ap: maximum acceptance probability
    max_p_0: location of maximum acceptance probability

    It is unclear whether we can rigorously solve for this using mixed monotonic programming
    accept_prob_2_sided is mixed monotonic in p_0, but I'm not sure if we can define the function
    with split p_0 inputs.
    """

    p_0s = np.linspace(0,1,num=n_grid, endpoint=True)
    aps = [accept_prob_2_sided(5, [p_0, alpha, n, p, p]) for p_0 in p_0s]

    argmax = np.argmax(aps)
    max_p_0 = p_0s[argmax]
    max_ap = aps[argmax]

    return max_ap, max_p_0