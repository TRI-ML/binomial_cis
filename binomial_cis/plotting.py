import matplotlib.pyplot as plt
import numpy as np
from binomial_cis.conf_intervals import llc_accept_prob, llc_accept_prob_2_sided
from binomial_cis.volume import expected_shortage, expected_shortage_mixed_monotonic, expected_shortage_cp, expected_shortage_mixed_monotonic_cp, max_expected_shortage
from binomial_cis.volume import expected_width, expected_width_mixed_monotonic, max_expected_width

def plot_expected_shortage_mixed_monotonic(alpha, n, p1s, p2s, verbose=True, randomized=True):
    """
    Contour plot of the mixed-monotonic form of expected shortage.

    Inputs
    alpha: miscoverage rate
    n: number of samples
    p1s: prob. of success parameter for limit of integration
    p2s: prob. of success parameter for CDF
    randomized: if False then use Clopper-Pearson
    
    Returns
    Contour plot showing expected shortage for each parameter combination.
    """
    # given grid of p1s and p2s
    p1s2D, p2s2D = np.meshgrid(p1s, p2s, indexing="ij")
    ns = np.zeros_like(p1s2D)
    print("ns.shape: ", ns.shape)

    # solve for expected shortage
    for i,p1 in enumerate(p1s):
        print("iter: ", i, "/", len(p1s)) if verbose else None
        for j,p2 in enumerate(p2s):
            if randomized:
                ns[i,j] = expected_shortage_mixed_monotonic(llc_accept_prob, alpha, n, p1, p2)
            else:
                ns[i,j] = expected_shortage_mixed_monotonic_cp(alpha, n, p1, p2)

    # plot results
    fig, ax = plt.subplots()
    c = ax.contourf(p1s2D, p2s2D, ns, cmap='RdBu', levels=20)
    ax.set_title('pcolormesh')
    ax.axis([p1s.min(), p1s.max(), p2s.min(), p2s.max()])
    fig.colorbar(c, ax=ax)
    plt.xlabel(r"$p_1$")
    plt.ylabel(r"$p_2$")
    plt.title( r"Expected Shortage as Mixed Monotonic Function" )


def plot_shortage_curve(alpha, n, num_p=21, verbose=True, randomized=True):
    """
    Plot expected shortage as a function of true prob of success p.

    Inputs
    alpha: miscoverage rate
    n: number of samples
    num_p: number of values of p to compute shortage for
    randomized: if False then use Clopper-Pearson

    Returns
    Plot of curve which visualizes expected shortage as a function of true prob of success p.
    """
    ps = np.linspace(0,1,num=num_p)
    exp_shortages = np.zeros(num_p)
    for i,p in enumerate(ps):
        print("p: ", p) if verbose else None
        if randomized:
            exp_shortages[i] = expected_shortage(llc_accept_prob, alpha, n, p)
        else:
            exp_shortages[i] = expected_shortage_cp(alpha, n, p)
    print("Finished computing grid of expected shortages.")

    print("\nBegginning computation of max expected shortage.")
    max_es, lb, p_lb, num_iters = max_expected_shortage(alpha, n, tol=1e-3, verbose=verbose, randomized=randomized)

    print("\nGlobal Maximum E[shortage]:  ", max_es)
    print("Sample Maximum E[shortage]:  ", max(exp_shortages))
    
    fig = plt.figure(figsize=(5,5))
    plt.plot(ps, exp_shortages, color="cornflowerblue", linewidth=3)
    plt.axhline(y = max_es, color = 'orchid', linestyle = '--', linewidth=3, label="Max E[shortage]")
    title = str(r"Visualizing Expected Shortage: $\alpha=$" + str(alpha) + r", $n=$" + str(n))
    plt.title(title)
    plt.xlabel(r'True Probability of Success, $p$')
    plt.ylabel(r'Expected Shortage')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.legend()









def plot_expected_width_mixed_monotonic(alpha, n, p1s, p2s, verbose=True):
    """
    Contour plot of the mixed-monotonic form of expected width.

    Inputs
    alpha: miscoverage rate
    n: number of samples
    p1s: prob. of success parameter for first CDF tern
    p2s: prob. of success parameter for second CDF term
    
    Returns
    Contour plot showing expected width for each parameter combination.
    """
    # given grid of taus and alphas
    p1s2D, p2s2D = np.meshgrid(p1s, p2s, indexing="ij")
    ns = np.zeros_like(p1s2D)
    print("ns.shape: ", ns.shape)

    # solve for ns
    for i,p1 in enumerate(p1s):
        print("iter: ", i, "/", len(p1s)) if verbose else None
        for j,p2 in enumerate(p2s):
            ns[i,j] = expected_width_mixed_monotonic(llc_accept_prob_2_sided, alpha, n, p1, p2)

    # plot results
    fig, ax = plt.subplots()
    c = ax.contourf(p1s2D, p2s2D, ns, cmap='RdBu', levels=20)
    ax.set_title('pcolormesh')
    ax.axis([p1s.min(), p1s.max(), p2s.min(), p2s.max()])
    fig.colorbar(c, ax=ax)
    plt.xlabel(r"$p_1$")
    plt.ylabel(r"$p_2$")
    plt.title( r"Expected Width as Mixed Monotonic Function" )


def plot_width_curve(alpha, n, num_p=21, verbose=True):
    """
    Plot expected width as a function of true prob of success p.

    Inputs
    alpha: miscoverage rate
    n: number of samples
    num_p: number of values of p to compute width for

    Returns
    Plot of curve which visualizes expected width as a function of true prob of success p.
    """
    ps = np.linspace(0,1,num=num_p)
    exp_widths = np.zeros(num_p)
    for i,p in enumerate(ps):
        print("p: ", p) if verbose else None
        exp_widths[i] = expected_width(llc_accept_prob_2_sided, alpha, n, p)
    print("Finished computing grid of expected width.")

    print("\nBegginning computation of max expected width.")
    max_ew, lb, p_lb, num_iters = max_expected_width(alpha, n, tol=1e-3, verbose=verbose)

    print("\nGlobal Maximum E[width]:  ", max_ew)
    print("Sample Maximum E[width]:  ", max(exp_widths))
    
    fig = plt.figure(figsize=(5,5))
    plt.plot(ps, exp_widths, color="cornflowerblue", linewidth=3)
    plt.axhline(y = max_ew, color = 'orchid', linestyle = '--', linewidth=3, label="Max E[width]")
    title = str(r"Visualizing Expected Width: $\alpha=$" + str(alpha) + r", $n=$" + str(n))
    plt.title(title)
    plt.xlabel(r'True Probability of Success, $p$')
    plt.ylabel(r'Expected Width')
    plt.xlim(0, 1)
    plt.ylim(0, 1)
    plt.legend()