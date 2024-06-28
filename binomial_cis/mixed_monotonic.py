import numpy as np
import warnings


class Interval:
    """
    Class for storing and performing common operations on intervals.

    Parameters
    ----------
    x_min: float
        Lower bound of interval
    x_max: float
        Upper bound of interval
    """
    def __init__(self, x_min, x_max):
        self.x_min = x_min
        self.x_max = x_max
    
    def get_midpoint(self):
        """
        Computes the midpoint of the interval.

        Returns
        -------
        float
            Midpoint of interval
        """
        return self.x_min + (self.x_max - self.x_min)/2
    
    def volume(self):
        """
        Computes the width of the interval.

        Returns
        -------
        float
            Width of interval
        """
        return self.x_max - self.x_min
    
    def get_feasible(self, F):
        """
        Finds a feasible point in the interval.

        Parameters
        ----------
        F: function
            Mixed monotonic function (increasing in first arg, decreasing in second arg)

        Returns
        -------
        x_feas: float
            Midpoint of interval
        x_feas_val: float
            Value of midpoint on the function `F`
        """
        x_feas = self.get_midpoint()
        x_feas_val = F(x_feas, x_feas)
        return x_feas, x_feas_val
    
    def get_ub(self, F):
        """
        Computes a certified upper bound of the function `F` over the interval.

        Parameters
        ----------
        F: function
            Mixed monotonic function (increasing in first arg, decreasing in second arg)
        
        Returns
        -------
        ub: float
            Certified upper bound of max(F) over Interval
        """
        ub = F(self.x_max, self.x_min)
        return ub
    
    def split(self):
        """
        Split along longest axis, else split along first axis.

        Returns
        -------
        I1, I2: Interval
            Intervals resulting from splitting the orignal in half.
        """
        midpoint = self.get_midpoint()
        I1 = Interval(self.x_min, midpoint)
        I2 = Interval(midpoint, self.x_max)
        return I1, I2


def mmp_solve(F, I, tol=1e-3, max_iters=1000, verbose=True):
    """
    Uses mixed-monotonic programming to solve for a certified maximum of the function F over the interval I.

    Parameters
    ----------
    F: function
            Mixed monotonic function (increasing in first arg, decreasing in second arg)
    I: Interval
        Domain to maximize over.
    
    Returns
    -------
    ub: float
        Certified upper bound of max `F` over `I` within `tol` of global max.
    lb: float
        Certified lower bound of max `F` over `I` within `tol` of global max.
    x_lb: float
        Input that achieves `lb`.
    num_iters: int
        Number of iterations taken for the solve.
    """
    # initialize set of intervals and associated ubs in each
    intervals = [I]
    ubs = np.array([I.get_ub(F)])

    # compute initial feasible point and value
    x_feas, x_feas_val = I.get_feasible(F)

    # branch and bound
    for i in range(max_iters):
        # split interval with highest upper bound
        idx_highest = np.argmax(ubs)
        I_popped = intervals[idx_highest]
        I1, I2 = I_popped.split()

        # remove old interval from set
        del intervals[idx_highest]
        ubs = np.delete(ubs, idx_highest)

        # get feasible values of new intervals
        x_feas_I1, x_feas_val_I1 = I1.get_feasible(F)
        x_feas_I2, x_feas_val_I2 = I2.get_feasible(F)

        # update best feasible value
        if x_feas_val_I1 > x_feas_val:
            x_feas = x_feas_I1
            x_feas_val = x_feas_val_I1
        
        if x_feas_val_I2 > x_feas_val:
            x_feas = x_feas_I2
            x_feas_val = x_feas_val_I2
        
        # only add new intervals if they can plausibly contain the maximum
        ub_I1 = I1.get_ub(F)
        ub_I2 = I2.get_ub(F)

        if ub_I1 >= x_feas_val:
            intervals.append(I1)
            ubs = np.append(ubs, ub_I1)
        
        if ub_I2 >= x_feas_val:
            intervals.append(I2)
            ubs = np.append(ubs, ub_I2)
        
        # sanity checks
        # print("\nI1:", I1.x_min, I1.x_max)
        # print("ub:", ub_I1)
        # print("x_feas:", x_feas_I1, x_feas_val_I1)
        assert len(intervals) == len(ubs)
        assert ub_I1 >= x_feas_val_I1
        assert ub_I2 >= x_feas_val_I2
        assert np.max(ubs) >= x_feas_val

        # prune all intervals with ub < lb
        keep_idxs = np.argwhere(ubs >= x_feas_val)[:,0]
        intervals = [intervals[idx] for idx in keep_idxs]
        ubs = ubs[keep_idxs]
        remaining_volume = sum(I.volume() for I in intervals)
        
        # print status
        if verbose:
            print("\n# of Splits: ", i+1)
            print("Current ub-lb gap: ", np.max(ubs) - x_feas_val)
            print("Remaining volume: ", remaining_volume)
        
        # return if there exists no interval for which ub > x_feas_val + tol
        if np.max(ubs) <= x_feas_val + tol:
            return np.max(ubs), x_feas_val, x_feas, i+1

    warnings.warn('Reached maximum number of iterations without achieving tolerance')
    ub, lb, x_lb = np.max(ubs), x_feas_val, x_feas
    return ub, lb, x_lb, max_iters