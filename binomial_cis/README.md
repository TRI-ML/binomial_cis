# File Organization
Below is a list of each source code file and its purpose.

- `binomial_helper.py`: This file defines helper functions for the binomial distribution. Specifically, computation of the binomial coefficient, binomial pmf, and binomial cdf are implemented. This implementation is used instead of SciPy. This is because without code compilation computing expected shortage (and expected excess and expected width) is too slow, and to compile the expected shortage (excess, width) function we have to implement and compile these binomial helper functions ourselves.

- `conf_intervals.py`: This file contains the main source code for the package. The most important function defined in this file is the `binom_ci()` function which is the primary function the user of this package will call.

- `mixed_monotonic.py`: This file defines the `mmp_solve()` function which implements a branch and bound procedure to solve mixed monotonic programs. For background on mixed monotonic programming refer to [*Mixed Monotonic Programming for Fast Global Optimization*](https://arxiv.org/pdf/1910.07853.pdf) by Matthiesen and Jorswieck. In our package, the computation of maximum expected shortage (excess, width) is framed as a mixed monotonic program.

- `plotting.py`: This file defines functions for visualizing the mixed-monotonic property of expected shortage (width) as well as visualizing expected shortage (width) as a function of the true probability of success.

- `volume.py`: This file contains functions for computing expected shortage (excess, width) and maximum expected shortage (excess, width).