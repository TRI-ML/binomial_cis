The tests are written in the form of interactive notebooks.

- `2_side_validation.ipynb`: This notebook can be used to more easily inspect any differences between our 2-sided bounds and those from the Blyth paper.
- `binom_ci_test.py`: This file can be used to test the correctness of the lower, upper, and 2-sided bounds. 
- `binom_helper_validation.ipynb`: This notebook can be used to validate the binomial helper functions from `binomial_cis/binomial_helper.py`. Accuracy and speed of our implementation is tested against the SciPy implementation.
- `conf_set_validation.ipynb`: This notebook can be used to validate the probabilistic guarantees of the confidence intervals using Monte Carlo simulation.