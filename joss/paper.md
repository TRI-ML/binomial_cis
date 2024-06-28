---
title: 'binomial_cis: A Python Package for Optimal Binomial Confidence Intervals'
tags:
  - Python
  - statistics
  - confidence intervals
  - binomial
authors:
  - name: Joseph A. Vincent
    orcid: 0000-0002-2270-7395
    affiliation: 1
affiliations:
 - name: Department of Aeronautics and Astronautics, Stanford University
   index: 1
date: 28 June 2024
bibliography: paper.bib
---


# Summary
[binomial_cis](https://github.com/TRI-ML/binomial_cis) is a Python package for computing confidence intervals for the probability of success parameter, $p$, of a binomial distribution. The binomial distribution represents the likelihood of observing $k$ successes in $n$ trials where the probability of success for each trial is $p$. For example, $p$ may be the probability of a coin flip landing on heads, and $k$ the number of heads we observe after $n$ flips. One often does not know the value of $p$ and wishes to estimate it. A confidence interval is a set, constructed based on $k, n$, that covers the unknown parameter $p$ with some user-specified probability. The binomial_cis package computes confidence intervals that lower and/or upper bound $p$ with a user-specified probability. 



# Statement of Need

Constructing confidence intervals for an unknown probability success given samples of successes and failures is one of the most fundamental problems in statistical inference.
Confidence intervals for success probabilities are used in many disciplines including medicine[@altman2013statistics], astronomy[@cameron2011estimation], and robotics[@vincent2024generalizable].
Research into this question dates back at least to the 1930s with the work of Clopper and Pearson [@clopper_pearson].
A foundational result for constructing binomial confidence intervals of minimal width was given by [@eudey1949] and is formalized in [@lehmann_textbook].
We refer to these intervals as *optimal binomial confidence intervals* and they have the property of being uniformly most accurate (UMA) and uniformly most accurate unbiased (UMAU).
Practically, these intervals can provide better inference of $p$ at small sample sizes ($n \le 50$).
The binomial_cis package is the first open-source implementation of these optimal binomial confidence intervals.
In addition, this package provides worst-case analysis of the tightness for the confidence intervals, a feature that is not present in other software for binomial confidence intervals.
Practically, this feature assists the user in understanding how many samples an experiment should have in order to meet a desired level of accuracy in inferring the value of $p$.


# Comparison to Existing Software
There are many existing software packages for computing binomial confidence intervals.
binomial_cis differs from the existing software by providing:

1. Open-source implementations for the optimal binomial confidence intervals given by [@eudey1949] and formalized in [@lehmann_textbook].

2. Functionality for worst-case analysis of the tightness for the confidence intervals, which helps guide users on selecting the sample size for their experiments.


# Research Usage

binomial_cis has been used to compute confidence intervals for the success rate of robots in simulated and real-world tasks [@vincent2024generalizable].



# Acknowledgements

Financial support was provided by Toyota Research Institute, where the author began development of the software during an internship.



# References
