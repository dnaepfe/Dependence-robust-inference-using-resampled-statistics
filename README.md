This repository contains the generating code for the semester thesis 'Dependence-robust inference
using resampled statistics' by Dominic Näpfer.

Standard inference procedures such as the t-test rely on the assumption of independence. If this assumption fails, e.g. when working with clustered data, these procedures can produce severely misleading results, primarily due to inflation of the Type
I error rate. This thesis studies methods that remain valid under general forms of
weak dependence, following the framework of (Leung, 2021). We introduce two test
statistics: a mean-type and a U-type statistic, both based on random resampling of
observations. We show that under the assumption of √n-consistency, both statistics
converge in distribution to a standard Gaussian. This holds without requiring any
knowledge of the underlying dependence structure and establishes inference procedures for hypothesis testing for any asymptotically linear estimator. We conduct a
Monte Carlo simulation to evaluate the empirical performance of the proposed methods in different scenarios for clustered data, comparing the statistics to the standard
t-test and cluster-robust methods, such as clubSandwich. We conclude with proposing
a practical guidance for applying the methods in empirical research.

