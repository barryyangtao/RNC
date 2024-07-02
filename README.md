# Random Number Converter (RNC) to Universal Distributions
A Simple Algorithm for Converting Random Number Generator Outputs to Universal Distributions to Aid Teaching and Research in Modern Physical Chemistry

- Contributors: Barry Y. Li, Tim Duong, Daniel Neuhauser, Anastassia N. Alexandrova, and Justin R. Caram
  
  (Department of Chemistry and Biochemistry, UCLA)

Abstract
Molecules and materials often display distributed heterogeneous properties.  Modeling the behavior of a single molecule with these properties (e.g. via Monte-Carlo simulations) requires a rapid way to sample distributions of interest and describe how they evolve. Here, we pedagogically introduce a simple algorithm that converts evenly distributed random number generator outputs to any defined probability density function (PDF) called the Random Number Converter (RNC). We demonstrate a numerical approach to obtain cumulative distribution functions (CDFs) and utilize a binary search algorithm to circumvent the need for analytical inverse CDFs. This simple method is demonstrated for various distributions, including single-exponential and Gaussian distributions and non-standard PDFs, for which neither the CDF nor its inverse are analytically solvable. We then apply this algorithm to the rate analysis of catalytic turnover cycles and Fourier spectroscopy, enabling the study of complex reaction kinetics and retrospective interferometry analysis. Our algorithm and examples can be used to train undergraduates on the tools that underlay Monte Carlo methods and connect physical chemistry to standard computer science and statistics. We provide a MATLAB and Python module in the Supporting Information with customizable parameters.

- The usage of the codes is detailed in the file: rnc_si_v6d0.pdf
- A preprint of the relevant article is in: https://doi.org/10.26434/chemrxiv-2024-4thlz
