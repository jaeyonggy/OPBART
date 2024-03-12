# Ordered probit Bayesian additive regression trees


## OPBART and semi-OPBART

This is the implementation of OPBART and semi-OPBART from the paper: 

**Lee, J., & Hwang, B. S. (2024). Ordered probit Bayesian additive regression trees for ordinal data. Stat, 13(1), e643. https://doi.org/10.1002/sta4.643**

Our work is heavily based on the R package 'SoftBart' by Linero and Yang (2018) as we used their BART implementation as the basis for our models.

OPBART is a very capable alternative to ordered probit regression for ordinal data. Semi-OPBART is a semiparametric version of OPBART which offers enhanced interpretability. Semiparametric BART follows the general BART framework proposed by Tan and Roy (2019). We also implemented Bayesian ordered probit regression from Albert and Chib (1993) to compare as a baseline model. The implementations are available in the R folder.

For a detailed example, check out demonstration.Rmd in the Demo folder or the link below.

[Demonstration](https://rawcdn.githack.com/jaeyonggy/OPBART/main/Demo/demonstration.nb.html)


## Reference

- Albert, J. H., and Chib, S. (1993). Bayesian analysis of binary and polychotomous response data. Journal of the American statistical Association, 88(422), 669-679.
- Linero, A. R., and Yang, Y. (2018). Bayesian regression tree ensembles that adapt to smoothness and sparsity. Journal of the Royal Statistical Society. Series B (Statistical Methodology), 80(5), 1087-1110.
- Tan, Y. V., and Roy, J. (2019). Bayesian additive regression trees and the General BART model. Statistics in medicine, 38(25), 5048-5069.


