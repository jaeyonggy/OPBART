# Ordered probit (semiparametric) Bayesian additive regression trees

## OPBART, OPSMBART

This is the implementation of OPBART and OPSMBART from the paper: OPBART (2023).

Our work is heavily based on the R package 'SoftBart' by Linero and Yang (2018) as we used their BART implementation as a basis for our models.

OPBART is a very capable alternative to ordered probit regression for ordinal data. OPSMBART is a semiparametric version of OPBART which offers enhanced interpretability. We also implemented Bayesian ordered probit regression from Albert and Chib (1993) to compare as a baseline model. The implementations are available in the R folder.

For a detailed example, check out demonstration.Rmd in the Demo folder or the link below.

[Demonstration](https://rawcdn.githack.com/jaeyonggy/OPBART/main/Demo/demonstration.nb.html)


## Warning

The calculation of marginal effects for dummy variables by OPSMBART is not implemented yet.

## Future work

- Implement marginal effects for dummy variables by OPSMBART.
- Optimize the codes in general.

## Reference

- Albert, J. H., and Chib, S. (1993). Bayesian analysis of binary and polychotomous response data. Journal of the American statistical Association, 88(422), 669-679.
- Linero, A. R., and Yang, Y. (2018). Bayesian regression tree ensembles that adapt to smoothness and sparsity. Journal of the Royal Statistical Society. Series B (Statistical Methodology), 80(5), 1087-1110.


