# Ordered probit (semiparametric) Bayesian additive regression trees

## OPBART, OPSMBART

This is the implementation of OPBART and OPSMBART from the paper

- OPBART (2023)

The implementation is heavily based on the R package “SoftBart” by Linero(2018) as we used their BART implementation as a basis for our model.

OPBART is a very capable alternative to ordered probit regression for ordinal data. OPSMBART is a semiparametric version of OPBART which has an added interpretability. We also implemented Bayesian ordered probit regression from Albert(1993) to compare as a baseline model.

For detailed example, check out demonstration.Rmd in the Demo folder or the link below.

[Demonstration](https://rawcdn.githack.com/jaeyonggy/OPBART/main/Demo/demonstration.nb.html)


## Warning

The calculation of marginal effects for dummy variables by OPSMBART is not implemented yet.

## Future work

- Implement marginal effects for dummy variables by OPSMBART.
- Optimize the codes in general.




