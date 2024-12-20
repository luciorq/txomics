---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# txomics

<!-- badges: start -->
<!-- badges: end -->

The goal of txomics is to ...

## Installation

You can install the development version of txomics from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("luciorq/txomics")
```

## Example

This is a basic example which shows you how to solve a common problem:

```{r example}
library(txomics)
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so:

```{r cars}
summary(cars)
```

You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this.

You can also embed plots, for example:

```{r pressure, echo = FALSE}
plot(pressure)
```

In that case, don't forget to commit and push the resulting figure files, so they display on GitHub and CRAN.
