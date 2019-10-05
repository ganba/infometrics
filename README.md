
<!-- README.md is generated from README.Rmd. Please edit that file -->

# infometrics

<!-- badges: start -->

<!-- badges: end -->

The goal of the infometrics package is to provide tools for estimation,
statistical inference of Linear Inverse Models with Noise, General
Linear Models with Noise, Matrix Balancing (Markov) Models, Multinomial,
Conditional and Mixed Choice Discrete Choice Models using the
infometrics framework.

## Installation

You can install the released version of infometrics from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("infometrics")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ganba/infometrics")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(infometrics)
## basic example code
```

An example of the Matrix Balancing Problem.

``` r
y <- c(0.3, 0.2, 0.5)
x <- c(0.1, 0.25, 0.3, 0.2, 0.15)
gce_matrix(y, x)
#> $lambda
#> [1] -0.02861574 -0.11657431  0.14582581
#> 
#> $p
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 0.3313380 0.3280720 0.3269115 0.3291968 0.3302855
#> [2,] 0.3255602 0.3139563 0.3101061 0.3178159 0.3216841
#> [3,] 0.3431018 0.3579717 0.3629824 0.3529874 0.3480304
#> 
#> $w
#>           [,1]      [,2]      [,3]      [,4]      [,5]
#> [1,] 0.2116069 0.2056374 0.1998363 0.1941989 0.1887205
#> [2,] 0.2491154 0.2217037 0.1973083 0.1755973 0.1562753
#> [3,] 0.1462764 0.1692410 0.1958109 0.2265522 0.2621196
#> 
#> $e
#> [1] -0.02860559 -0.11589334  0.14449885
#> 
#> $Sp
#> [1] 0.9990139
#> 
#> $S_p_j
#> [1] 0.9997824 0.9986300 0.9980228 0.9991252 0.9995091
#> 
#> $p_e_j
#> [1] -0.09882992 -0.09998231 -0.10058947 -0.09948708 -0.09910318
#> 
#> $H_p_w
#> [1] 0.02034069
#> 
#> $ER
#> [1] 0.01083344
#> 
#> $Pseudo_R2
#> [1] 0.0009861024
#> 
#> $conv
#> [1] 0
#> 
#> $hess
#>             [,1]       [,2]        [,3]
#> [1,]  1.09814930 -0.0464314 -0.05278199
#> [2,] -0.04643140  1.0795569 -0.05055660
#> [3,] -0.05278199 -0.0505566  1.07626690
```
