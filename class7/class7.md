Class 7: Functions and Packages
================
Nasha Luevit
January 29, 2019

Functions Revisit
-----------------

``` r
source("http://tinyurl.com/rescale-R")
```

Try the **rescale()** function

``` r
rescale( c(1, 5, 10) )
```

    ## [1] 0.0000000 0.4444444 1.0000000

Try **rescale2()** with the stop() function catch for non-numeric input

``` r
#rescale2( c(1:5, "string"))
#doesn't work
```

Functions practice

``` r
x <- c(3, 7, NA, 4, 8, NA)
which( is.na(x) )
```

    ## [1] 3 6

Define an example x and y
=========================

``` r
x <- c( 1, 2, NA, 3, NA)
y <- c(NA, 3, NA, 3, 4)
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
sum( is.na(x) )
```

    ## [1] 2

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

``` r
sum( is.na(y) )
```

    ## [1] 2

Put it together
===============

``` r
is.na(x)
```

    ## [1] FALSE FALSE  TRUE FALSE  TRUE

``` r
is.na(y)
```

    ## [1]  TRUE FALSE  TRUE FALSE FALSE

``` r
is.na(x) & is.na(y)
```

    ## [1] FALSE FALSE  TRUE FALSE FALSE

``` r
sum( is.na(x) & is.na(y))
```

    ## [1] 1

Take working snippet and make first function

``` r
both_na <- function(x, y) {
  sum( is.na(x) & is.na(y) ) 
}
```

``` r
both_na(x,y)
```

    ## [1] 1

``` r
x <- c(NA, NA, NA)
y1 <- c( 1, NA, NA)
y2 <- c( 1, NA, NA, NA)
y3 <- c(1, NA, NA, NA, NA)
```

``` r
both_na(x, y2)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 3

``` r
both_na(x, y3)
```

    ## Warning in is.na(x) & is.na(y): longer object length is not a multiple of
    ## shorter object length

    ## [1] 4

``` r
both_na2 <- function(x, y) {
 if(length(x) != length(y)) {
 stop("Input x and y should be the same length")
 }
 sum( is.na(x) & is.na(y) )
}
```

``` r
#both_na2(x, y2)
#error
```

``` r
both_na3 <- function(x, y) {
 if(length(x) != length(y)) {
 stop("Input x and y should be vectors of the same length")
 }

 na.in.both <- ( is.na(x) & is.na(y) )
 na.number <- sum(na.in.both)
 na.which <- which(na.in.both)
 message("Found ", na.number, " NA's at position(s):",
 paste(na.which, collapse=", ") )

 return( list(number=na.number, which=na.which) )
}
```

``` r
both_na3(x, y1)
```

    ## Found 2 NA's at position(s):2, 3

    ## $number
    ## [1] 2
    ## 
    ## $which
    ## [1] 2 3

Install package
===============

``` r
#install.packages("shiny")
```

``` r
library(shiny)
```

``` r
?shiny
```
