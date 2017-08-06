###########################################################################
##
##   Copyright 2013, 2014 Charles Danko and Minho Chae.
##
##   This program is part of the groHMM R package
##
##   groHMM is free software: you can redistribute it and/or modify it
##   under the terms of the GNU General Public License as published by
##   the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   This program is distributed in the hope that it will be useful, but
##   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
##   or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
##   for more details.
##
##   You should have received a copy of the GNU General Public License along
##   with this program.  If not, see <http://www.gnu.org/licenses/>.
##
##########################################################################

########################################################################
##
##  Distribution Fitting Functions -- MLEfit given X.
##  Date: 2014-03-27
##
##  Presently, fits:
##      -- Gamma distribution  (RgammaMLE)
##      -- Normal distribution (Rnorm)
##
##  Todo:
##  (1) Add fits for:
##      -- Normal by MLE
##      -- Poisson
##      -- Negative binomial (??)
##
########################################################################


#' RgammaMLE fits a gamma distribution to a specified data vector using maximum
#' likelihood.
#'
#' @param X A vector of observations, assumed to be real numbers in the
#' interval [0,+Inf).
#' @return Returns a list of parameters for the best-fit gamma
#' distribution (shape and scale).
#' @author Charles G. Danko
RgammaMLE <- function(X) {
    if(sum(X<0) > 0) message("Negative values not allowed!")
    N <- as.double(NROW(X))
    sumxis <- as.double(sum(X))
    sumlogxis <- as.double(sum(log(X)))
    Fit <- .Call("Rgamma", N, sumxis, sumlogxis, PACKAGE = "groHMM")
    return(Fit)
}

#' Rnorm fits a normal distribution to a specified data vector using maximum
#' likelihood.
#'
#' @param X A vector of observations, assumed to be real numbers in the
#' interval (-Inf,+Inf).
#' @return Returns a list of parameters for the best-fit normal distribution
#' (mean and variance).
#' @author Charles G. Danko
Rnorm <- function(X) {
    returnList      <- list()
    returnList$mean <- mean(X)
    returnList$var  <- var(X)

    return(returnList)
}

################################
##
##  R interface to MLEFit for alpha*Normal+(1-alpha)*Exponential hybrid
##  distribution.
##
################################

#' Rnorm.exp fits a normal+exponential distribution to a specified data
#' vector using maximum likelihood.
#'
#' Distrubution function defined by: alpha*Normal(mean, variance)+(1-alpha)
#' *Exponential(lambda).
#'
#' Fits nicely with data types that look normal overall, but have a long
#' tail starting for positive values.
#'
#' @param xi A vector of observations, assumed to be real numbers in the
#' interval (-Inf,+Inf).
#' @param wi A vector of weights.  Default: vector of repeating 1; indicating
#' all observations are weighted equally. (Are these normalized internally?!
#' Or do they have to be [0,1]?)
#' @param guess Initial guess for parameters.  Default: c(0.5, 0, 1, 1).
#' @param tol Convergence tolerance.  Default: sqrt(.Machine$double.eps).
#' @param maxit  Maximum number of iterations.  Default: 10,000.
#' @return Returns a list of parameters for the best-fit normal distribution
#' (alpha, mean, variance, and lambda).
#' @author Charles G. Danko
Rnorm.exp <- function(xi, wi=rep(1,NROW(xi)), guess=c(0.5, 0, 1, 1),
    tol=sqrt(.Machine$double.eps), maxit=10000) {
    Fit <- .Call("RNormExpMLE", xi, wi, guess, tol, as.integer(maxit),
        PACKAGE = "groHMM")
}
