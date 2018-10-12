# Import required libraries ----

# 1) Vectors ----
#a) Generate 100 rnorm, keep only ones larger than 1
x <- rnorm(100); x[x<1]

#b)
#' @param n integer
#' @param min numerical
#' @description Return n independent random variables from a standard normal distribution truncated below by 0
#' @example 
#' floor_rnorm(4, 0)
floor_rnorm <- function (n, min=0){
  pmax(rnorm(n), min)
}

floor_rnorm <- function (n, min=0){
  x <- rnorm(n)
  x[x>min]
}


# c) Generate 10,000 truncated normals with min at -1, plot histogram with adjusted bins
truncated_normal_rvs <- floor_rnorm(10^4,-1)
hist(truncated_normal_rvs, breaks = seq(-1,5,0.5))

# d) Improve floor_rnorm sampling of rnorm bigger than 4 when min is large
# CDF of normal, with restricted uniform sample, inversion
# Larger variance, and MH
# Exponential and rejection sampling
# Box Muller with restricted uniform
# Take both tails
# MCMC


# 2 Data ----
library(MASS)
data(hills)

# a) What object is 'hills'?
class(hills)
# data.frame

# b) How many cols in 'hills'?
dim(hills)[2]
# 3

# c) Change race 'Two Breweries' to 'Three Breweries'
rownames(hills) <- gsub(x = rownames(hills), pattern = 'Two Breweries', replacement = 'Three Breweries')

# d) Using 'with', find the mean time for races with a climb greater than 1000
with(hills, {mean(time[climb>1000])})

library(nlme)
data(Orthodont)
# e) What object is 'Orthodont'?
class(Orthodont)
# [1] "nfnGroupedData" "nfGroupedData"  "groupedData"    "data.frame"   
#f, g) Print method for grouped data
methods(print)
# print.groupedData*
nlme:::print.groupedData

# 3) Recursion
# a) Fibonacci Recursion
#' @param n integer
#' @description return nth fibonacci number
#' @example 
#' fib(4)
fib <- function(n)
  if(n<=2) { if(n>=0) 1 else 0 } else Recall(n-1) + Recall(n-2)

fib2 <- function(n)
  if(n<=2) { if(n>=0) 1 else 0 } else fib(n-1) + fib(n-2)

#b) Eval 20th fib
fib(20)

# c) Number of evals in 20: 2^19

# d) Fib loop
#' @param n
#' @description return nth fibonacci number
#' @example 
#' fib_loop(4)
fib_loop <- function(n){
  if(n<=2) { if(n>=0) 1 else 0 } else {
    fn_minus_1 <- 1
    fn_minus_2 <- 1
    steps_remaining <- n-2 
    while(steps_remaining > 0){
      steps_remaining <- steps_remaining-1
      
      fn <- fn_minus_1 + fn_minus_2
      fn_minus_2 <- fn_minus_1
      fn_minus_1 <- fn
    }
    fn
  }
}

system.time(fib_loop(1000))
system.time(fib(1000))


#4) MCMC ----
# a) 
#' @param x
#' @param alpha
#' @param beta
#' @description Calculate log-posterior of alpha, beta
log_posterior <- function(x, alpha, beta){
  dexp(alpha, rate=1, log=T) + # alpha prior density
  dexp(beta, rate=1, log=T) + # beta prior density
  sum(dgamma(x, shape = alpha, rate = beta, log=T)) # likelihood
}

# b) Single MH step
#' @param x
#' @param alpha
#' @param beta
#' @param sigma
single_mh_step <- function(x, alpha, beta, sigma){
  # Proposal
  prop_alpha <- rnorm(1, mean =alpha, sd=sigma)
  prop_beta <- rnorm(1, mean=beta, sd=sigma)
  if (any(c(prop_alpha, prop_beta)) < 0){
    return(list("alpha" = alpha,
         "beta"  = beta))
  }
  
  # Accept/ Reject
  numerator <- log_posterior(x, alpha = prop_alpha, beta = prop_beta) 
  denominator <-   log_posterior(x, alpha = alpha, beta = beta) 
  
  accept_prob <- exp(numerator-denominator)

  if(runif(1) < accept_prob){
    list("alpha" = prop_alpha, 
         "beta"  = prop_beta)
  } else 
    list("alpha" = alpha,
         "beta"  = beta)
}

# b) Multiple MH step
#' @param x
#' @param alpha
#' @param beta
#' @param sigma
#' @param n
mh_algo <- function(x, alpha, beta, sigma, n){
  latest_alpha <- alpha
  latest_beta <- beta
  
  result_store <- matrix(NA,nrow = n, ncol = 2)
  for (i in 1:n){
    mh_result <- single_mh_step(x, latest_alpha, latest_beta, sigma)
    latest_alpha <- mh_result$alpha
    latest_beta <- mh_result$beta
    result_store[i,] <- c(latest_alpha, latest_beta)
  }
  result_store
}


# d) Read file and plot hist
file_path <- file.path('../data', 'airpol.txt')
x <- scan(file_path)
hist(x, breaks = 100, freq = FALSE)


posterior <- mh_algo(x, alpha = 1, beta = 1, n=5000, sigma = 0.02)

# Estimates
alphahat <- colMeans(posterior)[1]
betahat <- colMeans(posterior)[2]

# Plots
hist(posterior[,1], main = 'Alpha Posterior')
hist(posterior[,2], main = 'Beta Posterior')

plot(posterior[,1], type = 'l')
points(posterior[,2], col= 'blue', type='l')

hist(x, breaks = 100, freq = FALSE) 

f <- function(y) 
  dgamma(y, alphahat, betahat) 
plot(f, 0, 70, add = TRUE, col = 2, lwd = 2)

# 5) Methods ----
# a) List of normal rvs and call is class biv
biv <- function(){
  value <- list("x"=rnorm(20), 'y'=rpois(20,5))
  attr(value,'class') <- 'biv'
  value
}

print.biv <- function(obj) { 
  n <- length(obj$x) 
  cat("Bivariate data,", n, "entries\n")
  len <- min(n, 6) 
  dots <- ifelse(n > 6, "...", "")
  cat("x : ", obj$x[1:len], dots, "\n") 
  cat("y : ", obj$y[1:len], dots, "\n")
  invisible(obj)
}

plot.biv <- function(obj,...){
  par(mfrow=c(1,3))
  plot(obj$x, obj$y, xlab = 'x', ylab= 'y')
  boxplot(obj$x, main = 'x') 
  boxplot(obj$y, main ='y')
  
  par(mfrow=c(1,1))
}

plot(obj)
biv_inst <- biv()
plot(1:10)


# low leverl constructor
.Biv <- setClass("Biv", 
                    slots = c(
                      x = "numeric", 
                      y = "numeric"
                    )
)


# helper constructor
Biv <- function(x=rnorm(20), y=rpois(20,5)){
  # validation here on arguments
  
  # constructor
  .Biv(x=x,y=y)
}

# print, is "show" in S4
setMethod("show", "Biv", function(object) {
  cat("Bivariate data \n")
  
  
  nx <- length(object@x)
  ny <- length(object@y)
  
  lenx <- min(nx, 6) 
  leny <- min(ny, 6) 
  
  dotsx <- ifelse(nx > 6, "...", "")
  dotsy <- ifelse(ny > 6, "...", "")
  
  cat(is(object)[[1]], "\n")
  cat("x : ", object@x[1:lenx], dotsx, "\n") 
  cat("y : ", object@y[1:leny], dotsy, "\n")
  
})


setMethod("plot", "Biv", function(object) {
  par(mfrow=c(1,3))
  plot(object@x, object@y, xlab = 'x', ylab= 'y')
  boxplot(object@x, main = 'x') 
  boxplot(object@y, main ='y')
  
  par(mfrow=c(1,1))
})

setGeneric("print", function(object, ...) standardGeneric("print"))

# print, is "show" in S4
setMethod("print", "Biv", function(object) {
  cat("Bivariate data \n")
  
  
  nx <- length(object@x)
  ny <- length(object@y)
  
  lenx <- min(nx, 6) 
  leny <- min(ny, 6) 
  
  dotsx <- ifelse(nx > 6, "...", "")
  dotsy <- ifelse(ny > 6, "...", "")
  
  cat(is(object)[[1]], "\n")
  cat("x : ", object@x[1:lenx], dotsx, "\n") 
  cat("y : ", object@y[1:leny], dotsy, "\n")
  
  invisible(object)
  
})

obj <- Biv()
plot(obj)

# 6) Functions ----
# a) Return matrix of specific form

#' @param x, numeric vector
#' @param z, numeric vector
#' @description return matrix with first column 1, 
#' second column x, third z 
#' and fourth the element-wise project 
#' @example gen_matric(1:10,11:20)
gen_matrix <- function(x, z) {
 init_matrix <- matrix(NA, nrow=min(length(x), length(z)), ncol = 4)
 init_matrix[]
}