library(assertthat)
library(testthat)
set.seed(42)

# Unit tests
# Test1. log_concave_check
print("Test1. log_concave_check")
test_that("log_concave_check", {
  # exponential
  f <- function(x) { dexp(x, rate=1) }
  log_f <- convert_log(f)
  expect_equal(log_concave_check(log_f, 0.1, 5),2)

  # inverse logit function
  f <- function(x) { exp(x)/(exp(x)+1) };
  log_f <- convert_log(f);
  expect_equal(log_concave_check(log_f, 0.1, 5),3)

  # x^a function
  f <- function(x) { x^5 }
  log_f <- convert_log(f)
  expect_equal(log_concave_check(log_f, 0.1, 5),3)
})


# Test2. initial_tk
print("Test2. initial_tk")

test_that("x^a function",
          {
            tk <- initial_Tk(dnorm, -5, 5)
            expect_true(all(tk >= -5))
            expect_true(all(tk <= 5))
          })



# Test3. update_u and zlist
print("Test3. update_u and zlist")
test_that("update_u and zlist",{

  left <- -5
  right <- 5
  f <- function(x) { dnorm(x, mean=0, sd=1) }
  log_f = convert_log(f)
  tk <- initial_Tk(dnorm, left, right)
  assert_that(all(tk == c(-5, 0, 5)) == TRUE)
  z <- initial_zlist(log_f, tk, left, right)
  # check zs are within bounds
  expect_true(all(z[[1]] >= left))
  expect_true(all(z[[1]] <= right))
  init_u <- update_u(log_f, tk, z[[1]])
  init_l <- update_l(log_f, tk)

  # Check that we built a proper upper and lower bound
  s <- seq(from=left, to=right, by=0.1)
  for (x in s) {
    u <- init_u(x)
    l <- init_l(x)
    expect_true(u >= l - .Machine$double.eps)
    expect_true(log_f(x) >= l - .Machine$double.eps)
    expect_true(log_f(x) <= u + .Machine$double.eps)
  }
})

# Test 4. update_Tk
print("Test4. update_Tk")
test_that("update_Tk",{
  left <- -5
  right <- 5
  f <- function(x) { dnorm(x, mean=0, sd=1) }
  log_f = convert_log(f)
  tk <- initial_Tk(dnorm, left, right)
  expect_true(all(update_Tk(tk, 3) == c(-5, 0, 3, 5)) == TRUE)
})

# Test5. samp_ars
print("Test5. samp_ars")
test_that("samp_ars",{

  left <- -5
  right <- 5
  f <- function(x) { dnorm(x, mean=0, sd=1) }
  log_f = convert_log(f)
  tk <- initial_Tk(dnorm, left, right)
  z <- initial_zlist(log_f, tk, left, right)
  xs <- replicate(1000, samp_ars(log_f, tk, start = left, end = right, zlist = z))
  # make sure all the sampled values are within bounds
  expect_true(all(xs <= right) == TRUE)
  expect_true(all(xs >= left) == TRUE)

})

# Test6. mode_finding
print("Test6. mode_finding")
test_that("mode_finding",{
  f <- function(x) { dnorm(x, mean=0, sd=1) }
  expect_true(abs(mode_finding(f, -10)) < 10e-1)
  f <- function(x) { dnorm(x, mean=10, sd=10) }
  expect_true(abs(mode_finding(f, -100) - 10) < 10e-1)

  # mode of beta is given by (shape1-1)/(shape1+shape2-2)
  f = function(x) {dbeta(x,shape1 = 2 ,shape2 = 3)}
  mode = (2 - 1) / (2 + 3 - 2)
  expect_true(abs(mode_finding(f, -1) - mode) < 10e-1)

})



# Test7. ars result ks.test
print("Test7. ars result ks.test")
test_that("ars result ks.test",{
  ## sample from normal
  f = function(x) {dnorm(x, mean=10, sd=1)}
  start = -Inf
  end = Inf
  sample = ars(f,start =start,end = end, N =1000, k=4)
  t <- ks.test(sample, rnorm(1000, mean=10, sd=1))
  expect_gt(t$p.value,0.05)

  ## Sample from beta
  shape1 <- 2
  shape2 <- 3
  f = function(x) {dbeta(x,shape1 = 2, shape2 = 3)}
  start = 0
  end = 1
  sample = ars(f,start =start,end = end,N = 1000,k = 4)
  t <- ks.test(sample, rbeta(1000, shape1=shape1, shape2=shape2))
  expect_gt(t$p.value,0.05)
})
