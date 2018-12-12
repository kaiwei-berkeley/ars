library(rlang)
library(numDeriv)
library(testthat)
### Checking function
convert_log <- function(f) {
  # constructing the log of the input function
  log_f <- function(x) {}
  old_f <- rlang::get_expr(body(f))
  body(log_f) <- rlang::get_expr(quo(log(!!(old_f))))
  return(log_f)
}

log_concave_check <- function(log_f, x1, xk) {
  # Checking the gradient of every consecutive cords
  x_val = seq(from = x1, to = xk, by = 0.1)
  x_val_initial <- x_val[1:length(x_val)-1]
  x_val_initial_plus_1 <- x_val[2:length(x_val)]

  grad_vec <- (log_f(x_val_initial_plus_1) - log_f(x_val_initial))/(x_val_initial_plus_1 - x_val_initial)
  dif_grad <- grad_vec[1:length(grad_vec)-1] - grad_vec[2:length(grad_vec)]
  dif_grad <- dif_grad[is.finite(dif_grad)]
  dif_grad <- round(dif_grad, digits = 8)

  if (all(grad_vec== 0)) {
    # case 1 is uniform distribution
    case = 1
  } else if (all(dif_grad == 0)) {
    # case 2 is exp distribution
    case = 2
  }  else if ((min(dif_grad) >= 0) == TRUE) {
    # case 3 is dsitribution that satisfy log-concavity check
    case = 3
  }  else {
    # case 4 is the distribution does not work
    case = 4
  }
  return(case)
}

domain_check <- function(log_f,D_lower, D_upper, x_lower, x_upper) {
  if ((D_lower > D_upper)==TRUE) stop("Invalid domain")
  else {

    # here we force the user to input x_lower (x_1) and x_upper (x_k) to be valid
    if (x_lower >= x_upper || is.infinite(x_lower) || is.infinite(x_upper)
        || x_lower < D_lower || x_upper > D_upper) stop("Invalid lower and upper values")
    else {
      if ((D_lower == -Inf) == TRUE) {
        deriv_val <- numDeriv::grad(func = log_f,
                                    x = x_lower, method = "Richardson")
        if ((deriv_val <= 0) == TRUE) stop("Derivative at x_1 must be positive")
      }
      if ((D_upper == Inf) == TRUE) {
        deriv_val <- numDeriv::grad(func = log_f, x = x_upper, method = "Richardson")
        if ((deriv_val >= 0) == TRUE) stop("Derivative at x_k must be negative")
      }
    }
  }
}

### Define initialize function
initial_Tk = function(f, x1, xk, k = 3) {
  # x1 and xk should be in domain
  Tk = seq(from = x1, to = xk, length.out = k)
  return(Tk)
}

### Define updating function
update_Tk = function(Tk, new_x){
  return (sort(c(Tk, new_x)))
}

##This takes updated Tk and updated z and return the updated u function
update_u = function(f, Tk, z) {
  f = f
  k = length(Tk)
  u = function(x) {
    v <- findInterval(x, z)
    index <- ifelse(v > k, k, v)
    df = numDeriv::grad(func = f,x = Tk[index], method = "simple")
    return(f(Tk[index]) + (x - Tk[index]) * df)
  }
  return(u)
}

##This takes updated Tk and returns the updated l function
update_l = function(f, Tk) {
  f = f
  l = function(x) {
    index = ifelse(x < min(Tk), 0,
                   ifelse(x == min(Tk), 1,
                          max(which(x>Tk))))
    out = ifelse(index == 0 || index == length(Tk),
                 -Inf,
                 ((Tk[index+1]-x) * f(Tk[index]) + (x-Tk[index]) * f(Tk[index+1]))/
                   (Tk[index+1]-Tk[index]))
    return(out)
  }
  return(l)
}

##############################12/04 update

### This function takes Tk and return a list contain:
# 1. z0 - zk (k+1 elements)
# 2. slopes evaluated at x1-xk (k elements)
initial_zlist = function(f,Tk,start,end){
  ## initialize z and slope
  k = length(Tk)
  z = numeric(k+1)
  # get all the slope
  slope = numDeriv::grad(func = f,x = Tk, method = "simple")
  ## Get z0 and zk
  z[1] = start
  z[k+1] = end

  ## get z_1 to z_k-1,using Tk[1] - Tk[k-1]
  for (j in 2:k) {
    df1 = slope[j-1]
    df2 = slope[j]
    z[j] = (f(Tk[j])-f(Tk[j-1])-Tk[j]*df2+Tk[j-1]*df1)/(df1-df2)
  }
  zlist = list(z,slope)
  return(zlist)
}

### call this before update Tk
### given previous z list, previous Tk and x_star; return a list of updated z, and updated slope
update_zlist = function(f,zlist,Tk,x){
  ind = ifelse(x>max(Tk),length(Tk)+1,min(which(x<Tk)))
  df_new = (f(x+0.0001)-f(x))/0.0001

  if (ind == length(Tk)+1) {
    z = zlist[[1]]; slope = zlist[[2]]
    Tk = update_Tk(Tk,new_x = x)
    j = length(Tk)

    slope[length(Tk)] = df_new

    df1 = slope[j-1]; df2 = slope[j]
    z[length(Tk)] = (f(Tk[j])-f(Tk[j-1])-Tk[j]*df2+Tk[j-1]*df1) / (df1-df2)

    zlist = list(z, slope)
    return(zlist)
  }

  Tk = update_Tk(Tk,new_x = x) ## room for improvement
  z = zlist[[1]]
  slope = zlist[[2]]
  ## update slope vector
  #df_new = numDeriv::grad(func = f,x = x, method = "simple")
  slope[(ind+1):(length(slope)+1)] = slope[ind:length(slope)]
  slope[ind] = df_new

  ## update z using new slope
  ## Get z0 and zk
  k = length(Tk)
  ## get z_1 to z_k-1,using Tk[1] - Tk[k-1]

  df1 <- slope[1:k-1]
  df2 <- slope[2:k]
  z[2:k] = (f(Tk[2:k])-f(Tk[1:k-1])-Tk[2:k]*df2+Tk[1:k-1]*df1)/(df1-df2)

  zlist = list(z,slope)
  return(zlist)
}

### sample from u function
# function to sample x, the function samp_ars is to get the area for each segment
# then obtain the cdf from this so that we can sample x.
# First generate one number to choose the line segment we will use to sample

samp_ars  = function(f,Tk,start,end,zlist){
  u1 = runif(1,0,1)
  # Initialize values we are going to use
  k = length(Tk); df = zlist[[2]];
  intercept = numeric(); area = numeric(); left_val = numeric(); right_val = numeric();
  slope = numeric()
  # left_val and right_val are vectors we use to identify the starting and ending point of the segment we want to sample from
  j = 1

  # If -Inf is the lower bound then use
  # the slope of Tk[1] and calculate the intercept
  if (is.infinite(start)) {
    intercept[1] = f(Tk[1]) - df[1]*Tk[1]

    # area obtained by integrate from -Inf to Tk[1] as well as the left and
    # right value for each segment
    area[1] = exp(intercept[1])/df[1]*(exp(df[1]*Tk[1])-0)
    left_val[1] = -Inf
    right_val[1] = Tk[1]

    slope[1] = df[1]
    # count variable for all the initial vector
    j = length(area) + 1
  } else {
    slope1 = df[1]
    intercept1 = f(Tk[1]) - slope1*Tk[1]
    area1 = exp(intercept1)/slope1*(exp(slope1*Tk[1]) - exp(slope1*zlist[[1]][1]))
    area[j] = area1
    left_val[j] = zlist[[1]][1]
    right_val[j] = Tk[1]
    slope[j] = slope1

    j = length(area) + 1
  }

  # Loops all over the interior segments
  for (i in 1:(length(Tk)-1)) {

    # The left slope and intercept of each upper hull
    df1 = df[i]
    intercept1 = f(Tk[i]) - df1*Tk[i]

    # The right slope and intercept of each upper hull
    df2 = df[i+1]
    intercept2 = f(Tk[i+1]) - df2*Tk[i+1]

    # intersection point
    z = zlist[[1]][i+1]

    # get the left area
    pr1 = exp(intercept1)/df1*(exp(df1*z) - exp(df1*Tk[i]))


    # storing all the values for the left upper hull
    area[j] = pr1; slope[j] = df1; intercept[j] = intercept1; left_val[j] = Tk[i]; right_val[j] = z

    # increase the count
    j = length(area) + 1

    # get the right area
    pr2 = exp(intercept2)/df2*(exp(df2*Tk[i+1]) - exp(df2*z))

    # storing all the values for the right upper hull
    area[j] = pr2; slope[j] = df2; intercept[j] = intercept2; left_val[j] = z; right_val[j] = Tk[i+1]

    j = length(area) + 1
  }



  # If -Inf is the upper bound then use
  # the slope of Tk[1] and calculate the intercept
  if (is.infinite(end)) {
    slope[j] = df[k]
    intercept[j] = f(Tk[k]) - slope[j]*Tk[k]
    area[j] = exp(intercept[j])/slope[j]*(0-exp(slope[j]*Tk[k]))
    left_val[j] = Tk[k]
    right_val[j] = Inf
  } else {
    slope_k = df[k]
    intercept_k = f(Tk[k]) - slope_k*Tk[k]
    area_k = exp(intercept_k)/slope_k*(exp(slope_k*end) - exp(slope_k*Tk[k]))
    area[j] = area_k
    left_val[j] = Tk[k]
    right_val[j] = end
    slope[j] = slope_k
  }

  # Summing all vector and create the probability vector of each segment. Then
  # create the cdf of the upper hull function
  T = sum(area); prob_vec = area/T; cdf = cumsum(prob_vec)
  len = length(prob_vec)
  # get the segment index of the segment we want to use


  if (length(which(u1<=cdf)) == 0) {
    ind = sample(c(1:len), 1)
  } else {
    ind = min(which(u1<=cdf))
  }

  # retrieve the values for slope, intercept, left, and right values of the segment chosen
  m = slope[ind]; b = intercept[ind]; left = left_val[ind]; right = right_val[ind]

  # generate the random value
  u2 = runif(1,0,1);
  x_star = log(u2*(exp(m*right) - exp(m*left)) + exp(m*left))/m
  return(x_star)
}

samp_ars3  = function(f,Tk,start,end,zlist){
  # Initialize values we are going to use
  k = length(Tk); df = zlist[[2]];
  intercept = numeric(); area = numeric(); left_val = numeric(); right_val = numeric();
  slope = numeric()
  # left_val and right_val are vectors we use to identify the starting and ending point of the segment we want to sample from
  j = 1

  # If -Inf is the lower bound then use
  # the slope of Tk[1] and calculate the intercept
  if (is.infinite(start)) {
    intercept[1] = f(Tk[1]) - df[1]*Tk[1]

    # area obtained by integrate from -Inf to Tk[1] as well as the left and
    # right value for each segment
    area[1] = exp(intercept[1])/df[1]*(exp(df[1]*Tk[1])-0)
    left_val[1] = -Inf
    right_val[1] = Tk[1]

    slope[1] = df[1]
    # count variable for all the initial vector
    j = length(area) + 1
  } else {
    slope1 = df[1]
    intercept1 = f(Tk[1]) - slope1*Tk[1]
    area1 = exp(intercept1)/slope1*(exp(slope1*Tk[1]) - exp(slope1*zlist[[1]][1]))
    area[j] = area1
    left_val[j] = zlist[[1]][1]
    right_val[j] = Tk[1]
    slope[j] = slope1

    j = length(area) + 1
  }

  # Loops all over the interior segments
  for (i in 1:(length(Tk)-1)) {

    # The left slope and intercept of each upper hull
    df1 = df[i]
    intercept1 = f(Tk[i]) - df1*Tk[i]

    # The right slope and intercept of each upper hull
    df2 = df[i+1]
    intercept2 = f(Tk[i+1]) - df2*Tk[i+1]

    # intersection point
    z = zlist[[1]][i+1]

    # get the left area
    pr1 = exp(intercept1)/df1*(exp(df1*z) - exp(df1*Tk[i]))


    # storing all the values for the left upper hull
    area[j] = pr1; slope[j] = df1; intercept[j] = intercept1; left_val[j] = Tk[i]; right_val[j] = z

    # increase the count
    j = length(area) + 1

    # get the right area
    pr2 = exp(intercept2)/df2*(exp(df2*Tk[i+1]) - exp(df2*z))

    # storing all the values for the right upper hull
    area[j] = pr2; slope[j] = df2; intercept[j] = intercept2; left_val[j] = z; right_val[j] = Tk[i+1]

    j = length(area) + 1
  }



  # If -Inf is the upper bound then use
  # the slope of Tk[1] and calculate the intercept
  if (is.infinite(end)) {
    slope[j] = df[k]
    intercept[j] = f(Tk[k]) - slope[j]*Tk[k]
    area[j] = exp(intercept[j])/slope[j]*(0-exp(slope[j]*Tk[k]))
    left_val[j] = Tk[k]
    right_val[j] = Inf
  } else {
    slope_k = df[k]
    intercept_k = f(Tk[k]) - slope_k*Tk[k]
    area_k = exp(intercept_k)/slope_k*(exp(slope_k*end) - exp(slope_k*Tk[k]))
    area[j] = area_k
    left_val[j] = Tk[k]
    right_val[j] = end
    slope[j] = slope_k
  }

  # Summing all vector and create the probability vector of each segment. Then
  # create the cdf of the upper hull function
  T = sum(area); prob_vec = area/T; cdf = cumsum(prob_vec)
  len = length(prob_vec)

  # get the segment index of the segment we want to use
  u1 = runif(1,0,1)
  if (length(which(u1<=cdf)) == 0) {
    ind = sample(c(1:len), 1)
  } else {
    ind = min(which(u1<=cdf))
  }

  # retrieve the values for slope, intercept, left, and right values of the segment chosen
  m = slope[ind]; b = intercept[ind]; left = left_val[ind]; right = right_val[ind]

  # generate the random value
  u2 = runif(1,0,1);

  # calculate the x_star values
  # if(m == 0){
  #   x_star = (right-left)*u2 + left
  # }else{
  # }
  x_star = log(u2*(exp(m*right) - exp(m*left)) + exp(m*left))/m
  if(is.infinite(x_star)){
    stop("Generated numbers that exceed machine maximum, run again or change the input h(x)")
  }
  return(x_star)
}



samp_ars2  = function(f,Tk,start,end,zlist){
  u1 = runif(1,0,1)
  # Initialize values we are going to use
  k = length(Tk);
  df = zlist[[2]];
  intercept = numeric(); area = numeric();
  left_val = numeric(); right_val = numeric(); slope = numeric()
  # left_val and right_val are vectors we use to identify the starting and ending point of the segment we want to sample from
  j = 1

  # If -Inf is the lower bound then use
  # the slope of Tk[1] and calculate the intercept
  if (is.infinite(start)) {
    intercept[1] = f(Tk[1]) - df[1]*Tk[1]

    # area obtained by integrate from -Inf to Tk[1] as well as the left and
    # right value for each segment
    area[1] = exp(intercept[1])/df[1]*(exp(df[1]*Tk[1])-0)
    left_val[1] = -Inf
    right_val[1] = Tk[1]

    slope[1] = df[1]
    # count variable for all the initial vector
    j = length(area) + 1

  } else {
    slope1 = df[1]
    intercept1 = f(Tk[1]) - slope1*Tk[1]
    area1 = exp(intercept1)/slope1*(exp(slope1*Tk[1]) - exp(slope1*zlist[[1]][1]))
    area[j] = area1
    left_val[j] = zlist[[1]][1]
    right_val[j] = Tk[1]
    slope[j] = slope1

    j = length(area) + 1
  }

  # For each interior segment
  index = 1:(length(Tk)-1)
  index_plus = index+1

  # The left slope and intercept of each upper hull
  df1 <- df[index]
  intercept_total <- f(Tk) - df*Tk
  intercept1 <- intercept_total[index]
  # The right and intercept of each upper hull
  df2 = df[index+1]
  intercept2 <- intercept_total[index_plus]
  # intersection point
  zs = zlist[[1]][index_plus]

  # get the left area
  pr1 = exp(intercept1)/df1*(exp(df1*zs) - exp(df1*Tk[index]))

  indices <- seq(from=1, to=2*length(pr1)-1, by=2)
  area[indices] = pr1;
  slope[indices] = df1;
  intercept[indices] = intercept1;
  left_val[indices] = Tk[index];
  right_val[indices] = zs

  pr2 = exp(intercept2)/df2*(exp(df2*Tk[index+1]) - exp(df2*zs))
  indices <- indices + 1
  area[indices] = pr2;
  slope[indices] = df2;
  intercept[indices] = intercept2;
  left_val[indices] = zs;
  right_val[indices] = Tk[index+1]

  # If -Inf is the upper bound then use
  # the slope of Tk[1] and calculate the intercept
  j = length(area) + 1
  if (is.infinite(end)) {
    slope[j] = df[k]
    intercept[j] = f(Tk[k]) - slope[j]*Tk[k]
    area[j] = exp(intercept[j])/slope[j]*(0-exp(slope[j]*Tk[k]))
    left_val[j] = Tk[k]
    right_val[j] = Inf
  } else {
    slope_k = df[k]
    intercept_k = f(Tk[k]) - slope_k*Tk[k]
    area_k = exp(intercept_k)/slope_k*(exp(slope_k*end) - exp(slope_k*Tk[k]))
    area[j] = area_k
    left_val[j] = Tk[k]
    right_val[j] = end
    slope[j] = slope_k
  }

  # Summing all vector and create the probability vector of each segment. Then
  # create the cdf of the upper hull function
  T = sum(area); prob_vec = area/T; cdf = cumsum(prob_vec)
  len = length(prob_vec)
  # get the segment index of the segment we want to use

  # here I tried to avoid the error Error in if (index[i] == 0 | index[i] == length(Tk)) { :
  # missing value where TRUE/FALSE needed
  # In addition: There were 11 warnings (use warnings() to see them)
  ##########################################
  # ---- but it seems like it does not work


  # if (is.infinite(min(which(u1<cdf)))) {
  #  ind = sample(c(1:len), 1)
  # } else {
  #  ind = min(which(u1<cdf))}
  m <- which(u1 <= cdf)
  if (length(m) == 0) {
    ind = sample(c(1:len), 1)
  } else {
    ind = min(m)
  }

  # retrieve the values for slope, intercept, left, and right values of the segment chosen
  m = slope[ind]; b = intercept[ind];
  left = left_val[ind]; right = right_val[ind]

  # generate the random value
  u2 = runif(1,0,1);

  x_star = log(u2*(exp(m*right) - exp(m*left)) + exp(m*left))/m
  return(x_star)
}

# this function is used when either of the lower bound or upper bound is bounded
# i.e. D =[-Inf,a] or [a,Inf] or [-Inf,Inf]
log_odd_transform = function(f) {
  # constructing the log of the input function
  log_odd_f = function(x) {}
  old_f = rlang::get_expr(body(f))
  body(log_odd_f) = rlang::get_expr(quo(log(!!(old_f)/(1-!!(old_f)))))
  return(log_odd_f)
}

###### use if unbounded on the left

optimal_point_left <- function(f,expected_val){
  # f : log_f function
  # expected_val: maximum value
  optimal = FALSE
  der_val = numDeriv::grad(func = f, x = expected_val, method = "Richardson")
  n = 10
  iteration = 100
  while (optimal == FALSE && iteration > 0) {
    if (der_val < 5 && der_val > 3) {
      der_val = der_val
      optimal = TRUE
    } else if (der_val >= 5) {
      expected_val = expected_val + n
      der_val = numDeriv::grad(func = f, x = expected_val, method = "Richardson")
    } else {
      expected_val = expected_val - n
      der_val = numDeriv::grad(func = f, x = expected_val, method = "Richardson")
    }
    iteration = iteration - 1
    n = n/2
  }
  return(expected_val)
}

###### use if unbounded on the right

optimal_point_right <- function(f,expected_val){
  # f : log_f function
  # expected_val: maximum value
  optimal = FALSE
  der_val = numDeriv::grad(func = f, x = expected_val, method = "Richardson")
  n = 10
  iteration = 1000
  while (optimal == FALSE && iteration > 0) {
    if (der_val < -3 && der_val > -5) {
      der_val = der_val
      optimal = TRUE
    } else if (der_val <= -5) {
      expected_val = expected_val - n
      der_val = numDeriv::grad(func = f, x = expected_val, method = "Richardson")
    } else {
      expected_val = expected_val + n
      der_val = numDeriv::grad(func = f, x = expected_val, method = "Richardson")
    }
    iteration = iteration - 1
    n = n/2
  }
  return(expected_val)
}





# this function is used when either of the lower bound or upper bound is bounded
# i.e. D =[-Inf,a] or [a,Inf] or [-Inf,Inf]
initial_point_sample = function(f,start,end) {
  # f should be the orginal function ####
  h = convert_log(f)

  # convert using log_odd transformation
  optimal_approx = mode_finding(f,start)


  # finding optimal left point
  if (is.infinite(start)) {
    x1 = optimal_point_left(h,optimal_approx)
  } else {x1 = start + 0.01}
  # finding optimal right point
  if (is.infinite(end)) {
    xk = optimal_point_right(h,optimal_approx)
  } else {xk = end - 0.01}

  result = c(x1,optimal_approx,xk)
  return(result)
}

#### optimize function
mode_finding = function(f,start) {
  if (is.infinite(start)) {
    start_new = -100000
  } else {start_new = start}
  stop1 = FALSE
  while (stop1 == FALSE) {
    a = seq(from = start_new, to = (start_new + 50.001), length.out = 500)
    fa = f(a)
    b = seq(from = (start_new + 50.001), to = (start_new + 100.002), length.out = 500)
    fb = f(b)
    max_val_a = max(fa)
    max_val_b = max(fb)
    max_val_ind_a = which.max(fa)
    max_val_ind_b = which.max(fb)
    if (max_val_a > max_val_b) {
      stop1 = TRUE
    } else {
      start_new = start_new + 100.002
    }
  }
  start_opt = a[max_val_ind_a]
  end_opt = b[max_val_ind_b]
  a = seq(from = start_opt, to = end_opt, by = 0.5)
  fa = f(a)

  max_val_a = max(fa)

  max_val_ind_a = which.max(fa)
  mode_point = a[max_val_ind_a]
  return(mode_point + 0.01)
}
