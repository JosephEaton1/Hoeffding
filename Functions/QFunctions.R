#no need to go past T = 0.5 we found the tail functions are 'computationally' zero after that
#can change h to bigger step size for quicker, will affect accuracy of final CDF
#these are values we used for paper, note the computation took a couple hours to get 
#all the way to 110, we found all central moments, except mean, were accurate after only 20 

h <- 0.0001    # note when we did BKR comparisons we changed this to 0.01/pi^4 so we could do direct comparisons 
T <- 0.5
vec <- seq(from = 0, to = T, by = h)
n <- length(vec)
k_vals <- 1:50
lower_bounds <- (2*k_vals - 1) * pi
upper_bounds <- 2*k_vals * pi

compute_Ck_j <- function(j, vec, lower_bounds, upper_bounds, n) {
  #function computes cdf function of C(x(pij)^2) for values in vec
  #note this is just our CvM function with added scale factor
  
  Ck <- numeric(n)
  Ck[1] <- 0
  scale_factor_base <- pi^2 * j^2 / 2
  
  for (i in 2:n) {
    scale_factor <- scale_factor_base * vec[i]
    terms <- numeric(50)
    
    for (k in 1:50) {
      g <- function(r) exp(-r^2 * scale_factor) / sqrt(-r * sin(r))
      int_result <- integrate(g, lower_bounds[k], upper_bounds[k],
                              subdivisions = 1000, rel.tol = 1e-7)
      terms[k] <- (-1)^(k-1) * int_result$value
    }
    
    Ck[i] <- 1 - (sum(terms) * 2 / pi)
  }
  #sometimes function gave negative outputs very close to zero due to the singularity
  return(pmax(Ck, 0))
}

#this is longest step, maybe can speed up code somewhere above? or compute distributions
#parallel as they are independent.
Ck_list <- list()
for (j in 1:110) {
  cat('computing', j, '\n')
  Ck_list[[j]] <- compute_Ck_j(j, vec, lower_bounds, upper_bounds, n)
}


#vec_list[[i]] is the convolution of the partial sum to i, as said before we went
# to 110 for most accurate results but maybe if using a bigger step size the 
#convolution error would build up so smaller partial sum more accurate? 

M <- length(vec)
vec_list <- list()
vec_list[[1]] <- Ck_list[[1]]

for (j in 2:110) {
  cat('convolution', j, '\n')
  
  vj <- numeric(M)
  vj[1] <- 0
  dF <- diff(Ck_list[[j]])
  
  weights_all <- numeric(M)
  weights_all[1] <- dF[1]/2
  if (M > 2) {
    weights_all[2:(M-1)] <- (dF[2:(M-1)] + dF[1:(M-2)])/2
  }
  weights_all[M] <- dF[M-1]/2
  
  for (l in 2:M) {
    F_prev <- vec_list[[j-1]][l:2]
    weights <- weights_all[1:(l-1)]
    vj[l] <- sum(F_prev * weights)
  }
  
  vec_list[[j]] <- vj
}

min_len <- M - 110   # NA builds up from convolution so need to chop off to integrate
vec_chopped <- vec[1:min_len]

for (i in 1:110) {
  vec_list[[i]] <- vec_list[[i]][1:min_len]
}

trapz <- function(x, y) {
  sum(diff(x) * (y[-1] + y[-length(y)]) / 2)
}

#use this to test accuracy of partial sum
moment_calc <-function(k){
  #calculates moments of Q for given partial sum to k
  Q_func <- 1 - vec_list[[k]]
  
  EQ <- trapz(vec_chopped, Q_func)
  
  xQ_func <- Q_func * 2 * vec_chopped
  EQ2 <- trapz(vec_chopped, xQ_func)
  
  q3_func <- Q_func * 3 * vec_chopped^2
  EQ3 <- trapz(vec_chopped, q3_func)
  
  q4_func <- Q_func * 4 * vec_chopped^3
  EQ4 <- trapz(vec_chopped, q4_func)
  
  VarQ <- EQ2 - EQ^2
  
  mu <- EQ
  mu2 <- VarQ
  mu3 <- EQ3 - 3*mu*EQ2 + 2*mu^3
  mu4 <- EQ4 - 4*mu*EQ3 + 6*mu^2*EQ2 - 3*mu^4
  
  Skewness <- mu3 / (mu2^(3/2))
  Kurtosis <- mu4 / (mu2^2)
  ExcessKurtosis <- Kurtosis - 3
  return(list(
    Mean = mu,
    Variance = mu2,
    Skewness = Skewness,
    Kurtosis = Kurtosis,
    ExcessKurtosis = ExcessKurtosis
  ))
}


true_values <- data.frame(
  k = "Theoretical",
  E_Q = 1/36,
  Var_Q = 2/(90^2),
  Skewness = 80*sqrt(2)/49,
  Kurtosis = 579/49,
  Excess_Kurtosis = 432/49
)
