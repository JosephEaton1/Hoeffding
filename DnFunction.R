Dn <- function(sample){
  #function calculates Dn (accounting fir ties) and zeta for given sample
  x <- sample[[1]]
  y <- sample[[2]]
  n <- length(x)
  
  dx <- outer(x, x, "-")
  dy <- outer(y, y, "-")
  
  #quadrant contributions and counting ties
  tie_x <- dx == 0
  tie_y <- dy == 0
  diag(tie_x) <- FALSE
  diag(tie_y) <- FALSE
  tie_both <- tie_x & tie_y
  tie_x_only <- tie_x & !tie_y
  tie_y_only <- tie_y & !tie_x
  
  a_mat <- ((dx > 0) & (dy < 0)) +
    (tie_x_only & (dy < 0)) * 0.5 +
    (tie_y_only & (dx > 0)) * 0.5 +
    tie_both * 0.25
  
  b_mat <- ((dx < 0) & (dy < 0)) +
    (tie_x_only & (dy < 0)) * 0.5 +
    (tie_y_only & (dx < 0)) * 0.5 +
    tie_both * 0.25
  
  c_mat <- ((dx > 0) & (dy > 0)) +
    (tie_x_only & (dy > 0)) * 0.5 +
    (tie_y_only & (dx > 0)) * 0.5 +
    tie_both * 0.25
  
  d_mat <- ((dx < 0) & (dy > 0)) +
    (tie_x_only & (dy > 0)) * 0.5 +
    (tie_y_only & (dx < 0)) * 0.5 +
    tie_both * 0.25
  
  a <- rowSums(a_mat)
  b <- rowSums(b_mat)
  c <- rowSums(c_mat)
  d <- rowSums(d_mat)
  
  sigma <- a*(a-1)*d*(d-1) + b*(b-1)*c*(c-1) - 2*a*b*c*d
  D <- sum(sigma) / (n*(n-1)*(n-2)*(n-3)*(n-4))
  
  #calculating zeta estimat
  w <- b*(c+1) - a*d
  g <- c + 1 + d
  f <- c + 1 + a
  
  one_mat   <- (dx <= 0) & (dy <= 0)
  two_mat   <- (dx <= 0)
  three_mat <- (dy <= 0)
  diag(one_mat)   <- TRUE
  diag(two_mat)   <- TRUE
  diag(three_mat) <- TRUE
  
  U <- 2 * sum(w * (c + 1)) / (5 * n^4)
  
  S <- w^2 / (5 * n^4)
  T_vec <- 2 * as.vector(one_mat %*% w) / (5 * n^3)
  V_vec <- -2 * as.vector(two_mat %*% (w * g)) / (5 * n^4)
  W_vec <- -2 * as.vector(three_mat %*% (w * f)) / (5 * n^4)
  
  zeta <- S + T_vec + V_vec + W_vec + U
  zeta_mean <- mean(zeta)
  zeta_variance <- var(zeta)
  Dn_variance <- 25*zeta_variance /n
  
  table <- data.frame(xi = x, yi = y, ai = a, bi = b, ci = c, di = d)
  
  return(list(
    details = table,
    Dn = D,
    zeta_mean = zeta_mean,
    zeta_variance = zeta_variance,
    Dn_variance = Dn_variance
  ))
}