#run all these functions before Dn_table as it uses them all
#we split steps into seperate functions for clarity

local_values <- function(P) {
  #function computes the local abcd values
  rows <- nrow(P)
  cols <- ncol(P)
  C <- matrix(0L, rows, cols)
  if (cols >= 2) {
    for (i in 2:cols) {
      for (j in 1:(i - 1)) {
        C[, i] <- C[, i] + as.integer(P[, j] < P[, i])
      }
    }
  }
  xrank_ <- matrix(rep(0:(cols - 1), each = rows), nrow = rows, ncol = cols)
  A <- xrank_ - C
  D <- (P - 1L) - C
  B <- (cols - 1L) - A - C - D
  list(a_i = A, b_i = B, c_i = C, d_i = D)
}

S_value <- function(a_i, b_i, c_i, d_i) {
  rowSums(a_i * (a_i - 1) * d_i * (d_i - 1) + b_i * (b_i - 1) * c_i * (c_i - 1) - 2 * a_i * b_i * c_i * d_i)
}

freq_table <- function(m, upper_correct, lower_correct, mode) {
  #function globalises bd/ac depending on mode and produces freq table of S_k(p) or compliment
  #works for A and Ac #magic#
  
  P <- arrangements::permutations(m)
  rows <- nrow(P)
  loc <- local_values(P)
  a_i <- loc$a_i
  b_i <- loc$b_i
  c_i <- loc$c_i
  d_i <- loc$d_i
  
  if (mode == "bd") {
    d_i <- d_i + matrix(upper_correct[P], rows, m)
    b_i <- b_i + matrix(lower_correct[P], rows, m)
  } else {
    a_i <- a_i + matrix(lower_correct[P], rows, m)
    c_i <- c_i + matrix(upper_correct[P], rows, m)
  }
  
  S_values <- S_value(a_i, b_i, c_i, d_i)
  table(S_values) 
}

cdf_table <- function(freq_table, total = sum(freq_table$Freq)) {
  #function takes frequency table (from freq_table function) and produces a CDF table
  freq_table <- freq_table[order(freq_table$S), ]
  Prob <- freq_table$Freq / total
  CDF  <- cumsum(Prob)
  
  data.frame(
    S = freq_table$S,
    CDF = CDF
  )
}

#note n=18 took us rougly 7 hours to compute, similarly n=17 took few hours too
#we added progress bar so you can see how long they will take

Dn_table <- function(n, k) {
  #function ties in all previous functions to produce exact cdf tables
  combs <- combn(n, k)
  n_comb <- ncol(combs)
  
  values_list   <- vector("list", n_comb)
  counts_list <- vector("list", n_comb)
  
  pb <- txtProgressBar(min = 0, max = n_comb, style = 3)
  
  for (i in seq_len(n_comb)) {
    A <- combs[, i] #fix A and Ac
    Ac <- setdiff(1:n, A)
    
    upper_correct_A <- A - (1:k) # corrections 
    lower_correct_A <- (n - k) - upper_correct_A 
    upper_correct_Ac <- Ac - (1:(n - k)) 
    lower_correct_Ac <- k - upper_correct_Ac 
    
    freq_table_1 <- freq_table(k,     upper_correct_A, lower_correct_A, "bd")
    freq_table_2 <- freq_table(n - k, upper_correct_Ac, lower_correct_Ac, "ac")
    
    v1 <- as.numeric(names(freq_table_1)); n1 <- as.numeric(freq_table_1) 
    v2 <- as.numeric(names(freq_table_2)); n2 <- as.numeric(freq_table_2)
    
    Smat <- outer(v1, v2, "+") 
    Wmat <- outer(n1, n2, "*")
    conv_i <- rowsum(as.vector(Wmat), group = as.vector(Smat))
    values_list[[i]]   <- as.numeric(rownames(conv_i))
    counts_list[[i]] <- as.numeric(conv_i[, 1])
    
    setTxtProgressBar(pb, i)
  }
  
  close(pb)
  
  all_values   <- unlist(values_list)
  all_counts <- unlist(counts_list)
  
  agg <- rowsum(all_counts, group = all_values)
  result <- data.frame(S = as.numeric(rownames(agg)), Freq = as.numeric(agg[, 1]))
  result <- result[order(result$S), ] #orders frequency table in increasing order
  rownames(result) <- NULL
  cdf <- cdf_table(result)
  # remember to convert S to Dn: S = n(n-1)(n-2)(n-3)(n-4) * Dn
  denom <- n * (n - 1) * (n - 2) * (n - 3) * (n - 4)
  cdf$Dn <- cdf$S / denom
  cdf <- cdf[, c("S", "Dn", "CDF")]
  return(cdf)
}

Dn_quantile <- function(cdf_table, p) {
  stats::approx(x = cdf_table$CDF, y = cdf_table$Dn, xout = p)$y
}
