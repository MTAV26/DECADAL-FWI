
library(s2dv)
library(multiApply)

CorrEno <- function(exp, obs, time_dim = 'year', member_dim = NULL, method = 'pearson', alpha = 0.05, test.type = 'two-sided', pval = FALSE, handle.na = 'return.na', ncores = 1){
  
  ## Ensemble mean
  if (!is.null(member_dim)){
    exp <- multiApply::Apply(data = exp, target_dims = member_dim, fun = mean, na.rm = FALSE, ncores = ncores)$output1
  }
  
  ## Correlation coefficient and its significance
  output <- multiApply::Apply(data = list(exp = exp, obs = obs), target_dims = time_dim, 
                              fun = .CorrEno, time_dim = time_dim, method = method, 
                              alpha = alpha, test.type = test.type, pval = pval, 
                              handle.na = handle.na, ncores = ncores)
  return(output)
}

.CorrEno <- function(exp, obs, time_dim, method, alpha, test.type, pval, handle.na){
  
  .correlation_eno <- function(exp, obs, time_dim, method, alpha, test.type, pval){
    
    cor <- NULL
    cor$r = cor(x = exp, y = obs, method = method) # Correlation coefficient
    
    n_eff = s2dv::Eno(data = obs, time_dim = time_dim, na.action = na.pass, ncores = 1)
    
    if (test.type == 'one-sided'){
      
      t_alpha_n2 = qt(p=alpha, df = n_eff-2, lower.tail = FALSE)
      t = cor$r * sqrt(n_eff-2) / sqrt(1-cor$r^2)
      
      if (anyNA(c(t,t_alpha_n2)) == FALSE & t >= t_alpha_n2 & cor$r > 0){
        cor$sign = TRUE
      } else {
        cor$sign = FALSE
      }
      
      if (isTRUE(pval)){
        cor$pval <- pt(q = t, df = n_eff-2, lower.tail = FALSE)
      }
      
    } else if (test.type == 'two-sided'){
      
      t_alpha2_n2 = qt(p=alpha/2, df = n_eff-2, lower.tail = FALSE)
      t = abs(cor$r) * sqrt(n_eff-2) / sqrt(1-cor$r^2)
      
      if (anyNA(c(t,t_alpha2_n2)) == FALSE & t >= t_alpha2_n2){
        cor$sign = TRUE
      } else {
        cor$sign = FALSE
      }
      
      # cor$n_eff <- n_eff
      
      if (isTRUE(pval)){
        cor$pval <- 2 * pt(q = t, df = n_eff-2, lower.tail = FALSE)
      }
      
    } else {stop('test.type not supported')}
    
    return(cor)
  }
  
  #==================================================
  
  if (anyNA(exp) | anyNA(obs)) { ## There are NAs 
    if (handle.na == 'only.complete.pairs') {
      nna <- is.na(exp) | is.na(obs) # A vector of T/F
      if (all(nna)) {
        # stop("There is no complete set of forecasts and observations.")
        output <- list(r = NA, sign = NA)
        if (pval) {
          output <- c(output, list(p.val = NA))
        }
      } else {
        # Remove the incomplete set
        exp <- exp[!nna]
        obs <- obs[!nna]
        output <- .correlation_eno(exp, obs, time_dim, method, alpha, test.type, pval)
      }
    } else if (handle.na == 'return.na') {
      # Data contain NA, return NAs directly without passing to .correlation_eno
      output <- list(r = NA, sign = NA)
      if (pval) {
        output <- c(output, list(p.val = NA))
      }
    }
    
  } else { ## There is no NA  
    output <- .correlation_eno(exp, obs, time_dim, method, alpha, test.type, pval)
  }
  
  return(output)
}
