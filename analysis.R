count_models = function(list, k = 1){
  list = base::lapply(X = list,
                      FUN = select_k,
                      k = k)
  lapply(X = list,
         FUN = base::NROW)
}

###############################################

IC = function(vector, data, glmFamily){
model = base::paste0("y~", vector[1]) %>%
  stats::as.formula()
glmFit = robustbase::glmrob(formula = model,
                             data = data,
                             family = glmFamily,
                             control = robustbase::glmrobMqle.control(maxit = 1000))
L1Norm = glmFit$coefficients %>%
  base::abs(x = .) %>%
  base::sum(x = .)

results = c(vector, L1Norm)

return(results)
}

###############################################

check_model_space = function(list, k = 1, data, glmFamily){
  list = base::lapply(X = list,
                FUN = select_k,
                k = k)
  output = base::lapply(X = list,
                  FUN = matrix_IC,
                  data = data,
                  glmFamily = glmFamily)

  for( i in 1:length(output)){
    output[[i]] = output[[i]] %>%
      base::t(.)
    output[[i]] = data.frame(output[[i]], stringsAsFactors = FALSE)
    base::colnames(output[[i]])[4] = "IC"
    rownames(output[[i]]) = paste0("Model ",1:NROW(output[[i]]))
  }

  output = base::lapply(X = output,
                        FUN = best_model)

  return(output)

}

###############################################

select_k = function(df, k){
  df$Count = df$Count %>%
    as.numeric()
  df %>%
    filter(Count >= k)
}

###############################################

matrix_IC = function(matrix, data, glmFamily){
  base::apply(X = matrix,
        MARGIN = 1,
        FUN = IC,
        data = data,
        glmFamily = glmFamily)
}


###############################################

best_model = function(matrix){
  lowest = matrix[,4] %>%
    as.numeric() %>%
    which.min()
  matrix$Best = FALSE
  matrix$Best[lowest] = TRUE
  return(matrix)
}

