model_space = function(data, B, m, nStrata = 5, glmFamily = binomial(link = "logit"), k = 1, resid = "pearson", coef = TRUE, wald = TRUE, dev = FALSE){

  # Parameters
  n = NROW(data)
  p = NCOL(data) - 1
  nMethods = coef + wald + dev

  # Variables
  bootstraps = base::matrix(data = NA_integer_,
                            nrow = B,
                            ncol = m)

  bootstrapModels = base::rep(x = list(NA),
                        times = B)

  varNames = base::colnames(x = x)

  coefValues = matrix(data = NA_integer_,
                      nrow = B,
                      ncol = p)
  colnames(coefValues) = varNames

  waldValues = matrix(data = NA_integer_,
                      nrow = B,
                      ncol = p)
  colnames(waldValues) = varNames

  devValues = matrix(data = NA_integer_,
                      nrow = B,
                      ncol = p)
  colnames(devValues) = varNames

  output = base::list()

  sVals = base::list()

  # Initial Fit

  glmFull = robustbase::glmrob(formula = y~.,
                               data = data,
                               family = glmFamily,
                               control = robustbase::glmrobMqle.control(maxit = 1000))
  glmFullResid = stats::residuals(glmFull,
                                  type = resid)
  residOrder = rank(x = glmFullResid)
  residStrata = ceiling(x = residOrder * nStrata / n)

  # Generate weights
  for(b in 1:B){
    for(k in 1:nStrata){
      bootstraps[b, ((k - 1) * m / nStrata + 1 : (m / nStrata))] = base::sample(x = which(residStrata == k),
                                                                                size = m / nStrata)
    }
  }

  # Fit Models

  sValues = plyr::alply(.data = bootstraps,
                        .margins = 1,
                        .fun = s_values,
                        data = data,
                        n= n,
                        p = p,
                        glmFamily = glmFamily,
                        coef = coef,
                        wald = wald,
                        dev = dev,
                        .parallel = FALSE)

  if(coef == TRUE){
    coefPath = lapply(X = sValues,
                      FUN = solution_path,
                      row = 1)

    coefSpace = reduced_space(list = coefPath, vector = varNames)
    output$coefSpace = coefSpace
    sVals$coef = base::matrix(data = NA_integer_,
                              nrow = B,
                              ncol = p)
    for (i in 1:B){
      sVals$coef[i,] = sValues[[i]][1,]
    }
    base::colnames(sVals$coef) = varNames
    base::rownames(sVals$coef) = base::paste0("Bootstrap",1:B)

  }

  if(wald == TRUE){
    waldPath = lapply(X = sValues,
                      FUN = solution_path,
                      row = 2)

    waldSpace = reduced_space(list = waldPath, vector = varNames)
    output$waldSpace = waldSpace
    sVals$wald = base::matrix(data = NA_integer_,
                              nrow = B,
                              ncol = p)
    for (i in 1:B){
      sVals$wald[i,] = sValues[[i]][2,]
    }
    base::colnames(sVals$wald) = varNames
    base::rownames(sVals$wald) = base::paste0("Bootstrap",1:B)
  }

  if(dev == TRUE){
    devPath = lapply(X = sValues,
                      FUN = solution_path,
                      row = 3)

    devSpace = reduced_space(list = devPath, vector = varNames)
    output$devSpace = devSpace
    sVals$dev = base::matrix(data = NA_integer_,
                              nrow = B,
                              ncol = p)
    for (i in 1:B){
      sVals$dev[i,] = sValues[[i]][3,]
    }
    base::colnames(sVals$dev) = varNames
    base::rownames(sVals$dev) = base::paste0("Bootstrap",1:B)
  }

  countModels = count_models(output,
                             k = k)

  output = check_model_space(output,
                             k = k,
                             data = data,
                             glmFamily = glmFamily)

  additional = countModels %>%
    base::names()

  addition = base::paste0(additional, "Count")

  for(i in 1:nMethods){
    output[[nMethods + i]] = countModels[[i]]
  }

  names(output)[nMethods + 1:nMethods] = addition

  for(i in 1:nMethods){
    output[[2*nMethods + i]] = sVals[[i]]
  }

  names(output)[2*nMethods + 1:nMethods] = base::paste0(additional, "SVal")

  return(output)

}




#####################################################################

s_values = function(weights, data, n, p, glmFamily, coef, wald, dev = FALSE, anovaTest = "QD"){
  glmBoot = robustbase::glmrob(formula = y~.,
                               data = data,
                               subset = weights,
                               family = glmFamily,
                               control = robustbase::glmrobMqle.control(maxit = 1000))

  glmBootSummary = summary(glmBoot)

  sVal = base::matrix(data = NA_integer_,
                      nrow = 3,
                      ncol = p)

  colnames(sVal) = base::names(glmBoot$coefficients[-1])
  rownames(sVal) = c("Coef", "Wald", "Dev")

  if(coef == TRUE){
    sVal[1,] = glmBootSummary$coefficients[-1, 1]
  }

  if(wald == TRUE){
    sVal[2,] = glmBootSummary$coefficients[-1, 3]
  }

  if(dev == TRUE){
    for(i in 1:p){

      variables = utils::head(x = names(data),
                        n = p)

      formulaDev = as.formula(paste( "y ~", paste(variables[-i], collapse = "+")))

      glmDev = robustbase::glmrob(formula = formulaDev,
                                    data = data,
                                    subset = weights,
                                    family = glmFamily,
                                    control = robustbase::glmrobMqle.control(maxit = 1000))
      anovaDev = anova(glmDev,
                         glmBoot,
                         test = anovaTest)

      sVal[3,] = anovaDev$Test.Stat[2]

    }
  }

  return(sVal)

}

#####################################################################

solution_path = function(matrix, row){

  vector = matrix[row, ]

  p = base::length(vector)

  absVector = base::abs(vector)
  ranks =  p + 1 - base::rank(absVector)

  var = base::names(ranks)
  n = base::length(ranks)

  matrix = base::matrix(data = NA_integer_,
                        nrow = (n + 1),
                        ncol = 2)
  df = base::data.frame(matrix)
  base::names(df) <- c("Variables","Dimension")

  for (i in 0:n){
    whichVar = ranks <= i
    string = base::paste(c("1", var[whichVar]), sep = "", collapse = "+")
    df$Variables[i + 1] = string
    df$Dimension[i + 1] = i
  }

  return(df)

}

#####################################################################

reduced_space = function(list, vector){

  allModels = list %>%
    plyr::ldply(., data.frame) %>%
    plyr::ddply(., c("Variables", "Dimension"), summarise,
          Count = length(Variables)) %>%
    dplyr::arrange(Dimension) %>%
    dplyr::mutate(Variables = paste(Variables, sep = "+"))

  return(allModels)

}

#####################################################################

stepwise_robust = function(data){

}

