step_glmrob = function(data, glmFamily, anovaTest = "QDapprox", pVal = 0.05){

  glmFull = robustbase::glmrob(y~.,
                       data = data,
                       family = glmFamily)
  glmSummary = list()

  variables = names(data)
  variables = head(x = variables,
                   n = p)

  vars = variables

  glmSteps = matrix(data = NA,
                    nrow = 3,
                    ncol = p)

  base::row.names(glmSteps) = c("Initial Model", "Remove variable", "p-value")

  continue = TRUE

  iter = 1

  while(continue == TRUE){

    pval = base::rep(x = 0,
                     times = length(vars))
    model = base::rep(x = NA,
                      times = length(vars))

    for(k in 1:length(vars)){

      formula = stats::as.formula(base::paste( "y ~", paste(vars[-k], collapse = "+")))

      glmFit = robustbase::glmrob(formula,
                                      data = data,
                                      family = glmFamily,
                                      control = robustbase::glmrobMqle.control(maxit = 1000))

      anovaFit = anova(glmFit,
                        glmFull,
                        test = anovaTest)

      model[k] = as.character(x = formula)
      pval[k] = anovaFit$`Pr(>chisq)`[2]

    }

    check = sum(pval >= pVal)
    if(check == 0){
      continue = FALSE
    }

    remove = base::which.max(x = pval)

    glmSteps[,iter] = c(base::paste( "y ~", base::paste(vars, collapse = "+")), vars[remove] , pval[remove])

    vars = vars[-remove]

    iter = iter + 1
  }
  glmSteps = glmSteps[,-c(iter, iter-1)]
  glmSteps = janitor::remove_empty(glmSteps, which = "cols")
  base::colnames(glmSteps) = base::paste0("Step", 1:base::NCOL(glmSteps))
  glmSteps = glmSteps %>%
    data.frame(., stringsAsFactors = FALSE)

  return(glmSteps)

}


############################################################################

step_ic = function(data, glmFamily, k = 2){
  glmFull = stats::glm(y~.,
                       data = data,
                       family = glmFamily)
  step = MASS::stepAIC(object = glmFull,
                       direction = "backward",
                       trace = 0,
                       k = k)
  return(step)
}

############################################################################

bestglm_ic = function(data, glmFamily, IC = "AIC"){
  search = bestglm::bestglm(Xy = data,
                            family = glmFamily,
                            IC = IC)
  return(search)
}

