run_everything = function(data, glmFamily, c = 0.5, saveName){
  # My method
  n = base::NROW(data)
  m = n*c
  myMethod = model_space(data = data,
                         B = 100,
                         m = m,
                         glmFamily = glmFamily)
  fileName = paste0("D:/robustSelectData/myMethod",saveName,".rds")
  base::saveRDS(myMethod, file = fileName)

  #Backwards Robust Dev

  stepGLMRob = step_glmrob(data = data,
                           glmFamily = glmFamily)
  fileName = paste0("D:/robustSelectData/stepGLMRob",saveName,".rds")
  base::saveRDS(stepGLMRob, file = fileName)

  #Backwards AIC

  stepAIC = step_ic(data = data,
                           glmFamily = glmFamily)
  fileName = paste0("D:/robustSelectData/stepAIC",saveName,".rds")
  base::saveRDS(stepAIC, file = fileName)

  #Backwards BIC

  n = base::nrow(data)
  stepBIC = step_ic(data = data,
                    glmFamily = glmFamily,
                    k = log(n))
  fileName = paste0("D:/robustSelectData/stepBIC",saveName,".rds")
  base::saveRDS(stepBIC, file = fileName)

  #Bestglm AIC
  bestAIC = bestglm_ic(data = data,
                    glmFamily = glmFamily,
                    IC = "AIC")
  fileName = paste0("D:/robustSelectData/bestglmAIC",saveName,".rds")
  base::saveRDS(bestAIC, file = fileName)

  #BestGLM BIC

  bestBIC = bestglm_ic(data = data,
                    glmFamily = glmFamily,
                    IC = "BIC")
  fileName = paste0("D:/robustSelectData/bestglmBIC",saveName,".rds")
  base::saveRDS(bestBIC, file = fileName)

  return(list(myMethod = myMethod,
              stepGLMRob = stepGLMRob,
              stepAIC = stepAIC,
              stepBIC = stepBIC,
              bestAIC = bestAIC,
              bestBIC = bestBIC))
}
