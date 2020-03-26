

mosaic_plot = function(list, measure = 'wald'){
  name = base::paste0(measure, "SpaceSVal")
  which = which(names(list) == name)

  S = list[[which]]

  p = base::NCOL(S)
  B = base::NROW(S)

  hMatrix = base::matrix(data = NA_real_,
                         nrow = p,
                         ncol = p)

  base::colnames(hMatrix) = base::colnames(S)
  base::rownames(hMatrix) = base::colnames(S)

  sAbs = base::t(base::apply(X = S,
                             MARGIN = 1,
                             FUN = base::abs))
  sRank = p + 1 - base::t(base::apply(X = sAbs,
                                      MARGIN = 1,
                                      FUN = rank))

  for(i in 1:p){
    for(j in 1:p){
      hMatrix[i,j] = base::sum(sRank[,i] < sRank[,j])
    }
  }

  hMatrix = hMatrix/B

  df = reshape::melt(hMatrix)
  base::colnames(df) = c("Variable Y", "Variable X", "Probability")

  order = base::order(base::abs(base::apply(X = S,
                                            MARGIN = 2,
                                            FUN = mean)))

  df$`Variable Y` = base::factor(df$`Variable Y`,
                                 levels = colnames(hMatrix)[order])
  df$`Variable X` = base::factor(df$`Variable X`,
                                 levels = base::rev(base::colnames(hMatrix)[order]))

  return(list(hMatrix = hMatrix,
              plotData = df))
}

###############################################


variable_inclusion_plot = function(list,
                                    measure = 'wald'){
  name = base::paste0(measure, "SpaceSVal")
  which = which(names(list) == name)

  S = list[[which]]

  p = base::NCOL(S)
  B = base::NROW(S)

  gMatrix = base::matrix(data = NA_real_,
                   nrow = p,
                   ncol = p + 1)

  colnames(gMatrix) = 0:p
  rownames(gMatrix) = colnames(S)

  sAbs = base::t(base::apply(X = S,
                             MARGIN = 1,
                             FUN = abs))
  sRank = p + 1 - t(apply(X = sAbs, MARGIN = 1, FUN = rank))

  for(i in 1:(p + 1)){
    for(j in 1:p){
      gMatrix[j, i] = base::sum(sRank[,j] <= (i - 1))
    }
  }

  gMatrix = gMatrix/B

  base::colnames(gMatrix) = base::paste0("Dim",0:p)

  df = reshape::melt(gMatrix)
  base::colnames(df) = c("Variable", "Dimension", "Probability")
  df$Dimension =  factor(df$Dimension, levels = colnames(gMatrix))

  return(list(gMatrix = gMatrix,
              plotData = df))
}
