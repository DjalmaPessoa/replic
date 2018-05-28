## Função para estimar se(b) e se(a) considerando desenho de replicação

svrepRMw <- function(design,...){
  
  require(RM.weights)
  dados <-design$variables
  
  dadosfies<- subset(dados,
                     select=c("v2103","v2105","v2107","v2109","v2113","v2115","v2117","v2121") )
  
  dadosfies <- transform (dadosfies, v2103 = as.numeric(v2103), v2105 = as.numeric(v2105),
                          v2107 = as.numeric(v2107), v2109 = as.numeric(v2109), v2113 = as.numeric(v2113), 
                          v2115 = as.numeric(v2115), v2117 = as.numeric(v2117), v2121 = as.numeric(v2121))
  
  dadosfies <- transform (dadosfies,
                          vv2103 = ifelse(is.na(v2103),0, ifelse(v2103==1,1,0)),
                          vv2105 = ifelse(is.na(v2105),0, ifelse(v2105==1,1,0)),
                          vv2107 = ifelse(is.na(v2107),0, ifelse(v2107==1,1,0)),
                          vv2109 = ifelse(is.na(v2109),0, ifelse(v2109==1,1,0)),
                          vv2113 = ifelse(is.na(v2113),0, ifelse(v2113==1,1,0)),
                          vv2115 = ifelse(is.na(v2115),0, ifelse(v2115==1,1,0)),
                          vv2117 = ifelse(is.na(v2117),0, ifelse(v2117==1,1,0)),
                          vv2121 = ifelse(is.na(v2121),0, ifelse(v2121==1,1,0))
  )
  
  data.FAO_Brasil<-plyr::rename(dadosfies, 
                                c("vv2103"="WORRIED","vv2105"="RUNOUT","vv2107"="HEALTHY","vv2109"="FEWFOOD"
                                  ,"vv2113"="SKIPPED","vv2115"='ATELESS',"vv2117"="HUNGRY","vv2121"="WHLDAY"))
  
  XX.Brasil = data.FAO_Brasil[,9:16]
  
  wsf <- weights(design, "sampling")
  
  
  full <- RM.w(XX.Brasil, wsf)
  
  rval_b <- full$b
  
  rval_a <- full$a
  
  wwf <- weights(design, "analysis")
  
  
  mat_a <- matrix(ncol = length(rval_a), nrow = ncol(wwf))
  mat_b <- matrix(ncol = length(rval_b), nrow = ncol(wwf))
  
  for(i in 1:ncol(wwf)){
    wi<- wwf[,i] 
    fulli <- RM.w(XX.Brasil, wi)
    mat_a[i,] <- fulli$a
    mat_b[i,] <- fulli$b
  }
  
  v_a <- svrVar(mat_a, design$scale, design$rscales, 
                mse = design$mse, coef = rval_a)
  
  v_b <- svrVar(mat_b, design$scale, design$rscales, 
                mse = design$mse, coef = rval_b)
  
  dimnames(v_b) <- list(names(rval_b),names(rval_b))
  
  se_a <- sqrt(diag(v_a))
  se_b <- sqrt(diag(v_b))
  names(se_b)<- names(rval_b)
  
  list(b = rval_b, covmat_b = v_b, e.p_b = se_b, a = rval_a, covmat_a = v_a, ep_a = se_a)
  
  
  
  ## Exemplo:
  
  result <- svrepRMw(pnad_dom_design)
  
  