## função para estimar erros padrão das estimativas dos parâmetros de posição b_i por 
## replicação de pesos.

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
}
  
# Exemplo de aplicação da função svrepMw aos dados PNAD2013:


### Leitura de dados da PNAD2013

# carrega library lodown

library(lodown)

# lista todos os arquivos de microdados da PNAD disponíveis
pnad_cat <-
  get_catalog( "pnad" ,
               output_dir = file.path( path.expand( "~" ) , "PNAD" ) )

# 2013 apenas
pnad_cat <- subset( pnad_cat , year == 2013)

# baixa os microdados para seu computador
pnad_cat <- lodown( "pnad" , pnad_cat )

# arquivo de microdados da pnad2013 (merge DOM e PES)

pnad_df <- readRDS( pnad_cat[ 1 , 'output_filename' ] )


## Carrega a library survey
library(survey)
options( survey.lonely.psu = "adjust" )


# cria réplicas de peso bootstrap

# réplicas de pesos antes de pós-estratificar
pnad_bootw <- bootweights(pnad_df$v4617, pnad_df$v4618, replicates = 80)

# cria desenho de replicação antes de pós-estratificação
pnad_boot_design <-
  survey::svrepdesign(
    weight = ~ pre_wgt ,
    repweights = pnad_bootw$repweights ,
    type = "bootstrap",
    combined.weights = FALSE ,
    scale = pnad_bootw$scale ,
    rscales = pnad_bootw$rscales ,
    data = pnad_df
  )

# totais de pós-estratificação
pop_types <-
  data.frame(
    v4609 = unique( pnad_df$v4609 ) ,
    Freq = unique( pnad_df$v4609 )
  )


# pós-estratificação do desenho de relicação
pnad_boot_design_post <- postStratify(design = pnad_boot_design ,
                                      strata = ~ v4609 ,
                                      population = pop_types )

# objeto de desenho para domicílios:

pnad_dom_design <-
  subset(pnad_boot_design_post , v0401==1)


# Aplica a função svrepMw

 result <- svrepMw (pnad_dom_design)
 result



