################################### Tese #############################################

##### Objetivo: Estimar o *erro padr?o* de alguma estimativa da minha disserta??o usando bootstrap

setwd("C:/Users/Etienne/Documents/DOUTORADO ENCE/TESE/Base de dados")
pnad2013DOM<-readRDS("pnad2013DOM.rds")

#Usando o survey para estimativas (levando em conta o peso amostral)

#install.packages("survey")
library(survey)
options(survey.lonely.psu = "adjust")

### escolhendo uma estimativa ###

# Base de dados FIES para o BRASIL (segundo PNAD 2013)
# Leva em considera??o somente os 8 itens referentes a domic?lios com maiores de 18 anos

dadosfies<- subset(pnad2013DOM, select=c("uf","v0302","v8005","v0404","v0602","v9001","v4728","v4743"
                                   ,"vv2103","vv2105","vv2107","vv2109","vv2113","vv2115","vv2117","vv2121"
                                   ,"v2138","v2139","v4618","v4617","v4609", "v4611","one","region","pre_wgt")) 

attach(dadosfies)

#install.packages("plyr")
library(plyr)

#renomeando os itens segundo o FIES
data.FAO_Brasil<-rename(dadosfies, c("vv2103"="WORRIED","vv2105"="RUNOUT","vv2107"="HEALTHY","vv2109"="FEWFOOD"
                                     ,"vv2113"="SKIPPED","vv2115"='ATELESS',"vv2117"="HUNGRY","vv2121"="WHLDAY"
                                     ,"v2138"="Atitude","v2139"="Outra atitude"))

head(pre_wgt) #peso b?sico do desenho (=v4610)
head(v4611) #peso do domic?lio segundo o dicion?rio
head(one) #peso 1 (teste)

#Fixando o banco de dados
attach(data.FAO_Brasil)
names(data.FAO_Brasil)

#Salvando o banco de dados do FIES e seus pesos correspondentes
XX = data.FAO_Brasil[,9:16]  #banco de dados
str(XX)
wt= v4611 #peso

# Calcula os escores por linha (N?mero de "sim" para as oito quest?es do FIES)
rv=rowSums(XX)

# N?mero de Itens (quest?es) do FIES
k = ncol(XX)

# An?lise Psicom?trica

#install.packages("RM.weights")
library(RM.weights) #pacote desenvolvido pela FAO
#?RM.w
#rr= RM.w(XX, wt) #fun??o constru?da o qual usa o Modelo Rasch para estima??o dos par?metros
  #Alguns resultados oriundos desta fun??o:
  # Gravidade do Item (par?metro de posi??o dos itens) dado por:
    #rr$b
  # Erro padr?o do Item dado por:
    #rr$se.b

### Estimando o par?metro de posi??o dos itens (rr$b) via bootstrap com base na observa??o 
    #de como a fun??o acima RM.w foi constru?da ### 


# cria r?plicas de peso bootstrap 

# r?plicas de pesos antes de p?s-estratificar
pnad_bootw <- bootweights(pnad2013DOM$v4617, pnad2013DOM$v4618, replicates = 80)

# cria desenho de replica??o antes da p?s estratifica??o
pnad_boot_design <-
  survey::svrepdesign(
    weight = ~ pre_wgt ,
    repweights = pnad_bootw$repweights ,
    type = "bootstrap",
    combined.weights = FALSE ,
    scale = pnad_bootw$scale ,
    rscales = pnad_bootw$rscales ,
    data = pnad2013DOM
  )

# totais de p?s estratifica??o
pop_types <- 
  data.frame( 
    v4609 = unique( pnad2013DOM$v4609 ) , 
    Freq = unique( pnad2013DOM$v4609 )
  )

# p?s estratifica??o do desenho de replica??o
pnad_boot_design_post <- postStratify(design = pnad_boot_design ,
                                      strata = ~ v4609 ,
                                      population = pop_types )

# pesos de desenho, usado para calcular estimativas pontuais 
wsf <- weights(pnad_boot_design_post, "sampling")

# r?plicas de pesos, usadas para estimar a vari?ncia:

wwf <- weights(pnad_boot_design_post, "analysis")

#Fun??o para calcular a estimativa de posi??o do itens

#install.packages("psychotools")
library(psychotools)

# Ao inv?s de usar wt, usar wsf
length(wt) #peso convencional (usado na disserta??o 2017)
length(wsf) #peso do desenho
length(wwf) #peso p/ c?lculo da vari?ncia

  #weighted likelihood estimation
  #?wle.rasch

wle.fit <- function(x, .data, .rawscores, .w) {
  .data <- XX
  .rawscores = rv
  ml.select <- which(!((.rawscores == 0) | (.rawscores == k)))
  .data = .data[ml.select,]
  .rawscores = .rawscores[ml.select]
  .w = wsf[ml.select]
  gamma.r.v = elementary_symmetric_functions(x)$"0"
  teil1a <- - as.matrix(.data) %*% x
  LL <- sum(.w * (teil1a - log(gamma.r.v[.rawscores+1])))
  return(-LL)
}

P.i.b <- function(x) {
  .data <- XX
  .rawscores = rv
  ml.select <- which(!((.rawscores == 0) | (.rawscores == k)))
  .data = .data[ml.select, ]
  .rawscores = .rawscores[ml.select]
  .w = wwf[ml.select]
  .w = (.w) * length(ml.select)/sum(.w)
  gamma.r.v = elementary_symmetric_functions(x)$"0"
  p = LL = matrix(NA, ncol=k, nrow=nrow(.data))
  for(i in 1:k){ 
    gamma.r.v1 = elementary_symmetric_functions(x[(-i)],ord=0)$"0"
    teil1a <- exp(- x[i]) * gamma.r.v1[.rawscores]
    p[,i] = teil1a / gamma.r.v[.rawscores+1] 
    LL[,i] <- (p[,i] * (1 - p[,i]))*.w
  } 
  LL = rowsum(LL, group = .rawscores)	
  return(LL)
}

# calcula estimativa pontual do par?metro de posi??o do dos itens

opt.w <- optim(seq(-3,3,length.out=k), wle.fit, method = "BFGS", hessian = T)
(b.w = opt.w$par) #usando wsf 
names(b.w) = colnames(XX)

p.i.b = colSums( P.i.b(b.w) )
(se.b.w = sqrt(1/p.i.b))  #erro padr?o usando wwf

    #d?vida: uso correto dos pesos? 



