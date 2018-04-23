
### Estimação original: Tabela 7.2 da dissertação 

#Limpar a memória do R
#rm(list=ls(all=TRUE))

#Leitura da PNAD2013

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

# carrega library survey

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

# objeto de desenho para domocílios: 

pnad_dom_design <- 
  subset(pnad_boot_design_post , v0401==1)



# pesos de desenho, usado para calcular estimativa pontual 
wsf <- weights(pnad_dom_design, "sampling")


# 80 réplicas de pesos, usadas para estimar variâncias:

wwf <- weights(pnad_dom_design, "analysis")

# verificar se peso de domicílios dados por v4611 coincide com
# wsf a menos de arredondamento:

sum(abs(pnad_dom_design$variables$v4611-wsf)>1)

## seleciona variáveis de interesse

dadosfies<- subset(pnad_dom_design$variables, select=c("uf","v0302","v8005","v0404","v0602","v9001","v4728","v4743"
                                         ,"v2103","v2105","v2107","v2109","v2113","v2115","v2117","v2121"
                                         ,"v2138","v2139","v4618","v4617","v4609", "v4611","one","region","pre_wgt")) 

## as variáveis "vvxxxx" não estão no dicionário da pnad e não foram lidas. Como obteve?

dadosfies <- transform (dadosfies, v2103 = as.numeric(v2103), v2105 = as.numeric(v2105),
v2107 = as.numeric(v2107), v2109 = as.numeric(v2109), v2113 = as.numeric(v2113), 
v2115 = as.numeric(v2115), v2117 = as.numeric(v2117), v2121 = as.numeric(v2121),
v2138 = as.numeric(v2138), v2139 = as.numeric(v2139))


#renomeando os itens segundo o FIES
#install.packages("plyr")
library(plyr)


data.FAO_Brasil<-rename(dadosfies, c("v2103"="WORRIED","v2105"="RUNOUT","v2107"="HEALTHY","v2109"="FEWFOOD"
                                     ,"v2113"="SKIPPED","v2115"='ATELESS',"v2117"="HUNGRY","v2121"="WHLDAY"
                                     ,"v2138"="Atitude","v2139"="Outra atitude"))



### Saving the FIES data and corresponding weights
XX.Brasil = data.FAO_Brasil[,9:16]


# calcula estimativa do item severity para nível WORRIED

library(RM.weights)

### Calculate raw scores (number of yes for each individual to the 8 questions)
rv.Brasil=rowSums(XX.Brasil)

### Number of items (questions) of the FIES
k = ncol(XX.Brasil)


# Section 2: Psychometric analysis ----------------------------------------------------------------
library(RM.weights)


#função para estimar os parâmetros do modelo Rasch de TRI.


# estima rval dever aproximadamente a da tabela 7.2

rval <- RM.w(XX.Brasil, wsf)$b["WORRIED"]

# calculas réplicas de estimativas do Gini, para estimar variância 

qq <- apply(wwf, 2, function(wi)RM.w(XX.Brasil, wi)$b["WORRIED"])

# Estimativa da variância a partir das réplicas de estimativas

variance <- survey::svrVar(qq, pnad_dom_design$scale, pnad_dom_design$rscales, 
                           mse = pnad_dom_design$mse, coef = rval)

# estimativa de cv:

100*sqrt(variance)/rval


# Compara valores

sqrt(variance)

RM.w(XX.Brasil, wsf)$se.b["WORRIED"]





