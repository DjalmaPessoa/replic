
# Estimativa da variancia do Indice de Gini
## Usa lodown para ler dados da PNAD

library(lodown)
pnad_cat <-
  get_catalog( "pnad" ,
               output_dir = file.path( path.expand( "~" ) , "PNAD" ) )

# 2011 only
pnad_cat <- subset( pnad_cat , year == 2011 )

# download the microdata to your local computer
pnad_cat <- lodown( "pnad" , pnad_cat )


## ler arquivo
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

# restringe população para calcular índice de Gini conforme IBGE:

sub_pnad_design <- 
  subset( 
    pnad_boot_design_post , 
    !is.na( v4720 ) & v4720 != 0 & v8005 >= 15
  )


# pesos de desenho, usado para calcular estimativas pontuais 
wsf <- weights(sub_pnad_design, "sampling")


# réplicas de pesos, usadas para estimar variância:

wwf <- weights(sub_pnad_design, "analysis")


## Função para calcular o índice de GINI:

ComputeGini <- function(x, w) {
  w <- w[order(x)]
  x <- x[order(x)]
  N <- sum(w)
  n <- length(x)
  big_t <- sum(x * w)
  r <- cumsum(w)
  Num <- sum((2 * r - 1) * x * w)
  Den <- N * big_t
  (Num/Den) - 1
}

# calcula estimativa do índice de Gini 

rval <- ComputeGini(sub_pnad_design$variables$v4720, wsf)

# calculas réplicas de estimativas do Gini, para estimar variância 

qq <- apply(wwf, 2, function(wi)ComputeGini(sub_pnad_design$variables$v4720,wi))

# Estimativa da variância a partir das réplicas de estimativas

variance <- survey::svrVar(qq, sub_pnad_design$scale, sub_pnad_design$rscales, 
                           mse = sub_pnad_design$mse, coef = rval)

# estimativa de cv:

100*sqrt(variance)/rval

#teste

