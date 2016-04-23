# Lista de Afazeres #

## Modelo GLM COM-Poisson Convencional ##

 - [x] Acomodar _ofsset_ no modelo COM-Poisson `cmp`;
 - [ ] Adotar uma classe para o objeto retornado pela função `cmp` e,
   baseando-se na `MASS::glm.nb`, criar as funções método:
       - [x] `print.compois`:
       - [ ] `summary`
       - [ ] `print.summary`
       - [ ] `vcov`
       - [x] `logLik`
       - [ ] `coef`
       - [ ] `extractAIC`
       - [ ] `anova`
       - [ ] `print.anova`
       - [ ] `family`
 - [ ] Elaborar _vignettes_ com a análise dos _datasets_ contidos no
   pacote, exemplificando o uso do modelo COM-Poisson convencional.

## Modelo GLM COM-Poisson com Inflação de Zeros ##

 - [ ] Criar funções para estimação de um modelo COM-Poisson para dados
   com excesso de zeros.

## Modelo GLM COM-Poisson com Efeitos Aleatórios ##

 - [ ] Criar funções para estimação de um modelo COM-Poisson para acomodar
   efeitos aleatórios, inicialmente normais.
