<img src = "./inst/img/tccPackage.svg" width=150px align="right" display="block">

## Extensões e Aplicações do Modelo COM-Poisson ##

> [Eduardo E. R. Junior](https://gitlab.c3sl.ufpr.br/u/eerj12) - Curso
de Estatística - Universidade Federal do Paraná  
> [Walmes M. Zeviani](https://gitlab.c3sl.ufpr.br/u/walmes) -
Departamento de Estatística - Universidade Federal do Paraná

Pacote R para ajuste e diagnóstico de modelos _Conway-Maxwell-Poisson
(COM-Poisson)_ para dados de contagem sub ou superdispersos. Modelo de
regressão tradicional, modelos para dados inflacionados de zeros e
modelos mistos são considerados. Este pacote faz parte de meu trabalho
de conclusão de curso e é desenvolvido em conjunto com o repositório
[tccDocument], que possui uma completa revisão e definição dos modelos
considerados no pacote.

***

### Introdução ###

Dados de contagem não raramente são analisados utilizando a distribuição
Poisson como distribuição associada a um modelo linear generalizado, no
caso de modelos de regressão. A distribuição Poisson possui somente um
parâmetro, que representa a média e variância. Essa propriedade, chamada
de equidispersão, é particularidade do modelo Poisson que pode não ser
adequada a diversas situações. Quando o modelo Poisson é aplicado sob
negligência desta suposição, o modelo apresenta erros padrões
inconsistentes para as estimativas dos parâmetros e, por consequência,
para toda função desses parâmetros.

A distribuição COM-Poisson é uma distribuição que contempla os casos de
sub e superdispersão devido a adição de mais um parâmetro. O pacote
`tccPackage` se propõem no ajuste e diagnóstico de modelos de regressão
que utilizam essa distribuição. Dados de contagem inflacionados de
zeros, e com estrutura para inclusão de efeitos mistos são considerados.

### Download e Instalação ###

O pacote é desenvolvido sob versionamento [Git] e mantido no serviço de
hospedagem [GitLab do C3SL]. Você poderá fazer o _download_ e instalação
do pacote de duas formas:

1. Usando o pacote `devtools` (disponível no [CRAN]), diretamente pelo
   endereço do repositório do GitLab.
```r
library(devtools)
install_git("https://gitlab.c3sl.ufpr.br/eerj12/tccPackage.git")
```

2. Baixando o código-fonte comprimido e instalando de forma usual. Aqui
   há uma diferença para o sistema operacional utilizado.
   - **Linux/Mac**
   Baixe o arquivo [tccPackage_0.0.1.tar.gz] e instale-o em uma sessão R
```r
install.packages("tccPackage_0.0.1.tar.gz", repos = NULL)
```

   - **Windows**
   Baixe o arquivo [tccPackage_0.0.1.zip] e instale-o em uma sessão R
```r
install.packages("tccPackage_0.0.1.zip", repos = NULL)
```