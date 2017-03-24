#--------------------------------------------------------------------------
# Séries Temporais - Modelagem SARIMA(p,d,q)x(P,D,Q)s
#--------------------------------------------------------------------------

# Definindo a biblioteca
.libPaths(
      "C:/Users/JoséRoberto/AppData/Roaming/SPB_16.6/R/win-library/3.2")

# Definindo pacotes não padrões a serem utilizados 
suppressMessages(library(gridExtra))

#--------------------------------------------------------------------------
#                       Fase de Identificação - 1ª Etapa
#--------------------------------------------------------------------------

# Estabelemento das funções de plotagem armazenadas em ggPlotSARIMA.R
source("ggPlotSARIMA.R", encoding = "UTF-8")

# Plotagem da série original
ts2ggplot(usmelec)

# Teste para avaliação da variância para Tranformação Box-Cox
ZbarW(usmelec)

# Definindo as séries de modelagem e de previsão, bem como o número de 
# meses a serem previstos
numprev <- 24
logUsmelec <- ts(log(usmelec), start = start(usmelec), frequency = 12)
modlogUsmelec <- logUsmelec[seq(length.out = length(logUsmelec)-numprev)]
prevlogUsmelec <- logUsmelec[-seq(length.out = length(logUsmelec)-numprev)]

# Verificação de acf para determinar parâmetro d
grid.arrange(ggACF(modlogUsmelec), nrow=3,
             ggACF(diff(modlogUsmelec),d=1) + labs(title=""),
             ggACF(diff(modlogUsmelec,1,2),d=2) + labs(title=""))
# Conclusão: d = 1

# Verificação de acf sobre o modelo d=1 para determinar parâmetro D
grid.arrange(ggACF(diff(modlogUsmelec),d=1), nrow=3,
             ggACF(diff(diff(modlogUsmelec),12),d=1,D=1)+labs(title=""),
             ggACF(diff(diff(modlogUsmelec),12,2),d=1,D=2)+labs(title=""))
# Conclusão: D = 1

# Série obtida é W[t]=(1-B)*(1-B^s)*Y[t] -- d=1 e D=1
Wt <- diff(diff(modlogUsmelec),12)

# Determinação dos valores de p, q, P e Q mediante plotagem da acf e da 
# pacf para a série previamente obtida
grid.arrange(ggACF(Wt,d=1,D=1), ggPACF(Wt,d=1,D=1), nrow=2)

# Visualização dos lags sazonais {12,24,36,...} da série W[t]
acfSlags <- acf(Wt, 240, plot = FALSE)
acfSlags[seq(from=12,to=240,by=12),]

pacfSlags <- pacf(Wt, 240, plot = FALSE)
pacfSlags[seq(from=12,to=240,by=12),]

#--------------------------------------------------------------------------
#                 Fase de Estimação e Diagnóstico - 2ª Etapa
#--------------------------------------------------------------------------

# Estabelemento da função de modelagem armazenada em modelSARIMA.R
source("modelSARIMA.R", encoding = "UTF-8")

# O 1º modelo tentativo é SARIMA(0,1,3)(0,1,1)[12]
SARIMA <- modelSARIMA(modlogUsmelec,0,1,3,0,1,1,12)

# Testes de Diagnóstico do 1º modelo tentativo
# Testes de 1ª ordem para os resíduos
mean(SARIMA$residuals)
ZbarW(SARIMA$residuals)

# Testes de 2ª ordem para os resíduos: autocorrelação residual
grid.arrange(ggACF(SARIMA,d=1,D=1,is.resd=T), nrow=2,
             ggPACF(SARIMA,d=1,D=1,is.resd=T))

#                                   ***
# O 2º modelo tentativo, com acréscimo de termos MA nos lags 10 e 15, é
# SARIMA(0,1,15)x(0,1,1)[12]
SARIMA <- modelSARIMA(modlogUsmelec,0,1,3,0,1,1,12,addLags.MA = c(10,15))

# Testes de Diagnóstico do 2º modelo tentativo
# Testes de 1ª ordem para os resíduos
mean(SARIMA$residuals)
ZbarW(SARIMA$residuals)

# Testes de 2ª ordem para os resíduos: autocorrelação residual
grid.arrange(ggACF(SARIMA$residuals,d=1,D=1,is.resd=T), nrow=2,
             ggPACF(SARIMA$residuals,d=1,D=1,is.resd=T))

# Teste de Normalidade para verificar distribuição gaussinana do resíduo
ggQQ(SARIMA)

#--------------------------------------------------------------------------
#                       Fase de Previsão - 3ª Etapa
#--------------------------------------------------------------------------

# Previsão com base no modelo final e preparação dos dados para plotagem
zt <- function(ts, model, logArg=T) { 
      if (logArg==T) exp(ts+var(model$residuals)/2) else ts
}
data <- zt(logUsmelec,SARIMA) %>% ts(start=start(logUsmelec),frequency=12)
fitted <- zt(modlogUsmelec - SARIMA$residuals, SARIMA) %>% 
      ts(start = start(data), frequency = 12)
forcst <- forecast::forecast(SARIMA,h=numprev)$mean %>% zt(SARIMA) %>%
      ts(end = end(data), frequency = 12)

# Estabelecimento do resultados gráficos
ggPlot.results(data, fitted, forcst, model = SARIMA)

# Testes para avaliação do modelo final
# 1) Variância do resíduo
ZbarW(SARIMA$residuals)
# 2) MSE = mean square error
MSE <- mean((graphset$Dados - graphset$Previsão)^2, na.rm = TRUE)
# 3) MAPE = mean absolut percentage error
MAPE <- 100*mean(abs(graphset$Dados-graphset$Previsão)/abs(graphset$Dados),
                 na.rm = TRUE)

#                                   ***
# Plotagem dos gráficos das séries analisadas
# Plotagem dos dados originais, dos ajustados e das previsões
ggplotALL()

# Plotagem dos dados de previsão, juntamente com seu intervalo de confiança
ggplotPrev()

# Plotagem dos resíduos da modelagem
ggplotResd()
