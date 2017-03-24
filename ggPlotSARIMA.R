#--------------------------------------------------------------------------
# Séries Temporais - Funções para plotagem dos dados de modelagem_sarima.R
#--------------------------------------------------------------------------

# Este arquivo .R produz as funções necessários para plotar as modelagem 
# SARIMA(p,d,q)x(P,D,Q)[s] realizada em modelagem_sarima.R, cuja série sob
# análise é "usmelec", do pacote "fpp". 

# Definindo a biblioteca
.libPaths(
      "C:/Users/JoséRoberto/AppData/Roaming/SPB_16.6/R/win-library/3.2")

# Definindo pacotes não padrões a serem utilizados 
suppressMessages(library(ggplot2))
suppressMessages(library(fpp))

# A função ts2ggplot() plota exclusivamente os dados de "ts", que deve ser
# a série original ou sua transformação logaritmica (Box-Cox). 
ts2ggplot <- function(ts) {
      data <- data.frame(Tempo=time(ts) %>% as.numeric,
                         Dados=as.vector(ts))
      ggplot(data, aes(x=Tempo,y=Dados)) + 
            theme(plot.title = element_text(hjust=0.5)) +
            geom_line(size = 1.1, color = "cyan3") +
            scale_x_continuous(
                  breaks=data$Tempo[t <- seq(1,dim(data)[1],by=36)],
                  labels=paste0(rep("Jan/",length(t)),
                                round(data$Tempo[t]))) +
            labs(title = "Geração total de energia elétrica - EUA",
                 y = "Energia elétrica (terawatt-hora)",
                 x = "Tempo (mês/ano)")
}

# A função ZbarW() realiza um teste para avaliação da variância de "ts" no
# intuito de averiguar a necessidade da Tranformação de Box-Cox
ZbarW <- function(ts, N = 10) {
      Zbar <- W <- vector()
      for (i in seq(1,length(ts),N)) {
            tsWindow <- ts[seq(i, length.out = N)]
            Zbar <- c(Zbar, mean(tsWindow))
            W <- c(W, max(tsWindow)-min(tsWindow))
      }
      data <- data.frame(Zbar = Zbar, W = W)
      ggplot(data, aes(x = Zbar, y = W)) + 
            geom_point(size = 2, color = "blue3", alpha = 0.7, na.rm = T) +
            labs(title = "Teste para Avaliação da Variância \n",
                 x = expression(paste("Média (", bar(Z),")")),
                 y = "Amplitude (W)") +
            theme(plot.title = element_text(hjust=0.5)) +
            stat_smooth(na.rm = T)
}

# A função ggQQ() reproduz o Teste de Normalidade, que pretende verificar 
# se a distribuição do resíduo é gaussiana.
ggQQ <- function(model) {
      y <- quantile(model$residuals, c(0.25, 0.75), na.rm = T)
      x <- qnorm(c(0.25, 0.75))
      slope <- diff(y)/diff(x)
      int <- y[1L] - slope * x[1L]
      d <- data.frame(resids = as.vector(model$residuals))
      p <- ggplot(d,aes(sample=resids)) + 
            geom_abline(slope=slope,intercept=int,color="blue") +
            theme(plot.title = element_text(hjust=0.5)) +
            stat_qq(alpha=0.5, size=2, color="darkslategray") +
            labs(title="Teste de Normalidade do Resíduo",
                 y = "Resíduo (tWh)",
                 x = "Distribuição (%)")
      return(p)
}

# As funções ggACF() e ggPACF() produzem os gráficos acf (autocorrelation 
# function) e pacf (partial autocorrelation function) de séries temporais. 
# -------------------- Descrição dos argumentos ---------------------------
# (1) "ts", objeto da classe ts (ou o modelo forecast::Arima) cuja acf ou 
# pacf (residuais) será calculada;
# (2) "d" e "D", parâmetros que indicam quantas diferenças simples e sazo-
# nais, resp., foram realizadas na série original para obter "ts";
# (3) "log.scale", valor lógico que indica se "ts" passou por transformação
# logaritmica (Box-Cox) ou não;
# (4) "is.resd", valor lógico que indica se "ts" é resíduo da modelagem ou
# não; e
# (5) "lag.max", parâmetro que indica quantos lags a acf ou a pacf plotará.
ggACF <- function(ts, d=0, D=0, log.scale=T, is.resd=F, lag.max=50) {
      p <- if (is.Arima(ts)) ggAcf(ts$residuals,lag.max) else ggAcf(ts,lag.max)
      p <- p + geom_segment(lineend="butt",color="red",size=2) +
            theme(plot.title=element_text(hjust=0.5),
                  plot.subtitle=element_text(hjust=0.5))
      p <- if (log.scale==T & is.resd==F) {
            p + ggtitle(expression(
                  "Autocorrelação de "*Delta[12]^D*Delta^d*Y[t]),
                  subtitle = paste("D =", D, ", d =", d, "e Y = ln(Z)"))
      } else {
            idx <- stringi::stri_join(ts$arma[c(1,3)], ts$arma[c(6,7)],
                                      ts$arma[c(2,4)], sep=",")
            if (log.scale==T) {
                  p + ggtitle(expression(
                        "Autocorrelação dos Resíduos de "*tilde(Y)[t]),
                        subtitle = paste0("SARIMA(", idx[1],")(",idx[2],
                                         ")[",ts$arma[5],"]"))
            } else {
                  p + ggtitle(expression(
                        "Autocorrelação dos Resíduos de "*tilde(Z)[t]),
                        subtitle = paste0("SARIMA(", idx[1],")(",idx[2],
                                         ")[",ts$arma[5],"]"))
            }
      }
      return(p)
}

ggPACF <- function(ts, d=0, D=0, log.scale=T, is.resd=F, lag.max=50) {
      p <- if (is.Arima(ts)) ggAcf(ts$residuals,lag.max) else ggAcf(ts,lag.max)
      p <- p + geom_segment(lineend="butt",color="red",size=2) +
            theme(plot.title=element_text(hjust=0.5),
                  plot.subtitle=element_text(hjust=0.5))
      p <- if (log.scale==T & is.resd==F) {
            p + ggtitle(expression(
                  "Autocorrelação Parcial de "*Delta[12]^D*Delta^d*Y[t]),
                  subtitle = paste("D =", D, ", d =", d, "e Y = ln(Z)"))
      } else {
            idx <- stringi::stri_join(ts$arma[c(1,3)], ts$arma[c(6,7)],
                                      ts$arma[c(2,4)], sep=",")
            if (log.scale==T) {
                  p + ggtitle(expression(
                        "Autocorrelação Parcial dos Resíduos de "*tilde(Y)[t]),
                        subtitle = paste0("SARIMA(", idx[1],")(",idx[2],
                                         ")[",ts$arma[5],"]"))
            } else {
                  p + ggtitle(expression(
                        "Autocorrelação Parcial dos Resíduos de "*tilde(Z)[t]),
                        subtitle = paste0("SARIMA(", idx[1],")(",idx[2],
                                         ")[",ts$arma[5],"]"))
            }
      }
      return(p)
}

# A função ggPlot.results() produz três funções de plotagem dos resultados
# da modelagem sarima de projSARIMA.R.
ggPlot.results <- function(data, fitted, forcst, model) {
      # ----------------- Descrição dos argumentos ------------------------
      # (1) "data", objeto da classe ts que representa a série original;
      # (2) "fitted", objeto da classe ts que representa a série estimada
      # pelo modelo do argumento "model";
      # (3) "forcst", objeto da classe ts que representa a série de previ-
      # são com base no modelo do argumento "model"; e
      # (4) "model", o modelo gerado na aplicação do algoritmo de 
      # projSARIMA.R.
      graphset <<- ts.union(Tempo=time(data),Dados=data, 
                           Ajuste=fitted,Previsão=forcst) %>% as.data.frame
      graphset.all <- reshape2::melt(graphset, id = "Tempo", 
                                     variable.name = "Tipo",
                                     value.name = "Valor")
      graphset.resd <- data.frame(Tempo = as.vector(time(model$residuals)),
                                  Residuo = as.vector(model$residuals))

      MSE <- mean((graphset$Dados - graphset$Prev)^2, na.rm = TRUE)
      
      # A função ggplotALL() realiza plotagem geral e completa dos dados 
      ggplotALL <<- function() {
            ggplot(graphset.all, aes(x=Tempo, y=Valor, color=Tipo)) +
                  theme(plot.title=element_text(hjust=0.5)) +
                  geom_line(size = 1.1, na.rm = TRUE) +
                  geom_vline(xintercept=max(graphset.all$Tempo)-numprev/12,
                             lty=2) +
                  scale_color_manual(
                        values=c("cyan3","darkorange","darkslategray4")) +
                  scale_x_continuous(
                        breaks=graphset$Tempo[
                              t <- seq(1,dim(graphset)[1],36)],
                        labels=paste0(rep("Jan/",length(t)),
                                      round(graphset$Tempo[t]))) +
                  labs(title="Geração de Energia Elétrica nos EUA",
                       y = "Energia Elétrica (terawatt-hora)",
                       x = "Tempo (mês/ano)")
      }
      # A função ggplotPrev() realiza plotagem dos dados de previsão, bem
      # como seu intervalo de confiança
      ggplotPrev <<- function() {
            ggplot(graphset, aes(x=Tempo)) +
                  theme(plot.title=element_text(hjust=0.5)) +
                  geom_line(aes(y=Dados, colour="Dados"), size=1.1, na.rm = T) +
                  geom_line(aes(y=Previsão, colour="Previsão"), size=1.1, na.rm=T) +
                  geom_ribbon(data=graphset,
                              aes(ymin=Previsão-1.96*sqrt(MSE), ymax=Previsão+1.96*sqrt(MSE)),
                              alpha=0.2) +
                  scale_color_manual("Tipo", breaks = c("Dados","Previsão"),
                                     values = c("cyan3", "darkslategray4")) +
                  scale_x_continuous(
                        breaks=graphset$Tempo[seq(1,dim(graphset)[1],by=6)],
                        labels=paste0(rep(c("Jan/","Jul/"),length(seq(1,dim(graphset)[1],by=12))),
                                      round(graphset$Tempo[seq(1,dim(graphset)[1],by=6)])),
                        limits = c(max(graphset$Tempo)-numprev/12,
                                   max(graphset$Tempo))) +
                  labs(title = "Geração de Energia Elétrica nos EUA",
                       x = "Tempo (mês/ano)", y = "Energia Elétrica (terawatt-hora)")
      }
      # A função ggplotResd() para plotagem dos resíduos da modelagem
      ggplotResd <<- function() {
            ggplot(graphset.resd, aes(x=Tempo, y=Residuo)) +
                  scale_color_manual(values = c("darkslategray")) +
                  geom_line(size = 1.1, na.rm = TRUE) +
                  geom_hline(yintercept=mean(model$residuals), size=1.2,
                             lty="dashed", color="blue3") +
                  theme(plot.title=element_text(hjust=0.5),
                        plot.subtitle=element_text(hjust=0.5)) +
                  labs(title="Geração de Energia Elétrica nos EUA",
                       y="Energia elétrica (terawatt-hora)",
                       subtitle="Resíduos da modelagem",
                       x="Tempo (mês/ano)")
      }
      
}




