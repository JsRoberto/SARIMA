#--------------------------------------------------------------------------
# Séries Temporais - Função para madelagem SARIMA(p,d,q)(P,D,Q)[s]
#--------------------------------------------------------------------------

# Este arquivo modelSARIMA.R produz a função responsável pela modelagem
# SARIMA(p,d,q)(P,D,Q)[s] realizada em projSARIMA.R. Sua vantagem é que
# facilmente podemos incluir termos AR, MA, SAR e SMA em lags específicos. 

# Definindo a biblioteca
.libPaths(
      "C:/Users/JoséRoberto/AppData/Roaming/SPB_16.6/R/win-library/3.2")

# Definindo pacotes não padrões a serem utilizados 
suppressMessages(library(forecast))

# A função modelSARIMA() produz um modelo mediante forecast::Arima(). 
# -------------------- Descrição dos argumentos ---------------------------
# (1) "ts", objeto que representa a serie temporal;
# (2) "d", "p" e "q", parâmetros que indicam a ordem das diferenças 
# simples, dos temos AR e dos termos MA, resp.;
# (3) "D", "P" e "Q", parâmetros que indicam a ordem das diferenças 
# sazonais, dos temos SAR e dos termos SMA, resp.;
# (4) "addLags.'tipo'", vetor numérico que indica em quais lags serão adi-
# cionados termos de 'tipo' definido - AR, MA, SAR ou SMA.
modelSARIMA <- function(ts, p, d, q, P=0, D=0, Q=0, S=NA,
                        addLags.AR = 0, addLags.MA = 0,
                        addLags.SAR = 0, addLags.SMA = 0) {
      sign1 <- sign2 <- sign3 <- sign4 <- F
      extra <- if (d==0 & D==0) NA else numeric()
      cond.AR <- if (max(addLags.AR)>p) {sign1 <- T; addLags.AR} else p
      cond.MA <- if (max(addLags.MA)>q) {sign2 <- T; addLags.MA} else q
      cond.SAR <- if (max(addLags.SAR)>P) {sign3 <- T; addLags.SAR} else P
      cond.SMA <- if (max(addLags.SMA)>Q) {sign4 <- T; addLags.SMA} else Q
      fixed.VCT <- if (sign1 | sign2 | sign3 | sign4) {
            vctAR <- c(rep(NA,p),rep(0,ifelse(max(addLags.AR)-p < 0, 0,
                                              max(addLags.AR)-p)))
            if (addLags.AR[1]!=0) vctAR[addLags.AR] <- NA
            vctMA <- c(rep(NA,q),rep(0,ifelse(max(addLags.MA)-q < 0, 0,
                                              max(addLags.MA)-q)))
            if (addLags.MA[1]!=0) vctMA[addLags.MA] <- NA
            vctSAR <- c(rep(NA,P),rep(0,ifelse(max(addLags.SAR)-P < 0, 0,
                                               max(addLags.SAR)-P)))
            if (addLags.SAR[1]!=0) vctSAR[addLags.SAR] <- NA
            vctSMA <- c(rep(NA,Q),rep(0,ifelse(max(addLags.SMA)-Q < 0, 0,
                                               max(addLags.SMA)-Q)))
            if (addLags.SMA[1]!=0) vctSMA[addLags.SMA] <- NA
            fixed.vct <- c(vctAR, vctMA, vctSAR, vctSMA, extra)
      }
      # OBS.: função Arima(), em alguns casos (d=D=0 e p,P,q,Q != 0), está
      # produzindo erro. Uma solução possível é aplicar argumento method =
      # = "CSS", isso modifica ligeiramente os termos obtidos.
      SARIMA <- Arima(ts, order = c(cond.AR[length(cond.AR)], d, cond.MA[length(cond.MA)]),
                      seasonal = list(order = c(cond.SAR[length(cond.SAR)], D, cond.SMA[length(cond.SMA)]), period=S),
                      fixed = fixed.VCT)
      return(SARIMA)
}