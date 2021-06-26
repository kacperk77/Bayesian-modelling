setwd('C:\\Users\\kacpe\\Desktop\\ekonometria bayesowska')

library(tidyverse)
library(corrplot)
library(caret)
library(mvtnorm)
library(manipulate)

remove(list = ls())



#1. Regresja wieloraka (OLS) i statystyki dostateczne 
dane_aktualne <- read.csv('dane_aktualne.csv', sep = ';')

colnames(dane_aktualne)[1] <- 'Stacja'
colnames(dane_aktualne)[2] <- 'MaxTemp'
colnames(dane_aktualne)[3] <- 'MinTemp'
colnames(dane_aktualne)[4] <- 'Rainfall'
colnames(dane_aktualne)[5] <- 'Snowfall'


dane_aktualne <- dane_aktualne %>% group_by(Stacja) %>%
  mutate(nextday_mintemp = lead(MinTemp)) 

dane_aktualne <- dane_aktualne %>% filter(nextday_mintemp != '-')
dane_aktualne <- dane_aktualne %>% filter(Snowfall != '-')

dane_aktualne$MinTemp <- as.numeric(dane_aktualne$MinTemp)
dane_aktualne$Snowfall <- as.numeric(dane_aktualne$Snowfall)
dane_aktualne$Rainfall <- as.numeric(dane_aktualne$Rainfall)
dane_aktualne$MaxTemp <- as.numeric(dane_aktualne$MaxTemp)
dane_aktualne$Data <- as.Date(dane_aktualne$Data)

dane_aktualne <- dane_aktualne %>% select(Data, nextday_mintemp, MinTemp, Rainfall, Snowfall, Stacja)


#losowanie danych 
set.seed(1234)
sample_data <- dane_aktualne %>% group_by(Data) %>% sample_n(1)

#cor <- cor(sample_data[,
 #                 c("Rainfall" ,       "MaxTemp" ,        "MinTemp"  ,          
 #                   "Snowfall")])

#corrplot.mixed(cor,  tl.col="black", tl.pos = "lt")


#regresja liniowa
OLS_results <- lm(nextday_mintemp~MinTemp+Rainfall+Snowfall, data = sample_data)
summary(OLS_results)


y <- as.matrix(sample_data[, c('nextday_mintemp')])
N.data <- length(y)
X <- cbind(as.matrix(rep(1, N.data)), 
           as.matrix(sample_data[, c('MinTemp','Rainfall', 'Snowfall')]))
Beta.ols.data <- OLS_results$coefficients
v.data <- OLS_results$df.residual
XTX.data <- t(X) %*% X
s2.data <- sum((OLS_results$residuals) ^ 2) / v.data



#2. Przygotowanie danych do parametrow a priori
dane_apriori_slownik <- read.csv('Weather Station Locations.csv')
dane_apriori_slownik <- dane_apriori_slownik %>% filter(Latitude >50)




dane_apriori <- read.csv('Summary of Weather.csv')
dane_apriori <- dane_apriori %>% right_join(dane_apriori_slownik, by = c('STA' = 'WBAN'))

excluded_vars <-  c('DR', 'SPD', 'SND', 'FT', 'FB', 'FTI', 'ITH', 'PGT', 'TSHDSBRSGF',
                    'SD3', 'RHX', 'RHN', 'RVG', 'WTE', 'PoorWeather', 'WindGustSpd')

dane_apriori_2 <- dane_apriori %>% select(-one_of(excluded_vars))
dane_apriori_2 <- na.omit(dane_apriori_2)


dane_apriori_2 <- dane_apriori_2 %>% filter(Precip != 'T')
dane_apriori_2 <- dane_apriori_2 %>% filter(Snowfall != '#VALUE!')


dane_apriori_2$Precip <- as.numeric(dane_apriori_2$Precip)
dane_apriori_2$Snowfall <- as.numeric(dane_apriori_2$Snowfall)
dane_apriori_2$Date <- as.Date(dane_apriori_2$Date) 


dane_apriori_2 <- dane_apriori_2 %>% group_by(STA) %>% mutate(nextday_mintemp = lead(MinTemp)) %>% 
  filter(YR %in% c(44,45))

dane_apriori_2 <- na.omit(dane_apriori_2)

dane_apriori_2 <- dane_apriori_2 %>% group_by(STA) %>%  sample_n(34)  

colnames(dane_apriori_2)[3] <- 'Rainfall'


dane_apriori_2 <- dane_apriori_2 %>% select(Date, nextday_mintemp, MinTemp, Rainfall, Snowfall, STA)


model <- lm(nextday_mintemp~+MinTemp+Rainfall+Snowfall, data = dane_apriori_2)
summary(model)



#3. Parametry a priori
Beta.prior <- model$coefficients
sm2.prior <- sd(c(model$residuals))^-2 
U.prior <- matrix(rep(0, 16), nrow = 4, ncol = 4)
U.prior[1,1]  <- vcov(model)[1,1]*1.5
U.prior[2,2]  <- vcov(model)[2,2]*1.5
U.prior[3,3]  <- vcov(model)[3,3]*1.5
U.prior[4,4]  <- vcov(model)[4,4]*1.5



v.prior <- nrow(dane_apriori_2)
k <- ncol(X)
N <- length(y)

vs2.prior <- v.prior / sm2.prior

#4. Parametry a posteriori
Beta.posterior <- solve(solve(U.prior) + XTX.data) %*% (solve(U.prior) %*% Beta.prior + XTX.data %*% Beta.ols.data)
U.posterior <- solve(solve(U.prior) + XTX.data)
v.posterior <- v.prior + N.data
vs2.posterior <- v.prior / sm2.prior + v.data * s2.data + t(Beta.ols.data-Beta.prior) %*% solve(U.prior + solve(XTX.data)) %*% (Beta.ols.data - Beta.prior)
sm2.posterior <- 1 / (vs2.posterior / v.posterior)


#5. Dane do wykresów - gêstoœæ a priori i a posteriori
beta.space <- seq(from = -1.5, to = 1.5, by = 0.01)
n_eval_points <- length(beta.space)
n_parameters <- length(Beta.posterior)
prior.marg.dens.beta <- matrix(NA, nrow = n_parameters, ncol = n_eval_points)
posterior.marg.dens.beta <- matrix(NA, nrow = n_parameters, ncol = n_eval_points)
for(ii in 1:length(Beta.posterior)) {
  prior.marg.dens.beta[ii, ] <- apply(as.matrix(beta.space), 1, dmvt,
                                      delta = Beta.prior[ii], sigma = as.matrix(U.prior[ii, ii] / sm2.prior), df = v.prior, log = FALSE)
  posterior.marg.dens.beta[ii, ] <- apply(as.matrix(beta.space), 1, dmvt,
                                          delta = Beta.posterior[ii], sigma = as.matrix(U.posterior[ii, ii] / sm2.posterior), df = v.posterior, log = FALSE)
}

grey_area <- rgb(160, 160, 160, 80, names = NULL, maxColorValue = 255)
grey_line <- rgb(80, 80, 80, 160, names = NULL, maxColorValue = 255)
green_area <- rgb(24, 121, 104, 80, names = NULL, maxColorValue = 255)
green_line <- rgb(13, 85, 72, 160, names = NULL, maxColorValue = 255)


#Polecenie, które tworzy wykres sk³adaj¹cy siê z 4 paneli (2x2)
par(mfrow = c(2, 2))
for(ii in 2:length(Beta.posterior)) {
  plot(beta.space, prior.marg.dens.beta[ii, ], las = 1, lwd = 2, bty = "n", col = grey_area,
       ylim = c(0, max(c(max(prior.marg.dens.beta[ii, ]),max(posterior.marg.dens.beta[ii, ]))) + 1), type = "l", ylab = "gêstoœæ", main = colnames(dane_aktualne)[ii + 1])
  polygon(c(beta.space, rev(beta.space)), c(prior.marg.dens.beta[ii, ], 
                                            rep(0, length(beta.space))), col = grey_area, border = NA)
  abline(v = Beta.prior[ii], col = grey_line, lwd = 3)
  text(Beta.prior[ii], max(prior.marg.dens.beta[ii, ]) + 0.4, paste("E(beta) a priori = ", Beta.prior[ii]), col = grey_line)
  abline(v = Beta.ols.data[ii], col = rgb(0, 0, 0, 1), lwd = 3)
  text(Beta.ols.data[ii], max(posterior.marg.dens.beta[ii, ]) + 0.2, paste("parametr OLS = ", round(Beta.ols.data[ii], 4)), col = rgb(0, 0, 0, 1))
  lines(beta.space, posterior.marg.dens.beta[ii, ], lwd = 2, col = green_line)
  polygon(c(beta.space, rev(beta.space)), c(posterior.marg.dens.beta[ii, ], 
                                            rep(0, length(beta.space))), col = green_area, border = NA)
  abline(v = Beta.posterior[ii], col = green_line, lwd = 3)
  text(Beta.posterior[ii], max(posterior.marg.dens.beta[ii, ]) + 0.6, paste("E(beta) a posteriori = ", round(Beta.posterior[ii], digits = 4)), col = green_line)
}

# 6. Szkicujemy HPDI wokó³ wybranego wspó³czynnika kierunkowego równania regresji
###  Wybierz:
###  2 <- Minimalna temperatura w ci¹gu dnia
###  3 <- Dobowa suma opadów deszczu (w mm)
###  4 <- Dobowa suma opadów œniegu (w mm)

ii <- 4

grey_area <- rgb(160, 160, 160, 80, names = NULL, maxColorValue = 255)
grey_line <- rgb(80, 80, 80, 160, names = NULL, maxColorValue = 255)
green_area <- rgb(24, 121, 104, 80, names = NULL, maxColorValue = 255)
green_line <- rgb(13, 85, 72, 160, names = NULL, maxColorValue = 255)
red_area <- rgb(255, 100, 123, 80, names = NULL, maxColorValue = 255)
red_line <- rgb(200, 0, 30, 160, names = NULL, maxColorValue = 255)

par(mfrow = c(1, 1))

manipulate( 
  {#Tworzymy zmienn¹ binarn¹ wskazuj¹c¹, gdzie bêdzie HPDI - tzn. o najwy¿ej gêstoœci a posteriori ponad zadany poziom
    highest <- data.frame(1:n_eval_points, posterior.marg.dens.beta[ii, ], stringsAsFactors = FALSE)
    highest <- highest[order(highest[, 2], decreasing = TRUE), ]
    highest <- data.frame(highest, cumsum(highest[, 2]) < conf_level, stringsAsFactors = FALSE)
    highest <- highest[order(highest[, 1], decreasing = FALSE), ]
    credible_set_indicator <- as.vector(as.integer(highest[, 3]))
    credible_set_begin <- match(1, credible_set_indicator)
    credible_set_end <- length(credible_set_indicator) - match(1, rev(credible_set_indicator))
    #Lewy i prawy brzeg HPDI
    x1 <- beta.space[credible_set_begin]
    x2 <- beta.space[credible_set_end]
    #Na potrzeby wykresu tworzymy wektor, który przyjmuje wartoœæ gêstoœci a posteriori w HPDI i zero poza nim
    posterior.cs <- posterior.marg.dens.beta[ii, ] * credible_set_indicator
    #Poziom ufnoœci
    HPDI_probab <- conf_level * 0.01
    #Wykres gêstoœci a posteriori
    plot(beta.space, posterior.marg.dens.beta[ii, ], las = 1, lwd = 2, bty = "n", col = green_line,
         ylim = c(0, max(posterior.marg.dens.beta[ii, ] + 1)), type = "l", ylab = "gêstoœæ", main = colnames(dane_aktualne)[ii + 1])
    polygon(c(beta.space, rev(beta.space)), 
            c(posterior.marg.dens.beta[ii, ], rep(0, length(beta.space))), 
            col = green_area, border = NA)
    text(Beta.posterior[ii], max(posterior.marg.dens.beta[ii, ]) + 0.6, paste("E(beta) a posteriori = ", round(Beta.posterior[ii], digits = 4)), col = green_line)
    abline(v = Beta.posterior[ii], col = green_line, lwd = 3)
    #Pole oznaczaj¹ce gêstoœæ a posteriori w przedziale ufnoœci HPDI
    polygon(c(beta.space, rev(beta.space)), 
            c(posterior.cs, rep(0, length(beta.space))), 
            col = red_area, border = NA)
    
    #Wyœwietl poziom ufnoœci i granice przedzia³u
    text(0, max(posterior.marg.dens.beta[ii, ]) + 0.2, paste(round(HPDI_probab * 100, digits = 1), "% przedzia³ HPDI: (", round(x1, digits = 2), " , ", round(x2, digits = 2), ")"), col = red_line)
  },
  conf_level = slider(0, 100, step = 1, initial = 95)
)


## Model 1: ze wszystkimi 4 zmiennymi
#obliczenia mamy ju¿ gotowe

# Niestety R nie poradzil sobie z obliczeniami - uzyty zostal kalkulatoro wysokiej precyzji

#P_y_M1 <- ((det(U.posterior) ^ 0.5) * gamma(v.posterior / 2) * ((vs2.posterior) ^ (- v.posterior / 2)))  /
 # ((pi ^ (N.data / 2)) * (det(U.prior) ^ 0.5) * gamma(v.prior / 2) * ((vs2.prior) ^ (- v.prior / 2)))


P_y_M1 = 2.359670293049659529759E-51



#Powtarzamy obliczenia z ograniczon¹ macierz¹ X...


ii = 2
X_2 <- X[, -c(ii)]
eval(parse(text = paste("OLS_results_2 <- lm(nextday_mintemp~ ", paste(colnames(X_2), collapse = "+"), ", data = dane_aktualne)", sep = "")))
Beta.ols.data_2 <- OLS_results_2$coefficients
v.data_2 <- OLS_results_2$df.residual
XTX.data_2 <- t(X_2) %*% X_2
s2.data_2 <- sum((OLS_results_2$residuals) ^ 2) / v.data_2

#Uwaga na inny rozk³ad a priori (ni¿sza liczba wymiarów wektora beta)
Beta.prior_2 <- Beta.prior[-c(ii)]
U.prior_2 <- U.prior[- c(ii), -c(ii)]
#Zak³adamy, ¿e rozk³ad brzegowy precyzji pozostanie bez zmian
sm2.prior_2 <- sm2.prior
v.prior_2 <- v.prior
vs2.prior_2 <- vs2.prior

Beta.posterior_2 <- solve(solve(U.prior_2) + XTX.data_2) %*% (solve(U.prior_2) %*% Beta.prior_2 + XTX.data_2 %*% Beta.ols.data_2)
U.posterior_2 <- solve(solve(U.prior_2) + XTX.data_2)
v.posterior_2 <- v.prior_2 + N.data
vs2.posterior_2 <- v.prior_2 / sm2.prior_2 + v.data_2 * s2.data_2 + t(Beta.ols.data_2 - Beta.prior_2) %*% solve(U.prior_2 + solve(XTX.data_2)) %*% (Beta.ols.data_2 - Beta.prior_2)
sm2.posterior_2 <- 1 / (vs2.posterior_2 / v.posterior_2)

#Tak jak w modelu 1 (pe³nym) obliczamy gêstoœæ brzegow¹:
#P_y_M2[ii - 1] <- ((det(U.posterior_2) ^ 0.5) * gamma(v.posterior_2 / 2) * ((vs2.posterior_2) ^ (- v.posterior_2 / 2))) /
  #((pi ^ (N.data / 2)) * (det(U.prior_2) ^ 0.5) * gamma(v.prior_2 / 2) * ((vs2.prior_2) ^ (- v.prior_2 / 2)))





ii = 3
X_2 <- X[, -c(ii)]
eval(parse(text = paste("OLS_results_2 <- lm(nextday_mintemp~ ", paste(colnames(X_2), collapse = "+"), ", data = dane_aktualne)", sep = "")))
Beta.ols.data_2 <- OLS_results_2$coefficients
v.data_2 <- OLS_results_2$df.residual
XTX.data_2 <- t(X_2) %*% X_2
s2.data_2 <- sum((OLS_results_2$residuals) ^ 2) / v.data_2

#Uwaga na inny rozk³ad a priori (ni¿sza liczba wymiarów wektora beta)
Beta.prior_2 <- Beta.prior[-c(ii)]
U.prior_2 <- U.prior[- c(ii), -c(ii)]
#Zak³adamy, ¿e rozk³ad brzegowy precyzji pozostanie bez zmian
sm2.prior_2 <- sm2.prior
v.prior_2 <- v.prior
vs2.prior_2 <- vs2.prior

Beta.posterior_2 <- solve(solve(U.prior_2) + XTX.data_2) %*% (solve(U.prior_2) %*% Beta.prior_2 + XTX.data_2 %*% Beta.ols.data_2)
U.posterior_2 <- solve(solve(U.prior_2) + XTX.data_2)
v.posterior_2 <- v.prior_2 + N.data
vs2.posterior_2 <- v.prior_2 / sm2.prior_2 + v.data_2 * s2.data_2 + t(Beta.ols.data_2 - Beta.prior_2) %*% solve(U.prior_2 + solve(XTX.data_2)) %*% (Beta.ols.data_2 - Beta.prior_2)
sm2.posterior_2 <- 1 / (vs2.posterior_2 / v.posterior_2)

#Tak jak w modelu 1 (pe³nym) obliczamy gêstoœæ brzegow¹:
#P_y_M3[ii - 1] <- ((det(U.posterior_2) ^ 0.5) * gamma(v.posterior_2 / 2) * ((vs2.posterior_2) ^ (- v.posterior_2 / 2))) /
 # ((pi ^ (N.data / 2)) * (det(U.prior_2) ^ 0.5) * gamma(v.prior_2 / 2) * ((vs2.prior_2) ^ (- v.prior_2 / 2)))


  ii = 4
  X_2 <- X[, -c(ii)]
  eval(parse(text = paste("OLS_results_2 <- lm(nextday_mintemp~ ", paste(colnames(X_2), collapse = "+"), ", data = dane_aktualne)", sep = "")))
  Beta.ols.data_2 <- OLS_results_2$coefficients
  v.data_2 <- OLS_results_2$df.residual
  XTX.data_2 <- t(X_2) %*% X_2
  s2.data_2 <- sum((OLS_results_2$residuals) ^ 2) / v.data_2
  
  #Uwaga na inny rozk³ad a priori (ni¿sza liczba wymiarów wektora beta)
  Beta.prior_2 <- Beta.prior[-c(ii)]
  U.prior_2 <- U.prior[- c(ii), -c(ii)]
  #Zak³adamy, ¿e rozk³ad brzegowy precyzji pozostanie bez zmian
  sm2.prior_2 <- sm2.prior
  v.prior_2 <- v.prior
  vs2.prior_2 <- vs2.prior
  
  Beta.posterior_2 <- solve(solve(U.prior_2) + XTX.data_2) %*% (solve(U.prior_2) %*% Beta.prior_2 + XTX.data_2 %*% Beta.ols.data_2)
  U.posterior_2 <- solve(solve(U.prior_2) + XTX.data_2)
  v.posterior_2 <- v.prior_2 + N.data
  vs2.posterior_2 <- v.prior_2 / sm2.prior_2 + v.data_2 * s2.data_2 + t(Beta.ols.data_2 - Beta.prior_2) %*% solve(U.prior_2 + solve(XTX.data_2)) %*% (Beta.ols.data_2 - Beta.prior_2)
  sm2.posterior_2 <- 1 / (vs2.posterior_2 / v.posterior_2)
  
  #Tak jak w modelu 1 (pe³nym) obliczamy gêstoœæ brzegow¹:
 # P_y_M4[ii - 1] <- ((det(U.posterior_2) ^ 0.5) * gamma(v.posterior_2 / 2) * ((vs2.posterior_2) ^ (- v.posterior_2 / 2))) /
   # ((pi ^ (N.data / 2)) * (det(U.prior_2) ^ 0.5) * gamma(v.prior_2 / 2) * ((vs2.prior_2) ^ (- v.prior_2 / 2)))

  
  

P_y_M2 <- 4.507864908397243158856E-227

P_y_M3 <- 1.054550892733039204098E-202

P_y_M4 <- 7.278606760086374987076E-204

# 2. Wyznaczamy czynniki Bayesa "dla poszczególnych zmiennych"
# (czyli, precyzyjniej, dla modeli z tymi zmiennymi i bez nich)
BF_1_2 <- P_y_M1 / P_y_M2

BF_1_3 <- P_y_M1 /  P_y_M3
  

BF_1_4 <- P_y_M1 / P_y_M4

BF_table <- data.frame(names(Beta.ols.data[2:4]), BF_1_2)

BF_table[2,2] <- BF_1_3
BF_table[3,2] <- BF_1_4
colnames(BF_table) <- c("zmienna", "czynnik Bayesa (analitycznie)")
(BF_table)


