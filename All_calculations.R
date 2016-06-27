############### Расчеты для MIDAS ##################

library('dplyr')
library("ggplot2")
library("midasr")
library("zoo")
library("formatR")
library("caret")
library("forecast") #для ARMA и AR
library("hydroGOF") #тут можно считать mse



# Наши данные с датами, 81 квартал или 243 месяца
lf_data <- read.csv("lf_data.csv") %>% tail(-2) %>% head(-1)
hf_data <- (read.csv("hf_data.csv") %>% head(-1) %>% tail(-2))[ ,-c(4,5)]
lf_data_inv <- sophisthse("INVFC_Q_I") %>% tail(-10) %>% head(-2)


#Все без дат
length(unemp)
gdp <- lf_data[, 2]
inv <- lf_data_inv[, 1] %>% as.vector()
unemp <- hf_data$unemp
cpi <- hf_data$cpi
bci <- hf_data$bci
oil <- hf_data$oil

trend <- 1:length(gdp)

MSE_GDP<- data.frame(nealmon_ind = rep(NA, 3), nbeta_ind = rep(NA, 3), nealmon_weight = rep(NA, 3), nbeta_weight = rep(NA, 3),
                     step_weight = rep(NA, 3), time_av = rep(NA, 3), ARIMA = rep(NA, 3), AR = rep(NA, 3))

MSE_inv<- data.frame(nealmon_ur = rep(NA, 3), nbeta_ur = rep(NA, 3), nealmon_weight = rep(NA, 3), nbeta_weight = rep(NA, 3),
                     step_weight = rep(NA, 3), time_av = rep(NA, 3), ARIMA = rep(NA, 3), AR = rep(NA, 3))


################## ВВП #####################
########## NEALMON_UR ##############

#мы в 3 месяце квартала 
set1 <- expand_weights_lags(
  weights = c("nealmon"), 
  from = 3, to = c(8,11), m = 1, 
  start = list(nealmon = rep(0,2)))


n_ur_1 <- midas_r_ic_table(gdp ~ trend + mls(gdp, 1:4, m = 1) + mls(unemp, 1, 3) + mls(cpi, 1, 3) + mls(bci,1, 3) + mls(oil, 1, 3), 
                           table = list(unemp = set1,
                                        cpi = set1,
                                        bci = set1,
                                        oil = set1))

avf_n_ur_1 <- average_forecast(n_ur_1$candlist,
                               data = list(trend = trend, gdp = gdp, unemp = unemp, cpi = cpi, bci = bci, oil = oil),
                               insample = 1:78, outsample = 79:81,
                               type = "recursive", 
                               measures = c("MSE"),
                               fweights = c("EW", "BICW", "MSFE", "DMSFE"),
                               show_progress = FALSE)

n <- which.min(avf_n_ur_1$accuracy$individual$MSE.in.sample)


MSE_GDP$nealmon_ur[1] <- avf_n_ur_1$accuracy$individual$MSE.out.of.sample[n]
MSE_GDP$nealmon_weight[1] <- avf_n_ur_1$accuracy$average[1,2]

#мы во 2 месяце квартала 
set2 <- expand_weights_lags(
  weights = c("nealmon"), 
  from = 4, to = c(9,12), m = 1, 
  start = list(nealmon = rep(0,2)))


n_ur_2 <- midas_r_ic_table(gdp ~ trend + mls(gdp, 1:4, m = 1) + mls(unemp, 1, 3) + mls(cpi, 1, 3) + mls(bci,1, 3) + mls(oil, 1, 3), 
                           table = list(unemp = set2,
                                        cpi = set2,
                                        bci = set2,
                                        oil = set2))

avf_n_ur_2<- average_forecast(n_ur_2$candlist,
                              data = list(trend = trend, gdp = gdp, unemp = unemp, cpi = cpi, bci = bci, oil = oil),
                              insample = 1:78, outsample = 79:81,
                              type = "recursive", 
                              measures = c("MSE"),
                              fweights = c("EW", "BICW", "MSFE", "DMSFE"),
                              show_progress = FALSE)

n2 <- which.min(avf_n_ur_2$accuracy$individual$MSE.in.sample)

MSE_GDP$nealmon_ur[2] <- avf_n_ur_2$accuracy$individual$MSE.out.of.sample[n2]
MSE_GDP$nealmon_weight[2] <- avf_n_ur_2$accuracy$average[1,2]


#мы в 1 месяце квартала
set3 <- expand_weights_lags(
  weights = c("nealmon"), 
  from = 5, to = c(10,13), m = 1, 
  start = list(nealmon = rep(0,2)))


n_ur_3 <- midas_r_ic_table(gdp ~ trend + mls(gdp, 1:4, m = 1) + mls(unemp, 1, 3) + mls(cpi, 1, 3) + mls(bci,1, 3) + mls(oil, 1, 3), 
                           table = list(unemp = set3,
                                        cpi = set3,
                                        bci = set3,
                                        oil = set3))

avf_n_ur_3<- average_forecast(n_ur_3$candlist,
                              data = list(trend = trend, gdp = gdp, unemp = unemp, cpi = cpi, bci = bci, oil = oil),
                              insample = 1:78, outsample = 79:81,
                              type = "recursive", 
                              measures = c("MSE"),
                              fweights = c("EW", "BICW", "MSFE", "DMSFE"),
                              show_progress = FALSE)

n3 <- which.min(avf_n_ur_3$accuracy$individual$MSE.in.sample)

MSE_GDP$nealmon_ur[3] <- avf_n_ur_3$accuracy$individual$MSE.out.of.sample[n3]
MSE_GDP$nealmon_weight[3] <- avf_n_ur_3$accuracy$average[1,2]

############################## NBETA ##############################


#мы в 3 месяце квартала
beta_set1 <- expand_weights_lags(
  weights = c("nbeta"), 
  from = 3, to = c(8,11), m = 1, 
  start = list(nbeta = c(0,1,0)))


b_ur_1 <- midas_r_ic_table(gdp ~ trend + mls(gdp, 1:4, m = 1) + mls(unemp, 1, 3) + mls(cpi, 1, 3) + mls(bci,1, 3) + mls(oil, 1, 3), 
                           table = list(unemp = beta_set1,
                                        cpi = beta_set1,
                                        bci = beta_set1,
                                        oil = beta_set1))

avf_b_ur_1 <- average_forecast(b_ur_1$candlist,
                               data = list(trend = trend, gdp = gdp, unemp = unemp, cpi = cpi, bci = bci, oil = oil),
                               insample = 1:78, outsample = 79:81,
                               type = "recursive", 
                               measures = c("MSE"),
                               fweights = c("EW", "BICW", "MSFE", "DMSFE"),
                               show_progress = FALSE)

beta_n <- which.min(avf_b_ur_1$accuracy$individual$MSE.in.sample)


MSE_GDP$nbeta_ur[1] <- avf_b_ur_1$accuracy$individual$MSE.out.of.sample[beta_n]
MSE_GDP$nbeta_weight[1] <- avf_b_ur_1$accuracy$average[1,2]

#мы во 2 месяце квартала
beta_set2 <- expand_weights_lags(
  weights = c("nbeta"), 
  from = 4, to = c(9,12), m = 1, 
  start = list(nbeta = c(0, 1, 0)))


b_ur_2 <- midas_r_ic_table(gdp ~ trend + mls(gdp, 1:4, m = 1) + mls(unemp, 1, 3) + mls(cpi, 1, 3) + mls(bci,1, 3) + mls(oil, 1, 3), 
                           table = list(unemp = beta_set2,
                                        cpi = beta_set2,
                                        bci = beta_set2,
                                        oil = beta_set2))

avf_b_ur_2<- average_forecast(b_ur_2$candlist,
                              data = list(trend = trend, gdp = gdp, unemp = unemp, cpi = cpi, bci = bci, oil = oil),
                              insample = 1:78, outsample = 79:81,
                              type = "recursive", 
                              measures = c("MSE"),
                              fweights = c("EW", "BICW", "MSFE", "DMSFE"),
                              show_progress = FALSE)

beta_n2 <- which.min(avf_b_ur_2$accuracy$individual$MSE.in.sample)

MSE_GDP$nbeta_ur[2] <- avf_b_ur_2$accuracy$individual$MSE.out.of.sample[beta_n2]
MSE_GDP$nbeta_weight[2] <- avf_b_ur_2$accuracy$average[1,2]


#мы в 1 месяце квартала 
beta_set3 <- expand_weights_lags(
  weights = c("nbeta"), 
  from = 5, to = c(10,13), m = 1, 
  start = list(nbeta = c(0,1,0)))


b_ur_3 <- midas_r_ic_table(gdp ~ trend + mls(gdp, 1:4, m = 1) + mls(unemp, 1, 3) + mls(cpi, 1, 3) + mls(bci,1, 3) + mls(oil, 1, 3), 
                           table = list(unemp = beta_set3,
                                        cpi = beta_set3,
                                        bci = beta_set3,
                                        oil = beta_set3))

avf_b_ur_3<- average_forecast(b_ur_3$candlist,
                              data = list(trend = trend, gdp = gdp, unemp = unemp, cpi = cpi, bci = bci, oil = oil),
                              insample = 1:78, outsample = 79:81,
                              type = "recursive", 
                              measures = c("MSE"),
                              fweights = c("EW", "BICW", "MSFE", "DMSFE"),
                              show_progress = FALSE)

beta_n3 <- which.min(avf_b_ur_3$accuracy$individual$MSE.in.sample)

MSE_GDP$nbeta_ur[3] <- avf_b_ur_3$accuracy$individual$MSE.out.of.sample[beta_n3]
MSE_GDP$nbeta_weight[3] <- avf_b_ur_3$accuracy$average[1,2]




############################### Unrestricted (step weightning) ###########################
gdp_tr = gdp[1:78]
unemp_tr = unemp[1:234]
bci_tr = bci[1:234]
cpi_tr = cpi[1:234]
oil_tr = oil[1:234]
trend_tr = 1:78

#мы в 3 месяце квартала 
ur_models1 <- list()
future1 <- list()
mse_in1 <- c()
mse_out1 <- c()

k=1
for (i1 in 1:4){
  for (i2 in 1:4){
    for (i3 in 1:4){
      for (i4 in 1:4){
        ur_models1[[k]] <- midas_r(gdp_tr ~ trend_tr + mls(gdp_tr, 3:6, 1) + mls(unemp_tr, 9:(i1+13), 3) + 
                                     mls(cpi_tr, 9:(i2+13), 3) + mls(bci_tr,9:(i3+13), 3) + mls(oil_tr, 9:(i4+13), 3), start = NULL)
        future1[[k]] <- forecast(ur_models1[[k]], newdata = list(trend_tr = 79:81, unemp_tr = rep(NA, 9), 
                                                                 cpi_tr = rep(NA, 9), bci_tr = rep(NA, 9), oil_tr = rep(NA, 9)))
        mse_in1[k] <- mse(future1[[k]]$fitted, gdp[7:78])
        mse_out1[k] <- mse(future1[[k]]$mean, gdp[79:81])
        k = k + 1
      }
    }
  }
}

MSE_GDP$step_weight[1] <- mse_out1[which.min(mse_in1)]

#мы вo 2 месяце квартала 
ur_models2 <- list()
future2 <- list()
mse_in2 <- c()
mse_out2 <- c()

k=1
for (i1 in 1:4){
  for (i2 in 1:4){
    for (i3 in 1:4){
      for (i4 in 1:4){
        ur_models2[[k]] <- midas_r(gdp_tr ~ trend_tr + mls(gdp_tr, 3:6, 1) + mls(unemp_tr, 10:(i1+14), 3) + 
                                     mls(cpi_tr, 10:(i2+14), 3) + mls(bci_tr,10:(i3+14), 3) + mls(oil_tr, 10:(i4+14), 3), start = NULL)
        future2[[k]] <- forecast(ur_models2[[k]], newdata = list(trend_tr = 79:81, unemp_tr = rep(NA, 9), 
                                                                 cpi_tr = rep(NA, 9), bci_tr = rep(NA, 9), oil_tr = rep(NA, 9)))
        mse_in2[k] <- mse(future2[[k]]$fitted, gdp[7:78])
        mse_out2[k] <- mse(future2[[k]]$mean, gdp[79:81])
        k = k + 1
      }
    }
  }
}
ur_models1[[34]]
ur_models2[[34]]
MSE_GDP$step_weight[2] <- mse_out2[which.min(mse_in2)]


#мы в 1 месяце квартала 

ur_models3 <- list()
future3 <- list()
mse_in3 <- c()
mse_out3 <- c()

k=1
for (i1 in 1:4){
  for (i2 in 1:4){
    for (i3 in 1:4){
      for (i4 in 1:4){
        ur_models3[[k]] <- midas_r(gdp_tr ~ trend_tr + mls(gdp_tr, 3:6, 1) + mls(unemp_tr, 11:(i1+15), 3) + 
                                     mls(cpi_tr, 11:(i2+15), 3) + mls(bci_tr,11:(i3+15), 3) + mls(oil_tr, 11:(i4+15), 3), start = NULL)
        future3[[k]] <- forecast(ur_models3[[k]], newdata = list(trend_tr = 79:81, unemp_tr = rep(NA, 9), 
                                                                 cpi_tr = rep(NA, 9), bci_tr = rep(NA, 9), oil_tr = rep(NA, 9)))
        mse_in3[k] <- mse(future3[[k]]$fitted, gdp[7:78])
        mse_out3[k] <- mse(future3[[k]]$mean, gdp[79:81])
        k = k + 1
      }
    }
  }
}

MSE_GDP$step_weight[3] <- mse_out3[which.min(mse_in3)]
################################## Time Averaging #################################



#Агрегируем данные
data_1 <- data.frame(gdp = gdp[1:78], unemp = rep(NA, 78), cpi = rep(NA, 78), bci = rep(NA, 78), oil = rep(NA, 78))
data_2 <- data.frame(gdp = gdp[1:78], unemp = rep(NA, 78), cpi = rep(NA, 78), bci = rep(NA, 78), oil = rep(NA, 78))
data_3 <- data.frame(gdp = gdp[1:78], unemp = rep(NA, 78), cpi = rep(NA, 78), bci = rep(NA, 78), oil = rep(NA, 78))

for (i in 2:5)
{
  k = 1
  for (j in seq(1, 231, 3)){
    data_3[k,i] <- mean(hf_data[ j:j+3,i])
    data_2[k,i] <- mean(hf_data[ j:j+3,i])
    data_1[k,i] <- mean(hf_data[ j:j+3,i])
    k = k+1
  }
}

for (i in 2:5){
  data_3[78,i] <- mean(hf_data[ 232:234,i])
  data_2[78,i] <- mean(hf_data[ 232:233,i])
  data_1[78,i] <- hf_data[ 232,i]
}

#мы в 3 месяце квартала
ta_models1 <- list()
ta_future1 <- list()
ta_mse_in1 <- c()
ta_mse_out1 <- c()

k=1
for (i1 in 1:2){
  for (i2 in 1:2){
    for (i3 in 1:2){
      for (i4 in 1:2){
        ta_models1[[k]] <- midas_r(gdp ~ trend_tr + mls(gdp, 3:6, 1) + mls(unemp, 3:(i1+3), 1) + 
                                     mls(cpi, 3:(i2+3), 1) + mls(bci,3:(i3+3), 1) + 
                                     mls(oil, 3:(i4+3), 1), data = list(data_3),start = NULL)
        ta_future1[[k]] <- forecast(ta_models1[[k]], newdata = list(trend_tr = 79:81, unemp = rep(NA, 3), 
                                                                    cpi = rep(NA, 3), bci = rep(NA, 3), oil = rep(NA, 3)))
        ta_mse_in1[k] <- mse(ta_future1[[k]]$fitted, gdp[7:78])
        ta_mse_out1[k] <- mse(ta_future1[[k]]$mean, gdp[79:81])
        k = k + 1
      }
    }
  }
}
MSE_GDP$time_av[1] <- ta_mse_out1[which.min(ta_mse_in1)]

#мы во 2 месяце квартала 

ta_models2 <- list()
ta_future2 <- list()
ta_mse_in2 <- c()
ta_mse_out2 <- c()

k=1
for (i1 in 1:2){
  for (i2 in 1:2){
    for (i3 in 1:2){
      for (i4 in 1:2){
        ta_models2[[k]] <- midas_r(gdp ~ trend_tr + mls(gdp, 3:6, 1) + mls(unemp, 3:(i1+3), 1) + 
                                     mls(cpi, 3:(i2+3), 1) + mls(bci,3:(i3+3), 1) + 
                                     mls(oil, 3:(i4+3), 1), data = list(data_2),start = NULL)
        ta_future2[[k]] <- forecast(ta_models2[[k]], newdata = list(trend_tr = 79:81, unemp = rep(NA, 3), 
                                                                    cpi = rep(NA, 3), bci = rep(NA, 3), oil = rep(NA, 3)))
        ta_mse_in2[k] <- mse(ta_future2[[k]]$fitted, gdp[7:78])
        ta_mse_out2[k] <- mse(ta_future2[[k]]$mean, gdp[79:81])
        k = k + 1
      }
    }
  }
}
MSE_GDP$time_av[2] <- ta_mse_out2[which.min(ta_mse_in2)]

#мы в 1 месяце квартала 


ta_models3 <- list()
ta_future3 <- list()
ta_mse_in3 <- c()
ta_mse_out3 <- c()

k=1
for (i1 in 1:2){
  for (i2 in 1:2){
    for (i3 in 1:2){
      for (i4 in 1:2){
        ta_models3[[k]] <- midas_r(gdp ~ trend_tr + mls(gdp, 3:6, 1) + mls(unemp, 3:(i1+3), 1) + 
                                     mls(cpi, 3:(i2+3), 1) + mls(bci,3:(i3+3), 1) + 
                                     mls(oil, 3:(i4+3), 1), data = list(data_1),start = NULL)
        ta_future3[[k]] <- forecast(ta_models3[[k]], newdata = list(trend_tr = 79:81, unemp = rep(NA, 3), 
                                                                    cpi = rep(NA, 3), bci = rep(NA, 3), oil = rep(NA, 3)))
        ta_mse_in3[k] <- mse(ta_future3[[k]]$fitted, gdp[7:78])
        ta_mse_out3[k] <- mse(ta_future3[[k]]$mean, gdp[79:81])
        k = k + 1
      }
    }
  }
}
MSE_GDP$time_av[3] <- ta_mse_out3[which.min(ta_mse_in3)]

####################### ARIMA ##############################
arima <- forecast(auto.arima(gdp[1:78]))$mean %>% as.vector()

MSE_GDP$ARIMA[1] <- mse(arima[1:3], gdp[79:81])
MSE_GDP$ARIMA[2] <- mse(arima[1:3], gdp[79:81])
MSE_GDP$ARIMA[3] <- mse(arima[1:3], gdp[79:81])


######################## AR ###########################
ar <- forecast(ar(gdp[1:78]))$mean %>% as.vector()

MSE_GDP$AR[1] <- mse(ar[1:3], gdp[79:81])
MSE_GDP$AR[2] <- mse(ar[1:3], gdp[79:81])
MSE_GDP$AR[3] <- mse(ar[1:3], gdp[79:81])


write.csv(MSE_GDP, "MSE_GDP.csv")


#################### Инвестиции ##################################
########## NEALMON_UR ##############

#мы в 3 месяце квартала 
inv_set1 <- expand_weights_lags(
  weights = c("nealmon"), 
  from = 3, to = c(8,11), m = 1, 
  start = list(nealmon = rep(0, 2)))


inv_n_ur_1 <- midas_r_ic_table(inv ~ trend + mls(inv, 1:4, m = 1) + mls(unemp, 1, 3) + mls(cpi, 1, 3) + mls(bci,1, 3) + mls(oil, 1, 3), 
                               table = list(unemp = inv_set1,
                                            cpi = inv_set1,
                                            bci = inv_set1,
                                            oil = inv_set1))

inv_avf_n_ur_1 <- average_forecast(inv_n_ur_1$candlist,
                                   data = list(trend = trend, inv = inv, unemp = unemp, cpi = cpi, bci = bci, oil = oil),
                                   insample = 1:78, outsample = 79:81,
                                   type = "recursive", 
                                   measures = c("MSE"),
                                   fweights = c("EW", "BICW", "MSFE", "DMSFE"),
                                   show_progress = FALSE)

inv_n <- which.min(inv_avf_n_ur_1$accuracy$individual$MSE.in.sample)


MSE_inv$nealmon_ur[1] <- inv_avf_n_ur_1$accuracy$individual$MSE.out.of.sample[inv_n]
MSE_inv$nealmon_weight[1] <- inv_avf_n_ur_1$accuracy$average[1,2]

#мы во 2 месяце квартала
set2 <- expand_weights_lags(
  weights = c("nealmon"), 
  from = 4, to = c(9,12), m = 1, 
  start = list(nealmon = rep(0,2)))


inv_n_ur_2 <- midas_r_ic_table(inv ~ trend + mls(inv, 1:4, m = 1) + mls(unemp, 1, 3) + mls(cpi, 1, 3) + mls(bci,1, 3) + mls(oil, 1, 3), 
                               table = list(unemp = set2,
                                            cpi = set2,
                                            bci = set2,
                                            oil = set2))

inv_avf_n_ur_2<- average_forecast(inv_n_ur_2$candlist,
                                  data = list(trend = trend, inv = inv, unemp = unemp, cpi = cpi, bci = bci, oil = oil),
                                  insample = 1:78, outsample = 79:81,
                                  type = "recursive", 
                                  measures = c("MSE"),
                                  fweights = c("EW", "BICW", "MSFE", "DMSFE"),
                                  show_progress = FALSE)

inv_n2 <- which.min(inv_avf_n_ur_2$accuracy$individual$MSE.in.sample)

MSE_inv$nealmon_ur[2] <- inv_avf_n_ur_2$accuracy$individual$MSE.out.of.sample[inv_n2]
MSE_inv$nealmon_weight[2] <- inv_avf_n_ur_2$accuracy$average[1,2]


#мы в 1 месяце квартала  
set3 <- expand_weights_lags(
  weights = c("nealmon"), 
  from = 5, to = c(10,13), m = 1, 
  start = list(nealmon = rep(0,2)))


inv_n_ur_3 <- midas_r_ic_table(inv ~ trend + mls(inv, 1:4, m = 1) + mls(unemp, 1, 3) + mls(cpi, 1, 3) + mls(bci,1, 3) + mls(oil, 1, 3), 
                               table = list(unemp = set3,
                                            cpi = set3,
                                            bci = set3,
                                            oil = set3))

inv_avf_n_ur_3<- average_forecast(inv_n_ur_3$candlist,
                                  data = list(trend = trend, inv = inv, unemp = unemp, cpi = cpi, bci = bci, oil = oil),
                                  insample = 1:78, outsample = 79:81,
                                  type = "recursive", 
                                  measures = c("MSE"),
                                  fweights = c("EW", "BICW", "MSFE", "DMSFE"),
                                  show_progress = FALSE)

inv_n3 <- which.min(inv_avf_n_ur_3$accuracy$individual$MSE.in.sample)

MSE_inv$nealmon_ur[3] <- inv_avf_n_ur_3$accuracy$individual$MSE.out.of.sample[inv_n3]
MSE_inv$nealmon_weight[3] <- inv_avf_n_ur_3$accuracy$average[1,2]

############################## NBETA ##############################


#мы в 3 месяце квартала 
beta_set1 <- expand_weights_lags(
  weights = c("nbeta"), 
  from = 3, to = c(8,11), m = 1, 
  start = list(nbeta = c(0,1,0)))


inv_b_ur_1 <- midas_r_ic_table(inv ~ trend + mls(inv, 1:4, m = 1) + mls(unemp, 1, 3) + mls(cpi, 1, 3) + mls(bci,1, 3) + mls(oil, 1, 3), 
                               table = list(unemp = beta_set1,
                                            cpi = beta_set1,
                                            bci = beta_set1,
                                            oil = beta_set1))

inv_avf_b_ur_1 <- average_forecast(inv_b_ur_1$candlist,
                                   data = list(trend = trend, inv = inv, unemp = unemp, cpi = cpi, bci = bci, oil = oil),
                                   insample = 1:78, outsample = 79:81,
                                   type = "recursive", 
                                   measures = c("MSE"),
                                   fweights = c("EW", "BICW", "MSFE", "DMSFE"),
                                   show_progress = FALSE)

inv_beta_n <- which.min(inv_avf_b_ur_1$accuracy$individual$MSE.in.sample)


MSE_inv$nbeta_ur[1] <- inv_avf_b_ur_1$accuracy$individual$MSE.out.of.sample[inv_beta_n]
MSE_inv$nbeta_weight[1] <- inv_avf_b_ur_1$accuracy$average[1,2]

#мы во 2 месяце квартала  
beta_set2 <- expand_weights_lags(
  weights = c("nbeta"), 
  from = 4, to = c(9,12), m = 1, 
  start = list(nbeta = c(0, 1, 0)))


inv_b_ur_2 <- midas_r_ic_table(inv ~ trend + mls(inv, 1:4, m = 1) + mls(unemp, 1, 3) + mls(cpi, 1, 3) + mls(bci,1, 3) + mls(oil, 1, 3), 
                               table = list(unemp = beta_set2,
                                            cpi = beta_set2,
                                            bci = beta_set2,
                                            oil = beta_set2))

inv_avf_b_ur_2<- average_forecast(inv_b_ur_2$candlist,
                                  data = list(trend = trend, inv = inv, unemp = unemp, cpi = cpi, bci = bci, oil = oil),
                                  insample = 1:78, outsample = 79:81,
                                  type = "recursive", 
                                  measures = c("MSE"),
                                  fweights = c("EW", "BICW", "MSFE", "DMSFE"),
                                  show_progress = FALSE)

inv_beta_n2 <- which.min(avf_b_ur_2$accuracy$individual$MSE.in.sample)

MSE_inv$nbeta_ur[2] <- avf_b_ur_2$accuracy$individual$MSE.out.of.sample[beta_n2]
MSE_inv$nbeta_weight[2] <- avf_b_ur_2$accuracy$average[1,2]


#мы в 1 месяце квартала 
beta_set3 <- expand_weights_lags(
  weights = c("nbeta"), 
  from = 5, to = c(10,13), m = 1, 
  start = list(nbeta = c(0,1,0)))


inv_b_ur_3 <- midas_r_ic_table(inv ~ trend + mls(inv, 1:4, m = 1) + mls(unemp, 1, 3) + mls(cpi, 1, 3) + mls(bci,1, 3) + mls(oil, 1, 3), 
                               table = list(unemp = beta_set3,
                                            cpi = beta_set3,
                                            bci = beta_set3,
                                            oil = beta_set3))

inv_avf_b_ur_3<- average_forecast(inv_b_ur_3$candlist,
                                  data = list(trend = trend, inv = inv, unemp = unemp, cpi = cpi, bci = bci, oil = oil),
                                  insample = 1:78, outsample = 79:81,
                                  type = "recursive", 
                                  measures = c("MSE"),
                                  fweights = c("EW", "BICW", "MSFE", "DMSFE"),
                                  show_progress = FALSE)

inv_beta_n3 <- which.min(inv_avf_b_ur_3$accuracy$individual$MSE.in.sample)

MSE_inv$nbeta_ur[3] <- inv_avf_b_ur_3$accuracy$individual$MSE.out.of.sample[inv_beta_n3]
MSE_inv$nbeta_weight[3] <- inv_avf_b_ur_3$accuracy$average[1,2]




############################### Unrestricted (step weightning) ###########################
inv_tr = inv[1:78]
unemp_tr = unemp[1:234]
bci_tr = bci[1:234]
cpi_tr = cpi[1:234]
oil_tr = oil[1:234]
trend_tr = 1:78

#мы в 3 месяце квартала 
inv_ur_models1 <- list()
inv_future1 <- list()
inv_mse_in1 <- c()
inv_mse_out1 <- c()

k=1
for (i1 in 1:4){
  for (i2 in 1:4){
    for (i3 in 1:4){
      for (i4 in 1:4){
        inv_ur_models1[[k]] <- midas_r(inv_tr ~ trend_tr + mls(inv_tr, 3:6, 1) + mls(unemp_tr, 9:(i1+13), 3) + 
                                         mls(cpi_tr, 9:(i2+13), 3) + mls(bci_tr,9:(i3+13), 3) + mls(oil_tr, 9:(i4+13), 3), start = NULL)
        inv_future1[[k]] <- forecast(inv_ur_models1[[k]], newdata = list(trend_tr = 79:81, unemp_tr = rep(NA, 9), 
                                                                         cpi_tr = rep(NA, 9), bci_tr = rep(NA, 9), oil_tr = rep(NA, 9)))
        inv_mse_in1[k] <- mse(inv_future1[[k]]$fitted, inv[7:78])
        inv_mse_out1[k] <- mse(inv_future1[[k]]$mean, inv[79:81])
        k = k + 1
      }
    }
  }
}

MSE_inv$step_weight[1] <- inv_mse_out1[which.min(inv_mse_in1)]

#мы вo 2 месяце квартала 
inv_ur_models2 <- list()
inv_future2 <- list()
inv_mse_in2 <- c()
inv_mse_out2 <- c()

k=1
for (i1 in 1:4){
  for (i2 in 1:4){
    for (i3 in 1:4){
      for (i4 in 1:4){
        inv_ur_models2[[k]] <- midas_r(inv_tr ~ trend_tr + mls(inv_tr, 3:6, 1) + mls(unemp_tr, 10:(i1+14), 3) + 
                                         mls(cpi_tr, 10:(i2+14), 3) + mls(bci_tr,10:(i3+14), 3) + mls(oil_tr, 10:(i4+14), 3), start = NULL)
        inv_future2[[k]] <- forecast(inv_ur_models2[[k]], newdata = list(trend_tr = 79:81, unemp_tr = rep(NA, 9), 
                                                                         cpi_tr = rep(NA, 9), bci_tr = rep(NA, 9), oil_tr = rep(NA, 9)))
        inv_mse_in2[k] <- mse(inv_future2[[k]]$fitted, inv[7:78])
        inv_mse_out2[k] <- mse(inv_future2[[k]]$mean, inv[79:81])
        k = k + 1
      }
    }
  }
}

MSE_inv$step_weight[2] <- inv_mse_out2[which.min(inv_mse_in2)]


#мы в 1 месяце квартала 

inv_ur_models3 <- list()
inv_future3 <- list()
inv_mse_in3 <- c()
inv_mse_out3 <- c()

k=1
for (i1 in 1:4){
  for (i2 in 1:4){
    for (i3 in 1:4){
      for (i4 in 1:4){
        inv_ur_models3[[k]] <- midas_r(inv_tr ~ trend_tr + mls(inv_tr, 3:6, 1) + mls(unemp_tr, 11:(i1+15), 3) + 
                                         mls(cpi_tr, 11:(i2+15), 3) + mls(bci_tr,11:(i3+15), 3) + mls(oil_tr, 11:(i4+15), 3), start = NULL)
        inv_future3[[k]] <- forecast(inv_ur_models3[[k]], newdata = list(trend_tr = 79:81, unemp_tr = rep(NA, 9), 
                                                                         cpi_tr = rep(NA, 9), bci_tr = rep(NA, 9), oil_tr = rep(NA, 9)))
        inv_mse_in3[k] <- mse(inv_future3[[k]]$fitted, inv[7:78])
        inv_mse_out3[k] <- mse(inv_future3[[k]]$mean, inv[79:81])
        k = k + 1
      }
    }
  }
}

MSE_inv$step_weight[3] <- inv_mse_out3[which.min(inv_mse_in3)]


################################## Time Averaging #################################



#Агрегируем данные
inv_data_1 <- data.frame(inv = inv[1:78], unemp = rep(NA, 78), cpi = rep(NA, 78), bci = rep(NA, 78), oil = rep(NA, 78))
inv_data_2 <- data.frame(inv = inv[1:78], unemp = rep(NA, 78), cpi = rep(NA, 78), bci = rep(NA, 78), oil = rep(NA, 78))
inv_data_3 <- data.frame(inv = inv[1:78], unemp = rep(NA, 78), cpi = rep(NA, 78), bci = rep(NA, 78), oil = rep(NA, 78))

for (i in 2:5)
{
  k = 1
  for (j in seq(1, 231, 3)){
    inv_data_3[k,i] <- mean(hf_data[ j:j+3,i])
    inv_data_2[k,i] <- mean(hf_data[ j:j+3,i])
    inv_data_1[k,i] <- mean(hf_data[ j:j+3,i])
    k = k+1
  }
}

for (i in 2:5){
  inv_data_3[78,i] <- mean(hf_data[ 232:234,i])
  inv_data_2[78,i] <- mean(hf_data[ 232:233,i])
  inv_data_1[78,i] <- hf_data[ 232,i]
}

#мы в 3 месяце квартала 
inv_ta_models1 <- list()
inv_ta_future1 <- list()
inv_ta_mse_in1 <- c()
inv_ta_mse_out1 <- c()

k=1
for (i1 in 1:2){
  for (i2 in 1:2){
    for (i3 in 1:2){
      for (i4 in 1:2){
        inv_ta_models1[[k]] <- midas_r(inv ~ trend_tr + mls(inv, 3:6, 1) + mls(unemp, 3:(i1+3), 1) + 
                                         mls(cpi, 3:(i2+3), 1) + mls(bci,3:(i3+3), 1) + 
                                         mls(oil, 3:(i4+3), 1), data = list(inv_data_3),start = NULL)
        inv_ta_future1[[k]] <- forecast(inv_ta_models1[[k]], newdata = list(trend_tr = 79:81, unemp = rep(NA, 3), 
                                                                            cpi = rep(NA, 3), bci = rep(NA, 3), oil = rep(NA, 3)))
        inv_ta_mse_in1[k] <- mse(inv_ta_future1[[k]]$fitted, inv[7:78])
        inv_ta_mse_out1[k] <- mse(inv_ta_future1[[k]]$mean, inv[79:81])
        k = k + 1
      }
    }
  }
}
MSE_inv$time_av[1] <- inv_ta_mse_out1[which.min(inv_ta_mse_in1)]

#мы во 2 месяце квартала

inv_ta_models2 <- list()
inv_ta_future2 <- list()
inv_ta_mse_in2 <- c()
inv_ta_mse_out2 <- c()

k=1
for (i1 in 1:2){
  for (i2 in 1:2){
    for (i3 in 1:2){
      for (i4 in 1:2){
        inv_ta_models2[[k]] <- midas_r(inv ~ trend_tr + mls(inv, 3:6, 1) + mls(unemp, 3:(i1+3), 1) + 
                                         mls(cpi, 3:(i2+3), 1) + mls(bci,3:(i3+3), 1) + 
                                         mls(oil, 3:(i4+3), 1), data = list(inv_data_2),start = NULL)
        inv_ta_future2[[k]] <- forecast(inv_ta_models2[[k]], newdata = list(trend_tr = 79:81, unemp = rep(NA, 3), 
                                                                            cpi = rep(NA, 3), bci = rep(NA, 3), oil = rep(NA, 3)))
        inv_ta_mse_in2[k] <- mse(inv_ta_future2[[k]]$fitted, inv[7:78])
        inv_ta_mse_out2[k] <- mse(inv_ta_future2[[k]]$mean, inv[79:81])
        k = k + 1
      }
    }
  }
}
MSE_inv$time_av[2] <- inv_ta_mse_out2[which.min(inv_ta_mse_in2)]

#мы в 1 месяце квартала 


inv_ta_models3 <- list()
inv_ta_future3 <- list()
inv_ta_mse_in3 <- c()
inv_ta_mse_out3 <- c()

k=1
for (i1 in 1:2){
  for (i2 in 1:2){
    for (i3 in 1:2){
      for (i4 in 1:2){
        inv_ta_models3[[k]] <- midas_r(inv ~ trend_tr + mls(inv, 3:6, 1) + mls(unemp, 3:(i1+3), 1) + 
                                         mls(cpi, 3:(i2+3), 1) + mls(bci,3:(i3+3), 1) + 
                                         mls(oil, 3:(i4+3), 1), data = list(inv_data_1),start = NULL)
        inv_ta_future3[[k]] <- forecast(inv_ta_models3[[k]], newdata = list(trend_tr = 79:81, unemp = rep(NA, 3), 
                                                                            cpi = rep(NA, 3), bci = rep(NA, 3), oil = rep(NA, 3)))
        inv_ta_mse_in3[k] <- mse(inv_ta_future3[[k]]$fitted, inv[7:78])
        inv_ta_mse_out3[k] <- mse(inv_ta_future3[[k]]$mean, inv[79:81])
        k = k + 1
      }
    }
  }
}
MSE_inv$time_av[3] <- inv_ta_mse_out3[which.min(inv_ta_mse_in3)]

####################### ARIMA ##############################
inv_arima <- forecast(auto.arima(inv[1:78]))$mean %>% as.vector()

MSE_inv$ARIMA[1] <- mse(inv_arima[1:3], inv[79:81])
MSE_inv$ARIMA[2] <- mse(inv_arima[1:3], inv[79:81])
MSE_inv$ARIMA[3] <- mse(inv_arima[1:3], inv[79:81])


######################## AR ###########################
inv_ar <- forecast(ar(inv[1:78]))$mean %>% as.vector()

MSE_inv$AR[1] <- mse(inv_ar[1:3], inv[79:81])
MSE_inv$AR[2] <- mse(inv_ar[1:3], inv[79:81])
MSE_inv$AR[3] <- mse(inv_ar[1:3], inv[79:81])

write.csv(MSE_inv, "MSE_inv.csv")

