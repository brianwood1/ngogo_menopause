library("car")
library("lmer")
library("lmerTest")
library("mgcv")
library("gamm4")

# Model fitting code for Demographic and hormonal evidence for menopause in wild chimpanzees

# The datasets that are available from the corresponding authors with reasonable request are
# LH_data 
# FSH_data_2016_2018
# FSH_data_2006_2007 
# estradiol_data 
# estrone_data 
# pregnanediol_data 

# Functions to fit the models
fit_categorical_LH_model <- function(d)
{
  x <- powerTransform(d$LH_SG,family="bcPower")
  p <- x$lambda[[1]]
  m <-lmer((LH_SG^p)/p~Reproductive_Status_LH_FSH + (1|Chimp.name),data=d)
  return(m)
}

fit_categorical_FSH_model_2016_2018 <- function(d)
{
  #d <- get_FSH_data_2016_2018()
  x <- powerTransform(d$FSH_SG,family="bcPower")
  p <- x$lambda[[1]]
  m <-lmer((FSH_SG^p)/p~Reproductive_Status_LH_FSH + (1|Chimp.name),data=d)
  return(m)
}

fit_categorical_FSH_model_2006_2007 <- function(d)
{
  #d <- get_FSH_data_2006_2007()
  x <- powerTransform(d$FSH_SG,family="bcPower")
  p <- x$lambda[[1]]
  m <-lmer((FSH_SG^p)/p~Reproductive_Status_LH_FSH + (1|Chimp.name),data=d)
  return(m)
}

fit_categorical_estradiol_model <- function(d)
{
  fit_model <-lmer((Estradiol_SG^0.0279222)/0.0279222~Reproductive_Status_Estradiol_Estrone_Pregnanediol_Cat + (1|Chimp.name),data=d)
  return(fit_model)
}

fit_categorical_estrone_model <- function(d)
{
  fit_model<-lmer((Estrone_SG^-0.3300254)/-0.3300254~Reproductive_Status_Estradiol_Estrone_Pregnanediol_Cat+(1|Chimp.name),data=d)
  return(fit_model)
}

fit_categorical_pregnanediol_model <- function(d)
{
  fit_model <-lmer((Pregnanediol_SG^-0.1911354)/-0.1911354~Reproductive_Status_Estradiol_Estrone_Pregnanediol_Cat+(1|Chimp.name),data=d)
  return(fit_model)
}

fit_age_linear_null_log_LH <- function(d)
{
  g <- gam(log_LH_SG ~ Age + s(Chimp.name,k=60, bs="re"), method="REML", family="gaussian", data=d)
  return(g)
}

fit_age_linear_null_log_FSH_2016_2018 <- function(d)
{
  #d <- get_FSH_data_2016_2018()
  m <-gam(Log_FSH_SG ~ Age + s(Chimp.name, k=60, bs="re"), method="REML", data=d)
  return(m)
}

fit_age_linear_null_log_FSH_2006_2007 <- function(d)
{
  #d <- get_FSH_data_2006_2007()
  m <-gam(Log_FSH_SG ~ Age + s(Chimp.name, k=20, bs="re"), method="REML", data=d)
  return(m)
}

fit_age_gam_log_LH <- function(d)
{
  #d <- get_LH_data()
  g <- gam(log_LH_SG ~ s(Age,k=10, bs="cr") + s(Chimp.name,k=60, bs="re"), method="REML", family="gaussian", data=d)
  return(g)  
}

fit_age_gam_log_FSH_2006_2007 <- function(d)
{
  #d <- get_FSH_data_2006_2007()
  m <-gam(Log_FSH_SG ~ s(Age,k=10,bs="cr") + s(Chimp.name, k=20, bs="re"), method="REML", data=d)
  return(m)
}

fit_age_gam_log_estradiol <- function(d)
{
  #d <- get_estradiol_data()
  g2 <- gam(log_Estradiol_SG ~ s(Age) + s(Chimp.name,k=60, bs="re"), method="REML", family="gaussian", data=d)
  return(g2)
}

fit_age_gam_log_estrone <- function(d)
{
  #d <- get_estrone_data()
  g2 <- gam(log_Estrone ~ s(Age) + s(Chimp.name,k=60, bs="re"), method="REML", family="gaussian", data=d)
  return(g2)
}

fit_age_gam_log_pregnanediol <- function(d)
{
  #d <- get_pregnanediol_data()
  g2 <- gam(log_Pregnanediol_SG ~ s(Age) + s(Chimp.name,k=60, bs="re"), method="REML", family="gaussian", data=d)
  return(g2)
}

# Model 1 - LH as a function of repro. category
model_1 <- fit_categorical_LH_model(LH_data)

# Model 2 - FSH 2016-2018 as a function of repro. category
model_2 <- fit_categorical_FSH_model_2016_2018(d=FSH_data_2016_2018)

# Model 3 - FSH 2006-2007 as a function of repro. category
model_3 <- fit_categorical_FSH_model_2006_2007(d=FSH_data_2006_2007)

# Model 4 - estradiol 2016-2018 as a function of repro. category
model_4 <- fit_categorical_estradiol_model(d=estradiol_data)

# Model 5 - estrone 2016-2018 as a function of repro. category
model_5 <- fit_categorical_estrone_model(d=estrone_data)

# Model 6 - pregnanediol 2016-2018 as a function of repro. category
model_6 <- fit_categorical_pregnanediol_model(d=pregnanediol_data)

# Model 7 - log LH 2016-2017 null linear for model comparison
model_7 <- fit_age_linear_null_log_LH(LH_data)

# Model 8 - FSH 2016-2017 null linear for model comparison
model_8 <- fit_age_linear_null_log_FSH_2016_2018(FSH_data_2016_2018)

# Model 9 - FSH 2006-2007 null linear for model comparison
model_9 <- fit_age_linear_null_log_FSH_2006_2007(FSH_data_2006_2007)

# Model 10 - Log LH 2016-2018 as a function of smoothed age
model_10 <- fit_age_gam_log_LH(LH_data)

# Model 11 - Log FSH 2016-2017 as a function of smoothed age
model_11 <- fit_age_gam_log_FSH_2016_2018(d=FSH_data_2016_2018)

# Model 12 - Log FSH 2006-2007 as a function of smoothed age
model_12 <- fit_age_gam_log_FSH_2006_2007(FSH_data_2006_2007)

# Model 13 - Log estradiol as a function of smoothed age
model_13<- fit_age_gam_log_estradiol(estradiol_data)

# Model 14 - Log estrone as a function of smoothed age
model_14<- fit_age_gam_log_estrone(estrone_data)

# Model 15 - Log pregnanediol as a function of smoothed age
model_15<- fit_age_gam_log_pregnanediol(pregnanediol_data)
