rm(list=ls())


# Load libraries - when fail to import, please run: install.packages("library-name")
library(readxl)
library(e1071)


# Define functions
colApply <- function(dat, cols = colnames(dat), func = as.factor) {
  dat[cols] <- lapply(dat[cols], func)
  return(dat)
}



mutual.information <- function(x, z) {
  p00 <- sum(x == 0 & z == 0)/length(z)
  p01 <- sum(x == 0 & z == 1)/length(z)
  p10 <- sum(x == 1 & z == 0)/length(z)
  p11 <- sum(x == 1 & z == 1)/length(z)
  p1. <- sum(x == 1)/length(z)
  p0. <- sum(x == 0)/length(z)
  p.1 <- sum(z == 1)/length(z)
  p.0 <- sum(z == 0)/length(z)
  out <- 0
  if (p00 > 0)
    out <- out + p00 * log(p00/p0./p.0)
  if (p01 > 0)
    out <- out + p01 * log(p01/p0./p.1)
  if (p10 > 0)
    out <- out + p10 * log(p10/p1./p.0)
  if (p11 > 0)
    out <- out + p11 * log(p11/p1./p.1)
  out
}



composite.xvalid <- function(keys, 
                             train, 
                             test, 
                             start.t=start.t, 
                             end.t=end.t,  
                             engraftment.aware=engraftment.aware, 
                             use.smoothing=TRUE, 
                             use.laplace=TRUE,
                             use.stationary=TRUE
) {
  
  
  if (FALSE) {
    "
      yy & y would be seen in pair in this fuction.
      yy stands for negative case 
      y stands for positive case
      "
  }
  
  
  ## Data manipulation
  all.patients <- complete.patients
  all.patient.data <- complete.patient.data
  all.severe <- complete.severe
  
  patients <- all.patients[train]
  severe <- intersect(all.severe, all.patients[train])
  
  
  ## Build Model - thresholds
  for (i in 1:length(train)) { 
    m <- train[i]
    if (df$event[m]==1) {
      onset.day <- df$onset.day[m]
      all.patient.data[[m]][1 : min(onset.day - post.window.size - 1, end.t) + offset.t, ] <- NA 
      
    }
  }
  
  
  thresholds <- array(NA, dim = c(day, length(keys)))
  for (k in 1:length(keys)) {
    key <- keys[k]
    
    yy <- c()
    for (i in 1:length(patients)) {
      patient <- patients[i]
      if (!(patient %in% severe)) {
        m <- (1:length(all.patients))[all.patients == patient]
        patient.data <- all.patient.data[[m]]
        j <- (1:ncol(patient.data))[names(patient.data) == key]
        yy <- cbind(yy, patient.data[, j])
      }
    }
    
    y <- c()
    for (i in 1:length(severe)) {
      patient <- severe[i]
      m <- (1:length(all.patients))[all.patients == patient]
      patient.data <- all.patient.data[[m]]
      j <- (1:ncol(patient.data))[names(patient.data) == key]
      y <- cbind(y, patient.data[, j])
    }
    
    
    for (t in start.t:end.t + offset.t) {
      if (sum(!is.na(y[t, ])) > 1 & sum(!is.na(yy[t, ])) > 1) {
        min.val <- min(y[t, ], yy[t, ], na.rm = T)
        max.val <- max(y[t, ], yy[t, ], na.rm = T)
        
        vals <- seq(from = min.val, to = max.val, length.out = 100)
        shannon <- array(NA, length(vals))
        for (iiii in 1:length(vals)) {
          thres <- vals[iiii]
          x <- as.numeric(c(y[t, ], yy[t, ]) > thres)
          z <- c(array(1, ncol(y)), array(0, ncol(yy)))
          I <- !is.na(x)
          x <- x[I]
          z <- z[I]
          shannon[iiii] <- mutual.information(x, z)
        }
        
        best.iiii <- (1:length(shannon))[shannon == max(shannon)][1]
        
        thres <- vals[best.iiii]
        p.y.over <- sum(y[t, ] > thres, na.rm = T)/sum(!is.na(y[t, ]))
        p.yy.over <- sum(yy[t, ] > thres, na.rm = T)/sum(!is.na(yy[t, ]))
        
        thresholds[t, k] <- thres
      }
      
    }
  }
  
  
  thresholds.high <- array(NA, dim = c(day, length(keys)))
  thresholds.low <- array(NA, dim = c(day, length(keys)))
  for (k in 1:length(keys)) {
    thresholds.high[, k] <- quantile(thresholds[, k], 3/4, na.rm = T)
    thresholds.low[, k] <- quantile(thresholds[, k], 1/4, na.rm = T)
  }
  
  
  ## Build Model - odds ratio
  p.y.high <- array(NA, dim = c(day, length(keys)))
  p.y.low <- array(NA, dim = c(day, length(keys)))
  p.yy.high <- array(NA, dim = c(day, length(keys)))
  p.yy.low <- array(NA, dim = c(day, length(keys)))
  for (k in 1:length(keys)) {
    key <- keys[k]
    
    yy <- c()
    for (i in 1:length(patients)) {
      patient <- patients[i]
      if (!(patient %in% severe)) {
        m <- (1:length(all.patients))[all.patients == patient]
        patient.data <- all.patient.data[[m]]
        j <- (1:ncol(patient.data))[names(patient.data) == key]
        yy <- cbind(yy, patient.data[, j])
      }
    }
    
    y <- c()
    for (i in 1:length(severe)) {
      patient <- severe[i]
      m <- (1:length(all.patients))[all.patients == patient]
      patient.data <- all.patient.data[[m]]
      j <- (1:ncol(patient.data))[names(patient.data) == key]
      y <- cbind(y, patient.data[, j])
    }
    
    for (t in start.t:end.t + offset.t) {
      if (!is.na(thresholds.high[t, k]) & sum(!is.na(y[t, ])) > 1 & sum(!is.na(yy[t, ])) > 1) {
        thres <- thresholds.high[t, k]
        p.y.high[t, k] <- sum(y[t, ] > thres, na.rm = T)/sum(!is.na(y[t, ]))
        p.yy.high[t, k] <- sum(yy[t, ] > thres, na.rm = T)/sum(!is.na(yy[t, ]))
      }
      if (!is.na(thresholds.low[t, k]) & sum(!is.na(y[t, ])) > 1 & sum(!is.na(yy[t, ])) > 1) {
        thres <- thresholds.low[t, k]
        p.y.low[t, k] <- sum(y[t, ] < thres, na.rm = T)/sum(!is.na(y[t, ]))
        p.yy.low[t, k] <- sum(yy[t, ] < thres, na.rm = T)/sum(!is.na(yy[t, ]))
      }
    }
  }
  
  
  odds.ratio.high <- array(NA, dim = c(day, length(keys)))
  odds.ratio.medium <- array(NA, dim = c(day, length(keys)))
  odds.ratio.low <- array(NA, dim = c(day, length(keys)))
  for (k in 1:length(keys)) {
    
    if (use.smoothing) {
      x <- start.t:end.t + offset.t
      y <- p.y.high[start.t:end.t + offset.t, k]
      I <- !is.na(y)
      if (sum(I) > 4) {
        model <- smooth.spline(x[I], y[I], all.knots = T, df = 4)
        p.y.high[start.t:end.t + offset.t, k] <- predict(model, start.t:end.t + offset.t)$y
        p.y.high[start.t:end.t + offset.t, k][p.y.high[start.t:end.t + offset.t, k] < 0] <- 0
      } else {
        p.y.high[start.t:end.t + offset.t, k] <- 0
      }
      
      x <- start.t:end.t + offset.t
      y <- p.yy.high[start.t:end.t + offset.t, k]
      I <- !is.na(y)
      if (sum(I) > 4) {
        model <- smooth.spline(x[I], y[I], all.knots = T, df = 4)
        p.yy.high[start.t:end.t + offset.t, k] <- predict(model, start.t:end.t + offset.t)$y
        p.yy.high[start.t:end.t + offset.t, k][p.yy.high[start.t:end.t + offset.t, k] < 0] <- 0
      } else {
        p.yy.high[start.t:end.t + offset.t, k] <- 0
      }
      
      x <- start.t:end.t + offset.t
      y <- p.y.low[start.t:end.t + offset.t, k]
      I <- !is.na(y)
      if (sum(I) > 4) {
        model <- smooth.spline(x[I], y[I], all.knots = T, df = 4)
        p.y.low[start.t:end.t + offset.t, k] <- predict(model, start.t:end.t + offset.t)$y
        p.y.low[start.t:end.t + offset.t, k][p.y.low[start.t:end.t + offset.t, k] < 0] <- 0
      } else {
        p.y.low[start.t:end.t + offset.t, k] <- 0
      }
      
      x <- start.t:end.t + offset.t
      y <- p.yy.low[start.t:end.t + offset.t, k]
      I <- !is.na(y)
      if (sum(I) > 4) {
        model <- smooth.spline(x[I], y[I], all.knots = T, df = 4)
        p.yy.low[start.t:end.t + offset.t, k] <- predict(model, start.t:end.t + offset.t)$y
        p.yy.low[start.t:end.t + offset.t, k][p.yy.low[start.t:end.t + offset.t, k] < 0] <- 0
      } else {
        p.yy.low[start.t:end.t + offset.t, k] <- 0
      }
    }
    
    if (use.laplace) {
      odds.ratio.high[, k] <- (p.y.high[, k] + 0.1)/(p.yy.high[, k] + 0.1)
      odds.ratio.low[, k] <- (p.y.low[, k] + 0.1)/(p.yy.low[, k] + 0.1)
      
      y.high.and.low <- p.y.high[, k] + p.y.low[, k]
      y.high.and.low[y.high.and.low > 1] <- 1
      
      yy.high.and.low <- p.yy.high[, k] + p.yy.low[, k]
      yy.high.and.low[yy.high.and.low > 1] <- 1
      
      odds.ratio.medium[, k] <- (1 - y.high.and.low + 0.1)/(1 - yy.high.and.low + 0.1)
    } else {
      odds.ratio.high[, k] <- (p.y.high[, k] + 0)/(p.yy.high[, k] + 0)
      odds.ratio.low[, k] <- (p.y.low[, k] + 0)/(p.yy.low[, k] + 0)
      
      odds.ratio.medium[, k] <- (1 - p.y.high[, k] - p.y.low[, k] + 0)/(1 - p.yy.high[, k] - p.yy.low[, k] + 0)
    }
    
    
    odds.ratio.high[start.t:end.t + offset.t, ][is.na(odds.ratio.high[start.t:end.t + offset.t, ])] <- 1
    odds.ratio.low[start.t:end.t + offset.t, ][is.na(odds.ratio.low[start.t:end.t + offset.t, ])] <- 1
    odds.ratio.medium[start.t:end.t + offset.t, ][is.na(odds.ratio.medium[start.t:end.t + offset.t, ])] <- 1
    
  }
  
  
  ## Test Model
  trials <- test
  
  out.odds <- array(NA, dim = c(day, length(trials)))
  is.severe <- array(NA, length(trials))
  data.points <- array(NA, length(trials))
  real.st <- array(NA, length(trials))
  real.et <- array(NA, length(trials))
  
  
  for (i in 1:length(trials)) {
    
    trial <- trials[i]
    
    patient <- all.patients[trial]
    m <- (1:length(all.patients))[all.patients == patient]
    patient.data <- all.patient.data[[m]]
    
    
    if (engraftment.aware) {
      real.start.t <- max(df$`Neutrophil.engraftment.days`[df$ID == patient], start.t)
    } else {
      real.start.t <- start.t
    }
    
    real.st[i] <- real.start.t
    
    real.end.t <- min(df$death.day[df$ID == patient], end.t)
    real.et[i] <- real.end.t
    
    
    odds <- array(NA, day) 
    if (sum(!is.na(patient.data)) > 0) {
      
      for (t in real.start.t:real.end.t) {
        
        if(sum(!is.na(patient.data[t + offset.t, ])) > 0){
          odds[t + offset.t] <- 0   
          for (k in 1:length(keys)) {
            key <- keys[k]
            val <- patient.data[t + offset.t, names(patient.data) == key]
            if (!is.na(val) & !is.na(thresholds.high[t + offset.t, k])) {
              if (val > thresholds.high[t + offset.t, k]) {
                odds[t + offset.t] <- odds[t + offset.t] + log(odds.ratio.high[t + offset.t, k])
              }
            }
            if (!is.na(val) & !is.na(thresholds.low[t + offset.t, k])) {
              if (val < thresholds.low[t + offset.t, k]) {
                odds[t + offset.t] <- odds[t + offset.t] + log(odds.ratio.low[t + offset.t, k])
              }
            }
            if (!is.na(val) & !is.na(thresholds.low[t + offset.t, k]) & !is.na(thresholds.high[t + offset.t, k])) {
              if (val >= thresholds.low[t + offset.t, k] & val <= thresholds.high[t + offset.t, k]) {
                odds[t + offset.t] <- odds[t + offset.t] + log(odds.ratio.medium[t + offset.t, k])
              }
            }
          }
        }
        
      }
      
      
      
      out.odds[, i] <- odds
      
      if (patient %in% all.severe) {
        is.severe[i] <- 1
      } else {
        is.severe[i] <- 0
      }
      
      data.points[i] <- data.density[all.patients == patient]
      
    }
    
  }
  
  
  
  ## Turn on use.stationary mode
  if (use.stationary) {
    
    severe.I.raw <- intersect(big.severe.I, train)
    mild.I.raw <- intersect(big.mild.I, train)
    
    
    for (t in start.t:end.t) {
      
      mild.I <- mild.I.raw

      severe.I <- severe.I.raw

      train.data <- subset(df, ID %in% complete.patients[c(mild.I, severe.I)])
      
      idx <- c("Age", "Sex")  
      
      model.nb <- naiveBayes(x=train.data[, idx],   
                             y=train.data$event, 
                             laplace=10)
      

      
      for (i in 1:length(trials)) {
        
        trial <- trials[i]
        
        patient <- all.patients[trial]
        m <- (1:length(all.patients))[all.patients == patient]
        patient.data <- df[m, ]

          
        prob.nb <- predict(model.nb, patient.data[, idx], type = c("raw"),)
        
        out.odds[t + offset.t, i] <- sum(c(0,
                                           log( (prob.nb[2]+0.01) / (prob.nb[1]+0.01) )
        ), na.rm=T)
        
      }
      
    }
    
  }
  
  
  
  ## Output
  out <- list(is.severe, out.odds, data.points, 
              odds.ratio.high, odds.ratio.medium, odds.ratio.low)
  
  out
  
}



engraftment.aware <- 0

start.t <- 0
end.t <- 25
post.window.size <- 14
offset.t <- 15
day <- end.t + offset.t
nfold <- 5 
grade <- '234'

root <- './mock_up_data/'     ######the absolute path of the mock_up_data folder under your project directory

encoded.name <- "ITP_encoded.csv"
keys.name <- "ITP_keys.csv"


keys <- read.csv(paste(root, keys.name, sep = ""))$keys


factor.col <- c("Confirmed_response", "Sex", "Death")


df <- read.csv(paste(root, encoded.name, sep = "")) 


df <- colApply(df, factor.col)


complete.patients <- df$ID
death.or.not <- df$Death
death.date <- df$Last.followup.time.for.death

df[!is.na(df$Confirmed_response_onset_day) & df$Confirmed_response_onset_day > 5, ]$Confirmed_response <- 0
df$event <- df$Confirmed_response 
df$onset.day <- df$Confirmed_response_onset_day 
df[df$event==0, "onset.day"] <- Inf


df$death.day <- df$Last.followup.time.for.death
df[df$Death==0, "death.day"] <- Inf

event <- df$event
fold <- c(df$fold)

complete.severe <- df$ID[df$event == 1]
complete.mild <- df$ID[df$event == 0]


data.density <- c()
complete.patient.data <- list()
for (i in 1:length(complete.patients)) {
  
  patient <- complete.patients[i]
  file.name <- paste(root, "patients/ITP_pt_", patient, ".csv", sep = "")
  patient.data <- read.csv(file.name, check.names = FALSE)
  patient.data <- patient.data[1:day, 2:ncol(patient.data)]
  
  # time-limited sample-and-hold
  for (w in 1:3) {
    hasValue <- list()
    for (t in offset.t:day) {
      hasValue[[t]] <- !is.na(patient.data[t, ])
    }
    
    for (t in day:offset.t) {
      patient.data[t, !hasValue[[t]]] <- patient.data[t - 1, !hasValue[[t]]]
    }
  }
  
  complete.patient.data[[i]] <- patient.data
  
  
  j <- (1:nrow(df))[df$ID == patient]
  
  complete.patient.data[[i]] <- complete.patient.data[[i]][1:day, ]
  data.density[i] <- sum(!is.na(complete.patient.data[[i]][, keys]))
  
  
}


big.severe.I <- c()
for (i in 1:length(complete.severe)) {
  big.severe.I <- c(big.severe.I, (1:length(complete.patients))[complete.patients == complete.severe[i]])
}
big.severe.I <- big.severe.I[!is.na(big.severe.I)]
big.mild.I <- setdiff(1:length(complete.patients), big.severe.I)


test.set <- (1:length(complete.patients))[df$fold=='Holdout2']
severe.I <- setdiff(big.severe.I, test.set)
mild.I <- setdiff(big.mild.I, test.set)

out_stationary <- composite.xvalid(keys=keys,
                        train=c(severe.I, mild.I),
                        test=(1:length(complete.patients)),
                        start.t=start.t,
                        end.t=end.t,
                        engraftment.aware=engraftment.aware,
                        use.smoothing=TRUE,
                        use.laplace=TRUE,
                        use.stationary=TRUE)

odds_stationary <- as.data.frame(out_stationary[[2]])
colnames(odds_stationary) <- complete.patients


out_dynamic <- composite.xvalid(keys=keys,  
                        train=c(severe.I, mild.I),
                        test=(1:length(complete.patients)),
                        start.t=start.t,
                        end.t=end.t,
                        engraftment.aware=engraftment.aware,
                        use.smoothing=TRUE,
                        use.laplace=TRUE,
                        use.stationary=FALSE)

odds_dynamic <- as.data.frame(out_dynamic[[2]])
colnames(odds_dynamic) <- complete.patients


if(sum(!(colnames(odds_stationary) == colnames(odds_dynamic))) != 0){
  stop("Error!")
}

add_elements <- function(x, y) {
  ifelse(is.na(x) & is.na(y), NA, ifelse(is.na(x), y, ifelse(is.na(y), x, x + y)))
}


eva.odds <- as.data.frame(mapply(add_elements, odds_stationary, odds_dynamic, SIMPLIFY = FALSE))


colnames(eva.odds) <- sub("^X", "", colnames(eva.odds))


outcome <- 'Confirmed_response'
start.day <- 1
end.day <- 5
ratio.low <- 0.5


day.seq <- seq(start.day, end.day)

time.colname <- "Confirmed_response_onset_day" 
status.colname <- "Confirmed_response"


data.II.IV <- t(eva.odds)[, -(1:15)]
thr.low.num.II.IV <- vector('numeric', length = length(day.seq))

for (i in 1:length(day.seq)) {
  day <- day.seq[i]
  
  pre.value.II.IV <- data.II.IV[1 : 85, day]
  pre.value.II.IV.temp <- pre.value.II.IV
  thr.low.II.IV <- quantile(pre.value.II.IV.temp, ratio.low, na.rm = T)         
  
  thr.low.num.II.IV[i] <- thr.low.II.IV
  
  
}

ind.list <- vector('list', length = length(day.seq))

for (i in 1:length(day.seq)) {
  day <- day.seq[i]

  pre.value.II.IV <- data.II.IV[, day]
  
  thr.low.II.IV <- thr.low.num.II.IV[i]
  
  I.high.II.IV <- pre.value.II.IV >= thr.low.II.IV
  I.low.II.IV <- pre.value.II.IV < thr.low.II.IV 
  
  ind.list[[i]] <- I.low.II.IV
  ind.list[[i]][I.high.II.IV & !is.na(I.high.II.IV)] <- 'H'
  ind.list[[i]][I.low.II.IV & !is.na(I.low.II.IV)] <- 'L'
  
}
ind.data <- as.matrix(as.data.frame(ind.list))
colnames(ind.data) <- 1:ncol(ind.data)


risk <- apply(ind.data, 1, function(x){
  if (sum(!is.na(x)) == 0) {
    NA
  } else {
    if(sum(x == "H") >= 4 & sum(x[3 : 5] == "H") == 3){
      "High-confidence"
    }else{
      "Low-confidence"
    }
  }
  
})

df_risk <- data.frame(ID = complete.patients, risk = risk, plt_30_time = NA)

plt_id <- which(colnames(complete.patient.data[[1]]) == "cbc_PLT")


for (i in 1:length(complete.patients)) {
  
  patient <- complete.patients[i]
  file.name <- paste(root, "patients/ITP_pt_", patient, ".csv", sep = "")
  patient.data <- read.csv(file.name, check.names = FALSE)
  patient.data <- patient.data[, 2:ncol(patient.data)]
  
  plt.data <- patient.data[-c(1 : 15), plt_id]
  if(min(which(plt.data >= 30)) <= 5){
    df_risk[i, "plt_30_time"] <- "fast"
  }
  else{
    df_risk[i, "plt_30_time"] <- "slow/no response"
  }
}


df_risk$risk_bind <- NA
df_risk[df_risk$plt_30_time == "slow/no response", ]$risk_bind <- "Slow/no response"
df_risk[df_risk$plt_30_time == "fast", ]$risk_bind <- paste(df_risk[df_risk$plt_30_time == "fast", ]$risk, 
                                                             df_risk[df_risk$plt_30_time == "fast", ]$plt_30_time)

