library(tseries)
library(forecast)

normalizeQuant <- function (actual, pred, quant){
  sumQuant = 0
  for (i in 1:length(actual)){
    if (pred[i] > actual[i]){
      multiple = -(1-quant)
    }
    else{
      multiple = quant
    }
    QuantRes= 2*(pred[i]-actual[i])*multiple
    sumQuant = sumQuant + QuantRes
  }
  QuantLoss = sumQuant/sum(actual)
  return(QuantLoss)
}

#import data from a csv file
data <- read.csv("LD2011_2014.txt", sep =";",header = TRUE, dec =",", stringsAsFactors = FALSE)



#convert all data except dates into numeric class
cols <- c(-1)
nclient <- 370
#convert all countss except dates into numeric
data[,cols] <- lapply(data[,cols], as.numeric)

#position of important data points
train.end <- 138432

#size of hourly observations
shortTrain <- 24*7
longTrain <- 24*7*2
train.start <- train.end - (longTrain*4) +1

train.size <- (train.end - train.start +1)/4
#size for hourly observation
test1d.size <- 24
test7d.size <- 24*7

test.size <- 24*7

#total number of data points to take to create hourly observation
data.sample.size <- (train.size+test.size)*4


#extract useful data from dataframe
dt <- data[train.start:(train.start+data.sample.size -1),]

#delete the column representing dates
dt[,1] <- NULL


hourly <- data.frame("time.hourly" = c(1:(data.sample.size/4)))

test.start <- train.end + 1


#create a data frame that stores sum of loss result of all clients
results <- data.frame("ARIMA-Q50.1step" = seq(1,370), 
                      "ARIMA-Q50.7step" =seq(1,370), 
                      "ARIMA-Q90.1step" = seq(1,370), 
                      "ARIMA-Q90.7step" = seq(1,370),
                      "ETS-Q50.1step" =seq(1,370), 
                      "ETS-Q50.7step" = seq(1,370), 
                      "ETS-Q90.1step" =seq(1,370), 
                      "ETS-Q90.7step" =seq(1,370))

counts <- 0
for (client in 1:370){
  nr <- 1
  for (i in seq(1,data.sample.size, by = 4))
  {
    #sum of 4 rows
    sum.hourly = sum(dt[i:(i+3),client])
    #insert into data frame "hourly"
    hourly[nr,client] = sum.hourly
    nr <- nr+1
  }
  print(client)
  
  #check whether the row contains only zeros
  sumrow = sum(hourly[1:(data.sample.size/4),client])
  if (sumrow == 0){
    next
  }
  else{
    
      #convert training data to time series format
      tshourly_short <- ts(hourly[(train.size - shortTrain +1):(data.sample.size/4),client], frequency = 7)
      tshourly_long <- ts(hourly[1:(data.sample.size/4),client], frequency = 7)
      
      #actual values of the forecast
      actual1d <- tshourly_short[(shortTrain +1):(shortTrain + test1d.size)]
      actual7d <-  tshourly_long[(longTrain+1):(longTrain + test7d.size)]
      
      #to avoid infinite loss result, eliminate actual values containing only zeros
      if(sum(actual1d)==0){
        next
      }
      
      else{
        
        counts <- counts+1 
        #find the first non zero value of the column
        rowNum_short<- 1
        rowVal <- tshourly_short[rowNum_short]
        while (rowVal ==0){
          rowNum_short = rowNum_short+1
          rowVal = tshourly_short[rowNum_short]
        }
        
        
        rowNum_long <- 1
        rowVal <- tshourly_long[rowNum_long]
        while (rowVal ==0){
          rowNum_long = rowNum_long+1
          rowVal = tshourly_long[rowNum_long]
        }
        
        
        ########################## Forecast with ARIMA######################
        
        ######ROLLING WINDOW FORECAST __ 1day forecast over 7 days
        #the length of training set remains the same; 
        #taking newly obtained forecast values as the newest data to the training set while deleting a same number of oldest data from the training set
        
        ARIMAfit_short <- auto.arima(tshourly_short[rowNum_short:shortTrain])
        #1 day forecast
        pd_s <- forecast(ARIMAfit_short, h = test1d.size)
        #quantile forecast results for 50 percentile
        pd50 <- pd_s$mean
        #quantile forecast results for 90 percentile = upper 80%
        pd90 <- pd_s$upper[,1]
        
        t.start = rowNum_short+24
        half1 <- tshourly_short[t.start:shortTrain]
        half2 <-  pd50
        
        #use rolling-day forecast for 7 days
        while (t.start < 169 ) {
          
          #forecast the vector c(half1,half2)
          ARIMAfit_short <- auto.arima(c(half,1,half2))
          
          #1 day forecast
          pd_s <- forecast(ARIMAfit_short, h = test1d.size)
          #quantile forecast results for 50 percentile
          pd50 <- pd_s$mean
          
          #quantile forecast results for 90 percentile = upper 80%
          pd90 <- pd_s$upper[,1]
          
          t.start <- t.start + 24
          half1 <- tshourly_short[t.start:shortTrain]
          half2 <- c(half2,pd50)
          
        }
        #calculate loss result with actual values and the last round of rolling-day forecast
        results[counts,1] <- abs(normalizeQuant(actual1d,pd50,0.5))
        results[counts,2] <- abs(normalizeQuant(actual1d,pd90,0.9))
        
        
        
        #########ROLLING WINDOW FORECAST __ 7days forecast
        
        
        #forecast directly for 7 days
        ARIMAfit_long <- auto.arima(tshourly_long[rowNum_long:train.size])
        #7 days forecast
        pd_l <- forecast(ARIMAfit_long, h = test7d.size)
        
        #quantile forecast results for 50 percentile
        pd50_l <- pd_l$mean
        #quantile forecast results for 90 percentile = upper 80%
        pd90_l <- pd_l$upper[,1]
        
        results[counts,3] <- abs(normalizeQuant(actual7d,pd50_l,0.5))
        results[counts,4] <- abs(normalizeQuant(actual7d,pd90_l,0.9))
      
        ########################## Forecast with ETS ########################
        
        #the length of training set remains the same; 
        #taking newly obtained forecast values as the newest data to the training set while deleting a same number of oldest data from the training set
        ETS_short <- ets(tshourly_short[rowNum_short:shortTrain])
        #1 day forecast
        pdETS_s <- forecast(ETS_short, h = test1d.size)
        #quantile forecast results for 50 percentile
        pdETS50 <- pdETS_s$mean
        #quantile forecast results for 90 percentile = upper 80%
        pdETS90 <- pdETS_s$upper[,1]
        
        
        
        t.start = rowNum_short +24
        half1 <- tshourly_short[t.start:shortTrain]
        half2 <- pdETS50
        
        while (t.start < 169 ) {
        
          ETSfit_short <- ets(c(half1,half2))
          
          #1 day forecast
          pdETS_s <- forecast(ETSfit_short, h = test1d.size)
          #quantile forecast results for 50 percentile
          pdETS50 <- pdETS_s$mean
          
          #quantile forecast results for 90 percentile = upper 80%
          pdETS90 <- pdETS_s$upper[,1]
          
          t.start <- t.start + 24
          half1 <- tshourly_short[t.start:shortTrain]
          half2 <- c(half2, pdETS50)
        }
        
        #calculate loss result with actual values and the last round of rolling-day forecast
        results[counts,5] <- abs(normalizeQuant(actual1d,pdETS50,0.5))
        results[counts,6] <- abs(normalizeQuant(actual1d,pdETS90,0.9))
        
        ############7 days forecast
        ETSfit_long <- ets(tshourly_long[rowNum_long:train.size])
        pdETS <- forecast(ETSfit_long, h = test7d.size)
        #quantile forecast results for 50 percentile
        pdETS50_l <- pdETS$mean
        #quantile forecast results for 90 percentile = upper 80%
        pdETS90_l <- pdETS$upper[,1]
        
        results[counts,7] <- abs(normalizeQuant(actual7d,pdETS50_l,0.5))
        results[counts,8] <- abs(normalizeQuant(actual7d,pdETS90_l,0.9))
        
      }
      
      
      
    
  }
}


comparison <- data.frame("ARIMA-1day-0.5" = c(mean(results[1:counts,1]),0.154),
                         "ARIMA-7day-0.5" = c(mean(results[1:counts,3]),0.283),
                         "ARIMA-1day-0.9" = c(mean(results[1:counts,2]),0.102),
                         "ARIMA-7day-0.9" = c(mean(results[1:counts,4]),0.109),
                         "ETS-1day-0.5" = c(mean(results[1:counts,5]),0.101),
                         "ETS-7day-0.5" = c(mean(results[1:counts,7]),0.121),
                         "ETS-1day-0.9" = c(mean(results[1:counts,6]),0.077),
                         "ETS-7day-0.9" = c(mean(results[1:counts,8]),0.101))
rownames(comparison) <- c("experimental", "Observed")