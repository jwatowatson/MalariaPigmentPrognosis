## Outstanding issues
## Multiplication by 2 or 4 for the MONO50 to MONO200 and

ff=list.files(path = 'Data/SMAC/')[1]
my_cols = colnames(read.csv(paste0('Data/SMAC/',ff)))

myMergedData <-
  do.call(rbind,
          lapply(list.files(path = 'Data/SMAC/'), function(x) {
            dat = read.csv(file = paste0('Data/SMAC/',x))
            print(x)
            country = unlist(strsplit(x, split = 'pig'))[1]
            dat = data.frame(dat[, my_cols])
            dat = cbind(dat, country=rep(country, nrow(dat)))
            }))
table(myMergedData$SITENO)
table(myMergedData$country)


## Make missing data NA
myMergedData$NEUTRO100[myMergedData$NEUTRO100==999]=NA
myMergedData$NEUTROUL[myMergedData$NEUTROUL>=999]=NA
myMergedData$POLYL200[myMergedData$POLYL200>200]=NA
myMergedData$MONOL200[myMergedData$MONOL200==999]=NA
myMergedData$MONO50[myMergedData$MONO50>50]=NA
myMergedData$MONOUL[myMergedData$MONOUL==9999]=NA
myMergedData$OUTCOME[myMergedData$OUTCOME>2]=NA
myMergedData$OUTCOME[myMergedData$OUTCOME==1]=0
myMergedData$OUTCOME[myMergedData$OUTCOME==2]=1
myMergedData$HB[myMergedData$HB==99.9]=NA
myMergedData$HCT[myMergedData$HCT==99.9]=NA
myMergedData$TEMP[myMergedData$TEMP==99.9]=NA
myMergedData$RESPRATE[myMergedData$RESPRATE==999]=NA
myMergedData$PH[myMergedData$PH==9.999]=NA
myMergedData$PCO2[myMergedData$PCO2==99.9]=NA
myMergedData$BE[myMergedData$BE==999.9]=NA
myMergedData$AGE[myMergedData$AGE==999]=NA
myMergedData$LACTATE[myMergedData$LACTATE==99.9]=NA
myMergedData$BMS[myMergedData$BMS==9]=NA
myMergedData$BVS[myMergedData$BVS==9]=NA
myMergedData$BES[myMergedData$BES==9]=NA
myMergedData$AGE[myMergedData$AGE>180]=NA
myMergedData$PARASIT[myMergedData$PARASIT==9999999]=NA


apply(myMergedData, 2, function(x) round(100*mean(!is.na(x))))

# parasite data to check scale
hist(log10(myMergedData$PARASIT+50))

ind = is.na(myMergedData$MONOL200) & !is.na(myMergedData$MONO50)
sum(ind)
myMergedData$MONOL200[ind] = myMergedData$MONO50[ind] #*4

ind = is.na(myMergedData$POLYL200) & !is.na(myMergedData$NEUTRO100)
sum(ind)
myMergedData$POLYL200[ind] = myMergedData$NEUTRO100[ind] #*2


table(myMergedData$POLYL200>100, myMergedData$country)
table(myMergedData$MONOL200>100, myMergedData$country)

myMergedData$HB[myMergedData$HB<1]=NA
plot(myMergedData$HB, myMergedData$HCT, pch='.')
abline(v=15, h=45)
myMergedData$HB[myMergedData$HB>15] = NA
myMergedData$HCT[myMergedData$HCT>45] = NA
hist(myMergedData$HCT/myMergedData$HB, breaks = 100)
ind = !is.na(myMergedData$HCT/myMergedData$HB) & myMergedData$HCT/myMergedData$HB>4
table(myMergedData$country[ind])
myMergedData$HCT[ind] = NA
myMergedData$HB[ind] = NA
hist(myMergedData$HCT/myMergedData$HB, breaks = 100)
ind = !is.na(myMergedData$HCT/myMergedData$HB) & myMergedData$HCT/myMergedData$HB<2
myMergedData$HCT[ind] = NA
myMergedData$HB[ind] = NA
hist(myMergedData$HCT/myMergedData$HB, breaks = 100)
myMergedData$mycols = RColorBrewer::brewer.pal(n = 6,name = 'Dark2')[as.numeric(as.factor(myMergedData$country))]
plot(myMergedData$HB, myMergedData$HCT, pch=20,
     col = myMergedData$mycols)
lines(0:15, 3*(0:15), col=1, lwd=3, lty=2)
lines(0:15, 3*(0:15)-5, col=1, lwd=3, lty=2)
sum(!is.na(myMergedData$HB/myMergedData$HCT))
legend('topleft', col = unique(myMergedData$mycols),
       legend = unique(myMergedData$country), pch=20)
title(paste('n=',sum(!is.na(myMergedData$HB/myMergedData$HCT)),sep = ''))

table(myMergedData$country[!is.na(myMergedData$HCT/myMergedData$HB)])


ind = !is.na(myMergedData$HCT/myMergedData$HB) & myMergedData$HB<=5&
  myMergedData$HB>1
table(myMergedData$country[ind])

myMergedData$country_names=plyr::mapvalues(myMergedData$country,
                from = c("gam" ,"gha", "ken" ,"lam" ,"lib", "mal"),
                to = c('The Gambia','Ghana','Kenya',
                       'Gabon (Lambarene)', 'Gabon (Libreville)','Malawi'))
par(las=1, family='serif')
plot(jitter(myMergedData$HB[ind],amount = .025),
     myMergedData$HCT[ind], pch=20,
     col = myMergedData$mycols[ind], panel.first = grid(),
     xlab='Haemoglobin (g/dL)', ylab='Haematocrit (%)')
legend('topleft', col = unique(myMergedData$mycols),inset = 0.03,
       legend = unique(myMergedData$country_names), pch=20)
title(paste('n=',sum(ind),sep = ''))

mod = (MASS::rlm(HCT ~ HB, myMergedData[ind, ]))
lines(1:5, predict(mod, data.frame(HB=1:5)),lwd=3, lty=2)
lines(0:15, 3*(0:15), col=1, lwd=3, lty=1)



ind = !is.na(myMergedData$HCT/myMergedData$HB) & myMergedData$HB>5
table(myMergedData$country[ind])

plot(myMergedData$HB[ind], myMergedData$HCT[ind], pch=20,
     col = myMergedData$mycols[ind], panel.first = grid(),
     xlab='Hb', ylab='HCT')
legend('topleft', col = unique(myMergedData$mycols),
       legend = unique(myMergedData$country), pch=20)
title(paste('n=',sum(ind),sep = ''))

mod = (MASS::rlm(HCT ~ HB, myMergedData[ind, ]))
lines(1:50, predict(mod, data.frame(HB=1:50)),lwd=3, lty=2)
lines(0:15, 3*(0:15), col=1, lwd=3, lty=1)
mod = (lm(HCT ~ HB, myMergedData[ind, ]))
lines(1:50, predict(mod, data.frame(HB=1:50)),lwd=3, lty=2,col='red')

aggregate(HCT/HB ~ country, myMergedData, mean)
aggregate(HCT/HB ~ country, myMergedData, sd)


sum(is.na(myMergedData$HB) & is.na(myMergedData$HCT))
ind = is.na(myMergedData$HB) & !is.na(myMergedData$HCT)
myMergedData$HB[ind] = (myMergedData$HCT[ind]-0.5)/3
# myMergedData = dplyr::filter(myMergedData,
#                              !(is.na(POLYL200)&
#                                  is.na(MONOL200)),
#                              !is.na(OUTCOME))

write.csv(x = myMergedData, file = '~/Downloads/compiled_Kremsner_data.csv')
save(myMergedData, file = 'RData/SMAC_data.RData')
