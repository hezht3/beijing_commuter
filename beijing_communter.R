
require(dplyr)
require(ggplot2)
require(rcompanion)
require(pairwiseCI)
require(readr)
require(Hmisc)
require(PerformanceAnalytics)
require(chron)
require(lubridate)

#####functions#################
filterMDv4<-function(fileTOfilter){
  #  fileTOfilter<-MDdat2#for testing
  tempdat0<-fileTOfilter
  tempdat<-tempdat0[-1,] #first row (1-sec) removed to match output from 'ratioV2'
  
  ##SIZE RANGE CHANGES (MING 17-06-2015 email)
  validSize<-tempdat$Size>10 & tempdat$Size<300
  validFlow<-tempdat$Flow>=0.9 & tempdat$Flow<=1.1
  #MD error breakdown:
  #the following converts the MDerror value into a binary number with 32 places.
  #this is used to figure out what errors are present by column...
  md<-as.data.frame(do.call(rbind, lapply(tempdat$Error,function(x){as.integer(intToBits(x))})))
  names(md)<-c('MDw1','MDw2','MDw3','MDw4','MDw5','gap','MDe1','MDe2','MDe3','MDe4','MDe5','MDe6','MDe7','MDe8','MDe9','MDe10')
  Errors <- md$MDw1|md$MDw2|md$MDw3|md$MDw4|md$MDw5|md$MDe1|md$MDe2|md$MDe3|md$MDe4|md$MDe5|md$MDe6|md$MDe7|md$MDe8|md$MDe9|md$MDe10 
  notErrors<-!Errors  
  validCount<-tempdat$Number>100
  Factor<-ratioV3(tempdat0$Number,10)  #the # of rows will now match those of other valid* vars.
  validFactor<-!Factor[,4]
  keepIndex<-validSize & validFlow & notErrors & validCount & validFactor
  cleanData<-tempdat[keepIndex,]
  percOf24hrs<-dim(cleanData)[1]/60/60/24*100
  stats<-data.frame(orig_n=length(validSize),total_drop=sum(!keepIndex),size_drop=sum(!validSize),flow_drop=sum(!validFlow),Num_drop=sum(!validCount), Factor_drop=sum(!validFactor), error_drop=sum(Errors),perc_total_drop=sum(!keepIndex)/length(validSize)*100,perc_of_24hrs=percOf24hrs,MDw1=sum(md$MDw1),MDw2=sum(md$MDw2),MDw3=sum(md$MDw3),MDw4=sum(md$MDw4),MDw5=sum(md$MDw5),MDe1=sum(md$MDe1),MDe2=sum(md$MDe2),MDe3=sum(md$MDe3),MDe4=sum(md$MDe4),MDe5=sum(md$MDe5),MDe6=sum(md$MDe6),MDe7=sum(md$MDe7),MDe8=sum(md$MDe8),MDe9=sum(md$MDe9),MDe10=sum(md$MDe10))
  #remove the first row (1-second) of cleanData & stats as
  return(list(cleanData,stats))
  
}




importMD_dm<-function(filename){
  #LOAD MINIDISC DATA INTO R:
  dat<-read.csv(file=filename,header=TRUE,skip=5,sep="\t")
  oldjava<-dim(dat)[2]==16
  #Following removes extra last column that shows up in dat thru javascript:
  if(dim(dat)[2]==16)rmcol<-16  #output file processed with earlier javascript 
  if(dim(dat)[2]==18)rmcol<-18  #output file processed with javascript v1.217.
  dat<-dat[,-rmcol] #remove extra \t column (turns out that somewhere between season 1&2, the javatool provides 2 additional columns (Surface & Mass); therefore, need to remove column 16 before, and column 18 after new columns introduced...
  timemd <- strptime(dat[,1],format="%d-%b-%Y %H:%M:%S") #CONVERT TIMESTAMP TO R time standard
  dat2<-cbind(PosixTime=timemd,dat)
  return(list(dat2,oldjava))
}





ratioV3<-function(v,x){ #(vector v, factor x)
  vectorLen<-length(v)
  fact<-x
  tempdat0<-v[-vectorLen] #drop the last observation
  tempdat1<-v[-1]         #drop the first observation
  ratio<-tempdat1/tempdat0
  
  ratio1<-as.data.frame(cbind(ratio,ratio>fact|ratio<1/fact)) #the length of this vector is 1 less than the input
  
  
  cumsum1<-cumsum(ratio1[,2]) #cumulative sum.
  rle1<-rle(cumsum1) #run length encoding--Compute the lengths and values of runs of equal values
  tossDF<-as.data.frame(cbind(rle1$values,rle1$lengths)) #dataframe with values and their lengths
  names(tossDF)<-c('values','lengths')
  tossDF<-tossDF[-1,] #remove value of '0' from being tossed.
  index<-which(tossDF$lengths==1) #if a value occurs only once, that means that it is surrounded on both sides by values that are a factor of >10 different.
  valsToRm<-tossDF$values[index] #using the index, pick out the cumulative sums to be removed.
  to_remove<-is.element(cumsum1,valsToRm)                     
  
  ratio1<-cbind(ratio1,cumsum1,to_remove)
  names(ratio1)[2]<-paste('GTorLTfact_',fact,sep='')
  return(ratio1)
}



gm_mean<-function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}



setwd("C:/Users/zy125/Box Sync/PhD/Beijing")

##############################################################################
#############load bus data###################################################
###############################################################################

path <- "C:/Users/zy125/Box Sync/PhD/Beijing/bus"
file.names <- dir(path, pattern =".txt")
setwd(path)
bus <- list()

for (f in file.names) {
  print(f)
  MDdat2 <- importMD_dm(f)
  MDdat2filt <- filterMDv4(MDdat2[[1]])
  bus[[f]] <- MDdat2filt[[1]][, c(1,4,5,6,7,8,16)]
}


for (n in 1:18){
  bus[[n]]['id'] = n
}
bus.all = do.call(rbind, bus)

bus_summary<-bus.all %>%
  group_by(id) %>%
  dplyr::summarize(number.amean = mean(Number), number.gmean=gm_mean(Number),number.median=median(Number), 
                   number.max=max(Number), number.min=min(Number), number.sd=sd(Number),number.q5=quantile(Number, probs = 0.05),
                   number.q25=quantile(Number, probs = 0.25), number.q75=quantile(Number, probs = 0.75),
                   number.q95=quantile(Number, probs = 0.95))



bus_summary_size<-bus.all %>%
  group_by(id) %>%
  dplyr::summarize(size.mean = mean(Size),size.median=median(Size), 
                   size.max=max(Size), size.min=min(Size))



bus_summary_time<-bus.all %>% 
                  group_by(id) %>% 
                  summarise(start_time = min(PosixTime), end_time = max(PosixTime), length=length(id))

bus_merge<-merge(bus_summary_time, bus_summary, by="id") %>% merge(bus_summary_size, by="id")

bus_merge$type<-"bus"

#################################################################################
###########################load car data######################################### 
################################################################################


path <- "C:/Users/zy125/Box Sync/PhD/Beijing/car"
file.names <- dir(path, pattern =".txt")
setwd(path)
car <- list()

for (f in file.names) {
  print(f)
  MDdat2 <- importMD_dm(f)
  MDdat2filt <- filterMDv4(MDdat2[[1]])
  car[[f]] <- MDdat2filt[[1]][, c(1,4,5,6,7,8,16)]
}


for (n in 1:16){
  car[[n]]['id'] = n
}
car.all = do.call(rbind, car)

car_summary<-car.all %>%
  group_by(id) %>%
  dplyr::summarize(number.amean = mean(Number), number.gmean=gm_mean(Number),number.median=median(Number), 
                   number.max=max(Number), number.min=min(Number), number.sd=sd(Number),number.q5=quantile(Number, probs = 0.05),
                   number.q25=quantile(Number, probs = 0.25), number.q75=quantile(Number, probs = 0.75),
                   number.q95=quantile(Number, probs = 0.95))

car_summary_size<-car.all %>%
                  group_by(id) %>%
                  dplyr::summarize(size.mean = mean(Size),size.median=median(Size), 
                  size.max=max(Size), size.min=min(Size))


car_summary_time<-car.all %>% 
                  group_by(id) %>% 
                  summarise(start_time = min(PosixTime), end_time = max(PosixTime), length=length(id))

car_merge<-merge(car_summary_time, car_summary, by="id") %>% merge(car_summary_size, by="id")

car_merge$type<-"car"

######################################################################################################
#################################load bicycle data#####################################################
########################################################################################################


path <- "C:/Users/zy125/Box Sync/PhD/Beijing/bicycle"
file.names <- dir(path, pattern =".txt")
setwd(path)
bicycle <- list()

for (f in file.names) {
  print(f)
  MDdat2 <- importMD_dm(f)
  MDdat2filt <- filterMDv4(MDdat2[[1]])
  bicycle[[f]] <- MDdat2filt[[1]][, c(1,4,5,6,7,8,16)]
}


for (n in 1:10){
  bicycle[[n]]['id'] = n
}
bicycle.all = do.call(rbind, bicycle)

bicycle_summary<-bicycle.all %>%
  group_by(id) %>%
  dplyr::summarize(number.amean = mean(Number), number.gmean=gm_mean(Number),number.median=median(Number), 
                   number.max=max(Number), number.min=min(Number), number.sd=sd(Number),number.q5=quantile(Number, probs = 0.05),
                   number.q25=quantile(Number, probs = 0.25), number.q75=quantile(Number, probs = 0.75),
                   number.q95=quantile(Number, probs = 0.95))


bicycle_summary_size<-bicycle.all %>%
                      group_by(id) %>%
                      dplyr::summarize(size.mean = mean(Size),size.median=median(Size), 
                      size.max=max(Size), size.min=min(Size))

bicycle_summary_time<-bicycle.all %>% 
                      group_by(id) %>% 
                      summarise(start_time = min(PosixTime), end_time = max(PosixTime), length=length(id))

bicycle_merge<-merge(bicycle_summary_time, bicycle_summary, by="id") %>% merge(bicycle_summary_size, by="id")

bicycle_merge$type<-"bicycle"

###################################################################################################
################load subway data###################################################################
###################################################################################################


path <- "C:/Users/zy125/Box Sync/PhD/Beijing/output/subway"
file.names <- dir(path, pattern =".txt")
setwd(path)
subway<- list()

for (f in file.names) {
  print(f)
  MDdat2 <- importMD_dm(f)
  MDdat2filt <- filterMDv4(MDdat2[[1]])
  subway[[f]] <- MDdat2filt[[1]][, c(1,4,5,6,7,8,16)]
}


for (n in 1:21){
  subway[[n]]['id'] = n
}
subway.all = do.call(rbind, subway)





##############################################################################################################################
#subway data id==7&12 have small amount of data available. Thus, deleted.id=14 has strange data inclusion. deleted as well.###
###############################################################################################################################
subway.all1<-subway.all[-which(subway.all$id==7 | subway.all$id==12),]

subway_summary<-subway.all1 %>%
  group_by(id) %>%
  dplyr::summarize(number.amean = mean(Number), number.gmean=gm_mean(Number),number.median=median(Number), 
                   number.max=max(Number), number.min=min(Number), number.sd=sd(Number),number.q5=quantile(Number, probs = 0.05),
                   number.q25=quantile(Number, probs = 0.25), number.q75=quantile(Number, probs = 0.75),
                   number.q95=quantile(Number, probs = 0.95))

subway_summary_size<-subway.all1 %>%
                     group_by(id) %>%
                     dplyr::summarize(size.mean = mean(Size),size.median=median(Size), 
                     size.max=max(Size), size.min=min(Size))

subway_summary_time<- subway.all1 %>% 
                      group_by(id) %>% 
                      summarise(start_time = min(PosixTime), end_time = max(PosixTime), length=length(id))


subway_merge<-merge(subway_summary_time, subway_summary, by="id") %>% merge(subway_summary_size, by="id")

subway_merge$type<-"subway"

############################################################################
###############load info data################################################
##########################################################################

commuter<-bind_rows(bus_merge, car_merge) %>% bind_rows(bicycle_merge) %>% bind_rows(subway_merge)


commuter<-commuter %>% mutate(peak = case_when(
  is.weekend(as.Date(start_time, tz="asia/shanghai")) ~ FALSE,
  (hour(start_time)>=7 & hour(start_time)<9) ~ TRUE,
  (hour(start_time)>=17 & hour(start_time)<20) ~ TRUE,
  TRUE ~ FALSE)
)

commuter$weekdays<-as.factor(ifelse(weekdays(commuter$start_time) %in% c("Saturday", "Sunday"), "weekend", "weekday"))

beijing<- read_csv("commuter.csv") ##manually add information


##############################################################################
##################combine travel modes#######################################
##################################################################################



beijing %>% group_by(type)%>% 
  dplyr::summarize(n=n())


beijing %>% group_by(period)%>% 
  dplyr::summarize(n=n())



beijing %>% group_by(pattern, type)%>% 
  dplyr::summarize(n=n())

beijing_summary_size<-beijing %>%
   group_by(type) %>%
   dplyr::summarize( mean = mean(size.mean),sd=sd(size.mean), max=max(size.mean), min=min(size.mean))



beijing_cor<-beijing[,c(13,14,15,16,18,20,27)]

bus_cor<-bus_info[,c(13,14,15,17,20,26)]

bicycle_cor<-bicycle[,c(11,12,13,14,15,17,20,26)]

res2 <- rcorr(as.matrix(beijing_cor))


outdoor<-rbind(bus, bicycle)

outdoor_cor<-outdoor[,c(11,12,13,14,15,17,20,26)]

res2_outdoor <- rcorr(as.matrix(outdoor_cor))


chart.Correlation(beijing_cor, histogram=TRUE)



bus_car<-beijing[which(beijing$type=="bus"| beijing$type=="car" ),]
t.test(g.mean~type, data=bus_car)


bus_bicycle<-beijing[which(beijing$type=="bus"| beijing$type=="bicycle" ),]
t.test(g.mean~type, data=bus_bicycle)

bus_subway<-beijing[which(beijing$type=="bus"| beijing$type=="subway" ),]

t.test(g.mean~type, data=bus_subway)

car_bicycle<-beijing[which(beijing$type=="bicycle"| beijing$type=="car" ),]

t.test(g.mean~type, data=car_bicycle)

car_subway<-beijing[which(beijing$type=="subway"| beijing$type=="car" ),]

t.test(g.mean~type, data=car_subway)

bicycle_subway<-beijing[which(beijing$type=="bicycle"| beijing$type=="subway" ),]

t.test(g.mean~type, data=bicycle_subway)

######plot####################################################################################

tiff("Figure1.tiff", units="in", width=8, height=6, res=300)

ggplot(beijing, aes( x=type, y=g.mean))+geom_boxplot(fill="grey")+theme_bw()+
  scale_x_discrete(labels= c("Bicycle (n=10)", "Bus (n=18)","Car (n=16)","Subway (n=19)"))+labs(x=NULL, y=expression(paste("UFP Number Concentration (Particle Number /  ", cm^3, ")")), fill="Averages Type")+
  theme(text = element_text(size=14),legend.position="bottom")
dev.off()



ggplot(beijing, aes( x=type, y=size.mean))+geom_boxplot(fill="grey")+theme_bw()+
  scale_x_discrete(labels= c("Bicycle (n=10)", "Bus (n=18)","Car (n=16)","Subway (n=19)"))+labs(x=NULL, y=expression(paste("UFP Number Concentration (Particle Number /  ", cm^3, ")")), fill="Averages Type")+
  theme(text = element_text(size=14),legend.position="bottom")


ggplot(beijing, aes( x=type, y=PM25_30))+geom_boxplot(fill="grey")+theme_bw()+
  scale_x_discrete(labels= c("Bicycle (n=10)", "Bus (n=18)","Car (n=16)","Subway (n=18)"))+labs(x=NULL, y=expression(paste("UFP Number Concentration (Particle Number /  ", cm^3, ")")), fill="Averages Type")+
  theme(text = element_text(size=14),legend.position="bottom")



ggplot(beijing, aes( x=type))+geom_boxplot(aes(y=PM25_30*1000),fill="grey")+geom_boxplot(aes(y=g.mean),fill="grey")+theme_bw()+
  scale_x_discrete(labels= c("Bicycle (n=10)", "Bus (n=18)","Car (n=16)","Subway (n=18)"))+labs(x=NULL, y=expression(paste("UFP Number Concentration (Particle Number /  ", cm^3, ")")), fill="Averages Type")+
  theme(text = element_text(size=14),legend.position="bottom")+scale_y_continuous(
    
    # Features of the first axis
    name = "Temperature (Celsius °)",
    
    # Add a second axis and specify its features
    sec.axis = sec_axis(~./100, name="Price ($)")
  )



beijing.peak<-beijing[which(beijing$pattern=="peak"),]

ggplot(beijing.peak, aes(type, a.mean))+geom_boxplot()+theme_bw()


beiing_count_peak<-beijing.peak %>%
  group_by(type) %>% 
  dplyr::summarize(n= n(), mean=mean(a.mean) )

ggplot(subway_vild, aes(a.mean,pattern))+geom_boxplot()


ggplot(beijing.peak, aes(direction, a.mean))+geom_boxplot()


beijing.peak %>%
  group_by(type) %>% 
  dplyr::summarize(n= n(), mean=mean(a.mean),con.sd=sd(a.mean),
                   min=min(a.mean), median=median(a.mean),
                   max=max(a.mean), 
                   len=mean(length), len.sd=sd(length))

beijing.peak %>%
  group_by(type) %>% 
  dplyr::summarize(n= n(), mean=mean(g.mean),con.sd=sd(g.mean),
                   min=min(g.mean), median=median(g.mean),
                   max=max(g.mean), 
                   len=mean(length), len.sd=sd(length))

beijing %>%
  group_by(type, period, pattern) %>% 
  dplyr::summarize(n= n(), total=sum(length), 
                   len=mean(length), len.sd=sd(length))

beijing%>%
  group_by(type) %>% 
  dplyr::summarize(n= n(), total=sum(length), 
                   len=mean(length), len.sd=sd(length))
beijing %>%
  group_by(type, date) %>% 
  dplyr::summarize(n=count(date))

beijing %>%
  group_by(period, type) %>% 
  dplyr::summarize(n= n())

subway.min.num<-subway.all1 %>%
  group_by(time = cut(PosixTime, breaks="1 min")) %>%
  dplyr::summarize(con = gm_mean(Number), id2=unique(id))


subway_micro_summary <- read_csv("C:/Users/zy125/Box Sync/PhD/Beijing/subway_micr.csv")
subway_micro_summary$environment<-factor(subway_micro_summary$environment, levels= c("zgc", "underground", "transfer","overground", "sd") )


subway_micro_summary$type<-NA

subway_micro_summary$type[subway_micro_summary$environment== "zgc"|subway_micro_summary$environment== "sd"]<-"Entrance"

subway_micro_summary$type[subway_micro_summary$environment== "underground"|subway_micro_summary$environment== "overground"]<-"Carriage"

subway_micro_summary$type[subway_micro_summary$environment== "transfer"]<-"Tunnel"


subway_micro_sum<-subway_micro_summary %>%
  group_by(environment, id2) %>%
  dplyr::summarize(n=n(), pmean = mean(con), pmax=max(con), pmin=min(con), pmedian=median(con))


subway_micro_sum<-subway_micro_summary %>%
  group_by(environment) %>%
  dplyr::summarize(n=n(), pmean = mean(con), psd=sd(con),pmax=max(con), pmin=min(con), pmedian=median(con))




zgc_sd<-subway_micro_summary[which(subway_micro_summary$environment=="zgc"| subway_micro_summary$environment=="sd" ),]
t.test(con~environment, data=zgc_sd)

zgc_over<-subway_micro_summary[which(subway_micro_summary$environment=="zgc"| subway_micro_summary$environment=="overground" ),]
t.test(con~environment, data=zgc_over)

zgc_under<-subway_micro_summary[which(subway_micro_summary$environment=="zgc"| subway_micro_summary$environment=="underground" ),]
t.test(con~environment, data=zgc_under)

zgc_transfer<-subway_micro_summary[which(subway_micro_summary$environment=="zgc"| subway_micro_summary$environment=="transfer" ),]
t.test(con~environment, data=zgc_transfer)

sd_over<-subway_micro_summary[which(subway_micro_summary$environment=="sd"| subway_micro_summary$environment=="overground" ),]
t.test(con~environment, data=sd_over)

sd_under<-subway_micro_summary[which(subway_micro_summary$environment=="sd"| subway_micro_summary$environment=="underground" ),]
t.test(con~environment, data=sd_under)

sd_transfer<-subway_micro_summary[which(subway_micro_summary$environment=="sd"| subway_micro_summary$environment=="transfer" ),]
t.test(con~environment, data=sd_transfer)

over_under<-subway_micro_summary[which(subway_micro_summary$environment=="overground"| subway_micro_summary$environment=="underground" ),]
t.test(con~environment, data=over_under)

over_transfer<-subway_micro_summary[which(subway_micro_summary$environment=="overground"| subway_micro_summary$environment=="transfer" ),]
t.test(con~environment, data=over_transfer)


under_transfer<-subway_micro_summary[which(subway_micro_summary$environment=="underground"| subway_micro_summary$environment=="transfer" ),]
t.test(con~environment, data=under_transfer)




tiff("Figure3.tiff", units="in", width=8, height=6, res=300)

ggplot(subway_micro_summary, aes(environment,con, fill=type))+geom_boxplot()+theme_bw()+
  labs(x=NULL, y=expression(paste("UFP Number Concentration (Particle Number / ", cm^3, ")")), fill="Micro-Environments Type")+
  scale_x_discrete(labels= c("ZGC (n=19)", "Underground (n=19)","Transfer (n=19)","Overground (n=19)", "SD (n=19)"))+
  theme(text = element_text(size=14),legend.position="bottom")
dev.off()





subway_micro_summary_envir<-subway_micro_summary %>%
  group_by(environment) %>%
  dplyr::summarize(length=n(), mean = mean(con),sd=sd(con), max=max(con),
            min=min(con))





tem_hum_bp_summary<-tem_hum_bp %>%
  summarize(tem.mean=mean(TEMP), tem.min=min(TEMP), tem.max=max(TEMP), tem.sd=sd(TEMP),
            rh.mean=mean(RH), rh.min=min(RH), rh.max=max(RH), rh.sd=sd(RH),
            bp.mean=mean(BP), bp.min=min(BP), bp.max=max(BP), bp.sd=sd(BP))




wd_ws_1$hour<-cut(wd_ws_1$date, breaks="hour")

wd_ws_hour<-wd_ws_1 %>% 
  group_by(hour) %>% 
  dplyr::summarize(wd=mean(WD), ws=mean(WS))



beijing_amean<-beijing[, c(1,3,4,5,6,7,8,9,10)]
beijing_gmean<-beijing[, c(1,3,4,5,6,7,8,9,11)]


colnames(beijing_amean)[9]<-"mean"
beijing_amean$average<-"arithmetic"
colnames(beijing_gmean)[9]<-"mean"
beijing_gmean$average<-"geometric"
beijing_1<-rbind(beijing_amean, beijing_gmean)


pairwiseCI(a.mean~type, data=commuter, method="Param.diff")

ci.ratio<-pairwiseCI(g.mean~type, data=beijing, method="Param.ratio")

write.csv(ci.ratio, "ci.ratio.csv")

micro_ci_ratio<-pairwiseCI(con~environment, data=subway_micro_summary, method="Param.ratio")

write.csv(micro_ci_ratio, "micro_ci_ratio.csv")
pd <- position_dodge(width = 1)


tiff("Figure3.tiff", units="in", width=8, height=6, res=300)

ggplot(CI, aes(y=estimate, x=comparison))+
  geom_pointrange(aes(ymax=upper, ymin=lower), position = pd, size=0.8)+
  geom_hline(yintercept = 1)+theme_classic()+
  labs(x=NULL, y="Ratio")+theme(text = element_text(size=14))
dev.off()


file.name<- dir(getwd(), pattern =".txt")
col_md<-importMD_dm(file="5519I422_output.txt")
col_md_filt<- filterMDv4(col_md[[1]])
col_md_dat <- col_md_filt[[1]][, c(1,4,5)]


col_md_min<-col_md_dat %>%
  group_by(time = cut(PosixTime, breaks="1 mins")) %>%
  summarize(con = mean(Number))

col_md_min$five<- as.POSIXct(ceiling(as.numeric(col_md_min$PosixTime)/(5 * 60))*(5*60), origin='1970-01-01')
col_md_five <- aggregate(con ~ five, col_md_min, mean)


PNSD <- read_csv("~/Box Sync/Beijing/bj_data/PNSD.csv", col_types = cols(`NaN` = col_datetime(format = "%d/%m/%Y %H:%M:%S")))

ref_dat<-PNSD[,c(1,2)]
colnames(ref_dat)<-c("time","number")

ref_col<-ref_dat[ref_dat$time>="2015-05-19 09:45:00" & ref_dat$time<= "2015-05-19 12:35:00",]

colnames(ref_dat)<-c("time","number")
col<-merge(ref_col, col_md_five, by="time")


ggplot(col, aes(number.x, number.y))+geom_point()+
  geom_smooth(method=lm, se=FALSE)


#10.7-600 nm
PNSD3 <- read_csv("~/Box Sync/Beijing/bj_data/PNSD3.csv", col_types = cols(`NaN` = col_datetime(format = "%d/%m/%Y %H:%M:%S")))

ref_dat3<-PNSD3[,c(1,46)]
colnames(ref_dat3)<-c("time","number")

ref_col3<-ref_dat3[ref_dat3$time>="2015-05-19 09:45:00" & ref_dat3$time<= "2015-05-19 12:35:00",]


col3<-merge(ref_col3, col_md_five, by="time")


#10.7-300 nm
PNSD4 <- read_csv("C:/Users/zy125/Box Sync/PhD/Beijing/bj_data/10_300.csv", col_types = cols(`NaN` = col_datetime(format = "%d/%m/%Y %H:%M:%S")))

ref_dat4<-PNSD4[,c(1,44)]
colnames(ref_dat4)<-c("time","number")

ref_col4<-ref_dat4[ref_dat4$time>="2015-05-19 09:45:00" & ref_dat4$time<= "2015-05-19 12:35:00",]


col4<-merge(ref_col4, col_md_five, by="time")


#08/31/2020 updated to include the 10nm-300nm data as reference data
ref_dat_hr<-ref_dat4 %>%
  group_by(time = cut(time, breaks="1 hour")) %>%
  dplyr::summarize(con = mean(number), gmean=gm_mean(number))

#journey data include all commuter types

journey_bicycle<-journey[which(journey$type=="bicycle"),]

cor(journey_bicycle$pek.a.mean, journey_bicycle$g.mean, use="complete.obs")


journey1<-journey[-c(61,62),] #remove nas





sum_ratio<-groupwiseMean (ratio ~ type, data=journey1, conf=0.95, digits=3)

ggplot(sum_ratio, aes(x = Mean, y = type))+geom_errorbar(aes(xmin = Trad.lower,
                                                            xmax = Trad.upper),
                                                        width = 0.05, size  = 1) +
  geom_point(shape = 15,size  = 4) + theme_bw() + theme(axis.title   = element_text(face  = "bold")) +
  
  ylab("Commuter Types")+geom_vline(xintercept=1)+theme_classic()



#########################temporal adjusted######################################
########################################################################################
ggplot(trip_temp, aes(type, a.mean.adj))+geom_boxplot()+theme_bw()+
  scale_fill_manual(values=c("#336699", "#cc9999"))+ 
  scale_x_discrete(labels= c("Bicycle (n=10)", "Bus (n=18)","Car (n=16)","Subway (n=19)"))+labs(x=NULL, y=expression(paste("Temporal Adjusted Arithmetic Averages of UFP Number Concentration (n/  ",  cm^3,")")))+
  theme( text = element_text(size=16), axis.title.y = element_text(size=14, vjust=0.5, margin = margin(0, 0, 0, 10)))+ theme_bw()+theme(plot.margin=unit(c(1,2,1,2),"cm"))


############################################################################################
######################################################################################


tiff("Figure2.tiff", units="in", width=8, height=6, res=300)
ggplot(temporal_adjusted, aes(type, number, fill=method))+geom_boxplot()+theme_bw()+
  scale_fill_manual(labels = c("Temporal Adjusted Measurements", "Real Time Measurements"), values=c("#336699", "#cc9999"))+ 
  scale_x_discrete(labels= c("Bicycle (n=10)", "Bus (n=18)","Car (n=16)","Subway (n=19)"))+labs(x=NULL, y=expression(paste("UFP Number Concentration (n/ ",  cm^3,")")))+
  theme( text = element_text(size=16), axis.title.y = element_text(size=14, vjust=0.5, margin = margin(0, 0, 0, 2)))+ theme_bw()+theme(plot.margin=unit(c(1,2,1,2),"cm"))+theme(legend.title = element_blank(),legend.position = "bottom")

dev.off()
