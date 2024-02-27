
rm(list=ls())
## ================ Load packages ==============================================
library("tidyverse")
library("data.table")
library("mgsub")
library("kableExtra")
library("nnet") #multinomial logistic regression
library("splines")
library("boot")
library("car")
library("effects")
library("meta")
library("gridExtra")
library("cowplot")
library("ggrepel")




# High TB burden country lists
country_HBC<-c('Angola','Bangladesh','Brazil','Central African Republic','China','Congo','Democratic People\'s Republic of Korea','Democratic Republic of the Congo','Ethiopia','Gabon','India','Indonesia','Kenya','Lesotho','Liberia','Mongolia','Mozambique','Myanmar','Namibia','Nigeria','Pakistan','Papua New Guinea','Philippines','Sierra Leone','South Africa','Thailand','Uganda','United Republic of Tanzania','Viet Nam','Zambia')
country_MDR<-c('Angola','Azerbaijan','Bangladesh','Belarus','China','Democratic People\'s Republic of Korea','Democratic Republic of the Congo','India','Indonesia','Kazakhstan','Kyrgyzstan','Mongolia','Mozambique','Myanmar','Nepal','Nigeria','Pakistan','Papua New Guinea','Peru','Philippines','Republic of Moldova','Russian Federation','Somalia','South Africa','Tajikistan','Ukraine','Uzbekistan','Viet Nam','Zambia','Zimbabwe')
country_HIV<-c('Botswana','Brazil','Cameroon','Central African Republic','China','Congo','Democratic Republic of the Congo','Eswatini','Ethiopia','Gabon','Guinea','Guinea-Bissau','India','Indonesia','Kenya','Lesotho','Liberia','Malawi','Mozambique','Myanmar','Namibia','Nigeria','Philippines','Russian Federation','South Africa','Thailand','Uganda','United Republic of Tanzania','Zambia','Zimbabwe')
country_list<-unique(c(country_HBC,country_MDR,country_HIV))

#load data
WHO_not<-fread("TB_outcomes_2023-12-04.csv")

# Use only HBC countries
WHO_not$g_whoregion<-mgsub(WHO_not$g_whoregion,c("AFR","AMR","EMR","EUR","SEA","WPR"),c("African Region","Region of the Americas","Eastern Mediterranean Region","European Region","South-East Asia Region","Western Pacific Region"))
WHO_not<-WHO_not[country%in%country_list]

WHO_not<-WHO_not[year%in%c(2012:2021)]

# Simplify data set
dat <- dplyr::select(WHO_not,country,iso3,g_whoregion,year,newrel_succ,newrel_fail,newrel_died,newrel_lost)



#calculate observed proportions
dat <- dat %>%
  mutate(
    success=newrel_succ/(newrel_succ+newrel_fail+newrel_died+newrel_lost),
    fail=newrel_fail/(newrel_succ+newrel_fail+newrel_died+newrel_lost),
    died=newrel_died/(newrel_succ+newrel_fail+newrel_died+newrel_lost),
    lost=newrel_lost/(newrel_succ+newrel_fail+newrel_died+newrel_lost)
  )






#prediction data frame for all countries
countries<-unique(dat$iso3)
gr<-expand.grid(2012:2021,countries)
pred_multi_sp_boot<-data.frame(iso3=gr[,2],year=gr[,1],success=NA,fail=NA,died=NA,lost=NA, L.success=NA,L.fail=NA,
                               L.died=NA, L.lost=NA, U.success=NA, U.fail=NA,  U.died=NA, U.lost=NA,
                               p.success=NA,p.fail=NA,p.died=NA,p.lost=NA,
                               ci.success=NA,ci.fail=NA,ci.died=NA,ci.lost=NA)


# Create an empty list to store the results
results_list_sp_boot <- list()


# Iterate over each country code
for (c in countries) {
  # Subset the data for the current country
  country_data <- dat%>% dplyr::filter(iso3 == c & !is.na(newrel_succ) & !is.na(newrel_fail) & !is.na(newrel_died) & !is.na(newrel_lost))
  country_data_before2020 <- country_data %>% dplyr::filter(year<2020)
  
  
  # Create the prediction data frame for the country
  gr_country <- expand.grid(2012:2021, c)
  pred_country <- data.frame(iso3 = gr_country[, 2], year = gr_country[, 1], 
                             success = NA, fail = NA, died = NA, lost = NA)
  
  # Define the function to return bootstrapped coefficients
  multiCoef <- function(data, indices, formula) {
    #while(length(unique(indices))<6){
    while(class(try(ns(data$year[indices], df=3)))[1]=="try-error"){
      indices<-sample(x=1:nrow(data),size=nrow(data),replace=TRUE)
    }
    d <- data[indices, ]
    fit <- multinom(formula, data = d)
    pred <- predict(fit, newdata = pred_country, se.fit = TRUE, type = "probs")
    return(pred)
  }
  
  # Set the seed for reproducibility
  set.seed(373)
  
  # Perform bootstrapping for the current country
  ## multinomial logistic regression model
  ##natural splines with 3 degrees of freedom
  coef.boot2 <- boot(data = country_data_before2020, statistic = multiCoef, R = 1e4, 
                     formula = cbind(newrel_succ, newrel_fail, newrel_died, newrel_lost) ~ ns(year, df=3), cl = NULL)
  
  #calculating pvalues and confidence intervals for the ratios between observed and predicted probailities for 2020 and 2021- column 9 (success), 19 (failure), 29 (death), 39 (LTFU) correspond to 2020 values and column 10 (success),20 (failure),30 (death),40 (LTFU) correspond to 2021 values
  idx<-9
  ff<-function(x){abs(dat$success[dat$iso3==c & dat$year==2020]/quantile(coef.boot2$t[,idx],probs=x, na.rm = TRUE)-1)}
  pTmp<-optimize(f=ff,interval=c(0,1))$minimum
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c & pred_multi_sp_boot$year==2020,"p.success"]<-2*ifelse(pTmp>0.5,1-pTmp,pTmp)
  ciPrctl<-paste(collapse="_",quantile(dat$success[dat$iso3==c & dat$year==2020]/coef.boot2$t[,idx],probs=c(0.025,0.975),na.rm = TRUE))
  ciHdi<-paste(collapse="_",HDInterval::hdi(dat$success[dat$iso3==c & dat$year==2020]/coef.boot2$t[,idx]))
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c & pred_multi_sp_boot$year==2020,"ci.success"]<-ciPrctl
 
  idx<-10
  ff<-function(x){abs(dat$success[dat$iso3==c & dat$year==2021]/quantile(coef.boot2$t[,idx],probs=x, na.rm = TRUE)-1)}
  pTmp<-optimize(f=ff,interval=c(0,1))$minimum
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c & pred_multi_sp_boot$year==2021,"p.success"]<-2*ifelse(pTmp>0.5,1-pTmp,pTmp)
  ciPrctl<-paste(collapse="_",quantile(dat$success[dat$iso3==c & dat$year==2021]/coef.boot2$t[,idx],probs=c(0.025,0.975),na.rm = TRUE))
  ciHdi<-paste(collapse="_",HDInterval::hdi(dat$success[dat$iso3==c & dat$year==2021]/coef.boot2$t[,idx]))
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c & pred_multi_sp_boot$year==2021,"ci.success"]<-ciPrctl
  
  idx<-19
  ff<-function(x){abs(dat$fail[dat$iso3==c & dat$year==2020]/quantile(coef.boot2$t[,idx],probs=x, na.rm = TRUE)-1)}
  pTmp<-optimize(f=ff,interval=c(0,1))$minimum
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c & pred_multi_sp_boot$year==2020,"p.fail"]<-2*ifelse(pTmp>0.5,1-pTmp,pTmp)
  ciPrctl<-paste(collapse="_",quantile(dat$fail[dat$iso3==c & dat$year==2020]/coef.boot2$t[,idx],probs=c(0.025,0.975),na.rm = TRUE))
  ciHdi<-paste(collapse="_",HDInterval::hdi(dat$fail[dat$iso3==c & dat$year==2020]/coef.boot2$t[,idx]))
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c & pred_multi_sp_boot$year==2020,"ci.fail"]<-ciPrctl
  
  idx<-20
  ff<-function(x){abs(dat$fail[dat$iso3==c & dat$year==2021]/quantile(coef.boot2$t[,idx],probs=x, na.rm = TRUE)-1)}
  pTmp<-optimize(f=ff,interval=c(0,1))$minimum
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c & pred_multi_sp_boot$year==2021,"p.fail"]<-2*ifelse(pTmp>0.5,1-pTmp,pTmp)
  ciPrctl<-paste(collapse="_",quantile(dat$fail[dat$iso3==c & dat$year==2021]/coef.boot2$t[,idx],probs=c(0.025,0.975),na.rm = TRUE))
  ciHdi<-paste(collapse="_",HDInterval::hdi(dat$fail[dat$iso3==c & dat$year==2021]/coef.boot2$t[,idx]))
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c & pred_multi_sp_boot$year==2021,"ci.fail"]<-ciPrctl
  
  
  
  idx<-29
  ff<-function(x){abs(dat$died[dat$iso3==c & dat$year==2020]/quantile(coef.boot2$t[,idx],probs=x, na.rm = TRUE)-1)}
  pTmp<-optimize(f=ff,interval=c(0,1))$minimum
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c & pred_multi_sp_boot$year==2020,"p.died"]<-2*ifelse(pTmp>0.5,1-pTmp,pTmp)
  ciPrctl<-paste(collapse="_",quantile(dat$died[dat$iso3==c & dat$year==2020]/coef.boot2$t[,idx],probs=c(0.025,0.975),na.rm = TRUE))
  ciHdi<-paste(collapse="_",HDInterval::hdi(dat$died[dat$iso3==c & dat$year==2020]/coef.boot2$t[,idx]))
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c & pred_multi_sp_boot$year==2020,"ci.died"]<-ciPrctl
 
  idx<-30
  ff<-function(x){abs(dat$died[dat$iso3==c & dat$year==2021]/quantile(coef.boot2$t[,idx],probs=x, na.rm = TRUE)-1)}
  pTmp<-optimize(f=ff,interval=c(0,1))$minimum
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c & pred_multi_sp_boot$year==2021,"p.died"]<-2*ifelse(pTmp>0.5,1-pTmp,pTmp)
  ciPrctl<-paste(collapse="_",quantile(dat$died[dat$iso3==c & dat$year==2021]/coef.boot2$t[,idx],probs=c(0.025,0.975),na.rm = TRUE))
  ciHdi<-paste(collapse="_",HDInterval::hdi(dat$died[dat$iso3==c & dat$year==2021]/coef.boot2$t[,idx]))
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c & pred_multi_sp_boot$year==2021,"ci.died"]<-ciPrctl
  
   
  idx<-39
  ff<-function(x){abs(dat$lost[dat$iso3==c & dat$year==2020]/quantile(coef.boot2$t[,idx],probs=x, na.rm = TRUE)-1)}
  pTmp<-optimize(f=ff,interval=c(0,1))$minimum
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c & pred_multi_sp_boot$year==2020,"p.lost"]<-2*ifelse(pTmp>0.5,1-pTmp,pTmp)
  ciPrctl<-paste(collapse="_",quantile(dat$lost[dat$iso3==c & dat$year==2020]/coef.boot2$t[,idx],probs=c(0.025,0.975),na.rm = TRUE))
  ciHdi<-paste(collapse="_",HDInterval::hdi(dat$lost[dat$iso3==c & dat$year==2020]/coef.boot2$t[,idx]))
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c & pred_multi_sp_boot$year==2020,"ci.lost"]<-ciPrctl
  
  idx<-40
  ff<-function(x){abs(dat$lost[dat$iso3==c & dat$year==2021]/quantile(coef.boot2$t[,idx],probs=x, na.rm = TRUE)-1)}
  pTmp<-optimize(f=ff,interval=c(0,1))$minimum
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c & pred_multi_sp_boot$year==2021,"p.lost"]<-2*ifelse(pTmp>0.5,1-pTmp,pTmp)
  ciPrctl<-paste(collapse="_",quantile(dat$lost[dat$iso3==c & dat$year==2021]/coef.boot2$t[,idx],probs=c(0.025,0.975),na.rm = TRUE ))
  ciHdi<-paste(collapse="_",HDInterval::hdi(dat$lost[dat$iso3==c & dat$year==2021]/coef.boot2$t[,idx]))
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c & pred_multi_sp_boot$year==2021,"ci.lost"]<-ciPrctl
  
  # Store the results in the list
  results_list_sp_boot[[c]] <- coef.boot2
  
  #extracting confidence intervals for the predicted probabilities
  CI <- Confint(results_list_sp_boot[[c]], type = "perc")

  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c,3:6] <- results_list_sp_boot[[c]]$t0[,1:4]
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c,7] <- CI[1:10]
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c,8] <- CI[11:20]
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c,9] <- CI[21:30]
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c,10] <- CI[31:40]
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c,11] <- CI[1:10,2]
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c,12] <- CI[11:20,2]
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c,13] <- CI[21:30,2]
  pred_multi_sp_boot[pred_multi_sp_boot$iso3==c,14] <- CI[31:40,2]
  
}


#show data
pred_multi_sp_boot %>%
  dplyr::mutate(
    totalProb=success+fail+died+lost
  ) %>%
  knitr::kable(col.names=c("Country code","Year","Succ_prob (Expected)","Fail_prob (expected)","Died_prob (expected)","Lost_prob (expected)", "L.succ_prob","L.fail_prob",
                           "L.died_prob","L.lost_prob","U.succ_prob","U.fail_prob",
                           "U.died_prob","U.lost_prob","Sum of probabilities")) %>%
  kableExtra::kable_styling(full_width = FALSE)


#plot data

datObsLong_sp_boot<-dat %>%
  dplyr::select(iso3,year,g_whoregion,country,success,fail,died,lost) %>%
  tidyr::pivot_longer(cols=c(success,fail,died,lost),names_to="group",values_to="observed")
datExpLong_sp_boot<-pred_multi_sp_boot %>%
  dplyr::select(iso3,year,success,fail,died,lost) %>%
  tidyr::pivot_longer(cols=c(success,fail,died,lost),names_to="group",values_to="expected")
datExpLLong_sp_boot<-pred_multi_sp_boot %>%
  dplyr::select(iso3,year,L.success,L.fail,L.died, L.lost) %>%
  tidyr::pivot_longer(cols=c(L.success,L.fail,L.died, L.lost),names_to="group",values_to="L")
datExpLLong_sp_boot$group<-gsub(pattern="L.",replacement="",datExpLLong_sp_boot$group)
datExpULong_sp_boot<-pred_multi_sp_boot %>%
  dplyr::select(iso3,year,U.success,U.fail,U.died, U.lost) %>%
  tidyr::pivot_longer(cols=c(U.success,U.fail,U.died, U.lost),names_to="group",values_to="U")
datExpULong_sp_boot$group<-gsub(pattern="U.",replacement="",datExpULong_sp_boot$group)

datRRpvalue_sp_boot<-pred_multi_sp_boot %>%
  dplyr::select(iso3,year,p.success,p.fail,p.died, p.lost) %>%
  tidyr::pivot_longer(cols=c(p.success,p.fail,p.died, p.lost),names_to="group",values_to="p")
datRRpvalue_sp_boot$group<-gsub(pattern="p.",replacement="",datRRpvalue_sp_boot$group)

datRRci_sp_boot<-pred_multi_sp_boot %>%
  dplyr::select(iso3,year,ci.success,ci.fail,ci.died, ci.lost) %>%
  tidyr::pivot_longer(cols=c(ci.success,ci.fail,ci.died, ci.lost),names_to="group",values_to="ci")
datRRci_sp_boot$group<-gsub(pattern="ci.",replacement="",datRRci_sp_boot$group)

datFullLong_sp_boot<-datObsLong_sp_boot %>%
  mutate(
    expected=datExpLong_sp_boot$expected[match(paste(sep="_",iso3,year,group),paste(sep="_",datExpLong_sp_boot$iso3,datExpLong_sp_boot$year,datExpLong_sp_boot$group))],
    expectedLow=datExpLLong_sp_boot$L[match(paste(sep="_",iso3,year,group),paste(sep="_",datExpLLong_sp_boot$iso3,datExpLLong_sp_boot$year,datExpLLong_sp_boot$group))],
    expectedUpp=datExpULong_sp_boot$U[match(paste(sep="_",iso3,year,group),paste(sep="_",datExpULong_sp_boot$iso3,datExpULong_sp_boot$year,datExpULong_sp_boot$group))],
    RRci=datRRci_sp_boot$ci[match(paste(sep = "_",iso3,year,group),paste(sep = "_",datRRci_sp_boot$iso3,datRRci_sp_boot$year,datRRci_sp_boot$group))],
    RRpvalue=datRRpvalue_sp_boot$p[match(paste(sep = "_",iso3,year,group),paste(sep = "_",datRRpvalue_sp_boot$iso3,datRRpvalue_sp_boot$year,datRRpvalue_sp_boot$group))]
    )

datFullLong_sp_boot %>%
  ggplot(mapping=aes(x=year)) +
  geom_ribbon(mapping=aes(ymin=expectedLow,ymax=expectedUpp),fill="orange",alpha=0.5) +
  geom_line(mapping=aes(y=expected),col="orange",lwd=1.25) +
  geom_point(mapping=aes(y=observed),col="steelblue",size=2) + 
  facet_wrap(~paste(sep=" ",iso3,group),ncol = 4,scales = "free_y")

ggsave("Plot_ns_splines_boot_df3_20231204.png",device="png",width=20,height=92,units=c("cm"))


###########################################################################################################################



#calculating Risk ratios and the confidence intervals
refDat2<-datFullLong_sp_boot
refDat2<-refDat2 %>%
  dplyr::mutate(
    RR=
      round(observed/expected,2),
    RR_Low=
      round(observed/expectedUpp,2),
    RR_Upp=
      round(observed/expectedLow,2)
  ) %>%
  mutate(
    RR_95CI=paste(sep="",format(big.mark=",",RR)," (",format(big.mark=",",RR_Low),", ",format(big.mark=",",RR_Upp),")")
    
  )
refDat2 %>%
  dplyr::select(c(iso3,group,observed,expected,RR_95CI)) %>%
  knitr::kable(col.names=c("Country code","Group","Observed","Expected","Risk Ratio (95% CI)"),caption = "Observed and expected proportions for 2020") %>%
  kableExtra::kable_styling(full_width = FALSE)


##################################################################################################
#Meta analysis



#Success

meta_succ2020_2 <- metagen(studlab = country, sm="RR", TE=log(RR), subgroup = g_whoregion, fixed = FALSE,data = refDat2%>%filter(group=="success", year==2020), lower = log(RR_Low), upper = log(RR_Upp))

summary(meta_succ2020_2)

#gave an error: Error in (function (yi, vi, sei, weights, ai, bi, ci, di, n1i, n2i, x1i,  : Fisher scoring algorithm did not converge. See 'help(rma)' for possible remedies.
#added the control part at the end
meta_succ2021_2 <- metagen(studlab = country, sm="RR", TE=log(RR), subgroup = g_whoregion, fixed = FALSE,data = refDat2%>%filter(group=="success", year==2021), lower = log(RR_Low), upper = log(RR_Upp), control = list(stepadj=0.5))

summary(meta_succ2021_2)


png("forest_success2020_3.png", width = 1500, height = 2000, res=120)
forest(x=meta_succ2020_2, print.subgroup.name=FALSE,text.random="Overall summary",text.random.w="Regional summary",ref=1,at=(c(0.1, 1, 10)),
             leftcols=c("studlab"),rightcols=c("effect","ci"),weight.study="same",weight.subgroup="same",test.subgroup.random=FALSE,xlim=c(0.1,10),
             lty.random = 0,sortvar=-TE,leftlabs = c("Country"),smlab=("Observed:Expected (Success)"),rightlabs = c("Risk Ratio","95% CI"),
             col.study="black",col.diamond = "white",col.diamond.lines ="black",col.square="black",col.inside="black",squaresize=0.5,col.by="black",allstudies=FALSE)
dev.off()

png("forest_success2021_3.png", width = 1500, height = 2000, res=120)
forest(x=meta_succ2021_2, print.subgroup.name=FALSE,text.random="Overall summary",text.random.w="Regional summary",ref=1,at=(c(0.1, 1, 10)),
       leftcols=c("studlab"),rightcols=c("effect","ci"),weight.study="same",weight.subgroup="same",test.subgroup.random=FALSE,xlim=c(0.1,10),
       lty.random = 0,sortvar=-TE,leftlabs = c("Country"),smlab=("Observed:Expected (Success)"),rightlabs = c("Risk Ratio","95% CI"),
       col.study="black",col.diamond = "white",col.diamond.lines ="black",col.square="black",col.inside="black",squaresize=0.5,col.by="black",allstudies=FALSE)
dev.off()

#Failure
meta_fail2020_2 <- metagen(studlab = country, sm="RR", TE=log(RR), subgroup = g_whoregion, fixed = FALSE,data = refDat2%>%filter(group=="fail", year==2020), lower = log(RR_Low), upper = log(RR_Upp))

summary(meta_fail2020_2)

meta_fail2021_2 <- metagen(studlab = country, sm="RR", TE=log(RR), subgroup = g_whoregion, fixed = FALSE,data = refDat2%>%filter(group=="fail", year==2021), lower = log(RR_Low), upper = log(RR_Upp))

summary(meta_fail2021_2)

png("forest_fail2020_2.png", width = 1500, height = 2000, res=120)
forest(x=meta_fail2020_2, print.subgroup.name=FALSE,text.random="Overall summary",text.random.w="Regional summary",ref=1,at=(c(0.1, 1, 10)),
             leftcols=c("studlab"),rightcols=c("effect","ci"),weight.study="same",weight.subgroup="same",test.subgroup.random=FALSE,xlim=c(0.1,10),
             lty.random = 0,sortvar=-TE,leftlabs = c("Country"),smlab=("Observed:Expected (Failure)"),rightlabs = c("Risk Ratio","95% CI"),
       col.study="black",col.diamond = "white",col.diamond.lines ="black",col.square="black",col.inside="black",squaresize=0.5,col.by="black",allstudies=FALSE)
dev.off()


png("forest_fail2021_2.png", width = 1500, height = 2000, res=120)
forest(x=meta_fail2021_2, print.subgroup.name=FALSE,text.random="Overall summary",text.random.w="Regional summary",ref=1,at=(c(0.1, 1, 10)),
       leftcols=c("studlab"),rightcols=c("effect","ci"),weight.study="same",weight.subgroup="same",test.subgroup.random=FALSE,xlim=c(0.1,10),
       lty.random = 0,sortvar=-TE,leftlabs = c("Country"),smlab=("Observed:Expected (Failure)"),rightlabs = c("Risk Ratio","95% CI"),
       col.study="black",col.diamond = "white",col.diamond.lines ="black",col.square="black",col.inside="black",squaresize=0.5,col.by="black",allstudies=FALSE)
dev.off()

#Died
meta_died2020_2 <- metagen(studlab = country, sm="RR", TE=log(RR), subgroup = g_whoregion, fixed = FALSE,data = refDat2%>%filter(group=="died", year==2020), lower = log(RR_Low), upper = log(RR_Upp))

summary(meta_died2020_2)

meta_died2021_2 <- metagen(studlab = country, sm="RR", TE=log(RR), subgroup = g_whoregion, fixed = FALSE,data = refDat2%>%filter(group=="died", year==2021), lower = log(RR_Low), upper = log(RR_Upp))

summary(meta_died2021_2)


png("forest_died2020_2.png", width = 1500, height = 2000, res=120)
forest(x=meta_died2020_2, print.subgroup.name=FALSE,text.random="Overall summary",text.random.w="Regional summary",ref=1,at=(c(0.1, 1, 10)),
       leftcols=c("studlab"),rightcols=c("effect","ci"),weight.study="same",weight.subgroup="same",test.subgroup.random=FALSE,xlim=c(0.1,10),
       lty.random = 0,sortvar=-TE,leftlabs = c("Country"),smlab=("Observed:Expected (Death)"),rightlabs = c("Risk Ratio","95% CI"),
       col.study="black",col.diamond = "white",col.diamond.lines ="black",col.square="black",col.inside="black",squaresize=0.5,col.by="black",allstudies=FALSE)

dev.off()


png("forest_died2021_2.png", width = 1500, height = 2000, res=120)
forest(x=meta_died2021_2, print.subgroup.name=FALSE,text.random="Overall summary",text.random.w="Regional summary",ref=1,at=(c(0.1, 1, 10)),
             leftcols=c("studlab"),rightcols=c("effect","ci"),weight.study="same",weight.subgroup="same",test.subgroup.random=FALSE,xlim=c(0.1,10),
             lty.random = 0,sortvar=-TE,leftlabs = c("Country"),smlab=("Observed:Expected (Death)"),rightlabs = c("Risk Ratio","95% CI"),
       col.study="black",col.diamond = "white",col.diamond.lines ="black",col.square="black",col.inside="black",squaresize=0.5,col.by="black",allstudies=FALSE)
dev.off()

#Lost
meta_lost2020_2 <- metagen(studlab = country, sm="RR", TE=log(RR), subgroup = g_whoregion, fixed = FALSE,data = refDat2%>%filter(group=="lost", year==2020), lower = log(RR_Low), upper = log(RR_Upp))

summary(meta_lost2020_2)

meta_lost2021_2 <- metagen(studlab = country, sm="RR", TE=log(RR), subgroup = g_whoregion, fixed = FALSE,data = refDat2%>%filter(group=="lost", year==2021), lower = log(RR_Low), upper = log(RR_Upp))

summary(meta_lost2021_2)

png("forest_lost2020_2.png", width = 1500, height = 2000, res=120)
forest(x=meta_lost2020_2, print.subgroup.name=FALSE,text.random="Overall summary",text.random.w="Regional summary",ref=1,at=(c(0.1, 1, 10)),
             leftcols=c("studlab"),rightcols=c("effect","ci"),weight.study="same",weight.subgroup="same",test.subgroup.random=FALSE,xlim=c(0.1,10),
             lty.random = 0,sortvar=-TE,leftlabs = c("Country"),smlab=("Observed:Expected (LTFU)"),rightlabs = c("Risk Ratio","95% CI"),
       col.study="black",col.diamond = "white",col.diamond.lines ="black",col.square="black",col.inside="black",squaresize=0.5,col.by="black",allstudies=FALSE)
dev.off()


png("forest_lost2021_2.png", width = 1500, height = 2000, res=120)
forest(x=meta_lost2021_2, print.subgroup.name=FALSE,text.random="Overall summary",text.random.w="Regional summary",ref=1,at=(c(0.1, 1, 10)),
             leftcols=c("studlab"),rightcols=c("effect","ci"),weight.study="same",weight.subgroup="same",test.subgroup.random=FALSE,xlim=c(0.1,10),
             lty.random = 0,sortvar=-TE,leftlabs = c("Country"),smlab=("Observed:Expected (LTFU)"),rightlabs = c("Risk Ratio","95% CI"),
       col.study="black",col.diamond = "white",col.diamond.lines ="black",col.square="black",col.inside="black",squaresize=0.5,col.by="black",allstudies=FALSE)
dev.off()
#############################################################################################3
#Creating year-outcome specific data frames and adding the p-value calculated from the meta analysis

#success

refDat_succ2020_2 <- refDat2 %>% 
  filter(group=="success", year==2020)


#refDat_succ2020_2$pvalue <- meta_succ2020_2$pval

refDat_succ2021_2 <- refDat2 %>% 
  filter(group=="success", year==2021)


#refDat_succ2021_2$pvalue <- meta_succ2021_2$pval


#fail

refDat_fail2020_2 <- refDat2 %>% 
  filter(group=="fail", year==2020)

#refDat_fail2020_2$pvalue <- meta_fail2020_2$pval

refDat_fail2021_2 <- refDat2 %>% 
  filter(group=="fail", year==2021)

#refDat_fail2021_2$pvalue <- meta_fail2021_2$pval

#died

refDat_died2020_2 <- refDat2 %>% 
  filter(group=="died", year==2020)

#refDat_died2020_2$pvalue <- meta_died2020_2$pval


refDat_died2021_2 <- refDat2 %>% 
  filter(group=="died", year==2021)

#refDat_died2021_2$pvalue <- meta_died2021_2$pval



#loss to follow up
refDat_lost2020_2 <- refDat2 %>% 
  filter(group=="lost", year==2020)

#refDat_lost2020_2$pvalue <- meta_lost2020_2$pval

refDat_lost2021_2 <- refDat2 %>% 
  filter(group=="lost", year==2021)

#refDat_lost2021_2$pvalue <- meta_lost2021_2$pval


#######################################################################################################################################################
#Table 1

dat <- dat %>%
  filter (iso3 %in% c("AGO", "AZE", "BGD", "BLR", "BWA", "BRA", "CMR", "CAF", "CHN", "COG", "PRK", "COD" ,"SWZ", "ETH" ,"GAB" ,"GIN" ,"GNB" ,"IND" ,"IDN" ,"KEN" ,"KGZ",
                      "LSO", "LBR", "MWI", "MNG", "MOZ", "MMR", "NAM", "NPL", "NGA", "PAK", "PNG", "PER" ,"PHL", "MDA", "RUS", "SLE", "SOM", "ZAF", "TJK" ,"THA" ,"UGA", "UKR",
                      "TZA", "VNM" ,"ZMB", "ZWE")) 

dat_2012 <- dat %>%
  filter(year==2012)

dat_2013 <- dat %>%
  filter(year==2013)

dat_2014 <- dat %>%
  filter(year==2014)

dat_2015 <- dat %>%
  filter(year==2015)

dat_2016 <- dat %>%
  filter(year==2016)

dat_2017 <- dat %>%
  filter(year==2017)

dat_2018 <- dat %>%
  filter(year==2018)

dat_2019 <- dat %>%
  filter(year==2019)

dat_2020 <- dat %>%
  filter(year==2020)

dat_2021 <- dat %>%
  filter(year==2021)

tab <- data.frame(Success = c(sum(dat_2012$newrel_succ, na.rm = T), sum(dat_2013$newrel_succ, na.rm = T), sum(dat_2014$newrel_succ, na.rm = T), sum(dat_2015$newrel_succ, na.rm = T), sum(dat_2016$newrel_succ, na.rm = T), sum(dat_2017$newrel_succ, na.rm = T), sum(dat_2018$newrel_succ, na.rm = T), sum(dat_2019$newrel_succ, na.rm = T), sum(dat_2020$newrel_succ, na.rm = T), sum(dat_2021$newrel_succ, na.rm = T)),
                  Failure = c(sum(dat_2012$newrel_fail, na.rm = T), sum(dat_2013$newrel_fail, na.rm = T), sum(dat_2014$newrel_fail, na.rm = T), sum(dat_2015$newrel_fail, na.rm = T), sum(dat_2016$newrel_died, na.rm = T), sum(dat_2017$newrel_fail, na.rm = T), sum(dat_2018$newrel_fail, na.rm = T), sum(dat_2019$newrel_fail, na.rm = T), sum(dat_2020$newrel_fail, na.rm = T), sum(dat_2021$newrel_fail, na.rm = T)),
                  Death = c(sum(dat_2012$newrel_died, na.rm = T), sum(dat_2013$newrel_died, na.rm = T), sum(dat_2014$newrel_died, na.rm = T), sum(dat_2015$newrel_died, na.rm = T), sum(dat_2016$newrel_died, na.rm = T), sum(dat_2017$newrel_died, na.rm = T), sum(dat_2018$newrel_died, na.rm = T), sum(dat_2019$newrel_died, na.rm = T), sum(dat_2020$newrel_died, na.rm = T), sum(dat_2021$newrel_died, na.rm = T)),
                  Loss = c(sum(dat_2012$newrel_lost, na.rm = T), sum(dat_2013$newrel_lost, na.rm = T), sum(dat_2014$newrel_lost, na.rm = T), sum(dat_2015$newrel_lost, na.rm = T), sum(dat_2016$newrel_lost, na.rm = T), sum(dat_2017$newrel_lost, na.rm = T), sum(dat_2018$newrel_lost, na.rm = T), sum(dat_2019$newrel_lost, na.rm = T), sum(dat_2020$newrel_lost, na.rm = T), sum(dat_2021$newrel_lost, na.rm = T)))

row.names(tab) <- c("2012", "2013", "2014", "2015", "2016", "2017", "2018", "2019", "2020", "2021") 

tab <- tab %>%
  mutate(
    prob_success=Success/(Success+Failure+Death+Loss),
    prob_fail=Failure/(Success+Failure+Death+Loss),
    prob_died=Death/(Success+Failure+Death+Loss),
    prob_lost=Loss/(Success+Failure+Death+Loss)
  )

tab_2 <- tab %>%
  dplyr::select(c(prob_success,prob_fail,prob_died,prob_lost))
tab_2 <- t(tab_2)
row.names(tab_2) <- c("Success", "Failure", "Death", "Loss to follow up")

tab_2 %>%
  knitr::kable(caption = "Probailities of TB treatment Outcomes by Year (2012-2021) aggregated across 49 high TB, TB/HIV and MDR-TB burden countries", digits = 4) %>%
  kableExtra::kable_styling(full_width = FALSE)


write.table(tab_2, file = "table1.csv", sep = ",", row.names = TRUE)

#####################################################################################################################################


#Figure 1 
datFullLong_sp_boot <- datFullLong_sp_boot %>%
  mutate(
    Outcome=case_when(
      group=="success" ~ "Success",
      group=="fail" ~ "Failure",
      group=="died" ~ "Death",
      group=="lost" ~ "Loss-to-follow-up"
    )
  )


#plot of observed vs expected proportions for 2020


# ggplot(datFullLong_sp_boot %>% filter(year==2020), aes(x=observed, y=expected, col=g_whoregion, label=iso3))+
#   geom_point() +
#   #geom_text(check_overlap = TRUE) +
#   geom_text_repel(aes(label=iso3), size=2) +
#   geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
#   facet_wrap(~paste(sep=" ",Outcome),ncol = 2, scales = "free") +
#   labs(color="WHO Region") +
#   xlab("Observed") +
#   ylab("Expected") +
#   theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10),
#         legend.title = element_text(size = 12), strip.text.x = element_text(size = 12))
# 
# 
# ggplot(refDat2020, aes(x=observed, y=expected, col=g_whoregion, label=iso3))+
#   geom_point() +
#   geom_text(label=ifelse(refDat2020$pvalue<0.05,refDat2020$iso3,""), hjust=0.75, vjust=0) +
#   geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
#   geom_hline(yintercept=0, linetype="solid", colour="black", size=0.1)+
#   geom_vline(xintercept=0, linetype="solid", colour="black", size=0.1)+
#   scale_x_continuous(expand = c(0, 0), limits = c(0, NA)) + 
#   scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
#   facet_wrap(~paste(sep=" ",Outcome),ncol = 2, scales = "free") +
#   labs(color="WHO Region") +
#   xlab("Observed") +
#   ylab("Expected") +
#   theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10),
#         legend.title = element_text(size = 12), strip.text.x = element_text(size = 12))

 f1 <- ggplot(refDat_succ2020_2, aes(x=observed, y=expected, col=g_whoregion, label=iso3))+
  geom_point(size=4) +
  #geom_text(label=ifelse(refDat_succ2020_2$pvalue<0.05,refDat_succ2020_2$iso3,""), hjust=0.75, vjust=0, size=7) +
  geom_text_repel(aes(label=ifelse(refDat_succ2020_2$RRpvalue<0.05,refDat_succ2020_2$iso3,"")), size=6)+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_hline(yintercept=0, linetype="solid", colour="black", size=0.1)+
  geom_vline(xintercept=0, linetype="solid", colour="black", size=0.1)+
  scale_x_continuous(expand = c(0, 0), limits = c(0.6, 1)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0.6, 1)) +
  labs(color="WHO Region") +
  ggtitle("Success")+
  xlab("Observed") +
  ylab("Expected") +
  theme(plot.title = element_text(hjust=0.5, size = 16), axis.title = element_text(size = 14),legend.position = "none")


f2 <-ggplot(refDat_fail2020_2, aes(x=observed, y=expected, col=g_whoregion, label=iso3))+
  geom_point(size=4) +
  #geom_text(label=ifelse(refDat_fail2020_2$pvalue<0.05,refDat_fail2020_2$iso3,""), hjust=0.75, vjust=0, size=7) +
  geom_text_repel(aes(label=ifelse(refDat_fail2020_2$RRpvalue<0.05,refDat_fail2020_2$iso3,"")), size=6, max.overlaps = 30)+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_hline(yintercept=0, linetype="solid", colour="black", size=0.1)+
  geom_vline(xintercept=0, linetype="solid", colour="black", size=0.1)+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.15)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.15)) +
  labs(color="WHO Region") +
  ggtitle("Failure")+
  xlab("Observed") +
  ylab("Expected") +
  theme(plot.title = element_text(hjust=0.5, size = 16), axis.title = element_text(size = 14),legend.position = "none")


f3 <- ggplot(refDat_died2020_2, aes(x=observed, y=expected, col=g_whoregion, label=iso3))+
  geom_point(size=4) +
  #geom_text(label=ifelse(refDat_died2020_2$pvalue<0.05,refDat_died2020_2$iso3,""), hjust=0.75, vjust=0, size=7) +
  geom_text_repel(aes(label=ifelse(refDat_died2020_2$RRpvalue<0.05,refDat_died2020_2$iso3,"")), size=6)+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_hline(yintercept=0, linetype="solid", colour="black", size=0.1)+
  geom_vline(xintercept=0, linetype="solid", colour="black", size=0.1)+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.18)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.18)) +
  labs(color="WHO Region") +
  ggtitle("Death")+
  xlab("Observed") +
  ylab("Expected") +
  theme(plot.title = element_text(hjust=0.5, size = 16), axis.title = element_text(size = 14), legend.position = "none")

f4 <- ggplot(refDat_lost2020_2, aes(x=observed, y=expected, col=g_whoregion, label=iso3))+
  geom_point(size=4) +
  #geom_text(label=ifelse(refDat_lost2020_2$pvalue<0.05,refDat_lost2020_2$iso3,""), hjust=0.75, vjust=0, size=7) +
  geom_text_repel(aes(label=ifelse(refDat_lost2020_2$RRpvalue<0.05,refDat_lost2020_2$iso3,"")), size=6, max.overlaps = 15)+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_hline(yintercept=0, linetype="solid", colour="black", size=0.1)+
  geom_vline(xintercept=0, linetype="solid", colour="black", size=0.1)+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.35)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.35)) +
  labs(color="WHO Region") +
  ggtitle("Loss to follow up")+
  xlab("Observed") +
  ylab("Expected") +
  theme(plot.title = element_text(hjust=0.5, size = 16), axis.title = element_text(size = 14),legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

legend <- get_legend(f4)

f4 <- f4 + theme(legend.position = "none")

g <- arrangeGrob(f1,f2,legend, f3,f4, ncol=3, widths=c(3.0,3.0,2.0))

ggsave("Exp_Obs_splines_boot_df3_2020_20240223.png",g,device="png",width=38,height=30,units=c("cm"))

###################################################################################################################################################

#plot of observed vs expected proportions for 2021



ggplot(datFullLong_sp_boot %>% filter(year==2021), aes(x=observed, y=expected, col=g_whoregion, label=iso3))+
  geom_point() +
  geom_text(check_overlap = TRUE) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  facet_wrap(~paste(sep=" ",Outcome),ncol = 2, scales = "free") +
  labs(color="WHO Region") +
  xlab("Observed") +
  ylab("Expected") +
  theme(axis.title = element_text(size = 12), legend.text = element_text(size = 10),
        legend.title = element_text(size = 12), strip.text.x = element_text(size = 12))



g1<- ggplot(refDat_succ2021_2, aes(x=observed, y=expected, col=g_whoregion, label=iso3))+
  geom_point(size=3) +
  #geom_text(label=ifelse(refDat_succ2021_2$pvalue<0.05,refDat_succ2021_2$iso3,""), hjust=0.75, vjust=0, size=5) +
  geom_text_repel(aes(label=ifelse(refDat_succ2021_2$RRpvalue<0.05,refDat_succ2021_2$iso3,"")), size=6)+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_hline(yintercept=0, linetype="solid", colour="black", size=0.1)+
  geom_vline(xintercept=0, linetype="solid", colour="black", size=0.1)+
  scale_x_continuous(expand = c(0, 0), limits = c(0.7, 1)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0.7, 1)) +
  labs(color="WHO Region") +
  ggtitle("Success")+
  xlab("Observed") +
  ylab("Expected") +
  theme(plot.title = element_text(hjust=0.5, size=16), axis.title = element_text(size = 14), legend.position = "none")
#removing Liberia for the plot - it has a significant pvalue but RR and CI are zeros
#removing UZB and KAZ - no data
refDat_fail2021_2_2<- refDat_fail2021_2 %>%
  filter(!iso3=="LBR")
refDat_fail2021_2_2 <- refDat_fail2021_2_2 %>%
  filter(!iso3==c("UZB", "KAZ"))

g2 <-ggplot(refDat_fail2021_2_2, aes(x=observed, y=expected, col=g_whoregion, label=iso3))+
  geom_point(size=3) +
  #geom_text(label=ifelse(refDat_fail2021_2$pvalue<0.05,refDat_fail2021_2$iso3,""), hjust=0.75, vjust=0, size=6) +
  geom_text_repel(aes(label=ifelse(refDat_fail2021_2_2$RRpvalue<0.05,refDat_fail2021_2_2$iso3,"")), size=6, max.overlaps = 40)+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_hline(yintercept=0, linetype="solid", colour="black", size=0.1)+
  geom_vline(xintercept=0, linetype="solid", colour="black", size=0.1)+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.15)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.15)) +
  labs(color="WHO Region") +
  ggtitle("Failure")+
  xlab("Observed") +
  ylab("Expected") +
  theme(plot.title = element_text(hjust=0.5, size=16), axis.title = element_text(size = 14), legend.position = "none")


g3 <- ggplot(refDat_died2021_2, aes(x=observed, y=expected, col=g_whoregion, label=iso3))+
  geom_point(size=3) +
  #geom_text(label=ifelse(refDat_died2021_2$pvalue<0.05,refDat_died2021_2$iso3,""), hjust=0.75, vjust=0, size=5) +
  geom_text_repel(aes(label=ifelse(refDat_died2021_2$RRpvalue<0.05,refDat_died2021_2$iso3,"")), size=6)+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_hline(yintercept=0, linetype="solid", colour="black", size=0.1)+
  geom_vline(xintercept=0, linetype="solid", colour="black", size=0.1)+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.18)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.18)) +
  labs(color="WHO Region") +
  ggtitle("Death")+
  xlab("Observed") +
  ylab("Expected") +
  theme(plot.title = element_text(hjust=0.5, size = 16), axis.title = element_text(size = 14), legend.position = "none")

g4 <- ggplot(refDat_lost2021_2, aes(x=observed, y=expected, col=g_whoregion, label=iso3))+
  geom_point(size=3) +
  #geom_text(label=ifelse(refDat_lost2021_2$pvalue<0.05,refDat_lost2021_2$iso3,""), hjust=0.75, vjust=0, size=6) +
  geom_text_repel(aes(label=ifelse(refDat_lost2021_2$RRpvalue<0.05,refDat_lost2021_2$iso3,"")), size=6)+
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  geom_hline(yintercept=0, linetype="solid", colour="black", size=0.1)+
  geom_vline(xintercept=0, linetype="solid", colour="black", size=0.1)+
  scale_x_continuous(expand = c(0, 0), limits = c(0, 0.35)) + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, 0.35)) +
  labs(color="WHO Region") +
  ggtitle("Loss to follow up")+
  xlab("Observed") +
  ylab("Expected") +
  theme(plot.title = element_text(hjust=0.5, size = 16), axis.title = element_text(size = 14), legend.text = element_text(size = 14),
        legend.title = element_text(size = 14))

legend2 <- get_legend(g4)

g4 <- g4 + theme(legend.position = "none")

g_2021 <- arrangeGrob(g1,g2,legend2, g3,g4, ncol=3, widths=c(3.0,3.0,2.0))

ggsave("Exp_Obs_splines_boot_df3_2021_20240223.png",g_2021,device="png",width=38,height=30,units=c("cm"))


#########################################################################################################################################################