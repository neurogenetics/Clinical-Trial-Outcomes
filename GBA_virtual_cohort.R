##### Part of a series of analyses associated with the manuscript "Genetic variability and potential effects on clinical
##### trial outcomes: persepctives in Parkinson's disease" 
##### Hampton Leonard 09/26/18


##### Create virtual GBA cohort with effect sizes from study


### Virtual Cohorts

### 60000 global pop

library(tidyr)
library(dplyr)
library(car)
library(ggplot2)
library(ggthemes)

set.seed(1234567)

cohort1 <- data.frame("ID" = numeric(60000), "rs76763715" = numeric(60000), "rs2230288" = numeric(60000), "rs75548401" = numeric(60000), "SCORE" = numeric(60000))

cohort1$ID <- 1:60000

## SNP rs76763715 aka N370S
## Assign variant frequencies according to ExAC Release 1 database and HWE

p <- 0.9966
q <- 0.0034

## HWE

TT_freq <- p^2

TC_freq <- 2*p*q

CC_freq <- q^2

TT_pop <- round(TT_freq * 60000)-1

TC_pop <- round(TC_freq * 60000)

CC_pop <- round(CC_freq * 60000)


samp1 <- sample_n(cohort1, TT_pop)
samp1$rs76763715 <- 0

samp2 <- sample_n(cohort1[cohort1$ID %in% samp1$ID == FALSE,], TC_pop)
samp2$rs76763715 <- 1

sampa <- sample_n(cohort1[cohort1$ID %in% samp1$ID == FALSE & cohort1$ID %in% samp2$ID == FALSE,], CC_pop)
sampa$rs76763715 <- 2


samp3 <- rbind(samp1, samp2, sampa)

samp3 <- samp3[,-c(3,4,5)]


## rs2230288 aka E326K or E365K

p <- 0.9852
q <- 0.0148


AA_freq <- p^2

AG_freq <- 2*p*q

GG_freq <- q^2

AA_pop <- round(AA_freq * 60000) 

AG_pop <- round(AG_freq * 60000)

GG_pop <- round(GG_freq * 60000)


samp4 <- sample_n(cohort1, AA_pop)
samp4$rs2230288 <- 0

samp5 <- sample_n(cohort1[cohort1$ID %in% samp4$ID == FALSE,], AG_pop)
samp5$rs2230288 <- 1

sampb <- sample_n(cohort1[cohort1$ID %in% samp4$ID == FALSE & cohort1$ID %in% samp5$ID == FALSE,], GG_pop)
sampb$rs2230288 <- 2


samp6 <- rbind(samp4, samp5, sampb)

samp6 <- samp6[, -c(2,4,5)]



## SNP rs75548401 aka K26R

p <- 0.9902
q <- 0.0098


AA_freq <- p^2

AG_freq <- 2*p*q

GG_freq <- q^2

AA_pop <- round(AA_freq * 60000) 

AG_pop <- round(AG_freq * 60000)

GG_pop <- round(GG_freq * 60000)


samp7 <- sample_n(cohort1, AA_pop)
samp7$rs75548401 <- 0

samp8 <- sample_n(cohort1[cohort1$ID %in% samp7$ID == FALSE,], AG_pop)
samp8$rs75548401 <- 1

sampc <- sample_n(cohort1[cohort1$ID %in% samp7$ID == FALSE & cohort1$ID %in% samp8$ID == FALSE,], GG_pop)
sampc$rs75548401 <- 2


samp9 <- rbind(samp7, samp8, sampc)


samp9 <- samp9[, -c(2,3)]



cohort_c <- merge(samp3, samp6)
cohort_c <- merge(cohort_c, samp9)



## Assign effects from study betas

for(i in 1:nrow(cohort_c)){
  cohort_c$SCORE[i] <- (cohort_c$rs76763715[i] * 0.7467) + (cohort_c$rs2230288[i] * 0.6357) +(cohort_c$rs75548401[i] * 0.3619)
}


## Remove any virtual patients that don't have one of the GBA mutations

cohort_c <- cohort_c[cohort_c$SCORE != 0,]

predictors <- names(cohort_c)[c(2,3,4)]


set.seed(1234567)

nIterations <- 1000 #test at 1000 iterations
# nIterations <- 10 #test at 1000 iterations
for(z in c(50,100,200,300,600,700,800,1000,2000))
#for(z in c(200))
{
  GBA_scores <- data.frame(variance_treated = numeric(nIterations), variance_untreated = numeric(nIterations), variance_sig = numeric(nIterations),
                         mean_treated_GRS = numeric(nIterations), mean_untreated_GRS = numeric(nIterations), z_score = numeric(nIterations))
                         
  
  output <- matrix(ncol=length(predictors), nrow=nIterations)
  for(y in 1:nIterations)
  {
    temp <- sample_n(cohort_c, z)
    
    ## randomly assign to trial arms for unbalanced cohorts, comment out if balancing
    
    temp$pheno <- rep(0:1, length.out = length(temp$ID))
    
    
    ## use below for balancing, comment out to line 211 if not balancing
    
    # temp1 <- temp[sample(nrow(temp), 1),]
    # temp1$pheno <- 0
    # 
    # temp <- temp[!(temp$ID %in% temp1$ID),]
    # 
    # temp2 <- temp[sample(nrow(temp), 1),]
    # temp2$pheno <- 1
    # 
    # temp3 <- temp[!(temp$ID %in% temp2$ID),]
    # 
    # temp3$pheno <- "NA"
    # 
    # 
    # for(i in 1:nrow(temp3)){
    #   test1 <- rbind(temp3[i,2:4],temp1[,2:4])
    #   sums1 <- colSums(test1[,1:3], na.rm = T)
    #   sumstemp1 <- colSums(temp2[2:4], na.rm = T)
    #   
    #   diff1 <- abs(sums1 - sumstemp1)
    #   final1 <- sum(diff1)
    #   
    #   test2 <- rbind(temp3[i,2:4],temp2[,2:4])
    #   sums2 <- colSums(test2[,1:3], na.rm = T)
    #   sumstemp2 <- colSums(temp1[2:4], na.rm = T)
    #   
    #   diff2 <- abs(sums2 - sumstemp2)
    #   final2 <- sum(diff2)
    #   
    #   if(final1 < final2){
    #     temp1 <- rbind(temp1, temp3[i,])
    #     temp1$pheno <- 0
    #   }else{
    #     temp2 <- rbind(temp2, temp3[i,])
    #     temp2$pheno <- 1
    #   }  
    # }
    # 
    # 
    # balanced_final <- rbind(temp1, temp2)
    # 
    # temp <- balanced_final
    
    
    
    
      print(paste(" at y iterations = ", y, " for sampling N =", z, sep = ""))
     
      
      ## test for equality of variance, testing for where trial arms have significantly different
      ## variance in GRS
      
      levene <- leveneTest(SCORE ~ factor(pheno), data = temp)
      lp <- levene$`Pr(>F)`[1]
      
      
      ## Split into simulated treatment and placebo arms
      
      untreated <- temp[temp$pheno == 0,]
      treated <- temp[temp$pheno == 1,]
      
      
      ut_mean <- mean(untreated$SCORE)
      t_mean <- mean(treated$SCORE)
      
      
      
      GBA_scores[y,1] <- var(treated$SCORE)
      GBA_scores[y,2] <- var(untreated$SCORE)
      GBA_scores[y,3] <- lp
      GBA_scores[y,4] <- mean(treated$SCORE)
      GBA_scores[y,5] <- mean(untreated$SCORE)
      GBA_scores[y,6] <- (ut_mean - t_mean)/sqrt(((sd(untreated$SCORE)^2))/nrow(untreated) + ((sd(treated$SCORE)^2))/nrow(treated))
      
    
    
  }
  write.table(GBA_scores, paste("GBA_iterations",nIterations,"sampling",z,".tab",sep = ""), quote = F, sep = "\t")  
}


### Analysis for signficant iterations - read in tables and subset by Z cut-off score  

gba_s50 <- read.table("GBA_iterations1000sampling50.tab")
gba_s50_sig <- gba_s50[gba_s50$z_score > 1.96 | gba_s50$z_score < -1.96,]

gba_s100 <- read.table("GBA_iterations1000sampling100.tab")
gba_s100_sig <- gba_s100[gba_s100$z_score > 1.96 | gba_s100$z_score < -1.96,]

gba_s200 <- read.table("GBA_iterations1000sampling200.tab")
gba_s200_sig <- gba_s200[gba_s200$z_score > 1.96 | gba_s200$z_score < -1.96,]

gba_s300 <- read.table("GBA_iterations1000sampling300.tab")
gba_s300_sig <- gba_s300[gba_s300$z_score > 1.96 | gba_s300$z_score < -1.96,]

gba_s600 <- read.table("GBA_iterations1000sampling600.tab")
gba_s600_sig <- gba_s600[gba_s600$z_score > 1.96 | gba_s600$z_score < -1.96,]

gba_s700 <- read.table("GBA_iterations1000sampling700.tab")
gba_s700_sig <- gba_s700[gba_s700$z_score > 1.96 | gba_s700$z_score < -1.96,]

gba_s800 <- read.table("GBA_iterations1000sampling800.tab")
gba_s800_sig <- gba_s800[gba_s800$z_score > 1.96 | gba_s800$z_score < -1.96,]

gba_s1000 <- read.table("GBA_iterations1000sampling1000.tab")
gba_s1000_sig <- gba_s1000[gba_s1000$z_score > 1.96 | gba_s1000$z_score < -1.96,]

gba_s2000 <- read.table("GBA_iterations1000sampling2000.tab")
gba_s2000_sig <- gba_s2000[gba_s2000$z_score > 1.96 | gba_s2000$z_score < -1.96,]



## Calculate percent difference in mean GRS between arms for significant iteratios

df_list <- list(gba_s50_sig =gba_s50_sig, gba_s100_sig = gba_s100_sig, gba_s200_sig = gba_s200_sig, gba_s300_sig = gba_s300_sig, 
                gba_s600_sig = gba_s600_sig, gba_s700_sig = gba_s700_sig, gba_s800_sig = gba_s800_sig, gba_s1000_sig = gba_s1000_sig,
                gba_s2000_sig = gba_s2000_sig)


## Function for calculating percent difference in mean between arms

change_fnc <- function(x){
  
  for(i in 1:nrow(x)){
    
    x$percent_diff_GRS[i] <- abs(((x[,5][i] - x[,4][i])/((x[,5][i] + x[,4][i])/2)) * 100)
  }
  return(x)
}



frames <- lapply(df_list, FUN = change_fnc)

gba_s50_sig <- as.data.frame(frames[1])
colnames(gba_s50_sig) <- sub(".*\\.", "", colnames(gba_s50_sig))

gba_s100_sig <- as.data.frame(frames[2])
colnames(gba_s100_sig) <- sub(".*\\.", "", colnames(gba_s100_sig))

gba_s200_sig <- as.data.frame(frames[3])
colnames(gba_s200_sig) <- sub(".*\\.", "", colnames(gba_s200_sig))

gba_s300_sig <- as.data.frame(frames[4])
colnames(gba_s300_sig) <- sub(".*\\.", "", colnames(gba_s300_sig))

gba_s600_sig <- as.data.frame(frames[5])
colnames(gba_s600_sig) <- sub(".*\\.", "", colnames(gba_s600_sig))

gba_s700_sig <- as.data.frame(frames[6])
colnames(gba_s700_sig) <- sub(".*\\.", "", colnames(gba_s700_sig))

gba_s800_sig <- as.data.frame(frames[7])
colnames(gba_s800_sig) <- sub(".*\\.", "", colnames(gba_s800_sig))

gba_s1000_sig <- as.data.frame(frames[8])
colnames(gba_s1000_sig) <- sub(".*\\.", "", colnames(gba_s1000_sig))

gba_s2000_sig <- as.data.frame(frames[9])
colnames(gba_s2000_sig) <- sub(".*\\.", "", colnames(gba_s2000_sig))




## create dataframe of average percent difference in GRS for each sample size

diff_df <- data.frame(Sample_Size = numeric(9), Avg_mean_GRS_Diff = numeric(9))
diff_df$Sample_Size <- c(50, 100, 200, 300, 600, 700, 800, 1000, 2000)
diff_df$Avg_mean_GRS_Diff <- c(mean(gba_s50_sig$percent_diff_GRS), mean(gba_s100_sig$percent_diff_GRS),
                               mean(gba_s200_sig$percent_diff_GRS), mean(gba_s300_sig$percent_diff_GRS),
                               mean(gba_s600_sig$percent_diff_GRS),mean(gba_s700_sig$percent_diff_GRS), 
                               mean(gba_s800_sig$percent_diff_GRS),mean(gba_s1000_sig$percent_diff_GRS),
                               mean(gba_s2000_sig$percent_diff_GRS))




diff_df$Avg_mean_GRS_Diff <- round(diff_df$Avg_mean_GRS_Diff, 2)

## Confidence Intervals

t <- t.test(gba_s50_sig$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df$Avg_mean_GRS_Diff[1] <- paste(diff_df$Avg_mean_GRS_Diff[1], "(", t1, "-", t2, ")" )


t <- t.test(gba_s100_sig$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df$Avg_mean_GRS_Diff[2] <- paste(diff_df$Avg_mean_GRS_Diff[2], "(", t1, "-", t2, ")" )


t <- t.test(gba_s200_sig$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df$Avg_mean_GRS_Diff[3] <- paste(diff_df$Avg_mean_GRS_Diff[3], "(", t1, "-", t2, ")" )


t <- t.test(gba_s300_sig$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df$Avg_mean_GRS_Diff[4] <- paste(diff_df$Avg_mean_GRS_Diff[4], "(", t1, "-", t2, ")" )


t <- t.test(gba_s600_sig$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df$Avg_mean_GRS_Diff[5] <- paste(diff_df$Avg_mean_GRS_Diff[5], "(", t1, "-", t2, ")" )


t <- t.test(gba_s700_sig$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df$Avg_mean_GRS_Diff[6] <- paste(diff_df$Avg_mean_GRS_Diff[6], "(", t1, "-", t2, ")" )


t <- t.test(gba_s800_sig$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df$Avg_mean_GRS_Diff[7] <- paste(diff_df$Avg_mean_GRS_Diff[7], "(", t1, "-", t2, ")" )


t <- t.test(gba_s1000_sig$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df$Avg_mean_GRS_Diff[8] <- paste(diff_df$Avg_mean_GRS_Diff[8], "(", t1, "-", t2, ")" )


t <- t.test(gba_s2000_sig$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df$Avg_mean_GRS_Diff[9] <- paste(diff_df$Avg_mean_GRS_Diff[9], "(", t1, "-", t2, ")" )



## Repeat for all iterations 

df_list_all <- list(gba_s50 =gba_s50, gba_s100 = gba_s100, gba_s200 = gba_s200, gba_s300 = gba_s300, 
                gba_s600 = gba_s600, gba_s700 = gba_s700, gba_s800 = gba_s800, gba_s1000 = gba_s1000,
                gba_s2000 = gba_s2000)



frames <- lapply(df_list_all, FUN = change_fnc)

gba_s50 <- as.data.frame(frames[1])
colnames(gba_s50) <- sub(".*\\.", "", colnames(gba_s50))

gba_s100 <- as.data.frame(frames[2])
colnames(gba_s100) <- sub(".*\\.", "", colnames(gba_s100))

gba_s200 <- as.data.frame(frames[3])
colnames(gba_s200) <- sub(".*\\.", "", colnames(gba_s200))

gba_s300 <- as.data.frame(frames[4])
colnames(gba_s300) <- sub(".*\\.", "", colnames(gba_s300))

gba_s600 <- as.data.frame(frames[5])
colnames(gba_s600) <- sub(".*\\.", "", colnames(gba_s600))

gba_s700 <- as.data.frame(frames[6])
colnames(gba_s700) <- sub(".*\\.", "", colnames(gba_s700))

gba_s800 <- as.data.frame(frames[7])
colnames(gba_s800) <- sub(".*\\.", "", colnames(gba_s800))

gba_s1000 <- as.data.frame(frames[8])
colnames(gba_s1000) <- sub(".*\\.", "", colnames(gba_s1000))

gba_s2000 <- as.data.frame(frames[9])
colnames(gba_s2000) <- sub(".*\\.", "", colnames(gba_s2000))


## create dataframe of average percent difference in GRS for each sample size

diff_df_all <- data.frame(Sample_Size = numeric(9), Avg_mean_GRS_Diff = numeric(9))
diff_df_all$Sample_Size <- c(50, 100, 200, 300, 600, 700, 800, 1000, 2000)
diff_df_all$Avg_mean_GRS_Diff <- c(mean(gba_s50$percent_diff_GRS), mean(gba_s100$percent_diff_GRS),
                               mean(gba_s200$percent_diff_GRS), mean(gba_s300$percent_diff_GRS),
                               mean(gba_s600$percent_diff_GRS),mean(gba_s700$percent_diff_GRS), 
                               mean(gba_s800$percent_diff_GRS),mean(gba_s1000$percent_diff_GRS),
                               mean(gba_s2000$percent_diff_GRS))



diff_df_all$Avg_mean_GRS_Diff <- round(diff_df_all$Avg_mean_GRS_Diff, 2)

## Confidence Intervals

t <- t.test(gba_s50$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_all$Avg_mean_GRS_Diff[1] <- paste(diff_df_all$Avg_mean_GRS_Diff[1], "(", t1, "-", t2, ")" )


t <- t.test(gba_s100$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_all$Avg_mean_GRS_Diff[2] <- paste(diff_df_all$Avg_mean_GRS_Diff[2], "(", t1, "-", t2, ")" )


t <- t.test(gba_s200$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_all$Avg_mean_GRS_Diff[3] <- paste(diff_df_all$Avg_mean_GRS_Diff[3], "(", t1, "-", t2, ")" )


t <- t.test(gba_s300$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_all$Avg_mean_GRS_Diff[4] <- paste(diff_df_all$Avg_mean_GRS_Diff[4], "(", t1, "-", t2, ")" )


t <- t.test(gba_s600$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_all$Avg_mean_GRS_Diff[5] <- paste(diff_df_all$Avg_mean_GRS_Diff[5], "(", t1, "-", t2, ")" )


t <- t.test(gba_s700$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_all$Avg_mean_GRS_Diff[6] <- paste(diff_df_all$Avg_mean_GRS_Diff[6], "(", t1, "-", t2, ")" )


t <- t.test(gba_s800$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_all$Avg_mean_GRS_Diff[7] <- paste(diff_df_all$Avg_mean_GRS_Diff[7], "(", t1, "-", t2, ")" )


t <- t.test(gba_s1000$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_all$Avg_mean_GRS_Diff[8] <- paste(diff_df_all$Avg_mean_GRS_Diff[8], "(", t1, "-", t2, ")" )


t <- t.test(gba_s2000$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_all$Avg_mean_GRS_Diff[9] <- paste(diff_df_all$Avg_mean_GRS_Diff[9], "(", t1, "-", t2, ")" )


## min and max for each sample size

diff_df_all$min[1] <- round(min(gba_s50$percent_diff_GRS), 4)
diff_df_all$max[1] <- round(max(gba_s50$percent_diff_GRS), 4)

diff_df_all$min[2] <- round(min(gba_s100$percent_diff_GRS), 4)
diff_df_all$max[2] <- round(max(gba_s100$percent_diff_GRS), 4)

diff_df_all$min[3] <- round(min(gba_s200$percent_diff_GRS), 4)
diff_df_all$max[3] <- round(max(gba_s200$percent_diff_GRS), 4)

diff_df_all$min[4] <- round(min(gba_s300$percent_diff_GRS), 4)
diff_df_all$max[4] <- round(max(gba_s300$percent_diff_GRS), 4)

diff_df_all$min[5] <- round(min(gba_s600$percent_diff_GRS), 4)
diff_df_all$max[5] <- round(max(gba_s600$percent_diff_GRS), 4)

diff_df_all$min[6] <- round(min(gba_s700$percent_diff_GRS), 4)
diff_df_all$max[6] <- round(max(gba_s700$percent_diff_GRS), 4)

diff_df_all$min[7] <- round(min(gba_s800$percent_diff_GRS), 4)
diff_df_all$max[7] <- round(max(gba_s800$percent_diff_GRS), 4)

diff_df_all$min[8] <- round(min(gba_s1000$percent_diff_GRS), 4)
diff_df_all$max[8] <- round(max(gba_s1000$percent_diff_GRS), 4)

diff_df$min[9] <- round(min(gba_s2000$percent_diff_GRS), 4)
diff_df$max[9] <- round(max(gba_s2000$percent_diff_GRS), 4)




## Fit multipolynomial regression line to data for information on relationship between sample size and percent difference in mean GRS
## LOESS for visual inspection of relationship and figures 



a <- c(mean(gba_s50_sig$percent_diff_GRS), mean(gba_s100_sig$percent_diff_GRS),
       mean(gba_s200_sig$percent_diff_GRS), mean(gba_s300_sig$percent_diff_GRS),
       mean(gba_s600_sig$percent_diff_GRS),mean(gba_s700_sig$percent_diff_GRS), 
       mean(gba_s800_sig$percent_diff_GRS),mean(gba_s1000_sig$percent_diff_GRS),
       mean(gba_s2000_sig$percent_diff_GRS))

b <- diff_df$Sample_Size


testlm <- lm (a ~ poly(b,2))

summary(testlm)

ggplot(diff_df, aes(x = b, y = a)) + geom_point(size = 3) + stat_smooth(method = "loess", span = 1.17, formula =y ~ x,  col = "red")+
  labs(x = "Sample Size", y = "Average % Difference", title = "Average Percent Difference in GRS Between GBA Treatment Groups (Significant Iterations)")+
  theme_tufte()+ geom_rangeframe()+theme(text = element_text(size = 15)) 









#########################################################################
## Supplemental analyses: variance in GRS between arms

gba_s50_v <- read.table("GBA_iterations1000sampling50.tab")
gba_s50_sig_v <- gba_s50_v[gba_s50_v$variance_sig < 0.05,]

gba_s100_v <- read.table("GBA_iterations1000sampling100.tab")
gba_s100_sig_v <- gba_s100_v[gba_s100_v$variance_sig < 0.05,]

gba_s200_v <- read.table("GBA_iterations1000sampling200.tab")
gba_s200_sig_v <- gba_s200_v[gba_s200_v$variance_sig < 0.05,]

gba_s300_v <- read.table("GBA_iterations1000sampling300.tab")
gba_s300_sig_v <- gba_s300_v[gba_s300_v$variance_sig < 0.05,]

gba_s600_v <- read.table("GBA_iterations1000sampling600.tab")
gba_s600_sig_v <- gba_s600_v[gba_s600_v$variance_sig < 0.05,]

gba_s700_v <- read.table("GBA_iterations1000sampling700.tab")
gba_s700_sig_v <- gba_s700_v[gba_s700_v$variance_sig < 0.05,]

gba_s800_v <- read.table("GBA_iterations1000sampling800.tab")
gba_s800_sig_v <- gba_s800_v[gba_s800_v$variance_sig < 0.05,]

gba_s1000_v <- read.table("GBA_iterations1000sampling1000.tab")
gba_s1000_sig_v <- gba_s1000_v[gba_s1000_v$variance_sig < 0.05,]

gba_s2000_v <- read.table("GBA_iterations1000sampling2000.tab")
gba_s2000_sig_v <- gba_s2000_v[gba_s2000_v$variance_sig < 0.05,]


## Significant Iterations

df_list_v <- list(gba_s50_sig_v =gba_s50_sig_v, gba_s100_sig_v = gba_s100_sig_v, gba_s200_sig_v = gba_s200_sig_v, gba_s300_sig_v = gba_s300_sig_v, 
                  gba_s600_sig_v = gba_s600_sig_v, gba_s700_sig_v = gba_s700_sig_v, gba_s800_sig_v = gba_s800_sig_v, gba_s1000_sig_v = gba_s1000_sig_v,
                  gba_s2000_sig_v = gba_s2000_sig_v)



## Function for calculating percent difference in variance

change_fnc_v <- function(x){
  
  
  for(i in 1:nrow(x)){
    
    x$percent_diff_GRS_Variance[i] <- abs(((x[,2][i] - x[,1][i])/((x[,2][i] + x[,1][i])/2)) * 100)
  }
  return(x)
}


frames <- lapply(df_list_v, FUN = change_fnc_v)

gba_s50_sig_v <- as.data.frame(frames[1])
colnames(gba_s50_sig_v) <- sub(".*\\.", "", colnames(gba_s50_sig_v))

gba_s100_sig_v <- as.data.frame(frames[2])
colnames(gba_s100_sig_v) <- sub(".*\\.", "", colnames(gba_s100_sig_v))

gba_s200_sig_v <- as.data.frame(frames[3])
colnames(gba_s200_sig_v) <- sub(".*\\.", "", colnames(gba_s200_sig_v))

gba_s300_sig_v <- as.data.frame(frames[4])
colnames(gba_s300_sig_v) <- sub(".*\\.", "", colnames(gba_s300_sig_v))

gba_s600_sig_v <- as.data.frame(frames[5])
colnames(gba_s600_sig_v) <- sub(".*\\.", "", colnames(gba_s600_sig_v))

gba_s700_sig_v <- as.data.frame(frames[6])
colnames(gba_s700_sig_v) <- sub(".*\\.", "", colnames(gba_s700_sig_v))

gba_s800_sig_v <- as.data.frame(frames[7])
colnames(gba_s800_sig_v) <- sub(".*\\.", "", colnames(gba_s800_sig_v))

gba_s1000_sig_v <- as.data.frame(frames[8])
colnames(gba_s1000_sig_v) <- sub(".*\\.", "", colnames(gba_s1000_sig_v))

gba_s2000_sig_v <- as.data.frame(frames[9])
colnames(gba_s2000_sig_v) <- sub(".*\\.", "", colnames(gba_s2000_sig_v))




diff_df_var <- data.frame(Sample_Size = numeric(9), Avg_var_GRS_Diff = numeric(9))
diff_df_var$Sample_Size <- c(50, 100, 200, 300, 600, 700, 800, 1000, 2000)
diff_df_var$Avg_var_GRS_Diff <- c(mean(gba_s50_sig_v$percent_diff_GRS_Variance), mean(gba_s100_sig_v$percent_diff_GRS_Variance),
                              mean(gba_s200_sig_v$percent_diff_GRS_Variance), mean(gba_s300_sig_v$percent_diff_GRS_Variance),
                              mean(gba_s600_sig_v$percent_diff_GRS_Variance), mean(gba_s700_sig_v$percent_diff_GRS_Variance),
                              mean(gba_s800_sig_v$percent_diff_GRS_Variance), mean(gba_s1000_sig_v$percent_diff_GRS_Variance),
                              mean(gba_s2000_sig_v$percent_diff_GRS_Variance))




diff_df_var$Avg_var_GRS_Diff <- round(diff_df_var$Avg_var_GRS_Diff, 2)

t <- t.test(gba_s50_sig_v$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_var$Avg_var_GRS_Diff[1] <- paste(diff_df_var$Avg_var_GRS_Diff[1], "(", t1, "-", t2, ")" )


t <- t.test(gba_s100_sig_v$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_var$Avg_var_GRS_Diff[2] <- paste(diff_df_var$Avg_var_GRS_Diff[2], "(", t1, "-", t2, ")" )


t <- t.test(gba_s200_sig_v$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_var$Avg_var_GRS_Diff[3] <- paste(diff_df_var$Avg_var_GRS_Diff[3], "(", t1, "-", t2, ")" )


t <- t.test(gba_s300_sig_v$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_var$Avg_var_GRS_Diff[4] <- paste(diff_df_var$Avg_var_GRS_Diff[4], "(", t1, "-", t2, ")" )


t <- t.test(gba_s600_sig_v$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_var$Avg_var_GRS_Diff[5] <- paste(diff_df_var$Avg_var_GRS_Diff[5], "(", t1, "-", t2, ")" )


t <- t.test(gba_s700_sig_v$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_var$Avg_var_GRS_Diff[6] <- paste(diff_df_var$Avg_var_GRS_Diff[6], "(", t1, "-", t2, ")" )


t <- t.test(gba_s800_sig_v$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_var$Avg_var_GRS_Diff[7] <- paste(diff_df_var$Avg_var_GRS_Diff[7], "(", t1, "-", t2, ")" )


t <- t.test(gba_s1000_sig_v$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_var$Avg_var_GRS_Diff[8] <- paste(diff_df_var$Avg_var_GRS_Diff[8], "(", t1, "-", t2, ")" )


t <- t.test(gba_s2000_sig_v$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_var$Avg_var_GRS_Diff[9] <- paste(diff_df_var$Avg_var_GRS_Diff[9], "(", t1, "-", t2, ")" )



## Calculate percent difference in variance for all iterations 

df_list_v_all <- list(gba_s50_v =gba_s50_v, gba_s100_v = gba_s100_v, gba_s200_v = gba_s200_v, gba_s300_v = gba_s300_v, 
                  gba_s600_v = gba_s600_v, gba_s700_v = gba_s700_v, gba_s800_v = gba_s800_v, gba_s1000_v = gba_s1000_v,
                  gba_s2000_v = gba_s2000_v)



frames <- lapply(df_list_v_all, FUN = change_fnc_v)

gba_s50_v <- as.data.frame(frames[1])
colnames(gba_s50_v) <- sub(".*\\.", "", colnames(gba_s50_v))

gba_s100_v <- as.data.frame(frames[2])
colnames(gba_s100_v) <- sub(".*\\.", "", colnames(gba_s100_v))

gba_s200_v <- as.data.frame(frames[3])
colnames(gba_s200_v) <- sub(".*\\.", "", colnames(gba_s200_v))

gba_s300_v <- as.data.frame(frames[4])
colnames(gba_s300_v) <- sub(".*\\.", "", colnames(gba_s300_v))

gba_s600_v <- as.data.frame(frames[5])
colnames(gba_s600_v) <- sub(".*\\.", "", colnames(gba_s600_v))

gba_s700_v <- as.data.frame(frames[6])
colnames(gba_s700_v) <- sub(".*\\.", "", colnames(gba_s700_v))

gba_s800_v <- as.data.frame(frames[7])
colnames(gba_s800_v) <- sub(".*\\.", "", colnames(gba_s800_v))

gba_s1000_v <- as.data.frame(frames[8])
colnames(gba_s1000_v) <- sub(".*\\.", "", colnames(gba_s1000_v))

gba_s2000_v <- as.data.frame(frames[9])
colnames(gba_s2000_v) <- sub(".*\\.", "", colnames(gba_s2000_v))




## create dataframe of average percent difference in variance in GRS for each sample size

diff_df_var_all <- data.frame(Sample_Size = numeric(9), Avg_var_GRS_Diff = numeric(9))
diff_df_var_all$Sample_Size <- c(50, 100, 200, 300, 600, 700, 800, 1000, 2000)
diff_df_var_all$Avg_var_GRS_Diff <- c(mean(gba_s50_v$percent_diff_GRS_Variance), mean(gba_s100_v$percent_diff_GRS_Variance),
                                  mean(gba_s200_v$percent_diff_GRS_Variance), mean(gba_s300_v$percent_diff_GRS_Variance),
                                  mean(gba_s600_v$percent_diff_GRS_Variance), mean(gba_s700_v$percent_diff_GRS_Variance),
                                  mean(gba_s800_v$percent_diff_GRS_Variance), mean(gba_s1000_v$percent_diff_GRS_Variance),
                                  mean(gba_s2000_v$percent_diff_GRS_Variance))




diff_df_var_all$Avg_var_GRS_Diff <- round(diff_df_var_all$Avg_var_GRS_Diff, 2)

t <- t.test(gba_s50_v$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_var_all$Avg_var_GRS_Diff[1] <- paste(diff_df_var_all$Avg_var_GRS_Diff[1], "(", t1, "-", t2, ")" )


t <- t.test(gba_s100_v$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_var_all$Avg_var_GRS_Diff[2] <- paste(diff_df_var_all$Avg_var_GRS_Diff[2], "(", t1, "-", t2, ")" )


t <- t.test(gba_s200_v$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_var_all$Avg_var_GRS_Diff[3] <- paste(diff_df_var_all$Avg_var_GRS_Diff[3], "(", t1, "-", t2, ")" )


t <- t.test(gba_s300_v$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_var_all$Avg_var_GRS_Diff[4] <- paste(diff_df_var_all$Avg_var_GRS_Diff[4], "(", t1, "-", t2, ")" )


t <- t.test(gba_s600_v$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_var_all$Avg_var_GRS_Diff[5] <- paste(diff_df_var_all$Avg_var_GRS_Diff[5], "(", t1, "-", t2, ")" )


t <- t.test(gba_s700_v$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_var_all$Avg_var_GRS_Diff[6] <- paste(diff_df_var_all$Avg_var_GRS_Diff[6], "(", t1, "-", t2, ")" )


t <- t.test(gba_s800_v$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_var_all$Avg_var_GRS_Diff[7] <- paste(diff_df_var_all$Avg_var_GRS_Diff[7], "(", t1, "-", t2, ")" )


t <- t.test(gba_s1000_v$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_var_all$Avg_var_GRS_Diff[8] <- paste(diff_df_var_all$Avg_var_GRS_Diff[8], "(", t1, "-", t2, ")" )


t <- t.test(gba_s2000_v$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_var_all$Avg_var_GRS_Diff[9] <- paste(diff_df_var_all$Avg_var_GRS_Diff[9], "(", t1, "-", t2, ")" )


## min and max for each sample size

diff_df_var_all$min[1] <- round(min(gba_s50_v$percent_diff_GRS), 4)
diff_df_var_all$max[1] <- round(max(gba_s50_v$percent_diff_GRS), 4)

diff_df_var_all$min[2] <- round(min(gba_s100_v$percent_diff_GRS), 4)
diff_df_var_all$max[2] <- round(max(gba_s100_v$percent_diff_GRS), 4)

diff_df_var_all$min[3] <- round(min(gba_s200_v$percent_diff_GRS), 4)
diff_df_var_all$max[3] <- round(max(gba_s200_v$percent_diff_GRS), 4)

diff_df_var_all$min[4] <- round(min(gba_s300_v$percent_diff_GRS), 4)
diff_df_var_all$max[4] <- round(max(gba_s300_v$percent_diff_GRS), 4)

diff_df_var_all$min[5] <- round(min(gba_s600_v$percent_diff_GRS), 4)
diff_df_var_all$max[5] <- round(max(gba_s600_v$percent_diff_GRS), 4)

diff_df_var_all$min[6] <- round(min(gba_s700_v$percent_diff_GRS), 4)
diff_df_var_all$max[6] <- round(max(gba_s700_v$percent_diff_GRS), 4)

diff_df_var_all$min[7] <- round(min(gba_s800_v$percent_diff_GRS), 4)
diff_df_var_all$max[7] <- round(max(gba_s800_v$percent_diff_GRS), 4)

diff_df_var_all$min[8] <- round(min(gba_s1000_v$percent_diff_GRS), 4)
diff_df_var_all$max[8] <- round(max(gba_s1000_v$percent_diff_GRS), 4)

diff_df_var_all$min[9] <- round(min(gba_s2000_v$percent_diff_GRS), 4)
diff_df_var_all$max[9] <- round(max(gba_s2000_v$percent_diff_GRS), 4)




## Fit multipolynomial regression line to data for information on relationship between sample size and percent difference in mean GRS
## LOESS for visual inspection of relationship and figures 



a <- c(mean(gba_s50_sig_v$percent_diff_GRS_Variance), mean(gba_s100_sig_v$percent_diff_GRS_Variance),
          mean(gba_s200_sig_v$percent_diff_GRS_Variance), mean(gba_s300_sig_v$percent_diff_GRS_Variance),
          mean(gba_s600_sig_v$percent_diff_GRS_Variance), mean(gba_s700_sig_v$percent_diff_GRS_Variance),
          mean(gba_s800_sig_v$percent_diff_GRS_Variance), mean(gba_s1000_sig_v$percent_diff_GRS_Variance),
          mean(gba_s2000_sig_v$percent_diff_GRS_Variance))


b <- diff_df_var$Sample_Size


testlm <- lm (a ~ poly(b,2))

summary(testlm)

ggplot(diff_df, aes(x = b, y = a)) + geom_point(size = 3) + stat_smooth(method = "loess", span = 1.17, formula =y ~ x,  col = "red")+
  labs(x = "Sample Size", y = "Average % Difference", title = "Average Percent Difference in Variance in GRS Between GBA Treatment Groups (Significant Iterations)")+
  theme_tufte()+ geom_rangeframe()+theme(text = element_text(size = 15)) 





## GBA iterations that have been balanced
## Either only 1 or 0 significant iterations so using all iterations


s50_gb <- read.table("GBA_iterations_balanced1000sampling50.tab")

s100_gb <- read.table("GBA_iterations_balanced1000sampling100.tab")

s200_gb <- read.table("GBA_iterations_balanced1000sampling200.tab")

s300_gb <- read.table("GBA_iterations_balanced1000sampling300.tab")

s600_gb <- read.table("GBA_iterations_balanced1000sampling600.tab")

s700_gb <- read.table("GBA_iterations_balanced1000sampling700.tab")

s800_gb <- read.table("GBA_iterations_balanced1000sampling800.tab")

s1000_gb <- read.table("GBA_iterations_balanced1000sampling1000.tab")

s2000_gb <- read.table("GBA_iterations_balanced1000sampling2000.tab")



df_list_bal <- list(s50_gb =s50_gb, s100_gb = s100_gb, s200_gb = s200_gb, s300_gb = s300_gb, 
                  s600_gb = s600_gb, s700_gb = s700_gb, s800_gb = s800_gb, s1000_gb = s1000_gb, s2000_gb = s2000_gb)



frames <- lapply(df_list_bal, FUN = change_fnc)

s50_gb <- as.data.frame(frames[1])
colnames(s50_gb) <- sub(".*\\.", "", colnames(s50_gb))

s100_gb <- as.data.frame(frames[2])
colnames(s100_gb) <- sub(".*\\.", "", colnames(s100_gb))

s200_gb <- as.data.frame(frames[3])
colnames(s200_gb) <- sub(".*\\.", "", colnames(s200_gb))

s300_gb <- as.data.frame(frames[4])
colnames(s300_gb) <- sub(".*\\.", "", colnames(s300_gb))

s600_gb <- as.data.frame(frames[5])
colnames(s600_gb) <- sub(".*\\.", "", colnames(s600_gb))

s700_gb <- as.data.frame(frames[6])
colnames(s700_gb) <- sub(".*\\.", "", colnames(s700_gb))

s800_gb <- as.data.frame(frames[7])
colnames(s800_gb) <- sub(".*\\.", "", colnames(s800_gb))

s1000_gb <- as.data.frame(frames[8])
colnames(s1000_gb) <- sub(".*\\.", "", colnames(s1000_gb))

s2000_gb <- as.data.frame(frames[9])
colnames(s2000_gb) <- sub(".*\\.", "", colnames(s2000_gb))



## create dataframe of average percent difference in variance in GRS for each sample size

diff_df_bal <- data.frame(Sample_Size = numeric(9), Avg_mean_GRS_Diff = numeric(9), min= numeric(9), max = numeric(9))
diff_df_bal$Sample_Size <- c(50, 100, 200, 300, 600, 700, 800, 1000, 2000)
diff_df_bal$Avg_mean_GRS_Diff <- c(mean(s50_gb$percent_diff_GRS), mean(s100_gb$percent_diff_GRS),
                                  mean(s200_gb$percent_diff_GRS), mean(s300_gb$percent_diff_GRS),
                                  mean(s600_gb$percent_diff_GRS), mean(s700_gb$percent_diff_GRS),
                                  mean(s800_gb$percent_diff_GRS), mean(s1000_gb$percent_diff_GRS),
                                  mean(s2000_gb$percent_diff_GRS))



diff_df_bal$Avg_mean_GRS_Diff <- round(diff_df_bal$Avg_mean_GRS_Diff, 2)

t <- t.test(s50_gb$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_bal$Avg_mean_GRS_Diff[1] <- paste(diff_df_bal$Avg_mean_GRS_Diff[1], "(", t1, "-", t2, ")" )


t <- t.test(s100_gb$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_bal$Avg_mean_GRS_Diff[2] <- paste(diff_df_bal$Avg_mean_GRS_Diff[2], "(", t1, "-", t2, ")" )


t <- t.test(s200_gb$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_bal$Avg_mean_GRS_Diff[3] <- paste(diff_df_bal$Avg_mean_GRS_Diff[3], "(", t1, "-", t2, ")" )


t <- t.test(s300_gb$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_bal$Avg_mean_GRS_Diff[4] <- paste(diff_df_bal$Avg_mean_GRS_Diff[4], "(", t1, "-", t2, ")" )


t <- t.test(s600_gb$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_bal$Avg_mean_GRS_Diff[5] <- paste(diff_df_bal$Avg_mean_GRS_Diff[5], "(", t1, "-", t2, ")" )


t <- t.test(s700_gb$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_bal$Avg_mean_GRS_Diff[6] <- paste(diff_df_bal$Avg_mean_GRS_Diff[6], "(", t1, "-", t2, ")" )


t <- t.test(s800_gb$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_bal$Avg_mean_GRS_Diff[7] <- paste(diff_df_bal$Avg_mean_GRS_Diff[7], "(", t1, "-", t2, ")" )


t <- t.test(s1000_gb$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_bal$Avg_mean_GRS_Diff[8] <- paste(diff_df_bal$Avg_mean_GRS_Diff[8], "(", t1, "-", t2, ")" )


t <- t.test(s2000_gb$percent_diff_GRS)
t1 <- round(t$conf.int[1], 2)
t2 <- round(t$conf.int[2], 2)
diff_df_bal$Avg_mean_GRS_Diff[9] <- paste(diff_df_bal$Avg_mean_GRS_Diff[9], "(", t1, "-", t2, ")" )


# min and max

diff_df_bal$min[1] <- round(min(s50_gb$percent_diff_GRS), 4)
diff_df_bal$max[1] <- round(max(s50_gb$percent_diff_GRS), 4)

diff_df_bal$min[2] <- round(min(s100_gb$percent_diff_GRS), 4)
diff_df_bal$max[2] <- round(max(s100_gb$percent_diff_GRS), 4)

diff_df_bal$min[3] <- round(min(s200_gb$percent_diff_GRS), 4)
diff_df_bal$max[3] <- round(max(s200_gb$percent_diff_GRS), 4)

diff_df_bal$min[4] <- round(min(s300_gb$percent_diff_GRS), 4)
diff_df_bal$max[4] <- round(max(s300_gb$percent_diff_GRS), 4)

diff_df_bal$min[5] <- round(min(s600_gb$percent_diff_GRS), 4)
diff_df_bal$max[5] <- round(max(s600_gb$percent_diff_GRS), 4)

diff_df_bal$min[6] <- round(min(s700_gb$percent_diff_GRS), 4)
diff_df_bal$max[6] <- round(max(s700_gb$percent_diff_GRS), 4)

diff_df_bal$min[7] <- round(min(s800_gb$percent_diff_GRS), 4)
diff_df_bal$max[7] <- round(max(s800_gb$percent_diff_GRS), 4)

diff_df_bal$min[8] <- round(min(s1000_gb$percent_diff_GRS), 4)
diff_df_bal$max[8] <- round(max(s1000_gb$percent_diff_GRS), 4)

diff_df_bal$min[9] <- round(min(s2000_gb$percent_diff_GRS), 4)
diff_df_bal$max[9] <- round(max(s2000_gb$percent_diff_GRS), 4)
