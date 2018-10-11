##### Part of a series of analyses associated with the manuscript "Genetic variability and potential effects on clinical
##### trial outcomes: perspectives in Parkinson's disease"
##### Hampton Leonard 10/09/18


##### Create virtual cohort for interaction of rs9298897 and rs17710829


library(dplyr)
library(car)

set.seed(1234567)


cohort1 <- data.frame("ID" = numeric(5000), "rs9298897" = numeric(5000), "rs17710829"=numeric(5000))

cohort1$ID <- 1:5000


## Assign variant frequencies according to ExAC Release 1 database and HWE

## snp rs9298897 
## G=0.3393	A=0.6607


p <- 0.3393
q <- 0.6607

GG_freq <- p^2

GA_freq <- 2*p*q

AA_freq <- q^2


GG_pop <- round(GG_freq *5000)

GA_pop <- round(GA_freq * 5000)

AA_pop <- round(AA_freq * 5000)-1

samp1 <- sample_n(cohort1, GG_pop)
samp1$rs9298897 <- 2

samp2 <- sample_n(cohort1[cohort1$ID %in% samp1$ID == FALSE,], GA_pop)
samp2$rs9298897 <- samp2$rs9298897 <- 1

samp3 <- sample_n(cohort1[cohort1$ID %in% samp1$ID == FALSE & cohort1$ID %in% samp2$ID == FALSE,], AA_pop)
samp3$rs9298897 <- 0


samp4 <- rbind(samp1, samp2, samp3)

cohort1$rs9298897 <- samp4$rs9298897


## rs17710829
## T = 0.9427 C = 0.0573


p <- 0.9427
q <- 0.0573

TT_freq <- p^2

TC_freq <- 2*p*q

CC_freq <- q^2


TT_pop <- round(TT_freq *5000)

TC_pop <- round(TC_freq * 5000)+1

CC_pop <- round(CC_freq * 5000)

samp5 <- sample_n(cohort1, TT_pop)
samp5$rs17710829 <- 0

samp6 <- sample_n(cohort1[cohort1$ID %in% samp5$ID == FALSE,], TC_pop)
samp6$rs17710829 <- 1

samp7 <- sample_n(cohort1[cohort1$ID %in% samp5$ID == FALSE & cohort1$ID %in% samp6$ID == FALSE,], CC_pop)
samp7$rs17710829 <- 2

cohortfinal <- rbind(samp5, samp6, samp7)




set.seed(1234567)

## Assign arbitrary range of moderate baseline UPDRS scores to all patients

cohortfinal$updrs_baseline <- sample(15:25, 5000, replace = T)

set.seed(1234567)
nIterations <- 1000 #test at 1000 iterations
# nIterations <- 10 #test at 1000 iterations
for(z in c(50, 100, 200, 300, 600, 700, 800))
#for(z in c(50))
{
  
  ## Dataframe: treated_mean_change = mean for change in updrs for treatement arm without snp effect; 
  ## untreated_mean_change = mean for change in updrs for placbeo arm without snp effect; treated_mean_change_wsnp = mean 
  ## for change in updrs for treatment arm with snp effect; untreated_mean_change_wsnp = mean for change in placebp arm with snp effect; z_score = z score for difference
  ## without snp effect; z_score_wsnp = z score for difference with snp effect; sum_treated = number of patients with the interaction in the treatment arm; sum_untreated =
  ## number of patients with the interaction in the placebo arm
  
  snp_scores <- data.frame(treated_mean_change = numeric(nIterations), untreated_mean_change = numeric(nIterations), 
                           treated_mean_change_wsnp = numeric(nIterations), untreated_mean_change_wsnp = numeric(nIterations),
                           z_score = numeric(nIterations), z_score_wsnp = numeric(nIterations), sum_treated = numeric(nIterations),
                           sum_untreated = numeric(nIterations))
  
  
  
  for(y in 1:nIterations){
    
    temp <- sample_n(cohortfinal, z)
    
    ## Randomly assign to treatment or placebo
    
    temp$pheno <- rep(0:1, length.out = length(temp$updrs_baseline))
    
    ## Randomly assign average yearly progression to all patients
    
    temp$updrs_final <- temp$updrs_baseline + sample(1:2, length(temp$updrs_baseline), replace = T)
    
    
    print(paste(" at y iterations = ", y, " for sampling N =", z, sep = ""))
    
    
    ## Split into treated and untreated  
    
    untreated <- temp[temp$pheno == 0,]
    treated <- temp[temp$pheno == 1,]
    
    
    ## Treatment effect: -0.5 updrs points per year; comment out if looking at false positives  
    
    #treated$updrs_final <- treated$updrs_final - 0.5
    
    ## Calculate change in updrs
    
    treated$change <- treated$updrs_final - treated$updrs_baseline
    untreated$change <- untreated$updrs_final - untreated$updrs_baseline
    
    t_change_mean <- mean(treated$change)
    ut_change_mean <- mean(untreated$change)
    
    
    
    snp_scores[y,1] <- mean(treated$change)
    snp_scores[y,2] <- mean(untreated$change)
    snp_scores[y,5] <- (t_change_mean - ut_change_mean)/sqrt(((sd(untreated$change)^2))/nrow(untreated) + ((sd(treated$change)^2))/nrow(treated))
    snp_scores[y,7] <- (nrow(treated[treated$rs9298897 >0 & treated$rs17710829 >0,]))
    snp_scores[y,8] <- (nrow(untreated[untreated$rs9298897 >0 & untreated$rs17710829 >0,]))
    
    
    
    ## Add SNP interaction effect 
    
    treated$updrs_final_wv <- treated$updrs_final
    untreated$updrs_final_wv <- untreated$updrs_final
    
    for(i in 1:nrow(treated)){
      if(treated$rs9298897[i] >0 & treated$rs17710829[i] >0){
        treated$updrs_final_wv[i] <- treated$updrs_final[i] + 2.374
      }
    }
    
    
    for(i in 1:nrow(untreated)){
      if(untreated$rs9298897[i] >0 & untreated$rs17710829[i] >0){
        untreated$updrs_final_wv[i] <- untreated$updrs_final[i] + 2.374
      }
    }
    
    
    ## Calculate change in updrs
    
    treated$change2 <- treated$updrs_final_wv - treated$updrs_baseline
    untreated$change2 <- untreated$updrs_final_wv - untreated$updrs_baseline
    
    t_change_mean2 <- mean(treated$change2)
    ut_change_mean2 <- mean(untreated$change2)
    
    
    
    snp_scores[y,3] <- mean(treated$change2)
    snp_scores[y,4] <- mean(untreated$change2)
    snp_scores[y,6] <- (t_change_mean2 - ut_change_mean2)/sqrt(((sd(untreated$change2)^2))/nrow(untreated) + ((sd(treated$change2)^2))/nrow(treated))
    
    
    
    
    ## Add a second year; for only analysis of the first year, comment below out until line 259
    
    # ## Randomly assign additional yearly progression in updrs
    # 
    # treated$updrs_final <- treated$updrs_final + sample(1:2, length(treated$updrs_final), replace = T)
    # untreated$updrs_final <- untreated$updrs_final + sample(1:2, length(untreated$updrs_baseline), replace = T)
    # 
    # 
    # ## Treatment effect: subtract an additional 0.5 updrs points for an additional year
    # 
    # treated$updrs_final <- treated$updrs_final - 0.5
    # 
    # 
    # ## Calculate change in updrs
    # 
    # treated$change <- treated$updrs_final - treated$updrs_baseline
    # untreated$change <- untreated$updrs_final - untreated$updrs_baseline
    # 
    # t_change_mean <- mean(treated$change)
    # ut_change_mean <- mean(untreated$change)
    # 
    # snp_scores[y,1] <- mean(treated$change)
    # snp_scores[y,2] <- mean(untreated$change)
    # snp_scores[y,5] <- (t_change_mean - ut_change_mean)/sqrt(((sd(untreated$change)^2))/nrow(untreated) + ((sd(treated$change)^2))/nrow(treated))
    # 
    # 
    # 
    # 
    # ## Year 2 adding SNP interaction effect
    # 
    # 
    # treated$updrs_final_wv <- treated$updrs_final
    # untreated$updrs_final_wv <- untreated$updrs_final
    # 
    # 
    # ## Adding 4.374 to account for 2 years of the SNP interaction effect
    # 
    # for(i in 1:nrow(treated)){
    #   if(treated$rs9298897[i] >0 & treated$rs17710829[i] >0){
    #     treated$updrs_final_wv[i] <- treated$updrs_final[i] + 4.374
    #   }
    # }
    # 
    # 
    # for(i in 1:nrow(untreated)){
    #   if(untreated$rs9298897[i] >0 & untreated$rs17710829[i] >0){
    #     untreated$updrs_final_wv[i] <- untreated$updrs_final[i] + 4.374
    #   }
    # }
    # 
    # 
    # ## Calculate change in updrs
    # 
    # treated$change2 <- treated$updrs_final_wv - treated$updrs_baseline
    # untreated$change2 <- untreated$updrs_final_wv - untreated$updrs_baseline
    # 
    # t_change_mean2 <- mean(treated$change2)
    # ut_change_mean2 <- mean(untreated$change2)
    # 
    # 
    # snp_scores[y,3] <- mean(treated$change2)
    # snp_scores[y,4] <- mean(untreated$change2)
    # snp_scores[y,6] <- (t_change_mean2 - ut_change_mean2)/sqrt(((sd(untreated$change2)^2))/nrow(untreated) + ((sd(treated$change2)^2))/nrow(treated))
    
    
  }
  
  write.table(snp_scores, paste("snp_interaction_year1_fp",nIterations,"sampling",z,".tab",sep = ""), quote = F, sep = "\t")  
}


## Replace file names for trial year 2 files as necessary
## False negatives

s50 <- read.table("snp_interaction_year11000sampling50.tab", header = TRUE)

## How many are false negatives without the SNP effect

s50_s <- s50[(s50$z_score < 1.96 & s50$z_score > -1.96),]

## How many are false negatives when SNP effect is added

s50_v <- s50[(s50$z_score_wsnp < 1.96 & s50$z_score_wsnp > -1.96),]

## How many false negatives were originally true positives before SNP effect was added? (How many FN due to SNP effect alone)

s50_a <- s50[s50$z_score > 1.96 | s50$z_score  < -1.96,] 

s50_b <- s50_a[s50_a$z_score_wsnp < 1.96 & s50_a$z_score_wsnp > -1.96,]

nrow(s50_b)/nrow(s50_v)


## Repeat for all sample sizes

s100 <- read.table("snp_interaction_year11000sampling100.tab", header = TRUE)


s100_s <- s100[(s100$z_score < 1.96 & s100$z_score > -1.96),]

s100_v <- s100[(s100$z_score_wsnp < 1.96 & s100$z_score_wsnp > -1.96),]

s100_a <- s100[s100$z_score > 1.96 | s100$z_score  < -1.96,] 

s100_b <- s100_a[s100_a$z_score_wsnp < 1.96 & s100_a$z_score_wsnp > -1.96,]

nrow(s100_b)/nrow(s100_v)



s200 <- read.table("snp_interaction_year11000sampling200.tab", header = TRUE)


s200_s <- s200[(s200$z_score < 1.96 & s200$z_score > -1.96),]

s200_v <- s200[(s200$z_score_wsnp < 1.96 & s200$z_score_wsnp > -1.96),]

s200_a <- s200[s200$z_score > 1.96 | s200$z_score  < -1.96,] 

s200_b <- s200_a[s200_a$z_score_wsnp < 1.96 & s200_a$z_score_wsnp > -1.96,]

nrow(s200_b)/nrow(s200_v)



s300 <- read.table("snp_interaction_year11000sampling300.tab", header = TRUE)


s300_s <- s300[(s300$z_score < 1.96 & s300$z_score > -1.96),]

s300_v <- s300[(s300$z_score_wsnp < 1.96 & s300$z_score_wsnp > -1.96),]

s300_a <- s300[s300$z_score > 1.96 | s300$z_score  < -1.96,] 

s300_b <- s300_a[s300_a$z_score_wsnp < 1.96 & s300_a$z_score_wsnp > -1.96,]

nrow(s300_b)/nrow(s300_v)



s600 <- read.table("snp_interaction_year11000sampling600.tab", header = TRUE)


s600_s <- s600[(s600$z_score < 1.96 & s600$z_score > -1.96),]

s600_v <- s600[(s600$z_score_wsnp < 1.96 & s600$z_score_wsnp > -1.96),]

s600_a <- s600[s600$z_score > 1.96 | s600$z_score  < -1.96,] 

s600_b <- s600_a[s600_a$z_score_wsnp < 1.96 & s600_a$z_score_wsnp > -1.96,]

nrow(s600_b)/nrow(s600_v)



s700 <- read.table("snp_interaction_year11000sampling700.tab", header = TRUE)


s700_s <- s700[(s700$z_score < 1.96 & s700$z_score > -1.96),]

s700_v <- s700[(s700$z_score_wsnp < 1.96 & s700$z_score_wsnp > -1.96),]

s700_a <- s700[s700$z_score > 1.96 | s700$z_score  < -1.96,] 

s700_b <- s700_a[s700_a$z_score_wsnp < 1.96 & s700_a$z_score_wsnp > -1.96,]

nrow(s700_b)/nrow(s700_v)



s800 <- read.table("snp_interaction_year11000sampling800.tab", header = TRUE)


s800_s <- s800[(s800$z_score < 1.96 & s800$z_score > -1.96),]

s800_v <- s800[(s800$z_score_wsnp < 1.96 & s800$z_score_wsnp > -1.96),]

s800_a <- s800[s800$z_score > 1.96 | s800$z_score  < -1.96,] 

s800_b <- s800_a[s800_a$z_score_wsnp < 1.96 & s800_a$z_score_wsnp > -1.96,]

nrow(s800_b)/nrow(s800_v)




## False positives


s50_fp <- read.table("snp_interaction_year1_fp1000sampling50.tab", header = T)

## How many false positives without SNP effect?

s50_fp_s <- s50_fp[s50_fp$z_score > 1.96 | s50_fp$z_score < -1.96,]

## How many false positives with SNP effect?

s50_fp_v <- s50_fp[s50_fp$z_score_wsnp > 1.96 | s50_fp$z_score_wsnp < -1.96,]

## How many false positives were true negatives before the SNP effect was added?

s50_fp_a <- s50_fp[s50_fp$z_score < 1.96 & s50_fp$z_score  > -1.96,] 

s50_fp_b <- s50_fp_a[s50_fp_a$z_score_wsnp > 1.96 | s50_fp_a$z_score_wsnp < -1.96,]

nrow(s50_fp_b)/nrow(s50_fp_v)


## Repeat for all sample sizes


s100_fp <- read.table("snp_interaction_year1_fp1000sampling100.tab", header = T)


s100_fp_s <- s100_fp[s100_fp$z_score > 1.96 | s100_fp$z_score < -1.96,]

s100_fp_v <- s100_fp[s100_fp$z_score_wsnp > 1.96 | s100_fp$z_score_wsnp < -1.96,]

s100_fp_a <- s100_fp[s100_fp$z_score < 1.96 & s100_fp$z_score  > -1.96,] 

s100_fp_b <- s100_fp_a[s100_fp_a$z_score_wsnp > 1.96 | s100_fp_a$z_score_wsnp < -1.96,]

nrow(s100_fp_b)/nrow(s100_fp_v)



s200_fp <- read.table("snp_interaction_year1_fp1000sampling200.tab", header = T)


s200_fp_s <- s200_fp[s200_fp$z_score > 1.96 | s200_fp$z_score < -1.96,]

s200_fp_v <- s200_fp[s200_fp$z_score_wsnp > 1.96 | s200_fp$z_score_wsnp < -1.96,]

s200_fp_a <- s200_fp[s200_fp$z_score < 1.96 & s200_fp$z_score  > -1.96,] 

s200_fp_b <- s200_fp_a[s200_fp_a$z_score_wsnp > 1.96 | s200_fp_a$z_score_wsnp < -1.96,]

nrow(s200_fp_b)/nrow(s200_fp_v)



s300_fp <- read.table("snp_interaction_year1_fp1000sampling300.tab", header = T)


s300_fp_s <- s300_fp[s300_fp$z_score > 1.96 | s300_fp$z_score < -1.96,]

s300_fp_v <- s300_fp[s300_fp$z_score_wsnp > 1.96 | s300_fp$z_score_wsnp < -1.96,]

s300_fp_a <- s300_fp[s300_fp$z_score < 1.96 & s300_fp$z_score  > -1.96,] 

s300_fp_b <- s300_fp_a[s300_fp_a$z_score_wsnp > 1.96 | s300_fp_a$z_score_wsnp < -1.96,]

nrow(s300_fp_b)/nrow(s300_fp_v)



s600_fp <- read.table("snp_interaction_year1_fp1000sampling600.tab", header = T)


s600_fp_s <- s600_fp[s600_fp$z_score > 1.96 | s600_fp$z_score < -1.96,]

s600_fp_v <- s600_fp[s600_fp$z_score_wsnp > 1.96 | s600_fp$z_score_wsnp < -1.96,]

s600_fp_a <- s600_fp[s600_fp$z_score < 1.96 & s600_fp$z_score  > -1.96,] 

s600_fp_b <- s600_fp_a[s600_fp_a$z_score_wsnp > 1.96 | s600_fp_a$z_score_wsnp < -1.96,]

nrow(s600_fp_b)/nrow(s600_fp_v)



s700_fp <- read.table("snp_interaction_year1_fp1000sampling700.tab", header = T)


s700_fp_s <- s700_fp[s700_fp$z_score > 1.96 | s700_fp$z_score < -1.96,]

s700_fp_v <- s700_fp[s700_fp$z_score_wsnp > 1.96 | s700_fp$z_score_wsnp < -1.96,]

s700_fp_a <- s700_fp[s700_fp$z_score < 1.96 & s700_fp$z_score  > -1.96,] 

s700_fp_b <- s700_fp_a[s700_fp_a$z_score_wsnp > 1.96 | s700_fp_a$z_score_wsnp < -1.96,]

nrow(s700_fp_b)/nrow(s700_fp_v)



s800_fp <- read.table("snp_interaction_year1_fp1000sampling800.tab", header = T)


s800_fp_s <- s800_fp[s800_fp$z_score > 1.96 | s800_fp$z_score < -1.96,]

s800_fp_v <- s800_fp[s800_fp$z_score_wsnp > 1.96 | s800_fp$z_score_wsnp < -1.96,]

s800_fp_a <- s800_fp[s800_fp$z_score < 1.96 & s800_fp$z_score  > -1.96,] 

s800_fp_b <- s800_fp_a[s800_fp_a$z_score_wsnp > 1.96 | s800_fp_a$z_score_wsnp < -1.96,]

nrow(s800_fp_b)/nrow(s800_fp_v)

