
#FinnDiane prospective survival analysis codes


########################## Step  1   ###########################
#Exctract M values from raw methylation data. 
#The resulted file "top_CpG_M_values.txt" and the file with clinical data (FD_EWAS_pheno.txt) 
#contain sensitive individual-level data and are therefore not publicly available. 
#Regarding further information on data availability we refer to the manuscript ("Data availability section",Smyth and Dahlstrom et al 2022)


## Read in data

  #M values 
  m <- read.table("top_CpG_M_values.txt", header=T, sep="\t", as.is=T, na.strings="NA") 
  #phenotype table 
  pheno <- read.table("FD_EWAS_pheno.txt", header=T, sep="\t", as.is=T, na.strings="NA") 

## Modify the M values file m
  m_t <- transpose(m)
  colnames(m_t) <- rownames(m)

#Add individual ids
  id <- colnames(m)
  m_t <- cbind(id, m_t)

#rename transposed m_t to m 
  m <- m_t 

#Combine pheno table and M values by “id”
  data <- left_join(pheno, m, by=c("id")) 

#select only individuals with macroalbuminuria
  data <- subset(data, CASE_CONT==c("case")) 

#Remove those with no progression information 
  data <- subset(data, !is.na(get_esrd) & !is.na(time_to_esrd_from_visit)) 
  nrow(data) #397

## Run Cox 
  
  #List CpGs (in our data CpG M values start at column 90)
  cgs <- colnames(data[90:ncol(data)]) 

  #create a table to collect Cox results from 
  collect <- data.frame(matrix(ncol=19, nrow=length(cgs)))
  colnames(collect) <- 
  c("CpG_name", "M1_variables", "M2_variables", "M2_variables", "M4_variables", 
  "M1_CpG_P", "M1_CpG_HR", "M1_CpG_HR95_low", "M1_CpG_HR95_up", 
  "M2_CpG_P", "M2_CpG_HR", "M2_CpG_HR95_low", "M2_CpG_HR95_up", 
  "M4_CpG_P", "M4_CpG_HR", "M4_CpG_HR95_low", "M4_CpG_HR95_up",
  "M3_C_index", "M4_C_index")
   
  collect$CG_name <- cgs 
 
##loop to analyse CpGs with model 1 and 2
  
for (i in 1:length(cgs)) {
CG <- cgs[i]
 
  #MODEL1
  model <- c("M1")
  covariates <- "CG, sex, age, 6 WCC"
  cox <- coxph(Surv(data$time_to_esrd_from_visit, data$get_esrd) ~ dataS[,CG] + SEX + AGE + CD8T + CD4T + Bcell + NK + Gran + Mono, data=data)  
  
  #put to collect
  collect[paste0(model, "_", c("variables"))][i,] <- covariates
  collect[paste0(model, "_", c("CpG_P"))][i,] <- coef(summary(cox))[1,5]
  collect[paste0(model, "_", c("CpG_HR"))][i,] <- round(coef(summary(cox))[1,2],3)
  collect[paste0(model, "_", c("CpG_HR95_low"))][i,]  <- round(exp(confint(cox))[1,1], 3)
  collect[paste0(model, "_", c("CpG_HR95_up"))][i,]   <- round(exp(confint(cox))[1,2], 3)  
  
  #MODEL 2
  model <- c("M2")
  covariates <- "CG, eGFR, male sex, age, ageonset, SBP, HbA1c, trigly, retinal laser, smoking, 6 WCC"
  cox <- coxph(Surv(data$time_to_esrd_from_visit, data$get_esrd) ~ dataS[,CG] + eGFR + SEX + AGE + CD8T + CD4T + Bcell + NK + Gran + Mono, data)  
 
  #put to collect
  collect[paste0(model, "_", c("variables"))][i,] <- covariates
  collect[paste0(model, "_", c("CpG_P"))][i,] <- coef(summary(cox))[1,5]
  collect[paste0(model, "_", c("CpG_HR"))][i,] <- round(coef(summary(cox))[1,2],3)
  collect[paste0(model, "_", c("CpG_HR95_low"))][i,]  <- round(exp(confint(cox))[1,1], 3)
  collect[paste0(model, "_", c("CpG_HR95_up"))][i,]   <- round(exp(confint(cox))[1,2], 3)  
  
}
 
##Extra filtering step to compare Model3 and Model4 C-indexes
  data < - subset(data, !is.na(data$eGFR_idms) & !is.na(data$SEX) & !is.na(data$AGE) & !is.na(data$AGEONSET) & !is.na( data$SBP)  & !is.na(data$hba1cprc) & !is.na(data$TG) & !is.na(data$LASER) & !is.na(data$CURRENTSMOKER))
  nrow(data) #388
 
##loop to analyse all CpGs with models 3 and 4 
  
for (i in 1:length(cgs)) {
CG <- cgs[i]
 
  
    #MODEL3
    model <- c("M3")
    covariates <- "sex, age, ageonset, SBP, HbA1c, triglycerides, retinal photocoagulation, smoking, 6 WCC"
    cox <- coxph(Surv(data$time_to_esrd_from_visit, data$get_esrd) ~ SEX + AGE + AGEONSET + SBP + hba1c + TG + LASER + CURRENTSMOKER + Bcell + Mono + CD4T + CD8T + Gran + NK, data = data)
   
    #put to collect
    collect[paste0(model, "_", c("variables"))][i,] <- covariates
    collect[paste0(model, "_", c("C_index"))][i,] <- concordance(cox)[1]
   
    #MODEL4 = MODEL3 + CpG
    model <- c("M4")
    covariates <- "CG, sex, age, ageonset, SBP, HbA1c, triglycerides, retinal photocoagulation, smoking, 6 WCC"
    cox <- coxph(Surv(data$time_to_esrd_from_visit, data$get_esrd) ~ dataS[,CG] + SEX + AGE + AGEONSET + SBP + hba1c + TG + LASER + CURRENTSMOKER + Bcell + Mono + CD4T + CD8T + Gran + NK, data = data)
    
    #put to collect
    collect[paste0(model, "_", c("covariates"))][i,] <- covariates
    collect[paste0(model, "_", c("C_index"))][i,] <- concordance(cox)[1] 
    collect[paste0(model, "_", c("CpG_P"))][i,] <- coef(summary(cox))[1,5]
    collect[paste0(model, "_", c("CpG_HR"))][i,] <- round(coef(summary(cox))[1,2],3)
    collect[paste0(model, "_", c("CpG_HR95_low"))][i,]  <- round(exp(confint(cox))[1,1], 3)
    collect[paste0(model, "_", c("CpG_HR95_up"))][i,]   <- round(exp(confint(cox))[1,2], 3)
  }
  
  
write.table(collect, file="cox_results_all_modes.txt", quote = FALSE, sep="\t", row.names=F, col.names=T) 

## Compare C-indexes of model3 and model4
# Instructions here https://cran.r-project.org/web/packages/survival/vignettes/concordance.pdf
# Calculate for all CpGs, in loop or separately

CG <- "this is the CpG studied"
fit3 <- coxph(Surv(data$time_to_esrd_from_visit, data$get_esrd) ~ SEX + AGE + AGEONSET + SBP + hba1c + TG + LASER + CURRENTSMOKER + Bcell + Mono + CD4T + CD8T + Gran + NK, data = data) 
fit4 <- coxph(Surv(data$time_to_esrd_from_visit, data$get_esrd) ~ dataS[,CG] + SEX + AGE + AGEONSET + SBP + hba1c + TG + LASER + CURRENTSMOKER + Bcell + Mono + CD4T + CD8T + Gran + NK, data = data)
ctest <- concordance(fit1, fit2)
contr <- c(-1, 1)
dtest <- contr %*% coef(ctest)
dtest <- contr %*% coef(ctest)
dvar <- contr %*% vcov(ctest) %*% contr
res <- c(contrast=dtest, sd=sqrt(dvar), z=dtest/sqrt(dvar))

#P value
p <- 2*pnorm(q=as.numeric(res[3]), lower.tail=FALSE)
p


