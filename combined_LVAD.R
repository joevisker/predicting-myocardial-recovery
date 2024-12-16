library(tidyverse)
library(data.table)
library(randomForest)
library(randomForestSRC)
library(glmnet)
library(pdp)
library(party)
library(permimp)
library(forcats)
library(caret)
library(pdp)
library(pROC)
library(e1071)
library(gbm)
library(performanceEstimation)
library(furrr)
library(nnet)
library(mice)
library(caret)
library(energy)
#setwd("XXX")

# Load and combine datasets
clin_dat_der <- read_csv("Clin_dat_Der.csv") %>% filter(!is.na(study_id)) %>% mutate(study_id=as.character(study_id))
clin_dat_val <- read_csv("Clin_dat_Val.csv") %>% filter(!is.na(study_id)) %>% select(colnames(clin_dat_der))


Deriv_idx <- read_csv("Derivation Cohort_KEY_INDEX.csv")
Valid_idx <- read_csv("Validation Cohort_KEY_INDEX.csv")

rna_deriv <- read_csv("Derivation Cohort Repository Analysis_R-Log Values (n=130).csv")
rna_valid <- read_csv("Validation Cohort Repository Analysis_R-Log Values (n=134).csv")


# Changing RNA data from wide to long format
newcols=rna_deriv$Gene
df=data.frame(t(rna_deriv %>% select(-`Gene`,-`Associated Gene Name`)))
colnames(df)=newcols
df_rna=data.frame(setDT(df, keep.rownames = TRUE)[])
df_rna[1:5,1:5]
###
newcols=rna_valid$ensembl_id
df_valid=data.frame(t(rna_valid %>% select(-gene_name,-ensembl_id)))
colnames(df_valid)=newcols
df_rna_valid=data.frame(setDT(df_valid, keep.rownames = TRUE)[])
df_rna_valid[1:5,1:5]
#They must share gene names 
df_rna = df_rna %>% select(rn,one_of(colnames(df_rna)[which(colnames(df_rna) %in% colnames(df_rna_valid))]))
df_rna_valid = df_rna_valid %>% select(rn,one_of(colnames(df_rna_valid)[which(colnames(df_rna_valid) %in% colnames(df_rna))]))
df_rna_valid = df_rna_valid %>% select(one_of(colnames(df_rna)))

# Filtering to only pre-LVAD info and certain data

Deriv_idx = Deriv_idx %>% rename(RNA_ID=`RNA Seq. ID`) %>% mutate(study_id=as.character(`LVAD ID`)) %>%  filter(`Pre/Post`=="PRE")#

Valid_idx=Valid_idx %>% rename(RNA_ID=`RNAseq ID`) %>% mutate(study_id=case_when((startsWith(`Patient Name`,"UU") & substr(`Patient Name`,3,3)==0) ~ substr(`Patient Name`,4,nchar(`Patient Name`)),
                                                  startsWith(`Patient Name`,"UU") ~ substr(`Patient Name`,3,nchar(`Patient Name`)),
                                                  TRUE~`Patient Name`)) 

IDlist=c(Deriv_idx$study_id[which(Deriv_idx$study_id %in% clin_dat_der$study_id)],Valid_idx$study_id[which(Valid_idx$study_id %in% clin_dat_val$study_id)])
IDlist=IDlist[-which(!(rbind(Deriv_idx %>% select(RNA_ID,study_id),Valid_idx %>% select(RNA_ID,study_id))$RNA_ID[match(IDlist,rbind(Deriv_idx %>% select(RNA_ID,study_id),Valid_idx %>% select(RNA_ID,study_id))$study_id)] %in% c(df_rna$rn,df_rna_valid$rn)))]
IDlist=IDlist[-82]
clin_dat=rbind(clin_dat_der,clin_dat_val)


IDresp=clin_dat %>% select(hospital,study_id,resp=`Response Status (Non-Responders & Partial Responders=0, Responders=1, to be used for the analysis)`) %>% filter(study_id %in% IDlist) 
table(IDresp$resp)

clin_dat = clin_dat %>% mutate(resp=`Response Status (Non-Responders & Partial Responders=0, Responders=1, to be used for the analysis)`)

clin_dat = clin_dat %>% mutate(ACEi_ARB=as.numeric(`ACE inhibitor pre-LVAD (0=no, 1=yes)`|`ARB pre-LVAD (0=no, 1=yes)`)) %>% select(-`ACE inhibitor pre-LVAD (0=no, 1=yes)`,
                                                                                                                                                -`ARB pre-LVAD (0=no, 1=yes)`)


RNA_study_ids=rbind(Deriv_idx %>% select(RNA_ID,study_id),Valid_idx %>% select(RNA_ID,study_id)) %>% filter(RNA_ID != "12422X1")
rna_dat = rbind(df_rna,df_rna_valid)[-c(268,269)] %>% rename(RNA_ID=rn)

rna_dat = clin_dat %>% left_join(RNA_study_ids,by="study_id")  %>% left_join(rna_dat,by="RNA_ID") %>% select(study_id,resp,starts_with("ENSG"))

idx_rem=which(apply(rna_dat,1,function(x) sum(is.na(x)))==22371)

rna_dat=rna_dat[-idx_rem,]
clin_dat=clin_dat[-idx_rem,]

# Selecting predictor and response columns only 

clin_dat = cbind(clin_dat %>% select(-contains('Post'),-contains('change'),-`Response Status (U-NOVA stages)`,-
                                         `Response Status (Non-Responders & Partial Responders=0, Responders=1, to be used for the analysis)`),
                       clin_dat %>% select(`Post-LVAD LVEF (max)`,`Post-LVAD LVEDD at max LVEF`))


clin_dat=clin_dat %>% select(-`smoking (0=no, 1=yes)`,-`substance_abuse (0=no, 1=yes)`,-`alcohol (0=no, 1=yes)`,
  -`race/ethnicity`,-height,-weight,-`heart failure etiology (expanded)`,-`NYHA class`,
  -`LVAD indication (1=BTT, 2=DT, 3=BTD, 4=BTR)`,#-`Device therapy pre-LVAD (0=None, 1=CRT-D, 2=ICD)`,
  -`Aspirin pre-LVAD (0=no, 1=yes)`,-`Clopidogrel pre-LVAD (0=no, 1=yes)`,-`Anticoagulation pre-LVAD (0=no, 1=yes)`,-`Antiarrhythmic pre-LVAD (0=no, 1=yes)`,-`systolic blood pressure pre-LVAD`,-
    `diastolic blood pressure pre-LVAD`,-`mean blood pressure pre-LVAD`,-`heart rate pre-LVAD`,-`RV systolic pressure pre-LVAD`,-`RV diastolic pressure pre-LVAD`,-`RV mean pressure pre-LVAD`,-
    `pvr pre-LVAD`,-
    `svr pre-LVAD`,-
    `papi pre-LVAD`,-
    `right atrial pressure/pcwp pre-LVAD`,-
    `hematoctrit pre-LVAD`,-
    `glucose pre-LVAD`,-
    `uric acid pre-LVAD`,-
    `alp pre-LVAD`,-
    `direct bilirubin pre-LVAD`,-
    `inr pre-LVAD`            ,-
    `hba1c pre-LVAD`          ,-
    `prealbumin pre-LVAD`     ,-
    `ldh pre-LVAD`            ,-
    heart_rate_echo           ,-
    rvot_echo                 ,-
    ivss_echo                 ,-
    lvesd_echo                ,-
    pws_echo                  ,-
    lads_echo                 ,-
    e_echo                    ,-
    a_echo                    ,-
    dt_echo                   ,-
    e_lateral_echo            ,-
    a_lateral_echo            ,-
    sm_lateral_echo           ,-
    e_septal_echo             ,-
    a_septal_echo             ,-
    sm_septal_echo            ,-
    e_rv_echo                 ,-
    a_rv_echo                 ,-
    sm_rv_echo                ,-
    lvedv_4ch_echo            ,-
    lvesv_4ch_echo            ,-
    rvedd_echo                ,-
    rasv_echo                 ,-
    ravm_echo                 ,-
    radv_echo                 ,-
    lasv_4ch_echo             ,-
    lavm_4ch_echo             ,-
    ladv_4ch_echo             ,-
    lvedv_2ch_echo            ,-
    lvesv_2ch_echo            ,-
    lasv_2ch_echo             ,-
    lavm_2ch_echo             ,-
    ladv_2ch_echo             ,-
    mr_echo                   ,-
    ar_echo                   ,-
    pr_echo                   ,-
    tr_echo                   ,-
    rvsp_echo                 ,-
    rvwmsi_a_echo             ,-
    rvwmsi_m_echo             ,-
    rvwmsi_b_echo             ,-
    ao_valve_echo             ,-
    `Baseline LVEF`)

# Remove columns with more than 39% missing values

clin_dat=clin_dat[,-as.numeric(which(t(rbind(clin_dat) %>% summarize_all(~sum(is.na(.))/n()>=.39))[,1]))]

# Fixing column  names 
clin_dat = clin_dat %>% mutate_at(vars(colnames(clin_dat)[which(grepl("1=",colnames(clin_dat)))]),~as.factor(.))
clin_dat = clin_dat %>% mutate_if(is.character,~as.factor(.))


clin_dat = clin_dat %>% select(-`Device therapy pre-LVAD (0=no, 1=yes)`)



colnames(clin_dat)=gsub("\\s*\\([^\\)]+\\)","",colnames(clin_dat))
colnames(clin_dat)=gsub("-","_",colnames(clin_dat))
colnames(clin_dat)=gsub("/","_",colnames(clin_dat))
colnames(clin_dat)=gsub(" ","",colnames(clin_dat))



clin_dat=clin_dat %>% mutate(LVADtype=fct_other(LVADtype,keep=c("HeartMate2","HeartMate3","HeartWare")))

clin_dat = clin_dat %>% mutate_at(vars(c("study_id","hospital")),~as.character(.))


# Dummy coding categorical variables

dmy <- dummyVars(" ~ .", data = clin_dat %>% select_if(is.factor),fullRank = T)
trsf <- data.frame(predict(dmy, newdata = clin_dat %>% select_if(is.factor)))
clin_dat=data.frame(clin_dat[,c(1,2)],cbind(trsf,clin_dat %>% select_if(is.numeric)))

clin_dat %>% summarize_all(~sum(is.na(.)))

### Association analysis
rna_reg=1:(dim(rna_dat)[2]-2) %>% purrr::map_df(function(x) data.frame(summary(glm(rna_dat$resp~rna_dat[[x+2]],family="binomial"))$coef)[2,])
write.csv(rna_reg %>% as_tibble() %>% mutate(RNV_Variable=colnames(rna_dat)[3:(dim(rna_dat)[2])]) %>% arrange(Pr...z..)  %>% select(RNV_Variable,Estimate,Std..Error,p.value=Pr...z..) %>% filter(p.value<.05),file="rna_resp.csv")

rna_LVEF=1:(dim(rna_dat)[2]-2) %>% purrr::map_df(function(x) data.frame(summary(lm(clin_dat$Post_LVADLVEF~rna_dat[[x+2]]))$coef)[2,])
write.csv(rna_LVEF %>% as_tibble() %>% mutate(RNV_Variable=colnames(rna_dat)[3:(dim(rna_dat)[2])]) %>% arrange(Pr...t..)  %>% select(RNV_Variable,Estimate,Std..Error,p.value=Pr...t..) %>% filter(p.value<.05),file="rna_lvef.csv")

rna_LVEDD=1:(dim(rna_dat)[2]-2) %>% purrr::map_df(function(x) data.frame(summary(lm(clin_dat$Post_LVADLVEDDatmaxLVEF~rna_dat[[x+2]]))$coef)[2,])
write.csv(rna_LVEDD %>% as_tibble() %>% mutate(RNV_Variable=colnames(rna_dat)[3:(dim(rna_dat)[2])]) %>% arrange(Pr...t..)  %>% select(RNV_Variable,Estimate,Std..Error,p.value=Pr...t..) %>% filter(p.value<.05),file="rna_lvedd.csv")

clin_reg=c(3:52,54) %>% purrr::map_df(function(x) data.frame(summary(glm(clin_dat$resp~clin_dat[[x]],family="binomial"))$coef)[2,])
write.csv(clin_reg %>% as_tibble() %>% mutate(Clin_Variable=colnames(clin_dat)[c(3:52,54)]) %>% arrange(Pr...z..)  %>% select(Clin_Variable,Estimate,Std..Error,p.value=Pr...z..),file="clin_resp.csv")

clin_LVEF=c(3:52,54) %>% purrr::map_df(function(x) data.frame(summary(lm(clin_dat$Post_LVADLVEF~clin_dat[[x]]))$coef)[2,])
write.csv(clin_LVEF %>% as_tibble() %>% mutate(Clin_Variable=colnames(clin_dat)[c(3:52,54)]) %>% arrange(Pr...t..)  %>% select(Clin_Variable,Estimate,Std..Error,p.value=Pr...t..),file="clin_lvef.csv")

clin_LVEDD=c(3:52,54) %>% purrr::map_df(function(x) data.frame(summary(lm(clin_dat$Post_LVADLVEDDatmaxLVEF~clin_dat[[x]]))$coef)[2,])
write.csv(clin_LVEDD %>% as_tibble() %>% mutate(Clin_Variable=colnames(clin_dat)[c(3:52,54)]) %>% arrange(Pr...t..)  %>% select(Clin_Variable,Estimate,Std..Error,p.value=Pr...t..),file="clin_lvedd.csv")
###

# Using multiple imputation with chained equations to impute missing values

MI_clin_dat=mice(clin_dat %>% select(-resp,-study_id,-hospital,-Post_LVADLVEF,-Post_LVADLVEDDatmaxLVEF),5)

MI_Clin=data.frame(resp=clin_dat$resp,study_id=clin_dat$study_id,
                   hospital=clin_dat$hospital,
                   Post_LVADLVEF=clin_dat$Post_LVADLVEF,
                   Post_LVADLVEDDatmaxLVEF=clin_dat$Post_LVADLVEDDatmaxLVEF,
                   complete(MI_clin_dat,1:5))


MI_RNA=rna_dat[match(MI_Clin$study_id,rna_dat$study_id),]

#Cross Validation 
MI_RNA=rbind(train_rna_dat,test_rna_dat)[match(MI_Clin$study_id,rbind(train_rna_dat,test_rna_dat)$study_id),]

# Splitting data into training and testing sets using proportional stratified sampling

folds=1:100 %>% purrr::map(~
                             (MI_Clin %>% select(resp,study_id,hospital) %>% distinct %>% group_by(resp,hospital) %>% dplyr::slice_sample(prop=.1))$study_id
                           )#createFolds(unique(MI_Clin$study_id),k=5) %>% purrr::map(~unique(MI_Clin$study_id)[.])



# Save folds and data for future use

saveRDS(list(MI_Clin,MI_RNA,folds),file="lvad_dats_folds_3_08_22.rds")


# test ICAR

datICAR=readRDS("datICAR.rds")
folds=readRDS("lvad_dats_folds_3_08_22.rds")[[3]]

datICAR=datICAR %>% mutate(age=as.numeric(age>=50),hfe=as.numeric(heartfailureetiology.1==1),hfsd=as.numeric(heartfailuresymptomsduration>=24),dt=as.numeric(Devicetherapypre_LVAD.1==1|Devicetherapypre_LVAD.2==1),
                   cpl=as.numeric(creatininepre_LVAD>1.2),le=as.numeric(lvedd_echo>=6.5))


datICAR=readRDS("lvad_dats_folds_3_08_22.rds")[[1]] %>% left_join(clin_dat %>% select(study_id,Devicetherapypre_LVAD.1,Devicetherapypre_LVAD.2)) %>%
  filter(!is.na(Devicetherapypre_LVAD.1),!is.na(Devicetherapypre_LVAD.2))

datICAR=datICAR %>% mutate(age=as.numeric(age>=50),hfe=as.numeric(heartfailureetiology.1==1),hfsd=as.numeric(heartfailuresymptomsduration>=24),dt=as.numeric(Devicetherapypre_LVAD.1==1|Devicetherapypre_LVAD.2==1),
                           cpl=as.numeric(creatininepre_LVAD>1.2),le=as.numeric(lvedd_echo>=6.5))

auc=folds %>% purrr::map(function(x){mod=glm(resp~hfe+hfsd+dt+cpl+le,data=datICAR %>% filter(!(study_id %in% x)),family="binomial")
            as.numeric(roc(datICAR$resp[which(datICAR$study_id %in% x)],predict(mod,datICAR %>% filter((study_id %in% x))))$auc)

roc=folds %>% purrr::map(function(x){mod=glm(resp~hfe+hfsd+dt+cpl+le,data=datICAR %>% filter(!(study_id %in% x)),family="binomial")
            roc=roc(datICAR$resp[which(datICAR$study_id %in% x)],predict(mod,datICAR %>% filter((study_id %in% x))))
            data.frame(se=roc$se,sp=roc$sp)})
roc_ICAR=bind_rows(roc,.id="iter")
                                                           })
mean(unlist(auc)) #0.735
sd(unlist(auc)) #0.134



auc_rf=(bind_rows(bin_auc) %>% filter(source==2,model=="rf",nvar==20))$auc
diff=data.frame(rf=auc_rf,icar=unlist(auc)) %>% mutate(diff=icar-rf) %>% dplyr::summarize(diff=mean(diff))

plan(multisession)
dist=1:10000 %>% future_map(~t(apply(data.frame(rf=auc_rf,icar=unlist(auc)),1,function(x) x[sample(1:2,2)])) %>% as.data.frame %>% mutate(diff=V2-V1) %>% dplyr::summarize(mean(diff)),.progress=T)

1-mean(diff$diff > abs(as.numeric(unlist(dist))))
qplot(as.numeric(unlist(dist)))

### Function for getting responder important variables and cross-validation 
get_imps=function(MI_Clin,MI_RNA){
MI_Clin=MI_Clin %>% dplyr::select(-study_id)
MI_RNA=MI_RNA %>% dplyr::select(-study_id)


out_clin=colnames(MI_Clin)[-(1:4)] %>% purrr::map(function(x) {
  y=summary(glm(MI_Clin$resp~MI_Clin[[x]],family="binomial"))
  y$coefficients[2,4]})

out_rna=colnames(MI_RNA)[-1] %>% purrr::map(function(x) {
  y=summary(glm(MI_RNA$resp~MI_RNA[[x]],family="binomial"))
  y$coefficients[2,4]})



df_imps_rna1=data.frame(var=colnames(MI_RNA)[-1],P=c(unlist(out_rna))) %>% arrange(P)
df_imps_clin1=data.frame(var=colnames(MI_Clin)[-(1:4)],
                         P=unlist(out_clin)) %>% arrange(P)

df_imps=rbind(df_imps_rna1,df_imps_clin1) %>% arrange(P)



rf_dat=cbind(MI_Clin %>% select(one_of(df_imps$var[1:500])), MI_RNA %>% select(one_of(df_imps$var[1:500])))
rf=cforest(resp ~ .,data=data.frame(resp=MI_Clin$resp,rf_dat),controls=cforest_unbiased(ntree=5000,mtry=floor((dim(rf_dat)[2]-1)/3)))
rf_imp=permimp(rf,conditional=T,AUC=T)            
df_imps2=data.frame(var=names(rf_imp$values),imps=as.numeric(rf_imp$values)) %>% arrange(desc(imps))

rf_dat=MI_Clin %>% select(one_of(df_imps_clin1$var[1:500]))
rf=cforest(resp ~ .,data=data.frame(resp=MI_Clin$resp,rf_dat),controls=cforest_unbiased(ntree=5000,mtry=floor((dim(rf_dat)[2]-1)/3)))
rf_imp=permimp(rf,conditional=T,AUC=T)            
df_imps_clin=data.frame(var=names(rf_imp$values),imps=as.numeric(rf_imp$values)) %>% arrange(desc(imps))

rf_dat=MI_RNA %>% select(one_of(df_imps_rna1$var[1:500]))
rf=cforest(resp ~ .,data=data.frame(resp=MI_Clin$resp,rf_dat),controls=cforest_unbiased(ntree=5000,mtry=floor((dim(rf_dat)[2]-1)/3)))
rf_imp=permimp(rf,conditional=T,AUC=T)            
df_imps_rna=data.frame(var=names(rf_imp$values),imps=as.numeric(rf_imp$values)) %>% arrange(desc(imps))

return(list(df_imps2,df_imps_clin,df_imps_rna))


cross_valid = function(imps,tstData,nvars){
out=
  bind_rows(imps %>% purrr::map(function(z){
  train=cbind(MI_Clin %>% select(study_id,resp,one_of(z$var[1:nvars])), MI_RNA %>% select(one_of(z$var[1:nvars]))) %>% filter(!(study_id %in% tstData[[fold]]))
  test=cbind(MI_Clin %>% select(study_id,resp,one_of(z$var[1:nvars])), MI_RNA %>% select(one_of(z$var[1:nvars]))) %>% filter((study_id %in% tstData[[fold]]))
  rf1=cforest(resp ~ .,data=data.frame(train %>% select(-study_id)),controls=cforest_unbiased(ntree=5000,mtry=floor((dim(train)[2]-1)/3)))
  glm1=glm(resp~.,data=train %>% select(-study_id),family="binomial")
  svm1=e1071::svm(resp~.,data=train %>% select(-study_id),probability=T)
  gbm1=gbm(resp~.,data=train %>% select(-study_id) %>% mutate(resp=as.character(resp)),n.trees=5000,interaction.depth = 3)
  pred_rf=matrix(unlist(predict(rf1,newdata=test,type="prob")),ncol=2,byrow=T)[,2]
  pred_glm=as.numeric(predict(glm1,newdata=test,type="response"))
  pred_svm=as.numeric(attr(predict(svm1,newdata=test,probability=T),"probabilities")[,1])
  pred_gbm=predict(gbm1,newdata=test %>% mutate(resp=as.character(resp)),n.trees=5000,type="response")
  bind_rows(map2(list(pred_rf,pred_glm,pred_svm,pred_gbm),c("rf","glm","svm","gbm"),function(x,y) data.frame(model=y,pred=x,nvar=nvars,fold=fold,source=z$num[1],truth=test$resp)))
  #data.frame(fold=fold,source=z$num[1],nvar=nvars,model=c("rf","glm","svm","gbm"),auc=unlist(list(pred_rf,pred_glm,pred_svm,pred_gbm) %>% purrr::map(~as.numeric(roc(test$resp,.)$auc))))
  }))
}

}

# Function for getting LVEF important variables and cross-validation
get_imps=function(MI_Clin,MI_RNA){
MI_Clin=MI_Clin %>% dplyr::select(-study_id)
MI_RNA=MI_RNA %>% dplyr::select(-study_id)


out_clin=colnames(MI_Clin)[-(1:4)] %>% purrr::map(function(x) {
  y=summary(glm(MI_Clin$Post_LVADLVEF~MI_Clin[[x]]))
  y$coefficients[2,4]})

out_rna=colnames(MI_RNA)[-1] %>% purrr::map(function(x) {
  y=summary(glm(MI_Clin$Post_LVADLVEF~MI_RNA[[x]]))
  y$coefficients[2,4]})



df_imps_rna1=data.frame(var=colnames(MI_RNA)[-1],P=c(unlist(out_rna))) %>% arrange(P)
df_imps_clin1=data.frame(var=colnames(MI_Clin)[-(1:4)],
                         P=unlist(out_clin)) %>% arrange(P)

df_imps=rbind(df_imps_rna1,df_imps_clin1) %>% arrange(P)



rf_dat=cbind(MI_Clin %>% select(one_of(df_imps$var[1:500])), MI_RNA %>% select(one_of(df_imps$var[1:500])))
rf=cforest(Post_LVADLVEF ~ .,data=data.frame(Post_LVADLVEF=MI_Clin$Post_LVADLVEF,rf_dat),controls=cforest_unbiased(ntree=5000,mtry=floor((dim(rf_dat)[2]-1)/3)))
rf_imp=permimp(rf,conditional=T)            
df_imps2=data.frame(var=names(rf_imp$values),imps=as.numeric(rf_imp$values)) %>% arrange(desc(imps))

rf_dat=MI_Clin %>% select(one_of(df_imps_clin1$var[1:500]))
rf=cforest(Post_LVADLVEF ~ .,data=data.frame(Post_LVADLVEF=MI_Clin$Post_LVADLVEF,rf_dat),controls=cforest_unbiased(ntree=5000,mtry=floor((dim(rf_dat)[2]-1)/3)))
rf_imp=permimp(rf,conditional=T)            
df_imps_clin=data.frame(var=names(rf_imp$values),imps=as.numeric(rf_imp$values)) %>% arrange(desc(imps))

rf_dat=MI_RNA %>% select(one_of(df_imps_rna1$var[1:500]))
rf=cforest(Post_LVADLVEF ~ .,data=data.frame(Post_LVADLVEF=MI_Clin$Post_LVADLVEF,rf_dat),controls=cforest_unbiased(ntree=5000,mtry=floor((dim(rf_dat)[2]-1)/3)))
rf_imp=permimp(rf,conditional=T)            
df_imps_rna=data.frame(var=names(rf_imp$values),imps=as.numeric(rf_imp$values)) %>% arrange(desc(imps))

return(list(df_imps2,df_imps_clin,df_imps_rna))
}


cross_valid = function(imps,tstData,nvars){
out=
  bind_rows(imps %>% purrr::map(function(z){
  train=cbind(MI_Clin %>% select(study_id,Post_LVADLVEF,one_of(z$var[1:nvars])), MI_RNA %>% select(one_of(z$var[1:nvars]))) %>% filter(!(study_id %in% tstData[[fold]]))
  test=cbind(MI_Clin %>% select(study_id,Post_LVADLVEF,one_of(z$var[1:nvars])), MI_RNA %>% select(one_of(z$var[1:nvars]))) %>% filter((study_id %in% tstData[[fold]]))
  rf1=cforest(Post_LVADLVEF ~ .,data=data.frame(train %>% select(-study_id)),controls=cforest_unbiased(ntree=5000,mtry=floor((dim(train)[2]-1)/3)))
  glm1=glm(Post_LVADLVEF~.,data=train %>% select(-study_id))
  svm1=e1071::svm(Post_LVADLVEF~.,data=train %>% select(-study_id))
  gbm1=gbm(Post_LVADLVEF~.,data=train %>% select(-study_id),n.trees=5000,interaction.depth = 3)
  pred_rf=as.numeric(predict(rf1,newdata=test))
  pred_glm=as.numeric(predict(glm1,newdata=test))
  pred_svm=as.numeric(predict(svm1,newdata=test))
  pred_gbm=predict(gbm1,newdata=test,n.trees=5000)
  bind_rows(map2(list(pred_rf,pred_glm,pred_svm,pred_gbm),c("rf","glm","svm","gbm"),function(x,y) data.frame(model=y,pred=x,nvar=nvars,fold=fold,source=z$num[1],truth=test$Post_LVADLVEF)))
  #data.frame(fold=fold,source=z$num[1],nvar=nvars,model=c("rf","glm","svm","gbm"),auc=unlist(list(pred_rf,pred_glm,pred_svm,pred_gbm) %>% purrr::map(~as.numeric(roc(test$Post_LVADLVEF,.)$auc))))
  }))
}


# Function for getting LVEDD important variables and cross-validation 
get_imps=function(MI_Clin,MI_RNA){
MI_Clin=MI_Clin %>% dplyr::select(-study_id)
MI_RNA=MI_RNA %>% dplyr::select(-study_id)

idx_rem=which(is.na(MI_Clin$Post_LVADLVEDDatmaxLVEF))
MI_Clin=MI_Clin[-idx_rem,]
MI_RNA=MI_RNA[-idx_rem,]

out_clin=colnames(MI_Clin)[-(1:4)] %>% purrr::map(function(x) {
  y=summary(glm(MI_Clin$Post_LVADLVEDDatmaxLVEF~MI_Clin[[x]]))
  y$coefficients[2,4]})

out_rna=colnames(MI_RNA)[-1] %>% purrr::map(function(x) {
  y=summary(glm(MI_Clin$Post_LVADLVEDDatmaxLVEF~MI_RNA[[x]]))
  y$coefficients[2,4]})



df_imps_rna1=data.frame(var=colnames(MI_RNA)[-1],P=c(unlist(out_rna))) %>% arrange(P)
df_imps_clin1=data.frame(var=colnames(MI_Clin)[-(1:4)],
                         P=unlist(out_clin)) %>% arrange(P)

df_imps=rbind(df_imps_rna1,df_imps_clin1) %>% arrange(P)



rf_dat=cbind(MI_Clin %>% select(one_of(df_imps$var[1:500])), MI_RNA %>% select(one_of(df_imps$var[1:500])))
rf=cforest(Post_LVADLVEDDatmaxLVEF ~ .,data=data.frame(Post_LVADLVEDDatmaxLVEF=MI_Clin$Post_LVADLVEDDatmaxLVEF,rf_dat),controls=cforest_unbiased(ntree=5000,mtry=floor((dim(rf_dat)[2]-1)/3)))
rf_imp=permimp(rf,conditional=T)            
df_imps2=data.frame(var=names(rf_imp$values),imps=as.numeric(rf_imp$values)) %>% arrange(desc(imps))

rf_dat=MI_Clin %>% select(one_of(df_imps_clin1$var[1:500]))
rf=cforest(Post_LVADLVEDDatmaxLVEF ~ .,data=data.frame(Post_LVADLVEDDatmaxLVEF=MI_Clin$Post_LVADLVEDDatmaxLVEF,rf_dat),controls=cforest_unbiased(ntree=5000,mtry=floor((dim(rf_dat)[2]-1)/3)))
rf_imp=permimp(rf,conditional=T)            
df_imps_clin=data.frame(var=names(rf_imp$values),imps=as.numeric(rf_imp$values)) %>% arrange(desc(imps))

rf_dat=MI_RNA %>% select(one_of(df_imps_rna1$var[1:500]))
rf=cforest(Post_LVADLVEDDatmaxLVEF ~ .,data=data.frame(Post_LVADLVEDDatmaxLVEF=MI_Clin$Post_LVADLVEDDatmaxLVEF,rf_dat),controls=cforest_unbiased(ntree=5000,mtry=floor((dim(rf_dat)[2]-1)/3)))
rf_imp=permimp(rf,conditional=T)            
df_imps_rna=data.frame(var=names(rf_imp$values),imps=as.numeric(rf_imp$values)) %>% arrange(desc(imps))

return(list(df_imps2,df_imps_clin,df_imps_rna))
}


cross_valid = function(imps,tstData,nvars){
out=
  bind_rows(imps %>% purrr::map(function(z){
  train=cbind(MI_Clin %>% select(study_id,Post_LVADLVEDDatmaxLVEF,one_of(z$var[1:nvars])), MI_RNA %>% select(one_of(z$var[1:nvars]))) %>% filter(!(study_id %in% tstData[[fold]]))
  test=cbind(MI_Clin %>% select(study_id,Post_LVADLVEDDatmaxLVEF,one_of(z$var[1:nvars])), MI_RNA %>% select(one_of(z$var[1:nvars]))) %>% filter((study_id %in% tstData[[fold]]))
  rf1=cforest(Post_LVADLVEDDatmaxLVEF ~ .,data=data.frame(train %>% select(-study_id)),controls=cforest_unbiased(ntree=5000,mtry=floor((dim(train)[2]-1)/3)))
  glm1=glm(Post_LVADLVEDDatmaxLVEF~.,data=train %>% select(-study_id))
  svm1=e1071::svm(Post_LVADLVEDDatmaxLVEF~.,data=train %>% select(-study_id))
  gbm1=gbm(Post_LVADLVEDDatmaxLVEF~.,data=train %>% select(-study_id),n.trees=5000,interaction.depth = 3)
  pred_rf=as.numeric(predict(rf1,newdata=test))
  pred_glm=as.numeric(predict(glm1,newdata=test))
  pred_svm=as.numeric(predict(svm1,newdata=test))
  pred_gbm=predict(gbm1,newdata=test,n.trees=5000)
  bind_rows(map2(list(pred_rf,pred_glm,pred_svm,pred_gbm),c("rf","glm","svm","gbm"),function(x,y) data.frame(model=y,pred=x,nvar=nvars,fold=fold,source=z$num[1],truth=test$Post_LVADLVEDDatmaxLVEF)))
  #data.frame(fold=fold,source=z$num[1],nvar=nvars,model=c("rf","glm","svm","gbm"),auc=unlist(list(pred_rf,pred_glm,pred_svm,pred_gbm) %>% purrr::map(~as.numeric(roc(test$Post_LVADLVEDDatmaxLVEF,.)$auc))))
  }))
}


### Summarizing results
library(tidyverse)
library(pROC)
library(furrr)
library(rms)
library(patchwork)
library(party)
library(pdp)
library(PRROC)

#Read in original data with splits
dats=readRDS("lvad_dats_folds_3_08_22.rds")

#Read in importance and CV results (usually run on cluster)
dat_bin = 1:100 %>% purrr::map(~readRDS(paste0("bin_fold100_v2/res100fold_",.,".rds")))

# Produce ROC curve using average values from each fold
roc_avg=dat_bin %>% purrr::map(function(x) {
  best=x %>% filter(model=="gbm",nvar==80)
  roccurve=with(best,roc(truth,pred))
  data.frame(se=roccurve$sensitivities,sp=roccurve$specificities)
})

roc_dat=bind_rows(roc_avg,.id="iter") %>% mutate(sp_min=1-sp)

roc_smth=roc_avg %>% purrr::map(function(x) smooth.spline(x$sp,x$se,df=4) %>% predict() %>% as.data.frame())

roc_ICAR=roc_ICAR %>% mutate(sp_min=1-sp)

png(file="roc_curve.png",width=800,height=800)
smooth.spline(roc_dat$sp_min,roc_dat$se,df=4) %>% predict() %>% as.data.frame() %>% ggplot(aes(x=x,y=y)) + geom_line(size=.8,color="blue") +
  geom_segment(aes(x=0,xend=1,y=0,yend=1)) + 
  #geom_line(data=bind_rows(roc_smth,.id="iter"),aes(x=1-x,y=y,group=iter),alpha=.1) + 
  geom_line(data=bind_rows(roc_avg,.id="iter"),aes(x=1-se,y=sp,group=iter),alpha=.15) + 
  geom_line(data=smooth.spline(roc_ICAR$sp_min,roc_ICAR$se,df=4) %>% predict() %>% as.data.frame(),aes(x=x,y=y),color="red") +

  theme_bw() + xlab("1-Specificity") + ylab("Sensitivity") + 
  scale_x_continuous(limits=c(0,1),breaks=seq(0,1,by=.1)) + scale_y_continuous(limits=c(0,1),breaks=seq(0,1,by=.1))
dev.off()

# Create figure showing AUC and PRAUC for each model for binary response outcome 
bin_auc=dat_bin %>% purrr::map(. %>% group_by(nvar,fold,source,model) %>% dplyr::summarize(auc=as.numeric(roc(truth,pred)$auc),
                                                                                           prauc=as.numeric(pr.curve(pred[which(truth==1)],pred[which(truth==0)])$auc.integral)))

bind_rows(bin_auc) %>% group_by(nvar,source,model) %>% dplyr::summarize(aucs=mean(auc),sdevs=sd(auc)) %>% group_by(source) %>% 
  filter(aucs==max(aucs))

Fig1=bind_rows(bin_auc) %>% group_by(nvar,source,model) %>% dplyr::summarize(mnauc=mean(auc),lower=quantile(auc,.25),upper=quantile(auc,.75),std=sd(auc)) %>% 
  mutate(Model=case_when(model=="gbm"~"GBM",model=="glm"~"LR",model=="rf"~"RF",model=="svm"~"SVM"),
    source=c("Both","Clinical Only","RNA Only")[source]) %>%
  ggplot(aes(x=nvar,y=mnauc,color=Model,group=model,ymin=lower,ymax=upper)) + geom_line() + geom_point() + #geom_errorbar(alpha=.25) + 
  facet_wrap(~source) + theme_bw() + 
  #scale_y_continuous(labels=seq(.55,.75,by=.05),breaks=seq(.55,.75,by=.05),limits=c(.55,.75)) + 
  ylab("AUC") + xlab("No. of Variables") + 
  scale_color_viridis_d(end=.95);Fig1

Fig1b=bind_rows(bin_auc) %>% group_by(nvar,source,model) %>% dplyr::summarize(mnprauc=mean(prauc),std=sd(prauc)) %>% 
  mutate(Model=case_when(model=="gbm"~"GBM",model=="glm"~"LR",model=="rf"~"RF",model=="svm"~"SVM"),
         source=c("Both","Clinical Only","RNA Only")[source]) %>%
  ggplot(aes(x=nvar,y=mnprauc,color=Model,group=model)) + geom_line() + geom_point() + 
  facet_wrap(~source) + theme_bw() + 
  scale_y_continuous(labels=seq(.2,.5,by=.05),breaks=seq(.2,.5,by=.05),limits=c(.2,.5)) + 
  ylab("AUC") + xlab("No. of Variables") + 
  scale_color_viridis_d(end=.95);Fig1b

ggsave("Fig1_v2.png",units="in",width=10,height=5)

# Assess calibration of models

cals=bind_rows(dat_bin) %>% filter(nvar==20,source==2,model=="rf") %>% dplyr::mutate(logit=qlogis(pred)) %>% filter(logit!=0,!is.infinite(logit)) %>%
  group_by(fold) %>% #filter(!all(truth==1) & !all(truth==0)) %>%
  dplyr::summarize(ints=summary(glm(truth~1,offset=logit,family="binomial"))$coef[1,1],
            slps=summary(glm(truth~logit,family="binomial"))$coef[2,1],
            intlow=as.numeric(confint.lm(glm(truth~1,offset=logit,family="binomial"))[1]),
            intup=as.numeric(confint.lm(glm(truth~1,offset=logit,family="binomial"))[2]),
            slplow=as.numeric(confint.lm(glm(truth~logit,family="binomial"))[2,1]),
            slpup=as.numeric(confint.lm(glm(truth~logit,family="binomial"))[2,2]))


tst=bind_rows(dat_bin) %>% filter(nvar==80,source==1,model=="gbm") %>% filter(fold==2) %>% mutate(logit=qlogis(pred)) %>% filter(logit!=0,!is.infinite(logit))
summary(glm(truth~1,offset=logit,family="binomial",data=tst))$coef[1,1]

bind_rows(dat_bin) %>% filter(nvar==80,source==1,model=="gbm") %>% ggplot(aes(x=pred,y=jitter(as.numeric(truth)))) + geom_point()

tst=bind_rows(dat_bin) %>% filter(nvar==80,source==1,model=="gbm") %>% dplyr::mutate(logit=qlogis(pred)) %>% filter(logit!=0,!is.infinite(logit)) %>%
  group_by(fold) %>% #filter(fold==2)
dplyr::summarize(intlow=as.numeric(confint.lm(glm(truth~1,offset=logit,family="binomial",data=tst))[1]))


apply(cals,
     2,mean)
ids=bind_rows(1:100 %>% purrr::map(function(x) dats[[1]] %>% filter(study_id %in% dats[[3]][[x]]) %>% select(study_id,lvef_echo,lvedd_echo)))

df=data.frame(bind_rows(dat_bin) %>% filter(nvar==20,source==2,model=="rf"),ids)
df=df %>% group_by(study_id) %>% 
  dplyr::summarize(pred=mean(pred)) %>% left_join(df %>% select(study_id,truth) %>% distinct)
df$truth=as.numeric(df$truth)-1
obs_prop=as.numeric(unlist(df$pred %>% purrr::map(function(z)  
  df %>% filter(pred>z-.05,pred<z+.05) %>% dplyr::summarize(mean(truth)))))

a=ggplot(df,aes(x=pred,y=obs_prop)) + geom_point() + geom_line(data=data.frame(x=c(0,1),y=c(0,1)),aes(x=x,y=y)) + theme_bw() + 
  xlab("RF Prediction") + ylab("Observed Proportion");a

print(val.prob(p=df$pred,y=df$truth, pl=TRUE, smooth=TRUE, logistic.cal=TRUE,statloc=F))
ggsave("Cal_Fig_v2.png",units="in",width=5,height=5)


# Summarize results for LVEF model 
dat_LVEF = 1:100 %>% purrr::map(function(x) readRDS(paste0("LVEF_fold100_v2/resLVEF70_100fold",x,".rds")))
#dat_LVEF2 = 1:100 %>% purrr::map(~readRDS(paste0("LVEF_fold100/resLVEF70_100fold",.,".rds")))


LVEF_MAE=c(dat_LVEF) %>% purrr::map(. %>% group_by(nvar,fold,source,model) %>% dplyr::summarize(mae=mean(abs(pred-truth))))

bind_rows(LVEF_MAE) %>% group_by(nvar,source,model) %>% dplyr::summarize(maes=mean(mae),sdevs=sd(mae)) %>% group_by(source) %>% 
  filter(maes==min(maes))

bind_rows(LVEF_MAE) %>% group_by(nvar,source,model) %>% dplyr::summarize(mae=mean(mae)) %>% 
  mutate(source=c("Both","Clinical Only","RNA Only")[source],
         Model=case_when(model=="gbm"~"GBM",model=="glm"~"LR",model=="rf"~"RF",model=="svm"~"SVM")) %>%
  ggplot(aes(x=nvar,y=mae,color=Model,group=Model)) + geom_line() + geom_point() +
  facet_wrap(~source) + theme_bw() + 
  ylab("MAE") + xlab("No. of Variables") + scale_color_viridis_d(end=.95) + 
  scale_y_continuous(limits=c(9,17),breaks=9:17,labels=9:17)

ggsave("Fig2_v2.png",units="in",width=10,height=5)


df=data.frame(bind_rows(c(dat_LVEF)) %>% filter(nvar==10,source==2,model=="rf"),ids)

df=df %>% group_by(study_id) %>% 
  dplyr::summarize(pred=mean(pred)) %>% left_join(df %>% select(study_id,truth,lvef_echo,lvedd_echo) %>% distinct)


a=df %>% arrange(truth) %>%
  ggplot(aes(x=lvef_echo,xend=truth,y=1:194,yend=1:194)) + geom_segment() + geom_point() + 
  geom_point(aes(x=truth),color="blue") +
  geom_point(aes(x=pred),color="red") + xlab("LVEF") + ylab("Test Patient") +
  theme_bw()
b=df %>% 
    ggplot(aes(x=truth,y=truth-pred)) + geom_point() + theme_bw() + ylab("Difference Between Truth and Predicted") + xlab("Post-LVAD LVEF")

a +b
ggsave("Fig3_v3.png",units="in",width=8,height=8)

# Summarize results for LVEDD model 
dat_LVEDD = 1:100 %>% purrr::map(~readRDS(paste0("LVEDD_fold100_v2/resLVEDD70_100fold",.,".rds")))
#dat_LVEDD2 = 1:100 %>% purrr::map(~readRDS(paste0("LVEDD_fold100/resLVEDD70_100fold",.,".rds")))
MI_Clin=dats[[1]] 
MI_RNA=dats[[2]]
tstData=dats[[3]]

idx_rem=which(is.na(MI_Clin$Post_LVADLVEDDatmaxLVEF))

LVEDD_MAE=cbind(dat_LVEDD) %>% purrr::map(. %>% group_by(nvar,fold,source,model) %>% dplyr::summarize(mae=mean(abs(pred-truth))))

bind_rows(LVEDD_MAE) %>% group_by(nvar,source,model) %>% dplyr::summarize(maes=mean(mae),sdevs=sd(mae)) %>% group_by(source) %>% 
  filter(maes==min(maes))

bind_rows(LVEDD_MAE) %>% group_by(nvar,source,model) %>% dplyr::summarize(mae=mean(mae)) %>% 
  mutate(source=c("Both","Clinical Only","RNA Only")[source],
         Model=case_when(model=="gbm"~"GBM",model=="glm"~"LR",model=="rf"~"RF",model=="svm"~"SVM")) %>%
  ggplot(aes(x=nvar,y=mae,color=Model,group=Model)) + geom_line() + geom_point() +
  facet_wrap(~source) + theme_bw() + 
  ylab("MAE") + xlab("No. of Variables") + scale_color_viridis_d(end=.95) + 
  scale_y_continuous(limits=c(.8,1.5),breaks=seq(.5,1.5,by=.1),labels=seq(.5,1.5,by=.1))

ggsave("Fig4_v2.png",units="in",width=10,height=5)

ids=bind_rows(1:100 %>% purrr::map(function(x) dats[[1]][-idx_rem,] %>% filter(study_id %in% dats[[3]][[x]]) %>% select(study_id,lvef_echo,lvedd_echo)))

df=data.frame(bind_rows(c(dat_LVEDD)) %>% filter(nvar==6,source==2,model=="rf"),ids)

df=df %>% group_by(study_id) %>% 
  dplyr::summarize(pred=mean(pred)) %>% left_join(df %>% select(study_id,truth,lvef_echo,lvedd_echo) %>% distinct)


a=df %>% arrange(truth) %>%
  ggplot(aes(x=lvedd_echo,xend=truth,y=1:189,yend=1:189)) + geom_segment() + geom_point() + 
  geom_point(aes(x=truth),color="blue") +
  geom_point(aes(x=pred),color="red") + xlab("LVEDD") + ylab("Test Patient") +
  theme_bw()
b=df %>% 
  ggplot(aes(x=truth,y=truth-pred)) + geom_point() + theme_bw() + ylab("Difference Between Truth and Predicted") + xlab("Post-LVAD LVEDD")

a +b
ggsave("Fig5_v3.png",units="in",width=8,height=8)

### Produce partial dependency plots for each important variable (full data) 
#imps=read_csv("orig_imps/bin_imps_list.csv")$var
imps=readRDS("fulldata_imps_bin.rds")[[2]]$var
MI_Clin=dats[[1]]
MI_RNA=dats[[2]]
rf1_dat=MI_Clin %>% select(one_of(imps[1:10]))
rf2_dat=MI_RNA %>% select(one_of(imps[11:20]))
rf1=cforest(resp ~ .,data=data.frame(resp=MI_Clin$resp,rf1_dat),controls=cforest_control(ntree=5000,mtry=floor((dim(rf1_dat)[2]-1)/3),replace=F))
rf2=cforest(resp ~ .,data=data.frame(resp=MI_Clin$resp,rf2_dat),controls=cforest_control(ntree=5000,mtry=floor((dim(rf2_dat)[2]-1)/3),replace=F))
partial_plots=as.character(imps[1:10]) %>% purrr::map(function(x) rf1 %>% pdp::partial(pred.var = x))
partial_plots2=as.character(imps[11:20]) %>% purrr::map(function(x) rf2 %>% pdp::partial(pred.var = x))

saveRDS(partial_plots,file="clin_pdp_bin_v2.rds")
saveRDS(partial_plots2,file="rna_pdp_bin.rds")


#imps=read_csv("orig_imps/LVEF_imps_list.csv")$var
imps=readRDS("fulldata_imps_lvef.rds")[[2]]$var
MI_Clin=dats[[1]]
MI_RNA=dats[[2]]
rf1_dat=MI_Clin %>% select(one_of(imps[1:10]))
rf2_dat=MI_RNA %>% select(one_of(imps[11:20]))
rf1=cforest(Post_LVADLVEF ~ .,data=data.frame(Post_LVADLVEF=MI_Clin$Post_LVADLVEF,rf1_dat),controls=cforest_control(ntree=5000,mtry=floor((dim(rf1_dat)[2]-1)/3),replace=F))
rf2=cforest(Post_LVADLVEF ~ .,data=data.frame(Post_LVADLVEF=MI_Clin$Post_LVADLVEF,rf2_dat),controls=cforest_control(ntree=5000,mtry=floor((dim(rf2_dat)[2]-1)/3),replace=F))

partial_plots=as.character(imps[1:10]) %>% purrr::map(function(x) rf1 %>% pdp::partial(pred.var = x))
saveRDS(partial_plots,file="clin_pdp_lvef_v2.rds")
partial_plots2=as.character(imps[11:20]) %>% purrr::map(function(x) rf2 %>% pdp::partial(pred.var = x))
saveRDS(partial_plots2,file="rna_pdp_lvef.rds")

#imps=read_csv("orig_imps/LVEDD_imps_list.csv")$var
imps=readRDS("fulldata_imps_lvedd.rds")[[2]]$var

MI_Clin=dats[[1]]
MI_RNA=dats[[2]]
idx_rem=which(is.na(MI_Clin$Post_LVADLVEDDatmaxLVEF))
MI_Clin=MI_Clin[-idx_rem,]
MI_RNA=MI_RNA[-idx_rem,]
rf1_dat=MI_Clin %>% select(one_of(imps[1:10]))
rf2_dat=MI_RNA %>% select(one_of(imps[11:20]))
rf1=cforest(Post_LVADLVEDDatmaxLVEF ~ .,data=data.frame(Post_LVADLVEDDatmaxLVEF=MI_Clin$Post_LVADLVEDDatmaxLVEF,rf1_dat),controls=cforest_control(ntree=500,mtry=floor((dim(rf1_dat)[2]-1)/3),replace=F))
rf2=cforest(Post_LVADLVEDDatmaxLVEF ~ .,data=data.frame(Post_LVADLVEDDatmaxLVEF=MI_Clin$Post_LVADLVEDDatmaxLVEF,rf2_dat),controls=cforest_control(ntree=5000,mtry=floor((dim(rf2_dat)[2]-1)/3),replace=F))

partial_plots=as.character(imps[1:10]) %>% purrr::map(function(x) rf1 %>% pdp::partial(pred.var = x))
saveRDS(partial_plots,file="clin_pdp_lvedd_v2.rds")
partial_plots2=as.character(imps[11:20]) %>% purrr::map(function(x) rf2 %>% pdp::partial(pred.var = x))
saveRDS(partial_plots2,file="rna_pdp_lvedd.rds")

#Get correlations between variables
LVEF_imps=MI_RNA %>% select(
ENSG00000134597,
ENSG00000114654,
ENSG00000140830,
ENSG00000177363,
ENSG00000137804,
ENSG00000205221,
ENSG00000109794,
ENSG00000151690,
ENSG00000162772,
ENSG00000076770)

LVEDD_imps=MI_RNA %>% select(
ENSG00000176907,
ENSG00000107821,
ENSG00000078053,
ENSG00000130037,
ENSG00000159399,
ENSG00000229807,
ENSG00000099725,
ENSG00000138376,
ENSG00000177363,
ENSG00000109794)


rm=which(is.na(MI_Clin %>% select(Post_LVADLVEDDatmaxLVEF)))

write.csv(data.frame(LVEF_cor=cor(cbind(MI_Clin %>% select(Post_LVADLVEF),LVEF_imps))[-1,1]),file="LVEF_cor.csv")
data.frame(LVEDD_cor=cor(cbind(MI_Clin %>% select(Post_LVADLVEDDatmaxLVEF),LVEDD_imps)[-rm,])[-1,1]) %>% write.csv(file="LVEDD_cor.csv")

# Produce partial dependency plots for supplement 
nms=read_csv("nms_v2.csv",col_names = T)$Bin
nms[3]="LVAD Configuration"
tst1=c(readRDS("clin_pdp_bin_v2.rds"),readRDS("rna_pdp_bin.rds"))
tst1[[2]]$MRApre_LVAD.1=c("No","Yes")
tst1[[3]]$LVADmode.2=c("Axial","Centrifugal")
tst1[[4]]$LVADtype.HeartWare=c("Other","HeartWare")
tst1[[6]]$hypetension.1=c("No","Yes")
tst1[[10]]$rh_type.1=c("Negative","Positive")


tst1[11:20]  = tst1[11:20] %>% purrr::map(function(x) {
  x$yhat=1-x$yhat
  return(x)
  })

plts1=map2(tst1,nms, function(x,y){
  #if(is.factor(x[[1]])) x[[1]]=recode(x[[1]],'-1' = "Missing")
  ggplot() + geom_line(aes(x=x[,1],y=x[,2])) + geom_point(aes(x=x[,1],y=x[,2])) + 
    xlab(str_wrap(y,width=15)) + ylab("Marginal Predicted Probability") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust=1))
}
)
wrap_plots(plts1)
ggsave("SupFig1_v2.png",units="in",width=12,height=12)

nms=read_csv("nms_v2.csv",col_names = T)$LVEF
nms[3]="LVAD Configuration"
nms[1]="Sex"
tst2=c(readRDS("clin_pdp_lvef_v2.rds"),readRDS("rna_pdp_lvef.rds"))

tst2[[1]]$gender.2=c("Male","Female")
tst2[[3]]$LVADmode.2=c("Axial","Centrifugal")
tst2[[4]]$rh_type.1=c("Negative","Positive")
tst2[[5]]$LVADtype.HeartWare=c("Other","HeartWare")
tst2[[8]]$MRApre_LVAD.1=c("No","Yes")


plts2=map2(tst2,nms, function(x,y){
  #if(is.factor(x[[1]])) x[[1]]=recode(x[[1]],'-1' = "Missing")
  ggplot() + geom_point(aes(x=x[,1],y=x[,2])) + geom_line(aes(x=x[,1],y=x[,2]))  + 
    xlab(str_wrap(y,width=15)) + ylab("Marginal Predicted Value") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
}
)
wrap_plots(plts2)
ggsave("SupFig2_v2.png",units="in",width=12,height=12)




nms=read_csv("nms_v2.csv",col_names = T)$LVEDD
nms[1]="Sex"
nms[9]="LVAD Configuration"

tst3=c(readRDS("clin_pdp_lvedd_v2.rds"),readRDS("rna_pdp_lvedd.rds"))
tst3[[1]]$gender.2=c("Male","Female")
tst3[[3]]$LVADtype.HeartWare=c("Other","HeartWare")
tst3[[6]]$rh_type.1=c("Negative","Positive")
tst3[[7]]$ACEi_ARB=c("No","Yes")
tst3[[9]]$LVADmode.2=c("Axial","Centrifugal")

tst3[[10]]$diabetes.1=c("No","Yes")


plts3=map2(tst3,nms, function(x,y){
  #if(is.factor(x[[1]])) x[[1]]=recode(x[[1]],'-1' = "Missing")
  ggplot() + geom_line(aes(x=x[,1],y=x[,2])) + geom_point(aes(x=x[,1],y=x[,2])) + 
    xlab(str_wrap(y,width=15)) + ylab("Marginal Predicted Value") + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
}
)
wrap_plots(plts3)
ggsave("SupFig3_v2.png",units="in",width=12,height=12)
