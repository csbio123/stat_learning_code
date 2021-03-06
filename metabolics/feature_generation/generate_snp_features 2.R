library(foreign)
dat_1 <- read.dta("/Users/ti1/Google Drive/study2_snps_as\ predictors/additional_variables.dta")
dat_2 <- read.dta("/Users/ti1/Google Drive/study2_snps_as\ predictors/Psychiatric diagnoses in SELCoH.dta")
dat_3 <- read.dta("/Users/ti1/Google Drive/study2_snps_as\ predictors/SELCoH ancestry scores.dta")
dat_4 <- read.dta("/Users/ti1/Google\ Drive/study2_snps_as\ predictors/selcoh_akt1pathway.dta")



#write.table(dat_1, file="/Users/ti1/Google\ Drive/study2_snps_as\ predictors/additional_variables.csv")
#write.table(dat_2, file="/Users/ti1/Google\ Drive/study2_snps_as\ predictors/Psychiatric diagnoses in SELCoH.csv")
#write.table(dat_3, file="/Users/ti1/Google\ Drive/study2_snps_as\ predictors/SELCoH ancestry scores.csv")
#write.table(dat_4, file="/Users/ti1/Google\ Drive/study2_snps_as\ predictors/selcoh&akt1pathway.csv")


columns = colnames(dat_4)

delete = c("am1_intvra", "cd1_hhprsn", "cd1_countborn", "cd1_yrsuk", "cd1_eligible", "sd1_age", "sd1_genhealth", "sd1_agesmoke", "sd1_smokenow",
           "sd1_nosmoke", "sd1_stopsmoke", "sd1_yrcannab", "sd1_agecanna", "sd1_nocanna", "sd1_mthcanna", "sd1_frqcann",
           "c1_depscr", "c1_depcomm", "c1_depidea", "c1_phoblevl", "c1_worrylvl", "c1_anxscr", "c1_panicscr",
           "c1_cmplscr",  "c1_obsvthts",  "c1_social_impair",  "c1_total_score",  "c1_primdiag",  "c1_seconddiag", "c1_suiciderisk",
           "dvd1_hhcluster", "cm1_vretentionscrper", "phouse", "new_psq_mania_1", "new_psq_mania_2", "new_psq_mania_3", "new_psq_thought_1",
           "new_psq_thought_2", "new_psq_paranoia_1", "new_psq_paranoia_2", "new_psq_paranoia_3", "new_psq_strange_1", "new_psq_strange_2", "new_psq_seehear_1", "new_psq_seehear_2",
          "psq_level_1_man", "psq_level_2_man", "psq_level_3_man", "psq_level_1_no_man", "psq_level_2_no_man",
           "canuse_times", "canuse_times_2", "canuse_times_3", "contact_1_rels", "practsupport", "emotsupport", "non_threat_events_ever",
           "age", "edu", "alone", "single", "psq_sum", "emp", "relat_stat", "housing_ten", "emp2", "social_class", "social_class_2", "index_depriv", "overall", "threat_1", "threat_2", "threat_3",
           "threat_4", "threat_5", "threat_lstyr", "social_class_3", "country", "Dvd1_country", "country2", "Dvd1_country2", "contact_2_rels", "area_threat") 

cluster_identification = c("cd1_hhid")

features = c("ID", "cd1_sex", 
        "cd1_age", "sd1_ethnicity", "sd1_smoke", "sd1_cannab", "sd1_amphet", "sd1_cocaine", "sd1_ecstasy", "sd1_lsd", "sd1_tranq", "sd1_crack", "sd1_heroin",
        "dvd1_socialclass", "dvd1_edu", "pw1", "ethnic", "sexabuse", "anynonpx", "canuse", "sexphys", "sexphys_2", "total_events_ever", "total_events_ever_2", 
        "threat_events_ever", "threat_events_ever_2", "threat_events_ever_3", "sexorphys", "violence")

target = "psq_level_3_no_man"

snps = columns[!columns %in% c(target, features, cluster_identification, delete)]


target_data = dat_4[, c(target)]
features_data = dat_4[, c(features)]
snps_data = dat_4[, c(snps)]

complete_data = list(target=target_data, individual_factors=features_data, snps=snps_data)
