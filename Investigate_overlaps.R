setwd("C:/Users/spjtcoi/Google Drive/top_hits_trainingset_bmi_chip")
top200_training<-read.csv("remove_top_200.csv", header= F, stringsAsFactors = F)
top100_training<-read.csv("remove_top_100.csv", header= F, stringsAsFactors = F)
top50_training<-read.csv("remove_top_50.csv", header= F, stringsAsFactors = F)
top20_training<-read.csv("remove_top_20.csv", header= F, stringsAsFactors = F)
top10_training<-read.csv("remove_top_10.csv", header= F, stringsAsFactors = F)

setwd("C:/Users/spjtcoi/Google Drive/top_hits_validation_bmi_chip")
top200_test<-read.csv("remove_top_200.csv", header= F, stringsAsFactors = F)
top100_test<-read.csv("remove_top_100.csv", header= F, stringsAsFactors = F)
top50_test<-read.csv("remove_top_50.csv", header= F, stringsAsFactors = F)
top20_test<-read.csv("remove_top_20.csv", header= F, stringsAsFactors = F)
top10_test<-read.csv("remove_top_10.csv", header= F, stringsAsFactors = F)
#with chip
retained200 <- top200[which(top200$V1 %in% same_columns),]
retained10 <- top10[which(top10$V1 %in% same_columns),]
retained20<- top20[which(top20$V1 %in% same_columns),]
