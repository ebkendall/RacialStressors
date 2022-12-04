rsa_data = read.csv('Data/_FinalDataforRSASecondsStatesCovariates.csv', na.strings = "")

data_format = rsa_data[, c("ID..", "Time", "State", "RSA")]

na_patients = data_format[is.na(data_format$State), ]
na_patients = unique(na_patients$ID..)

data_format = data_format[!(data_format$ID.. %in% na_patients), ]

length(unique(data_format$ID..))

data_format$Time = data_format$Time / 1000

data_format$State[data_format$State == "Baseline"] = 1
data_format$State[data_format$State == "Stress"] = 2
data_format$State[data_format$State == "Recovery"] = 3
data_format$State[data_format$State == "recovery"] = 3

data_format$State = as.numeric(data_format$State)
data_format$Time = as.numeric(data_format$Time)
data_format$RSA = as.numeric(data_format$RSA)

save(data_format, file = 'Data/data_format.rda')


# mean_RSA = matrix(0, nrow = length(unique(data_format$ID..)), ncol = 3)
# 
# for(i in unique(data_format$ID..)) {
#     sub_dat = data_format[data_format$ID.. == i, ]
# 
#     mean_RSA[which(unique(data_format$ID..) == i), 1] = mean(sub_dat$RSA[sub_dat$State == 1])
#     mean_RSA[which(unique(data_format$ID..) == i), 2] = mean(sub_dat$RSA[sub_dat$State == 2])
#     mean_RSA[which(unique(data_format$ID..) == i), 3] = mean(sub_dat$RSA[sub_dat$State == 3])
#     # mean_RSA[which(unique(data_format$ID..) == i), 4] = i
# }
# colMeans(mean_RSA, na.rm = T)
# 
# temp = t(apply(mean_RSA, 1, diff))
# colMeans(temp, na.rm = T)
# 
# hist(mean_RSA[,2], breaks = sqrt(nrow(mean_RSA)), col = 2)
# hist(mean_RSA[,1], breaks = sqrt(nrow(mean_RSA)), col = 3, add = T)
# hist(mean_RSA[,3], breaks = sqrt(nrow(mean_RSA)), col = 4, add = T)