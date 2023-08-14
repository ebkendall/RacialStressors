# Format the data into workable conditions
formatting_data = function(data_time, t) {
    data_format = NULL
    colnums = NULL
    if(t == 15)  colnums = 6:26
    if(t == 30) colnums = 6:15
    for(i in 1:nrow(data_time)) {
        split_id = strsplit(data_time$StudySubjectTask[i], split="_", fixed = T)
        id = as.numeric(split_id[[1]][2])
        state = split_id[[1]][3]
        
        RSA = data_time[i, colnums]
        if(!is.na(RSA[1])) {
            RSA = RSA[!is.na(RSA)]
            
            temp = cbind(id, state, as.numeric(RSA))
            rownames(temp) = NULL
            colnames(temp) = NULL
            
            data_format = rbind(data_format, temp)
        } else {
            print(split_id[[1]])
            print("has no data")
        }
    }
    colnames(data_format) = c("id", "state", "rsa")
    
    # Re-order the timeline
    for(i in 1:length(unique(data_format[,'id']))) {
        id = unique(data_format[,'id'])[i]
        sub_dat = data_format[data_format[,'id'] == id, ]
        
        s1 = sub_dat[sub_dat[,'state'] == "Baseline", ,drop=F]
        s2 = sub_dat[sub_dat[,'state'] == "Stress", ,drop=F]
        s3 = sub_dat[sub_dat[,'state'] == "Recovery", ,drop=F]
        if(nrow(sub_dat) != nrow(s1) + nrow(s2) + nrow(s3)) print("error")
        
        reordered_dat = rbind(s1, s2, s3)
        
        data_format[data_format[,'id'] == id, ] = reordered_dat
    }
    data_format = as.data.frame(data_format)
    data_format[,'rsa'] = as.numeric(data_format[,'rsa'])
    
    return(data_format)
}

# three different time frames of data
data_30 = read.csv('Data/Raw IBI files_HRV_summary_30secepocs_032323.csv')
data_30 = data_30[-1,]
# colnames(data_30)
data_15  = read.csv('Data/Raw IBI files_HRV_summary_15secepochs_032323.csv')
data_15 = data_15[-1,]
# colnames(data_15)


data_format_30 = formatting_data(data_30, 30)
data_format_15  = formatting_data(data_15, 15)

# checking the mean rsa values 
print(paste0("30s --> S1: ", round(mean(data_format_30[data_format_30[,'state'] == "Baseline", 'rsa']), 3),
                    ", S2: ", round(mean(data_format_30[data_format_30[,'state'] == 'Stress', 'rsa']), 3),
                    ", S3: ", round(mean(data_format_30[data_format_30[,'state'] == 'Recovery', 'rsa']), 3)))
print(paste0("15s --> S1: ", round(mean(data_format_15[data_format_15[,'state'] == "Baseline", 'rsa']), 3),
             ", S2: ", round(mean(data_format_15[data_format_15[,'state'] == 'Stress', 'rsa']), 3),
             ", S3: ", round(mean(data_format_15[data_format_15[,'state'] == 'Recovery', 'rsa']), 3)))


# Chaning the state labels
data_format_15[data_format_15[,"state"] == "Baseline", "state"] = 1
data_format_15[data_format_15[,"state"] == "Stress", "state"] = 2
data_format_15[data_format_15[,"state"] == "Recovery", "state"] = 3
data_format_15 = cbind(data_format_15, 0)
colnames(data_format_15)[4] = "time"
for(i in unique(data_format_15[,"id"])) {
    sub_data = data_format_15[data_format_15[,'id'] == i, ]
    sub_data[,"time"] = 0:(nrow(sub_data)-1)
    
    data_format_15[data_format_15[,'id'] == i, ] = sub_data
}

data_format_30[data_format_30[,"state"] == "Baseline", "state"] = 1
data_format_30[data_format_30[,"state"] == "Stress", "state"] = 2
data_format_30[data_format_30[,"state"] == "Recovery", "state"] = 3
data_format_30 = cbind(data_format_30, 0)
colnames(data_format_30)[4] = "time"
for(i in unique(data_format_30[,"id"])) {
    sub_data = data_format_30[data_format_30[,'id'] == i, ]
    sub_data[,"time"] = 0:(nrow(sub_data)-1)
    
    data_format_30[data_format_30[,'id'] == i, ] = sub_data
}
colnames(data_format_15) = colnames(data_format_30) = c('ID..', 'State', 'RSA', 'Time')

save(data_format_30, file = 'Data/data_format_30.rda')
save(data_format_15, file = 'Data/data_format_15.rda')

# Adding Baseline covariates to the model
rsa_covariates = read.csv('Data/_FinalDataforRSASecondsStatesCovariates.csv', na.strings = "")
rsa_covariates = as.matrix(rsa_covariates)
load('Data/data_format_15.rda'); print(length(unique(data_format_15$ID..)))
load('Data/data_format_30.rda'); print(length(unique(data_format_30$ID..)))
cov_names = colnames(rsa_covariates)[8:16]

ids = unique(data_format_15$ID..)
cov_df = matrix(nrow = length(ids), ncol = length(cov_names)+1)
cov_df[,1] = ids
colnames(cov_df) = c("ID", cov_names)

for(i in 1:length(ids)) {
    id_index = which(rsa_covariates[,"ID.."] == ids[i])[1]
    temp = c(rsa_covariates[id_index, 8:16])
    cov_df[i, 2:ncol(cov_df)] = temp
}
save(cov_df, file = "Data/cov_df.rda")
# plot(data_format_1[data_format_1[,'id'] == 25897, 'rsa'])
# plot(data_format_5[data_format_5[,'id'] == 25897, 'rsa'])
# plot(data_format_30[data_format_30[,'id'] == 25897, 'rsa'])
# 
# load('Data/data_format_1.rda')
# load('Data/data_format_30.rda')
# load('Data/data_format_5.rda')
# mean(data_format_1[data_format_1[,"state"] == "Baseline", "rsa"], na.rm = T)
# mean(data_format_1[data_format_1[,"state"] == "Stress", "rsa"])
# mean(data_format_1[data_format_1[,"state"] == "Recovery", "rsa"])

# sd = 1.5
# mean baseline = 6.5


# Empirical estimates and information -----------------------------------------
# load('Data/data_format_15.rda')
# data_format = data_format_15
load('Data/data_format_30.rda')
data_format = data_format_30
# load('Data/sim_data_1_a.rda')
# data_format = sim_data
# data_format = as.data.frame(data_format)

unique_id = unique(data_format$ID..)

info_criteria = matrix(ncol = 6, nrow = length(unique_id))
colnames(info_criteria) = c("m_mu", "m_alpha", "m_beta", "sd_1", "sd_2", "sd_3")

for(i in 1:length(unique_id)) {
    sub_i = unique_id[i]
    data_sub = data_format[data_format$ID.. == sub_i, ]
    data_sub$State = as.numeric(data_sub$State)
    data_sub$RSA = as.numeric(data_sub$RSA)
    
    mean_s1 = mean(data_sub$RSA[data_sub$State == 1])
    mean_s2 = mean(data_sub$RSA[data_sub$State == 2])
    mean_s3 = mean(data_sub$RSA[data_sub$State == 3])
    
    sd_s1 = sd(data_sub$RSA[data_sub$State == 1])
    sd_s2 = sd(data_sub$RSA[data_sub$State == 2])
    sd_s3 = sd(data_sub$RSA[data_sub$State == 3])
    
    info_criteria[i, ] = c(mean_s1, mean_s1 - mean_s2, mean_s1 - mean_s3,
                           sd_s1, sd_s2, sd_s3)
    
}

colMeans(info_criteria, na.rm = T)

miss_info = NULL
miss_info = c(miss_info, unique_id[is.nan(info_criteria[,"m_alpha"]) | is.nan(info_criteria[,"m_beta"])])
miss_info = unique(miss_info)



