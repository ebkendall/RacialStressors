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

# Adding covariates to the data ------------------------------------------------
rsa_covariates = read.csv('Data/SEEL_Covariates_updated.csv', na.strings = "")
rsa_covariates = as.matrix(rsa_covariates)

load('Data/data_format_15.rda'); print(length(unique(data_format_15$ID..)))
rsa_covariates_sub = rsa_covariates[rsa_covariates[,"ID"] %in% data_format_15[,"ID.."], ]
rsa_covariates_sub = rsa_covariates_sub[which(apply(is.na(rsa_covariates_sub), 1, sum) == 0), ]

final_rsa_cov = matrix(0, nrow = nrow(rsa_covariates_sub), ncol = 6)
colnames(final_rsa_cov) = c("ID", "Age", "sex1", "edu_yes", "DLER_avg", "DLER_dis")
final_rsa_cov[,"ID"] = rsa_covariates_sub[,"ID"]
final_rsa_cov[,"Age"] = as.numeric(rsa_covariates_sub[,"Age"])
final_rsa_cov[,"sex1"] = as.numeric(rsa_covariates_sub[,"Sex"] == 1)
final_rsa_cov[,"edu_yes"] = as.numeric(rsa_covariates_sub[,"R_PEdu"] == 1)
final_rsa_cov[,"DLER_avg"] = as.numeric(rsa_covariates_sub[,"ODLERavg"])
final_rsa_cov[,"DLER_dis"] = as.numeric(rsa_covariates_sub[,"DLER_Dis_NoDis"] == 1)

data_format_15 = data_format_15[data_format_15[,"ID.."] %in% final_rsa_cov[,"ID"], ]
print(length(unique(data_format_15$ID..)))
data_format_add_on = matrix(nrow = nrow(data_format_15), ncol = 5)
colnames(data_format_add_on) = c("Age", "sex1", "edu_yes", "DLER_avg", "DLER_dis")
for(i in 1:nrow(data_format_15)) {
    dat_i = c(final_rsa_cov[final_rsa_cov[,"ID"] == data_format_15[i,"ID.."], ])
    data_format_add_on[i,] = dat_i[2:6]
}
data_format_15 = cbind(data_format_15, data_format_add_on)
data_format_15[,"ID.."] = as.numeric(data_format_15[,"ID.."])
data_format_15[,"State"] = as.numeric(data_format_15[,"State"])
save(data_format_15, file = 'Data/data_format_15.rda')

load('Data/data_format_30.rda'); print(length(unique(data_format_30$ID..)))
data_format_30 = data_format_30[data_format_30[,"ID.."] %in% final_rsa_cov[,"ID"], ]
print(length(unique(data_format_30$ID..)))
data_format_add_on = matrix(nrow = nrow(data_format_30), ncol = 5)
colnames(data_format_add_on) = c("Age", "sex1", "edu_yes", "DLER_avg", "DLER_dis")
for(i in 1:nrow(data_format_30)) {
    dat_i = c(final_rsa_cov[final_rsa_cov[,"ID"] == data_format_30[i,"ID.."], ])
    data_format_add_on[i,] = dat_i[2:6]
}
data_format_30 = cbind(data_format_30, data_format_add_on)
data_format_30[,"ID.."] = as.numeric(data_format_30[,"ID.."])
data_format_30[,"State"] = as.numeric(data_format_30[,"State"])
save(data_format_30, file = 'Data/data_format_30.rda')