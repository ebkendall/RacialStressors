# Format the data into workable conditions
formatting_data = function(data_time, t) {
    data_format = NULL
    colnums = NULL
    if(t == 1)  colnums = 6:326
    if(t == 5)  colnums = 6:69
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
        print(nrow(sub_dat) == nrow(s1) + nrow(s2) + nrow(s3))
        
        reordered_dat = rbind(s1, s2, s3)
        
        data_format[data_format[,'id'] == id, ] = reordered_dat
    }
    data_format = as.data.frame(data_format)
    data_format[,'rsa'] = as.numeric(data_format[,'rsa'])
    
    return(data_format)
}

# three different time frames of data
data_30 = read.csv('Data/Raw IBI files_HRV_summary_30secondepochs.csv')
colnames(data_30)
data_5  = read.csv('Data/Raw IBI files_HRV_summary_5secondepochs.csv')
colnames(data_5)
data_1  = read.csv('Data/Raw IBI files_HRV_1secondepochs.csv')
colnames(data_1)


data_format_30 = formatting_data(data_30, 30)
data_format_5  = formatting_data(data_5, 5)
data_format_1  = formatting_data(data_1, 1)

# checking the mean rsa values 
print(paste0("30s --> S1: ", round(mean(data_format_30[data_format_30[,'state'] == "Baseline", 'rsa']), 3),
                    ", S2: ", round(mean(data_format_30[data_format_30[,'state'] == 'Stress', 'rsa']), 3),
                    ", S3: ", round(mean(data_format_30[data_format_30[,'state'] == 'Recovery', 'rsa']), 3)))

print(paste0("5s ---> S1: ", round(mean(data_format_5[data_format_5[,'state'] == "Baseline", 'rsa']), 3),
             ", S2: ", round(mean(data_format_5[data_format_5[,'state'] == 'Stress', 'rsa']), 3),
             ", S3: ", round(mean(data_format_5[data_format_5[,'state'] == 'Recovery', 'rsa']), 3)))

print(paste0("1s ---> S1: ", round(mean(data_format_1[data_format_1[,'state'] == "Baseline", 'rsa']), 3),
             ", S2: ", round(mean(data_format_1[data_format_1[,'state'] == 'Stress', 'rsa']), 3),
             ", S3: ", round(mean(data_format_1[data_format_1[,'state'] == 'Recovery', 'rsa']), 3)))

save(data_format_30, file = 'Data/data_format_30.rda')
save(data_format_5 , file = 'Data/data_format_5.rda')
save(data_format_1 , file = 'Data/data_format_1.rda')
