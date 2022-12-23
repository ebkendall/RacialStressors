library(msm)

seedInd = 5
set.seed(seedInd)
print(seedInd)

new_data_gen <- function(d, dat_fr) {
    dat_fr$disc_time = dat_fr$Time - (dat_fr$Time %% d)
    
    newDat = NULL
    num = 0
    
    # Adding censored rows
    for(i in unique(dat_fr$ptnum)){
        
        current <- NULL
        subject <- dat_fr[dat_fr$ptnum==i,,drop=FALSE]
        
        #------------------------------------
        censoredAges <- seq(0, max(subject$Time), d)  
        yrs_round = round(subject$Time, digits = 5)
        disc_round = round(subject$disc_time, digits = 5)
        for(t in censoredAges){
            tempRow = NULL
            
            t_round = round(t, digits = 5)
            
            # If 't' corresponds to an observed age, then the next row will 
            # include the observed clinical visit data.
            if(t_round %in% yrs_round){	
                current <- rbind( current, subject[subject$disc_time==t_round,]) 
            } else{
                
                # Create a CENSORED row for each subject at each discrete unit of time.
                tempRow['ptnum'] <- i
                tempRow['Time'] <- t
                tempRow['State'] <- 99
                tempRow['RSA'] <- NA
                tempRow['obstrue'] <- 1  
                tempRow['disc_time'] <- t
                
                current <- rbind( current, tempRow)
                
                # If 't' corresponds to an observed discrete grid time, then the 
                # subject was observed some time during this time interval.  
                # Hence, the next row will include the observed clinical visit 
                # data.
                if(t_round %in% disc_round){ current <- rbind( current, subject[disc_round==t_round,]) }
            }
            
        }
        #------------------------------------
        
        newDat <- rbind( newDat, current)
        rownames(newDat) = NULL
        num <- num+1
    }
    
    dat_fr = newDat
    return(dat_fr)
}

load('Data/data_format.rda')
colnames(data_format)[1] = 'ptnum'

obstrue = 0
data_format = cbind(data_format, obstrue)

data_format = new_data_gen(1, data_format)

# Initializing the rate matrix
qmat <- matrix(c(       0,exp(-2),exp(-2),
                        0,      0,exp(-2),
                        0,      0,      0), ncol=3, byrow=TRUE)
dimnames(qmat) <- list( c('Well', 'Mild','Severe'), c('Well', 'Mild','Severe'))

#----------------------------------------------------------------------------------------------------------------
# Run the msm implementation ------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------------

# Misclassification response matrix on the logit scale
emat = matrix(c(    1, exp(-3), exp(-3),
       		  exp(-3),       1, exp(-3),
			  exp(-3), exp(-3),       1), ncol=3, byrow=TRUE)
emat = emat / rowSums(emat)
dimnames(emat) <- list( c('Well','Mild','Severe'), c('Well','Mild','Severe'))

Output_msm <- msm(State ~ Time, subject=ptnum, data=data_format, qmatrix=qmat, covariates= ~ 1 + disc_time, 
				  center=FALSE, covinits=list(disc_time=c(0,0,0)), obstrue=obstrue, 
				  ematrix=emat, initprobs=c(1, 1, 1), est.initprobs=TRUE, deathexact=NULL, 
				  censor=99, censor.states=1:3, method='BFGS', control=list(fnscale=4000, maxit=10000))   

save(Output_msm,file=paste0('Model_out/Output_msm',seedInd,'.rda'))