library(matrixStats)
library(plotrix)

# Information defining which approach to take ----------------------------------
trial_num = 1
simulation = F
thirty = F
use_labels = F
case_b = T
s1 = T
# ------------------------------------------------------------------------------

Dir = 'Model_out/'

file_name = NULL
if(simulation) {
    if(thirty) {
        if(case_b) {
            file_name = paste0("Model_out/B_chain_", trial_num, "_sim_30b.rda")
        } else {
            if(use_labels) {
                file_name = paste0("Model_out/B_chain_", trial_num, "_sim_30.rda")
            } else {
                file_name = paste0("Model_out/B_chain_", trial_num, "_sim_30_nl.rda")
            }   
        }
    } else {
        if(case_b) {
            file_name = paste0("Model_out/B_chain_", trial_num, "_sim_15b.rda")
        } else {
            if(use_labels) {
                file_name = paste0("Model_out/B_chain_", trial_num, "_sim_15.rda")
            } else {
                file_name = paste0("Model_out/B_chain_", trial_num, "_sim_15_nl.rda")   
            }   
        }
    }
} else {
    if(thirty) {
        if(case_b) {
            if(s1) {
                file_name = paste0("Model_out/B_chain_", trial_num, "_30b_s1.rda")
            } else {
                file_name = paste0("Model_out/B_chain_", trial_num, "_30b_old.rda")   
            }
        } else {
            if(use_labels) {
                file_name = paste0("Model_out/B_chain_", trial_num, "_30.rda")
            } else {
                if(s1) {
                    file_name = paste0("Model_out/B_chain_", trial_num, "_30_nl_s1.rda")
                } else {
                    file_name = paste0("Model_out/B_chain_", trial_num, "_30_nl.rda")   
                }
            }   
        }
    } else {
        if(case_b) {
            if(s1) {
                file_name = paste0("Model_out/B_chain_", trial_num, "_15b_s1.rda")
            } else {
                file_name = paste0("Model_out/B_chain_", trial_num, "_15b_old.rda")   
            }
        } else {
            if(use_labels) {
                file_name = paste0("Model_out/B_chain_", trial_num, "_15.rda")
            } else {
                if(s1) {
                    file_name = paste0("Model_out/B_chain_", trial_num, "_15_nl_s1.rda")   
                } else {
                    file_name = paste0("Model_out/B_chain_", trial_num, "_15_nl.rda")      
                }
            }   
        }
    }
}

load(file_name)

if(simulation) {
    # Simulation
    if(thirty) {
        load('Data/sim_data_1_30.rda')
        data_format = sim_data
    } else {
        load('Data/sim_data_1_15.rda')
        data_format = sim_data  
    } 
} else {
    # Real data analysis
    if(thirty) {
	    load('Data/data_format_30.rda')
    	data_format = data_format_30
    } else {
	    load('Data/data_format_15.rda')
    	data_format = data_format_15   
    }
    
    miss_info = c(26296, 29698, 30625, 401, 423, 419, 457)
    data_format = data_format[!(data_format[,"ID.."] %in% miss_info), ]
}

EIDs = unique(data_format[,"ID.."])

# New patients ---------------------------------------------------------------
if(simulation) {
    if(thirty) {
        if(case_b) {
            pdf_title = paste0('Plots/chart_plot_', trial_num, '_sim_30b.pdf')
        } else {
            if(use_labels) {
                pdf_title = paste0('Plots/chart_plot_', trial_num, '_sim_30.pdf')
            } else {
                pdf_title = paste0('Plots/chart_plot_', trial_num, '_sim_30_nl.pdf')   
            }   
        }
    } else {
        if(case_b) {
            pdf_title = paste0('Plots/chart_plot_', trial_num, '_sim_15b.pdf')
        } else {
            if(use_labels) {
                pdf_title = paste0('Plots/chart_plot_', trial_num, '_sim_15.pdf')
            } else {
                pdf_title = paste0('Plots/chart_plot_', trial_num, '_sim_15_nl.pdf')   
            }
        }
    }
} else {
	if(thirty) {
	    if(case_b) {
	        if(s1) {
	            pdf_title = paste0('Plots/chart_plot_', trial_num, '_30b_s1.pdf')
	        } else {
	            pdf_title = paste0('Plots/chart_plot_', trial_num, '_30b_old.pdf')   
	        }
	    } else {
	        if(use_labels) {
	            pdf_title = paste0('Plots/chart_plot_', trial_num, '_30.pdf')
	        } else {
	            if(s1) {
	                pdf_title = paste0('Plots/chart_plot_', trial_num, '_30_nl_s1.pdf')   
	            } else {
	                pdf_title = paste0('Plots/chart_plot_', trial_num, '_30_nl.pdf')      
	            }
	        }
	    }
	} else {
	    if(case_b) {
	        if(s1) {
	            pdf_title = paste0('Plots/chart_plot_', trial_num, '_15b_s1.pdf')
	        } else {
	            pdf_title = paste0('Plots/chart_plot_', trial_num, '_15b_old.pdf')   
	        }
	    } else {
	        if(use_labels) {
	            pdf_title = paste0('Plots/chart_plot_', trial_num, '_15.pdf')
	        } else {
	            if(s1) {
	                pdf_title = paste0('Plots/chart_plot_', trial_num, '_15_nl_s1.pdf')
	            } else {
	                pdf_title = paste0('Plots/chart_plot_', trial_num, '_15_nl.pdf')      
	            }
	        }   
	    }
	}
}

pdf(pdf_title)

panels = c(4, 1)
par(mfrow=panels, mar=c(2,4,2,4), bg='black', fg='green')
for(i in EIDs){
	print(which(EIDs == i))
	indices_i = (data_format[,'ID..']==i)
	sub_dat = data_format[data_format[,"ID.."] == i, ]
	n_i = sum(indices_i)

	t_grid = data_format[indices_i, "Time"]
	
	if(simulation){
	    b_i = as.numeric(data_format[ indices_i,"True_state"])
	    to_s1 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==1]
	    to_s2 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==2]
	    to_s3 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==3]
	}

	color_choice = c('deeppink', 'deeppink', 'deeppink')

	plot(t_grid, data_format[indices_i, "RSA"], xlab='time', ylab = 'RSA', 
		col.main='green', main = paste0('Participant: ', i))
	axis( side=1, at=t_grid, col.axis='green', labels=t_grid)
	axis( side=2, at=seq(min(data_format[indices_i, "RSA"]), max(data_format[indices_i, "RSA"])), col.axis='green')
	
	if(simulation){
	    abline( v=t_grid[to_s1], col='deeppink', lwd=2)
	    abline( v=t_grid[to_s2], col='deeppink', lwd=2)
	    abline( v=t_grid[to_s3], col='deeppink', lwd=2)
	    col_choice = c('deeppink', 'deeppink', 'deeppink')
	    abline( v= t_grid[1], col = col_choice[b_i[1]], lwd = 2)
	} else {
	    s = diff(as.numeric(data_format$State[indices_i]))
	    abline(v = t_grid[1], col = color_choice[as.numeric(sub_dat$State[1])], lwd = 2)
	    for(j in 1:sum(s!=0)){
	        t_switch = which(s!=0)[j]+1
	        abline(v = t_grid[t_switch], col = color_choice[as.numeric(sub_dat$State[t_switch])], lwd = 2)
	    }
	}

	b_chain_ind = 20000:40000
	barplot( rbind(   colMeans(B_chain[b_chain_ind, indices_i] == 1),
				colMeans(B_chain[b_chain_ind, indices_i] == 2),
				colMeans(B_chain[b_chain_ind, indices_i] == 3)), 
				col=c( 'dodgerblue', 'firebrick1', 'yellow2'), 
				xlab='time', xaxt='n', space=0, 
				col.main='green', border=NA) 
	grid( nx=NA, NULL, col='white')
	legend( 'topright', inset=c(0,-.28), xpd=T, horiz=T, bty='n', x.intersp=.75,
			legend=c( 'Nominal', 'Stress', 'Recovery'), pch=15, pt.cex=1.5, 
					col=c( 'dodgerblue', 'firebrick1', 'yellow2'))
	axis( side=1, at=t_grid, col.axis='green', labels=t_grid)
	axis( side=2, at=0:1, col.axis='green')
}
dev.off()

print(pdf_title)