library(matrixStats)
library(plotrix)

args <- commandArgs(TRUE)
set.seed(args[1])

simulation = args[2]

trialNum = 6 # CHANGE EVERY TIME ******************

Dir = 'Model_out/'

# load(paste0( Dir, 'post_mcmc_out_dev',args[1],'_', trialNum, '.rda'))
load(paste0(Dir,'mcmc_out_',toString(args[1]),'_', trialNum,'.rda'))

if(simulation) {
    load('Data/Simulation/sim_data_1_c.rda')
    data_format = sim_data
} else {
    load('Data/data_format_15.rda')
    data_format = data_format_15
    # load('Data/data_format_30.rda')
    # data_format = data_format_30
}

EIDs = unique(data_format[,"ID.."])

# New patients ---------------------------------------------------------------
pdf(paste0('Plots/chart_plot_', trialNum, '_seed',toString(args[1]), '.pdf'))

panels = c(4, 1)
par(mfrow=panels, mar=c(2,4,2,4), bg='black', fg='green')
for(i in EIDs){
	print(which(EIDs == i))
	indices_i = (data_format[,'ID..']==i)
	sub_dat = data_format[data_format[,"ID.."] == i, ]
	n_i = sum(indices_i)

	t_grid = data_format[indices_i, "Time"]
	
	if(simulation){
	    b_i = data_format[ indices_i,'True_state']
	    to_s1 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==1]
	    to_s2 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==2]
	    to_s3 = (2:n_i)[diff(b_i)!=0 & b_i[-1]==3]
	}
	
	# bar_grid = seq( 0, n_i, by=5)[-1]

	color_choice = c('dodgerblue', 'firebrick1', 'yellow2')

	plot(t_grid, data_format[indices_i, "RSA"], xlab='time', ylab = 'RSA', 
		col.main='green', main = paste0('Participant: ', i)) #, " ", mean(data_format[indices_i, 'changed'])
	axis( side=1, at=t_grid, col.axis='green', labels=t_grid)
	axis( side=2, at=seq(min(data_format[indices_i, "RSA"]), max(data_format[indices_i, "RSA"])), col.axis='green')
	
	if(simulation){
	    abline( v=t_grid[to_s1], col='dodgerblue', lwd=2)
	    abline( v=t_grid[to_s2], col='firebrick1', lwd=2)
	    abline( v=t_grid[to_s3], col='yellow2', lwd=2)
	    col_choice = c('dodgerblue', 'firebrick1', 'yellow2')
	    abline( v= t_grid[1], col = col_choice[b_i[1]], lwd = 2)
	} else {
	    s = diff(as.numeric(data_format$State[indices_i]))
	    abline(v = t_grid[1], col = color_choice[as.numeric(sub_dat$State[1])], lwd = 2)
	    for(j in 1:sum(s)){
	        t_switch = which(s==1)[j]+1
	        abline(v = t_grid[t_switch], col = color_choice[as.numeric(sub_dat$State[t_switch])], lwd = 2)
	    }
	}

	barplot( rbind(   colMeans(mcmc_out$B_chain[, indices_i] == 1),
				colMeans(mcmc_out$B_chain[, indices_i] == 2),
				colMeans(mcmc_out$B_chain[, indices_i] == 3)), 
				col=c( 'dodgerblue', 'firebrick1', 'yellow2'), 
				xlab='time', xaxt='n', space=0, 
				col.main='green', border=NA) 
	grid( nx=NA, NULL, col='white')
	legend( 'topright', inset=c(0,-.28), xpd=T, horiz=T, bty='n', x.intersp=.75,
			legend=c( 'Nominal', 'Stress', 'Recovery'), pch=15, pt.cex=1.5, 
					col=c( 'dodgerblue', 'firebrick1', 'yellow2'))
	axis( side=1, at=t_grid, col.axis='green', labels=t_grid)
	axis( side=2, at=0:1, col.axis='green')
	
	# abline(v = t_grid[1], col = color_choice[sub_dat$State[1]], lwd = 2)
	# for(j in 1:sum(s)){
	# 	t_switch = which(s==1)[j]+1
	# 	abline(v = t_grid[t_switch], 
	# 		col = color_choice[sub_dat$State[t_switch]],
	# 		lwd = 2)
	# }
}
dev.off()
