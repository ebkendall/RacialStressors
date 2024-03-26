# This workflow file will run everything needed to reproduce the results in 
# the manuscript.
# Much of this workflow.sh is done in an embarrassingly parallel fashion. 

####################
# To simulate data
for c in {1..3}
do
Rscript simulate_data.r $c
done

# To run the MCMC code for the 100 simulated data sets, use the following code.
# Note that it is feasible to parallelize across data sets to improve computation
# time.
for seed in {1..100}
do
    for c in {1..3}
    do
    Rscript mcmc_runfile.r $seed $c
    done
done

# To obtain the maximum likelihood state sequence, run the following:
for c in {1..3}
do
Rscript state_sampler.r $c
done

# To visualize the trace plots for the MCMC, run the following:
for c in {1..3}
do
Rscript mcmc_outfile.r $c
done

# To produce the posterior probability plots, run the following:
for c in {1..3}
do
Rscript out_file_chart.r $c
done