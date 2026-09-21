#!/usr/bin/env bash


#### ------------- ####
#### FIT ADHERENCE ####
#### ------------- ####
# first, fit lognormal distribution to reported quantiles in Treasure et al. 2016
lognorm_params=$(python3 scripts/fit_pre_post_art_distr.py)
# then, generate a fit params file based on the parameters estimated above
cat config/fit_params.csv > config/tmp_fit_params.csv
while read i; do
	label=$(echo $i | cut -f1 -d":")
	val=$(echo $i | cut -f2 -d":")
	sed -i '' "s/${label}_HERE/${val}/g" config/tmp_fit_params.csv
done <<(echo $lognorm_params)

# fit param file uses the parameters for NNRTI-based regimens as most individuals during sampling were on NNRTIs
python3 -u scripts/forward_simulate_fit_inference.py \
	--config config/fit_config.csv\
	--params config/tmp_fit_params.csv \
	&> output/RCCSdata_R001_R020_tx_vl_pairs_uniq_llv_fit.log



#### --------------------- ####
#### SIMUALTE TRAJECTORIES ####
#### --------------------- ####
# get last 50% of trace
tail \
	-n "$(echo "$(wc -l < "output/est_rebound/RCCSdata_R001_R020_tx_vl_pairs_uniq_llv_mcmc.tsv" | awk '{print $1}') / 2" | bc)" \
	output/est_rebound/RCCSdata_R001_R020_tx_vl_pairs_uniq_llv_mcmc.tsv | \
	awk -F'\t' '{print $NF}' \
	> output/est_rebound/RCCSdata_R001_R020_tx_vl_pairs_uniq_llv_mcmc_50_all.tsv

# randomly sample values from last 50% of trace
shuf output/est_rebound/RCCSdata_R001_R020_tx_vl_pairs_uniq_llv_mcmc_50_all.tsv | \
	head -n 1000 | \
	awk -F'\t' '{print $NF}' \
	> output/est_rebound/RCCSdata_R001_R020_tx_vl_pairs_uniq_llv_mcmc_50_sampled.tsv

# get sim_params with lognormal parameters
cat config/sim_params.csv > config/tmp_sim_params.csv
while read i; do
	label=$(echo $i | cut -f1 -d":")
	val=$(echo $i | cut -f2 -d":")
	sed -i '' "s/${label}_HERE/${val}/g" config/tmp_sim_params.csv
done <<(echo $lognorm_params)

# for each simulate 1000 trajectories
# first 10 with figure
mkdir -p output/simulations
rm -rf output/simulations/sampled
rm -rf figures/simulations/sampled
mkdir -p output/simulations/sampled
mkdir -p figures/simulations/sampled
idx=0
while read p; do
  	cat \
  		<(sed "s/HERE_OMEGA_IDX/${idx}_HERE_OMEGA/g" config/sampled_sim_config.csv | \
  			sed "s/HERE_OMEGA/${p}/g") \
  		<(echo "figure,'true'") \
		> config/tmp_sampled_sim_config.csv
	python3 scripts/forward_simulate_fit_inference.py \
		--config config/tmp_sampled_sim_config.csv \
		--params config/tmp_sim_params.csv
	idx=$((idx+1))
done < <(head -n 10 output/est_rebound/RCCSdata_R001_R020_tx_vl_pairs_uniq_llv_mcmc_50_sampled.tsv)

# remainder with no fig
while read p; do
  	cat \
  		<(sed "s/HERE_OMEGA_IDX/${idx}_HERE_OMEGA/g" config/sampled_sim_config.csv | \
  			sed "s/HERE_OMEGA/${p}/g") \
		> config/tmp_sampled_sim_config.csv
	python3 scripts/forward_simulate_fit_inference.py \
		--config config/tmp_sampled_sim_config.csv \
		--params config/tmp_sim_params.csv
	idx=$((idx+1))
done < <(tail -n 990 output/est_rebound/RCCSdata_R001_R020_tx_vl_pairs_uniq_llv_mcmc_50_sampled.tsv)

# finally, caclulate transmission risk from those simulated trajectories 
python3 scripts/calc_sim_trans_risk.py \
	--transFit 'output/ve_constrained' \
	--adhFit 'output/simulations/sampled/*_tstates.tsv.gz'

# simulated across fixed rebound rates
rm -rf output/simulations/fixed
mkdir -p output/simulations/fixed
# we want annual probability of cessation from 
# 25% to 60%: 25 30 35 40 45 50 55 60 
# p = -ln(1-P)/365 
for p in 0.0007881700615 0.00097719 0.00118023 0.00139952 0.00163791 0.00189903 0.00218769 0.00251039; do
	cat \
  		<(sed "s/HERE_OMEGA_IDX/HERE_OMEGA/g" config/fixed_sim_config.csv | \
  			sed "s/HERE_OMEGA/${p}/g") \
		> config/tmp_fixed_sim_config.csv
	python3 scripts/forward_simulate_fit_inference.py \
		--config config/tmp_fixed_sim_config.csv \
		--params config/tmp_sim_params.csv
done

# consolidate trans risk from these simulations
python3 scripts/consolidate_sim_trans_risk.py \
	--transFit 'output/ve_constrained/*.csv' \
	--adhSim 'output/simulations/fixed/*_tstates.tsv.gz' \
	--labelFmt 'lambda k: k.split("/")[-1].split("_")[1]'

