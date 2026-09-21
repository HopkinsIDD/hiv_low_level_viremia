#!/usr/bin/env bash

python3 scripts/plot_participants.py \
	--copiesDat 'data/RCCSdata_R016_R020_vl.tsv' \
	--copiesPrevDat 'output/RCCSdata_R016_R020_prev_copies_cat.tsv' \
	--periodDat 'data/not_shared/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d.tsv'

python3 scripts/plot_linkage_sum.py \
	--periodDat 'data/not_shared/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d.tsv' \
	--p24Tree 'output/genetic_linkage/p24_phylo/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d_p24_focal_nn.fasta.treefile' \
	--gp41Tree 'output/genetic_linkage/gp41_phylo/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d_gp41_focal_nn.fasta.treefile' \
	--linkageSum 'output/genetic_linkage/couples_linkage_sum.tsv'

python3 scripts/plot_vl_trans_risk.py \
	--mergedPeriodDat 'data/not_shared/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d.tsv' \
	--fit 'output/ve_constrained'

python3 scripts/plot_vl_traj.py \
	--llvDat 'data/RCCSdata_R001_R020_tx_vl_pairs_uniq_llv.tsv' \
	--suprDat 'data/RCCSdata_R001_R020_tx_vl_pairs_uniq_supr.tsv' \
	--simDat $(ls output/simulations/sampled/*_tstates.tsv.gz | grep "_1_") \
	--sampledRisk "output/simulated_trajectory_risk_sum.tsv" \
	--fixRisk "output/fixed_simulated_trajectory_risk_sum.tsv" \
	--sampledRate "output/est_rebound/RCCSdata_R001_R020_tx_vl_pairs_uniq_llv_mcmc_50_all.tsv"

