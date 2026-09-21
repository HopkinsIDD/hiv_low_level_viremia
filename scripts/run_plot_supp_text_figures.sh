#!/usr/bin/env bash

# viral load trajectories
python3 scripts/plot_select_vl_traj.py \
	--periodCopiesDat 'data/not_shared/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d.tsv' \
	--rccsDat 'data/not_shared/RCCSdata_R001_R019_VOIs_clean.tsv' \
	--filter 'lambda k: (((k.copies1 > 200) & (k.copies1 < 1000))|((k.copies2 > 200) & (k.copies2 < 1000))) & (k.conversion == True)'


# big trees
python3 scripts/plot_big_trees.py \
	--periodDat 'data/not_shared/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d.tsv' \
	--p24Tree 'output/genetic_linkage/p24_phylo/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d_p24_focal_nn.fasta.treefile' \
	--gp41Tree 'output/genetic_linkage/gp41_phylo/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d_gp41_focal_nn.fasta.treefile'


# risk ratio 
python3 scripts/plot_intra_extra_rr.py \
	--fit output/ve_constrained


# trans risk x sex
python3 scripts/plot_trans_risk.py \
		--mergedPeriodDat 'data/not_shared/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d.tsv' \
		--fit output/ve_constrained_sex \
		--label ve_constrained_sex
 

# viral loads
# using r1-r19 data here as relevant to transmission analysis
# couple data only available through r19
python3 scripts/plot_copies.py \
	--copiesDat 'data/not_shared/RCCSdata_R001_R019_vl.tsv' \
	&> output/RCCSdata_R001_R019_copies.log


# omega^cease trace
# plot trace
python3 scripts/plot_trace.py \
	--log output/est_rebound/RCCSdata_R001_R020_tx_vl_pairs_uniq_llv_mcmc.tsv \
	--label rccs_omega_fit
