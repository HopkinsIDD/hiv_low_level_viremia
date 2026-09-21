#!/usr/bin/env bash
mkdir -p output


#### ------------------ ####
#### DO GENETIC LINKAGE ####
#### ------------------ ####
# Flag: will re-run trees
bash scripts/run_do_linkage.sh


#### ------------- ####
#### CALC LLV PREV ####
#### ------------- ####
python3 -u scripts/calc_prev.py \
	--dat data/RCCSdata_R016_R020_vl.tsv \
	&> output/RCCSdata_R016_R020_est_prev_llv.log



#### ------------------ ####
#### FIT VL-TRANS MODEL ###
#### ------------------ ####
bash scripts/run_vl_trans_models.sh \
	data/not_shared/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_d.tsv \
	data/not_shared/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d.tsv \
	output/genetic_linkage/couples_linkage_sum_simple.tsv

# with mock data
python3 -u scripts/full_bayesian_model.py \
	--copiesDat data/mock/MOCK_copies.tsv \
	--periodDat data/mock/MOCK_periods.tsv \
	--stan stan/vl_trans.stan \
	--transVL prox \
	--transForm ve \
	--linkedPeriods data/mock/MOCK_linkage_sum_simple.tsv \
	--label MOCK \
	&> output/MOCK.log


#### ---------------------- ####
#### SIMULATED TRAJECTORIES ####
#### ---------------------- ####
bash scripts/run_simulate_trajectories.sh


#### ------------ ####
#### PLOT FIGURES ####
#### ------------ ####
bash scripts/run_plot_main_text_figures.sh

#### ------------ ####
#### PLOT SFIGURES ####
#### ------------ ####
bash scripts/run_plot_plot_supp_text_figures.sh

#### --------------- ####
#### CALC STATISTICS ####
#### --------------- ####
python3 scripts/calc_stats.py \
	--dat data/not_shared/RCCSdata_R001_R020_VOIs_clean.tsv \
	--llvPrevDat output/RCCSdata_R016_R020_est_prev_llv.tsv \
	--periodDat 'data/not_shared/donor-all_recipient-all.tsv' 'data/not_shared/donor-all_recipient-exclusivelymonogamous.tsv' \
	--periodDatLabels 'donor-all_recipient-all' 'donor-all_recipient-exclusivelymonogamous' \
	--mergedPeriodDat = 'data/not_shared/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d.tsv' \
	--mergedPeriodDatLabels 'merged_donor-noacute_recipient-exclusivelymonogamous' \
	--linkageDat output/couples_linkage_sum.tsv \
	--fitDat \
		'output/ve_constrained' \
		'output/ve_constrained_sex' \
		'output/ve_constrained_estC' \
		'output/ve_artificially_constrained' \
		'output/ve_unconstrained' \
		'output/ve_constrained_0_1' \
		'output/ve_constrained_39_40' \
		'output/ve_nosaturation' \
		'output/ve_pretreatment' \
		'output/sigmoid_constrained' \
	--fitDatLabels \
		'main_analysis' \
		'sex_predictor' \
		'estC' \
		'artificially_constrained' \
		'no_genetic_constraints' \
		'low_art_tn_assumption' \
		'high_art_tn_assumption' \
		'no_saturation' \
		'pre-treatment' \
		'sigmoid' \
	--fitCompare output/transmission_model_selection.tsv \
	--sampledRiskDat output/simulated_trajectory_risk.tsv \
	--fixedRiskDat output/fixed_simulated_trajectory_risk_sum.tsv \
	--sampledReboundRate output/rccs_r1_r20_vl_pairs_treated_uniq_llv_fit_mcmc_all.tsv \

/Applications/datatooltk/bin/datatooltk --csv output/statistics.csv --output output/statistics.dbtex --csv-encoding UTF-8
