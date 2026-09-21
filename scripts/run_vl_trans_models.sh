#!/usr/bin/env bash
#1=data/not_shared/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_d.tsv
#2=data/not_shared/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d.tsv
#3=output/genetic_linkage/couples_linkage_sum_simple.tsv

#### ------------------ ####
#### FIT VL-TRANS MODEL ###
#### ------------------ ####
# all transmission as a function of proximal VL with constrained intra-partner ratio
python3 -u scripts/full_bayesian_model.py \
	--copiesDat $1 \
	--periodDat $2 \
	--stan stan/vl_trans.stan \
	--linkedPeriods $3 \
	--label $4ve_constrained \
	&> output/$4ve_constrained.log
scp figures/ve_constrained_pairs.pdf figures/pdf/ve_constrained_pairs.pdf

# with sex as a predictor
python3 -u scripts/full_bayesian_model.py \
	--copiesDat $1 \
	--periodDat $2 \
	--stan stan/vl_trans_pred.stan \
	--linkedPeriods $3 \
	--label $4ve_constrained_sex \
	&> output/$4ve_constrained_sex.log

# same as main analysis but we try to est c
python3 -u scripts/full_bayesian_model.py \
	--copiesDat $1 \
	--periodDat $2 \
	--stan stan/vl_trans.stan \
	--linkedPeriods $3 \
	--label $4ve_constrained_estC \
	--c -1 \
	&> output/$4ve_constrained_estC.log
scp figures/ve_constrained_estC_pairs.pdf figures/pdf/ve_constrained_estC_pairs.pdf


# artificially setting only the 2 llv couples that
# are linked in both regions as linked
# and forcing the remainder to be unlinked
# NOT robust to changes in IDs
cat \
	<(awk -F'\t' 'NR==FNR{a[$1] = True; next}{if (!($1 in a)) print $0}' \
		output/genetic_linkage/artificial_llv_linkage_sum_simple.tsv \
		$3) \
	output/genetic_linkage/artificial_llv_linkage_sum_simple.tsv \
	> output/genetic_linkage/$4artificial_linkage_sum_simple.tsv

python3 -u scripts/full_bayesian_model.py \
	--copiesDat $1 \
	--periodDat $2 \
	--stan stan/vl_trans.stan \
	--linkedPeriods output/genetic_linkage/artificial_linkage_sum_simple.tsv \
	--label $4ve_artificially_constrained \
	&> output/$4ve_artificially_constrained.log

# without constraining by genetic linkage data
python3 -u scripts/full_bayesian_model.py \
	--copiesDat $1 \
	--periodDat $2 \
	--stan stan/vl_trans.stan \
	--label $4ve_unconstrained \
	&> output/$4ve_unconstrained.log

# assuming BD viral loads are all between 0 and 1
python3 -u scripts/full_bayesian_model.py \
	--copiesDat $1 \
	--periodDat $2 \
	--stan stan/vl_trans.stan \
	--linkedPeriods $3 \
	--vSuprMin 0 \
	--vSuprMax 1 \
	--label $4ve_constrained_0_1 \
	&> output/$4ve_constrained_0_1.log

# assuming BD viral loads are all between 39 and 40
python3 -u scripts/full_bayesian_model.py \
	--copiesDat $1 \
	--periodDat $2 \
	--stan stan/vl_trans.stan \
	--linkedPeriods $3 \
	--vSuprMin 39 \
	--vSuprMax 40 \
	--label $4ve_constrained_39_40 \
	&> output/$4ve_constrained_39_40.log
	
# assuming no curve saturation 
python3 -u scripts/full_bayesian_model.py \
	--copiesDat $1 \
	--periodDat $2 \
	--stan stan/vl_trans.stan \
	--linkedPeriods $3 \
	--allowSaturation 0 \
	--label $4ve_nosaturation \
	&> output/$4ve_nosaturation.log

# only fitting to pre-treatment data
python3 -u scripts/full_bayesian_model.py \
	--copiesDat $1 \
	--periodDat $2 \
	--stan stan/vl_trans.stan \
	--transVL prox_pt \
	--linkedPeriods $3 \
	--label $4ve_pretreatment \
	&> output/$4ve_pretreatment.log

# sigmoid functional form
python3 -u scripts/full_bayesian_model.py \
	--copiesDat $1 \
	--periodDat $2 \
	--stan stan/vl_trans.stan \
	--transVL prox \
	--transForm sigmoid \
	--linkedPeriods $3 \
	--label $4sigmoid_constrained \
	&> output/$4sigmoid_constrained.log


# Run model selection on select models
Rscript scripts/model_selection.R