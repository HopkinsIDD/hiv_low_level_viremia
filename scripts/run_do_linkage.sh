#!/usr/bin/env bash

# linkage based distance and get closely related sequences
python3 scripts/couples_genetic_linkage.py \
	--p24 data/not_shared/rccs_all_p24_aln_dg.fasta \
	--gp41 data/not_shared/rccs_all_gp41_aln_dg.fasta \
	--periodCopiesDat data/not_shared/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d.tsv \
	--allCouplesDat data/not_shared/couples_AN_20220904_clean.tsv.gz \
	&> output/genetic_linkage.log

# build trees of de-duplicated sequences
mkdir output/genetic_linakge/p24_phylo
mv \
	data/not_shared/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d_p24_focal_nn.fasta \
	output/genetic_linakge/p24_phylo/.
iqtree2 \
	-T AUTO \
	-o K03455.1 \
	-B 1000 \
	-s output/genetic_linakge/p24_phylo/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d_p24_focal_nn.fasta \
	-redo

mkdir output/genetic_linakge/gp41_phylo
mv \
	data/not_shared/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d_gp41_focal_nn.fasta \
	output/genetic_linakge/gp41_phylo/.
iqtree2 \
	-T AUTO \
	-o K03455.1 \
	-B 1000 \
	-s output/genetic_linakge/gp41_phylo/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d_gp41_focal_nn.fasta  \
	-redo

# do phylo linkage
python3 scripts/couples_phylo_linkage.py \
	--label p24 \
	--tree output/genetic_linakge/p24_phylo/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d_p24_focal_nn.fasta.treefile \
	--periodCopiesDat data/not_shared/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d.tsv \
	&> output/phylo_linkage_p24.log
python3 scripts/couples_phylo_linkage.py \
	--label gp41 \
	--tree output/genetic_linakge/gp41_phylo/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d_gp41_focal_nn.fasta.treefile\
	--periodCopiesDat data/not_shared/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d.tsv \
	&> output/phylo_linkage_gp41.log

# summarize linkage dat
cat \
	<(echo "study_idDonor\tstudy_idRecipient\tcopies1\tcopies2\tp24_j\tp24_d_min\tp24_genolinked\tgp41_j\tgp41_d_min\tgp41_genolinked\tp24_phylolinked\tgp41_phylolinked") \
	<(awk -F'\t' 'FNR==NR{a[$1"_"$2] = "true"; b[$1"_"$2] = $3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10; next}{a[$1] = "true"; c[$1] = $2"\t"$3}END{for (i in a) print i"\t"b[i]"\t"c[i]}' \
		<(awk -F'\t' 'FNR==NR{for(i=1; i<=NF; i++) cols[$i] = i;next}
				{print $cols["study_idDonor"]"\t"$cols["study_idRecipient"]"\t"$cols["copies1"]"\t"$cols["copies2"]"\t"$cols["p24_j"]"\t"$cols["p24_d_min"]"\t"$cols["p24_linked"]"\t"$cols["gp41_j"]"\t"$cols["gp41_d_min"]"\t"$cols["gp41_linked"]}' \
			<(head -n 1 output/genetic_linkage/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d_dist.tsv) \
			<(tail -n +2 output/genetic_linkage/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d_dist.tsv)) \
		<(awk -F'\t' 'FNR==NR{a[$1"_"$2] = "true"; b[$1"_"$2] = $3; next}{a[$1"_"$2] = "true"; c[$1"_"$2] = $3}END{for (i in a) print i"\t"b[i]"\t"c[i]}' \
			<(tail -n +2 output/genetic_linkage/p24_phylo/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d_p24_phylo_dist.tsv) \
			<(tail -n +2 output/genetic_linkage/gp41_phylo/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d_gp41_phylo_dist.tsv)) | \
		sed 's/_/\t/g') > output/genetic_linkage/couples_linkage_sum.tsv

# summarize it even further
cat \
	<(echo "study_idDonor\tstudy_idRecipient\tlinked") \
	<(awk -F'\t' '{print $1"\t"$2"\t"$7"\t"$10"\t"$11"\t"$12}' output/genetic_linkage/couples_linkage_sum.tsv | grep -E "True|False" | \
		awk -F'\t' '{if (((($3 == "True")&&($5=="True")))||((($4 == "True")&&($6=="True")))){print$1"\t"$2"\tTrue"}else{print $1"\t"$2"\tFalse"}}') | \
	> output/genetic_linkage/couples_linkage_sum_simple.tsv
