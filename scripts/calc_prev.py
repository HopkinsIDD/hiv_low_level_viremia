import os
import pandas as pd
import numpy as np
import argparse
import cmdstanpy
import glob
from cmdstanpy import CmdStanModel
try:
	from scripts.utils import import_config
except:
	from utils import import_config


def run():
	parser = argparse.ArgumentParser()
	parser.add_argument('--dat', 
	    help='path to input data')
	parser.add_argument('--config',
		help='path to config csv', default='config/config.csv')
	args = parser.parse_args()
	config=import_config(args.config)
	#args.dat = 'data/RCCSdata_R016_R020_vl.tsv'
	#args.pred= ['round', 'round_mid_date']
	#args.strata = ['age_cat', 'sex', 'comm_type']
	#args.outcome = 'lambda k: np.array([pd.NA,"hiv-","bd_copies","40_200_copies","200_1000_copies","1000_inf_copies"])[np.where((k.finalhiv == "N")|k.finalhiv.isnull(),1,np.where((k.finalhiv == "P") & k.copies.isnull(),0,np.digitize(k.copies,[0, 40, 200, 1000])+1))]'
	#args.filter = 'lambda k: (k["round"] >= 16)'
	
	dat = pd.read_csv(args.dat, sep='\t').query('round >= 16')
	
	# replace missing finalhiv with string
	dat = dat.assign(finalhiv = lambda k: np.where(k.finalhiv.isnull(), "NA", k.finalhiv))
	# bin viral load data
	copies_cat_labels = np.array(['', 
		'0_200_copies',
		'200_1000_copies',
		'1000_inf_copies'])
	dat = dat.assign(copies_cat = lambda k: 
		np.where(k.copies.isnull(),
			-1,
			np.digitize(k.copies, config['breaks_copies'])))

	strata = ['age_cat', 'sex', 'comm_type']
	
	group_vars = ['round']+strata
	# group by survey round and strata
	# and get:
	# Number of participants (N)
	# Number of participants with HIV viral load data (N_finalhiv)
	# Number of participants with HIV seropositivity (N_finalhivP)
	# number of participants with VL data (N_copies_cat)
	# number of participants with VL in each category (N_copies_catX)
	# get an index for each strata
	strata_idx = dat.groupby(group_vars).size().\
		reset_index(name='n_strata').\
		assign(
			strata_idx = lambda k: np.arange(0,k.shape[0])+1,
			p_strata = lambda k: k.groupby(['round']).n_strata.transform(lambda g: g/g.sum()))

	dat_sum = dat.groupby(group_vars + ['finalhiv', 'copies_cat']).size().\
		reset_index(name='n').\
		merge(strata_idx)

	# add strata-level variables
	dat_sum = dat_sum.\
		merge(
			dat_sum.query('finalhiv != "NA"').groupby(group_vars).n.sum().reset_index(name='n_finalhiv_strata')).\
		merge(
			dat_sum.query('finalhiv == "P"').groupby(group_vars).n.sum().reset_index(name='n_finalhivP_strata')).\
		merge(
			dat_sum.query('copies_cat > 0').groupby(group_vars).n.sum().reset_index(name='n_copies_strata')).\
		merge(
			dat_sum.query('copies_cat >= 3').groupby(group_vars).n.sum().reset_index(name='n_viremic_strata')).\
		query('copies_cat > 0')
	
	dat_sum = dat_sum.merge(
		dat_sum.drop(['copies_cat', 'n'], axis=1).drop_duplicates().\
			reset_index().\
			assign(
				x = lambda k: 
					(k.n_viremic_strata / k.n_copies_strata) * 
						k.n_finalhivP_strata/k.n_finalhiv_strata *
						k.p_strata).\
			assign(p_strata_given_viremic = lambda k: k.groupby('round').x.transform(lambda g: g/g.sum()))\
			[group_vars+['p_strata_given_viremic']],
		how='left',
		on=group_vars)
			

	prev_among_par = dat_sum.assign(x = lambda k: 
		(k.n / k.n_copies_strata) * (k.n_finalhivP_strata/k.n_finalhiv_strata)*k.p_strata).\
		groupby(['round', 'copies_cat']).x.sum().\
		reset_index(name='value').\
		assign(label = 'prev_among_par')

	prev_among_viremic = dat_sum.query('copies_cat >= 2').\
		assign(x = lambda k: 
			(k.n / k.n_viremic_strata) * \
				k.p_strata_given_viremic).\
		groupby(['round', 'copies_cat']).x.sum().\
		reset_index(name='value').\
		assign(label = 'prev_among_viremic')

	# save point estimates
	dat_sum.to_csv('output/RCCSdata_R016_R020_strata_copies_sum.tsv', sep='\t', index=None)
	pd.concat([prev_among_par, prev_among_viremic]).\
		assign(copies_cat_label = lambda k: copies_cat_labels[k.copies_cat]).\
		to_csv('output/RCCSdata_R016_R020_prev_copies_cat.tsv', sep='\t', index=None)

	# add strata idx to individual level data
	dat = dat.\
		merge(dat_sum[['strata_idx', 'round']+strata].drop_duplicates(),
			how='left',
			on=['round']+strata)
	
	m = CmdStanModel(stan_file='stan/llv_prev.stan')
	# make output dir
	os.makedirs('output', exist_ok=True)
	
	# calculate prob HIV given finalhiv
	# calcualte prob viremic given copies_dat
	# calculate prob LLV given viremic given copies_dat
	est_probs = {}
	'''
	label, filter_str, y_str = ('p_finalhivP_given_finalhiv', 'lambda k: k.finalhiv != "NA"', 'lambda k: k.finalhiv == "P"')
	label, filter_str, y_str = ('p_viremic_given_copies_cat', 'lambda k: k.copies_cat > 0', 'lambda k: k.copies_cat >= 3')
	label, filter_str, y_str = ('p_llv_given_viremic_copies_cat', 'lambda k: k.copies_cat >= 3', 'lambda k: k.copies_cat == 3')
	'''
	for (label, filter_str, y_str) in \
			[('p_finalhivP_given_finalhiv', 'lambda k: k.finalhiv != "NA"', 'lambda k: k.finalhiv == "P"'),
				('p_viremic_given_copies_cat', 'lambda k: k.copies_cat > 0', 'lambda k: k.copies_cat >= 2'),
				('p_llv_given_viremic_copies_cat', 'lambda k: k.copies_cat >= 2', 'lambda k: k.copies_cat == 2')]:
		out_dir = f'output/{label}'
		os.makedirs(out_dir, exist_ok=True)	
		# removes old fits
		_ = [os.remove(i) for i in glob.glob(out_dir+'/*')]
		i_dat = dat[eval(filter_str)(dat)]
		# add individual IDs
		i_dat = i_dat.merge(
			i_dat[['study_id']].drop_duplicates().\
				assign(ind_idx = lambda k: np.arange(k.shape[0])+1),
			how='left',
			on='study_id')
		# prepare data
		stan_dat = {
			'N': i_dat.shape[0],
			'N_ind': i_dat.ind_idx.max(),
			'N_strata': i_dat.strata_idx.max(),
			'y': eval(y_str)(i_dat),
			'x1': i_dat.ind_idx,
			'x2': i_dat.strata_idx
		}
		# fit model
		f = m.sample(data = stan_dat, chains=4, iter_warmup=2500, iter_sampling=2500, show_console=True,
			output_dir=out_dir)
		# get draws
		dr = f.draws_pd()
		dr[['alpha', 'alpha_i_sd'] + [i for i in dr.columns if 'beta' in i]].to_csv(out_dir+'/select_draws.tsv', sep='\t', index=None)
		# dr = cmdstanpy.from_csv(glob.glob(out_dir + '/*.csv')).draws_pd()
		# get strata probs
		est_probs[label] = dr[[i for i in dr.columns if 'strata_prev' in i]].\
			assign(citer = lambda k: np.arange(k.shape[0])).\
			melt(id_vars = 'citer', value_name=label).\
			assign(strata_idx = lambda k: [int(j.split('[')[1].replace(']', '')) for j in k.variable])



	# calculate the prevalence of stil
	est_dat_sum = dat_sum.drop(['n', 'copies_cat'], axis=1).drop_duplicates()
	for key, value in est_probs.items():
		if key != 'label':
			est_dat_sum = est_dat_sum.merge(value.drop(['variable'], axis=1),
				how='left',
				on=['strata_idx'] if 'citer' not in est_dat_sum.columns else ['strata_idx', 'citer'])

	est_prev_among_par = est_dat_sum.\
		assign(x = lambda k: 
			k.p_llv_given_viremic_copies_cat * \
				k.p_viremic_given_copies_cat * \
				k.p_finalhivP_given_finalhiv * \
				k.p_strata).\
			groupby(['round', 'citer']).x.agg(lambda g: g.sum()).\
			reset_index(name='value').\
			groupby(['round']).value.quantile([0.025, 0.25, 0.5, 0.75, 0.975]).\
			reset_index().pivot(index='round', columns='level_1', values='value').\
			assign(label = 'prev_among_par')

	est_dat_sum = est_dat_sum.\
		reset_index().\
		assign(
			x = lambda k: 
				(k.p_viremic_given_copies_cat) * 
					k.p_finalhivP_given_finalhiv *
					k.n_strata).\
		assign(p_strata_given_viremic = lambda k: k.groupby(['round', 'citer']).x.transform(lambda g: (g/g.sum())))

	est_prev_among_viremic = est_dat_sum.\
			assign(x = lambda k: 
				(k.p_llv_given_viremic_copies_cat) * \
					k.p_strata_given_viremic).\
			groupby(['round', 'citer']).x.agg(lambda g: g.sum()).\
			reset_index(name='value').\
			groupby(['round']).value.quantile([0.025, 0.25, 0.5, 0.75, 0.975]).\
			reset_index().pivot(index='round', columns='level_1', values='value').\
			assign(label = 'prev_among_viremic')

	pd.concat([est_prev_among_par, est_prev_among_viremic]).\
		assign(copies_cat = 2, copies_cat_label = copies_cat_labels[2]).\
		reset_index().to_csv('output/RCCSdata_R016_R020_est_prev_llv.tsv', sep='\t', index=None)
		


if __name__ == '__main__':
    run()



