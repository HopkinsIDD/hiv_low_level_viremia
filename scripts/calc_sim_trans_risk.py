import argparse
import pandas as pd
import numpy as np
import cmdstanpy
import glob
import matplotlib.pyplot as plt
try:
	from scripts.utils import plot_style, import_config, get_ind_n_trans_bd, summarize_draws
except:
	from utils import plot_style, import_config, get_ind_n_trans_bd, summarize_draws


def run():
	parser = argparse.ArgumentParser()
	parser.add_argument('--transFit',
		help='file with transmission risk mcmc samples')
	parser.add_argument('--adhFit',
		help='file with adherence fits')
	parser.add_argument('--out',
		help='output file path',
		default='simulated_trajectory_risk')
	args = parser.parse_args()
	#args.transFit = 'output/ve_constrained/*.csv'
	#args.adhFit = 'output/simulations/sampled/*_tstates.tsv.gz'
	#args.out = 'output/simulated_risk'
	#out_name = '.'.join('_'.join(args.transFit.replace('*', '').split('/')[1:]).split('.')[:-1]) + \
	#	'.'.join('_'.join(args.adhFit.replace('*', '').split('/')[1:]).split('.')[:-2])
	out_name = args.out
	# read in posterior draws
	#f = cmdstanpy.from_csv(glob.glob(args.transFit))
	#dr = f.draws_pd() 
	dr = pd.read_csv(args.transFit+'/select_draws.tsv', sep='\t')

	
	risks_backup = []
	rebounded_risks_backup = []
	for idx, i in enumerate(glob.glob(args.adhFit)[:-1]):
		i_sim_dat = pd.read_csv(i, compression='gzip', sep='\t')
		i_sim_dat = i_sim_dat.\
			merge(
				i_sim_dat.query('t == 365.0')[['idx', 'v']].\
					rename(columns={'v': 'v_365'}),
				how='left',
				on='idx')
		# each simulation is paired with one posterior draw
		# all risk
		risks_backup.append(
			i_sim_dat[['t', 'v']].assign(
				trans_risk = lambda k: get_ind_n_trans_bd(dr.iloc[[idx]], t=1/365, v=k.v.values)).\
				groupby('t')[['trans_risk']].\
				agg(lambda k: k.sum()).cumsum().\
				reset_index().\
				assign(idx = idx))
		# risk contributed by those with VL > 1,000 copies/mL at 1yr
		rebounded_risks_backup.append(
			i_sim_dat.query('(v_365 > 3)&(t<=365)')[['t', 'v']].assign(
				trans_risk = lambda k: get_ind_n_trans_bd(dr.iloc[[idx]], t=1/365, v=k.v.values)).\
				groupby('t')[['trans_risk']].\
				agg(lambda k: k.sum()).cumsum().\
				reset_index().\
				assign(
					idx = idx,
					p_rebound = (i_sim_dat.v_365 > 3).sum()/i_sim_dat.shape[0]))
		# add rel risk
		rebounded_risks_backup[-1]['rel_trans_risk'] = \
			rebounded_risks_backup[-1].trans_risk / \
				risks_backup[-1].trans_risk
		


	cum_risk_sum = pd.concat(risks_backup)
	rebounded_cum_risk_sum = pd.concat(rebounded_risks_backup)

	cum_risk_sum_hpd = cum_risk_sum.drop(columns=['idx']).groupby('t').\
		trans_risk.apply(lambda k: summarize_draws([k])).\
		reset_index()
	
	# save
	cum_risk_sum.to_csv(f'output/{out_name}.tsv', sep='\t', index=None)
	cum_risk_sum_hpd.to_csv(f'output/{out_name}_sum.tsv', sep='\t', index=None)
	rebounded_cum_risk_sum.to_csv(f'output/{out_name}_rebounded_1yr.tsv', sep='\t', index=None)


if __name__ == "__main__":
    run()
