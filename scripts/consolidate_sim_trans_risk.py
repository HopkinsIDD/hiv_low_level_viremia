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
	parser.add_argument('--adhSim',
		help='file with adherence simulations')
	parser.add_argument('--labelFmt',
		default = 'lambda k: k')
	parser.add_argument('--out',
		help='output file path',
		default='fixed_simulated_trajectory_risk')
	args = parser.parse_args()
	#args.transFit = 'output/pt_vl_fnr_tx_all_bdMM_constrainedprox_all4Constrained_RCCSdata_R001_R019_all_copies_donor-noacute_recipient-exclusivelymonogamous/*.csv'
	#args.adhSim = 'output/simulations/fixed/*_tstates.tsv.gz'
	#args.labelFmt = 'lambda k: k.split("/")[-1].split("_")[1]'
	out_name = args.out

	# read in posterior draws
	f = cmdstanpy.from_csv(glob.glob(args.transFit))
	dr = f.draws_pd()[['use_param1', 'use_param2', 'use_param3']]


	risks_backup = []
	dat_backup = []
	for idx, i in enumerate(glob.glob(args.adhSim)):
		dat_backup.append(pd.read_csv(i, compression='gzip', sep='\t'))
		risks_backup.append(
			dat_backup[-1].\
				assign(gdx = lambda k: np.floor(k.idx/100).astype(int)).\
				set_index('gdx').\
				assign(trans_per_dt = \
					lambda k: 
						k.groupby('gdx')['v'].\
							transform(lambda z: get_ind_n_trans_bd(dr.iloc[[z.index[0]]], t=1/365.25, v=z.values.flatten()).flatten())).\
				reset_index().\
				groupby(['gdx', 't'])[['trans_per_dt']].agg(lambda k: k.sum()).reset_index().\
				assign(
					label = lambda k: eval(args.labelFmt)(i),
					trans_risk = lambda k: k.groupby('gdx').trans_per_dt.cumsum()))


	'''
	output/simulations/fixed/fixed_0.0010_tstates.tsv.gz
	output/simulations/fixed/fixed_0.00150_tstates.tsv.gz
	output/simulations/fixed/fixed_0.00250_tstates.tsv.gz
	output/simulations/fixed/fixed_0.0020_tstates.tsv.gz
	output/simulations/fixed/fixed_0.0030_tstates.tsv.gz
	'''


	risk_sum = pd.concat(risks_backup)


	# get 95% HDs
	risk_sum_hpd = risk_sum.groupby(['label', 't']).\
			trans_risk.apply(lambda k: summarize_draws([k])).\
			reset_index()

	# save
	risk_sum.to_csv(f'output/{out_name}.tsv', sep='\t', index=None)
	risk_sum_hpd.to_csv(f'output/{out_name}_sum.tsv', sep='\t', index=None)


if __name__ == "__main__":
    run()



