import sys
import argparse
import glob
import pandas as pd
import numpy as np
import cmdstanpy
try:
	from scripts.utils import get_ind_n_trans_bd, get_couple_prob_trans_bd, summarize_draws
except:
	from utils import get_ind_n_trans_bd, get_couple_prob_trans_bd, summarize_draws



def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--fit', 
	    help='path to model fit')
	args = parser.add_argument('--v',
		help='viral loads in copies/mL',
		nargs='+',
		type=float)
	args = parser.parse_args()
	#args.fit = 'output/pt_vl_fnr_tx_all_bdMM_constrainedprox_all4Constrained_RCCSdata_R001_R019_all_copies_donor-noacute_recipient-exclusivelymonogamous'
	#args.v = [50, 200, 1000]
	dr = cmdstanpy.from_csv(
			glob.glob(args.fit + '/*.csv')).draws_pd()
	# probability of intra-couple transmission
	prob_intra = summarize_draws(get_couple_prob_trans_bd(dr, t=1, v=np.log10(args.v))).\
		assign(v=args.v)
	print('#### probability of intra-couple transmission over 1 year ####', file=sys.stdout)
	print(prob_intra.to_string(), file=sys.stdout)
	# transmission rate
	trans_rate = summarize_draws(get_ind_n_trans_bd(dr, t=100, v=np.log10(args.v))).\
		assign(v=args.v)
	print('#### rate of transmission per 100 person-years (does not account for saturation) ####', file=sys.stdout)
	print(trans_rate.to_string(), file=sys.stdout)




if __name__ == '__main__':
    run()
