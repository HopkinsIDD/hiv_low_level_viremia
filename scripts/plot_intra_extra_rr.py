import os
import argparse
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
import string
from matplotlib.gridspec import GridSpec,GridSpecFromSubplotSpec
try:
	from scripts.utils import plot_style, split_col, get_ind_n_trans_bd, get_couple_prob_trans_bd, summarize_draws, import_config
except:
	from utils import plot_style, split_col, get_ind_n_trans_bd, get_couple_prob_trans_bd, summarize_draws, import_config


def run():
	parser = argparse.ArgumentParser()
	parser.add_argument('--fit', 
	    help='path to stan fit')
	parser.add_argument('--config',
		help='path to config file',
		default='config/config.csv')
	args = parser.parse_args()

	config = import_config(args.config) 

	dr = pd.read_csv(args.fit + '/select_draws.tsv', sep='\t')
	# get transmission risk across VL range	
	vmax = 7
	v = np.linspace(0, 7, 500)
	r = get_ind_n_trans_bd(dr, t=100, v=v)
	# summarize
	r_sum = summarize_draws(r).\
		assign(v = v)

	rr_sum = summarize_draws(
			get_couple_prob_trans_bd(dr, t=1, v=v) / \
				(1 - np.exp(-dr['beta_extra'].values*1))).\
		assign(v = v)

	plot_style()
	vl_ticks = np.log10([1, 200, 1000, 10000, 100000, 10**6, 10**7])
	vl_labels = [f'{int(10**k):,}' if 10**k < 1000000 else f'{int(10**k/1000000)}M'for k in vl_ticks]

	#config = {i['var']: i.value for idx, i in pd.read_csv('config/config.tsv', sep='\t').iterrows()}
	fig, ax = plt.subplots(constrained_layout=True)
	v_range = rr_sum.v.max() - rr_sum.v.min()
	ax.grid(axis='both', color='#eaeaea', zorder=1)
	ax.fill_between(
		rr_sum.v,
		rr_sum[0.025] - 0.025*(rr_sum[0.975] - rr_sum[0.025]),
		rr_sum[0.975] + 0.025*(rr_sum[0.975] - rr_sum[0.025]),
		color='white',
		linewidth=0,
		zorder=3)
	ax.fill_between(rr_sum.v, rr_sum[0.025], rr_sum[0.975], color=config['color_rr_intra_extra_95'], linewidth=0, zorder=3)
	ax.fill_between(rr_sum.v, rr_sum[0.25], rr_sum[0.75], color=config['color_rr_intra_extra_50'],  linewidth=0, zorder=4)
	ax.plot(rr_sum.v, rr_sum[0.5], color=config['color_rr_intra_extra_mid'], zorder=5)
	ax.set_xlim(1, np.floor(rr_sum.v.max()))
	#axs.set_xlabel(r'viral load (log$_{10}$ copies/mL)')
	ax.set_xlabel(r'viral load (copies/mL)', size=14)
	ax.set_ylabel('risk ratio\n(intra-couple / extra-couple)', size=14)
	ax.tick_params(axis='y', which='major', labelsize=10)
	_ = [ax.spines[i].set_visible(False) for i in ['top', 'right']]
	ax.set_xlim(2,4)
	ax.set_ylim(0, 
		np.ceil(rr_sum.iloc[np.argmin(np.abs(4 - rr_sum.v)),:][0.5]))
	ax.set_xticks(np.log10([200, 500, 1000, 5000, 10000]),
		labels=[200, 500, 1000, 5000, '10,000'],
		size=9)
	ax.axhline(1, ls='--', color='#333333', zorder=5)
	os.makedirs("figures/pdf", exist_ok=True)
	fig.savefig('figures/pdf/intra_extra_rr.pdf')
	os.makedirs("figures/eps", exist_ok=True)
	fig.savefig('figures/eps/intra_extra_rr.eps')
	plt.close()


if __name__ == '__main__':
    run()


