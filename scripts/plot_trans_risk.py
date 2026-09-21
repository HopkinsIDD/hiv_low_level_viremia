import argparse
import os
import ast
import string
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.colors as mcolors
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cmdstanpy
import glob
try:
	from scripts.utils import plot_style, split_pad_col, import_config, get_ind_n_trans_bd, summarize_draws, get_trans_form
	from scripts.plot_vl_trans_risk import plot_trans_risk, get_donor_copies_dat
except:
	from utils import plot_style, split_pad_col, import_config, get_ind_n_trans_bd, summarize_draws, get_trans_form
	from plot_vl_trans_risk import plot_trans_risk, get_donor_copies_dat




def run():
	parser = argparse.ArgumentParser()
	parser.add_argument('--mergedPeriodDat',
		help='path to rccs contact period data merged with viral load data')
	parser.add_argument('--fit', 
	    help='path to stan fit')
	parser.add_argument('--config',
		help='path to config file',
		default='config/config.csv')
	parser.add_argument('--label', default='')
	args = parser.parse_args()
	#args.mergedPeriodDat = 'output/ve_constrained_sex/periods.tsv'
	#args.fit = 'output/ve_constrained'
	config = import_config(args.config)

	# format merged period dat
	'''
	mp_d = pd.read_csv(args.mergedPeriodDat, sep='\t').\
		assign(
			ever_conversion = lambda k: 
				k.groupby(['study_idDonor', 'study_idRecipient']).conversion.transform('any'),
			min_lag_int_dateRecipient = lambda k: 
				k.groupby(['study_idDonor', 'study_idRecipient']).lag_int_dateRecipient.transform('min'),
			neg_max_int_dateRecipient = lambda k: 
				-1 * k.groupby(['study_idDonor', 'study_idRecipient']).int_dateRecipient.transform('max'))
	# subset to just seroconversions and assign y value
	sc_mp_d = mp_d.\
		query('ever_conversion').\
		merge(
			mp_d.query('ever_conversion')[['study_idDonor', 'study_idRecipient', 'min_lag_int_dateRecipient', 'neg_max_int_dateRecipient']].\
				drop_duplicates().\
				sort_values(by=['min_lag_int_dateRecipient', 'neg_max_int_dateRecipient']).\
				assign(y = lambda k: np.arange(k.shape[0]))\
				[['study_idDonor', 'study_idRecipient', 'y']],
			how='left',
			on=['study_idDonor', 'study_idRecipient'])
	# get donor copies dat
	donor_copies_dat = get_donor_copies_dat(sc_mp_d)
	# rounded maximum viral load
	vmax = 0.25*np.ceil(donor_copies_dat.donor_copies.max()/0.25)
	'''
	#vmax = np.ceil(np.log10(np.nanmax(np.hstack(split_pad_col(pd.read_csv(args.mergedPeriodDat, sep='\t').\
	#	query('conversion')['donor_copies'])).astype(float))))
	vmax=7
	# read in draws for transmission risk plot
	#f = cmdstanpy.from_csv(
	#		glob.glob(args.fit + '/*.csv'))
	dr = pd.read_csv(args.fit + '/select_draws.tsv', sep='\t').assign(trans_form = get_trans_form(args.fit))
	# get transmission risk across VL range	
	if 'mean_use_param1' in dr.columns:
		dr['use_param1'] = dr.mean_use_param1
		dr['use_param2'] = dr.mean_use_param2
		dr['use_param3'] = dr.mean_use_param3
	
	#
	np.quantile(
		get_ind_n_trans_bd(dr, t=100, v=np.array(np.log10([200])), form=dr.trans_form.iloc[0]),
		[0.025, 0.5, 0.975])
	np.quantile(
		get_ind_n_trans_bd(dr, t=100, v=np.array(np.log10([50])), form=dr.trans_form.iloc[0]),
		[0.025, 0.5, 0.975])
	
	v = np.linspace(0, vmax, 500)
	r = get_ind_n_trans_bd(dr, t=100, v=v, form=dr.trans_form.iloc[0])
	# summarize
	r_sum = summarize_draws(r).\
		assign(v = v)
	plot_style()
	vl_ticks = np.log10([1, 200, 1000, 10000, 100000, 10**6, 10**7])
	vl_labels = [f'{int(10**k):,}' if 10**k < 1000000 else f'{int(10**k/1000000)}M'for k in vl_ticks]

	#config = {i['var']: i.value for idx, i in pd.read_csv('config/config.tsv', sep='\t').iterrows()}
	fig,ax = plt.subplots(figsize=(4.8,3.6), layout="constrained")
	ax = plot_trans_risk(ax, r_sum, config)
	ax.set_xticks(vl_ticks, labels=vl_labels, size=9)
	ax.set_ylim(0, 25)
	os.makedirs("figures/pdf", exist_ok=True)
	fig.savefig(f'figures/pdf/intra_couple_trans_risk_{args.label}.pdf')
	os.makedirs("figures/eps", exist_ok=True)
	fig.savefig(f'figures/eps/intra_couple_trans_risk_{args.label}.eps')
	plt.close()

	if 'mean_use_param1' in dr.columns:
		if 'sex' in args.fit:
			v = np.linspace(0,7,100)
			r1 = pd.DataFrame(np.quantile(
				100*108 * dr.sex1_use_param1.values * (1 - (1 - dr.sex1_use_param2.values)**(dr.sex1_use_param3.values*10**v[:,np.newaxis])),
				[0.025, 0.25, 0.5, 0.75, 0.975],
				axis=1).T,
				columns=[0.025, 0.25, 0.5, 0.75, 0.975]).\
				assign(v = v)
			r2 = pd.DataFrame(np.quantile(
				100*108 * dr.sex2_use_param1.values * (1 - (1 - dr.sex2_use_param2.values)**(dr.sex2_use_param3.values*10**v[:,np.newaxis])),
				[0.025, 0.25, 0.5, 0.75, 0.975],
				axis=1).T,
				columns=[0.025, 0.25, 0.5, 0.75, 0.975]).\
				assign(v = v)
			fig,ax = plt.subplots(figsize=(4.8,3.6), layout="constrained")
			ax.grid(axis='both', color='#eaeaea', zorder=1)
			ax.fill_between(
				r1.v,
				r1[0.025], #- 0.025*(r_sum[0.975] - r_sum[0.025]),
				r1[0.975], #+ 0.025*(r_sum[0.975] - r_sum[0.025]),
				color='white',
				linewidth=2, zorder=2)
			ax.fill_between(
				r2.v,
				r2[0.025], #- 0.025*(r_sum[0.975] - r_sum[0.025]),
				r2[0.975], #+ 0.025*(r_sum[0.975] - r_sum[0.025]),
				color='white',
				linewidth=2, zorder=2)
			ax.fill_between(r1.v, r1[0.025], r1[0.975], color='indianred', linewidth=0, zorder=3, alpha=0.25)
			ax.fill_between(r2.v, r2[0.025], r2[0.975], color='steelblue', linewidth=0, zorder=3, alpha=0.25)
			ax.fill_between(r1.v, r1[0.25], r1[0.75], color='indianred', linewidth=0, zorder=4, alpha=0.5)
			ax.fill_between(r2.v, r2[0.25], r2[0.75], color='steelblue', linewidth=0, zorder=4, alpha=0.5)
			ax.plot(r1.v, r1[0.5], color='indianred', zorder=5, label='female')
			ax.plot(r2.v, r2[0.5], color='steelblue', zorder=5, label='male')
			ax.legend(title='index partner sex', fontsize=12, title_fontsize=12)
			ymax = 25
			ax.set_ylim(0, ymax)
			ax.set_yticks(np.arange(0,ymax+5,5))
			ax.set_xlim(1, np.floor(r_sum.v.max()))
			#axs.set_xlabel(r'viral load (log$_{10}$ copies/mL)')
			ax.tick_params(axis='y', which='major', labelsize=10)
			ax.set_xlabel(r'viral load (copies/mL)', size=14)
			ax.set_ylabel('transmissions/\n100 person-years', size=14)
			_ = [ax.spines[i].set_visible(False) for i in ['top', 'right']]
			os.makedirs("figures/pdf", exist_ok=True)
			fig.savefig(f'figures/pdf/intra_couple_trans_risk_{args.label}_sex.pdf')
			os.makedirs("figures/eps", exist_ok=True)
			fig.savefig(f'figures/eps/intra_couple_trans_risk_{args.label}_sex.eps')
			plt.close()





if __name__ == '__main__':
    run()
