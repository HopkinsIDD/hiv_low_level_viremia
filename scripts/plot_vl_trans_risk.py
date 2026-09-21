import argparse
import os
import ast
import string
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import cmdstanpy
import glob
try:
	from scripts.utils import plot_style, split_col, get_ind_n_trans_bd, get_couple_prob_trans_bd, summarize_draws, import_config
except:
	from utils import plot_style, split_col, get_ind_n_trans_bd, get_couple_prob_trans_bd, summarize_draws, import_config


def get_donor_copies_dat(sc_mp_d):
	# convert list to array
	#donor_copies = format_list_col(sc_mp_d.donor_copies)
	n_col = [i.shape[0] for i in sc_mp_d.donor_copies]
	donor_copies = np.hstack(sc_mp_d.donor_copies)
	# convert to log10, replace BD values (0) with 1 log10 copies/mL and flatten
	donor_copies = np.log10(np.where(
		donor_copies == 0., 1, donor_copies))
	# get dates
	#donor_copies_dates = format_list_col(sc_mp_d.donor_copies_date).flatten()
	donor_copies_dates = np.hstack(sc_mp_d.donor_copies_date)
	# get y values
	donor_copies_y = np.repeat(sc_mp_d.y.values, n_col).flatten()
	# return
	out = pd.DataFrame([donor_copies, donor_copies_dates, donor_copies_y],
				index=['donor_copies', 'donor_copies_date', 'donor_copies_y']).T.\
			drop_duplicates().dropna().\
			assign(donor_copies_y = lambda k: k.donor_copies_y.astype(float),
				donor_copies_date = lambda k: k.donor_copies_date.astype(float))
	return(out)
	

#ax = axs[0]
def plot_couples_dat(ax, period_dat, donor_copies_dat, config=None):
	from matplotlib.lines import Line2D
	from matplotlib import collections  as mc
	# first plot coupling periods
	ax.add_collection(mc.LineCollection([[(i[1], i[0]), (i[2], i[0])] for i in period_dat.values],
		colors='#333333',
		linewidths=0.5,
		zorder=5))
	# then add horizontal lines to delineate each coupling period
	ax.add_collection(mc.LineCollection([[(i[1], i[0]-0.4), (i[1], i[0]+0.4)] for i in np.unique(np.vstack([
			period_dat.values[:,[0,1]],
			period_dat.values[:,[0,2]]]), axis=1)],
		colors='#333333',
		linewidth=0.5,
		zorder=5))
	cols = np.array([
		config['color_0_1000_copies'],
		config['color_1000_inf_copies']])
	labels = np.array([
		config['label_0_1000_copies'],
		config['label_1000_inf_copies']])
	gt_1000 = donor_copies_dat.query('donor_copies > 3')
	lt_1000 = donor_copies_dat.query('donor_copies < 3')
	im = ax.scatter(
			gt_1000.donor_copies_date.values,
			gt_1000.donor_copies_y.values,
			facecolor=cols[1],
			marker = '^',
			zorder=3,
			s=25,
			edgecolor='#333333',
			linewidth=0.25)
	im = ax.scatter(
			lt_1000.donor_copies_date.values,
			lt_1000.donor_copies_y.values,
			facecolor=cols[0],
			marker = 'v',
			zorder=3,
			s=25,
			edgecolor='#333333',
			linewidth=0.25)
	# add legend
	ax.legend(handles=
		[Line2D([0], [0], marker=np.array(['v', '^'])[idx], 
				linestyle='None', 
				markeredgecolor='#333333', label=labels[idx], 
				markerfacecolor=cols[idx], markersize=10) for 
			idx in range(cols.shape[0]-1, -1, -1)],
		loc='upper left',
		fontsize=10,
		handletextpad=0.25,
		title='index partner',
		title_fontsize=10,
		frameon=True,
		edgecolor='#eaeaea')
	ax.grid(axis='x',zorder=1,color='#eaeaea')
	ax.set_xlabel('date', size=14)
	ax.set_ylabel('couple', size=14)
	_ = [ax.spines[i].set_visible(False) for i in ['top', 'right']]
	ax.autoscale()
	ax.set_ylim(-donor_copies_dat.donor_copies_y.max()*0.01, donor_copies_dat.donor_copies_y.max()*1.01)
	ax.set_yticks([])
	ax.set_xticks([1995, 2000, 2005, 2010, 2015, 2020],
		labels=[1995, 2000, 2005, 2010, 2015, 2020],
		size=10)
	return(ax)


def plot_trans_risk(ax, r_sum, config, annotation=True):
	ax.grid(axis='both', color='#eaeaea', zorder=1)
	ax.fill_between(
		r_sum.v,
		r_sum[0.025], #- 0.025*(r_sum[0.975] - r_sum[0.025]),
		r_sum[0.975], #+ 0.025*(r_sum[0.975] - r_sum[0.025]),
		color='white',
		linewidth=2, zorder=2)
	ax.fill_between(r_sum.v, r_sum[0.025], r_sum[0.975], color=config['color_trans_rate_95'], linewidth=0, zorder=3)
	ax.fill_between(r_sum.v, r_sum[0.25], r_sum[0.75], color=config['color_trans_rate_50'], linewidth=0, zorder=4)
	ax.plot(r_sum.v, r_sum[0.5], color=config['color_trans_rate_mid'], zorder=5)
	trans_rate200 = r_sum.iloc[np.argmin(np.abs(np.log10(200) - r_sum.v)),:]
	ax.plot([0, np.log10(200)],
		[trans_rate200[0.5],
			trans_rate200[0.5]],
		color='#333333',
		ls='--',
		zorder=6)
	ax.plot([np.log10(200), np.log10(200)],
		[0, trans_rate200[0.5]],
		color='#333333',
		ls='--',
		zorder=6)
	trans_rate1000 = r_sum.iloc[np.argmin(np.abs(3 - r_sum.v)),:]
	ax.plot([0, 3],
		[trans_rate1000[0.5],
			trans_rate1000[0.5]],
		color='#333333',
		ls='--', zorder=6)
	ax.plot([3, 3],
		[0, trans_rate1000[0.5]],
		color='#333333',
		ls='--', zorder=6)
	if annotation == True:
		# text labels
		ax.text(0.05, trans_rate1000[0.5],
			f'{np.round(trans_rate1000[0.5],1)} ({np.round(trans_rate1000[0.025],1)}, {np.round(trans_rate1000[0.975],1)})',
			color='#707070',
			va='bottom',
			ha='left')
		ax.text(0.05, trans_rate200[0.5],
			f'{np.round(trans_rate200[0.5],1)} ({np.round(trans_rate200[0.025],1)}, {np.round(trans_rate200[0.975],1)})',
			color='#707070',
			va='bottom',
			ha='left')
	ymax = 5*np.ceil(r_sum[0.975].max()/5)
	ax.set_ylim(0, ymax)
	ax.set_yticks(np.arange(0,ymax+5,5))
	ax.set_xlim(1, np.floor(r_sum.v.max()))
	#axs.set_xlabel(r'viral load (log$_{10}$ copies/mL)')
	ax.tick_params(axis='y', which='major', labelsize=10)
	ax.set_xlabel(r'viral load (copies/mL)', size=14)
	ax.set_ylabel('transmissions/\n100 person-years', size=14)
	_ = [ax.spines[i].set_visible(False) for i in ['top', 'right']]
	return(ax)


def plot_rr(ax, rr_sum, config):
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
	ax.set_xlim(rr_sum.v.min() - 0.025*v_range, rr_sum.v.max() + 0.025*v_range)
	#axs.set_xlabel(r'viral load (log$_{10}$ copies/mL)')
	ax.set_xlabel(r'viral load (copies/mL)', size=14)
	ax.set_ylabel(f'rate ratio \n v. {int(np.round(rr_sum.ref.iloc[0]/1000)*1000):,}' + r' log$_{10}$ copies/mL', size=14)
	ax.tick_params(axis='y', which='major', labelsize=10)
	ax.set_xticks(np.log10([200, 400, 600, 800, 1000]),
		labels = ['200', '400', '600', '800', '1,000'])
	ax.set_ylim(0, 1)
	_ = [ax.spines[i].set_visible(False) for i in ['top', 'right']]
	return(ax)




def run():
	parser = argparse.ArgumentParser()
	parser.add_argument('--mergedPeriodDat',
		help='path to rccs contact period data merged with viral load data')
	parser.add_argument('--fit', 
	    help='path to stan fit')
	parser.add_argument('--config',
		help='path to config file',
		default='config/config.csv')
	args = parser.parse_args()
	#args.mergedPeriodDat = 'data/not_shared/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d.tsv'
	#args.fit = 'output/ve_constrained'
	#args.config = 'config/config.csv'
	config = import_config(args.config) 

	# format merged period dat
	mp_d = pd.read_csv(args.mergedPeriodDat, sep='\t')

	for col in ['dt', 'donor_copies', 'donor_copies_date']:
		mp_d[col] = split_col(mp_d[col])

	for col in ['donor_copies_round', 'donor_copies_obs_idx']:
		mp_d[col] = split_col(mp_d[col], int)

	
	mp_d = mp_d.\
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
	#vmax = 0.25*np.ceil(donor_copies_dat.donor_copies.max()/0.25)
	vmax=7
	# get coupling periods dat
	period_dat = sc_mp_d[['y', 'int_dateRecipient', 'lag_int_dateRecipient']]
	

	# read in draws for transmission risk plot
	dr = pd.read_csv(args.fit + '/select_draws.tsv', sep='\t')
	# get transmission risk across VL range	
	v = np.linspace(0, vmax, 500)
	r = get_ind_n_trans_bd(dr, t=100, v=v)
	# summarize
	r_sum = summarize_draws(r).\
		assign(v = v)
	# finally, calculate risk ratio of intra-to-extra partner transmission
	# ratio of probability of transmission occurring within a couple over one year
	# to probability of extra-transmission occuring over one year
	#rr_intra_extra_sum = summarize_draws(
	#		get_couple_prob_trans_bd(dr, t=1, v=v) / \
	#			(1 - np.exp(-dr['beta_extra'].values*1))).\
	#	assign(v = v)
	# finally, calculate rate ratio of transmission from 200-1,000 copies/mL compared to median set-point viral load
	llv = np.hstack([np.linspace(np.log10(200),3,1000),3])
	rr_spvl_llv_sum = summarize_draws(
		get_ind_n_trans_bd(dr, t=1, v=llv) / \
			get_ind_n_trans_bd(dr, t=1, v=np.array([dr['mu_0'].median()]))).\
		assign(v = llv,
			ref = 10**dr['mu_0'].median())

	#fig, axs = plt.subplots(1,2,figsize=(6.4*1.25, 4.8*1.25), constrained_layout=True)
	plot_style()
	vl_ticks = np.log10([1, 200, 1000, 10000, 100000, 10**6, 10**7])
	vl_labels = [f'{int(10**k):,}' if 10**k < 1000000 else f'{int(10**k/1000000)}M'for k in vl_ticks]

	#config = {i['var']: i.value for idx, i in pd.read_csv('config/config.tsv', sep='\t').iterrows()}
	fig = plt.figure(figsize=(6.4*1.25, 4.8*1.25), layout="constrained")
	gs = GridSpec(2, 2, figure=fig)
	axs = [fig.add_subplot(gs[:, 0]), fig.add_subplot(gs[0,1]), fig.add_subplot(gs[1,1])]
	axs[0] = plot_couples_dat(axs[0], period_dat, donor_copies_dat, config)
	axs[1] = plot_trans_risk(axs[1], r_sum, config)
	axs[1].set_xticks(vl_ticks, labels=vl_labels, size=9)
	axs[1].legend(handles=
		 [Line2D([0], [0], color=config['color_trans_rate_mid'], lw=1, label='median'),
		 	mpatches.Patch(color=config['color_trans_rate_50'], label='50% CI'),
		 	mpatches.Patch(color=config['color_trans_rate_95'], label='95% CI')],
		 fontsize=10,
		 loc='upper left')
	axs[2] = plot_rr(axs[2], rr_spvl_llv_sum, config)
	axs[2].legend(handles=
		 [Line2D([0], [0], color=config['color_rr_intra_extra_mid'], lw=1, label='median'),
		 	mpatches.Patch(color=config['color_rr_intra_extra_50'], label='50% CI'),
		 	mpatches.Patch(color=config['color_rr_intra_extra_95'], label='95% CI')],
		 fontsize=10,
		 loc='upper left')
	#axs[2].set_xlim(2,4)
	#axs[2].set_ylim(0, 
	#	np.ceil(rr_intra_extra_sum.iloc[np.argmin(np.abs(4 - rr_intra_extra_sum.v)),:][0.5]))
	#axs[2].set_xticks(np.log10([200, 500, 1000, 5000, 10000]),
	#	labels=[200, 500, 1000, 5000, '10,000'],
	#	size=9)
	#axs[2].axhline(1, ls='--', color='#333333', zorder=5)
	x_adj = [-0.1, -0.25, -0.25]
	for adx,ax in enumerate(axs):
		ax.text(x_adj[adx], 1.05, list(string.ascii_uppercase)[adx], transform=ax.transAxes, 
		            fontsize=10, fontweight='bold', va='top', ha='left')


	os.makedirs("figures/pdf", exist_ok=True)
	fig.savefig('figures/pdf/intra_couple_transmission.pdf')
	os.makedirs("figures/eps", exist_ok=True)
	fig.savefig('figures/eps/intra_couple_transmission.eps')
	plt.close()



if __name__ == '__main__':
    run()


