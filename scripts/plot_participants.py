import os
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
import matplotlib.patches as mpatches
import string
from collections import defaultdict
try:
	from scripts.utils import plot_style, import_config, numeric_from_datetime
except:
	from utils import plot_style, import_config, numeric_from_datetime


def plot_participant_vl(ax, copies_dat, config):
	# tabulate
	# N participants
	# N participants w/ HIV
	# N participants with LLV
	copies_dat_sum = \
		(lambda k: pd.DataFrame(k.tolist(), index=k.index, columns=['n_par', 'n_hiv']).reset_index())(
				copies_dat.groupby(['round', 'round_label', 'round_mid_date']).finalhiv.agg(lambda k: [k.shape[0], (k == "P").sum()])).\
			merge(
				copies_dat.groupby(['round']).copies.agg(lambda k: ((k > 200) & (k < 1000)).sum()).\
					reset_index(name='n_llv'),
				how='left',
				on='round')
	copies_dat_sum['round_mid_date'] = [numeric_from_datetime(i) for i in pd.to_datetime(copies_dat_sum['round_mid_date'])]
	tot_width = 0.75
	spacing = 0.01
	n_cat = 3
	bw = (tot_width - (n_cat-1)*spacing)/n_cat
	rel_loc = [-(bw*1.5+spacing/2)*tot_width, 0, (bw*1.5+spacing/2)*tot_width]
	_ = ax.bar(rel_loc[0] + copies_dat_sum.round_mid_date, copies_dat_sum['n_par'],
		facecolor='white', edgecolor='white', width=bw,
		linewidth=4, zorder=2)
	_ = ax.bar(rel_loc[1] + copies_dat_sum.round_mid_date, copies_dat_sum['n_hiv'],
		facecolor='white', edgecolor='white', width=bw,
		linewidth=4, zorder=2)
	_ = ax.bar(rel_loc[2] + copies_dat_sum.round_mid_date, copies_dat_sum['n_llv'],
		facecolor='white', edgecolor='white', width=bw,
		linewidth=4, zorder=2)
	_ = ax.bar(rel_loc[0] + copies_dat_sum.round_mid_date, copies_dat_sum['n_par'],
		label='all',
		facecolor=config['color_base'], edgecolor='#333333', width=bw, zorder=3)
	_ = ax.bar(rel_loc[1] + copies_dat_sum.round_mid_date, copies_dat_sum['n_hiv'],
		label='living with HIV',
		facecolor=config['color_hiv+'], edgecolor='#333333', width=bw, zorder=3)
	_ = ax.bar(rel_loc[2] + copies_dat_sum.round_mid_date, copies_dat_sum['n_llv'],
		label='with LLV', 
		facecolor=config['color_200_1000_copies'], edgecolor='#333333', width=bw, zorder=3)
	_ = ax.set_yscale('log')
	_ = ax.grid(axis='y', color='#eaeaea', zorder=1)
	_ = ax.legend(ncol=3, fontsize=13, handletextpad=0.25, columnspacing=1,
			frameon=True, edgecolor='#eaeaea')
	_ = ax.set_title('Participants', fontsize=18)
	#_ = ax.set_xticks(copies_dat_sum.round_mid_date, labels=copies_dat_sum.round_label, rotation=45, ha='right', va='top')
	_ = ax.set_xticks(copies_dat_sum.round_mid_date, labels=copies_dat_sum.round_label)
	_ = ax.set_ylim(1, 100000)
	_ = ax.set_yticks([1, 10, 100, 1000, 10000, 100000], labels=['1', '10', '100', '1,000', '10,000', '100,000'])
	_ = ax.tick_params(axis='both', labelsize=16)
	_ = ax.set_ylabel('count '+r'(log$_{10}$ scale)', fontsize=18)
	#_ = ax.set_xlabel('date')
	return(ax)

#fig, ax = plt.subplots()
#prev_dat = prev_dat.query('label == "llv_prev_among_viremic"')
'''
prev_dat = prev_dat.query('label == "prev_among_par"')[['round_mid_date', 'copies_cat', 'value']].\
			pivot(index='round_mid_date', columns='copies_cat', values='value').\
			reset_index().\
			assign(HIVNeg = lambda k: 1 - k.values[:,1:].sum(axis=1)).\
			melt(id_vars='round_mid_date', var_name='copies_cat').\
			merge(prev_dat[['copies_cat', 'copies_cat_label']].drop_duplicates(),
				how='left',
				on='copies_cat').fillna('hiv-')
'''
def plot_prob_vl(ax, prev_dat, config, truncate=False):
	round_labels = prev_dat[['round_mid_date', 'round_label']].drop_duplicates()
	# drop missing VL 
	prev_dat_wide = prev_dat[['round_mid_date', 'copies_cat_label', 'value']].\
		pivot(index='round_mid_date', columns='copies_cat_label', values='value')
	cats = np.array(['hiv-', '0_200_copies', '200_1000_copies', '1000_inf_copies'])
	prev_dat_wide = prev_dat_wide[[i for i in cats if i in prev_dat_wide.columns]]
	valid_cat_idx = np.array([idx for idx, i in enumerate(cats) if i in prev_dat_wide.columns])
	ax.stackplot(pd.to_datetime(prev_dat_wide.index), 100*prev_dat_wide.values.T.astype(float),
		edgecolor='#333333',
		colors=np.array([config[f'color_{i}'] for i in cats])[valid_cat_idx],
		labels=np.array([config[f'label_{i}'] for i in cats])[valid_cat_idx])
	ymin=0
	if truncate:
		ymin = 5*np.floor((100*prev_dat_wide.values[:,0].max() - 5)/5)
	ax.set_ylim(ymin, 100)
	ax.set_yticks(np.arange(ymin, 100+5, 5))
	ax.set_xlim(*pd.to_datetime(prev_dat_wide.index[[0,-1]]))
	ax.set_xticks(pd.to_datetime(prev_dat_wide.index), labels=round_labels['round_label'], size=16)
	#	rotation=45, ha='right', va='top')
	ax.legend(loc='upper left', reverse=True, handletextpad=0.5, fontsize=12, frameon=True,
		edgecolor='#eaeaea')
	#ax.set_xlabel('date')
	ax.set_ylabel('proportion (%)', fontsize=18)
	#fig.savefig('test.pdf')
	return(ax)


def plot_couple_periods_old(ax, period_dat, config):
	period_dat_sum = pd.concat(
		[i.groupby('round_mid_year').size() for i in period_dat],
		axis=1).rename(
			columns={
				0: 'n_couple', 1: 'n_monogamous',
				2:'n_monogamous_vl',
				3:'n_monogamous_vl_seroconversion'}).fillna(0.0)
	tot_width = 0.75
	spacing = 0.05
	n_cat = 4
	bw = (tot_width - (n_cat-1)*spacing)/n_cat
	rel_loc = [
		-(bw*1.5+1.5*spacing),
		-(bw*0.5 + spacing*0.5),
		(bw*0.5 + spacing*0.5),
		(bw*1.5+1.5*spacing)]
	_ = ax.bar(rel_loc[0] + period_dat_sum.index, period_dat_sum['n_couple'],
		facecolor='white', edgecolor='white', width=bw,
		linewidth=4, zorder=2)
	_ = ax.bar(rel_loc[1] + period_dat_sum.index, period_dat_sum['n_monogamous'],
		facecolor='white', edgecolor='white', width=bw,
		linewidth=4, zorder=2)
	_ = ax.bar(rel_loc[2] + period_dat_sum.index, period_dat_sum['n_monogamous_vl'],
		facecolor='white', edgecolor='white', width=bw,
		linewidth=4, zorder=2)
	_ = ax.bar(rel_loc[3] + period_dat_sum.index, period_dat_sum['n_monogamous_vl_seroconversion'],
		facecolor='white', edgecolor='white', width=bw,
		linewidth=4, zorder=2)
	_ = ax.bar(rel_loc[0] + period_dat_sum.index, period_dat_sum['n_couple'],
		label='all',
		facecolor=config['color_base'], edgecolor='#333333', width=bw, zorder=3)
	_ = ax.bar(rel_loc[1] + period_dat_sum.index, period_dat_sum['n_monogamous'],
		label='monogamous',
		facecolor=config['color_monogamous'], edgecolor='#333333', width=bw, zorder=3)
	_ = ax.bar(rel_loc[2] + period_dat_sum.index, period_dat_sum['n_monogamous_vl'],
		label='monogamous w/ VL',
		facecolor=config['color_vl'], edgecolor='#333333', width=bw, zorder=3)
	_ = ax.bar(rel_loc[3] + period_dat_sum.index, period_dat_sum['n_monogamous_vl_seroconversion'],
		label='monogamous w/ VL & seroconv.',
		facecolor=config['color_seroconversion'], edgecolor='#333333', width=bw, zorder=3)
	#for bc in b:
	#	bc._hatch_color=(0.8, 0.36, 0.36, 1.0)
	#	bc.stale=True
	#_ = ax.set_yscale('log')
	_ = ax.grid(axis='y', color='#eaeaea', zorder=1)
	_ = ax.legend(ncol=2, fontsize=12, handletextpad=0.25, columnspacing=1,
			frameon=True, edgecolor='#eaeaea', title_fontsize=12 )
	_ = ax.set_title('Serodifferent couples')
	_ = ax.set_xticks(period_dat_sum.index, labels=period_dat_sum.index, size=12)
	_ = ax.set_ylim(0, 50*np.ceil(period_dat_sum.values.max()/50)+50)
	#_ = ax.set_yticks([1, 10, 100, 1000, 10000, 100000], labels=['1', '10', '100', '1,000', '10,000', '100,000'])
	_ = ax.set_ylabel('count')
	#_ = ax.set_xlabel('date')
	return(ax)


def plot_periods(axs, period_dat, config):
	# for each couple get earliest lag date, that is y axis value
	period_dat = period_dat.merge(
		period_dat.groupby(['study_idDonor', 'study_idRecipient']).lag_int_dateRecipient.min().\
			reset_index().\
			sort_values(by='lag_int_dateRecipient').\
			assign(y = lambda k: np.arange(k.shape[0])+1).\
			drop(['lag_int_dateRecipient'], axis=1),
		how='left',
		on=['study_idDonor', 'study_idRecipient']).\
		merge(
			period_dat.groupby(['study_idDonor', 'study_idRecipient']).conversion.any().reset_index(),
			how='left', on=['study_idDonor', 'study_idRecipient'])
	for idx, i in period_dat.iterrows():
		_ = axs[1*i.conversion_y].plot([i.lag_int_dateRecipient, i.int_dateRecipient],
			[i.y, i.y],
			color=config['color_seroconversion'] if i.conversion_y else '#333333',
			lw=0.75, zorder=2)
	xmin = np.floor(period_dat.lag_int_dateRecipient.min())
	xmax = np.ceil(period_dat.int_dateRecipient.max())
	_ = [ax.grid(axis='x', color='#eaeaea', zorder=1) for ax in axs]
	_ = [ax.set_xlim(xmin, xmax) for ax in axs]
	_ = [ax.set_ylim(1, period_dat.y.max()) for ax in axs]
	_ = [ax.set_yticks([1, *np.arange(100, period_dat.y.max(), 100), period_dat.y.max()]) for ax in axs]
	_ = [ax.set_xticks(np.arange(1995, 2025, 5), 
		labels = np.arange(1995, 2025, 5).astype(str),
		size=14, 
		#rotation=45,
		#ha='right', va='top'
		) for ax in axs]
	axs[0].set_ylabel('couple')
	axs[1].set_ylabel('couple')
	#axs[0].set_xlabel('date')
	#axs[1].set_xlabel('date')
	axs[0].set_title('Non-seroconverting\ncouples')
	axs[1].set_title('Seroconverting\ncouples')
	return(axs)


def plot_periods(ax, period_dat, config):
	# for each couple get earliest lag date, that is y axis value
	period_dat = period_dat.merge(
		period_dat.groupby(['study_idDonor', 'study_idRecipient']).lag_int_dateRecipient.min().\
			reset_index().\
			sort_values(by='lag_int_dateRecipient').\
			assign(y = lambda k: np.arange(k.shape[0])+1).\
			drop(['lag_int_dateRecipient'], axis=1),
		how='left',
		on=['study_idDonor', 'study_idRecipient']).\
		merge(
			period_dat.groupby(['study_idDonor', 'study_idRecipient']).conversion.any().reset_index(),
			how='left', on=['study_idDonor', 'study_idRecipient'])
	#for idx, i in period_dat.iterrows():
	#	_ = ax.plot([i.lag_int_dateRecipient, i.int_dateRecipient],
	#		[i.y, i.y],
	#		color='#7860c9' if i.conversion_y else '#adadad',
	#		lw=0.75, zorder=3+1*i.conversion_y)
	grey = '#999999'
	for idx, i in period_dat.iterrows():
		_ = ax.plot([i.lag_int_dateRecipient, i.int_dateRecipient],
			[i.y, i.y],
			color=grey,
			lw=0.75, zorder=2)
	ax.scatter(
		period_dat.query("conversion_x").int_dateRecipient,
		period_dat.query('conversion_x').y,
		edgecolor='white',
		lw=0.50,
		facecolor='white',
		s=50, zorder=3)
	ax.scatter(
		period_dat.query("conversion_x").int_dateRecipient,
		period_dat.query('conversion_x').y,
		edgecolor=grey,
		lw=0.25,
		facecolor=config['color_seroconversion'],
		s=50, zorder=4)
	ax.legend(handles=[Line2D([0], [0], marker='o', linestyle='None',
		label='partner seropositivity',
    	markerfacecolor=config['color_seroconversion'], 
    	markeredgecolor=grey,
    	markersize=10,
    	markeredgewidth=1)],
		fontsize=16, handletextpad=0.25, columnspacing=1,
		frameon=True, edgecolor='#eaeaea',
		title_fontsize=12)
	xmin = np.floor(period_dat.lag_int_dateRecipient.min())
	xmax = np.ceil(period_dat.int_dateRecipient.max())
	_ = ax.grid(axis='x', color='#eaeaea', zorder=1)
	_ = ax.set_xlim(xmin, xmax)
	_ = ax.set_ylim(1, period_dat.y.max())
	_ = ax.set_yticks([1, *np.arange(100, period_dat.y.max(), 100), period_dat.y.max()])
	_ = ax.set_xticks(np.arange(1995, 2025, 5), 
		labels = np.arange(1995, 2025, 5).astype(str),
		size=14, 
		#rotation=45,
		#ha='right', va='top'
		)
	_ = ax.tick_params(axis='both', labelsize=16)
	ax.set_ylabel('couple', fontsize=18)
	#ax.set_xlabel('date')
	return(ax)


def run():
	parser = argparse.ArgumentParser()
	parser.add_argument('--copiesDat', 
	    help='path to rccs viral load data')
	parser.add_argument('--copiesPrevDat',
		help='path to estimated prevalence of LLV')
	parser.add_argument('--periodDat',
		help='path to rccs couples data')
	parser.add_argument('--config',
		help='path to config csv', default='config/config.csv')
	args = parser.parse_args()
	#config = {i[0]:i[1] for idx, i in pd.read_csv(args.config, header=None).iterrows()}
	config = import_config(args.config)
	'''
	args.copiesDat = 'data/RCCSdata_R016_R020_vl.tsv'
	args.copiesPrevDat = 'output/rccs_r16_r20_prev_copies_cat.tsv'
	args.periodDat = 'output/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d.tsv'
	'''
	copies_dat = pd.read_csv(args.copiesDat, sep='\t').query('round >= 16')
	
	round_mid_date = pd.DataFrame({
		'round': [16, 17, 18, 19,20], 
		'round_label': ['2013-\n2015', '2015-\n2016', '2016-\n2018', '2018-\n2020', '2021-\n2023'],
		'round_mid_year': [2014, 2015, 2017, 2019, 2022],
		'round_mid_date': ['2014-04-07', '2015-11-03', '2017-07-17', '2019-06-20', '2022-03-04']}).\
		assign(round_mid_date = lambda k: pd.to_datetime(k.round_mid_date))
	
	copies_dat = copies_dat.merge(round_mid_date, how='left', on='round')
	#.\
	#	assign(int_date = lambda k: pd.to_datetime(k.int_date))
	#copies_dat = copies_dat.merge(
	#	copies_dat.groupby('round').int_date.agg(lambda k: k.median()).\
	#		reset_index(name='round_mid_year').\
	#		assign(round_mid_year = lambda k: k.round_mid_year.dt.year),
	#	how='left', on='round')
	#round_mid_year = pd.read_csv('data/rccs_r1_r20_cohort_dat.tsv',
	#	sep='\t').\
	#	assign(int_date = lambda k: pd.to_datetime(k.int_date)).\
	#	groupby('round').int_date.agg(lambda k: k.median()).\
	#		reset_index(name='round_mid_date').\
	#		assign(round_mid_year = lambda k: k.round_mid_date.dt.year)
	#print(round_mid_year)
	

	prev_dat = pd.read_csv(args.copiesPrevDat, sep='\t').\
		merge(round_mid_date, how='left', on='round')

	# read in period dat, get just consecutive rounds
	period_dat = pd.read_csv(args.periodDat, sep='\t').\
		rename(columns={'round_x': 'round'}).\
		merge(round_mid_date, how='left', on='round')


	plot_style()


	fig = plt.figure(figsize=(6.4*2.5, 4.8*2.2), constrained_layout=True)
	#gs = GridSpec(10,8, figure=fig)
	#axs = [fig.add_subplot(gs[:4,:4]), fig.add_subplot(gs[1:4,4:6]), fig.add_subplot(gs[1:4,6:]), 
	#	fig.add_subplot(gs[4:,:])]
	#lg_ax = fig.add_subplot(gs[0,4:8])
	gs = GridSpec(12,12, figure=fig)
	axs = [fig.add_subplot(gs[:5,:4]), fig.add_subplot(gs[1:5,4:8]), fig.add_subplot(gs[1:5,8:]), 
		fig.add_subplot(gs[5:,:])]
	lg_ax = fig.add_subplot(gs[0,4:])
	axs[0] = plot_participant_vl(axs[0], copies_dat, config)
	axs[2] = plot_prob_vl(axs[2], 
		prev_dat.query('label == "prev_among_viremic"')\
			[['round_mid_date', 'round_label', 'copies_cat', 'value', 'copies_cat_label']],
			config)
	# todo not manually
	axs[2].set_ylim(0, 100)
	axs[2].set_yticks([0, 25, 50, 75, 100])
	axs[1] = plot_prob_vl(axs[1], 
		prev_dat.query('label == "prev_among_par"')[['round_mid_date', 'copies_cat', 'value']].\
			pivot(index='round_mid_date', columns='copies_cat', values='value').\
			reset_index().\
			assign(HIVNeg = lambda k: 1 - k.values[:,1:].sum(axis=1)).\
			melt(id_vars='round_mid_date', var_name='copies_cat').\
			merge(prev_dat[['copies_cat', 'copies_cat_label']].drop_duplicates(),
				how='left',
				on='copies_cat').\
			assign(copies_cat_label = lambda k: np.where(k.copies_cat_label.isnull(), 'hiv-', k.copies_cat_label)).\
			merge(round_mid_date, how='left', on='round_mid_date'),
		config,
		truncate=True)
	axs[2].set_title('Participants with HIV\n>200 copies/mL', fontsize=18)
	axs[1].set_title('Participants', fontsize=18)	
	handles, labels = axs[1].get_legend_handles_labels()
	lg_ax.legend(handles, labels, ncol=4, loc='center',
		fontsize=16, handletextpad=0.25, columnspacing=1,
		frameon=True, edgecolor='#eaeaea',
		title_fontsize=12,
		handlelength=2)
	lg_ax.axis('off')
	for ax in [axs[1], axs[2]]:
		ax.get_legend().remove()


	axs[3] = plot_periods(axs[3], period_dat, config)
	in_ax = axs[3].inset_axes(
		[0.00, 0.5, 0.39, 0.39],
		transform=axs[3].transAxes, zorder=10)
	img = plt.imread('figures/model_fig.jpg')
	in_ax.imshow(img, interpolation='none', zorder=1)
	in_ax.set_xlim(-img.shape[1]*0.025, img.shape[1]*1.025)
	in_ax.set_ylim(img.shape[0]*1.025, -img.shape[0]*0.025)
	#in_ax.axis('off')
	in_ax.set_xticks([])
	in_ax.set_yticks([])
	_ = [in_ax.spines[i].set_color('#eaeaea') for i in ['left', 'right', 'top', 'bottom']]
	axs[3].set_title('\nSerodifferent monogamous-partner couples', fontsize=18)
	#axs[3] = plot_couple_periods(axs[3], period_dat, config)
	#x_adj = [-0.275, -0.5, -0.5, -0.4, -0.4]
	x_adj = [-0.35, -0.22, -0.22, -0.07]
	for adx,ax in enumerate(axs):
		ax.text(x_adj[adx], 1.05, list(string.ascii_uppercase)[adx], transform=ax.transAxes, 
			fontsize=12, fontweight='bold', va='top', ha='left')




	#fig.tight_layout()
	os.makedirs("figures/pdf", exist_ok=True)
	fig.savefig('figures/pdf/study_participants.pdf')
	os.makedirs("figures/eps", exist_ok=True)
	fig.savefig('figures/eps/study_participants.eps')
	plt.close()




if __name__ == '__main__':
    run()

