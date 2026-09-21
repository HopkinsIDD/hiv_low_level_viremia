import os
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.collections import LineCollection
import string
from collections import defaultdict
try:
	from scripts.utils import plot_style, split_pad_col, numeric_from_datetime, datetime_from_numeric, import_config
except:
	from utils import plot_style, split_pad_col, numeric_from_datetime, datetime_from_numeric, import_config


def plot_couple_vl_traj(fig, gs0, couple, period_dat, rccs_dat, config):
	#import matplotlib.dates as mdates
	from matplotlib.ticker import MaxNLocator
	from matplotlib.gridspec import GridSpecFromSubplotSpec
	from matplotlib.lines import Line2D
	# get periods of couple
	couple_period_dat = period_dat.\
		query('(study_idDonor == @couple.study_idDonor)&(study_idRecipient == @couple.study_idRecipient)')
	# get VL data over range of coupling periods
	# flag: read in 4 as argument
	cols = np.array([
		config['color_0_200_copies'],
		config['color_200_1000_copies'],
		config['color_1000_inf_copies']])
	donor_vl_dat = rccs_dat.\
		query('~copies.isnull()').\
		query('study_id == @couple.study_idDonor').\
		query('numeric_int_date >= @couple_period_dat.lag_int_dateRecipient.min()-4').\
		query('numeric_int_date <= @couple_period_dat.int_dateRecipient.max()+4').\
		sort_values(by='round').\
		assign(copies_color = lambda k: cols[np.digitize(k.copies, config['breaks_copies'])-1]).\
		assign(copies = lambda k: np.where(k.copies < 200, np.where(k.lod >= 200, 20, 10), k.copies))
	recipient_hiv_dat = rccs_dat.\
		query('~finalhiv.isnull()').\
		query('study_id == @couple.study_idRecipient').\
		query('numeric_int_date >= @couple_period_dat.lag_int_dateRecipient.min()-0.25').\
		query('round <= first_finalhivP').\
		sort_values(by='round')
	print(donor_vl_dat[['study_id', 'copies', 'lod']])
	# define axes
	gs00 = GridSpecFromSubplotSpec(4,1, subplot_spec=gs0, wspace=0.0, hspace=0.0)
	ax1 = fig.add_subplot(gs00[0, :])
	ax2 = fig.add_subplot(gs00[1:, :])
	# plot periods of coupling
	for (col, lw, z) in [('white', 2, 3), ('#333333', 1, 4)]:
		ax1.scatter(
			recipient_hiv_dat.numeric_int_date,
			[1]*recipient_hiv_dat.shape[0],
			edgecolor=col,
			linewidth=lw,
			zorder=z+1,
			facecolor=recipient_hiv_dat.finalhiv.map({'N': '#eaeaea', 'P': config['color_seroconversion']}))
		for rdx, row in couple_period_dat[['int_dateRecipient', 'lag_int_dateRecipient']].iterrows():
			ax1.plot([row.int_dateRecipient, 
					row.lag_int_dateRecipient], [1,1],
				color=col,
				lw=lw,
				zorder=z)
	labels=['seronegative', 'seropositive']
	cols=['#eaeaea', config['color_seroconversion']]
	ax1.legend(handles=
			[Line2D([0], [0], marker='o', 
					color='#333333',
					markeredgecolor='#333333', label='monogamy', 
					markersize=0)] + [Line2D([0], [0], marker='o', 
					linestyle='None', 
					markeredgecolor='#333333', label=labels[idx], 
					markerfacecolor=cols[idx], markersize=10) for 
				idx in range(2)],
			loc='upper right',
			fontsize=10,
			handletextpad=0.25,
			title_fontsize=10,
			frameon=True,
			edgecolor='#eaeaea',
			ncol=3)
	ax1.tick_params(axis='both',length=0)
	ax1.set_yticks([1], labels=['partner'])
	# plot donor VL values
	for (col, lw, z) in [('white', 2, 3), ('#333333', 1, 4)]:
		if (donor_vl_dat.pre_treatment == True).any():
			ax2.plot(
				#donor_vl_dat.query('pre_treatment == True')['int_date'], 
				np.hstack([donor_vl_dat.query('pre_treatment == True')['numeric_int_date'].values,
						[donor_vl_dat.query('pre_treatment == False')['numeric_int_date'].iloc[0]]]) if \
					(donor_vl_dat.pre_treatment == False).any() else 
						donor_vl_dat.query('pre_treatment == True')['numeric_int_date'].values, 
				np.log10(np.maximum(
					np.hstack([donor_vl_dat.query('pre_treatment == True')['copies'].values,
						[donor_vl_dat.query('pre_treatment == False')['copies'].iloc[0]]]) if \
					(donor_vl_dat.pre_treatment == False).any() else 
						donor_vl_dat.query('pre_treatment == True')['copies'].values, 
					1)),
				color=col, lw=lw, zorder=z,
				label='pre-treatment' if col != 'white' else None)
		if (donor_vl_dat.pre_treatment == False).any():
			ax2.plot(
				#donor_vl_dat.query('pre_treatment == False')['int_date'], 
				donor_vl_dat.query('pre_treatment == False')['numeric_int_date'], 
				np.log10(np.maximum(donor_vl_dat.query('pre_treatment == False')['copies'].values, 1)),
				color=col, lw=lw, zorder=z, ls='--' if col != 'white' else '-',
				label='treatment-experienced' if col != 'white' else None)
		ax2.scatter(
			#donor_vl_dat['int_date'],
			donor_vl_dat['numeric_int_date'],
			np.log10(np.maximum(donor_vl_dat['copies'].values, 1)),
			linewidth=lw, zorder=z,
			facecolor=donor_vl_dat['copies_color'],
			edgecolor=col)
	ax2.legend(fontsize=10,
			handletextpad=0.25,
			title_fontsize=10,
			frameon=True,
			edgecolor='#eaeaea')
	#_  = [ax2.spines[i].set_visible(False) for i in [ 'top', 'right']]
	ax2.axhspan(np.log10(200), np.log10(1000), color='#ecdada', zorder=1)
	ax2.set_ylabel('\nindex viral load')
	vl_ticks = np.log10([10, 200, 1000, 10000, 100000, 10**6, 10**7])
	ax2.set_ylim(0, vl_ticks[-1])
	vl_labels = ['<200'] + [f'{np.round(10**k,0).astype(int):,}' if 10**k < 1000000 else f'{int(10**k/1000000)}M'for k in vl_ticks][1:]
	ax2.set_yticks(
			vl_ticks, labels=vl_labels, size=12)
	ax2.tick_params(axis='x', labelsize=12)
	ax2.set_xlabel('date\n')
	xmin = np.floor(min([donor_vl_dat.numeric_int_date.min(), recipient_hiv_dat.numeric_int_date.min()]))
	xmax = np.ceil(max([donor_vl_dat.numeric_int_date.max(), recipient_hiv_dat.numeric_int_date.max()]))
	_ = [ax.grid(axis='x', color='#eaeaea', zorder=2) for ax in [ax1, ax2]]
	#_ = [ax.set_xlim(
	#	pd.to_datetime(datetime_from_numeric(xmin - (xmax-xmin)*0.025)), 
	#	pd.to_datetime(datetime_from_numeric(xmax+(xmax-xmin)*0.025))) for ax in [ax1, ax2]]
	#_ = [ax.set_xlim(
	#	(xmin - (xmax-xmin)*0.025), 
	#	(xmax+(xmax-xmin)*0.025)) for ax in [ax1, ax2]]
	_ = [ax.set_xlim(xmin, xmax) for ax in [ax1, ax2]]
	if (xmax-xmin) <= 3:
		_ = [ax.xaxis.set_major_locator(MaxNLocator(integer=True)) for ax in [ax1, ax2]]
	#ax2.xaxis.set_major_formatter(mdates.DateFormatter('%b\n%Y'))
	ax1.set_xticklabels([])
	ax1.set_title(f'{couple.study_idDonorLabel} ' + r'$\rightarrow$' + f' {couple.study_idRecipientLabel}')
	return(ax1, ax2)



def run():
	parser = argparse.ArgumentParser()
	parser.add_argument('--periodCopiesDat',
		help='path to rccs couples data')
	parser.add_argument('--rccsDat',
		help='path to rccs VL data')
	parser.add_argument('--filter',
		default='lambda k: True',
		help='filter to select couples')
	parser.add_argument('--config',
		help='path to config csv', default='config/config.csv')
	args = parser.parse_args()
	config = import_config(args.config)
	#args.periodCopiesDat = 'data/not_shared/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d.tsv'
	#args.rccsDat = 'data/not_shared/RCCSdata_R001_R019_VOIs_clean.tsv' 
	#args.filter = 'lambda k: (((k.copies1 > 200) & (k.copies1 < 1000))|((k.copies2 > 200) & (k.copies2 < 1000))) & (k.conversion == True)'

	period_dat = pd.read_csv(args.periodCopiesDat, sep='\t')
	period_dat[['copies1', 'copies2']] = split_pad_col(period_dat.donor_copies).astype(float)

	plot_couples = period_dat[eval(args.filter)(period_dat)].sort_values(by='int_dateRecipient')\
		[['study_idDonor', 'study_idRecipient']].drop_duplicates().reset_index().\
		assign(
			study_idDonorLabel = lambda k: 'LLV' + np.arange(k.shape[0]).astype(str)+'-I',
			study_idRecipientLabel = lambda k: 'LLV' + np.arange(k.shape[0]).astype(str)+'-R')

	rccs_dat = pd.read_csv(args.rccsDat, sep='\t').\
		query("~int_date.isnull()").\
		assign(int_date = lambda k: pd.to_datetime(k.int_date)).\
		assign(numeric_int_date = lambda k: k.int_date.apply(lambda r: numeric_from_datetime(r)))

	plot_style()

	fig = plt.figure(figsize=(5.4*3, 4.8*np.ceil(plot_couples.shape[0]/3)), layout="constrained")
	gs = GridSpec(3, np.ceil(plot_couples.shape[0]/3).astype(int), figure=fig)
	axs = []
	for cdx, couple in plot_couples.iterrows():
		axs.append(plot_couple_vl_traj(
			fig,
			gs[np.floor(cdx/3).astype(int), cdx % 3],
			couple, period_dat, rccs_dat, config))

	os.makedirs("figures/pdf", exist_ok=True)
	fig.savefig('figures/pdf/select_vl_traj.pdf')
	os.makedirs("figures/eps", exist_ok=True)
	fig.savefig('figures/eps/select_vl_traj.eps')
	plt.close()



if __name__ == '__main__':
    run()

