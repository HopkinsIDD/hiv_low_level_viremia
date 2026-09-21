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
	from scripts.utils import plot_style, import_config
except:
	from utils import plot_style, import_config


def plot_vl_traj(fig, gs0, paired_vl, config, max_y=10**7):
	if 'use_color' not in config.keys():
		config['use_color'] = '#333333'
	log_max_y = np.log10(max_y)
	paired_vl_logged = paired_vl[['copies', 'lag_copies']].values
	bd_val = 10
	paired_vl_logged = np.log10(np.where(paired_vl_logged < 200, bd_val, paired_vl_logged))
	from matplotlib.gridspec import GridSpecFromSubplotSpec
	gs00 = GridSpecFromSubplotSpec(1,4, subplot_spec=gs0, wspace=0.0, hspace=0.0)
	ax1 = fig.add_subplot(gs00[0, 0])
	ax2 = fig.add_subplot(gs00[0, 1:3])
	ax3 = fig.add_subplot(gs00[0,3])
	#max_y = np.log10(paired_vl[['copies', 'lag_copies']].values.max())
	cuts = np.hstack([-1, np.log10(config['breaks_copies'][1:])])
	colors = ['', config['color_0_200_copies'],
		config['color_200_1000_copies'], config['color_1000_inf_copies']]
	for cdx in range(1, len(cuts)):
		cut = cuts[cdx]
		_ = ax1.hist(
			paired_vl_logged[:,1][(paired_vl_logged[:,1] < cut) & (paired_vl_logged[:,1] > cuts[cdx-1])], 
			bins=np.arange(0, log_max_y, 0.1),
			facecolor=colors[cdx], edgecolor='#333333',
			orientation='horizontal', linewidth=0.5)
		_ = ax3.hist(
			paired_vl_logged[:,0][(paired_vl_logged[:,0] < cut) & (paired_vl_logged[:,0] > cuts[cdx-1])], 
			bins=np.arange(0, log_max_y, 0.1),
			facecolor=colors[cdx], edgecolor='#333333',
			orientation='horizontal', linewidth=0.5)
	ax1.invert_xaxis()
	ax1.set_ylabel('\nviral load\n(copies/mL)', size=14)
	vl_ticks = np.log10([bd_val, 200, 1000, 10000, 100000, 10**6, 10**7])
	vl_labels = ['<200'] + \
		[f'{np.round(10**k,0).astype(int):,}' if 10**k < 1000000 else f'{int(10**k/1000000)}M'for k in vl_ticks[1:]]
	ax1.set_yticks(
		vl_ticks, labels=vl_labels, size=12)
	_ = [ax.set_xticks([]) for ax in [ax1, ax3]]
	_ = [ax.set_yticks([]) for ax in [ax2, ax3]]
	_ = [ax.set_ylim(0, log_max_y) for ax in [ax1, ax2, ax3]]
	_ = [ax1.spines[i].set_visible(False) for i in ['bottom',  'top', 'right']]
	_ = [ax3.spines[i].set_visible(False) for i in ['bottom', 'left', 'top', 'right']]
	for i in paired_vl_logged:
		_ = ax2.plot([0,1], i[::-1],
			color=config['use_color'], lw=0.25, zorder=3, dashes=[10,10])
	ax2.set_xlim(0,1)
	ax2.set_xticks([0,1], labels=['$t_0$', '$t_1$'])
	# text annotate ax3
	binned_fu_vl = np.digitize(10**(paired_vl_logged[:,0]), config['breaks_copies'])
	bin_counts = np.unique(binned_fu_vl, return_counts=True)
	bin_props = np.zeros(binned_fu_vl.max()+1)
	bin_props[bin_counts[0]] = bin_counts[1]/paired_vl_logged.shape[0]
	bin_prop_labels = np.round(100*bin_props[1:], 0).astype(int).astype(str) + '%'
	breaks_copies = config['breaks_copies'].copy()
	breaks_copies[0] = bd_val
	bin_prop_locs = \
		np.nan_to_num(np.log10(breaks_copies[:-1])) + \
			(np.log10(np.hstack([breaks_copies[1:-1], max_y])) - np.nan_to_num(np.log10(breaks_copies[:-1])))/2
	#bin_prop_locs = \
	#	np.nan_to_num(np.log10(config['breaks_copies'][:-1])) + \
	#		(np.log10(np.hstack([config['breaks_copies'][1:-1], max_y])) - np.nan_to_num(np.log10(config['breaks_copies'][:-1])))/2
	for ldx, label in enumerate(bin_prop_labels):
		ax3.text(paired_vl_logged.shape[0]*0.05, bin_prop_locs[ldx], label,  ha='left', va='center')
	for cut in np.log10(config['breaks_copies'][1:]):
		_ = [ax.axhline(cut, color='#eaeaea', zorder=0, lw=0.5) for ax in [ax1, ax2, ax3]]
	return([ax1, ax2, ax3])

#fig
#gs0 = gs[1,0:2]
#sim_vl = sim_dat.query('t <= 365.0')
def plot_sim_vl_traj(fig, gs0, sim_vl, config, max_y=10**7):
	if 'use_color' not in config.keys():
		config['use_color'] = '#333333'
	log_max_y = np.log10(max_y)
	paired_vl_logged = sim_vl.query('(t == 0) | (t == t.max())')[['idx', 't', 'v']].\
		pivot(index='idx', columns='t', values='v').reset_index().drop(columns=['idx']).\
		values[:,[1,0]]
	#paired_vl_logged = np.log10(np.where(paired_vl_logged == 0, 1, paired_vl_logged))
	#fig = plt.figure(figsize=(6.4*1.5, 4.8*1.25))
	#gs0 = GridSpec(2, 2, figure=fig)[0,0]
	from matplotlib.gridspec import GridSpecFromSubplotSpec
	gs00 = GridSpecFromSubplotSpec(1,4, subplot_spec=gs0, wspace=0.0, hspace=0.0)
	ax1 = fig.add_subplot(gs00[0, 0])
	ax2 = fig.add_subplot(gs00[0, 1:3])
	ax3 = fig.add_subplot(gs00[0,3])
	#max_y = np.log10(paired_vl[['copies', 'lag_copies']].values.max())
	cuts = np.hstack([-1, np.log10(config['breaks_copies'][1:])])
	colors = ['', config['color_0_200_copies'],
		config['color_200_1000_copies'], config['color_1000_inf_copies']]
	for cdx in range(1, len(cuts)):
		cut = cuts[cdx]
		_ = ax1.hist(
			paired_vl_logged[:,1][(paired_vl_logged[:,1] < cut) & (paired_vl_logged[:,1] > cuts[cdx-1])], 
			bins=np.arange(0, log_max_y, 0.1),
			facecolor=colors[cdx], edgecolor='#333333',
			orientation='horizontal', linewidth=0.5)
		_ = ax3.hist(
			paired_vl_logged[:,0][(paired_vl_logged[:,0] < cut) & (paired_vl_logged[:,0] > cuts[cdx-1])], 
			bins=np.arange(0, log_max_y, 0.1),
			facecolor=colors[cdx], edgecolor='#333333',
			orientation='horizontal', linewidth=0.5)
	ax1.invert_xaxis()
	ax1.set_ylabel('viral load\n(copies/mL)', size=14)
	vl_ticks = np.log10([1, 200, 1000, 10000, 100000, 10**6, 10**7])
	vl_labels = ['1'] + \
		[f'{np.round(10**k,0).astype(int):,}' if 10**k < 1000000 else f'{int(10**k/1000000)}M'for k in vl_ticks[1:]]
	ax1.set_yticks(
		vl_ticks, labels=vl_labels, size=12)
	_ = [ax.set_xticks([]) for ax in [ax1, ax3]]
	_ = [ax.set_yticks([]) for ax in [ax2, ax3]]
	_ = [ax.set_ylim(0, log_max_y) for ax in [ax1, ax2, ax3]]
	_ = [ax1.spines[i].set_visible(False) for i in ['bottom',  'top', 'right']]
	_ = [ax3.spines[i].set_visible(False) for i in ['bottom', 'left', 'top', 'right']]
	cols=[config['color_0_200_copies'], config['color_200_1000_copies'], config['color_1000_inf_copies']]
	for idx, i in sim_vl.groupby('idx'):
		x_interp = np.linspace(i.t.min(), i.t.max(), 10000)
		y_interp = np.interp(x_interp, i.t, i.v)
		cat_breaks = np.hstack([
			0, 
			np.where(
				np.digitize(10**y_interp[0:-1], config['breaks_copies']) != 
				np.digitize(10**y_interp[1:], config['breaks_copies']))[0]+1,
			y_interp.shape[0]])
		for jdx, j in enumerate(cat_breaks):
			if (jdx > 0):
				'''
				print(' ')
				print(j)
				print(i.t[cat_breaks[jdx-1]:j])
				print(i.v[cat_breaks[jdx-1]:j])
				print(np.digitize(10**i.v.values[j-1], config['breaks_copies'])-1)
				'''
				_ = ax2.plot(x_interp[cat_breaks[jdx-1]:j], y_interp[cat_breaks[jdx-1]:j],
				color=cols[np.digitize(10**y_interp[j-1], config['breaks_copies'])-1], lw=0.25, zorder=3)
			
	ax2.set_xlim(0,365)
	ax2.set_xticks([0, 90, 180, 270, 365])
	ax2.set_xlabel('days')
	#ax2.set_xticks([0,1], labels=['$t_0$', '$t_1$'])
	# text annotate ax3
	binned_fu_vl = np.digitize(10**(paired_vl_logged[:,0]), config['breaks_copies'])
	bin_counts = np.unique(binned_fu_vl, return_counts=True)
	bin_props = np.zeros(binned_fu_vl.max()+1)
	bin_props[bin_counts[0]] = bin_counts[1]/paired_vl_logged.shape[0]
	bin_prop_labels = np.round(100*bin_props[1:], 0).astype(int).astype(str) + '%'
	bin_prop_locs = \
			np.nan_to_num(np.log10(config['breaks_copies'][:-1])) + \
				(np.log10(np.hstack([config['breaks_copies'][1:-1], max_y])) - np.nan_to_num(np.log10(config['breaks_copies'][:-1])))/2
	for ldx, label in enumerate(bin_prop_labels):
		ax3.text(paired_vl_logged.shape[0]*0.05, bin_prop_locs[ldx], label,  ha='left', va='center')
	for cut in np.log10(config['breaks_copies'][1:]):
		_ = [ax.axhline(cut, color='#eaeaea', zorder=0, lw=0.5) for ax in [ax1, ax2, ax3]]
	return(ax1,ax2,ax3)


def plot_sampled_trans_risk(ax, sim_trans_risk, config):
	ax.grid(axis='both', color='#eaeaea', zorder=1)
	ax.fill_between(
		sim_trans_risk.index,
		sim_trans_risk['0.025'],
		sim_trans_risk['0.975'],
		color='white',
		linewidth=2, zorder=2)
	ax.fill_between(
		sim_trans_risk.index,
		sim_trans_risk['0.025'],
		sim_trans_risk['0.975'],
		color=config['color_trans_rate_95'],
		linewidth=0, zorder=2)
	ax.fill_between(
		sim_trans_risk.index,
		sim_trans_risk['0.25'],
		sim_trans_risk['0.75'],
		color=config['color_trans_rate_50'],
		linewidth=0, zorder=2)
	ax.plot(
		sim_trans_risk.index,
		sim_trans_risk['0.5'],
		color=config['color_trans_rate_mid'], zorder=5)
	ax.set_xlabel('simulated years')
	ax.set_xticks([0, 180, 365, 547.5, 730], labels=[0, 0.5, 1.0, 1.5, 2.0], size=10)
	ax.set_xlim([0,730])
	ax.set_ylim(0, 8)
	ax.set_yticks([0,2,4,6, 8], labels=['0', '2', '4', '6', '8'], size=10)
	ax.set_ylabel('cumulative\ntransmissions/\n100 persons')
	_ = [ax.spines[i].set_visible(False) for i in ['top', 'right']]
	return(ax)


def plot_fixed_trans_risk(fig, gs0, fixed_trans_risk, sampled_rate, config):
	from scipy.stats import gaussian_kde
	from matplotlib.gridspec import GridSpecFromSubplotSpec
	yr_fixed_trans_risk = fixed_trans_risk.query('t == 365.0')
	sampled_prob = 100*(1 - np.exp(-sampled_rate.values[:,0]*365))
	kde = gaussian_kde(sampled_prob)
	yr_fixed_trans_risk['label'] = 100*(1 - np.exp(-yr_fixed_trans_risk.label*365))
	# make sub-grids
	gs00 = GridSpecFromSubplotSpec(5,1, subplot_spec=gs0, wspace=0.0, hspace=0.0)
	ax0 = fig.add_subplot(gs00[0, 0])
	ax = fig.add_subplot(gs00[1:, 0])
	# ensure x limits are shared
	#xmin = sampled_rate.values.min()
	#xmax = yr_fixed_trans_risk.label.values.max() + yr_fixed_trans_risk.label.values.min() - xmin
	#xmin = sampled_rate.values.min()
	#xmax = yr_fixed_trans_risk.label.values.max()*1.05
	xmin = sampled_prob.min()
	xmax = 65
	xval = np.linspace(xmin, xmax, 100000)
	# plot density
	ax0.plot(xval, kde(xval), color=config['color_trans_rate_mid'])
	ax0.fill_between(xval, kde(xval), 0, color=config['color_trans_rate_95'], alpha=1)
	# plot estimated transmission risks
	ax.grid(axis='y', color='#eaeaea', zorder=1)
	ax.scatter(
		yr_fixed_trans_risk.label,
		yr_fixed_trans_risk['0.5'],
		facecolor='white',
		edgecolor='white',
		s=75,
		zorder=2)
	ax.scatter(
		yr_fixed_trans_risk.label,
		yr_fixed_trans_risk['0.5'],
		facecolor='#eaeaea',
		edgecolor=config['color_trans_rate_mid'],
		s=50,
		zorder=5)
	for rdx, row in yr_fixed_trans_risk.iterrows():
		for level in [('0.025', '0.975', 1, config['color_trans_rate_95']), ('0.25', '0.75', 5, config['color_trans_rate_50'])]:
			_ = ax.plot(
				[yr_fixed_trans_risk.label,yr_fixed_trans_risk.label],
				[yr_fixed_trans_risk[level[0]], yr_fixed_trans_risk[level[1]]],
				color='white',
				lw=level[2]+1,
				zorder=2)
			_ = ax.plot(
				[yr_fixed_trans_risk.label,yr_fixed_trans_risk.label],
				[yr_fixed_trans_risk[level[0]], yr_fixed_trans_risk[level[1]]],
				color=level[3],
				lw=level[2],
				zorder=4)
	ax0.set_yticks([])
	ax0.set_xticks([])
	ax0.set_ylabel(' est.', size=10)
	#ax.set_xlabel('treatment cessation\n' + r'rate ($\times10^{-3}$ year$^{-1}$)', size=12)
	ax.set_xlabel('annual risk of\nvirological failure (%)', size=12)
	ax.set_ylabel('transmissions/\n100 person-years', size=12)
	# convert units to Probabaility of cessation (year^-1)
	#ax.set_xticks(yr_fixed_trans_risk.label, labels=yr_fixed_trans_risk.label*1000, size=10)
	ax.set_xticks(yr_fixed_trans_risk.label,
		labels=np.round(yr_fixed_trans_risk.label, 0).astype(int),
		size=9)
	ax.set_yticks(np.arange(0, np.ceil(yr_fixed_trans_risk['0.975'].max())+1, 1),
		labels=np.arange(0, np.ceil(yr_fixed_trans_risk['0.975'].max())+1, 1).astype(int), size=10)
	_ = [iax.set_xlim(xmin, xmax) for iax in [ax0, ax]]
	_ = [[iax.spines[i].set_visible(False) for i in ['top', 'right']] for iax in [ax0, ax]]
	return(ax, ax0)


def run():
	parser = argparse.ArgumentParser()
	parser.add_argument('--llvDat',
		help='path to rccs viral load data')
	parser.add_argument('--suprDat')
	parser.add_argument('--simDat',
		help='path to simulated viral load data')
	parser.add_argument('--sampledRate',
		help='path to mcmc trace')
	parser.add_argument('--sampledRisk',
		help='path to simulated transmission risk')
	parser.add_argument('--fixRisk',
		help='path to simulated transmission risk under fixed rebound values'	)
	parser.add_argument('--config',
		help='path to config csv', default='config/config.csv')
	args = parser.parse_args()
	#args.vlDat = 'data/rccs_r1_r20_cohort_dat.tsv'
	#args.simDat = 'output/simulations/sampled/rccs_r1_r20_vl_pairs_treated_uniq_llv_fit_mcmc_sampled_0.000556853168868114_40_tstates.tsv.gz'
	#args.sampledRisk = "output/simulated_trajectory_risk_sum.tsv"
	#args.fixRisk = "output/fixed_simulated_trajectory_risk_sum.tsv"
	#args.sampledRate = "output/rccs_r1_r20_vl_pairs_treated_uniq_llv_fit_mcmc_sampled.tsv"
	config = import_config(args.config)

	paired_llv = pd.read_csv(args.llvDat, sep='\t')
	paired_supr = pd.read_csv(args.suprDat, sep='\t')
	
	sim_dat = pd.read_csv(args.simDat, sep='\t')
	
	# simulated risk under fit rebound values
	sampled_trans_risk = pd.read_csv(args.sampledRisk, sep='\t')

	# simualted risk under fixed rebound values
	fixed_trans_risk = pd.read_csv(args.fixRisk, sep='\t')

	# mcmc sampled rate values
	sampled_rate = pd.read_csv(args.sampledRate, header=None)


	plot_style()
	fig = plt.figure(figsize=(6.4*2, 4.8*1.25))
	gs = GridSpec(2, 4, figure=fig, wspace=0.75, hspace=0.5)
	#fig.add_subplot(gs[0,:2])
	#fig.add_subplot(gs[0,2:])
	#fig.add_subplot(gs[1,:2])
	#fig.add_subplot(gs[1,2])
	#fig.add_subplot(gs[1,3])
	config.update({'use_color': config['color_200_1000_copies']})
	# only first per person within each category
	ax1, ax2, ax3 = plot_vl_traj(fig, gs[0,0:2], paired_llv, config)
	ax7,ax8,ax9 = plot_sim_vl_traj(fig, gs[1,0:2],
		sim_dat.query('t <= 365.0'), config)
	config.update({'use_color': config['color_0_200_copies']})
	ax4, ax5, ax6 = plot_vl_traj(fig, gs[0,2:],
		paired_supr, config)

	ax10 = plot_sampled_trans_risk(fig.add_subplot(gs[1,2]), sampled_trans_risk, config)
	ax10.legend(handles=
		 [Line2D([0], [0], color=config['color_trans_rate_mid'], lw=1, label='median'),
		 	mpatches.Patch(color=config['color_trans_rate_50'], label='50% CI'),
		 	mpatches.Patch(color=config['color_trans_rate_95'], label='95% CI')],
		 fontsize=10,
		 loc='upper left',
		 )
	ax11, ax12 = plot_fixed_trans_risk(fig, gs[1,3], fixed_trans_risk, sampled_rate, config)
	ax11.legend(handles=
		 [Line2D([0], [0], marker='o', color=None, linewidth=0,
		 	markerfacecolor='#eaeaea', markeredgecolor=config['color_trans_rate_mid'], markersize=10,
		 	label='median'),
		 Line2D([0], [0], color=config['color_trans_rate_50'], linewidth=5,
		 	label='50% CI'),
		 Line2D([0], [0], color=config['color_trans_rate_95'], linewidth=1,
		 	label='95% CI')],
		 fontsize=10,
		 loc='upper left')
	# add titles
	ax13 = fig.add_subplot(gs[0,:])
	ax13.axis('off')
	ax13.set_title('On self-reported treatment\n', fontweight='bold')
	ax14 = fig.add_subplot(gs[1,:])
	ax14.axis('off')
	ax14.set_title('Simulation\n', fontweight='bold')
	ax2.set_title('200-1,000 copies/mL')
	ax5.set_title('<200 copies/mL')
	ax8.set_title('200-1,000 copies/mL')
	# add subplot labels
	x_adj = [-1.15, -1.15, -1.15, -0.6, -0.45]
	for adx,ax in enumerate([ax1, ax4, ax7, ax10, ax12]):
		ax.text(x_adj[adx], 1.05, list(string.ascii_uppercase)[adx], transform=ax.transAxes, 
			fontsize=12, fontweight='bold', va='top', ha='left')


	#fig.tight_layout()
	os.makedirs("figures/pdf", exist_ok=True)
	fig.savefig('figures/pdf/vl_trajs.pdf', bbox_inches='tight', pad_inches=0)
	os.makedirs("figures/eps", exist_ok=True)
	# fill between breaks EPS export, so use pdftops to convert pdf to eps
	subprocess.run(["pdftops", "-eps", "figures/pdf/vl_trajs.pdf", "figures/eps/vl_trajs.eps"])
	#fig.savefig('figures/eps/vl_trajs.eps')
	#fig.savefig('figures/eps/vl_trajs.pdf', bbox_inches='tight', pad_inches=0)
	#fig.savefig('figures/eps/vl_trajs.png', bbox_inches='tight', pad_inches=0)
	
	plt.close()




if __name__ == '__main__':
    run()




