import os
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import string
try:
	from scripts.utils import plot_style, import_config
except:
	from utils import plot_style, import_config


def run():
	parser = argparse.ArgumentParser()
	parser.add_argument('--copiesDat', 
	    help='path to rccs viral load data')
	parser.add_argument('--config',
		help='path to config csv', default='config/config.csv')
	args = parser.parse_args()
	config = import_config(args.config)
	
	#args.copiesDat = 'data/RCCSdata_R001_R019_vl.tsv'

	d = pd.read_csv(args.copiesDat, sep='\t').\
		query('~copies.isnull()')

	log_max_y = np.ceil(np.log10(d.copies.max()))
	# just pre-treatment
	pt_d = d.query('round <= last_arvmedFalse')
	print(f'{(pt_d.query('round <= 10').copies == 0).sum()/pt_d.query('round <= 10').shape[0]}% of pre-treatment viral loads prior to or during round 10 are BD')
	print(f'{(pt_d.query('round > 10').copies == 0).sum()/pt_d.query('round > 10').shape[0]}% of pre-treatment viral loads after round 10 are BD')
	print(f'{(pt_d.query('round <= 10').copies < 400).sum()/pt_d.query('round <= 10').shape[0]}% of pre-treatment viral loads prior to or during round 10 are <400')
	print(f'{(pt_d.query('round > 10').copies < 400).sum()/pt_d.query('round > 10').shape[0]}% of pre-treatment viral loads after round 10 are <400')
	plot_style()
	bd_val = 100
	fig, axs = plt.subplots(1,2, figsize=(6.4*2, 4.8), constrained_layout=True)
	b0 = axs[0].hist(np.log10(np.where(pt_d.query('round <= 10').copies <400, bd_val, pt_d.query("round <= 10").copies)),
		bins=np.arange(np.log10(bd_val), log_max_y, 0.1),
		facecolor='#eaeaea',
		edgecolor='#333333', zorder=2)
	b1 = axs[1].hist(np.log10(np.where(pt_d.query('round > 10').copies <400, bd_val, pt_d.query("round > 10").copies)),
		bins=np.arange(np.log10(bd_val), log_max_y, 0.1),
		facecolor='#eaeaea',
		edgecolor='#333333', zorder=2)
	vl_ticks = np.hstack([[np.log10(bd_val) + 0.05], np.log10([400, 1000, 10000, 100000, 10**6, 10**7])])
	vl_labels = ['<400'] + \
		[f'{np.round(10**k,0).astype(int):,}' if 10**k < 1000000 else f'{int(10**k/1000000)}M'for k in vl_ticks[1:]]
	_ = [ax.set_xticks(vl_ticks, labels=vl_labels, fontsize=10) for ax in axs]
	_ = [ax.set_xlabel('viral load\n(copies/mL)') for ax in axs]
	_ = [ax.grid(axis='y', color='#eaeaea', zorder=1) for ax in axs]
	axs[0].set_ylabel('participant-visits')
	axs[1].set_ylabel(' ')
	axs[0].set_ylim(0, np.ceil(b0[0].max()/25+1)*25)
	axs[1].set_ylim(0, np.ceil(b1[0].max()/200+1)*200)
	axs[0].set_title('1994 through 2004')
	axs[1].set_title('2005 through 2020')
	x_adj = [-0.15, -0.15]
	for adx,ax in enumerate(axs):
		ax.text(x_adj[adx], 1.05, list(string.ascii_lowercase)[adx], transform=ax.transAxes, 
			fontsize=12, fontweight='bold', va='top', ha='left')


	fig.suptitle('People with HIV\nself-reporting no ART use')
	fig.savefig("figures/pdf/pre_treatment_copies.pdf")
	fig.savefig("figures/eps/pre_treatment_copies.eps")
	plt.close()

	tx_d = d.query('round > last_arvmedFalse').query('lod == 150')
	print(f'{(tx_d.copies <=150).sum()/tx_d.shape[0]}% of Abbott post-treatment viral loads are <150')
	
		
	fig, axs = plt.subplots(1,1, figsize=(6.4*1, 4.8), constrained_layout=True)
	bd_val = 5
	b0 = axs.hist(np.log10(np.where(tx_d.query('round > 10').copies <=150, 5, tx_d.query("round > 10").copies)),
		bins=np.hstack([[np.log10(5), np.log10(5)+0.1], np.arange(np.log10(150), log_max_y, 0.1)]),
		facecolor='#eaeaea',
		edgecolor='#333333', zorder=2)
	axs.set_yscale('log')
	vl_ticks = np.hstack([[np.log10(bd_val) + 0.05], np.log10([150, 1000, 10000, 100000, 10**6, 10**7])])
	vl_labels = ['<150'] + \
		[f'{np.round(10**k,0).astype(int):,}' if 10**k < 1000000 else f'{int(10**k/1000000)}M'for k in vl_ticks[1:]]
	_ = [ax.set_xticks(vl_ticks, labels=vl_labels) for ax in [axs]]
	_ = [ax.set_xlabel('viral load\n(copies/mL)') for ax in [axs]]
	_ = [ax.grid(axis='y', color='#eaeaea', zorder=1) for ax in [axs]]
	axs.set_ylabel('participant-visits')
	axs.set_yscale('log')
	axs.set_ylim(0.01, 10**5)
	axs.set_title('People with HIV\nreporting ART use\n2005 through 2020')
	fig.savefig("figures/pdf/treatment_copies.pdf")
	fig.savefig("figures/eps/treatment_copies.eps")
	plt.close()


if __name__ == "__main__":
    run()
