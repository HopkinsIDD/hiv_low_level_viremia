import os
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import baltic as bt
import string
from collections import defaultdict
from matplotlib.lines import Line2D
try:
	from scripts.utils import plot_style, import_config, split_pad_col
except:
	from utils import plot_style, import_config, split_pad_col



def run():
	parser = argparse.ArgumentParser()
	parser.add_argument('--periodDat',
		help='path to rccs couples data')
	parser.add_argument('--p24Tree', 
	    help='link to p24 tree')
	parser.add_argument('--gp41Tree', 
	    help='link to gp41 tree')
	parser.add_argument('--config',
		help='path to config csv', default='config/config.csv')
	args = parser.parse_args()
	#args.periodDat = 'output/ve_constrained/periods.tsv'
	#args.p24Tree = 'output/p24_phylo/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d_p24_focal_nn.fasta.treefile'
	#args.gp41Tree = 'output/gp41_phylo/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d_gp41_focal_nn.fasta.treefile'
	config = import_config(args.config)

	p_d = pd.read_csv(args.periodDat, sep='\t')
	p_d[['donor_copies1', 'donor_copies2']] = split_pad_col(p_d['donor_copies']).astype(float)

	llv_sc_couples = p_d.sort_values(by='int_dateRecipient').\
		query('conversion == 1')\
		[((p_d[['donor_copies1', 'donor_copies2']] > 200) & 
			((p_d[['donor_copies1', 'donor_copies2']] < 1000))).any(axis=1) ].sort_values(by='int_dateRecipient')\
		[['study_idDonor', 'study_idRecipient']].drop_duplicates().reset_index().\
		assign(couple_index = lambda k: (np.arange(k.shape[0])+1).astype(str)).\
		assign(study_idDonorLabel = lambda k: 'LLV' + k.couple_index+'-I',
			study_idRecipientLabel = lambda k: 'LLV' + k.couple_index+'-P')
	print(llv_sc_couples.to_string())
	'''
	Base colors:
	indianred
	steelblue
	mediumseagreen
	rebeccapurple
	sienna
	orchid
	teal
	coral
	darkmagenta
	'''
	#colors = np.vstack([
	#	['#a44949', '#d77c7c'],
	#	['#386890', '#6a9bc3'],
	#	['#308f5a', '#62c28d'],
	#	['#51287a', '#845bad'],
	#	['#804124', '#b37456'],
	#	['#ae59ab', '#e18cde'],
	#	['#006666', '#329999'],
	#	['#cc6540', '#ff9872'],
	#	['#6f006f', '#a232a2']])

	#colors = np.vstack([
	#	['#3776b7', '#9bbadb'],
	#	['#25704c', '#92b7a5'],
	#	['#663399', '#b299cc'],
	#	['#a56813', '#f2b460'],
	#	['#209497', '#6ce0e4'],
	#	['#a4900e', '#f1dd5a'],
	#	['#51525e', '#9e9fab'],
	#	['#883d3e', '#d58a8a'],
	#	['#333333', '#eaeaea']])

	colors = np.vstack([
		['#1f77b4', '#8fbbd9'],
		['#ff7f0e', '#ffbf86'],
		['#2ca02c', '#95cf95'],
		['#d62728', '#ea9393'],
		['#9467bd', '#c9b3de'],
		['#e7b800', '#f3db7f'],
		['#264653', '#92a2a9'],
		['#2a9d8f', '#94cec7'],
		['#333333', '#eaeaea']])

	trees = {'p24': bt.loadNewick(args.p24Tree), 'gp41': bt.loadNewick(args.gp41Tree)}

	plot_style()
	fig = plt.figure(figsize=(6.4*2, 4.8*2), layout='constrained')
	gs = GridSpec(5, 2, figure=fig)
	axs = [
		fig.add_subplot(gs[:4,0]),
		fig.add_subplot(gs[:4,1]),
		fig.add_subplot(gs[4,:])]

	from collections import defaultdict
	tips = defaultdict(lambda: set([]))
	# plot trees
	for rdx, region in enumerate(['p24', 'gp41']):
		axs[rdx].grid(axis='x', color='#eaeaea', zorder=1)
		trees[region].plotTree(axs[rdx], zorder=2, colour='white', linewidth=0.5*2)
		trees[region].plotTree(axs[rdx], zorder=3, colour='#333333', linewidth=0.5)
		axs[rdx].autoscale()
		_ = [axs[rdx].spines[i].set_visible(False) for i in ['top', 'right', 'left']]
		axs[rdx].set_yticks([])
		axs[rdx].set_xlabel('divergence (substitutions/site)')
		axs[rdx].set_title(region, size=18)
		# add donor tips
		for ddx, donor in enumerate(llv_sc_couples.study_idDonor):
			donor_tips = trees[region].getExternal(lambda k: donor in k.name)
			if len(donor_tips) > 0:
				tips[donor].add(region)
			axs[rdx].scatter(
				[i.x for i in donor_tips],
				[i.y for i in donor_tips],
				edgecolor='white',
				facecolor='white',
				marker='>',
				linewidth=2,
				s=150, zorder=4)
			axs[rdx].scatter(
				[i.x for i in donor_tips],
				[i.y for i in donor_tips],
				edgecolor='#333333',
				facecolor=colors[ddx,0],
				marker='>',
				s=150, zorder=5)
		for redx, recipient in enumerate(llv_sc_couples.study_idRecipient):
			recipient_tips = trees[region].getExternal(lambda k: recipient in k.name)
			if len(recipient_tips) > 0:
				tips[recipient].add(region)
			axs[rdx].scatter(
				[i.x for i in recipient_tips],
				[i.y for i in recipient_tips],
				edgecolor='white',
				facecolor='white',
				marker='<',
				linewidth=2,
				s=150, zorder=4)
			axs[rdx].scatter(
				[i.x for i in recipient_tips],
				[i.y for i in recipient_tips],
				edgecolor='#333333',
				facecolor=colors[redx,1],
				marker='<',
				s=150, zorder=6)




	# add legend
	def add_label(x):
		if len(x) == 0:
			return('missing')
		elif x == set(['p24']):
			return('p24')
		elif x == set(['gp41']):
			return('gp41')
		elif x == set(['p24', 'gp41']):
			return('p24, gp41')


	axs[-1].legend(
		handles=
			list(np.hstack([
				[
					Line2D([0], [0], marker='>', color='None', label=i.study_idDonorLabel + f' ({add_label(tips[i.study_idDonor])})',
						markerfacecolor= colors[idx,0] if len(tips[i.study_idDonor]) > 0 else 'white',
						markeredgecolor='#333333' if len(tips[i.study_idDonor]) > 0  else '#eaeaea', 
						markersize=14),
					Line2D([0], [0], marker='<', color='None', label=i.study_idRecipientLabel + f' ({add_label(tips[i.study_idRecipient])})',
						markerfacecolor=colors[idx,1] if len(tips[i.study_idRecipient]) > 0  else 'white', 
						markeredgecolor='#333333' if len(tips[i.study_idRecipient]) > 0 else '#eaeaea', 
						markersize=14)] for 
				idx, i in llv_sc_couples.iterrows()] + \
				[Line2D([0], [0], marker='<', color='None', label='',
						markerfacecolor='white', 
						markeredgecolor='white', 
						markersize=0)]*2)),
		loc='center', ncols=np.ceil(llv_sc_couples.shape[0]/2), 
		title='Seroconverting low-level viremia couples\n', fontsize=11, title_fontsize=14, labelspacing=1.5, handletextpad=0.4)


	axs[-1].axis('off')
	axs[-1].plot([0.035, 0.9], [0.32, 0.32], color='#eaeaea')
	axs[-1].plot([0.26+0.01, 0.26+0.01], [-0.1, 0.75], color='#eaeaea', clip_on=False)
	axs[-1].plot([0.42+0.01, 0.42+0.01], [-0.1, 0.75], color='#eaeaea', clip_on=False)
	axs[-1].plot([0.58+0.01, 0.58+0.01], [-0.1, 0.75], color='#eaeaea', clip_on=False)
	axs[-1].plot([0.74+0.01, 0.74+0.01], [-0.1, 0.75], color='#eaeaea', clip_on=False)
	
	axs[-1].set_ylim(0,1)
	axs[-1].set_xlim(0,1)
	axs[-1].text(0.045, 0.62, 'Index:', fontsize=14, transform=axs[-1].transAxes)
	axs[-1].text(0.045, 0.4, 'Partner:', fontsize=14, transform=axs[-1].transAxes)
	axs[-1].text(0.045, 0.17, 'Index:', fontsize=14, transform=axs[-1].transAxes)
	axs[-1].text(0.045, 0.17-0.22, 'Partner:', fontsize=14, transform=axs[-1].transAxes)
	os.makedirs("figures/pdf", exist_ok=True)
	fig.savefig('figures/pdf/full_trees.pdf')
	os.makedirs("figures/eps", exist_ok=True)
	fig.savefig('figures/eps/full_trees.eps')
	plt.close()



if __name__ == '__main__':
    run()

