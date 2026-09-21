import os
import pandas as pd
import numpy as np
import argparse
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.collections import LineCollection
import matplotlib.patheffects as path_effects
from matplotlib.patches import Patch
import baltic as bt
import string
from collections import defaultdict
import subprocess
try:
	from scripts.utils import plot_style, split_col
except:
	from utils import plot_style,  split_col


def plot_linkage_sum(ax, linkage_sum, llv_couples, config):
	regions = ['p24', 'gp41']
	# add summary columns
	for region in regions:
		linkage_sum[f'{region}_linked'] = \
			np.where(
				linkage_sum[f'{region}_genolinked'].isnull().values,
				'Missing',
				(linkage_sum[[f'{region}_genolinked', f'{region}_phylolinked']] == True).values.all(axis=1).astype(str))
	def summarize_linkage(df):
		return(df.groupby([f'{i}_linked' for i in regions]).\
			size().\
			sort_values(ascending=False).\
			reset_index(name='n').\
			assign(consensus = lambda k: 
				np.where(
					(k[[f'{i}_linked' for i in regions]] == 'True').any(axis=1),
					"True",
					np.where(
						(k[[f'{i}_linked' for i in regions]] == 'False').any(axis=1),
						'False',
						'Missing'))))
	# group and get counts
	linkage_sum_grouped = pd.concat([
		summarize_linkage(linkage_sum).\
			assign(type='all'),
		summarize_linkage(linkage_sum.merge(llv_couples, how='inner',on=['study_idDonor', 'study_idRecipient'])).\
			assign(type='llv')])
	linkage_sum_grouped = linkage_sum_grouped.\
		merge(
			linkage_sum_grouped.query('type == "all"')[['p24_linked', 'gp41_linked']].\
				assign(x = lambda k: np.arange(k.shape[0])),
			how='left',
			on=['p24_linked', 'gp41_linked']).\
		assign(x = lambda k: np.where(k.type == "all", k.x-0.175, k.x+0.175))
	linkage_sum_grouped['fill_cat'] = \
		linkage_sum_grouped['type'] + '_' + (linkage_sum_grouped['consensus'] == 'True').astype(str)
	# plot
	ax.axis('off')
	axin1 = ax.inset_axes([0, 0.25, 1, 0.75])
	axin2 = ax.inset_axes([0, 0, 1, 0.23])
	axin1.bar(linkage_sum_grouped.x, linkage_sum_grouped.n,
		edgecolor='#333333',
		color= linkage_sum_grouped.fill_cat.map({
			'all_True': config['color_seroconversion'], 
			'llv_True': config['color_seroconversion'],
			'all_False': config['color_seroconversion_unlinked'], 
			'llv_False': config['color_seroconversion_unlinked']}),
		hatch = linkage_sum_grouped.type.map({
			'all': '', 
			'llv': '//'}),
		hatchcolor = config['color_200_1000_copies'],
		width=0.3, zorder=2)
	axin1.grid(axis='y', zorder=1)
	axin1.set_ylabel('couples', size=20)
	axin1.set_yticks([0, 10, 20, 30, 40, 50],
		labels=[0,10,20,30,40,50], fontsize=20)
	axin1.set_ylim(0, 50)
	axin1.set_xticks([])
	_ = [axin1.spines[i].set_visible(False) for i in ['top', 'right']]
	# make legends by hand
	lg1 = axin1.legend(handles=
		[
	    	Patch(facecolor=config['color_seroconversion'], edgecolor='#333333'),
	    	Patch(facecolor=config['color_seroconversion'], edgecolor='#333333', 
	    		hatch='//', hatchcolor=config['color_200_1000_copies'])], 
	    labels=['all', 'index LLV'], 
	    loc='upper right', bbox_to_anchor=(0.5, 1.0),
	    title='verdict:\nlinked\n',
	    frameon=True,
	    fontsize=18,
	    title_fontproperties={'weight': 'bold', 'size': 18})
	#lg1._legend_box.align='left'
	lg1.get_title().set_multialignment('center')
	axin1.add_artist(lg1)
	lg2 = axin1.legend(handles=
		[
	    	Patch(facecolor=config['color_seroconversion_unlinked'], edgecolor='#333333'),
	    	Patch(facecolor=config['color_seroconversion_unlinked'], edgecolor='#333333', 
	    		hatch='//', hatchcolor=config['color_200_1000_copies'])], 
	    labels=['all', 'index LLV'], loc='upper right', 
	    title='verdict:\nambiguous, unlinked,\nor missing',
	    frameon=True,
	    fontsize=18,
	    title_fontproperties={'weight': 'bold', 'size': 18})
	lg2.get_title().set_multialignment('center')
	#axin1.legend(fontsize=20, handletextpad=0.25, columnspacing=1,
	#			frameon=True, edgecolor='#eaeaea')
	xmin,xmax = (-0.5, linkage_sum_grouped.x.max()+0.15+0.175)
	ylocs = {'Missing': 0, 'False': 1, 'True': 2}
	ymin,ymax = (-0.5, 3*len(region)-0.5)
	for idx in np.arange(0, np.ceil(ymax),2):
		axin2.fill_between([xmin, xmax], y1=idx-0.5, y2=idx+0.5, facecolor='#eaeaea', zorder=1)
	axin2.axhline(y=2.5, color='#d6d6d6', zorder=2)
	axin2.scatter(
		np.tile(
			linkage_sum_grouped.query('type == "all"').x+0.175,3*len(regions)),
		np.repeat(np.arange(0,3*len(regions)),
			linkage_sum_grouped.query('type == "all"').shape[0]),
		facecolor='#d6d6d6',
		edgecolor='#d6d6d6',
		linewidth='5',
		s=100,
		zorder=2)
	axin2.add_collection(LineCollection(
		[[(i.x+0.175, ylocs[i.p24_linked]+len(regions)+1), (i.x+0.175, ylocs[i.gp41_linked])] for 
			idx,i in linkage_sum_grouped.query('type == "all"').iterrows()],
			color='#333333', linewidth=2,
		zorder=3))
	axin2.scatter(np.tile(linkage_sum_grouped.query('type == "all"').x.values+0.175, 2),
		np.hstack([
			linkage_sum_grouped.query('type == "all"').p24_linked.map(ylocs).values + len(regions)+1,
			linkage_sum_grouped.query('type == "all"').gp41_linked.map(ylocs).values]),
		edgecolor='#333333',
		facecolor='#333333',
		zorder=4,
		s=100)
	axin2.set_xticks([])
	axin2.set_yticks(np.arange(0, len(regions)*3),
		labels=np.tile(['missing', 'unlinked', 'linked'], len(regions)),
		size=20)	
	for rdx, region in enumerate(regions[::-1]):
		axin2.text(-2.2, 0.9 + rdx*3, '  '+region+'  ', ha='right', va='center', color='#333333', size=20)
	axin2.set_ylim(-0.5, 3*len(regions)-0.5)
	_ = [axin2.spines[i].set_visible(False) for i in ['top', 'right', 'left', 'bottom']]
	_ = [ax.set_xlim(xmin, xmax) for ax in [axin1, axin2]]
	ax.set_title('Couples with\nseroconverting partners', size=24)
	return(ax)


def plot_subtree(ax, tree, donor_tips, recipient_tips, config):
	tips = donor_tips + recipient_tips
	tip_names = [i.name for i in tips]
	# get common ancestor of tips
	ca = tree.commonAncestor(tips)
	#if len(ca.parent.parent.parent.leaves) > 100:
	#	print('here')
	#	st = tree.subtree(ca.parent)
	#else:
	#	st = tree.subtree(ca.parent.parent.parent)
	st = tree.subtree(ca.parent.parent.parent)
	st.drawTree()
	if len(st.getExternal()) > 100:
		# which node has the most tips none of which are donor or recipient
		tmp = 0
		to_collapse =[i for i in 
				sorted(
					[i for i in st.getInternal() if len(set(tip_names) & set(i.leaves)) == 0],
					key = lambda k: -len(k.leaves)) if
				len(i.leaves) > 100]
		while len(to_collapse) > 0:
			st.collapseSubtree(to_collapse[0], f'collapsed_node_{tmp}', widthFunction=lambda x:len(x.leaves)/10)
			tmp += 1
			to_collapse =[i for i in 
				sorted(
					[i for i in st.getInternal() if len(set(tip_names) & set(i.leaves)) == 0],
					key = lambda k: -len(k.leaves)) if
				len(i.leaves) > 100]
		#st.traverse_tree()
		#st.sortBranches()
		st.drawTree()
	donor_tips = st.getExternal(lambda k: k.name in [i.name for i in donor_tips])
	recipient_tips = st.getExternal(lambda k: k.name in [i.name for i in recipient_tips])
	x_span = max([i.x for i in st.Objects]) - min([i.x for i in st.Objects])
	# plot
	'''
	ax.text(0, 1, 
		f'index: {recipient_tips[0].name.split("_")[0]}',
		transform=ax.transAxes, 
		ha='left',
		size=14,
		color=config['color_seroconversion'],
		path_effects=
		[path_effects.withStroke(linewidth=1, foreground="#595959")])
	ax.text(0, 0.94, 
		f'partner: {donor_tips[0].name.split("_")[0]}',
		transform=ax.transAxes, 
		ha='left',
		size=12,
		color=config['color_200_1000_copies'],
		path_effects=
		[path_effects.withStroke(linewidth=1, foreground="#595959")])
	'''
	st.plotTree(ax, zorder=2, colour='#595959')
	ssize=500
	ax.scatter([i.x for i in donor_tips],
		[i.y for i in donor_tips],
		edgecolor='white',
		facecolor='white',
		marker='>',
		linewidth=4,
		s=ssize, zorder=3,
		clip_on=False)
	ax.scatter([i.x for i in donor_tips],
		[i.y for i in donor_tips],
		edgecolor='#595959',
		facecolor=config['color_200_1000_copies'],
		marker='>',
		s=ssize, zorder=4,
		clip_on=False)
	ax.scatter([i.x for i in recipient_tips],
		[i.y for i in recipient_tips],
		edgecolor='white',
		facecolor='white',
		marker='<',
		linewidth=4,
		s=ssize, zorder=3,
		clip_on=False)
	ax.scatter([i.x for i in recipient_tips],
		[i.y for i in recipient_tips],
		edgecolor='#595959',
		facecolor=config['color_seroconversion'],
		marker='<',
		s=ssize, zorder=4,
		clip_on=False)
	# collapsed clades
	clades = [i for i in st.Objects if isinstance(i,bt.clade)]
	for clade in clades:
		ax.add_patch(
			plt.Polygon(
				([clade.x,clade.y-0.001*len(st.Objects)],
					[clade.x,clade.y+0.001*len(st.Objects)],
					[clade.x*2,clade.y+clade.width/2.0],
					[clade.x*2,clade.y-clade.width/2.0]),
				facecolor='#eaeaea',edgecolor='#595959',zorder=3))
	# scale bar
	span = 0.02
	y_adj = st.ySpan*0.015
	ax.plot([min([i.x for i in st.Objects]), min([i.x for i in st.Objects])],
		[-2*y_adj,0],
		color='#595959')
	ax.plot([min([i.x for i in st.Objects]), min([i.x for i in st.Objects])+span],
		[-y_adj, -y_adj],
		color='#595959')
	ax.plot([min([i.x for i in st.Objects])+span, min([i.x for i in st.Objects])+span],
		[-2*y_adj, 0],
		color='#595959')
	ax.text(min([i.x for i in st.Objects])+span/2, -st.ySpan*0.075,
		f'{span} subs./site',
		size=14,
		va='center',
		ha='center')
	ax.set_xticks([])
	ax.set_yticks([])
	_ = [ax.spines[i].set_visible(False) for i in ['left', 'right', 'top', 'bottom']]
	ax.set_ylim(-st.ySpan*0.025, st.ySpan*1.1)
	return(ax)

	


def run():
	parser = argparse.ArgumentParser()
	parser.add_argument('--periodDat',
		help='path to rccs couples data')
	parser.add_argument('--linkageSum', 
	    help='path to linkage summary data')
	parser.add_argument('--p24Tree', 
	    help='link to p24 tree')
	parser.add_argument('--gp41Tree', 
	    help='link to gp41 tree')
	parser.add_argument('--config',
		help='path to config csv', default='config/config.csv')
	args = parser.parse_args()
	# old:
	#args.periodDat = '../hiv-transmission-parameters/FreshRCCS_new/output/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d.tsv'
	#args.p24Tree = '../genetic_linkage/data/prox4.0_donor-all_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d_p24_focal_nn.fasta.treefile'
	#args.gp41Tree = '../genetic_linkage/data/prox4.0_donor-all_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d_gp41_focal_nn.fasta.treefile'
	#args.linkageSum = 'output/couples_linkage_sum.tsv'

	#args.periodDat = 'output/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d.tsv'
	#args.p24Tree = 'output/p24_phylo/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d_p24_focal_nn.fasta.treefile'
	#args.gp41Tree = 'output/gp41_phylo/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d_gp41_focal_nn.fasta.treefile'
	#args.linkageSum = 'output/couples_linkage_sum.tsv'
	config = {i[0]:i[1] for idx, i in pd.read_csv(args.config, header=None).iterrows()}
	linkage_sum = pd.read_csv(args.linkageSum, sep='\t')
	period_dat = pd.read_csv(args.periodDat, sep='\t').query('conversion == True')
	donor_copies = split_col(period_dat.donor_copies)
	donor_copies = [np.pad(i, (0, 2-len(i)), constant_values=np.nan) for i in donor_copies]
	period_dat[['donor_copies1', 'donor_copies2']] = donor_copies

	llv_couples = period_dat.sort_values(by='int_dateRecipient')\
		[((period_dat[['donor_copies1', 'donor_copies2']] > 200) & 
			((period_dat[['donor_copies1', 'donor_copies2']] < 1000))).any(axis=1) ].\
		assign(couple_idx = lambda k: 'LLV' + (np.arange(k.shape[0])+1).astype(str)).\
		assign(idDonor = lambda k: k.couple_idx + '-I', idRecipient = lambda k: k.couple_idx+'-P')\
		[['study_idDonor', 'study_idRecipient', 'idDonor', 'idRecipient']]
	# | ((study_idDonor == "R78014") & (study_idRecipient == "R04130"))').\	
	labels = np.array(['unlinked', 'unlinked', 'linked', 'missing'])
	linked_llv_couples = llv_couples.merge(linkage_sum, on=['study_idDonor', 'study_idRecipient'], how='inner').\
		assign(
			p24 = lambda k: k.p24_genolinked + k.p24_phylolinked,
			gp41 = lambda k: k.gp41_genolinked + k.gp41_phylolinked).\
		assign(
			p24 = lambda k: labels[np.nan_to_num(k.p24.values.astype(float), nan=3).astype(int)],
			gp41 = lambda k: labels[np.nan_to_num(k.gp41.values.astype(float), nan=3).astype(int)]).\
		query('(p24 == "linked") | (gp41 == "linked") | ((study_idDonor == "R78014") & (study_idRecipient == "R04130"))').\
		assign(gp41 = lambda k: np.where(k.study_idDonor ==  "R78014", 'likely linked', k.gp41)).\
		assign(linked = lambda k: 
			np.where(
				(k.gp41 == "linked") | (k.p24 == "linked"),
				'linked',
				np.where(
					(k.gp41 == "likely linked") | (k.p24 == "likely linked"),
					'likely linked',
					'unlinked')))


	trees = {'p24': bt.loadNewick(args.p24Tree), 'gp41': bt.loadNewick(args.gp41Tree)}
			


	plot_style()
	mpl.rcParams['hatch.linewidth'] = 3.0 
	fig = plt.figure(figsize=(6.4*2 + (6.4*0.3*linked_llv_couples.shape[0]), 4.8*1.75), layout="constrained")
	gs = GridSpec(2, 4 + linked_llv_couples.shape[0], figure=fig)
	axs = [
		fig.add_subplot(gs[:,:4]), 
		*[fig.add_subplot(gs[0,4+i]) for i in np.arange(linked_llv_couples.shape[0])],
		*[fig.add_subplot(gs[1,4+i]) for i in np.arange(linked_llv_couples.shape[0])]]

	axs[0] = plot_linkage_sum(axs[0], linkage_sum, llv_couples, config)


	ax_to_use = 1
	for region in ['p24', 'gp41']:
		for cdx, couple in linked_llv_couples.reset_index(drop=True).iterrows():
			if ~np.isnan(couple[f'{region}_phylolinked']):
				axs[ax_to_use] = plot_subtree(
					axs[ax_to_use], 
					trees[region],
					[i for i in trees[region].getExternal() if couple.study_idDonor in i.name],
					[i for i in trees[region].getExternal() if couple.study_idRecipient in i.name],
					config)
			else:
				axs[ax_to_use].axis('off')
				axs[ax_to_use].set_ylabel(' ')
			if region == 'p24':
				axs[ax_to_use].text(0.42, 1.1, 
					f'{couple.idDonor}',
					transform=axs[ax_to_use].transAxes, 
					ha='right',
					size=24,
					color=config['color_200_1000_copies'],
					path_effects=
					[path_effects.withStroke(linewidth=1, foreground="#595959")])
				axs[ax_to_use].text(0.5, 1.1, 
					r'$\rightarrow$',
					transform=axs[ax_to_use].transAxes, 
					ha='center',
					size=24,
					color='#eaeaea',
					path_effects=
					[path_effects.withStroke(linewidth=1, foreground="#595959")])
				axs[ax_to_use].text(0.58, 1.1, 
					f'{couple.idRecipient}',
					transform=axs[ax_to_use].transAxes, 
					ha='left',
					size=24,
					color=config['color_seroconversion'],
					path_effects=
					[path_effects.withStroke(linewidth=1, foreground="#595959")])
				axs[ax_to_use].text(0.5, 1.025, 
					f'verdict: {couple.linked}',
					transform=axs[ax_to_use].transAxes, 
					ha='center',
					size=20,
					color='#333333',
					fontweight='bold')
			axs[ax_to_use].text(0.5, 0.95, 
				f'{region}: {couple[region]}',
				transform=axs[ax_to_use].transAxes, 
				ha='center',
				size=20,
				color='#333333')
			ax_to_use += 1
			

	
	x_adj = [-0.2] + [-0.45]*13
	for adx,ax in enumerate(axs):
		ax.text(x_adj[adx], 1, list(string.ascii_uppercase)[adx], transform=ax.transAxes, 
			fontsize=20, fontweight='bold', va='top', ha='left')

	os.makedirs("figures/pdf", exist_ok=True)
	fig.savefig('figures/pdf/genetic_linkage.pdf')
	#os.makedirs("figures/eps", exist_ok=True)
	#fig.savefig('figures/eps/genetic_linkage.eps')
	subprocess.run(["pdftops", "-eps", "figures/pdf/genetic_linkage.pdf", "figures/eps/genetic_linkage.eps"])
	plt.close()




if __name__ == "__main__":
    run()

