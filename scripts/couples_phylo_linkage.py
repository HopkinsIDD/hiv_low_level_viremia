import sys
import argparse
import pandas as pd
import numpy as np
import baltic as bt
import matplotlib.pyplot as plt
from itertools import product


def label_nodes(tre):
	for idx, item in enumerate(tre.getInternal()):
		item.name = 'node_' + str(idx)
	return(tre)



def get_within_ind_nodes(tree, min_bl):
	from collections import defaultdict
	# list of tips associated with each individual
	ind_tips = defaultdict(lambda: [])
	for tip in tree.getExternal():
		ind_tips[tip.name.split('_')[0]].append(tip)
	# list of monophyletic nodes for each individual
	within_ind_nodes = defaultdict(lambda: set())
	#ind = 'R34729'
	#tips = ind_tips[ind]
	for ind, tips in ind_tips.items():
		if len(tips) > 1:
			for tip in tips:
				go = True
				test_node = tip
				while go:
					# if this node is monophyletic for this individual
					# or if this is a leaf and has a branch length  ~0
					# then we iterate one node up the tree
					# code written this way so that
					# R30131 in p24 is not identified as putative superinfection
					if ((len(set([i.split('_')[0] for i in test_node.parent.leaves])) == 1)  | 
							((test_node.branchType == 'leaf') & (test_node.length < min_bl))):
						test_node = test_node.parent
					else:
						within_ind_nodes[ind].add(test_node)
						go = False
		elif tips[0].length < min_bl:
			within_ind_nodes[ind].add(tips[0].parent)
		else:
			within_ind_nodes[ind].add(tips[0])
	return(within_ind_nodes)


def get_anc_ind_nodes(within_ind_nodes, tree, min_bl):
	from collections import defaultdict
	anc_ind_nodes = defaultdict(lambda: [])
	for ind, mono_nodes in within_ind_nodes.items():
		if len(mono_nodes) > 1:
			# if more than one within individual node, need to cluster closely related ones
			ind_tips = np.unique(np.hstack(
				[
					[j for j in i.leaves if ind in j] if i.branchType == 'node' else 
						i.name for i in mono_nodes]))
			all_anc_nodes = set([i.parent for i in mono_nodes])
			# sort by x value to get most ancestral
			all_anc_nodes = sorted(list(all_anc_nodes), key=lambda k: k.x)
			used_tips = []
			for anc_node in all_anc_nodes:
				if len(used_tips) < len(ind_tips):
					anc_ind_nodes[ind].append(anc_node)
					used_tips.extend([i for i in anc_node.leaves if ind in i])
		else:
			# if only one individual node, then just get ancestor
			anc_ind_nodes[ind].append(list(mono_nodes)[0].parent)
	return(anc_ind_nodes)


def format_list_col(x):
	fmt_x = [np.array(i.replace('[', '').lstrip().replace(']', '').split()).astype(np.dtypes.StringDType()) for i in x]
	return(np.array([np.pad(i, (0, 2-len(i)), constant_values=np.nan) for i in fmt_x]))


def get_paths(tree):
	import itertools
	# x-vals of all objects
	x_vals = {i.name: i.x for i in tree.Objects}
	parent_x_vals = {i.name: i.parent.x if i.parent else np.nan for i in tree.Objects}
	paths = {}
	for idx,i in enumerate(tree.Objects):
		paths[i.name] = set([(i.name, i.x)])
		# iterate up the tree
		x = i
		while hasattr(x.parent, 'name'):
			paths[i.name].add((x.parent.name, x.parent.x))
			x = x.parent
	# get MRCAs
	mrcas = []
	for pair in itertools.combinations(paths.keys(), 2):
		mrcas.append([*pair, *sorted(paths[pair[0]] & paths[pair[1]], key=lambda tup: tup[1])[-1]])
	# format out df
	mrcas = pd.DataFrame(mrcas,
		columns=['obj1', 'obj2', 'mrca', 'mrca_x']).\
		assign(
			obj1_x = lambda k: k.obj1.map(x_vals),
			obj1_parent_x = lambda k: k.obj1.map(parent_x_vals),
			obj2_x = lambda k: k.obj2.map(x_vals),
			obj2_parent_x = lambda k: k.obj2.map(parent_x_vals))
	return(mrcas)


def get_dist(obj1, obj2, cpd_dict):
	if obj1 == obj2:
		return(0)
	return(cpd_dict[tuple(sorted((obj1, obj2)))])


def get_null_distr(use_couples, within_ind_nodes, cpd_dict):
	from itertools import combinations, product
	# all tips in tree
	all_ids = [i for i in within_ind_nodes.keys() if 'K03455' not in i]
	# all true couples
	true_couples = set([tuple(sorted([i.study_idDonor, i.study_idRecipient])) for idx, i in use_couples.iterrows()])
	# iterate over potential couples and get distance
	null_distr = []
	for putative_null_couple in combinations(all_ids, 2):
		if tuple(sorted(putative_null_couple))  not in true_couples:
			null_distr.append(
				min([
					get_dist(i[0].name, i[1].name, cpd_dict) for 
						i in product(
							within_ind_nodes[putative_null_couple[0]],
							within_ind_nodes[putative_null_couple[1]])]))
	return(np.array(null_distr))


def run():
	parser = argparse.ArgumentParser()
	parser.add_argument('--tree')
	parser.add_argument('--label')
	parser.add_argument('--periodCopiesDat')
	args = parser.parse_args()
	#args.periodCopiesDat = 'output/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d.tsv'
	#args.tree = 'output/p24_phylo/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d_p24_focal_nn.fasta.treefile'
	#args.label = 'p24'
	if not args.label:
		args.label = args.tree
	print(f'Asessing topological clustering using treefile {args.tree}.', file=sys.stderr)
	# read in couples data
	periods = pd.read_csv(args.periodCopiesDat, sep='\t')
	periods[['copies1', 'copies2']] = format_list_col(periods.donor_copies).astype(float)
	
	# unique couples
	couples = periods[['study_idDonor', 'study_idRecipient']].drop_duplicates()
	# serodiscordant
	seroconverting_couples = periods.sort_values(by='round_x').groupby(['study_idDonor', 'study_idRecipient']).last().\
		reset_index().\
		query('conversion == True')

	# read in tree
	min_bl = 5E-4
	with open(args.tree, 'r') as fp:
		tree = bt.make_tree(fp.readlines()[0].strip())
		_ = tree.traverse_tree()
		# collapse nearly 0 branch lengths
		tree = label_nodes(tree.collapseBranches(collapseIf=lambda x: x.length < min_bl))

	
	# get pairwise genetic distance between all objects in the tree
	paths = get_paths(tree).\
		assign(cpd = lambda k: (k.obj1_x - k.mrca_x) + (k.obj2_x - k.mrca_x))

	cpd_dict = {tuple(sorted((i.obj1, i.obj2))): i.cpd for idx, i in paths.iterrows()}
	# for each individual, get list of their monophyletic nodes
	# also includes nodes on which an individual's tip falls on with BL = min_bl
	within_ind_nodes = get_within_ind_nodes(tree, min_bl)
	# consolidate into ancestral nodes
	anc_ind_nodes = get_anc_ind_nodes(within_ind_nodes, tree, min_bl)


	# for each linked couple where both individuals are in the tree, evaluate evidence of linkage
	tree_seroconverting_couples = seroconverting_couples[
		np.isin(seroconverting_couples.study_idDonor.values, list(anc_ind_nodes.keys())) & 
			np.isin(seroconverting_couples.study_idRecipient.values, list(anc_ind_nodes.keys()))][['study_idDonor', 'study_idRecipient']].copy()
	print(f'{tree_seroconverting_couples.shape[0]} couples have tips in the tree for both donor and recipient',
		file=sys.stderr)


	# are donors putatively superinfected?
	superinfected_donors = \
		seroconverting_couples.study_idDonor.values[
			np.isin(seroconverting_couples.study_idDonor, 
				[key for key, val in anc_ind_nodes.items() if len(val) > 1])]
	print(f'{superinfected_donors.shape[0]} donors show evidence of putative superinfection: {', '.join(superinfected_donors)}',
		file=sys.stderr)

	# add phylogenetic linkage
	tree_seroconverting_couples['phylo_linked'] = tree_seroconverting_couples.apply(lambda k: 
		any(np.hstack([[k.study_idRecipient in j for j in i.leaves] for i in anc_ind_nodes[k.study_idDonor]])),
		axis=1)
	print(f'{tree_seroconverting_couples["phylo_linked"].sum()} couples are phylogenetically linked')

	# save to output file
	tree_seroconverting_couples.to_csv('.'.join(args.tree.split('.')[:-1]) + f'_{args.label}_phylo_dist.tsv', sep='\t', index=None)

	# among LLV couples how many are phylo linked
	llv_tree_seroconverting_couples = tree_seroconverting_couples.merge(
		seroconverting_couples.query('((copies1 < 1000)&(copies1 > 200)) | ((copies2 < 1000)&(copies2 > 200))')\
			[['study_idDonor', 'study_idRecipient']],
			on=['study_idDonor', 'study_idRecipient'],
			how='inner')
	print(f'{llv_tree_seroconverting_couples.shape[0]} couples with index LLV have index and partner tips in tree')
	print(f'{llv_tree_seroconverting_couples["phylo_linked"].sum()} couples with index LLV are phylogenetically linked')

	# minimal cophenetic dsitance between each transmission pair
	tree_seroconverting_couples['phylo_min_d'] = tree_seroconverting_couples.apply(lambda k:
		np.hstack([get_dist(j[0], j[1], cpd_dict) for 
				j in product([i.name for i in within_ind_nodes[k.study_idDonor]],
			[i.name for i in within_ind_nodes[k.study_idRecipient]])]).min(),
		axis=1)


	# null distirbution of cophenetic distances
	null_distr = get_null_distr(couples, within_ind_nodes, cpd_dict)

	# print some basic statistics
	# plot bar charts
	bins = np.arange(0, null_distr.max()+0.005, 0.005)
	fig, axs = plt.subplots(3,1,figsize=(6.4, 4.8*1.5), constrained_layout=True,
			sharex=True)
	axs[0].hist(null_distr,
		bins=bins,
		facecolor='#eaeaea',
		edgecolor='#333333',
		label='random pairs')
	axs[1].hist(tree_seroconverting_couples.query('phylo_linked == False').phylo_min_d,
		bins=bins,
		facecolor='indianred',
		edgecolor='#333333',
		label='couples (non-phylolinked)')
	axs[2].hist(tree_seroconverting_couples.query('phylo_linked == True').phylo_min_d,
		bins=bins,
		facecolor='steelblue',
		edgecolor='#333333',
		label='couples (phylolinked)')
	_ = [ax.legend() for ax in axs]
	_ = [ax.set_ylabel('pairs') for ax in axs]
	axs[-1].set_xlabel('cophenetic distance (subs/site)')
	fig.suptitle(args.label)
	fig.savefig(f'figures/{".".join(args.label.split("/")[-1])}_phylo_linkage.pdf')
	plt.close()



if __name__ == "__main__":
    run()




