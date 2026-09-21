import argparse
import pandas as pd
import numpy as np
import os
from datetime import datetime
import matplotlib.pyplot as plt
from itertools import product
try:
	from scripts.utils import import_aln, map_arr, numeric_from_datetime, plot_style
except:
	from utils import import_aln, map_arr, numeric_from_datetime, plot_style


def read_fasta(fp):
    name, seq = None, []
    for line in fp:
        line = line.rstrip()
        if line.startswith(">"):
            if name: yield (name, ''.join(seq))
            name, seq = line[1:], []
        else:
            seq.append(line)
    if name: yield (name, ''.join(seq))


def import_aln(fh):
    s_names = []
    all_s = ''
    #fh = open(fasta_path, 'rt')
    with fh as fasta:
        for h,s in read_fasta(fasta):
            s_names.append(h)
            all_s += s
    fh.close()
    n_seqs = len(s_names)
    size = int(len(all_s)/n_seqs)
    s_arr = np.array(list(all_s.lower())).reshape((n_seqs, size))
    return(np.array(s_names), s_arr)


def k2p_arr_d(arr1, arr2):
	gt1 = np.isin(arr1, ['A', 'C', 'G', 'T', 'a', 'c', 't', 'g'])
	gt2 = np.isin(arr2, ['A', 'C', 'G', 'T', 'a', 'c', 't', 'g'])
	arr1_clean = np.where(~gt1, 'n', arr1)
	arr2_clean = np.where(~gt2, 'n', arr2)
	# for each item in seq_arr is it purine or pyridime 
	from collections import defaultdict
	n_type_dict = defaultdict(lambda: np.nan)
	n_type_dict.update({
		'a': 0,
		'c': 1,
		'g': 0,
		't': 1})
	n_type1 = map_arr(arr1_clean, n_type_dict)
	n_type2 = map_arr(arr2_clean, n_type_dict)
	# empty distance array
	d_arr = np.zeros((arr1.shape[0], arr2.shape[0]))
	ol_arr = np.zeros((arr1.shape[0], arr2.shape[0]))
	for sdx,s in enumerate(arr1_clean):
		# white sites are genotyped in both
		n = (gt1[sdx] & gt2)
		n_sum = n.sum(axis=1)
		# which sites have differences
		d = np.where((s != arr2_clean) & n)
		# for each difference
		# is this difference a transition (FALSE) 
		# or transversion (TRUE)
		d_type = n_type1[sdx][d[1]] == n_type2[d[0], d[1]]
		p = np.zeros(arr2.shape[0])
		q = np.zeros(arr2.shape[0])
		np.add.at(p, d[0], d_type)
		np.add.at(q, d[0], ~d_type)
		p /= n_sum
		q /= n_sum
		d_arr[sdx,:] = -(1/2)*np.log((1-2*p-q)*np.sqrt(1-2*q))
		ol_arr[sdx,:] = n.sum(axis=1)
	return(ol_arr, d_arr)


def get_study_idx(names):
		from collections import defaultdict
		out = defaultdict(lambda: np.array([]).astype('int64'))
		out.update({idx: np.array(i).astype('int64') for idx, i in pd.DataFrame(names, columns=['name']).\
			assign(idx = lambda k: np.arange(k.shape[0]),
				study_id = lambda k: [i.split('_')[0] for i in k.name]).\
				groupby('study_id').idx.agg(list).items()})
		return(out)

'''
def deduplicate_seqs(names, seqs):
	from collections import defaultdict
	name_idx = defaultdict(lambda: [])
	for idx, i in enumerate(names):
		name_idx[i].append(idx)
	keep = []
	# for each sequence name
	# get the one with the greatest number of ACTG values
	for key, val in name_idx.items():
		keep.append(
			val[np.argmax(np.isin(seqs[np.array(val),:], 
				['A', 'C', 'G', 'T', 'a', 'c', 'g', 't']).sum(axis=1))])
	return(names[keep], seqs[keep,:])
'''

def format_list_col(x):
	fmt_x = [np.array(i.replace('[', '').lstrip().replace(']', '').split()).astype(np.dtypes.StringDType()) for i in x]
	return(np.array([np.pad(i, (0, 2-len(i)), constant_values=np.nan) for i in fmt_x]))


def write_seqs(names, seqs, out):
	with open(out, 'w') as fp:
		for idx, i in enumerate(names):
			fp.write('>'+i+'\n'+''.join(seqs[idx])+'\n')


def get_null_distr(all_couples, region_dat, validity_threshold):
	from itertools import combinations
	# all available IDs wtih sequence data 
	all_ids = [key for key, value in region_dat['idx'].items() if 
		value.shape[0] > 0 and 'K03455' not in key]
	true_couples = set([tuple(sorted([i.study_id1, i.study_id2])) for idx, i in all_couples.iterrows()])
	# iterate over potential couples and get distance
	null_d = []
	null_overlap = []
	# if input to combinations is sorted, so are output tuples
	putative_null_couples = set(combinations(sorted(all_ids), 2))
	# back out true couples
	putative_null_couples = list(putative_null_couples - true_couples)
	# then select 1 million random pairs
	rng = np.random.default_rng(seed=111)
	use_idx = rng.choice(np.arange(len(putative_null_couples)), 1000000,replace=False)
	for putative_null_couple_idx in use_idx:
		putative_null_couple = putative_null_couples[putative_null_couple_idx]
		overlap, d = k2p_arr_d(
			region_dat['seqs'][region_dat['idx'][putative_null_couple[0]],:],
			region_dat['seqs'][region_dat['idx'][putative_null_couple[1]],:])
		null_overlap.append(overlap.flatten())
		null_d.append(d.flatten())
	null_d = np.hstack(null_d)
	null_overlap = np.hstack(null_overlap)
	return(null_d[null_overlap > validity_threshold])


def get_win_distr(region_dat, validity_threshold):
	win_distr =[]
	# get all individuals with >1 sequence
	for i in [key for key, value in region_dat['idx'].items() if value.shape[0] > 1]:
		overlap, d = k2p_arr_d(
				region_dat['seqs'][region_dat['idx'][i],:],
				region_dat['seqs'][region_dat['idx'][i],:])
		d[overlap <= validity_threshold] = np.nan
		win_distr.append(d[np.triu_indices(region_dat['idx'][i].shape[0], k=1)])
	win_distr = np.hstack(win_distr)
	return(win_distr[~np.isnan(win_distr)])


def get_j(arr1, arr2):
	t = np.linspace(0, max(arr1.max(), arr2.max()), 10000)[:,np.newaxis]
	j = (arr1 <= t).sum(axis=1)/arr1.shape[0] + (arr2 > t).sum(axis=1)/arr2.shape[0] - 1
	return(t[np.argmax(j)][0])


def get_nearest_neighbors(query_names, query_seqs, db_names, db_seqs, valid, n=5):
	ol, d = k2p_arr_d(query_seqs, db_seqs)
	# mask like-like and invalid comparisons
	keep = (query_names[:,np.newaxis] != db_names) & (ol > valid)
	d[~keep] = np.inf
	nn = np.unique(np.argsort(d, axis=1)[:,:n].flatten())
	return(nn)


def add_dists(seroconverting_couples, regions, validity_thresholds):
	for region in regions:
		seroconverting_couples[f'{region}_overlap'] = np.nan
		seroconverting_couples[f'{region}_overlap'] = \
			seroconverting_couples[f'{region}_overlap'].astype('object')
		seroconverting_couples[f'{region}_d'] = np.nan
		seroconverting_couples[f'{region}_d'] = \
			seroconverting_couples[f'{region}_d'].astype('object')
		all_overlap = []
		all_d = []
		all_min_d = []
		#all_first_d = []
		all_donor_dates = []
		all_recipient_dates = []
		for idx, i in seroconverting_couples.iterrows():
			if (i[f'{region}Donor'].shape[0] > 0) and (i[f'{region}Recipient'].shape[0] > 0):
				overlap,d = k2p_arr_d(seroconverting_couples.loc[idx,f'{region}Donor'], seroconverting_couples.loc[idx,f'{region}Recipient'])
				overlap = overlap.flatten()
				d = d.flatten()
				# get just valid comparisons
				valid = overlap > validity_thresholds[region]
				all_overlap.append(overlap[valid])
				all_d.append(d[valid])
				valid_donor_dates = np.tile(i[f'{region}DonorDate'], i[f'{region}RecipientDate'].shape[0])[valid]
				valid_recipient_dates = np.repeat(i[f'{region}RecipientDate'], i[f'{region}DonorDate'].shape[0])[valid]
				all_donor_dates.append(valid_donor_dates)
				all_recipient_dates.append(valid_recipient_dates)
				all_min_d.append(d[valid].min())
				#all_first_d.append(
				#	d[valid][
				#		(valid_donor_dates == valid_donor_dates.min()) & \
				#			(valid_recipient_dates == valid_recipient_dates.min())].\
				#		min())
			else:
				all_overlap.append(np.nan)
				all_d.append(np.nan)
				all_min_d.append(np.nan)
				#all_first_d.append(np.nan)
				all_donor_dates.append(np.nan)
				all_recipient_dates.append(np.nan)
		seroconverting_couples[f'{region}_d_overlap'] = all_overlap
		seroconverting_couples[f'{region}_d'] = all_d
		seroconverting_couples[f'{region}_d_donorDate'] = all_donor_dates
		seroconverting_couples[f'{region}_d_recipientDate'] = all_recipient_dates
		seroconverting_couples[f'{region}_d_min'] = all_min_d
	return(seroconverting_couples)


#couples = seroconverting_couples

def get_save_nearest_neighbors(couples, dat, regions, validity_thresholds, out_base):
	# skip WG 
	for region in [i for i in regions if i.lower()!='wg']:
		use = \
			pd.concat([
					couples[[i[f'{region}Donor'].shape[0] > 0 for idx, i in couples.iterrows()]]\
						[['study_idDonor',  f'{region}DonorID', f'{region}Donor']].\
						groupby('study_idDonor').first().reset_index().\
						rename(columns={
							'study_idDonor': 'study_id',
							f'{region}DonorID': 'id',
							f'{region}Donor': region}),
					couples[[i[f'{region}Recipient'].shape[0] > 0 for idx, i in couples.iterrows()]]\
						[['study_idRecipient',f'{region}RecipientID', f'{region}Recipient']].\
						groupby('study_idRecipient').first().reset_index().\
						rename(columns={
							'study_idRecipient': 'study_id',
							f'{region}RecipientID': 'id',
							f'{region}Recipient': region})]).\
				groupby('study_id').first().reset_index()
		use_names = np.hstack([np.repeat(i.study_id, i[region].shape[0]) for idx,i in use.iterrows()])
		use_ids = np.array([i for i in np.hstack(use.id)])
		use_seqs = np.vstack(use[region].values)
		nearest_neighbors_idx = \
			get_nearest_neighbors(
				use_names,
				use_seqs,
				np.array([i.split('_')[0] for i in dat[region]['names']]),
				dat[region]['seqs'],
				validity_thresholds[region])
		# and save
		save_seqs = pd.DataFrame({
				'name': np.hstack([
					[f'{use_ids[idx]}' for idx in range(use_ids.shape[0])],
					dat[region]['names'][nearest_neighbors_idx]]),
				'seq': np.hstack([
					[''.join(i) for i in use_seqs],
					[''.join(i) for i in dat[region]['seqs'][nearest_neighbors_idx]]])}).\
			drop_duplicates()
		with open(out_base+f'_{region}_focal_nn.fasta', 'w') as fp:
			# add ref 
			ref_idx = [idx for idx, i in enumerate(dat[region]['names']) if 'K03455.1' in i][0]
			_ = fp.write(f'>{dat[region]["names"][ref_idx]}\n{"".join(dat[region]["seqs"][ref_idx])}\n')
			for idx, i in save_seqs.iterrows():
				_ = fp.write(f'>{i['name']}\n{i.seq}\n')


def run():
	parser = argparse.ArgumentParser()
	#parser.add_argument('--wg')
	parser.add_argument('--p24')
	parser.add_argument('--gp41')
	parser.add_argument('--periodCopiesDat')
	parser.add_argument('--allCouplesDat')
	args = parser.parse_args()
	#args.periodCopiesDat = 'output/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d.tsv'
	#args.allCouplesDat = 'data/couples_AN_20220904_clean.tsv.gz'
	#args.p24 = 'data/rccs_all_p24_aln_dg.fasta'
	#args.gp41 = 'data/rccs_all_gp41_aln_dg.fasta'
	# read in sequence data

	out_dir = 'output/genetic_linkage'
	os.makedirs(out_dir, exist_ok=True)

	validity_thresholds = {
		'wg': 2500,
		'p24': 387*0.95,
		'gp41': 402*0.95}
	validity = {
		'wg': lambda k: np.isin(k, ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't']).sum(axis=1) > validity_thresholds['wg'],
		'p24': lambda k: np.isin(k, ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't']).sum(axis=1) > validity_thresholds['p24'],
		'gp41': lambda k: np.isin(k, ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't']).sum(axis=1) > validity_thresholds['gp41']}

	dat = {}
	regions = ['p24', 'gp41']
	for input_dat_idx, input_dat in enumerate([args.p24, args.gp41]):
		names, seqs = import_aln(open(input_dat, 'r'))
		#names, seqs = deduplicate_seqs(names, seqs)
		valid = validity[regions[input_dat_idx]](seqs)
		names = names[valid]
		seqs = seqs[valid,:]
		#dates = np.array([numeric_from_datetime(datetime.strptime(i.split('_')[1], '%Y-%m-%d'))
		#	if 'K03455.1' not in i else np.nan for i in names])
		dates = np.array([datetime.strptime(i.split('_')[1], '%Y-%m-%d')
			if 'K03455.1' not in i else np.nan for i in names])
		# save valid deduplicated alignment
		write_seqs(names, seqs, '.'.join(input_dat.split('.')[:-1])+'_dedup.fasta')
		idx = get_study_idx(names)
		dat[regions[input_dat_idx]] = {'names': names, 'dates': dates, 'seqs': seqs, 'idx': idx}

	'''
	Checking date counts
	x = pd.DataFrame(np.unique([int(i.split('_')[1].split('-')[0]) for i in dat['p24']['names'] if 'K03455.1' not in i], return_counts=True), index=['year', 'n']).T
	x.query('year < 2010').n.sum()
	x.query('year >= 2010').n.sum()

	y = pd.DataFrame(np.unique([int(i.split('_')[1].split('-')[0]) for i in dat['gp41']['names'] if 'K03455.1' not in i], return_counts=True), index=['year', 'n']).T
	y.query('year < 2010').n.sum()
	y.query('year >= 2010').n.sum()
	'''
	# coupling periods
	periods = pd.read_csv(args.periodCopiesDat, sep='\t')
	periods[['copies1', 'copies2']] = format_list_col(periods.donor_copies).astype(float)
	
	# unique couples
	couples = periods.sort_values(by='round').groupby(['study_idDonor', 'study_idRecipient']).last().\
		reset_index()

	# add sequence data
	# filter for only validity
	# add couple sequence arrays and valid comparison columns
	# validity:
	# whole genome: at least 2,500 nt of overlap
	# p24: 95% overlap
	# gp41: 95% overlap
	nucs = ['A', 'C', 'G', 'T', 'a', 'c', 'g', 't']
	couples = couples.\
		assign(
			p24Donor = lambda k: [dat['p24']['seqs'][dat['p24']['idx'][i],:] for i in k.study_idDonor],
			gp41Donor = lambda k: [dat['gp41']['seqs'][dat['gp41']['idx'][i],:] for i in k.study_idDonor],
			p24DonorDate = lambda k: [dat['p24']['dates'][dat['p24']['idx'][i]] for i in k.study_idDonor],
			gp41DonorDate = lambda k: [dat['gp41']['dates'][dat['gp41']['idx'][i]] for i in k.study_idDonor],
			p24DonorID = lambda k: [dat['p24']['names'][dat['p24']['idx'][i]] for i in k.study_idDonor],
			gp41DonorID = lambda k: [dat['gp41']['names'][dat['gp41']['idx'][i]] for i in k.study_idDonor],
			p24Recipient = lambda k: [dat['p24']['seqs'][dat['p24']['idx'][i],:] for i in k.study_idRecipient],
			gp41Recipient = lambda k: [dat['gp41']['seqs'][dat['gp41']['idx'][i],:] for i in k.study_idRecipient],
			p24RecipientDate = lambda k: [dat['p24']['dates'][dat['p24']['idx'][i]] for i in k.study_idRecipient],
			gp41RecipientDate = lambda k: [dat['gp41']['dates'][dat['gp41']['idx'][i]] for i in k.study_idRecipient],
			p24RecipientID = lambda k: [dat['p24']['names'][dat['p24']['idx'][i]] for i in k.study_idRecipient],
			gp41RecipientID = lambda k: [dat['gp41']['names'][dat['gp41']['idx'][i]] for i in k.study_idRecipient]
			)

	# subset to just seroconverters 
	keep_cols = [
		'study_idDonor', 'study_idRecipient',
		'copies1', 'copies2',
		'p24Donor', 'gp41Donor', #'wgDonor',
		'p24DonorDate', 'gp41DonorDate', #'wgDonorDate',
		'p24DonorID', 'gp41DonorID',
		'p24Recipient', 'gp41Recipient', #'wgRecipient',
		'p24RecipientDate', 'gp41RecipientDate', #'wgRecipientDate'
		'p24RecipientID', 'gp41RecipientID'
		]
	seroconverting_couples = couples.query('conversion == True').copy()[keep_cols]

	# get a list of nearest neighbors for all seroconverting couple sequences
	get_save_nearest_neighbors(seroconverting_couples, dat, regions, validity_thresholds,
		out_dir + f'/{".".join(args.periodCopiesDat.split("/")[-1].split(".")[:-1])}')

	# get dists
	seroconverting_couples = add_dists(seroconverting_couples, regions, validity_thresholds)

	# get a null distribution for non-couples in each region
	# read in a list of all couples ever reported
	all_couples = pd.read_csv(args.allCouplesDat, sep='\t')\
		[['study_id1', 'study_id2']].drop_duplicates()
	null_distrs = {region: 
		get_null_distr(all_couples, dat[region], validity_thresholds[region]) for 
			region in regions}

	win_distrs = {region:
		get_win_distr(dat[region], validity_thresholds[region]) for region in regions}

	# define cut-off by Youden's J. 
	cutoff = {region: 
		get_j(
				win_distrs[region],
				null_distrs[region]) for 
			region in regions}

	# print some basic statistics
	print(f'{couples.shape[0]} unique couples')
	print(f'{seroconverting_couples.shape[0]} seroconversions')
	for region in regions:
		# add cutoff to output file
		seroconverting_couples[f'{region}_j'] = cutoff[region]
		seroconverting_couples[f'{region}_linked'] = \
			np.where(
				pd.isna(seroconverting_couples[f'{region}_d_min']),
				'',
				seroconverting_couples[f'{region}_d_min'] < cutoff[region])
		# save null distribution
		with open(out_dir + '/'+'.'.join(args.periodCopiesDat.split('/').split('.')[:-1]) + f'_{region}_null.tsv', 'w') as fp:
			fp.write('\n'.join(null_distrs[region].astype(str)))
		# save win distr
		with open(out_dir + '/'+'.'.join(args.periodCopiesDat.split('/').split('.')[:-1]) + f'_{region}_win.tsv', 'w') as fp:
			fp.write('\n'.join(win_distrs[region].astype(str)))
		
		print(f'{sum([i.shape[0] > 0 for i in seroconverting_couples[f'{region}Donor']])}' + \
			f' index partners have valid {region} sequence data')
		print(f'{sum([i.shape[0] > 0 for i in seroconverting_couples[f'{region}Donor']])}' + \
			f' index partners with seroconverting partners have valid {region} sequence data')
		print(f'{sum(~seroconverting_couples[f'{region}_d_overlap'].isnull())}' + \
			f' couples have valid {region} sequence data for both partners')
		print(f'{np.nansum(seroconverting_couples[f'{region}_d_min'] < 0.05)}' + \
			f' couples with valid {region} sequence data for both partners have a genetic distance < 0.05 subs/site')
		print(f'{np.nansum(seroconverting_couples[f'{region}_d_min'] < cutoff[region])}' + \
			f' couples with valid {region} sequence data for both partners have a genetic distance < {cutoff[region]} subs/site')

	
	# statistics specific to llv couples
	llv_couples = seroconverting_couples.query('((copies1 < 1000)&(copies1 > 200)) | ((copies2 < 1000)&(copies2 > 200))')
	print(f'{llv_couples.shape[0]} seroconversions with index LLV')
	for region in regions:
		print(f'{sum([i.shape[0] > 0 for i in llv_couples[f'{region}Donor']])}' + \
			f' index partners with LLV have valid {region} sequence data')
		print(f'{sum([i.shape[0] > 0 for i in llv_couples[f'{region}Donor']])}' + \
			f' index partners with LLV with seroconverting partners have valid {region} sequence data')
		print(f'{sum(~llv_couples[f'{region}_d_overlap'].isnull())}' + \
			f' couples with index LLV have valid {region} sequence data for both partners')
		print(f'{np.nansum(llv_couples[f'{region}_d_min'] < 0.05)}' + \
			f' couples with index LLV with valid {region} sequence data for both partners have a genetic distance < 0.05 subs/site')
		print(f'{np.nansum(llv_couples[f'{region}_d_min'] < cutoff[region])}' + \
			f' couples with index LLV with valid {region} sequence data for both partners have a genetic distance < {cutoff[region]} subs/site')



	'''
	llv_couples[['study_idDonor', 'study_idRecipient', 'p24_min_d', 'gp41_min_d']]
	R16249 -> R04315: gp41 mono
	R39840 -> R23733: gp41 mono, p24 linked by another sequence co-clusters
	R69651 -> R03822: gp41 mono, p24 mono
	R71227 -> R12968: Do not co-cluster in p24
	R84146 -> R28479: Do not co-cluster in p24
	'''

	supr_couples = seroconverting_couples.query('(copies1 < 200) & (copies2 < 200)')
	print(f'{supr_couples.shape[0]} seroconversions with index with VL<200')
	for region in regions:
		print(f'{sum([i.shape[0] > 0 for i in supr_couples[f'{region}Donor']])}' + \
			f' index partners with VL<200 have valid {region} sequence data')
		print(f'{sum([i.shape[0] > 0 for i in supr_couples[f'{region}Donor']])}' + \
			f' index partners with VL<200 with seroconverting partners have valid {region} sequence data')
		print(f'{sum(~supr_couples[f'{region}_d_overlap'].isnull())}' + \
			f' couples with index with VL<200 have valid {region} sequence data for both partners')
		print(f'{np.nansum(supr_couples[f'{region}_d_min'] < 0.05)}' + \
			f' couples with index with VL<200 with valid {region} sequence data for both partners have a genetic distance < 0.05 subs/site')
		print(f'{np.nansum(supr_couples[f'{region}_d_min'] < cutoff[region])}' + \
			f' couples with index with VL<200 with valid {region} sequence data for both partners have a genetic distance < {cutoff[region]} subs/site')

	# format and save to output file
	for col in [i for i in seroconverting_couples.columns if 'date' in i.lower()]:
		seroconverting_couples[col] = \
			[';'.join([datetime.strftime(j, '%Y-%m-%d') for j in i]) if 
					isinstance(i, np.ndarray) else '' for 
				i in seroconverting_couples[col] ]
	for col in [i for i in seroconverting_couples.columns if (('_d' in i.lower()) & ('date' not in i.lower()) & ('_min' not in i.lower())) | ('ID' in i)]:
		seroconverting_couples[col] = \
			[';'.join(i.astype(str)) if 
					isinstance(i, np.ndarray) else '' for 
				i in seroconverting_couples[col] ]
	for col in list(product(regions, ['Donor', 'Recipient'])):
		seroconverting_couples[''.join(col)] = \
			[';'.join([''.join(j) for j in i]) if 
					i.shape[0] > 0 else '' for 
				i in seroconverting_couples[''.join(col)]]
	
	seroconverting_couples.to_csv(
		out_dir + '/'+
		'.'.join(args.periodCopiesDat.split('/')[-1].split('.')[:-1]) + '_dist.tsv',
		sep='\t',
		index=None)


	# plot bar charts
	plot_style()
	bins = np.arange(0, max([i.max() for i in null_distrs.values()])+0.005, 0.005)
	fig, axs = plt.subplots(3,2,figsize=(6.4*1.5, 4.8*1.5), constrained_layout=True, sharex='all', sharey='row')
	for idx, region in enumerate(regions):
		labels=[None, None]
		if idx == 0:
			labels=['intra-individual', 'random pairs', 'couples']
		axs[0, idx].hist(win_distrs[region].astype(float),
			bins=bins,
			facecolor='#848484',
			edgecolor='#333333',
			label=labels[0])
		axs[1, idx].hist(null_distrs[region],
			bins=bins,
			facecolor='#eaeaea',
			edgecolor='#333333',
			label=labels[1])
		x = seroconverting_couples[f'{region}_d_min'].dropna().values
		axs[1, idx].hist(x,
			bins=bins,
			facecolor='steelblue',
			edgecolor='#333333')
		axs[2, idx].hist(x[x >= cutoff[region]],
			bins=bins,
			facecolor='#eaeaea',
			edgecolor='#333333',
			label=labels[2]+' (non-genolinked)' if labels[1] else labels[1])
		axs[2, idx].hist(x[x < cutoff[region]],
			bins=bins,
			facecolor='#848484',
			edgecolor='#333333',
			label=labels[2] + ' (genolinked)' if labels[1] else labels[1])
		axs[0,idx].set_title(region)
		axs[1,idx].set_title(' ')
		axs[0,idx].set_xlabel(' ')
		axs[2,idx].set_xlabel('genetic distance (subs/site)')
		axs[0,idx].set_ylabel(' ')
		axs[1,idx].set_ylabel(' ')
		_ = [axs[i,idx].axvline(cutoff[region],
			color='#333333',
			ls='--') for i in [0,1,2]]
		


	_ = [axs[i,0].legend() for i in [0,1,2]]
	axs[0,0].set_ylabel('pairs')
	axs[1,0].set_ylabel('pairs')
	fig.savefig('test.pdf')
	fig.savefig('figures/pdf/genetic_linkage_distr.pdf')
	fig.savefig('figures/eps/genetic_linkage_distr.eps')
	plt.close()




if __name__ == "__main__":
    run()


