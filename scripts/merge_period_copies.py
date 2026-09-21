import sys
import os
import pandas as pd
import numpy as np
import cmdstanpy
from cmdstanpy import CmdStanModel
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
import glob


def numeric_from_datetime(dt):
	from calendar import isleap
	import pandas as pd
	days = np.array([366 if isleap(k) else 365 for k in dt.dt.year])
	res = dt.dt.year + dt.dt.dayofyear/days
	return(res)



def format_p_d_prox(p_d, v_win):
	# have to do a bunch of data wrangling if we're using most proximal viral load
	# first convert interview dates to numerics for easier handling
	# get just the donor VL observations that have some overlap with each coupling period
	# sort by observation date
	# round x = follow-up recipient round
	# round y = source round for viral load measurment
	# flag: is this necessary?
	# int_date: date of donor viral load measurement
	p_d = p_d.\
		assign(int_date = lambda k: numeric_from_datetime(pd.to_datetime(k.int_date))).\
		assign(
			min_valid_int_date = lambda k: k.int_date - v_win/2,
			max_valid_int_date = lambda k: k.int_date + v_win/2,
			int_dateRecipient = lambda k: numeric_from_datetime(pd.to_datetime(k.int_dateRecipient)),
			lag_int_dateRecipient = lambda k: numeric_from_datetime(pd.to_datetime(k.lag_int_dateRecipient))).\
		query('((max_valid_int_date >= lag_int_dateRecipient) & (min_valid_int_date <= int_dateRecipient))').\
		sort_values(by=['couple_idx', 'round', 'round_y'])
	# for each coupling period
	# ensure there is a VL measurement within args.vWin/2 of entire period
	p_d = p_d
	p_d = p_d.merge(
		p_d.groupby(['couple_idx', 'round'])\
			[['dt', 'int_date', 'lag_int_dateRecipient', 'int_dateRecipient']].\
			apply(lambda k: 
				(np.abs(
					k.int_date.values[:,np.newaxis] - \
					np.linspace(
						k.lag_int_dateRecipient.iloc[0],
						k.int_dateRecipient.iloc[0],
						1000)) < v_win/2).any(axis=0).all()).\
			reset_index(name='valid_vl_obs').query('valid_vl_obs == True'),
		how='inner', on=['couple_idx', 'round'])
	# for each coupling period, get most proximal viral load
	# measurement at grid points across period
	prox_vl = p_d.groupby(['couple_idx', 'round', 'dt'])\
				[['int_date', 'round_y', 'copies', 'lag_int_dateRecipient', 'int_dateRecipient']].\
				apply(lambda k: 
					[k.int_date.values, k.round_y.values,  k.copies.values, np.argmin(
						np.abs(
							k.int_date.values[:,np.newaxis] - \
							np.linspace(
								k.lag_int_dateRecipient.iloc[0],
								k.int_dateRecipient.iloc[0],
								1000)),
					axis=0)]).\
				reset_index(name='prox_vl')
	def trim_zeros(arr):
		return(arr[arr > 0])
	# prox vl column has values:
	# interview_date,round, copies, interview_date argmin
	prox_vl = prox_vl.assign(
			dt = lambda k: [trim_zeros(i['dt']*np.bincount(i.prox_vl[-1])/1000) for idx,i in k.iterrows()],
			donor_copies = lambda k: [i.prox_vl[2][np.unique(i.prox_vl[-1])] for idx,i in k.iterrows()],
			donor_copies_date = lambda k: [i.prox_vl[0][np.unique(i.prox_vl[-1])] for idx,i in k.iterrows()],
			donor_copies_round = lambda k: [i.prox_vl[1][np.unique(i.prox_vl[-1])] for idx,i in k.iterrows()]).\
		drop(['prox_vl'], axis=1)
	# merge both of these in
	# do some clean up of the Vl and date columns
	# and drop invalid observations
	p_d = p_d.\
		drop(['round_y', 'copies', 'int_date', 'min_valid_int_date', 'max_valid_int_date', 'valid_vl_obs'], axis=1).\
		rename(columns={'dt': 'dt_sum'}).\
		drop_duplicates().\
		merge(
			prox_vl,
			how='left',
			on=['couple_idx', 'round'])
	return(p_d)


'''
copies_dat = 'data/internal/RCCSdata_R001_R019_VOIs_clean.tsv'
period_dat = 'data/internal/donor-noacute_recipient-exclusivelymonogamous.tsv'
v_win = 4
'''
def format_d_p_d(copies_dat, period_dat, v_win):
	# read in viral load data, sort by individual and round
	d = pd.read_csv(copies_dat, sep='\t').\
		query('~copies.isnull()').\
		sort_values(by=['study_id', 'round'])
	# unique LODs
	LODs = d.lod.drop_duplicates().sort_values().values
	# add indicator for tx availability
	# add indicator for tx use
	# tabulate number of previous individual BD observations for each observation
	d = d.\
		assign(
			tx_avail = lambda k: 1*(k['round'] > 10),
			tx = lambda k: 1*(k['round'] > k['last_arvmedFalse']))
	LODs = d['lod'].drop_duplicates().sort_values().values
	prev_putative_tx = \
		k = pd.concat([
			d.query('round >= 11').query('copies == 0.0').query('lod == @lod').assign(x=1).pivot(index='study_id', columns='round', values='x').\
				reindex(columns=d['round'].drop_duplicates().sort_values()).\
					fillna(0).\
					cumsum(axis=1).\
					reset_index().\
					melt(id_vars='study_id').\
					sort_values(by='round').\
					assign(lod=lod, value = lambda k: k.groupby('study_id').value.transform(lambda g: g.shift())).\
					fillna(0) for lod in LODs]).\
		pivot(index=['study_id', 'round'], columns='lod', values='value').\
		sort_values(by='round').\
		groupby(['study_id', 'round']).transform(lambda g: g.ffill().fillna(0)).\
		reindex(columns=LODs).fillna(0).\
		assign(
			prev_putative_tx = lambda k: [np.hstack(i) for idx,i in k.iterrows()],
			lod = lambda k: ';'.join(k.columns[:-1].astype(str))).reset_index()
	d = d.\
		merge(prev_putative_tx[['study_id', 'round', 'prev_putative_tx']],
			how='left',
			on=['study_id', 'round']).\
		assign(prev_putative_tx = lambda k:
			[i.prev_putative_tx if isinstance(i.prev_putative_tx, (list, np.ndarray)) else np.repeat(0.0, LODs.shape[0]) for 
				idx, i in k.iterrows()])
	# read in coupling period data
	# drop those with negative or 0 duration
	# flag: why do these exist?
	# merge in all donor viral load data
	# int_date: date of donor viral load measurement
	p_d = pd.read_csv(period_dat, sep='\t').\
		query("dt > 0").\
		merge(d[['study_id', 'round', 'copies', 'int_date']].drop_duplicates(),
			left_on='study_idDonor',
			right_on='study_id',
			how='inner').\
		drop('study_id', axis=1).\
		rename(columns={'round_x': 'round'})
	p_d = format_p_d_prox(p_d, v_win)
	q = max([i.shape[0] for i in p_d.donor_copies])
	if q > 2:
		raise Exception('code currently supports at most 2 viral loads per coupling period')
	'''
	# We want to retain the following viral load observations
	# 1) Pre-treatment
	# 2) Those implicated in partnership intervals
	p_d_unique_vl = p_d[['study_idDonor', 'sexDonor', 'donor_copies_round']].\
		assign(
			study_idDonor = lambda k: [[i.study_idDonor]*len(i.donor_copies_round) for idx, i in k.iterrows()],
			sexDonor = lambda k: [[i.sexDonor]*len(i.donor_copies_round) for idx, i in k.iterrows()])

	keep_vl = pd.concat([
			d.query('round <= last_arvmedFalse')[['study_id', 'round', 'sex']],
			pd.DataFrame({
				'study_id': np.hstack(p_d_unique_vl.study_idDonor.values), 
				'round': np.hstack(p_d_unique_vl['donor_copies_round'].values),
				'sex': np.hstack(p_d_unique_vl.sexDonor.values)}).drop_duplicates()]).\
		drop_duplicates().\
		sort_values(by=['study_id', 'round']).\
		assign(obs_idx = lambda k: np.arange(k.shape[0])+1)
	'''
	keep_vl = d[['study_id', 'round', 'sex', 'tx_avail', 'tx', 'prev_putative_tx', 'lod', 'copies']].\
		sort_values(by=['study_id', 'round']).\
		assign(obs_idx = lambda k: np.arange(k.shape[0])+1)
	# add individual level identifiers
	keep_vl = keep_vl.merge(
		keep_vl[['study_id']].drop_duplicates().assign(ind_idx = lambda k: np.arange(k.shape[0])+1),
		how='left',
		on='study_id')
	# need to add observation level identifiers to p_d
	# first set up a dictionary
	from collections import defaultdict
	obs_idx_dict = {i: defaultdict(lambda: dict()) for i in keep_vl.study_id.drop_duplicates()}
	for idx, i in keep_vl.iterrows():
		obs_idx_dict[i.study_id][i['round']] = i['obs_idx']
	p_d = p_d.assign(donor_copies_obs_idx = lambda k: \
		[[obs_idx_dict[i.study_idDonor][j] for j in i.donor_copies_round] for idx, i in k.iterrows()])
	return(keep_vl, p_d, q, keep_vl)


def run():
	parser = argparse.ArgumentParser()
	parser.add_argument('--copiesDat', 
	    help='path to rccs viral load data')
	parser.add_argument('--periodDat',
		help='path to rccs contact period data')
	parser.add_argument('--vWin',
		type=float,
		default=3,
		help='size of viral load window when using proximal viral load measurements')
	parser.add_argument('--out', 
	    help='path to output file')
	args = parser.parse_args()
	#args.copiesDat = 'data/internal/RCCSdata_R001_R019_VOIs_clean.tsv'
	#args.periodDat = 'data/internal/donor-noacute_recipient-exclusivelymonogamous.tsv'
	#args.vWin = 4


	keep_vl, p_d, q, keep_vl  = format_d_p_d(args.copiesDat, args.periodDat, args.vWin)
	
	print(f'Merged periods include {p_d.conversion.sum()} seroconversions accumulated over {np.hstack(p_d.dt.values).sum()} couple-years of follow-up during {p_d.shape[0]} coupling periods among {p_d[['study_idDonor', 'study_idRecipient']].drop_duplicates().shape[0]} couples',
		file=sys.stderr)

	# Some basic checks on the prev_putative_tx col
	keep_vl = keep_vl.assign(n_prev_putative_tx = lambda k: [i.prev_putative_tx.sum() for idx, i in k.iterrows()])

	if keep_vl.query('tx_avail==0').query('n_prev_putative_tx > 0').shape[0] > 0:
		raise Exception('no prev putative tx prior to treatment availability')

	ever_prev_putative_tx = keep_vl.query("n_prev_putative_tx > 0").study_id.drop_duplicates()
	

	print(f'Individuals with putative tx', file=sys.stderr)
	for study_id in ever_prev_putative_tx :
		print(keep_vl.query('study_id == @study_id')[['study_id', 'round', 'tx_avail', 'tx', 'copies', 'lod','prev_putative_tx']],
			file=sys.stderr)


	# format list cols
	for col in ['dt', 'donor_copies', 'donor_copies_date', 'donor_copies_round', 'donor_copies_obs_idx']:
		p_d[col] = [';'.join(np.array(i).astype(str)) for i in p_d[col]]


	p_d.to_csv(args.out +'_p_d.tsv', sep='\t', index=None)


	keep_vl['prev_putative_tx'] = [';'.join(np.array(i).astype(str)) for i in keep_vl['prev_putative_tx']]
	keep_vl.to_csv(args.out + '_d.tsv', sep='\t', index=None)

	
if __name__ == "__main__":
    run()

	

