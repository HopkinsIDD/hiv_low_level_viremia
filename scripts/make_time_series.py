import pandas as pd
import numpy as np
import argparse
from io import StringIO
import os
import ast
import sys


# flag do we need this
def set_dates(x):
		for i in [j for j in x.columns if 'date' in j]:
			x[i] = pd.to_datetime(x[i])
		return(x)


def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--dat', 
	    help='path to selected rccs data')
	args = parser.add_argument('--couplesDat',
		help='path to rccs couples data')
	args = parser.add_argument('--filter',
		default='lambda k: (k.monogamousRecipient == True) & (k.exclusively_partneredRecipient == True) & (k["round"] <= k.last_arvmedFalseDonor) & (k["last_finalhivN_dateDonor"].isnull() | k["int_dateRecipient"].isnull() | ((k["int_dateRecipient"] - k[["first_finalhivP_dateDonor", "last_finalhivN_dateDonor"]].median(axis=1)) >= pd.Timedelta(180, unit="days")))',
		help='lambda function filter to indicate valid donor observations')
	args = parser.add_argument('--out',
		default='out',
		help='output file name')
	args = parser.parse_args()
	#args.dat = 'data/RCCSdata_R001_R019_VOIs_clean.tsv.gz'
	#args.couplesDat = 'data/couples_AN_20220904_clean.tsv.gz'
	#args.donorFilter = 'lambda k: (k["last_finalhivN_dateDonor"].isnull() | k["int_dateRecipient"].isnull() | (k["int_dateRecipient"] - k[["first_finalhivP_dateDonor", "last_finalhivN_dateDonor"]].median(axis=1) >= pd.Timedelta(180, unit="days")))'
	#args.filter = 'lambda k: (k.monogamousRecipient == True) & (k.exclusively_partneredRecipient == True)'
	#test = lambda k: (k["round"] <= k.last_arvmedFalseDonor) & (k["last_finalhivN_dateDonor"].isnull() | k["int_dateRecipient"].isnull() | ((k["int_dateRecipient"] - k[["first_finalhivP_dateDonor", "last_finalhivN_dateDonor"]].median(axis=1)) >= pd.Timedelta(180, unit="days")))


	# read in data
	d = pd.read_csv(args.dat, sep='\t').\
		assign(int_date = lambda k: pd.to_datetime(k.int_date))
	couples = pd.read_csv(args.couplesDat, sep='\t', compression='gzip')

	# decorate couples with metadata columns
	# and fill missing HIV statuses
	ind_vars = ['study_id', 'sex', 'last_finalhivN', 'first_finalhivP', 'last_arvmedFalse', 'first_finalhivP_date', 'last_finalhivN_date']
	ind_round_vars = ['study_id', 'round', 'int_date', 'finalhiv', 'monogamous', 'exclusively_partnered']
	all_vars = ind_vars + ind_round_vars
	for j in range(1,3):
		couples = couples.\
			merge(
				d[ind_vars].\
					drop_duplicates().\
					rename(columns={
						i: f'{i}{j}' for i in ind_vars}),
				how='left',
				on=f'study_id{j}').\
			merge(
				d[ind_round_vars].\
					rename(columns={
						i: f'{i}{j}' for i in ind_round_vars if i != 'round'}),
				how='left',
				on=[f'study_id{j}', 'round'])
		# flag: could maybe use bfil and ffil here? 
		couples[f'finalhiv{j}'] = \
			np.where(
					couples['round'] >= couples[f'first_finalhivP{j}'],
					'P',
					couples[f'finalhiv{j}'])
		couples[f'finalhiv{j}'] = \
			np.where(
					couples['round'] <= couples[f'last_finalhivN{j}'],
					'N',
					couples[f'finalhiv{j}'])

	# subset to heterosexual serodiscordant couples
	# keep only the first round at which a couple becomes seroconcordant
	couples = couples.merge(
		couples.assign(sd = lambda k: (k.finalhiv1 != k.finalhiv2) & ~k.finalhiv1.isnull() & ~k.finalhiv2.isnull()).\
			groupby('couple')['sd'].any().reset_index().query('sd == True')[['couple']],
		how='inner',
		on='couple').\
		query('sex1 != sex2')

	couples.to_csv('.'.join(args.couplesDat.split('.')[:-1])+'_fmt.tsv', sep='\t')
	
	couples = couples.\
		query('~((finalhiv1 == "N") & (finalhiv2 == "N"))').\
		query('~(((finalhiv1 == "N") & (finalhiv2.isnull())) | ((finalhiv2 == "N") & (finalhiv1.isnull())))').\
		query('(finalhiv1 != finalhiv2) | ((first_finalhivP1 == round) | (first_finalhivP2 == round))')


	if (couples[['first_finalhivP1', 'first_finalhivP2']].isnull().sum(axis=1) == 2).any():
		raise Exception('all couples should have at least one HIV+ partner')

	# identifies donor and renames columns accordingly
	donor_idx = np.nan_to_num(couples[['first_finalhivP1', 'first_finalhivP2']].values, nan=np.inf).argmin(axis=1)
	get_vars = [i for i in all_vars if i !='round']
	couples = pd.concat([
		couples.iloc[np.where(donor_idx == 0)[0],:].\
			rename(columns=
					{i+'1': i+'Donor' for i in get_vars} | 
					{i+'2': i+'Recipient' for i in get_vars}),
		couples.iloc[np.where(donor_idx == 1)[0],:].\
			rename(columns=
					{i+'2': i+'Donor' for i in get_vars} | 
					{i+'1': i+'Recipient' for i in get_vars})]).\
		reset_index().\
		pipe(set_dates)
		



	# filters observations
	keep = (eval(args.filter)(couples))

	if np.isnan(keep).any():
		raise Exception('donor or recipient filter returned NA')
	
	couples = couples[keep]

	# and now, accumulate our time intervals
	# diff_round = difference between rounds, to ensure that
	# valid coupling persists across rounds
	# lag_int_dateRecipient = last recipient interview date
	# in this coupling with a valid HIV testing result
	# this allows for periods across recipient missed interim rounds
	couples = couples.\
		sort_values(by=['couple', 'round']).\
		reset_index(drop=True).\
		assign(
			diff_round = lambda k: 
				np.where(
					k['couple'].shift() == k['couple'],
					k['round'].diff(),
					np.nan),
			lag_int_dateRecipient = lambda k:		
				pd.to_datetime(pd.Series(np.where(k['couple'].shift() == k['couple'],
					pd.Series(np.where(k.finalhivRecipient.isnull() | k.int_dateRecipient.isnull(),
						np.nan,
						k.int_dateRecipient.astype(str))).shift().ffill(),
					np.nan))))

	# filter to only valid time periods
	# have to observe couple in each round
	# recipient has to be observed and have valid HIV test
	# recipient has to have been observed at start of follow-up to get time
	keep = (~couples.diff_round.isnull()) & \
		(couples.diff_round == 1.0) & \
		(~couples.int_dateRecipient.isnull() & ~couples.finalhivRecipient.isnull()) & \
		(~couples.lag_int_dateRecipient.isnull())
		
	couples = couples[keep]

	# add source round for lag_int_dateRecipient
	couples = couples.\
		merge(
			d[['study_id', 'round', 'int_date', 'exclusively_partnered', 'monogamous']].\
				rename(columns={
					'study_id': 'study_idRecipient',
					'round': 'lag_int_dateRecipient_round',
					'int_date': 'lag_int_dateRecipient'}),
			how='left',
			on=['study_idRecipient', 'lag_int_dateRecipient'])

	# and finally, add time duration and event indicator
	couples = couples.assign(
		dt = lambda k:
			(k.int_dateRecipient - k.lag_int_dateRecipient).dt.days/365.25,
		conversion = lambda k:
			k.finalhivRecipient == 'P')

	dirs = args.out.split('/')[:-1]
	mkdir = dirs[0]
	os.makedirs(mkdir, exist_ok=True)
	for sub_dir in dirs[1:]:
		mkdir += f'/{sub_dir}'
		os.makedirs(mkdir, exist_ok=True)
	
	couples.to_csv(args.out + '.tsv', sep='\t', index=None)

	print(f'{couples.conversion.sum()} seroconversions observed over {couples.dt.sum()} years of follow up contributed by {couples.couple.drop_duplicates().shape[0]} couples in {couples.shape[0]} discrete periods',
			file=sys.stderr)


if __name__ == '__main__':
    run()
