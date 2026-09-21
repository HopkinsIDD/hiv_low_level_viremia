import pandas as pd
import numpy as np
import argparse
import cmdstanpy
import glob
try:
	from scripts.utils import numeric_from_datetime, datetime_from_numeric, split_pad_col, import_config, get_ind_n_trans_bd, get_couple_prob_trans_bd, summarize_draws, get_trans_form
except:
	from utils import numeric_from_datetime, datetime_from_numeric, split_pad_col, import_config, get_ind_n_trans_bd, get_couple_prob_trans_bd, summarize_draws, get_trans_form


def format_small(x):
		x_fmt = f'{x:.1e}'.split('e')
		return(x_fmt[0] + '$\\\\times10^{' + str(int(x_fmt[1])) + '}$')


def num_to_word(nums):
		num_words = np.array(['zero', 'one', 'two', 'three', 'four', 'five', 'six', 'seven', 'eight', 'nine'])
		out = nums.astype(str)
		out[nums < 10] = num_words[nums[nums < 10]]
		return(out)


def get_rccs_stats(dat_file, config, round_digits=1):
	i_output = []
	dat = pd.read_csv(dat_file, sep='\t')
	# survey dates
	i_output.append(dat[['round', 'int_date']].dropna().\
		assign(int_date = lambda k: pd.to_datetime(k.int_date)).\
		groupby('round').int_date.quantile([0.001, 0.5, 0.999]).\
		reset_index().\
		assign(
			label = lambda k: k.level_1.map({0.001: 'min', 0.500: 'median', 0.999: 'max'}),
			int_date = lambda k: pd.to_datetime(k.int_date).dt.strftime('%B %-d, %Y')).\
		assign(
			var = lambda k: 'r'+np.where(k['round'] == 15.1, '15.1', k['round'].astype(int).astype(str))+'_'+k.label+'Date')\
		[['var', 'int_date']].values)
	# survey years
	i_output.append(dat[['round', 'int_date']].dropna().\
		assign(int_date = lambda k: pd.to_datetime(k.int_date)).\
		groupby('round').int_date.quantile([0.001, 0.5, 0.999]).\
		reset_index().\
		assign(
			label = lambda k: k.level_1.map({0.001: 'min', 0.500: 'median', 0.999: 'max'}),
			int_date = lambda k: pd.to_datetime(k.int_date).dt.strftime('%Y')).\
		assign(
			var = lambda k: 'r'+np.where(k['round'] == 15.1, '15.1', k['round'].astype(int).astype(str))+'_'+k.label+'Year')\
		[['var', 'int_date']].values)
	# number of participants
	i_output.append(['n_par', f"{dat.study_id.drop_duplicates().shape[0]:,}"])
	i_output.append(['n_par_r16_plus', f"{dat.query("round >= 16").study_id.drop_duplicates().shape[0]:,}"])	
	n_par_per_round = dat.groupby('round').size().reset_index().\
		assign(label = lambda k: 'n_r'+np.where(k['round'] == 15.1, '15.1', k['round'].astype(int).astype(str))+'_par')
	i_output.append(n_par_per_round.\
			assign(val = lambda k: [f'{i:,}' for i in k[0]])[['label', 'val']].values)
	# numer of participants with HIV
	i_output.append(['n_plhiv', f"{dat.query('finalhiv == "P"').study_id.drop_duplicates().shape[0]:,}"])
	n_plhiv_r16_plus = dat.query("round >= 16").query('finalhiv == "P"').study_id.drop_duplicates().shape[0]
	i_output.append(['n_plhiv_r16_plus', f"{n_plhiv_r16_plus:,}"])
	n_plhiv_per_round = dat.query('finalhiv == "P"').groupby('round').size().reset_index().\
		assign(label = lambda k: 'n_r'+np.where(k['round'] == 15.1, '15.1', k['round'].astype(int).astype(str))+'_plhiv')
	i_output.append(n_plhiv_per_round.\
			assign(val = lambda k: [f'{i:,}' for i in k[0]])[['label', 'val']].values)
	# number of participants with viremia
	n_viremic_per_round = dat.query('finalhiv == "P"').query('~copies.isnull() & (copies > 200)').\
		groupby('round').size().reset_index().\
		assign(label = lambda k: 'n_r'+np.where(k['round'] == 15.1, '15.1', k['round'].astype(int).astype(str))+'_viremic')
	i_output.append(n_viremic_per_round.\
			assign(val = lambda k: [f'{i:,}' for i in k[0]])[['label', 'val']].values)
	# number of participants with LLV
	n_llv = dat.query('finalhiv == "P"').query('~copies.isnull()').\
		query('(copies > 200) & (copies < 1000)').shape[0]
	i_output.append(['n_llv', f"{n_llv:,}"])
	n_llv_r16_plus = dat.query("round >= 16").query('finalhiv == "P"').query('~copies.isnull()').\
		query('(copies > 200) & (copies < 1000)').shape[0]
	i_output.append(['n_llv_r16_plus', f"{n_llv_r16_plus:,}"])
	p_llv_r16_plus = 100*n_llv_r16_plus/n_plhiv_r16_plus
	i_output.append(['p_llv_r16_plus', f"{round(p_llv_r16_plus,2)}\\\\%"])
	n_llv_per_round = dat.query('finalhiv == "P"').query('~copies.isnull() & (copies > 200) & (copies < 1000)').\
		groupby('round').size().reset_index().\
		assign(label = lambda k: 'n_r'+np.where(k['round'] == 15.1, '15.1', k['round'].astype(int).astype(str))+'_llv')
	i_output.append(n_llv_per_round.\
			assign(val = lambda k: [f'{i:,}' for i in k[0]])[['label', 'val']].values)
	n_missing_vl_r16_plus = dat.query('round >= 16').query('finalhiv == "P"').query('copies.isnull()').shape[0]
	p_missing_vl_r16_plus = n_missing_vl_r16_plus / dat.query('round >= 16').query('finalhiv == "P"').shape[0]
	i_output.append(['n_missing_finalhiv_r16_plus', dat.query('round >= 16').query('finalhiv.isnull()').shape[0]])
	i_output.append(['n_missing_vl_r16_plus', n_missing_vl_r16_plus])
	i_output.append(['p_missing_vl_r16_plus', f'{round(100*p_missing_vl_r16_plus, 1)}\\\\%'])
	'''
	p_llv_among_par_per_round = n_llv_per_round.merge(n_par_per_round, on='round', how='left').\
		assign(val = lambda k: np.round(100*k['0_x']/k['0_y'],2).astype(str) + '\\\\%',
			label = lambda k: k.label_x.str.replace('n_', 'p_').str.replace('_llv', '_llv_among_par'))
	i_output.append(p_llv_among_par_per_round[['label', 'val']].values)	
	p_llv_among_plhiv_per_round = n_llv_per_round.merge(n_plhiv_per_round, on='round', how='left').\
		assign(val = lambda k: round(100*k['0_x']/k['0_y'], 1).astype(str) + '\\\\%',
			label = lambda k: k.label_x.str.replace('n_', 'p_').str.replace('_llv', '_llv_among_plhiv'))
	i_output.append(p_llv_among_par_per_round[['label', 'val']].values)	
	p_llv_among_viremic_per_round = n_llv_per_round.merge(n_viremic_per_round, on='round', how='left').\
		assign(val = lambda k: round(100*k['0_x']/k['0_y'], 1).astype(str) + '\\\\%',
			label = lambda k: k.label_x.str.replace('n_', 'p_').str.replace('_llv', '_llv_among_viremic'))
	i_output.append(p_llv_among_viremic_per_round[['label', 'val']].values)	
	'''
	# number unique paired tx llv viral load measurements
	# only consecutive rounds after treatment, HIV positive, with VL data
	paired_vl = dat.query('(round > last_arvmedFalse)&(finalhiv == "P")&(~copies.isnull())').\
		sort_values(by=['study_id', 'round']).\
		assign(
			lag_int_date = lambda k: k.int_date.shift(),
			lag_round = lambda k: k['round'].shift(),
			lag_copies = lambda k: k.copies.shift()).\
		query('(study_id == study_id.shift()) & (round == (round.shift()+1))')\
		[['study_id', 'round', 'lag_round', 'int_date', 'lag_int_date', 'copies', 'lag_copies']]
	# number with Vls in each bin
	vl_bins = config['breaks_copies']
	def get_bin_props(vl_dat, vl_bins):
		vl_binned = np.digitize(vl_dat, vl_bins)
		bin_counts = np.unique(vl_binned, return_counts=True)
		bin_props = np.zeros(vl_binned.max()+1)
		bin_props[bin_counts[0]] = bin_counts[1]/vl_dat.shape[0]
		bin_props_labels = np.round(100*bin_props[1:], round_digits).astype(str) + '\\\\%'
		return(bin_counts, bin_props, bin_props_labels)
	for label, use_dat in [
			('llv', paired_vl.query('(lag_copies > 200) & (lag_copies < 1000)').\
						assign(min_round = lambda k: k.groupby('study_id')['round'].transform('min')).\
						query('round == min_round')),
			('supr', paired_vl.query('(lag_copies < 40)').\
						assign(min_round = lambda k: k.groupby('study_id')['round'].transform('min')).\
						query('round == min_round'))]:
		i_output.append([f'n_uniq_tx_{label}_vl_pairs', f'{use_dat.shape[0]:,}'])
		# median and 95th percentile of follow-up duration
		time_between_pairs_q = np.round(np.quantile(
			(pd.to_datetime(use_dat.int_date) - pd.to_datetime(use_dat.lag_int_date)).dt.days/365.25, 
			[0.025, 0.5, 0.975]), 2)
		i_output.append([f'q050_uniq_tx_{label}_vl_pairs_years',
			f"{time_between_pairs_q[1]:,}"])
		i_output.append([f'q025_uniq_tx_{label}_vl_pairs_years',
			f"{time_between_pairs_q[0]:,}"])
		i_output.append([f'q975_uniq_tx_{label}_vl_pairs_years',
			f"{time_between_pairs_q[2]:,}"])
		dat_bin_counts, dat_bin_props, dat_bin_props_labels = get_bin_props(use_dat.copies, vl_bins)
		i_output.extend([[
			f'n_uniq_tx_{label}_vl_pairs_{np.where(vl_bins[idx-1]<0, 0, vl_bins[idx-1]).astype(int)}_to_{'inf' if vl_bins[idx] > 1000 else str(int(vl_bins[idx]))}_copies',
				f'{dat_bin_counts[1][idx-1]:,}'] for idx in np.arange(1, dat_bin_counts[0].shape[0]+1)])
		i_output.extend([[
			f'p_uniq_tx_{label}_vl_pairs_{np.where(vl_bins[idx-1]<0, 0, vl_bins[idx-1]).astype(int)}_to_{'inf' if vl_bins[idx] > 1000 else str(int(vl_bins[idx]))}_copies',
				dat_bin_props_labels[idx-1]] for idx in np.arange(1, dat_bin_counts[0].shape[0]+1)])
	return(i_output)


'''
def get_period_stats(period_file, period_label, config, round_digits=2):
	i_output = []
	period_dat = pd.read_csv(period_file, sep='\t')
	# remove any with negative duration of follow-up
	# flag: why exist?
	period_dat = period_dat.query('dt> 0')
	# number of couples
	i_output.append([f'n_{period_label}_couples', 
		f"{period_dat[['study_idDonor', 'study_idRecipient']].drop_duplicates().shape[0]:,}"])
	# number of years of follow-up
	i_output.append([f'n_{period_label}_years', 
		f"{round(period_dat.dt.sum(), 0).astype(int):,}"])
	# median and 95th percentile of follow-up duration
	follow_up_q = np.round(np.quantile(period_dat.dt, [0.025, 0.5, 0.975]), round_digits)
	i_output.append([f'q050_{period_label}_years',
		f"{follow_up_q[1]:,}"])
	i_output.append([f'q025_{period_label}_years',
		f"{follow_up_q[0]:,}"])
	i_output.append([f'q975_{period_label}_years',
		f"{follow_up_q[2]:,}"])
	# number of seroconversions
	i_output.append([f'n_{period_label}_seroconversion', 
		f"{period_dat.conversion.sum():,}"])
	return(i_output)
'''


#get_merged_period_stats(args.mergedPeriodDat[idx], args.mergedPeriodDatLabels[idx], config, round_digits=1)
#period_file = args.mergedPeriodDat[0]
#period_label = args.mergedPeriodDatLabels[0]
def get_period_stats(period_file, period_label, config, round_digits=2):
	i_output = []
	period_dat = pd.read_csv(period_file, sep='\t')
	if 'dt_sum' not in period_dat.columns:
		# not merged with viral load data
		period_dat['dt_sum'] = period_dat.dt
		period_dat['int_dateRecipient'] = [numeric_from_datetime(i) for i in pd.to_datetime(period_dat['int_dateRecipient'])]
		period_dat['lag_int_dateRecipient'] = [numeric_from_datetime(i) for i in pd.to_datetime(period_dat['int_dateRecipient'])]
	# remove any with negative duration of follow-up
	# flag: why exist?
	period_dat = period_dat.query('dt_sum > 0')
	# number of couples
	i_output.append([f'n_{period_label}_couples', 
		f"{period_dat[['study_idDonor', 'study_idRecipient']].drop_duplicates().shape[0]:,}"])
	# total number of years of follow-up
	i_output.append([f'n_{period_label}_years', 
		f"{round(period_dat.dt_sum.sum(), 0).astype(int):,}"])
	# median and 95th percentile of follow-up duration per couple
	follow_up_q = np.round(np.quantile(
		period_dat.groupby(['study_idDonor', 'study_idRecipient']).dt_sum.sum(), [0, 0.5, 1.0]), round_digits)
	i_output.append([f'q050_{period_label}_years',
		f"{follow_up_q[1]:,}"])
	i_output.append([f'q000_{period_label}_years',
		f"{follow_up_q[0]:,}"])
	i_output.append([f'q100_{period_label}_years',
		f"{follow_up_q[2]:,}"])
	# median and 95th percentile of time between tests
	interval_q = np.round(np.quantile(
			(period_dat.int_dateRecipient - 
				period_dat.lag_int_dateRecipient).values,
			[0.000, 0.5, 1.00]),
		round_digits)
	i_output.append([f'q050_{period_label}_interval',
		f"{interval_q[1]:,}"])
	i_output.append([f'q000_{period_label}_interval',
		f"{interval_q[0]:,}"])
	i_output.append([f'q100_{period_label}_interval',
		f"{interval_q[2]:,}"])
	# number of seroconversions
	i_output.append([f'n_{period_label}_seroconversion', 
		f"{period_dat.conversion.sum():,}"])
	return(i_output)



def get_merged_period_stats(period_file, period_label, config, round_digits=2):
	i_output = []
	period_dat = pd.read_csv(period_file, sep='\t')
	# remove any with negative duration of follow-up
	# flag: why exist?
	period_dat = period_dat.query('dt_sum > 0')
	thresh = 7
	# how many intervals have index partner viral load w/in thresh of follow-up test
	donor_copies_date = split_pad_col(period_dat.donor_copies_date.values)
	n_within_thresh = (np.nanmin(np.abs(period_dat.int_dateRecipient.values[:,np.newaxis] - donor_copies_date), axis=1) <= thresh/365.25).sum()
	i_output.append([f'n_{period_label}_follow-up_vl_testing_interval_lte_{thresh}days', n_within_thresh])
	i_output.append([f'p_{period_label}_follow-up_vl_testing_interval_lte_{thresh}days', np.round(100*n_within_thresh/period_dat.shape[0], round_digits)])
	n_within_year = (np.nanmin(np.abs(period_dat.int_dateRecipient.values[:,np.newaxis] - donor_copies_date), axis=1) <= 1).sum()
	i_output.append([f'n_{period_label}_follow-up_vl_testing_interval_lte_1yr', n_within_year])
	i_output.append([f'p_{period_label}_follow-up_vl_testing_interval_lte_1yr', np.round(100*n_within_year/period_dat.shape[0], round_digits)])
	# maximum duration between index parter viral load and follow-up test
	i_output.append([f'q100_{period_label}_follow_up_vl_testing_interval', 
		np.round(
			np.nanmin(np.abs(period_dat.int_dateRecipient.values[:,np.newaxis] - donor_copies_date), axis=1).max(),
			round_digits)])
	period_dat[['copies1', 'copies2']] = split_pad_col(period_dat.donor_copies.values)
	period_dat[['copies1_round', 'copies2_round']] = split_pad_col(period_dat.donor_copies_round.values)
	llv = period_dat.query('((copies1 >200) & (copies1 < 1000)) | ((copies2 >200) & (copies2 < 1000))')
	# number of llv couples
	i_output.append([f'n_{period_label}_llv_couples', 
		f"{llv[['study_idDonor', 'study_idRecipient']].drop_duplicates().shape[0]:,}"])
	# number of llv years
	i_output.append([f'n_{period_label}_llv_years', 
		f"{round(llv.dt_sum.sum(), 0).astype(int):,}"])
	# number of llv seroconversions
	i_output.append([f'n_{period_label}_llv_seroconversion', 
		llv.conversion.sum()])
	# on treatment and with LLV at initiatial measurement
	tx = {True: 'treated', False: 'pre-treatment'}
	init_llv_tx = llv.assign(tx=lambda k: (k['round'] > k.last_arvmedFalseDonor).map(tx)).\
			query('(tx == "treated") & (copies1 <1000) & (copies1_round == lag_int_dateRecipient_round)')
	i_output.extend([
		[f'n_{period_label}_llv_tx_200-1000_to_any_couples', init_llv_tx.shape[0]],
		[f'n_{period_label}_llv_tx_200-1000_to_any_seroconversion', 
			init_llv_tx.conversion.sum()]])
	return(i_output)


# cohort
# all serodifferent couples
# all serodifferent monogamous, post-acute couples
# all serodifferent monogamous, post-acute seroconverting
#get_epi_table(args.dat, args.periodMetadat, args.mergedPeriodDat[0], config, round_digits=1)
#cohort_file = args.dat
#merged_period_file = args.mergedPeriodDat[0]
#all_period_file = args.periodDat[0]
#mono_period_file = args.periodDat[1]
def get_epi_table(
		cohort_file, 
		all_period_file, 
		mono_period_file,
		merged_period_file, 
		config,
		round_digits=1):
	def summarize_demo(df, group_cols):
		out = []
		# for categorical variables get most common response
		for col in ['sex', 'comm_type']:
			out.append(df[[*group_cols, col]].\
				groupby(group_cols)[col].value_counts().\
				sort_values(ascending=False).reset_index().groupby(group_cols).first().drop('count', axis=1))
		# for numeric get median
		for col in ['ageyrs']:
			out.append(df[[*group_cols, col]].\
				groupby(group_cols)[col].agg('median'))
		out = pd.concat(out, axis=1).reset_index()
		return(out)
	def merge_period_cohort_dat(period_dat, cohort_dat):
		# rough check of format of int_dateRecipient
		if period_dat.int_dateRecipient.dtype != np.float64:
			period_dat =  period_dat.\
			assign(int_dateRecipient = lambda k: 
				[numeric_from_datetime(i) for i in pd.to_datetime(k.int_dateRecipient)])
		# need to get age at end of interval
		out = period_dat.\
			rename(columns={'study_idDonor': 'study_id', 'int_dateRecipient':'int_date'}).\
			merge(cohort_dat[['study_id', 'int_date', 'sex', 'comm_type', 'ageyrs']].\
					query('~sex.isnull() & ~comm_type.isnull() & ~ageyrs.isnull()').\
					assign(int_date = lambda k: [numeric_from_datetime(i) for i in pd.to_datetime(k.int_date)]).\
					rename(columns={'int_date':'int_date_y'}),
				how='left',
				on='study_id').\
			assign(dt_dif = lambda k: k['int_date'] - k['int_date_y']).\
			assign(dt_dif_abs = lambda k: np.abs(k.dt_dif)).\
			assign(dt_dif_min = lambda k: k.groupby(['study_id', 'study_idRecipient', 'round']).\
				dt_dif_abs.transform('min')).\
			query('dt_dif_abs == dt_dif_min').\
			assign(ageyrs = lambda k: k.ageyrs + k.dt_dif)
		return(out)
	# get cohort stats
	cohort_dat = pd.read_csv(cohort_file, sep='\t').\
		query('finalhiv == "P"')
	dat_sum = [summarize_demo(cohort_dat, ['study_id']).\
		assign(type='cohort')]
	for period_dat, label in [
				(pd.read_csv(all_period_file,sep='\t'),
					'all_couples'), 
				(pd.read_csv(mono_period_file,sep='\t'),
					'mono_couples'), 
			 	(pd.read_csv(merged_period_file,sep='\t'), 
			 		'mono_noacute_vl_couples'),
			 	(pd.read_csv(merged_period_file,sep='\t').query('conversion == True'), 
			 		'mono_noacute_vl_converting_couples')]:
		period_dat = period_dat.\
			rename(columns={'round': 'round'})
		if 'int_date' in period_dat.columns:
			period_dat = period_dat.drop(['int_date'], axis=1)
		if 'dt_sum' in period_dat.columns:
			period_dat = period_dat.assign(dt_sum = lambda k: k.dt_sum.astype('float')).query('dt_sum > 0')
		elif 'dt' in period_dat.columns:
			period_dat = period_dat.assign(dt = lambda k: k.dt.astype('float')).query('dt > 0')
		dat_sum.append(
			summarize_demo(
				merge_period_cohort_dat(period_dat, cohort_dat),
				['study_id', 'study_idRecipient']
				).assign(type=label).drop(['study_idRecipient'], axis=1))
	# concatenate data types, assign age cat
	breaks = [14,24,34,60]
	labels=np.array(['15-24', '25-34', '35-60'])
	dat = pd.concat(dat_sum).\
		assign(age_cat = lambda k: labels[np.digitize(k.ageyrs, breaks, right=True)-1]).\
		drop(['ageyrs'], axis=1)
	# tabulate
	tab = []
	# howmany in each group total
	totals = {idx: i for idx, i in dat.groupby('type').size().items()}
	labels = {'sex': 'Sex', 'comm_type': 'Community type', 'age_cat': 'Age (years)',
		'cohort': 'Participants with HIV', 
		'all couples': 'SD couples',
		'mono_couples': 'MSD couples', 
		'mono_noacute_vl_couples': 'MSDVL couples',
		'mono_noacute_vl_converting_couples': 'MSDVLC couples',
		'var': ''}
	tab = [dat.groupby(['type']).size().reset_index(name='n').assign(var='', n=lambda k: k['n'].map("{:,d}".format)).\
				pivot(columns='type', values='n', index='var').\
				reset_index().rename(columns={'var': ''}).fillna('')]
	for var in ['sex', 'comm_type', 'age_cat']:
		tab.append(pd.concat([
			pd.DataFrame([[labels[var], '', '', '', '', '']], columns=['var', *totals.keys()]),
			dat.groupby(['type', var]).size().reset_index(name='n').\
				assign(n = lambda k: k.n.map("{:,d}".format) + ' (' + round(100*k['n']/np.array([totals[i] for i in k.type]), round_digits).astype(str) + '%)').\
				pivot(columns='type', values='n', index=var).\
				reset_index().rename(columns={var: 'var'}).fillna('').assign(var = lambda k: '    '+k['var'])]))
	tab = pd.concat(tab)[[
		'var',
		'cohort',
		'all_couples',
		'mono_couples',
		'mono_noacute_vl_couples',
		'mono_noacute_vl_converting_couples']].rename(columns=labels).fillna('')
	ftab = tab.copy()
	ftab.iloc[:,0] = ftab.iloc[:,0].str.replace('    ', '\\\\hspace{3mm}')
	# bold headings
	ftab.iloc[:,0] = np.where(~ftab.iloc[:,0].str.contains('hspace'), 
		'\\\\textbf{' + ftab.iloc[:,0] + '}', 
		ftab.iloc[:,0])
	tex = '\\\\' + 'begin{tabular}{p{1.8cm}C{2.2cm}C{3cm}C{3cm}C{3cm}C{3cm}}\n'
	tex += '\\\\toprule\n'
	tex += '&\\\\multicolumn{5}{c}{\\\\textbf{Participants with HIV}}' + '\\\\\\\\\n'
	tex += '\\\\cline{2-6}\n'
	tex += '&\\\\textbf{All}&\\\\textbf{Index in serodifferent couple}&\\\\textbf{Index in serodifferent, monogamous-partner couple}&\\\\textbf{Index in serodifferent, monogamous-partner couple in transmission analysis}&\\\\textbf{Index in serodifferent, monogamous, seroconverting-partner couple}' + '\\\\\\\\\n'
	#tex += '&&&\\\\multicolumn{3}{c}{\\\\textbf{w/ monogamous partners}}' + '\\\\\\\\\n'
	#tex += '\\\\cline{4-6}\n'
	#tex += '&&&&\\\\multicolumn{2}{c}{\\\\textbf{in transmisison analysis}}' + '\\\\\\\\\n'
	#tex += '\\\\cline{5-6}\n'
	#tex += '&\\\\textbf{Participants w/ HIV}&\\\\textbf{All}&\\\\textbf{All}&\\\\textbf{All}&\\\\textbf{w/ seroconverting partners}' + '\\\\\\\\\n'
	tex += '\\\\midrule\n'
	tex += '\\\\\\\\'.join(['&'.join(i) for i in ftab.values]) + '\\\\\\\\'
	tex += '\n\\\\bottomrule\n'
	tex += '\\\\end{tabular}'
	tex = tex.replace('%', '\\\\%')
	return(tab, tex)


#args.mergedPeriodDat[0], args.linkageDat
#period_file = args.mergedPeriodDat[0]
#linkage_file = args.linkageDat
def get_llv_table(period_file, linkage_file, config):
	period_dat = pd.read_csv(period_file, sep='\t')
	period_dat[['copies1', 'copies2']] = split_pad_col(period_dat.donor_copies).astype(float)
	llv = period_dat.query('((copies1 >200) & (copies1 < 1000)) | ((copies2 >200) & (copies2 < 1000))')
	llv[['copies_date1', 'copies_date2']] = \
		split_pad_col(llv.donor_copies_date).astype(float)
	# merge in linkage data
	linkage_dat = pd.read_csv(linkage_file, sep='\t').\
		assign(
			p24 = lambda k: np.where(
				k.p24_genolinked.isnull(),
				'unavailable', 
				np.where(
					k.p24_genolinked & k.p24_phylolinked,
					'linked',
					'unlinked')),
			gp41 = lambda k: np.where(
				k.gp41_genolinked.isnull(),
				'unavailable', 
				np.where(
					k.gp41_genolinked & k.gp41_phylolinked,
					'linked',
					'unlinked')))
	# manually change this linkage data for table
	linkage_dat = linkage_dat.assign(gp41 = lambda k: np.where(
		(k.study_idDonor == "R78014") & (k.study_idRecipient == "R04130"),
		"likely linked",
		k.gp41))
	llv = llv.merge(linkage_dat[['study_idDonor', 'study_idRecipient', 'p24', 'gp41']], 
		how='left', on=['study_idDonor', 'study_idRecipient'])
	# couple #donortx #donorcopies1 #donorcopies2 #partnerhiv1 #partnerhiv2 
	tex = '\\\\' + 'begin{tabular}{l' + "c"*9 + '}\n'
	tex += '\\\\toprule\n'
	tex += '&\\\\multicolumn{3}{c}{\\\\textbf{Index}}&&\\\\multicolumn{2}{c}{\\\\textbf{Partner}}&&\\\\multicolumn{2}{c}{\\\\textbf{Linkage}}' + '\\\\\\\\\n'
	tex += '\\\\cline{2-4}\\\\cline{9-10}\n'
	tex += '\\\\textbf{Couple}&\\\\textbf{Treatment}&\\\\multicolumn{2}{c}{\\\\textbf{Copies/mL (date)}}&&\\\\multicolumn{2}{c}{\\\\textbf{HIV seropositivity}}&&\\\\textbf{p24}&\\\\textbf{gp41}' + '\\\\\\\\\n'
	tex += '\\\\midrule\n'
	tx = {True: 'treated', False: 'pre-treatment'}
	for cdx, couple in llv.query('conversion == True').sort_values(by='int_dateRecipient').reset_index().iterrows():
		#tex += f'{couple.study_idDonor}$\\\\rightarrow${couple.study_idRecipient}&{tx[couple['round'] > couple.last_arvmedFalseDonor]}&'
		tex += f'LLV{int(cdx+1)}&{tx[couple['round'] > couple.last_arvmedFalseDonor]}&'
		tex += f'{"{:,}".format(int(np.round(couple.copies1, -2))) if couple.copies1 > 0 else "$<40$"} ({format(datetime_from_numeric(couple.copies_date1), "%Y-%m")})&'
		if ~np.isnan(couple.copies2):
			tex += f'{"{:,}".format(int(np.round(couple.copies2, -2))) if couple.copies2 > 0 else "$<40$"} ({format(datetime_from_numeric(couple.copies_date2), "%Y-%m")})'
		else:
			tex += 'unavailable'
		#else:
		#	tex += '&'
		tex += f'&&N ({format(datetime_from_numeric(couple.lag_int_dateRecipient), "%Y-%m")})'
		tex += f'&{couple.finalhivRecipient} ({format(datetime_from_numeric(couple.int_dateRecipient), "%Y-%m")})' 
		tex += f'&&{couple.p24}&{couple.gp41}' + '\\\\\\\\\n'
	tex += '\\\\bottomrule\n'
	tex += '\\\\end{tabular}'
	tex = tex.replace('%', '\\\\%')
	return(tex)
	

#linkage_dat_file = 'output/couples_linkage_sum.tsv'
def get_linkage_stats(linkage_dat_file, config, round_digits=2):
	linkage_dat = pd.read_csv(linkage_dat_file, sep='\t').\
		assign(p24_linked = lambda k: 
				np.where(
					k.p24_genolinked.isnull(),
					np.nan,
					k.p24_genolinked & k.p24_phylolinked),
			gp41_linked = lambda k: 
				np.where(
					k.gp41_genolinked.isnull(),
					np.nan,
					k.gp41_genolinked & k.gp41_phylolinked))
	i_output = []
	# genetic distance thresholds
	i_output.extend([[f'{region}_cutoff', 
		round(linkage_dat[f'{region}_j'].iloc[0], round_digits+1)] for region in ['p24', 'gp41']])
	# number with data in at least one segment
	i_output.append([
		f'n_couple_linkage_dat',
		linkage_dat.query('~p24_linked.isnull() | ~gp41_linked.isnull()').shape[0]])
	i_output.append([
		f'p_couple_linkage_dat',
		f'{round(100*\
				linkage_dat.query('~p24_linked.isnull() | ~gp41_linked.isnull()').shape[0]/linkage_dat.shape[0],
			round_digits)}\\\\%'])
	# number linked in at least one segment
	i_output.append([
		f'n_couple_linked',
		linkage_dat.query('(p24_linked == 1) | (gp41_linked == 1)').shape[0]])
	# proportion linked in at least one segment
	i_output.append([
		f'p_couple_linked',
		f'{round(100*\
			linkage_dat.query('(p24_linked == 1) | (gp41_linked == 1)').shape[0] / \
				linkage_dat.query('~p24_linked.isnull() | ~gp41_linked.isnull()').shape[0],
			round_digits)}\\\\%'])
	# number with data in both segments
	i_output.append([
		f'n_couple_both_linkage_dat',
		linkage_dat.query('~p24_linked.isnull() & ~gp41_linked.isnull()').shape[0]])
	# number with discordant results
	i_output.append([
		f'n_couple_discordant_linkage_dat',
		linkage_dat.query('~p24_linked.isnull() & ~gp41_linked.isnull()').query('p24_linked != gp41_linked').shape[0]])
	# proportion with discordant results
	i_output.append([
		f'p_couple_discordant_linkage_dat',
		100*(np.round(linkage_dat.query('~p24_linked.isnull() & ~gp41_linked.isnull()').query('p24_linked != gp41_linked').shape[0] / \
			linkage_dat.query('~p24_linked.isnull() & ~gp41_linked.isnull()').shape[0],1))])
	# number llv with data in at least one segment
	llv_linkage_dat = linkage_dat.query('((copies1 > 200) & (copies1 < 1000)) | ((copies2 > 200) & (copies2 < 1000))')
	i_output.append([
		f'n_llv_linkage_dat',
		num_to_word(np.array(llv_linkage_dat.query('~p24_linked.isnull() | ~gp41_linked.isnull()').shape[0]))])
	# number llv linked in at least one segment
	i_output.append([
		f'n_llv_linked',
		num_to_word(np.array(llv_linkage_dat.query('(p24_linked == 1) | (gp41_linked == 1)').shape[0]))])
	i_output.append([
		f'n_llv_likely_linked',
		num_to_word(np.array(llv_linkage_dat.\
			query('((p24_linked == 1) | (gp41_linked == 1)) | ((study_idDonor == "R78014")&(study_idRecipient=="R04130"))').\
			shape[0]))])
	# proportion llv linked in at least one segment
	i_output.append([
			f'p_llv_linked',
			f'{round(100*\
				llv_linkage_dat.query('(p24_linked == 1) | (gp41_linked == 1)').shape[0] / \
					llv_linkage_dat.query('~p24_linked.isnull() | ~gp41_linked.isnull()').shape[0],
				round_digits)}\\\\%'])
	# number llv linked in both segments
	i_output.append([
		f'n_llv_linked_both',
		num_to_word(np.array(llv_linkage_dat.query('(p24_linked == 1) & (gp41_linked == 1)').shape[0]))])
	# number llv linked in gp41 but not p24
	i_output.append([
		f'n_llv_linked_only_gp41',
		num_to_word(np.array(llv_linkage_dat.query('((p24_linked ==0) | p24_linked.isnull()) & (gp41_linked == 1)').shape[0]))])
	# genetic distance of closely linked llv couple
	get_donor = 'R78014'
	get_recipient = 'R04130'
	i_output.append([
		f'{get_donor}_{get_recipient}_gp41_d_min',
		round(llv_linkage_dat.query('(study_idDonor == @get_donor) & (study_idRecipient == @get_recipient)').gp41_d_min.iloc[0],
			round_digits+1)])
	return(i_output)



#idx = 5
#dr = fitDr[idx]
#label = args.fitDatLabels[idx]
#get_fit_stats(dr, label, config)
def get_fit_stats(dr, label, config, round_digits=2):
	# first value taken as reference for risk ratios
	focal_vl = np.array([dr['mu_0'].median(), 7, np.log10(1000), np.log10(200)])
	focal_vl_labels = np.array(['mu', 'max', '1000', '200'])
	r = get_ind_n_trans_bd(dr, t=100, v=focal_vl, form=dr.trans_form.iloc[0])
	r_sum = summarize_draws(r).\
		assign(
			v = focal_vl,
			label=focal_vl_labels)
	rr_sum = summarize_draws(r[1:,:] / r[0,:]).\
		assign(v = focal_vl[1:],
			label = focal_vl_labels[1:] + f'_v_{focal_vl_labels[0]}')
	i_output = []
	i_output.append([f'q50_{label}_mu0', f'{int(round((10**dr['mu_0']).median(), 0)):,}'])
	i_output.extend([[f'q{i}_{label}_{row.label}_trans_rate_per_100_year', 
		round(row[i], 1 if row[i] > 1 else 2)] for i in r_sum.columns[:-2] for rdx, row in r_sum.iterrows()])
	i_output.extend([[f'q{i}_{label}_{row.label}_trans_rr', 
		round(row[i],  1 if row[i] > 1 else 2)] for i in rr_sum.columns[:-2] for rdx, row in rr_sum.iterrows()])
	i_output.extend([[f'q{i}_{label}_{row.label}_trans_rr_percent', 
		f'{int(round(row[i]*100,  0))}\\\\%'] for i in [0.5] for rdx, row in rr_sum.iterrows()])
	# INTRA V. EXTRA RATES
	# extra partner rate
	beta_extra = dr[['beta_extra']].values.T*100
	beta_extra_sum = summarize_draws(beta_extra)
	i_output.extend(
		[[f'q{i}_{label}_extra_trans_rate_per_100_year', 
				round(beta_extra_sum[i].iloc[0], 1 if beta_extra_sum[i].iloc[0] > 1 else 2)] for 
			i in beta_extra_sum.columns])
	# point of parity
	# use intra-couple rate as extra-partner rate is only valid for couples
	v_full = np.linspace(0, 7, 500)
	rr_intra_extra_sum = summarize_draws(
			get_couple_prob_trans_bd(dr, t=1, v=v_full, form=dr.trans_form.iloc[0]) / \
				(1 - np.exp(-dr['beta_extra'].values*1))).\
		assign(v = v_full)
	i_output.append(
		[f'q0.5_{label}_intra_extra_rr_parity',
		(10**v_full[np.argmin(np.abs(1  - rr_intra_extra_sum[0.5]))]).astype(int)])
	return(i_output)


def get_sampled_risk_stats(sampled_risk_dat_file, config, round_digits=2):
	i_output = []
	sampled_risk_dat = pd.read_csv(sampled_risk_dat_file, sep='\t')
	# risk after one year
	i_output.extend(
		[[f'q{i['index']}_sampled_trans_risk_1yr', round(i[0], 1)] for idx, i in \
			summarize_draws([sampled_risk_dat.query("t == 365").trans_risk.values]).iloc[0].reset_index().iterrows()])
	# risk after one year
	i_output.extend(
		[[f'q{i['index']}_sampled_trans_risk_2yr', round(i[0], 1)] for idx, i in \
			summarize_draws([sampled_risk_dat.query("t == 730").trans_risk.values]).iloc[0].reset_index().iterrows()])
	# how much cumulative risk is accumulated in first 90 days
	t_compare = 90
	i_output.extend(
		[[f'q{i['index']}_prop_trans_risk_90', str(round(i[0], 1)) + '\\\\%'] for idx, i in \
		summarize_draws([sampled_risk_dat.query('(t == 365) | (t == @t_compare)').\
			pivot(index='idx', columns='t', values='trans_risk').\
			reset_index().\
			assign(rel_risk = lambda k: 100*k[90.0]/k[365.0]).rel_risk.values]).iloc[0].reset_index().iterrows()])
	# how much risk contributed by rebounded individuals at 1 yr
	rebounded_risk_dat = pd.\
		read_csv(sampled_risk_dat_file.replace('.tsv', '_rebounded_1yr.tsv'), sep='\t').query("t == 365")
	i_output.extend(
		[[f'q{i['index']}_prop_rebounded_1yr', str(int(round(i[0],0))) + '\\\\%'] for idx, i in \
		summarize_draws([100*rebounded_risk_dat.query("t == 365.0").p_rebound.values]).iloc[0].reset_index().iterrows()])
	i_output.extend(
		[[f'q{i['index']}_prop_trans_risk_rebounded_1yr', str(int(round(i[0],0))) + '\\\\%'] for idx, i in \
		summarize_draws([100*rebounded_risk_dat.query("t == t.max()").rel_trans_risk.values]).iloc[0].reset_index().iterrows()])
	return(i_output)


def get_fixed_risk_stats(fixed_risk_dat_file, config, round_digits=2):
	i_output = []
	for rdx, row in pd.read_csv(fixed_risk_dat_file, sep='\t').\
			query('t == t.max()').iterrows():
		i_output.extend(
			[[f'q{i['index']}_{row.label}_fixed_trans_risk', str(round(i[0],round_digits))] for idx, i in 
				row[['0.025', '0.25', '0.5', '0.75', '0.975']].reset_index(name=0).iterrows()])
	return(i_output)


#llv_prev_dat_file = args.llvPrevDat
def get_llv_prev_stats(llv_prev_dat_file, config):
	i_output = []
	llv_prev_dat = pd.read_csv(llv_prev_dat_file, sep='\t').\
		drop(['copies_cat'], axis=1).\
		melt(id_vars=['round', 'label', 'copies_cat_label'])
	i_output.extend([
		[f'q{i["variable"]}_r{int(i["round"])}_{i.copies_cat_label}_{i.label}', 
			f'{round(100*i.value, 2 if i.value < 0.01 else 1)}\\\\%'] for \
			idx, i in llv_prev_dat.iterrows()])
	return(i_output)


def get_vl_param_table(
		fits, labels, 
		params=[
			'fnr', 
			'mu_0', 'omega_mu', 
			'sigma_0', 'omega_sigma',],
		param_labels=[
			'$\\\\pi$', 
			'$\\\\mu_0$', '$\\\\Omega_\\\\mu$',
			'$\\\\sigma_0$', '$\\\\Omega_\\\\sigma$'],
		round_digits=2):
	def format_table_val(x, round_digits=2):
		if x >= 0.01:
			return(np.round(x, round_digits))
		else:
			return(format_small(x))
	print(params)
	tex = '\\\\' + 'begin{tabular}{l' + "C{1cm}"*(len(params)+1) + '}\n'
	tex += '\\\\toprule\n'
	tex += '&\\\\multicolumn{'+str(len(params))+'}{c}{\\\\textbf{Parameter (95\\\\% HPD)}}' + '\\\\\\\\\n'
	tex += '\\\\cline{2-'+str(len(params)+1)+'}\n'
	tex += '\\\\textbf{Model}&' + '&'.join(['\\\\textbf{'+p+'}' for p in param_labels]) + '\\\\\\\\\n'
	tex += '\\\\midrule\n'
	tex += ''.join([labels[jdx]+'&' + '&'.join([f'{format_table_val(i[0.5])} ({format_table_val(i[0.025])}, {format_table_val(i[0.975])})' for idx,i in summarize_draws(j[params].values.T).iterrows()]) + '\\\\\\\\\n' for jdx, j in enumerate(fits)])
	tex += '\\\\bottomrule\n'
	tex += '\\\\end{tabular}'
	return(tex)


def get_trans_param_table(fits, labels, round_digits=2):
	def format_table_val(x, round_digits=2):
		if x >= 0.01:
			return(np.round(x, round_digits))
		else:
			return(format_small(x))
	# assumes hill is last in fits list
	params = [
		'use_param3', 'use_param2', 'use_param1', 'beta_extra']
	param_labels = [
		'$c$', '$\\\\rho_\\\\text{est.}$', '$\\\\rho_\\\\text{prod.}$', '$\\\\beta_\\\\text{extra}$']
	sigmoid_param_labels = [
		'$v_{50}$', '$k$', '$\\\\beta_\\\\text{max}$']
	bd_fits = [i for idx, i in enumerate(fits) if 'sigmoid' not in labels[idx].lower()]
	bd_labels = [i for i in labels if 'sigmoid' not in i.lower()]
	sigmoid_fits = [i for idx, i in enumerate(fits) if 'sigmoid' in labels[idx].lower()]
	sigmoid_labels = [i for i in labels if 'sigmoid' in i.lower()]
	tex = '\\\\' + 'begin{tabular}{l' + "C{2.5cm}" + "C{2cm}"*(len(params)-1) + '}\n'
	tex += '\\\\toprule\n'
	tex += '&\\\\multicolumn{'+str(len(params))+'}{c}{\\\\textbf{Parameter (95\\\\% HPD)}}' + '\\\\\\\\\n'
	tex += '\\\\cline{2-'+str(len(params)+1)+'}\n'
	tex += '\\\\textbf{Model}&' + '&'.join(['\\\\textbf{'+p+'}' for p in param_labels]) + '\\\\\\\\\n'
	tex += '\\\\midrule\n'
	tex += ''.join([bd_labels[jdx]+'&' + '&'.join([f'{format_table_val(i[0.5])} ({format_table_val(i[0.025])}, {format_table_val(i[0.975])})' for 
			idx,i in summarize_draws(j[params].values.T).iterrows()]) + '\\\\\\\\\n' for jdx, j in enumerate(bd_fits)])
	tex += '\\\\midrule\n'
	tex += '&' + '&'.join(['\\\\textbf{'+p+'}' for p in sigmoid_param_labels]) + '\\\\\\\\\n'
	tex += '\\\\midrule\n'
	tex += ''.join([sigmoid_labels[jdx]+'&' + '&'.join([f'{format_table_val(i[0.5])} ({format_table_val(i[0.025])}, {format_table_val(i[0.975])})' for 
			idx,i in summarize_draws(j[params].values.T).iterrows()]) + '\\\\\\\\\n' for jdx, j in enumerate(sigmoid_fits)])
	tex += '\\\\bottomrule\n'
	tex += '\\\\end{tabular}'
	return(tex)



def get_trans_rate_table(fits, labels, v=[200, 1000], round_digits=2):
	def format_table_val(x, round_digits=2):
		if x >= 0.01:
			return(np.round(x, round_digits))
		else:
			return(format_small(x))
	# calculate risk per 100 person years across different 
	r_sum = [
		summarize_draws(
			get_ind_n_trans_bd(fit, t=100, v=np.log10(v), form=fit.trans_form.iloc[0])) for 
		fit in fits]
	tex = '\\\\' + 'begin{tabular}{l' + "c"*len(v) + '}\n'
	tex += '\\\\toprule\n'
	tex += '&\\\\multicolumn{'+str(len(v))+'}{c}{\\\\textbf{Viral load (copies/mL)}}' + '\\\\\\\\\n'
	tex += '\\\\cline{2-'+str(len(v)+1)+'}\n'
	tex += '\\\\textbf{Model}&' + '&'.join(['\\\\textbf{'+f'{iv:,}'+'}' for iv in v]) + '\\\\\\\\\n'
	tex += '\\\\midrule\n'
	tex += ''.join([labels[jdx]+'&' + '&'.join([f'{format_table_val(i[0.5])} ({format_table_val(i[0.025])}, {format_table_val(i[0.975])})' for idx, i in j.iterrows()]) + '\\\\\\\\\n' for jdx, j in enumerate(r_sum)])
	tex += '\\\\bottomrule\n'
	tex += '\\\\end{tabular}'
	return(tex)


def get_trans_rate_compare_table(fit_compare_file, labels, round_digits=2):
	fit_compare=pd.read_csv(fit_compare_file, sep='\t')\
		[['model', 'elpd_loo', 'elpd_diff', 'se_diff']]
	label_idx = [labels.index(j)+1 for j in fit_compare.model]
	fit_compare.model = \
		[f'{label_idx[idx]}. {i}' for idx, i in enumerate(
			fit_compare.model.\
			str.replace('_', ' ').\
			str.replace('_bd_', '_BD_').\
			str.replace('stC', 'stimated C').\
			str.replace('art tn', 'ART-TN').\
			str.replace('re-treatment', 're-treatment participants only').\
			str.replace('igmoid', 'igmoidal transmission').\
			str.replace('_', ' ').\
			str.capitalize())]
	fit_compare.columns = fit_compare.columns.str.replace("_loo", '').str.replace('_', ' ').str.upper().\
		str.replace('DIFF', '$\\\\Delta$').str.replace("VAR1", 'Model')
	tex = '\\\\' + 'begin{tabular}{l' + "c"*(fit_compare.shape[1]-1) + '}\n'
	tex += '\\\\toprule\n'
	tex += '&'.join(['\\\\textbf{'+i+'}' for i in fit_compare.columns]) + '\\\\\\\\\n'
	tex += '\\\\midrule\n'
	tex += '\\\\\\\\\n'.join([
		f'{i.iloc[0]}&' + '&'.join([str(round(j, round_digits)) for j in i.iloc[1:]]) for 
			idx,i in fit_compare.iterrows()])
	tex += '\\\\\\\\\n'
	tex += '\\\\bottomrule\n'
	tex += '\\\\end{tabular}'
	return(tex)



def run():
	parser = argparse.ArgumentParser()
	# input files
	parser.add_argument('--dat', 
	    help='path to selected rccs data')
	parser.add_argument('--llvPrevDat')
	args = parser.add_argument('--periodDat', nargs='+')
	args = parser.add_argument('--periodDatLabels', nargs='+')
	args = parser.add_argument('--mergedPeriodDat', nargs='+')
	args = parser.add_argument('--mergedPeriodDatLabels', nargs='+')
	args = parser.add_argument('--linkageDat')
	args = parser.add_argument('--fitDat', nargs='+')
	args = parser.add_argument('--fitDatLabels', nargs='+')
	args = parser.add_argument('--fitCompare')
	args = parser.add_argument('--sampledRiskDat')
	args = parser.add_argument('--fixedRiskDat')
	args = parser.add_argument('--sampledReboundRate')

	args = parser.add_argument('--config', default='config/config.csv')
	args = parser.parse_args()

	'''
	args.dat = 'data/not_shared/RCCSdata_R001_R020_VOIs_clean.tsv'
	args.llvPrevDat = 'output/RCCSdata_R016_R020_est_prev_llv.tsv'
	args.periodDat = ['data/not_shared/donor-all_recipient-all.tsv', 'data/not_shared/donor-all_recipient-exclusivelymonogamous.tsv']
	args.periodDatLabels = ['donor-all_recipient-all', 'donor-all_recipient-exclusivelymonogamous']
	args.mergedPeriodDat = ['data/not_shared/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d.tsv']
	args.mergedPeriodDatLabels = ['merged_donor-noacute_recipient-exclusivelymonogamous']

	args.linkageDat = 'output/couples_linkage_sum.tsv'
	args.fitDat = [
		'output/ve_constrained',
		'output/ve_constrained_sex',
		'output/ve_constrained_estC',
		'output/ve_artificially_constrained',
		'output/ve_unconstrained',
		'output/ve_constrained_0_1',
		'output/ve_constrained_39_40',
		'output/ve_nosaturation',
		'output/ve_pretreatment',
		'output/sigmoid_constrained']
	args.fitDatLabels = [
		'main_analysis',
		'sex_predictor',
		'estC',
		'artificially_constrained',
		'no_genetic_constraints',
		'low_art_tn_assumption',
		'high_art_tn_assumption',
		'no_saturation',
		'pre-treatment',
		'sigmoid']

	args.fitCompare = 'output/transmission_model_selection.tsv'
	args.sampledRiskDat = 'output/simulated_trajectory_risk.tsv'
	args.fixedRiskDat = 'output/fixed_simulated_trajectory_risk_sum.tsv'
	args.sampledReboundRate = 'output/rccs_r1_r20_vl_pairs_treated_uniq_llv_fit_mcmc_all.tsv'
	'''
	
	config = import_config(args.config)

	fitDr = [pd.read_csv(i+'/select_draws.tsv', sep='\t').assign(trans_form = get_trans_form(i)) for i in args.fitDat]
	
	output = []

	#### STUDY POPULATION STATS ####
	output.extend(get_rccs_stats(args.dat, config))

	_ = [output.extend(get_period_stats(args.periodDat[idx], args.periodDatLabels[idx], config, round_digits=1)) for 
		idx in range(len(args.periodDat))]

	_ = [output.extend(
			get_period_stats(args.mergedPeriodDat[idx], args.mergedPeriodDatLabels[idx], config, round_digits=1) + \
				get_merged_period_stats(args.mergedPeriodDat[idx], args.mergedPeriodDatLabels[idx], config, round_digits=1))
			for idx in range(len(args.mergedPeriodDat))]


	output.append([
		'demographic_table',
		get_epi_table(args.dat, 
			args.periodDat[[idx for idx, i in enumerate(args.periodDatLabels) if i == "donor-all_recipient-all"][0]],
			args.periodDat[[idx for idx, i in enumerate(args.periodDatLabels) if i == "donor-all_recipient-exclusivelymonogamous"][0]],
			args.mergedPeriodDat[0], config, round_digits=1)[1]])
		


	#### DESCRIPTIVE LLV STATS ####
	output.extend(get_llv_prev_stats(args.llvPrevDat, config))
	_ = output.extend(get_linkage_stats(args.linkageDat, config, round_digits=1))

	output.append([
		'llv_table',
		get_llv_table(args.mergedPeriodDat[0], args.linkageDat, config)])



	#### BAYESIAN MODEL FIT STATS ###
	fit_param_table = []
	fit_param_table_labels = []
	for idx, i in enumerate(fitDr):
		if 'mean_use_param1' in i.columns:
			for lab, col in [('Mean', 'mean'), ('Female', 'sex1'), ('Male', 'sex2')]:
				fit_param_table.append(
					i.drop([j for j in ['use_param1', 'use_param2', 'use_param3'] if j in i.columns], axis=1).\
						rename(columns={
							f'{col}_use_param1': 'use_param1', 
							f'{col}_use_param2': 'use_param2', 
							f'{col}_use_param3': 'use_param3'}))
				fit_param_table_labels.append(
					str(idx+1)+'. '+
						args.fitDatLabels[idx][0].upper() + 
						args.fitDatLabels[idx][1:].\
						replace('_bd_', '_BD_').\
						replace('stC', 'stimated C').\
						replace('art_tn', 'ART-TN').\
						replace('re-treatment', 're-treatment participants only').\
						replace('igmoid', 'igmoidal transmission').\
						replace('_', ' ') + ' ' + 
						lab)
				output.extend(
					get_fit_stats(fit_param_table[-1], args.fitDatLabels[idx] + '_' + lab.lower(), config))
		else:
			fit_param_table.append(i)
			fit_param_table_labels.append(
				str(idx+1)+'. '+
						args.fitDatLabels[idx][0].upper() + 
						args.fitDatLabels[idx][1:].\
						replace('_bd_', '_BD_').\
						replace('stC', 'stimated c').\
						replace('art_tn', 'ART-TN').\
						replace('re-treatment', 're-treatment participants only').\
						replace('igmoid', 'igmoidal transmission').\
						replace('_', ' '))
			output.extend(
				get_fit_stats(
					fit_param_table[-1], args.fitDatLabels[idx], config))


	output.append([
		'fit_no_art_vl_param_table',
		get_vl_param_table(
			[i for idx,i in enumerate(fit_param_table) if 'male' not in fit_param_table_labels[idx].lower()], 
			[i.replace('mean', '') for i in fit_param_table_labels if 'male' not in i.lower()])])


	output.append([
		'fit_art_vl_param_table',
		get_vl_param_table(
			[i for idx,i in enumerate(fit_param_table) if 'male' not in fit_param_table_labels[idx].lower()], 
			[i.replace('mean', '') for i in fit_param_table_labels if 'male' not in i.lower()],
			params=['txr', 'theta_un_tx', 'mu_un_tx', 'sigma_un_tx', 'theta_tx', 'mu_tx', 'sigma_tx'],
			param_labels=['$\\\\gamma$', '$\\\\theta_\\\\text{unART}$', '$\\\\mu_\\\\text{unART}$', '$\\\\sigma_\\\\text{unART}$',
				'$\\\\theta_\\\\text{ART}$', '$\\\\mu_\\\\text{ART}$', '$\\\\sigma_\\\\text{ART}$'])])


	output.append([
		'fit_trans_param_table',
		get_trans_param_table(
			fit_param_table,
			fit_param_table_labels)])


	output.append([
		'fit_trans_rate_table',
		get_trans_rate_table(
			fit_param_table,
			fit_param_table_labels,
			np.array([200, 1000, np.round(10**fitDr[0]['mu_0'].median()).astype(int)]))])


	output.append([
		'fit_trans_rate_compare_table',
		get_trans_rate_compare_table(args.fitCompare, args.fitDatLabels)])


	#### SIMULATION MODEL FIT STATS ###
	_ = output.extend(get_sampled_risk_stats(args.sampledRiskDat, config, round_digits=2))


	_ = output.extend(get_fixed_risk_stats(args.fixedRiskDat, config, round_digits=2))


	_ = output.extend([[f'q{i['index']}_omega_cease_per_year', np.round(i.values[1], 2)] for idx, i in \
		summarize_draws([pd.read_csv(args.sampledReboundRate, header=None).values[:,0]*365]).\
			iloc[0].reset_index().iterrows()])


	_ = output.extend([[f'q{i['index']}_risk_failure_per_year', np.round(i.values[1], 0).astype(int)] for idx, i in \
		summarize_draws([100*(1 - np.exp(-pd.read_csv(args.sampledReboundRate, header=None).values[:,0]*365))]).\
			iloc[0].reset_index().iterrows()])


	# save
	pd.DataFrame(np.vstack(output),
		columns=['var', 'val']).\
		to_csv('output/statistics.csv', sep=',', index=None)

	pd.DataFrame(np.vstack(output),
		columns=['var', 'val']).\
		to_csv('output/statistics.scsv', sep=';', index=None)



if __name__ == '__main__':
    run()


'''
def get_bd_param_table(fits, labels, round_digits=2):
	def format_table_val(x, round_digits=2):
		if x >= 0.01:
			return(np.round(x, round_digits))
		else:
			return(format_small(x))
	params = [
		'fnr', 'txr',
		'mu_0', 'sigma_mu', 
		'sigma_0', 'sigma_sigma',
		'prev_tx_theta', 'mu_tx', 'sigma_tx',
		'use_param3', 'use_param2', 'use_param3']
	param_labels = [
		'$\\\\pi$', '$\\\\gamma$',
		'$\\\\mu_0$', '$\\\\sigma_\\\\mu$',
		'$\\\\sigma_0$', '$\\\\sigma_\\\\sigma$',
		'$\\\\phi$', '$\\\\mu_q$', '$\\\\sigma_q$',
		'$c$', '$\\\\rho_\\\\text{est.}$', '$\\\\rho_\\\\text{prod.}$']
	param_sum = [summarize_draws(j[params].values.T) for j in fits if j.trans_form.iloc[0] == 'bd']
	fit_labels = [i for idx, i in enumerate(labels) if fits[idx].trans_form.iloc[0] == 'bd']
	tex = '\\\\' + 'begin{tabular}{l' + "C{3cm}"*len(fit_labels) + '}\n'
	tex += '\\\\toprule\n'
	tex += '&\\\\multicolumn{'+str(len(fit_labels))+'}{c}{\\\\textbf{Model}}' + '\\\\\\\\\n'
	tex += '\\\\cline{2-'+str(len(fit_labels)+1)+'}\n'
	tex += '\\\\textbf{Parameter}&' + '&'.join(['\\\\textbf{'+p+'}' for p in fit_labels]) + '\\\\\\\\\n'
	tex += '\\\\midrule\n'
	tex += ''.join([param_labels[i]+'&' + '&'.join([f'{format_table_val(j.loc[i,0.5])} ({format_table_val(j.loc[i,0.025])}, \
		{format_table_val(j.loc[i,0.975])})' for jdx, j in enumerate(param_sum)]) + '\\\\\\\\\n' for i in param_sum[0].index])
	tex += '\\\\bottomrule\n'
	tex += '\\\\end{tabular}'
	return(tex)
'''
'''
d = pd.read_csv('output/couples_linkage_sum.tsv', sep='\t')
d[['study_idDonor', 'study_idRecipient', 'copies1', 'copies2', 'p24_genolinked', 'p24_phylolinked', 'gp41_genolinked', 'gp41_phylolinked']].\
	sort_values(by=['copies1', 'copies2']).\
	query('((copies1 > 200) & (copies1 < 1000)) | ((copies2 > 200)&(copies2 < 1000))').\
	assign(p24 = lambda k: k.p24_genolinked * k.p24_phylolinked).\
	assign(gp41 = lambda k: k.gp41_genolinked*k.gp41_phylolinked)\
	[['study_idDonor', 'study_idRecipient', 'copies1', 'copies2', 'p24', 'gp41']]

	.\
	query('(~p24.isnull()) & (~gp41.isnull())').\
	assign(tot = lambda k: k.p24 + k.gp41).\
	query('tot == 1')
	assign(gp41 = lambda k: (k.gp41_genolinked == True) & (k.gp41_phylolinked == True)).\
	query('(p24 & !)')
'''

'''

d = pd.read_csv('data/internal/RCCSdata_R001_R019_VOIs_clean_INTERNAL.tsv.gz', sep='\t', compression='gzip')
d.query('study_id == "R41593"')[['study_id', 'study_id_orig', 'round', 'finalhiv', 'copies']]
d.query('study_id == "R86742"')[['study_id', 'study_id_orig', 'round', 'finalhiv', 'copies']]


p_d = pd.read_csv('output/donor-all_recipient-all.tsv', sep='\t')
'''
