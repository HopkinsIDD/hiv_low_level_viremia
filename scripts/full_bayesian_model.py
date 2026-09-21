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
try:
	from scripts.utils import split_col
except:
	from utils import split_col


def run():
	parser = argparse.ArgumentParser()
	parser.add_argument('--copiesDat', 
	    help='path to rccs viral load data')
	parser.add_argument('--periodDat',
		help='path to rccs contact period data')
	parser.add_argument('--stan', 
	    help='path to stan file')
	parser.add_argument('--transVL',
		choices=['prox', 'prox_pt'],
		default='prox')
	parser.add_argument('--transForm',
		choices=['ve', 'sigmoid'],
		default='ve')
	parser.add_argument('--vSuprMin',
		help='imputed BD lower limit',
		type=float,
		default=0.0)
	parser.add_argument('--vSuprMax',
		help='imputed BD upper limit',
		type=float,
		default=40.0)
	parser.add_argument('--vWin',
		type=float,
		default=4,
		help='size of viral load window when using proximal viral load measurements')
	parser.add_argument('--linkedPeriods',
		help='optional input file outlining which couples are linked')
	parser.add_argument('--allowSaturation',
		help='whether to allow saturation of the transmission risk',
		type=int,
		default=1,
		choices=[0, 1])
	parser.add_argument('--c', type=float, default=1)
	parser.add_argument('--label', default="full_bayesian_model")
	args = parser.parse_args()
	print(args, file=sys.stderr)
	#args.copiesDat = 'data/internal/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_d.tsv'
	#args.periodDat = 'data/internal/prox4.0_donor-noacute_recipient-exclusivelymonogamous_RCCSdata_R001_R019_all_copies_merged_p_d.tsv'
	#args.copiesDat = 'data/MOCK_copies.tsv'
	#args.periodDat = 'data/MOCK_periods.tsv'
	#args.label = 'MOCK'
	#args.stan = 'stan/vl_trans.stan'
	#args.linkedPeriods = 'output/couples_linkage_sum_simple.tsv'
	#args.label = 'test'
	keep_d = pd.read_csv(args.copiesDat, sep='\t')	
	p_d = pd.read_csv(args.periodDat, sep='\t')
	# convert columns to arrays
	keep_d['prev_putative_tx'] = split_col(keep_d['prev_putative_tx'], int)

	for col in ['dt', 'donor_copies']:
		if col in p_d.columns:
			p_d[col] = split_col(p_d[col])

	for col in ['donor_copies_round', 'donor_copies_obs_idx']:
		if col in p_d.columns:
			p_d[col] = split_col(p_d[col], int)
		

	q = max(i.shape[0] for i in p_d['donor_copies_obs_idx'])
	print(f'Fitting VL model to {keep_d.query('tx == 0').shape[0]} pre-treatment viral load measurements taken from {keep_d.query('tx == 0').ind_idx.drop_duplicates().shape[0]} individuals of which {keep_d.query('tx == 0').query('copies > 0.0').shape[0]} are above the LOD.',
		file=sys.stderr)
	if 'couple_idx' in p_d.columns:
		print(f'Fitting transmission rate model to {p_d.conversion.sum()} seroconversions accumulated over {np.hstack(p_d.dt.values).sum()} couple-years of follow-up during {p_d.shape[0]} coupling periods among {p_d.couple_idx.drop_duplicates().shape[0]} couples',
			file=sys.stderr)
	else:
		print(f'Fitting transmission rate model to {p_d.conversion.sum()} seroconversions accumulated over {np.hstack(p_d.dt.values).sum()} couple-years of follow-up during {p_d.shape[0]} coupling periods among {p_d[['study_idDonor', 'study_idRecipient']].drop_duplicates().shape[0]} couples',
			file=sys.stderr)
	print(f'Using the proximal model for the relationship between viral load and transmission rate',
		file=sys.stderr)
	print(f'Assuming that viral load measurements are valid for a {args.vWin} year window around the observed date',
		file=sys.stderr)

	if args.linkedPeriods:
		constrained = pd.read_csv(args.linkedPeriods, sep='\t')
		if ('couple_idx' in p_d.columns) & ('couple_idx' in constrained.columns):
			p_d = p_d.merge(constrained, how='left', on='couple_idx')
		else:
			p_d = p_d.merge(constrained, how='left', on=['study_idDonor', 'study_idRecipient'])

	# list of unique LODs and assign index
	LODs = {i:idx+1 for idx, i in enumerate(keep_d['lod'].drop_duplicates().sort_values().values)}
	# we pad unobserved period viral loads with 0.1 to distinguish them
	# from BD observations which are indicated with a 0
	# since t = 0 they drop out of likelihood regardless of viral load
	stan_dat = {
		"N_ind": keep_d.ind_idx.max(),
		"N": keep_d.shape[0],
		"ind_idx": keep_d.ind_idx,
		"v": np.where(keep_d.copies == 0, 1, keep_d.copies),
		"L": len(LODs.keys()),
		"lod": list(LODs.keys()),
		"lod_idx": keep_d.lod.map(LODs),
		"tx_avail": keep_d.tx_avail.values,
		"prev_putative_tx": np.vstack([np.pad(i, (0, len(LODs)-len(i)), constant_values=0) for i in keep_d.prev_putative_tx]),
		"tx": keep_d.tx,
		"P": p_d.shape[0],
		"Q": q,
		"q": [len(i.donor_copies_obs_idx) for idx, i in p_d.iterrows()],
		"period_donor_v_idx": np.vstack([np.pad(i, (0, q-len(i)), constant_values=-1) for i in p_d.donor_copies_obs_idx]),
		"period_t": np.vstack([np.pad(i, (0, q-len(i)), constant_values=-1) for i in p_d.dt]),
		"sex": (1 + 1*(keep_d.sex.values == 'M')),
		"period_y": p_d.conversion.astype(int),
		"sample_from_posterior": 1,
		"v_supr_min": args.vSuprMin,
		"v_supr_max": args.vSuprMax,
		"trans_vl": 1 if args.transVL == 'prox' else 2,
		"trans_form": {'ve': 1, 'sigmoid': 2, 'gamma': 3}[args.transForm],
		"constrain_prob_intra_given_transmission": 1 if args.linkedPeriods else 0,
		"N_constrained": np.where(~p_d.linked.isnull() & p_d.conversion)[0].shape[0] if args.linkedPeriods else 0,
		"period_constrain_idx": np.where(~p_d.linked.isnull() & p_d.conversion)[0]+1 if args.linkedPeriods else [],
		"N_constrained_linked": p_d.query('~linked.isnull() & conversion == True').linked.sum() if args.linkedPeriods else 1,
		"allow_saturation": args.allowSaturation,
		"c": args.c
		}

	# compile model
	m = CmdStanModel(stan_file=args.stan)

	#m = CmdStanModel(stan_file='stan/pt_vl_fnr_tx_all_bdMM_constrained.stan')
	# make output dir
	os.makedirs('output', exist_ok=True)
	#out_dir = f'output/{".".join(args.stan.split("/")[-1].split(".")[:-1])}{args.transVL}{args.vWin}{args.imputedBdLl}{args.imputedBdUl}{"Constrained" if args.linkedPeriods else ""}_{".".join(args.copiesDat.split("/")[-1].split(".")[:-1])}_{".".join(args.periodDat.split("/")[-1].split(".")[:-1])}'
	out_dir = 'output/' + args.label
	os.makedirs(out_dir, exist_ok=True)	
	# removes old fits
	_ = [os.remove(i) for i in glob.glob(out_dir+'/*')]
	# save simplified period data
	#print(out_dir + '/periods.tsv')
	#p_d.\
	#	to_csv(out_dir + '/periods.tsv', sep='\t', index=False)
	# fit model
	init_mu = np.log10(stan_dat['v'][stan_dat['v'] > 1]).mean()
	init_std = np.log10(stan_dat['v'][stan_dat['v'] > 1]).std()
	#iter_warmup=3000,
	#iter_sampling=3000,
	'''
	f = m.sample(data = stan_dat, 
		chains=4, 
		parallel_chains=1,
		iter_warmup=2000,
		iter_sampling=2000,
		save_warmup=False,
		show_console=True,
		output_dir=out_dir,
		inits={
		'mu_0': init_mu,
		'omega_mu': init_mu/10,
		'eta_mu_i': np.repeat(0, stan_dat['N_ind']),
		'sigma_0': 2,
		'omega_sigma': 2/10,
		'eta_mu_i': np.repeat(0, stan_dat['N_ind']),
		'mu_tx': 4,
		'sigma_tx': 1,
		'mu_un_tx': 4,
		'sigma_un_tx': 2,
		'logit_fnr': np.log(0.05/(1-0.05)),
		'logit_txr': np.log(0.05/(1-0.05)),
		'logit_theta_un_tx': np.log(0.75/(1-0.75)),
		'logit_theta_tx': np.log(0.9/(1-0.9)),
		"beta_extra": 0.005,
		"beta_tx": 0.005,
		'est_param1': {'ve':-5, 'sigmoid': 0.1, 'gamma':-10}[args.transForm],
		'est_param2': {'ve':-5, 'sigmoid': 1, 'gamma':-10}[args.transForm],
		'est_param3': {'ve':-5, 'sigmoid': 4, 'gamma':-5}[args.transForm]})
		'''
	f = m.sample(data = stan_dat, 
		chains=4, 
		parallel_chains=4,
		iter_warmup=2500,
		iter_sampling=2500,
		save_warmup=False,
		show_console=True,
		output_dir=out_dir,
		inits={
		'mu_0': init_mu,
		'omega_mu': init_mu/10,
		'eta_mu_i': np.repeat(0, stan_dat['N_ind']),
		'sigma_0': 2,
		'omega_sigma': 2/10,
		'eta_mu_i': np.repeat(0, stan_dat['N_ind']),
		'mu_tx': 4,
		'sigma_tx': 1,
		'mu_un_tx': 4,
		'sigma_un_tx': 2,
		'logit_fnr': np.log(0.05/(1-0.05)),
		'logit_txr': np.log(0.05/(1-0.05)),
		'logit_theta_un_tx': np.log(0.75/(1-0.75)),
		'logit_theta_tx': np.log(0.9/(1-0.9)),
		"beta_extra": 0.005,
		"beta_tx": 0.005,
		'est_param1': {'ve':-5, 'sigmoid': 0.1, 'gamma':-10}[args.transForm],
		'est_param2': {'ve':-5, 'sigmoid': 1, 'gamma':-10}[args.transForm],
		'est_param3': {'ve':-5, 'sigmoid': 4, 'gamma':-5}[args.transForm]})

	'''
	out_dir = 'output/ve_constrained_oldGOOD'
	f = cmdstanpy.from_csv(
			glob.glob(out_dir + '/*.csv'))
	'''
	dr = f.draws_pd()	
	if args.transForm == 'gamma':
		dr['use_param1_inv'] = 1/dr['use_param1']

	# get select columns
	select_columns = [
		'mu_0', 'omega_mu', 'sigma_0', 'omega_sigma',
		'logit_fnr', 'fnr',
		'logit_txr', 'txr',
		'theta_tx', 'mu_tx', 'sigma_tx',
		'theta_un_tx', 'mu_un_tx', 'sigma_un_tx',
		'est_param1', 'est_param2', 'est_param3',
		'use_param1', 'use_param2', 'use_param3',
		'mean_use_param1', 'mean_use_param2', 'mean_use_param3',
		'sex1_use_param1', 'sex1_use_param2', 'sex1_use_param3',
		'sex2_use_param1', 'sex2_use_param2', 'sex2_use_param3',
		'est_param1_coeffs', 'est_param2_coeffs',
		'beta_extra', 'beta_tx']
	dr_select = dr[[i for i in select_columns if i in dr.columns]]
	dr_select.to_csv(out_dir + '/select_draws.tsv', sep='\t', index=False)

	dr[[i for i in dr.columns if 'log_prob_no_transmission' in i]].\
		to_csv(out_dir + '/log_prob_no_transmission.tsv', sep='\t', index=False)
	#### SUMMARIZE PARAMETERS ####
	params = [
		'lp__',
		'mu_0', 'omega_mu', 'sigma_0', 'omega_sigma',
		'fnr',
		'txr',
		'theta_un_tx', 'mu_un_tx', 'sigma_un_tx',
		'theta_tx', 'mu_tx', 'sigma_tx',
		'use_param1', 'use_param2', 'use_param3', 'use_param1_inv',
		'est_param1_coeffs', 'est_param2_coeffs',
		'beta_extra', 'beta_tx']


	dr_sum = dr[[i for i in params if i in dr.columns]].\
		assign(it = lambda k: np.arange(k.shape[0])).\
		melt(id_vars='it').\
		drop('it', axis=1).\
		groupby('variable').\
		quantile(q=[0.025, 0.25, 0.5, 0.75, 0.975]).\
		reset_index().\
		pivot(index='variable', columns='level_1', values='value').\
		loc[[i for i in params if i in dr.columns],:]

	dr_sum.to_csv(out_dir + '/sum.tsv',
		sep='\t')

	# print to log
	print(dr_sum.to_string(),
		file=sys.stderr)

	#### SUMMARIZE INDIVIDUAL-LEVEL VIRAL LOADS ####
	(lambda z: z.agg(mu='mean',sigma=lambda k: np.std(k)).\
		merge(z.quantile(q=[0.025, 0.25, 0.5, 0.75, 0.975]),
			left_index=True,
			right_index=True,
			how='left'))(
			dr[[i for i in dr.columns if ((('mu_i' in i) | ('sigma_i' in i)) & ('eta' not in i))]].\
				assign(iter = lambda k: range(k.shape[0])).\
				melt(id_vars='iter').\
				groupby('variable')['value']).\
				reset_index().\
				pivot(index=['variable', 'mu', 'sigma'], columns='level_1', values='value').\
				reset_index().\
				assign(ind_idx = lambda k: [int(i.split('[')[1].replace(']', '')) for i in k.variable]).\
				merge(keep_d[[i for i in keep_d.columns if i in ['ind_idx', 'study_id']]].drop_duplicates(),
					how='left',
					on='ind_idx').\
				sort_values(by='ind_idx').\
				to_csv(out_dir + '/vl_sum.tsv',
					sep='\t',
					index=None)
			
	#### PAIRS PLOT ####
	'''
	pairs_params = [
		'logit_fnr',
		'mu_0', 'omega_mu', 'sigma_0', 'omega_sigma',
		'logit_txr',
		'logit_prev_tx_bdr', 'mu_tx', 'sigma_tx',
		'est_param1', 'est_param2', 'est_param3' if args.c == -1 else None, 
		'est_param1_coeffs[1]', 'est_param1_coeffs[2]',
		'est_param2_coeffs[1]', 'est_param2_coeffs[2]',
		'beta_extra',
		'beta_tx' if args.transVL == 'prox_pt' else None]
	'''
	pairs_params = [
		'mu_0', 'omega_mu', 'sigma_0', 'omega_sigma',
		'logit_fnr',
		'logit_txr',
		'logit_theta_un_tx', 'mu_un_tx', 'sigma_un_tx',
		'logit_theta_tx', 'mu_tx', 'sigma_tx',
		'est_param1', 'est_param2', 'est_param3' if args.c == -1 else None, 
		'est_param1_coeffs[1]', 'est_param1_coeffs[2]',
		'est_param2_coeffs[1]', 'est_param2_coeffs[2]',
		'beta_extra',
		'beta_tx' if args.transVL == 'prox_pt' else None]

	use_draws = dr[[i for i in pairs_params if i in dr.columns]]
	if args.transForm == 've':
		use_draws = use_draws.rename(columns={
			'logit_txr': r'$\gamma$',
			'logit_fnr': r'$\pi$',
			'mu_0': r'$\mu_0$',
			'omega_mu': r'$\Omega_\mu$',
			'sigma_0': r'$\sigma_0$',
			'omega_sigma': r'$\Omega_\sigma$',
			'est_param1': r'logit($\rho_\text{prod.}$)',
			'est_param2': r'logit($\rho_\text{est.}$)',
			'est_param3': 'log(c)',
			'est_param1_coeffs[1]': r'$\beta^{\rho_\text{prod.}}_\text{FEMALE}$',
			'est_param1_coeffs[2]': r'$\beta^{\rho_\text{prod.}}_\text{MALE}$',
			'est_param2_coeffs[1]': r'$\beta^{\rho_\text{est.}}_\text{FEMALE}$',
			'est_param2_coeffs[2]': r'$\beta^{\rho_\text{est.}}_\text{MALE}$',
			'logit_theta_un_tx': r'logit($\theta_\text{unART}$)',
			'logit_theta_tx': r'logit($\theta_\text{ART}$)',
			'beta_extra': r'$\beta_\text{extra}$',
			'mu_tx': r'$\mu_\text{unART}$',
			'sigma_tx': r'$\sigma_\text{unART}$'})
	use_draws_r2 = use_draws.corr()**2
	p = sns.pairplot(use_draws.sample(1000))
	for row_idx in np.arange(p.axes.shape[0]):
		for col_idx in np.arange(p.axes.shape[1]):
			if row_idx != col_idx:
				_ = p.axes[row_idx, col_idx].text(.05, 0.95, 
					r'$r^2 = $' + f'{use_draws_r2.iloc[row_idx, col_idx]:.2g}', 
					transform=p.axes[row_idx, col_idx].transAxes, fontsize=12,
					 bbox=dict(facecolor='white', edgecolor='none', alpha=0.75))
	# add correlation coefficient text
	p.savefig(f'figures/{out_dir.split('/')[1]}_pairs.pdf')
	p.savefig(f'figures/{out_dir.split('/')[1]}_pairs.png')
	plt.close()

	# transmission plot
	if args.transForm == 'sigmoid':
		for test_vl in [np.log10(200), 3]:
			print(f'Estimated transmission rate with a viral load of {round(test_vl, 2)} log10 copies/mL:',
				file=sys.stderr)
			print(pd.DataFrame(np.quantile(
				dr.use_param1.values / (1 + 10**(-dr.use_param2.values*(np.array([test_vl])[:,np.newaxis] - dr.use_param3.values))),
				[0.025, 0.25, 0.5, 0.75, 0.975],
				axis=1).T,
				columns=[0.025, 0.25, 0.5, 0.75, 0.975]))
		v = np.linspace(0,7,100)
		r = pd.DataFrame(np.quantile(
			dr.use_param1.values / (1 + 10**(-dr.use_param2.values*(v[:,np.newaxis] - dr.use_param3.values))),
			[0.025, 0.25, 0.5, 0.75, 0.975],
			axis=1).T,
			columns=[0.025, 0.25, 0.5, 0.75, 0.975]).\
			assign(v = v)
	elif args.transForm == 've':
		# if we have predictors, use the mean value
		if 'pred' in args.stan:
			dr['use_param1'] = dr.mean_use_param1
			dr['use_param2'] = dr.mean_use_param2
			dr['use_param3'] = dr.mean_use_param3
		for test_vl in [np.log10(200), 3]:
			print(f'Estimated transmission rate with a viral load of {round(test_vl, 2)} log10 copies/mL:',
				file=sys.stderr)
			print(pd.DataFrame(np.quantile(
				108 * dr.use_param1.values * (1 - (1 - dr.use_param2.values)**(dr.use_param3.values * 10**np.array([test_vl])[:,np.newaxis])),
				[0.025, 0.25, 0.5, 0.75, 0.975],
				axis=1).T,
				columns=[0.025, 0.25, 0.5, 0.75, 0.975]),
				file=sys.stderr)
		v = np.linspace(0,7,100)
		r = pd.DataFrame(np.quantile(
			108 * dr.use_param1.values * (1 - (1 - dr.use_param2.values)**(dr.use_param3.values*10**v[:,np.newaxis])),
			[0.025, 0.25, 0.5, 0.75, 0.975],
			axis=1).T,
			columns=[0.025, 0.25, 0.5, 0.75, 0.975]).\
			assign(v = v)
	elif trans_form == 'gamma':
		for test_vl in [np.log10(200), 3]:
			print(f'Estimated transmission rate with a viral load of {round(test_vl, 2)} log10 copies/mL:',
				file=sys.stderr)
			print(pd.DataFrame(np.quantile(
				108*(1 - ((1/dr.use_param1.values) / ((1/dr.use_param1.values) + 10**np.array([test_vl])[:,np.newaxis]))**(dr.use_param2.values)),
				[0.025, 0.25, 0.5, 0.75, 0.975],
				axis=1).T,
				columns=[0.025, 0.25, 0.5, 0.75, 0.975]))
		v = np.linspace(0,7,100)
		r = pd.DataFrame(np.quantile(
			108*(1 - ((1/dr.use_param1.values) / ((1/dr.use_param1.values) + 10**v[:,np.newaxis]))**(dr.use_param2.values)),
			[0.025, 0.25, 0.5, 0.75, 0.975],
			axis=1).T,
			columns=[0.025, 0.25, 0.5, 0.75, 0.975]).\
			assign(v = v)
	fig, axs = plt.subplots(figsize=(4.8, 3.6), constrained_layout=True)
	axs.fill_between(r.v, r[0.025], r[0.975], color='#333333', alpha=0.25, linewidth=0, zorder=1)
	axs.fill_between(r.v, r[0.25], r[0.75], color='#333333', alpha=0.5, linewidth=0, zorder=2)
	axs.plot(r.v, r[0.5], color='#333333', zorder=3)
	axs.set_ylim(-0.025, 0.25)
	axs.set_xlabel(r'viral load (log$_{10}$ copies/mL)')
	axs.set_ylabel(r'transmission rate (year$^{-1}$)')
	fig.savefig(f'figures/{out_dir.split('/')[1]}_vl_trans_rate.pdf')
	fig.savefig(f'figures/{out_dir.split('/')[1]}_vl_trans_rate.png')
	plt.close()
	
	if 'pred' in args.stan:
		v = np.linspace(0,7,100)
		r1 = pd.DataFrame(np.quantile(
			108 * dr.sex1_use_param1.values * (1 - (1 - dr.sex1_use_param2.values)**(dr.sex1_use_param3.values*10**v[:,np.newaxis])),
			[0.025, 0.25, 0.5, 0.75, 0.975],
			axis=1).T,
			columns=[0.025, 0.25, 0.5, 0.75, 0.975]).\
			assign(v = v)
		r2 = pd.DataFrame(np.quantile(
			108 * dr.sex2_use_param1.values * (1 - (1 - dr.sex2_use_param2.values)**(dr.sex2_use_param3.values*10**v[:,np.newaxis])),
			[0.025, 0.25, 0.5, 0.75, 0.975],
			axis=1).T,
			columns=[0.025, 0.25, 0.5, 0.75, 0.975]).\
			assign(v = v)
		fig, axs = plt.subplots(figsize=(4.8, 3.6), constrained_layout=True)
		axs.fill_between(r1.v, r1[0.025], r1[0.975], color='indianred', alpha=0.25, linewidth=0, zorder=1)
		axs.fill_between(r2.v, r2[0.025], r2[0.975], color='steelblue', alpha=0.25, linewidth=0, zorder=1)
		axs.fill_between(r1.v, r1[0.25], r1[0.75], color='indianred', alpha=0.5, linewidth=0, zorder=2)
		axs.fill_between(r2.v, r2[0.25], r2[0.75], color='steelblue', alpha=0.5, linewidth=0, zorder=2)
		axs.plot(r1.v, r1[0.5], color='indianred', label='category1', zorder=3)
		axs.plot(r2.v, r2[0.5], color='steelblue', label='category2', zorder=3)
		axs.legend()
		axs.set_ylim(-0.025, 0.25)
		axs.set_xlabel(r'viral load (log$_{10}$ copies/mL)')
		axs.set_ylabel(r'transmission rate (year$^{-1}$)')
		fig.savefig(f'figures/{out_dir.split('/')[1]}_vl_trans_rate_stratified.pdf')
		fig.savefig(f'figures/{out_dir.split('/')[1]}_vl_trans_stratified.png')
		plt.close()



if __name__ == "__main__":
    run()

