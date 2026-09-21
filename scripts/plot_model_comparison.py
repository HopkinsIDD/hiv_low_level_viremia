import pandas as pd
import numpy as np
import argparse
import cmdstanpy
import matplotlib.pyplot as plt
import glob
import os
try:
	from scripts.utils import import_config, get_ind_n_trans_bd, summarize_draws, plot_style, get_trans_form, format_list_col
	from scripts.plot_vl_trans_risk import plot_trans_risk
except:
	from utils import import_config, get_ind_n_trans_bd, summarize_draws, plot_style, get_trans_form, format_list_col
	from plot_vl_trans_risk import plot_trans_risk


def run():
	parser = argparse.ArgumentParser()
	# input files
	args = parser.add_argument('--llvPrevDat')
	args = parser.add_argument('--mergedPeriodDat')
	args = parser.add_argument('--fitDat', nargs='+')
	args = parser.add_argument('--fitDatLabels', nargs='+')
	args = parser.add_argument('--config', default='config/config.csv')
	args = parser.parse_args()
	'''
	args.mergedPeriodDat = 'output/ve_constrained/periods.tsv'
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
		'Main_analysis',
		'Sex_as_predictor',
		'Estimated_$c$',
		'Artificial_genetic_constraints',
		'No_genetic_constraints',
		'Low_ART-TN_assumption',
		'High_ART-TN_assumption',
		'No_saturation_of_transmission',
		'Pre-treatment_participants_only',
		'Sigmoidal transmission']
	'''
	config = import_config(args.config)

	#fitDr = [cmdstanpy.from_csv(
	#			glob.glob(i + '/*.csv')).draws_pd().assign(trans_form = get_trans_form(i)) for i in args.fitDat]

	# format merged period dat
	#vmax = np.ceil(np.log10(np.nanmax(np.hstack(format_list_col(pd.read_csv(args.mergedPeriodDat, sep='\t').\
	#	query('conversion')['donor_copies'])).astype(float))))
	vmax=7
	v = np.linspace(1, vmax, 500)
		
	vl_ticks = np.log10([1, 200, 1000, 10000, 100000, 10**6, 10**7])
	vl_labels = [f'{int(10**k):,}' if 10**k < 1000000 else f'{int(10**k/1000000)}M'for k in vl_ticks]


	plot_style()
	fig, axs = plt.subplots(int(np.ceil(len(args.fitDat)/5)), 5,
		figsize=(6.4*3, 4.8*(np.ceil(len(args.fitDat)/5))), constrained_layout=True)
	axs_flat = axs.flatten()
	for drdx, drf in enumerate(args.fitDat):
		dr = pd.read_csv(drf+'/select_draws.tsv', sep='\t').assign(trans_form = get_trans_form(drf))
		#cmdstanpy.from_csv(
		#		glob.glob(drf + '/*.csv')).draws_pd().assign(trans_form = get_trans_form(drf))
		# get transmission risk across VL range	
		if 'mean_use_param1' in dr.columns:
			dr['use_param1'] = dr.mean_use_param1
			dr['use_param2'] = dr.mean_use_param2
			dr['use_param3'] = dr.mean_use_param3
		r = get_ind_n_trans_bd(dr, t=100, v=v, form=dr.trans_form.iloc[0])
		# summarize
		r_sum = summarize_draws(r).\
			assign(v = v)
		vl_ticks = np.log10([1, 200, 1000, 10000, 100000, 10**6, 10**7])
		vl_labels = [f'{int(10**k):,}' if 10**k < 1000000 else f'{int(10**k/1000000)}M'for k in vl_ticks]
		#config = {i['var']: i.value for idx, i in pd.read_csv('config/config.tsv', sep='\t').iterrows()}
		axs_flat[drdx] = plot_trans_risk(axs_flat[drdx], r_sum, config, annotation=True if 'nosaturation' not in drf else False)
		#axs_flat[drdx].set_xticks(vl_ticks, labels=vl_labels, size=9)
		axs_flat[drdx].set_ylim(0, 25)
		axs_flat[drdx].set_xlim(1, np.floor(r_sum.v.max()))
		axs_flat[drdx].set_xticks(vl_ticks, labels=vl_labels, size=9)
		axs_flat[drdx].set_title(f'{drdx+1}. ' + args.fitDatLabels[drdx].replace('_', ' '))
		axs_flat[drdx].set_xlabel('viral load (copies/mL)\n')
		axs_flat[drdx].set_ylabel('\ntarnsmissions/\n100 person-years')


	os.makedirs("figures/pdf", exist_ok=True)
	fig.savefig('figures/pdf/intra_couple_trans_risk_all.pdf')
	os.makedirs("figures/eps", exist_ok=True)
	fig.savefig('figures/eps/intra_couple_trans_risk_all.eps')
	plt.close()




if __name__ == '__main__':
    run()
