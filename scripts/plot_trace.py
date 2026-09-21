import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import seaborn as sns
try:
	from scripts.utils import plot_style
except:
	from utils import plot_style



def run():
	parser = argparse.ArgumentParser()
	parser.add_argument('--log')
	parser.add_argument('--config')
	parser.add_argument('--label', default='trace')
	args = parser.parse_args()
	#args.log = 'base_params_simulated_fit_mcmc.tsv'
	#args.log = 'base_params_simulated_fit_omega_mcmc.tsv'
	#args.config = 'scripts/inference_params/config_sim.csv'
	plot_style()

	if args.config:
		# flag: clean this up
		config = pd.read_csv(args.config)
		print(config)
		print(config.query('param == "theta_labels"'))
		theta_labels = eval(config.query('param == "theta_labels"')['definition'].values[0])
		theta_vals = eval(config.query('param == "theta"')['definition'].values[0])
		true_vals = {}
		for idx, i in enumerate(theta_vals):
			if isinstance(i, np.ndarray):
				for jdx, j in enumerate(i):
					if jdx == 0:
						true_vals[theta_labels[idx] + f'_{jdx}'] = j
					# FLAG: NOT ROBUST
					elif jdx == 1:
						true_vals[theta_labels[idx] + f'_{jdx}'] = j / (1 - i[0])
			else:
				true_vals[theta_labels[idx]] = i
				if 'omega_scale' in theta_labels[idx]:
					true_vals[f'log10_{theta_labels[idx]}'] = np.log10(i)
		

	log = pd.read_csv(args.log, sep='\t')
	# chck for multi dimensional columns and split those that are
	multi = [(type(i) == type('a')) and ('[' in i) for i in log.iloc[0,:].values]
	log_values = log.iloc[0,:].values.copy()
	for idx, i in enumerate(log_values):
		# flag: clean this up
		if (type(i) == type('a')) and ('[' in i):
			i_split = i.replace('[', '').replace(']', '').split()
			log[['_'.join(i) for i in zip([log.columns[idx]]*len(i_split), np.arange(len(i_split)).astype(str))]] = \
				pd.DataFrame([i.split() for i in log.iloc[:,idx].str.replace('[', '').str.replace(']', '')]).astype(float)


	log = log.drop(log.columns[:len(multi)][multi], axis=1)

	labels_dict = {}
	labels_dict['omega1'] = r'$\omega^\text{init.}$'
	labels_dict['omega2'] = r'$\omega^\text{cease}$'
	

	os.makedirs("figures/simulations", exist_ok=True)
	fig, axs = plt.subplots(log.shape[1]-1, 1, figsize=(8.4, 3.4*(log.shape[1]-1)), constrained_layout=True)
	for idx, i in enumerate(log.columns[1:]):
		_ = axs[idx].plot(log['iter'], log[i].astype(float), color='#333333')
		_ = axs[idx].set_xlabel('\n')
		_ = axs[idx].set_ylabel(labels_dict[i] if i in labels_dict.keys() else i)
		

	axs[-1].set_xlabel('iteration')
	fig.savefig(f'figures/pdf/{args.label}_trace.pdf')
	fig.savefig(f'figures/eps/{args.label}_trace.eps')
	
	plt.close()

	### PAIRS PLOT ##
	trim_log = log.iloc[int(log.shape[0]*0.50):,:]
	print('median parameter values')
	print(trim_log.iloc[:,4:].median())
	print('maximum posterior parameter values')
	print(trim_log.sort_values(by='posterior').iloc[-1,4:])
	trim_log = trim_log.iloc[:,4:]
	rng = np.random.default_rng(seed = 111)
	trim_log = trim_log.iloc[rng.choice(trim_log.shape[0], replace=False, size=np.minimum(1000, trim_log.shape[0])),:]
	
	fig = sns.pairplot(trim_log)
	if args.config:
		for idx in range(fig.axes.shape[0]):
			fig.axes[idx,idx].axvline(true_vals[fig.axes[-1,idx].get_xlabel()], ls='--', color='#333333')
	
	fig.savefig(f'figures/{args.label}_pairs.pdf')
	plt.close()

	trim_trans_log = trim_log
	if 'omega_scale' in trim_trans_log.columns:
		trim_trans_log.\
		assign(log10_omega_scale = lambda k: np.log10(k.omega_scale)).\
		drop(['omega_scale'], axis=1)


	from scipy.stats import pearsonr
	print(trim_trans_log.shape)
	fig = sns.pairplot(trim_trans_log)
	for idx in range(fig.axes.shape[0]):
		if args.config:
			fig.axes[idx,idx].axvline(true_vals[fig.axes[-1,idx].get_xlabel()], ls='--', color='#333333')
		for jdx in range(idx+1, fig.axes.shape[0]):
			i_j_dat = trim_trans_log.loc[:,[fig.axes[idx,jdx].get_xlabel(), fig.axes[idx,jdx].get_ylabel()]]
			r = np.round(pearsonr(i_j_dat.iloc[:,0], i_j_dat.iloc[:,1]).statistic, 2)
			fig.axes[idx,jdx].text(0.05, 0.95, f'r = {r}', transform=fig.axes[idx,jdx].transAxes,
				bbox=dict(facecolor='white', alpha=0.5))
			fig.axes[jdx,idx].text(0.05, 0.95, f'r = {r}', transform=fig.axes[jdx,idx].transAxes,
				bbox=dict(facecolor='white', alpha=0.5))



	fig.savefig(f'figures/{args.label}_trans_pairs.pdf')
	plt.close()


if __name__ == "__main__":
    run()


