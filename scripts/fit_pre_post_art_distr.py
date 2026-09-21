import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
from scipy.stats import gamma, lognorm
import matplotlib.pyplot as plt
import os
try:
	from scripts.utils import plot_style
except:
	from utils import plot_style


def calc_gamma_se(theta):
	from scipy.stats import gamma
	# params from theta 
	alpha, scale = theta
	# Reported quantiles in Treasure 2016 10.1097/QAI.0000000000000964
	x = np.array([0.02, 0.19, 0.50, 1.00, 1.85, 11.87])
	x_delta = np.hstack([x[0], np.diff(x)])
	q = np.array([0.05, 0.25, 0.50, 0.68, 0.75, 0.95])
	# how much does q change time x-step?
	q_delta = np.hstack([q[0], np.diff(q)])
	# estimated q based on theta
	est_q = gamma.cdf(x, a=alpha, scale=scale)
	est_q_delta = np.hstack([est_q[0], np.diff(est_q)])
	# mse normalized to size of interval
	return(((est_q_delta - q_delta)**2).sum())


def calc_lognorm_se(theta):
	from scipy.stats import lognorm
	# params from theta 
	alpha, scale = theta
	# Reported quantiles in Treasure 2016 10.1097/QAI.0000000000000964
	x = np.array([0.02, 0.19, 0.50, 1.00, 1.85, 11.87])
	x_delta = np.hstack([x[0], np.diff(x)])
	q = np.array([0.05, 0.25, 0.50, 0.68, 0.75, 0.95])
	# how much does q change time x-step?
	q_delta = np.hstack([q[0], np.diff(q)])
	# estimated q based on theta
	est_q = lognorm.cdf(x, s=alpha, scale=scale)
	est_q_delta = np.hstack([est_q[0], np.diff(est_q)])
	# mse normalized to size of interval
	return(((est_q_delta - q_delta)**2).sum())


result = minimize(calc_lognorm_se, [2, 2])

x = np.array([0.02, 0.19, 0.50, 1.00, 1.85, 11.87])
x_range = np.linspace(0, np.ceil(x[-1]), 1000)
cdf_est = lognorm.cdf(x_range, s=result.x[0], scale=result.x[1])
pdf_est = lognorm.pdf(x_range, s=result.x[0], scale=result.x[1])
q = np.array([0.05, 0.25, 0.50, 0.68, 0.75, 0.95])

plot_style()
fig, axs = plt.subplots(1, 2, figsize=(6.4*1.5, 4.8), constrained_layout=True)
axs[0].grid(axis='both', zorder=1, color='#eaeaea')
axs[0].plot(x_range, cdf_est, color='white', linewidth=3, zorder=2)
axs[0].plot(x_range, cdf_est, color='#333333', linewidth=2, zorder=3)
axs[0].plot(x_range, cdf_est, color='#333333', zorder=4, label='estimated')
axs[0].scatter(x, q, facecolor='white', edgecolor='white', linewidth=2, s=50, zorder=5)
axs[0].scatter(x, q, facecolor='#eaeaea', edgecolor='#333333', s=50, zorder=6, 
	label='reported (Treasure et al., JAIDS, 2016)')
axs[0].legend(fontsize=10, frameon=True)
axs[0].set_ylabel('cumulative density')
#axs[0].set_xlim(x_range[0], x_range[-1])
axs[1].grid(axis='both', zorder=1, color='#eaeaea')
axs[1].plot(x_range, pdf_est, color='white', linewidth=3, zorder=2)
axs[1].plot(x_range, pdf_est, color='#333333', linewidth=2, zorder=3)
axs[1].plot(x_range, pdf_est, color='#333333', zorder=4, label='estimated')
axs[1].set_ylabel('density')
axs[0].set_ylim(0,1)
axs[1].set_ylim(0,1.6)
for ax in axs:
	#ax.set_xticks([0, 1, 2, 3, 4, 5, 10, 15])
	ax.set_xticks(np.arange(0,13))
	ax.set_xlim(-12*0.025, 12*1.025)


fig.supxlabel('pre-ART to post-rebound\nset-point fold-change')
#axs[1].set_xlim(x_range[0], x_range[-1])
os.makedirs("figures/pdf", exist_ok=True)
fig.savefig('figures/pdf/pre_post_art_spvl.pdf')
os.makedirs("figures/eps", exist_ok=True)
fig.savefig('figures/eps/pre_post_art_spvl.eps')
plt.close()

print(f'LOGNORM_S:{result.x[0]}\nLOGNORM_SCALE:{result.x[1]}')

