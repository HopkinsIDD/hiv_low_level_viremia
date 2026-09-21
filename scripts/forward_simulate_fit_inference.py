import os
import sys
import argparse
import pandas as pd
import numpy as np
import scipy.stats
#from line_profiler import profile


#@profile
def bincount_2d(a, minlength):
	return(np.bincount(
			--			(a + (np.arange(a.shape[0])*(minlength))[:,np.newaxis]).flatten(),
					minlength=a.shape[0] * (minlength)).\
				reshape((a.shape[0], (minlength))))


def ffill(arr):
	mask = np.isnan(arr)
	out = np.maximum.accumulate(np.where(mask, 0, arr))
	out[(out == 0) & mask] = np.nan
	return out


#@profile
def remove_neg(a):
	return(a[a < 1000])


#@profile
def stickbreak(i):
	return(np.cumprod(np.hstack([1, 1-i]))*np.hstack([i, 1]))


#v = tstates[0,tcols['v'], adherent]
#s = states[cols['s'],adherent]
#d1 = states[cols['d1'],adherent]
#d2 = states[cols['d2'],adherent]
#a = states[cols['a'],adherent]
#t = tstates[:,tcols['t'],0]

#@profile
def adherent_viral_decay(v, cp, d1, d2, t):
	import numpy as np
	# biphasic exponential viral decay
	# for adherent individuals only
	# they have initiated treatment from ~set point and 
	# are on their way to suppression
	# v = current log10 viral load
	# a = proportion of viral load subject to second-phase decay
	# d1 = intial phase decay rate
	# d2 = second phase decay rate
	# t = time steps 
	# assumes that viral loads started decaying from set point and follow
	# V(t) = \max(v - s A, 0) e^{-d_1 t} + s a e^{-d_2 t}
	# output array
	t_col = t[:,np.newaxis]
	pow10_v = 10**v
	pow10_cp = 10**cp
	#vt = np.log10(np.maximum(pow10_v - pow10_cp, 0) * np.exp(-d1*t_col) + np.minimum(pow10_v, pow10_cp)*np.exp(-d2*t_col))
	vt = np.log10(np.maximum(
		np.maximum(pow10_v - pow10_cp, 0) * np.exp(-d1*t_col) + np.minimum(pow10_v, pow10_cp)*np.exp(-d2*t_col),
		1))
	return(vt)


#@profile
def dt_on_tx(v, cp, d1_dt, d2_dt):
	import numpy as np
	# biphasic exponential viral decay
	# for adherent individuals only
	# they have initiated treatment from ~set point and 
	# are on their way to suppression
	# v = current viral load
	# s = set point viral load
	# d1 = intial phase decay rate
	# d2 = second phase decay rate
	# a = proportion of set point viral load subject to second phase
	# dt = time step 
	pow10_v = 10**v
	pow10_cp = 10**cp
	#print(" ")
	#print(pow10_v)
	#print(pow10_cp)
	#print(np.maximum(pow10_v - pow10_cp, 0))
	#print(np.minimum(pow10_v, pow10_cp))
	vt = np.log10(
		np.maximum(
			np.maximum(pow10_v - pow10_cp, 0) * np.exp(-d1_dt) + np.minimum(pow10_v, pow10_cp)*np.exp(-d2_dt),
			1))
	#vt = np.log10(
	#		np.maximum(pow10_v - pow10_cp, 0) * np.exp(-d1_dt) + np.minimum(pow10_v, pow10_cp)*np.exp(-d2_dt))
	return(vt)




#@profile
def non_adherent_viral_growth(v, rp, pp, zr, zd3, t):
	# flag: need to track viral load peak
	t_col = t[:,np.newaxis]
	# output array
	vt = np.zeros((t.shape[0], v.shape[0]))
	vt[0,:] = v
	# current viral load is greater than rebound point
	# then decay to rebound point
	v_gt_rp = v > rp
	vt[1:,v_gt_rp] = np.maximum(v[v_gt_rp] - zd3[v_gt_rp]*t_col[1:,:], rp[v_gt_rp])
	# current viral load is less then rebound point
	# then grow to acute viral load peak
	v_lt_rp = v < rp
	if any(v_lt_rp):
		# exact time to viral load peak
		#t_to_peak = ((np.log((10**peak_log10_v) / v[v_lt_rp]) / r[v_lt_rp]) / (t[1] - t[0]))
		t_to_peak = (pp[v_lt_rp] - v[v_lt_rp]) / zr[v_lt_rp]
		# time step just before reaching viral load peak
		dt_before_peak = t_to_peak.astype(int)
		max_dt_before_peak = dt_before_peak.max()+1
		vt[1:max_dt_before_peak,v_lt_rp] = \
			np.minimum(v[v_lt_rp] + zr[v_lt_rp]*t_col[1:max_dt_before_peak,:], pp[v_lt_rp])
		vt[max_dt_before_peak:, v_lt_rp] = pp[v_lt_rp]
		# exact time from viral load peak to viral load set point
		# delta between p and s
		delta = pp - rp
		t_peak_to_rp = (delta[v_lt_rp]) / zr[v_lt_rp]
		t_to_rp = t_to_peak + t_peak_to_rp
		# time steps to just before set point
		dt_to_rp = t_to_rp.astype(int)
		decay_from_peak = \
			np.clip(-zd3[v_lt_rp] * (t_col - t_to_peak),
				a_min=-delta[v_lt_rp], a_max=0)
		vt[:,v_lt_rp] = vt[:,v_lt_rp]+ decay_from_peak
	# finally, if equal to set point, stay there
	vt[:,v == rp] = rp[v==rp]
	return(vt)


#s = states[cols['s'], rebounding]
#p = states[cols['p'], rebounding]
#zr_dt = states[cols['zr_dt'], rebounding]
#zd3_dt = states[cols['zd3_dt'], rebounding]
#r_rel_d3 = states[cols['r_rel_d3'], rebounding]
#pre_w = tstates[t-1,tcols['w'], rebounding]
#pre_rho = tstates[t-1,tcols['rho'], rebounding]
#pre_v = tstates[t-1,tcols['v'], rebounding]
#pre_pre_v = tstates[t-2,tcols['v'], rebounding]
#dt = config['dt']
#rand_draws = new_off_tx_rnd[t,:]

#@profile
def dt_off_tx(rp,pp,zr_dt,zd3_dt,r_rel_d3,pre_w,pre_rho,pre_v,pre_pre_v,dt,rand_draws):
	import numpy as np
	pre_w_bool = pre_w.astype(bool)
	pre_rho_bool = pre_rho.astype(bool)
	vt = np.zeros((pre_v.shape[0],2))
	#### CASES DURING WHICH VIRAL POPULATION GROWS ####
	# viral load less then set point
	pre_v_lt_rp = pre_v < rp
	# first time point off treatment and viral load greater than or equal to set-point
	# flag: any way to do this without creating new array of 0s? 
	new_off_tx_pre_v_gt_rp = pre_w_bool & ~pre_v_lt_rp
	new_off_tx_pre_v_gt_rp_grow = np.zeros((pre_v.shape[0])).astype(bool)
	# generate a random number as we randomly assign growth or decay for new_off_tx_v_gt_s
	# p(grow) = d3 / (r + d3)
	# because faster growth relative to decay decreases time in growth phase
	new_off_tx_pre_v_gt_rp_grow[new_off_tx_pre_v_gt_rp] = \
		rand_draws[:new_off_tx_pre_v_gt_rp.sum()] > r_rel_d3[new_off_tx_pre_v_gt_rp]
	# finally, if vl grew at last time step and hasn't yet peaked
	grow_not_yet_peaked = ~pre_w_bool & (pre_v > pre_pre_v) & ~pre_rho_bool
	to_grow = pre_v_lt_rp | new_off_tx_pre_v_gt_rp_grow | grow_not_yet_peaked
	grew = pre_v[to_grow] + zr_dt[to_grow]
	# check if peaked in the middle of the time frame
	# need to ensure does not drop below 
	peaked = grew > pp[to_grow]
	if any(peaked):
		vt[to_grow,0] = peaked
		grew[peaked] = np.maximum(
			pp[to_grow][peaked] - \
				zd3_dt[to_grow][peaked] * \
					(1 - (pp[to_grow][peaked] - pre_v[to_grow][peaked]) / zr_dt[to_grow][peaked]),
			rp[to_grow][peaked])
	vt[to_grow,1] = grew
	#### CASES WHERE VIRAL POPULATION DECAYS ###
	to_decay = (new_off_tx_pre_v_gt_rp & ~new_off_tx_pre_v_gt_rp_grow) | ((pre_v < pre_pre_v) & (pre_v >= rp)) | pre_rho_bool
	vt[to_decay,1] = np.maximum(pre_v[to_decay] - zd3_dt[to_decay], rp[to_decay])
	#### CASES WHERE VIRAL POPULATION IS STABLE ####
	at_rp = (pre_v == pre_pre_v) & (pre_v == rp)
	vt[at_rp,1] = rp[at_rp]
	return(vt)


#@profile
def treatment_changes(x, ind_tau, params, config):
	x = x.astype(int)
	max_t_index = config['t_max'] - 1
	# pulling out the first column in x as an integer
	#x0 = x[:,0,:].astype(int)
	# first, get random values
	change_rand_draws = config['rng'].uniform(0,1,(x.shape[0]-1, x.shape[2]))
	# calculate transition (y/n) conditional on current state
	omega = params['omega']
	scaled_rate = omega * config['dt']
	change = (change_rand_draws <= scaled_rate[:, np.newaxis, np.newaxis]).astype(int)
	# set treatment status, first at t0
	on_tx = x[0,0,:].astype(bool)
	# then iterate through remaining time steps
	ind_idx = np.arange(x.shape[2])
	for t in np.arange(1,x.shape[0]):
		x[t,0,:] = np.abs(x[t-1,0,:] - change[x[t-1,0,:],t-1,ind_idx])
	# which changed treatment?
	n_tx_dt = x[:,0,:].sum(axis=0)
	delta_tx = (n_tx_dt > 0) & (n_tx_dt < x.shape[0])
	# set wash out status
	# by default not in wash-out
	x[:,1,:] = 0 
	# tabulate where on treatment
	delta_tx_on_tx = x[:,0,:].astype(bool)
	# when last on treatment
	delta_tx_last_tx = np.where(delta_tx_on_tx, np.arange(delta_tx_on_tx.shape[0])[:,np.newaxis], np.nan)
	# then those not on treatment at time = 0 
	# who may have been on treatment at at time point just before time = 0
	# if not on tx at t=0, then approximate p(wash_out) as:
	# tau / (1/omega_on)
	# flag: confirm correct
	delta_tx_off_tx_0 = delta_tx & ~x[0,0,:].astype(bool)
	p_washout = ind_tau[delta_tx_off_tx_0] * params['omega'][0] * config['dt']
	delta_tx_last_tx[0,delta_tx_off_tx_0] = \
		np.where(p_washout > config['rng'].uniform(size=p_washout.shape[0]), -1, delta_tx_last_tx[0,delta_tx_off_tx_0])
	# forward fill and then replace nans with -inf for those never on treatment
	delta_tx_last_tx[:,delta_tx] = np.nan_to_num(ffill(delta_tx_last_tx[:,delta_tx]), nan=-np.inf)
	# finally, determine whether delta_tx people are in washout
	# get washout duration 
	delta_tx_ind_dt_tau = np.round(ind_tau[delta_tx]/config['dt']).astype(int)
	x[:,1,delta_tx] = (np.arange(x.shape[0])[:,np.newaxis] - delta_tx_last_tx[:,delta_tx]) <= delta_tx_ind_dt_tau
	return(x, delta_tx)


#@profile
def initiate_states(params, params_generator, config):
	# flag: below comment is out of date, to update
	# tstates (time inhomogeneous) = 
	#	0: time (t),
	#	1: on treatment (x_i,t),
	#	2: in wash out period (w_i,t),
	#	3. reached viral load peak (rho_i,t) 
	#	4: viral load (v_i,t)
	# states (time homogeneous) =
	# flag: needs to be updated
	#	0: category (g_i), 0: non-adherers, 1: transient-adherers, 2: adherers
	#	1: pre-treatment set-point (sp_i)
	#	2: post-rebound set-point (rb_i)
	#	2: first phase viral decay rate (d1_i),
	#	3: second phase viral decay rate (d2_i),
	#	4. fraction of viral load subject to second phase decay (a_i),
	#	5: viral growth rate (r_i),
	#	6: aucte peak above set point (delta__i)
	#	7. viral decay rate from acute peak (d3_i)
	#	8. duration of washout period (tau_i)
	#### TIME HOMOGENEOUS STATES ####
	# flag: clean this up, make it more flexible
	state_cols = ['g', 'sp', 'delta_rp', 'd1', 'd2', 'a', 'r', 'delta_pp', 'd3', 'tau']
	len_state_cols = len(state_cols)
	transformed_rate_cols = ['d1', 'd2', 'r', 'd3']
	relative_rate_cols = [('r', 'd3')]
	# add thing for composite params
	# pp = peak-point
	# rp = rebound-point
	# cp = change-point from fast to slow-phase decay on the log10 copies/mL scale
	composite_params = ['rp', 'pp', 'cp']
	rp_idx, pp_idx, cp_idx = np.arange(len(composite_params))
	composite_funcs = [
		lambda states, cols: states[cols['sp']] + np.log10(states[cols['delta_rp']]),
		lambda states, cols: states[cols['sp']] + np.log10(1 + states[cols['delta_pp']]),
		lambda states, cols: states[cols['sp']] + np.log10(states[cols['a']])]
	# define states
	states = np.full((len_state_cols + len(transformed_rate_cols)*3 + len(relative_rate_cols) + len(composite_params) + 1,
			config['N']), 
		np.nan)
	cols = {col:cdx for cdx, col in enumerate(state_cols)}
	states[:len(state_cols)] = np.vstack([params[i] for i in state_cols])
	# add transform of some rate parameters to be on the log10 viral load scale
	to_transform_rate_indices = np.array([cols[i] for i in transformed_rate_cols])
	range_transformed_rate_cols = np.arange(len(transformed_rate_cols))
	for idx, i in enumerate(transformed_rate_cols):
		cols.update({
		i+'_dt': len_state_cols + idx,
		'z'+i: len_state_cols + 4 + idx,
		'z'+i+'_dt': len_state_cols + 8 + idx})
	states[len_state_cols + range_transformed_rate_cols] = states[to_transform_rate_indices] * config['dt']
	states[len_state_cols + range_transformed_rate_cols + 4] = states[to_transform_rate_indices] / np.log(10)
	states[len_state_cols + range_transformed_rate_cols + 8] = states[len_state_cols + range_transformed_rate_cols + 4] * config['dt']
	# add relative rate parameters
	for cdx, col in enumerate(relative_rate_cols):
		cols[col[0] + '_rel_' + col[1]] = len_state_cols + 12 + cdx
		states[cols[col[0] + '_rel_' + col[1]]] = states[cols[col[0]]] / (states[cols[col[0]]] + states[cols[col[1]]])
	# add composite parameters
	for idx, i in enumerate(composite_params):
		cols[i] = states.shape[0] - len(composite_params) + idx - 1
		states[cols[i]] = composite_funcs[idx](states, cols)
	# if peak viral load is < current viral load,
	# redraw set-point, delta_rp, and delta_pp
	to_redraw = params['v0'] > states[cols['pp']]
	n_to_redraw = to_redraw.sum()
	while n_to_redraw > 0:
		new_config = {
			'N': n_to_redraw,
			'random_draws': config['rng'].uniform(size=n_to_redraw*2)}
		states[cols['sp'], to_redraw] = params_generator['sp']('-', params, new_config)
		states[cols['delta_pp'], to_redraw] = params_generator['delta_pp']('-', params, new_config)
		states[cols['delta_rp'], to_redraw] = params_generator['delta_rp']('-', params, new_config)
		states[cols['pp'], to_redraw] = \
			composite_funcs[pp_idx](states[:,to_redraw], cols)
		states[cols['rp'], to_redraw] = \
			composite_funcs[rp_idx](states[:,to_redraw], cols)
		to_redraw = params['v0'] > states[cols['pp']]
		n_to_redraw = to_redraw.sum()
	# and actually finally, add indicator for whether treatment status changs
	cols['delta_tx'] = states.shape[0] - 1
	states[cols['delta_tx']] = states[cols['g'],:] == 1
	delta_tx = states[cols['delta_tx']].astype(bool)
	##### TIME INHOMOGENEOUS STATES #####
	tstates = np.full(
		(int(config['t_max']/config['dt']),
			5,
			config['N']), 
		0).astype(float)
	tcols = {'t': 0, 'x': 1, 'w': 2, 'rho': 3, 'v': 4}
	# time column
	# flag don't repeat this each time
	tstates[:,tcols['t'],:] = np.arange(0,config['t_max'], config['dt'])[:,np.newaxis]
	# initial viral load
	tstates[0,tcols['v'],:] = params['v0']
	# t0 treatment status
	# x = 0 = off treatment
	# x = 1 = on treatment
	# if non-adherer then alwayss 0
	# if adherer then always one
	# for both, same for all time points
	tstates[:,tcols['x'],~delta_tx] = states[cols['g'],~delta_tx].astype(bool)
	# for transient adherers 
	tstates[0,tcols['x'],delta_tx] = params['x0'][delta_tx]
	return((states,cols),(tstates,tcols))


#@profile
def run_simulation(params, params_generator, config):
	#### TIME HOMOGENEOUS STATES ####
	(states,cols), (tstates,tcols) = initiate_states(params, params_generator, config)
	delta_tx = states[cols['delta_tx']].astype(bool)
	# for those who do change, need to simulate gillespie transitions
	# probability of change is state dependent
	tstates[:,tcols['x']:(tcols['w']+1),delta_tx], delta_tx[delta_tx] = treatment_changes(
			tstates[:,:,delta_tx][:,[tcols['x'], tcols['w']],:], 
			states[cols['tau'],delta_tx],
			params,
			config)
	# draw random numbers for treatment changes in case > set point
	# draws an excess of ranodm numbers but does it all at once 
	# so potentially more efficient
	# flag: to confirm
	new_off_tx_rnd = config['rng'].uniform(size=(tstates.shape[0], delta_tx.sum()))
	# whether in washout period or not for those who do not 
	# change treatment is the same as treatment status
	tstates[:,tcols['w'],~delta_tx] = tstates[:,tcols['x'],~delta_tx]
	# saving this for later
	on_tx = tstates[:,tcols['x'],:].astype(bool)
	delta_tx_on_tx = (delta_tx & on_tx).astype(bool)
	# simulate viral load dynamics conditional on
	# treatment status and growth and decay rates 
	# for those with no change, can do all at once
	# start with not on treatment
	# non_adherrent_viral_growth(v, s, p, zr, zd3, t)
	non_adherent = (~delta_tx & ~on_tx[0,:])
	if any(non_adherent):
		tstates[:,tcols['v'],non_adherent] = \
			non_adherent_viral_growth(
				tstates[0,tcols['v'], non_adherent],
				states[cols['rp'],non_adherent],
				states[cols['pp'],non_adherent],
				states[cols['zr'],non_adherent],
				states[cols['zd3'],non_adherent],
				tstates[:,tcols['t'],0])
	# then go on treatment
	# adherent_viral_decay(v, cp, d1, d2, a, t)
	adherent = (~delta_tx & on_tx[0,:])
	if any(adherent):
		tstates[:,tcols['v'],adherent] = \
			adherent_viral_decay(
				tstates[0,tcols['v'], adherent],
				states[cols['cp'],adherent],
				states[cols['d1'],adherent],
				states[cols['d2'],adherent],
				tstates[:,tcols['t'],0])
	if any(delta_tx):
		#washout_rand_draws = rng.uniform(size=p_washout.shape[0])
		#off_tx_to_set_w1 = remove_neg(ind_dt_tau_arr[:,~on_tx,:][:,washout_rand_draws < p_washout,:])
		for t in np.arange(1,config['t_max'], config['dt']):
			#if (t > 100) & (tstates[t-1, tcols['v'], 0] < tstates[t-2, tcols['v'], 0]):
			#	break
			# 1. Where on treatment, so declining
			tstates[t,tcols['v'],delta_tx & on_tx[t,:]] = \
				dt_on_tx(
					tstates[t-1,tcols['v'], delta_tx_on_tx[t,:]],
					states[cols['cp'],delta_tx_on_tx[t,:]],
					states[cols['d1_dt'],delta_tx_on_tx[t,:]],
					states[cols['d2_dt'],delta_tx_on_tx[t,:]])
			# 2. When not on treatment
			# 2.a. When not on treatment and in wash out, so stable
			wash_out = delta_tx & \
						~on_tx[t,:] & \
						tstates[t,tcols['w'],:].astype(bool)			
			tstates[t, tcols['v'], wash_out] = \
				tstates[t-1, tcols['v'], wash_out]
			# 2.b. When not on treatment and not in wash out, so rebounding
			rebounding = delta_tx & ~on_tx[t,:] & ~tstates[t,tcols['w'],:].astype(bool)
			# flag: assumes rho and v are consecutive columns, fix
			# do we need to reseed?
			#np.maximum(
					#	tstates[t-1,tcols['v'], rebounding],
					#	0),
			tstates[t, tcols['rho']:(tcols['v']+1), rebounding] =\
				dt_off_tx(
					states[cols['rp'], rebounding],
					states[cols['pp'], rebounding],
					states[cols['zr_dt'], rebounding],
					states[cols['zd3_dt'], rebounding],
					states[cols['r_rel_d3'], rebounding],
					tstates[t-1,tcols['w'], rebounding],
					tstates[t-1,tcols['rho'], rebounding],
					tstates[t-1,tcols['v'], rebounding],
					tstates[t-2,tcols['v'], rebounding] if t>1 else np.repeat(np.nan, (rebounding.sum())),
					config['dt'],
					new_off_tx_rnd[t,:])
	return(tstates, tcols, states, cols)


#@profile
def run_iteration(full_theta,data,params_generator, config):
	rng = config['rng']
	# set random draws to use in this iteration
	config['random_draws'] = config['rng'].uniform(size=config['n_random_draws'])
	params = {config['theta_labels'][idx]: i for idx, i in enumerate(full_theta)}
	for k, v in params_generator.items():
		params[k] = v(data, params, config)
	# flag match output to previous formats
	tstates, tcols, states, cols = run_simulation(params, params_generator, config)
	# get bin counts for each sampled time point
	# viral loads of 0 assigned to first bin
	simulated_bin_counts = \
		bincount_2d(
			np.digitize(tstates[data.unique_dt[0],tcols['v'],:], data.bin_edges[0], right=True),
			data.bin_edges[0].shape[0]+1)
	# add 0.5 to account for 0s 
	# this means that log-likelihood is not directly comparable when config['N'] changes
	simulated_bin_probs = (simulated_bin_counts + 0.5) / \
		(tstates.shape[2] + simulated_bin_counts.shape[1]*0.5)
	# data llh
	llh = np.log(simulated_bin_probs[data['unique_dt_idx'], data['v1_bin']]).sum()
	return(llh, simulated_bin_probs[-1,:], config)


#@profile
def mcmc_update(current_lposterior, current_lprior, current_llh, theta, data,params_generator, config):
	rng = config['rng']
	# first, randomly sample whcih parameters to update, each has 50% chance
	# can help break out of weird correlated parts of parameter space
	if len(theta) > 1:
		to_sample = config['rng'].uniform(size=len(theta)) > 0.50
	else:
		to_sample=np.array([True])
	proposed_theta = [i(theta[idx], config) if to_sample[idx] else theta[idx] for idx, i in enumerate(config['proposals'])]
	# expand simplex parameeters
	proposed_full_theta = [stickbreak(i) if 
			isinstance(i, np.ndarray) else i for i in proposed_theta]
	proposed_llh, _, config = run_iteration(proposed_full_theta, data, params_generator, config)
	proposed_lprior = sum([i.logpdf(proposed_full_theta[idx]) for idx, i in enumerate(config['priors'])])
	proposed_lposterior = proposed_lprior + proposed_llh
	if np.exp(proposed_lposterior - current_lposterior) > config['accept_threshold']:
		return(proposed_lposterior, proposed_lprior, proposed_llh, proposed_theta, data,params_generator, config)
	else:
		#if proposed_lposterior > current_lposterior:
		#	print("error")
		#	print(current_lposterior)
		#	print(proposed_lposterior)
		#	print(theta)
		#	print(proposed_theta)
		return(current_lposterior, current_lprior, current_llh, theta,data,params_generator, config)



#@profile
def mcmc(data,params_generator, config):
	rng = config['rng']
	# draw initial parameters
	# flag: this is not right for multi dimensional parameters!
	#theta = [i.rvs() for idx, i in enumerate(config['priors'])]
	theta = [i[:-1] if isinstance(i, np.ndarray) else i for i in config['theta']]
	print(f'running inference for {config["n_iter"]} MCMC iterations', file=sys.stderr)
	with open(config['out_file'] + '_mcmc.tsv', 'w') as fp:
		_ = fp.write('iter\tposterior\tprior\tlikelihood\t' + '\t'.join(config['theta_labels']) + '\n')
		i = 0
		full_theta = [stickbreak(i) if 
				isinstance(i, np.ndarray) else i for i in theta]
		config['accept_threshold'] = config['rng'].uniform()
		llh, _, config = run_iteration(full_theta, data, params_generator, config)
		lprior = sum([i.logpdf(full_theta[idx]) for idx, i in enumerate(config['priors'])])
		lposterior = lprior + llh
		print(f'{i}\t{lposterior}\t{lprior}\t{llh}', file=sys.stderr)
		# flag: don't waste time doing the string replacement?
		_ = fp.write(f'{0}\t{lposterior}\t{lprior}\t{llh}\t{'\t'.join([str(i) for i in theta])}\n')
		for i in np.arange(1, config['n_iter']):
			config['accept_threshold'] = config['rng'].uniform()
			lposterior, lprior, llh,theta, _, _, config = \
				mcmc_update(lposterior, lprior, llh, theta, data, params_generator, config)
			if i % 1000 == 0:
				print(f'{i}\t{lposterior}\t{lprior}\t{llh}', file=sys.stderr)
			_ = fp.write(f'{i}\t{lposterior}\t{lprior}\t{llh}\t{'\t'.join([str(i) for i in theta])}\n')


def plot_vl_trajectories(in_states):
	import matplotlib.pyplot as plt
	import matplotlib as mpl
	states = in_states.copy()
	states[:,1,:] = np.where(states[:,1,:] < 0, 0, states[:,1,:])
	max_y = states[:,1,:].max()*1.1
	fig = plt.figure(figsize=(12.8, 6.4))
	gs = fig.add_gridspec(1, 4)
	ax1 = fig.add_subplot(gs[0, 0])
	ax2 = fig.add_subplot(gs[0, 1:3])
	ax3 = fig.add_subplot(gs[0,3])
	_ = ax1.hist(states[0,1,:], 
		bins=np.arange(0, max_y, 0.1),
		facecolor='#eaeaea', edgecolor='#333333',
		orientation='horizontal')
	ax1.invert_xaxis()
	ax1.set_ylabel('viral load (log$_{10}$ copies/mL)', size=18)
	ax1.set_xticks([])
	ax1.set_ylim(0, max_y)
	[ax1.spines[i].set_visible(False) for i in ['bottom',  'top', 'right']]
	for i in np.arange(states.shape[2]):
		_ = ax2.plot(states[:,0,i], states[:,1,i],
			color='#333333', alpha=0.25)
	ax2.axhline(3, ls='--', color='steelblue')
	ax2.axhline(np.log10(200), ls='--', color='steelblue')
	ax2.axhline(np.log10(40), ls='--', color='indianred')
	ax2.set_ylim(0, max_y)
	ax2.set_xlim(0,states[-1,0,0])
	ax2.set_yticks([])
	#ax2.set_xticks(np.arange(0,states[-1,0,0]+0.25,0.25))
	ax2.set_xlabel('time (days)')
	ax3.hist(states[-1,1,:], 
		bins=np.arange(0, max_y, 0.1),
		facecolor='#eaeaea', edgecolor='#333333',
		orientation='horizontal')
	ax3.set_xticks([])
	ax3.set_yticks([])
	ax3.set_ylim(0, max_y)
	[ax3.spines[i].set_visible(False) for i in ['bottom', 'left', 'top', 'right']]
	fig.tight_layout()
	fig.subplots_adjust(hspace=0, wspace=0)
	plt.close()
	return(fig)


def simulate(params_generator, config):
	rng = config['rng']
	# set random draws to use in this iteration
	config['random_draws'] = config['rng'].uniform(size=config['n_random_draws'])	
	#config['theta'][1] = 0.08
	#config['theta'][2] = 0.02
	params = {config['theta_labels'][idx]: i for idx, i in enumerate(config['theta'])}
	# read in data or make mock data
	if 'data' in config.keys():
		data = pd.read_csv(config['data'], sep='\t')
	else:
		data = pd.DataFrame({'v0': config['rng'].uniform(np.log10(200),np.log10(1000),size=config['N'])})
	for k, v in params_generator.items():
		params[k] = v(data, params, config)
	# draw time between samples
	sample_t = rng.normal(config['t_mean'], config['t_scale'], size=config['N'])
	config['t_max'] = np.ceil(sample_t.max()).astype(int)
	tstates, tcols, states, cols = run_simulation(params, params_generator, config)
	tdat = pd.DataFrame(
			np.vstack([tstates[:,:,i] for i in np.arange(0,tstates.shape[2])]),
			columns=list(tcols.keys())).\
		assign(idx = np.repeat(np.arange(0,tstates.shape[2]), tstates.shape[0]))
	dat = pd.DataFrame(states.T,
			columns=list(cols.keys())).\
		assign(idx = np.arange(0,states.shape[1]))
	simple_dat = tdat[['idx', 't', 'v']]
	simple_dat = pd.concat(
		[simple_dat.query('t == 0').\
			assign(t = 0),	
		simple_dat.merge(
					pd.DataFrame({
						'idx': np.arange(config['N']),
						't': np.round(sample_t).astype(int)}),
					how='inner', on=['idx', 't']).\
			assign(t = 1)]).\
		pivot(index='idx', columns='t', values='v').\
		reset_index().\
		sort_values(by='idx').\
		assign(t = np.round(sample_t).astype(int))\
		[['t', 0, 1]]
	simple_dat.columns=['t', 'v0', 'v1']
	if not os.path.exists('output'):
		os.makedirs('output')
	os.makedirs("output/simulations", exist_ok=True)
	tdat.to_csv(f"output/simulations/{config['out_file']}_tstates.tsv.gz",
		sep='\t', index=None, compression='gzip')
	dat.to_csv(f"output/simulations/{config['out_file']}_states.tsv.gz",
		sep='\t', index=None, compression='gzip')
	simple_dat.to_csv(f"output/simulations/{config['out_file']}_simulation_clean.tsv",
		sep='\t', index=None)
	if 'figure' in config.keys() and config['figure'].lower() == "true":
		# finally, plot simulation results
		fig = plot_vl_trajectories(tstates[:,[tcols['t'],tcols['v']],:])
		os.makedirs("figures/simulations", exist_ok=True)
		fig = fig.savefig(f"figures/simulations/{config['out_file']}.pdf")
		


#@profile
def run():
	parser = argparse.ArgumentParser()
	parser.add_argument('--config',
		help='config file')
	parser.add_argument('--params',
		help='params file')
	args = parser.parse_args()
	#args.config = 'config/tmp_sampled_sim_config.csv'
	#args.params = 'config/sim_params.csv'
	#args.config = 'config/fit_config.csv'
	#args.params = 'config/tmp_fit_params.csv'
	
	#### SET UP ANALYSIS ####
	config = {i[0]: eval(i[1]) for 
		i in pd.read_csv(args.config, comment='#').values}

	import os
	os.makedirs('/'.join(config['out_file'].split('/')[:-1]), exist_ok=True)
	#config['theta'] = [0.003]
	# create functions to generate parameters
	params_generator = {i[0]: eval(i[1]) for 
		i in pd.read_csv(args.params, comment='#').values}

	if config['mode'] == 'inference':
		#### READ IN DATA ####
		data = pd.read_csv(config['data'], sep='\t' if config['data'][-4:] == '.tsv' else ',')
		print(data)
		data['dt'] = (data.t / config['dt']).astype(int)
		data['unique_dt_idx'] = np.digitize(data['dt'], np.unique(data['dt']))-1
		# get data deciles
		bd_vl = np.log10(150)
		data_quantiles = np.hstack([
			bd_vl,
			np.quantile(data.query('v1 > @bd_vl').v1,
				np.arange(0,1.1,0.1))])
		data['v1_bin'] = np.digitize(data.v1, data_quantiles, right=True)
		# wasteful to repeat this for each row
		# but saves sending multiple objects to run_iteration function
		# shrug
		data['bin_edges'] = pd.Series([data_quantiles]*data.shape[0])
		data['unique_dt'] = pd.Series([np.unique(data.dt)]*data.shape[0])
		config['t_max'] = data.t.max().astype(int) + 1
		mcmc(data,params_generator, config)
	elif config['mode'] == 'calculation':
		# flag: combine data code with above
		#### READ IN DATA ####
		data = pd.read_csv(config['data'], sep='\t' if config['data'][-4:] == '.tsv' else ',')
		data['dt'] = (data.t / config['dt']).astype(int)
		data['unique_dt_idx'] = np.digitize(data['dt'], np.unique(data['dt']))-1
		# flag: move this to data processing file
		#data.v0 = np.log10(data.v0)
		#data.v1 = np.log10(data.v1)
		# get data deciles
		bd_vl = np.log10(150)
		data_quantiles = np.hstack([
			bd_vl,
			np.quantile(data.query('v1 > @bd_vl').v1,
				np.arange(0,1.1,0.1))])
		# flag: why do I do this? 
		#data_quantiles[0] = -np.inf
		data['v1_bin'] = np.digitize(data.v1, data_quantiles, right=True)
		# wasteful to repeat this for each row
		# but saves sending multiple objects to run_iteration function
		# shrug
		data['bin_edges'] = pd.Series([data_quantiles]*data.shape[0])
		data['unique_dt'] = pd.Series([np.unique(data.dt)]*data.shape[0])
		config['t_max'] = data.t.max().astype(int) + 1
		config['random_draws'] = config['rng'].uniform(size=config['n_random_draws'])
		params = {config['theta_labels'][idx]: i for idx, i in enumerate(config['theta'])}
		for k, v in params_generator.items():
			params[k] = v(data, params, config)
		full_theta = [stickbreak(i) if 
				isinstance(i, np.ndarray) else i for i in config['theta']]
		llh, simulated_bin_probs, config = run_iteration(full_theta, data, params_generator, config)
		with open(config['out_file'] + '_llh.tsv', 'w') as fp:
			_ = fp.write('likelihood\t' + '\t'.join(config['theta_labels']) + '\tbin_edges\tbin_probs\n')
			_ = fp.write(f'{llh}\t{'\t'.join([str(i) for i in config['theta']])}\t{';'.join([str(i) for i in data['bin_edges'].iloc[0]])}\t{';'.join([str(i) for i in simulated_bin_probs])}\n')
	elif config['mode'] == 'simulation':
		out_file_base = config['out_file']
		for idx in range(config['n_iter']):
			config['out_file'] = out_file_base + str(idx)
			simulate(params_generator, config)


'''
fig, axs = plt.subplots()
axs.hist(data.t, bins=10)
axs.set_xlabel('time between samples')
axs.set_ylabel('participants')
fig.savefig('figures/unique_llv_pairs_delta_t.pdf')
plt.close()

data = pd.read_csv(args.data, sep='\t')

data.query('v1 == 0.0').shape

data = pd.read_csv('data/unique_llv_pairs.csv')

'''



if __name__ == "__main__":
    run()


