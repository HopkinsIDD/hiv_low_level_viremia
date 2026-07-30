/// NEED TO REVIEW VIRAL LOAD MODEL AND ENUSRE THAT THERE ARE NO FALSE NEGATIVES FOR TREATMENT DISTRIBUTIONS!!!!!
functions{
  real calc_log_prob_no_intra(int trans_form, real use_param1, real use_param2, real use_param3, real t, real log10_v){
    if (trans_form == 1){
      return -108*t*use_param1*(1 - pow(1 - use_param2, use_param3*pow(10, log10_v)));
    }else if (trans_form == 2){
      return t .* (use_param1 / (1 + pow(10, -use_param2 * (log10_v - use_param3))));
    }else{
      return 0;
    }
  }
  vector normalize_log_probs(vector x){
    real tot_log_prob = log_sum_exp(x);
    return x - tot_log_prob;
  }
}


data {
  // Individual-level viral load model
  int<lower=1> N_ind;   // number of individuals with viral load data (including pre-treatment viral load model and transmission rate model)
  int<lower=1> N; // number of viral load measurements
  array[N] int ind_idx; // index of individual assigned to each viral load measurement
  array[N] real<lower = 0> v; // viral load measurements (copies/mL), below detection = 1
  array[N] int<lower = 0, upper = 1> tx_avail; // boolean indicating if v is sampled during treatment use
  vector[N] prev_putative_tx; // individual level cumulative number of below detection viral loads that could have been due to unreported treatment
  array[N] int<lower=0, upper=1> tx; // boolean indicating self-reported treatment
  // Transmission model
  int<lower=1> P; // number of time periods over which transmission could have occurred
  int<lower=1> Q; // maximum number of viral loads per time period (code at present only tested with Q=2)
  array[P] int<lower=1, upper=Q> q; // number of observations per period
  array[P,Q] int<lower=-1> period_donor_v_idx; // donor viral load index in the viral load arrays above
  array[P,Q] real period_t; // duration of each time period
  array[N] int sex; // 1 = Female, 2 = male
  array[P] int<lower = 0, upper=1> period_y; // whether transmission occurred or didn't occurr in that period
  // Run settings
  int<lower = 0, upper = 1> sample_from_posterior; // whether to sample from the posterior or the prior (0 = prior, 1 = posterior)
  real lod; // limit of detection
  real imputed_bd_ll; // lower limit of imputed below-detection viral laods
  real imputed_bd_ul; // upper limit of imputed below-detection viral loads for each coupling period
  int<lower=1, upper=2> trans_vl; // whether to use all proximal viral loads (1) or just pre-treatment proximal viral loads (2)
  int<lower=0, upper=1> constrain_prob_intra_given_transmission; // whether to constrain probability of intra couple given transmission
  int<lower=0> N_constrained; // how many transmission events to constrain
  array[N_constrained] int period_constrain_idx; // period indices to use for constraint
  real N_constrained_linked; // of those couples with linkage how many are linked
  int<lower=0, upper=1> allow_saturation; // whether to all saturation of transmission risk
  int<lower=1, upper=2> trans_form; // functional form of the VL v. transmission risk relationship. 1: Virion-establishment, 2: Hill
  real<lower=-1> c; // fixed value for scaling factor c in birth-death form, if <0 then estimated
}


transformed data{
  // ---------------- //
  // PARAMETER LIMITS //
  // ---------------- //
  real est_param1_lower = negative_infinity();
  real est_param1_upper = positive_infinity();
  real est_param2_lower = negative_infinity();
  real est_param2_upper = positive_infinity();
  real est_param3_lower = negative_infinity();
  real est_param3_upper = positive_infinity();
  //if (trans_form == 1){
  //  est_param3_lower = 0;
  //}
  if (trans_form == 2){
    est_param1_lower = 0;
    est_param2_lower = 0;
  }
  // ---------------- //
  // VIRAL LOAD MODEL //
  // ---------------- //
  // 1. Categorize observations
  // --> 1. above LOD
  // --> 2. below LOD
  int N_no_tx_non_bd = 0;
  int N_no_tx_bd = 0;
  for (i in 1:N){
    if (tx[i] == 0){
      if (v[i] > lod){
        N_no_tx_non_bd += 1;
      }
      else if (v[i] <= lod){
        N_no_tx_bd += 1;
      }
    }
  }

  array[N_no_tx_non_bd] int no_tx_non_bd_idx;
  array[N_no_tx_bd] int no_tx_bd_idx;
  int tmp1 = 1;
  int tmp2 = 1;
  for (i in 1:N){
    if (tx[i] == 0){
      if (v[i] > lod){
        no_tx_non_bd_idx[tmp1] = i;
        tmp1 += 1;
      }
      else if (v[i] <= lod){
        no_tx_bd_idx[tmp2] = i;
        tmp2 += 1;
      }
    }
  }

  if (
    (N_no_tx_non_bd + 
      N_no_tx_bd) != (N - sum(tx))){
    print("ERROR: TABULATED VIRAL LOADS DOES NOT MATCH EXPECTED SIZE");
    print("EXPECTED:");
    print(N - sum(tx));
    print("OBSERVED:");
    print((N_no_tx_non_bd + 
      N_no_tx_bd));
  }

  
  // 2. Index of individuals who ever have a below-detection observation with no treatment 
  array[N_ind] int ever_bd_non_tx_ind_bool = rep_array(0, N_ind);
  for (i in 1:N){
    if ((v[i] < lod) && (tx[i] == 0)){
      ever_bd_non_tx_ind_bool[ind_idx[i]] = 1;
    }
  }
  int N_ever_bd_non_tx_ind = sum(ever_bd_non_tx_ind_bool);
  array[N_ever_bd_non_tx_ind] int ever_bd_non_tx_ind_idx;
  tmp1 = 1;
  for (i in 1:N_ind){
    if (ever_bd_non_tx_ind_bool[i] == 1){
      ever_bd_non_tx_ind_idx[tmp1] = i;
      tmp1 += 1;
    }
  }

  // 3. Index all individuals with previous or no prvious putative tx
  int N_no_prev_putative_tx = 0;
  int N_any_prev_putative_tx = 0;
  for (i in 1:N){
    if (prev_putative_tx[i] == 0){
      N_no_prev_putative_tx  += 1;
    }else if (prev_putative_tx[i] > 0){
      N_any_prev_putative_tx += 1;
    }
  }
  array[N] int<lower=0, upper=1> no_prev_putative_tx = rep_array(0,N);
  array[N_no_prev_putative_tx] int no_prev_putative_tx_idx;
  array[N] int<lower=0, upper=1> any_prev_putative_tx = rep_array(0,N);
  array[N_any_prev_putative_tx] int any_prev_putative_tx_idx;
  tmp1 = 1;
  tmp2 = 1;
  for (i in 1:N){
    if (prev_putative_tx[i] == 0){
      no_prev_putative_tx[i] = 1;
      no_prev_putative_tx_idx[tmp1] = i;
      tmp1 += 1;
    }else if (prev_putative_tx[i] > 0){
      any_prev_putative_tx[i] = 1;
      any_prev_putative_tx_idx[tmp2] = i;
      tmp2 += 1;
    }
  }


  // 4. Convert viral load to log scale
  real log10_lod = log10(lod); // limit of detection 
  array[N] real log10_v = log10(v); // pre-treatment viral load measurements (log10 copies/mL)
  
  // ----------------------- //
  // TRANSMISSION RATE MODEL //
  // ----------------------- //
  // 1. Categorize below-detection donor viral loads according to three categories
  // --> 2. Below LOD, no reported treatment, donor_v_no_tx_bd
  // --> 4. Below LOD, tx reported, donor_v_tx_bd
  int N_donor_v_no_tx_bd = 0;
  int N_donor_v_tx_bd = 0;
  for (i in 1:P){
    for (j in 1:q[i]){
      if (v[period_donor_v_idx[i,j]] <= lod){
        if (tx[period_donor_v_idx[i,j]] == 0){
          N_donor_v_no_tx_bd += 1;
        }
        else if (tx[period_donor_v_idx[i,j]] == 1){
          N_donor_v_tx_bd += 1;
        }
      }
    }
  }
  // idx maps positions in the N_donor_v_non_tx_bd-length array to the N-length array
  // rev_idx maps positions in the N-length array to positions in the N_donor_v_non_tx_bd array
  array[N] int donor_v_no_tx_bd_rev_idx;
  array[N_donor_v_no_tx_bd] int donor_v_no_tx_bd_idx;
  array[N] int donor_v_tx_bd_rev_idx;
  array[N_donor_v_tx_bd] int donor_v_tx_bd_idx;
  tmp1 = 1;
  tmp2 = 1;
  for (i in 1:P){
    for (j in 1:q[i]){
      if (v[period_donor_v_idx[i,j]] <= lod){
        if (tx[period_donor_v_idx[i,j]] == 0){
          donor_v_no_tx_bd_rev_idx[period_donor_v_idx[i,j]] = tmp1;
          donor_v_no_tx_bd_idx[tmp1] = period_donor_v_idx[i,j];
          tmp1 += 1;
        }
        else if (tx[period_donor_v_idx[i,j]] == 1){
          donor_v_tx_bd_rev_idx[period_donor_v_idx[i,j]] = tmp2;
          donor_v_tx_bd_idx[tmp2] = period_donor_v_idx[i,j];
          tmp2 += 1;
        }
      }
    }
  }


  // 2. Categorizes periods based on overlapping donor viral load categories
  // --> 1. above detection (uses VL as observed): period_v_ad
  // --> 3. below detection, prior to self-reported treatment (impute from intra-individual VL distribution or TX-BD distribution): period_v_no_tx_bd
  // --> 4. below detection, after self-reported treatment initiation (impute from TX-BD distribution): period_v_tx_bd
  int P_period_v_ad = 0;
  int P_period_v_no_tx_bd = 0;
  int P_period_v_tx_bd = 0;
  for (i in 1:P){
    for (j in 1:q[i]){
      if (v[period_donor_v_idx[i,j]] > lod){
        P_period_v_ad += 1;
      }
      else if (v[period_donor_v_idx[i,j]] <= lod){
        if (tx[period_donor_v_idx[i,j]] == 0){
          P_period_v_no_tx_bd += 1;
        }else if (tx[period_donor_v_idx[i,j]] == 1){
          P_period_v_tx_bd += 1;
        }
      }
    }
  }

  array[P_period_v_ad,2] int period_v_ad_idx; 
  array[P_period_v_no_tx_bd,3] int period_v_no_tx_bd_idx; 
  array[P_period_v_tx_bd,2] int period_v_tx_bd_idx; 
  tmp1 = 1;
  tmp2 = 1;
  int tmp3 = 1;
  for (i in 1:P){
    for (j in 1:q[i]){
      if (v[period_donor_v_idx[i,j]] > lod){
        period_v_ad_idx[tmp1,1] = i;
        period_v_ad_idx[tmp1,2] = j;
        tmp1 += 1;
      }
      else if (v[period_donor_v_idx[i,j]] <= lod){
        if (tx[period_donor_v_idx[i,j]] == 0){
          period_v_no_tx_bd_idx[tmp2,1] = i;
          period_v_no_tx_bd_idx[tmp2,2] = j;
          period_v_no_tx_bd_idx[tmp2,3] = tmp2;
          tmp2 += 1;
        }
        else if (tx[period_donor_v_idx[i,j]] == 1){
          period_v_tx_bd_idx[tmp3,1] = i;
          period_v_tx_bd_idx[tmp3,2] = j;
          tmp3 += 1;
        }
      }
    }
  }

  // 3. Categorize periods based on any previous putative unreported treatment observations
  array[P,Q] int<lower=0, upper=1> period_any_prev_putative_tx = rep_array(0, P, Q);
  array[P,Q] int<lower=0, upper=1> period_tx_avail = rep_array(0, P, Q);
  array[P,Q] int<lower=0, upper=1> period_tx = rep_array(0, P, Q);
  for (i in 1:P){
    for (j in 1:q[i]){
      if (prev_putative_tx[period_donor_v_idx[i,j]] > 0){
        period_any_prev_putative_tx[i,j] = 1;
      }else{
        period_any_prev_putative_tx[i,j] = 0;
      }
      if (tx_avail[period_donor_v_idx[i,j]] == 1){
        period_tx_avail[i,j] = 1;
      }else{
        period_tx_avail[i,j] = 0;
      }
      if (tx[period_donor_v_idx[i,j]] == 1){
        period_tx[i,j] = 1;
      }else{
        period_tx[i,j] = 0;
      }
    }
  }

  if ((P_period_v_ad + P_period_v_no_tx_bd + P_period_v_tx_bd) != sum(q)){
    print("ERROR: TABULATED PERIODS DOES NOT MATCH EXPECTED SIZE");
    print("EXPECTED:");
    print(P);
    print("OBSERVED (IND):");
    print(P_period_v_ad);
    print(P_period_v_no_tx_bd);
    print(P_period_v_tx_bd);
    print("OBSERVED (TOTAL):");
    print((P_period_v_ad + P_period_v_no_tx_bd + P_period_v_tx_bd));
  }


  // 4. Categorize periods by outcome
  int P_y1 = sum(period_y); // total number of seroconversions
  array[P_y1] int y1_idx; // index where seroconversiosn occurred
  array[P-P_y1] int y0_idx; // index where seroconversions did not occur
  tmp1 = 1;
  tmp2 = 1;
  for (i in 1:P){
    if (period_y[i] == 1){
      y1_idx[tmp1] = i;
      tmp1 += 1;
    }else if (period_y[i] == 0){
      y0_idx[tmp2] = i;
      tmp2 += 1;
    }
  }


  // 5. Sum total time across periods
  vector[P] tot_period_t = rep_vector(0, P);
  for (i in 1:P){
    for (j in 1:q[i]){
      tot_period_t[i] += period_t[i,j];
    }
  }


  // 6. Summarize linkage constraints
  real N_constrained_unlinked_plus_one = N_constrained - N_constrained_linked + 1;
  real N_constrained_linked_plus_one = N_constrained_linked + 1;

}



parameters {
  // viral load model
  real<lower=1> mu_0; // mean of individual mean VLs
  real<lower=0> omega_mu; // individual random effects standard deviation for mean VL
  vector[N_ind] eta_mu_i; // individual random effects for mean VL
  real<lower=0> sigma_0; // mean of individual VL standard deviation
  real<lower=0> omega_sigma; // individual random effects standard deviation for VL standard deviation
  vector[N_ind] eta_sigma_i; // individual random effects for VL standard deviation  
  real logit_fnr; // logit_false negative rate
  real logit_txr; // logit probability of unreported treatment initialization 
  real<lower=log10_lod> mu_tx; // mean log10 viral load of non-suppressed treatment-experienced
  real<lower=0> sigma_tx; // variation of non-suppressed treatment-experienced
  real logit_prev_tx_bdr; // mixing parameter between BD and N for those previously on tx
  // transmission model parameters
  real<lower=est_param1_lower, upper=est_param1_upper> est_param1;
  real<lower=est_param2_lower, upper=est_param2_upper> est_param2;
  real<lower=est_param3_lower, upper=est_param3_upper> est_param3;
  sum_to_zero_vector[2] est_param1_coeffs;
  sum_to_zero_vector[2] est_param2_coeffs;
  real<lower=0> beta_extra; // rate of extra-couple transmission
  real<lower=0> beta_tx; // VL-independent transmission rate post-treatment (only if trans_vl == 2)
  vector[N_donor_v_no_tx_bd] 
    donor_no_tx_bd_sampled_log10_v; // sampled log10 viral load associated with BD observations with tx avail, truncated at LOD
  vector<upper=log10_lod>[N_donor_v_no_tx_bd] 
    donor_no_tx_bd_trunc_sampled_log10_v; // sampled log10 viral load associated with BD observations with tx avail, truncated at LOD
  vector<lower=imputed_bd_ll, upper=imputed_bd_ul>[N_donor_v_no_tx_bd] 
    donor_no_tx_bd_imputed_v;  // imputed log10 viral load associated with BD observations with tx_avail, truncated at LOD
  vector<lower=imputed_bd_ll, upper=imputed_bd_ul>[N_donor_v_tx_bd] 
    donor_tx_bd_imputed_v; // imputed log10 viral load associated with BD observations with reported tx, truncated at LOD
    
}


transformed parameters{
  // ----------------------------- //
  // SHARED TRANSFORMED PARAMETERS //
  // ----------------------------- //
  real log_prob_fn = log_inv_logit(logit_fnr); // log probability of false negative 
  real log_prob_not_fn = log1m_exp(log_prob_fn); // log probability of not having a false negative 
  
  real log_prob_tx = log_inv_logit(logit_txr); // log probability of initiating treatment 
  real log_prob_not_tx = log1m_exp(log_prob_tx); // log probability of not initiating threatment

  real log_prob_bd_given_prev_tx = log_inv_logit(logit_prev_tx_bdr); // log probability of BD given previous unreported treatment
  real log_prob_ad_given_prev_tx = 
    log1m_exp(log_prob_bd_given_prev_tx) - 
      log1m_exp(normal_lcdf(log10_lod | mu_tx, sigma_tx)); // log probability of AD given previous unreported treatment, with normalization constant

  real log_prob_tx_or_fn = log_sum_exp([
    log_prob_tx,
    log_prob_fn,
    log_prob_tx + log_prob_fn]); // log probability of treatment initiation or false negative
  real log_prob_not_tx_not_fn = log_prob_not_tx + log_prob_not_fn; // log probability of not treatment initiatiation or false negative
  real log_prob_not_tx_fn = log_prob_not_tx + log_prob_fn; // logit probability of not treatment initiation but false negative

  vector[N_ind] mu_i = mu_0 + omega_mu .* eta_mu_i; // individual-level VL mean
  vector[N_ind] sigma_i = sigma_0*exp(omega_sigma .* eta_sigma_i); // individual-level VL standard deviation
  vector[N_ind] ind_log_prob_bd_given_no_tx_no_fn; // individual-level log-probability of true below-detection without treatment
  vector[N_ind] ind_log_prob_tx_given_bd;
  for (i in ever_bd_non_tx_ind_idx){
    ind_log_prob_bd_given_no_tx_no_fn[i] = normal_lcdf(log10_lod | mu_i[i], sigma_i[i]);
    ind_log_prob_tx_given_bd[i] = log_prob_tx - 
      log_sum_exp([
        log_prob_tx,
        log_prob_not_tx_fn,
        log_prob_not_tx_not_fn + ind_log_prob_bd_given_no_tx_no_fn[i]]);
  }

  vector[N] log_prob_ever_tx;
  vector[N] log_prob_never_tx;
  log_prob_ever_tx[any_prev_putative_tx_idx] = 
    log1m_exp(
      prev_putative_tx[any_prev_putative_tx_idx] .* 
        log1m_exp(ind_log_prob_tx_given_bd[ind_idx][any_prev_putative_tx_idx]));
  
  log_prob_never_tx[no_prev_putative_tx_idx] = rep_vector(0, N_no_prev_putative_tx); // no prev putative tx, thus prob_never_tx = 1
  log_prob_never_tx[any_prev_putative_tx_idx] = log1m_exp(log_prob_ever_tx[any_prev_putative_tx_idx]); // if some prev puptative tx, prob_never_tx <= 1
  
  // ----------------------- //
  // TRANSMISSION RATE MODEL //
  // ----------------------- //
  // DEFINE TRANSFORMED PARAMETERS:
  // 1. Parameter transformations
  vector[N] use_param1 = est_param1 + est_param1_coeffs[sex];
  vector[N] use_param2 = est_param2 + est_param2_coeffs[sex];
  real use_param3 = est_param3;
  if (trans_form == 1){
    use_param1 = inv_logit(use_param1);
    use_param2 = inv_logit(use_param2);
    use_param3 = exp(use_param3);
    if (allow_saturation == 0){
      use_param1 = rep_vector(1, N);
    }
    if (c > 0){
      use_param3 = c;
    }
  }

  // 2. Transform imputed values
  vector[N_donor_v_no_tx_bd] donor_no_tx_bd_imputed_log10_v = 
    log10(donor_no_tx_bd_imputed_v);
  vector[N_donor_v_tx_bd] donor_tx_bd_imputed_log10_v = 
    log10(donor_tx_bd_imputed_v);

  // 3. Temporary values
  matrix[3, P_period_v_no_tx_bd] norm_log_prob_bd; 

  // 4. Probability of no intra-couple transmission
  vector[P] log_prob_no_intra_transmission = rep_vector(0, P);
  vector[P] log_prob_no_transmission = rep_vector(0,P);
  vector[P] prob_intra_given_transmission; 
  real cum_prob_intra_given_transmission = 0; // cumulative probability of intra-couple transmission given transmission occurred
  for (idx in period_v_ad_idx){
    // if above detection, use VL as observed
    // period_v is on the copies/ml scale, so no need for 10^ transformation
    log_prob_no_intra_transmission[idx[1]] += 
      ((trans_vl == 2) && (period_tx[idx[1], idx[2]] == 1)) ? 
        -beta_tx*period_t[idx[1], idx[2]] :
        calc_log_prob_no_intra(
          trans_form,
          use_param1[period_donor_v_idx[idx[1], idx[2]]], use_param2[period_donor_v_idx[idx[1], idx[2]]], use_param3,
          period_t[idx[1], idx[2]],
          log10_v[period_donor_v_idx[idx[1], idx[2]]]);
  }
  for (idx in period_v_no_tx_bd_idx){
    // below-detection with on reported ART
    // --> prior to unreported treatment:
      // --> false negative, then impute from intra-individual VL distirbution
      // --> true negative, then impute from intra-individual VL distribution truncated at LOD
    // --> after unreported treatment
      // --> then impute ffrom tx-bd distribution
    // does this observed donor viral load occur after any previous putative tx
    // normalized probabilities of three potetial data-generating processes
    // dropping treatment initiatialization rate prior to round 10
    // false negatives
    norm_log_prob_bd[1, idx[3]] =
      log_prob_never_tx[period_donor_v_idx[idx[1],idx[2]]] + 
        [log_prob_fn, log_prob_not_tx_fn][(1 + period_tx_avail[idx[1], idx[2]])];
    // true negative with no treatment
    norm_log_prob_bd[2, idx[3]] = 
      log_prob_never_tx[period_donor_v_idx[idx[1],idx[2]]] + 
        [log_prob_not_fn, log_prob_not_tx_not_fn][(1 + period_tx_avail[idx[1], idx[2]])] +
        ind_log_prob_bd_given_no_tx_no_fn[ind_idx[period_donor_v_idx[idx[1], idx[2]]]];
    // true negatives with unreported treatment
    if (period_tx_avail[idx[1], idx[2]] == 1){
      if (period_any_prev_putative_tx[idx[1], idx[2]] == 1){
        norm_log_prob_bd[3, idx[3]] = log_sum_exp([
          log_prob_never_tx[period_donor_v_idx[idx[1],idx[2]]] + 
            log_prob_tx,
          log_prob_ever_tx[period_donor_v_idx[idx[1],idx[2]]] + 
            log_prob_bd_given_prev_tx]);
      }else{
        norm_log_prob_bd[3, idx[3]] = log_prob_tx;
      }
    }
    //norm_log_prob_bd[idx[3],:2] = normalize_log_probs(norm_log_prob_bd[idx[3],:2]')';
    norm_log_prob_bd[:(2 + (period_tx_avail[idx[1], idx[2]] == 1)), idx[3]] = 
      normalize_log_probs(norm_log_prob_bd[:(2 + (period_tx_avail[idx[1], idx[2]] == 1)), idx[3]]);
    // calculated weighted sum of potential probabilities of no intra-couple transmission
    // drop term describing suppression on unreported treatment when that is not possible
    if (trans_vl == 1){
      log_prob_no_intra_transmission[idx[1]] += 
        log_sum_exp([
          norm_log_prob_bd[1, idx[3]] + 
            calc_log_prob_no_intra(
              trans_form,
              use_param1[period_donor_v_idx[idx[1], idx[2]]], use_param2[period_donor_v_idx[idx[1], idx[2]]], use_param3,
              period_t[idx[1], idx[2]],
              donor_no_tx_bd_sampled_log10_v[donor_v_no_tx_bd_rev_idx[period_donor_v_idx[idx[1], idx[2]]]]),
           norm_log_prob_bd[2, idx[3]] +
            calc_log_prob_no_intra(
              trans_form,
              use_param1[period_donor_v_idx[idx[1], idx[2]]], use_param2[period_donor_v_idx[idx[1], idx[2]]], use_param3,
              period_t[idx[1], idx[2]],
              donor_no_tx_bd_trunc_sampled_log10_v[donor_v_no_tx_bd_rev_idx[period_donor_v_idx[idx[1], idx[2]]]]),
            ((period_tx_avail[idx[1], idx[2]] == 1) ? 
              norm_log_prob_bd[3, idx[3]] + 
                calc_log_prob_no_intra(
                  trans_form,
                  use_param1[period_donor_v_idx[idx[1], idx[2]]], use_param2[period_donor_v_idx[idx[1], idx[2]]], use_param3,
                  period_t[idx[1], idx[2]],
                  donor_no_tx_bd_imputed_log10_v[donor_v_no_tx_bd_rev_idx[period_donor_v_idx[idx[1], idx[2]]]]) : 
              0)][1:(2 + period_tx_avail[idx[1], idx[2]])]);
      }else if (trans_vl == 2){
        log_prob_no_intra_transmission[idx[1]] += 
          log_sum_exp([
            norm_log_prob_bd[1, idx[3]] + 
              calc_log_prob_no_intra(
                trans_form,
                use_param1[period_donor_v_idx[idx[1], idx[2]]], use_param2[period_donor_v_idx[idx[1], idx[2]]], use_param3,
                period_t[idx[1], idx[2]],
                donor_no_tx_bd_sampled_log10_v[donor_v_no_tx_bd_rev_idx[period_donor_v_idx[idx[1], idx[2]]]]),
             norm_log_prob_bd[2, idx[3]] +
              calc_log_prob_no_intra(
                trans_form,
                use_param1[period_donor_v_idx[idx[1], idx[2]]], use_param2[period_donor_v_idx[idx[1], idx[2]]], use_param3,
                period_t[idx[1], idx[2]],
                donor_no_tx_bd_trunc_sampled_log10_v[donor_v_no_tx_bd_rev_idx[period_donor_v_idx[idx[1], idx[2]]]]),
            ((period_tx_avail[idx[1], idx[2]]>0) ? 
              norm_log_prob_bd[3, idx[3]] + 
                -beta_tx*period_t[idx[1], idx[2]] : 
              0)][1:(2 + period_tx_avail[idx[1], idx[2]])]);
      }
  }
  for (idx in period_v_tx_bd_idx){
    // below-detection during self reported treatment
      // --> then impute ffrom tx-bd distribution
    log_prob_no_intra_transmission[idx[1]] += 
      (trans_vl == 2) ?
        -beta_tx*period_t[idx[1], idx[2]] :
        calc_log_prob_no_intra(
          trans_form,
          use_param1[period_donor_v_idx[idx[1], idx[2]]], use_param2[period_donor_v_idx[idx[1], idx[2]]], use_param3,
          period_t[idx[1], idx[2]],
          donor_tx_bd_imputed_log10_v[donor_v_tx_bd_rev_idx[period_donor_v_idx[idx[1], idx[2]]]]);
  }
  log_prob_no_transmission = log_prob_no_intra_transmission + -1*tot_period_t*beta_extra;
  if (constrain_prob_intra_given_transmission == 1){
    for (i in period_constrain_idx){
      if (period_y[i] == 1){
        prob_intra_given_transmission[i] = log1m_exp(log_prob_no_intra_transmission[i]) - log1m_exp(log_prob_no_transmission[i]);
        cum_prob_intra_given_transmission += exp(prob_intra_given_transmission[i]);
      }
   }
  }
}


model {  
  // ------ // 
  // PRIORS //
  // ------ //
  mu_0 ~ normal(4,2); // population mean viral load
  omega_mu ~ normal(0, 1); // mean viral load random effects std. dev
  eta_mu_i ~ normal(0, 1); // mean viral load random effects. non-centered parameterization so must be standard-normal
  sigma_0 ~ normal(0,1); // population mean intra-individual variability
  omega_sigma ~ normal(0,1); // mean intra-individual variability random effects std. dev
  eta_sigma_i ~ normal(0,1); // viral load variability random effects. non-centered parameterization so must be standard-normal
  logit_fnr ~ normal(0,1); // false negative rate
  logit_txr ~ normal(0,1); // logit probability of unreported treatment initialization 
  mu_tx ~ normal(4,2); // mean log10 viral load of non-suppressed treatment-experienced
  sigma_tx ~ normal(0,1); // variation of non-suppressed treatment-experienced
  logit_prev_tx_bdr ~ normal(0,1); // mixing parameter between BD and N for those previously on tx
  
  // ----------------------- //
  // TRANSMISSION RATE MODEL //
  // ----------------------- //
  est_param1 ~ normal(0,5);
  est_param2 ~ normal(0,5);
  est_param3 ~ normal(0,1);
  est_param1_coeffs ~ normal(0,1);
  est_param2_coeffs ~ normal(0,1);

  beta_tx ~ normal(0, 0.01);
  if (constrain_prob_intra_given_transmission == 1){
    // Wiki: expected value of the posterior distribution over p, namely Beta(s+1, n−s+1)
    (cum_prob_intra_given_transmission/N_constrained) ~ beta(N_constrained_linked_plus_one, N_constrained_unlinked_plus_one);
  }else{
    beta_extra ~ normal(0, 0.01);
  }

  // prior for imputed BD values
  donor_no_tx_bd_sampled_log10_v ~ 
    normal(
      mu_i[ind_idx[donor_v_no_tx_bd_idx]],
      sigma_i[ind_idx[donor_v_no_tx_bd_idx]]);
  donor_no_tx_bd_trunc_sampled_log10_v ~ 
    normal(
      mu_i[ind_idx[donor_v_no_tx_bd_idx]],
      sigma_i[ind_idx[donor_v_no_tx_bd_idx]]);
  donor_no_tx_bd_imputed_v ~ uniform(imputed_bd_ll, imputed_bd_ul);
  donor_tx_bd_imputed_v ~ uniform(imputed_bd_ll, imputed_bd_ul);


  // LIKELIHOOD
  if (sample_from_posterior == 1){
    // ---------------- //
    // VIRAL LOAD MODEL //
    // ---------------- //
    // 1. above LOD
    // -> Can arise from
    // --> 1. No false negative, no treatment initiation, no prior treatment, intra-individual VL variation
    // --> 2. Intra-individual previous unreported treatment distribution
    for (idx in no_tx_non_bd_idx){
      target += log_sum_exp(
        [
          log_prob_never_tx[idx] + 
            [log_prob_not_fn, log_prob_not_tx_not_fn][1 + tx_avail[idx]] + 
            normal_lpdf(log10_v[idx] | mu_i[ind_idx[idx]], sigma_i[ind_idx[idx]]),
          ((any_prev_putative_tx[idx]==1) ? 
            log_prob_ever_tx[idx] + 
              log_prob_ad_given_prev_tx + 
                normal_lpdf(log10_v[idx] | mu_tx, sigma_tx) :
            0)
        ][1:(1 + (any_prev_putative_tx[idx] == 1))]);
    }
    // 2. below LOD
    // -> Can arise from:
    // --> 1. False negative or treatment initiation
    // --> 2. No false negative, no treatment initiation, intra-individual VL variation
    // --> 3. Intra-individual previous unreported treatment VL variation
    for (idx in no_tx_bd_idx){
      target += log_sum_exp(
        [
          [log_prob_fn, log_prob_tx_or_fn][1 + tx_avail[idx]],
          [log_prob_not_fn, log_prob_not_tx_not_fn][1 + tx_avail[idx]] +
            ind_log_prob_bd_given_no_tx_no_fn[ind_idx[idx]],
          ((any_prev_putative_tx[idx]==1) ? 
            log_prob_ever_tx[idx] + 
              log_prob_bd_given_prev_tx :
            0)
        ][1:(2 + (any_prev_putative_tx[idx] == 1))]);
    }

    // transmission rate model
    for (i in 1:P){
      if (period_y[i] == 1){
        target += log1m_exp(log_prob_no_transmission[i]);
      }else if (period_y[i] == 0){
        target += log_prob_no_transmission[i];
      }
    }
    
  }
}

generated quantities {
  real prob_bd_given_prev_tx = inv_logit(logit_prev_tx_bdr);
  real prob_tx = inv_logit(logit_txr);
  real prob_fn = inv_logit(logit_fnr);
  real mean_use_param1 = 0;
  real mean_use_param2 = 0;
  real mean_use_param3 = 0;

  real sex1_use_param1 = inv_logit(est_param1 + est_param1_coeffs[1]);
  real sex2_use_param1 = inv_logit(est_param1 + est_param1_coeffs[2]);
  real sex1_use_param2 = inv_logit(est_param2 + est_param2_coeffs[1]);
  real sex2_use_param2 = inv_logit(est_param2 + est_param2_coeffs[2]);
  real sex1_use_param3 = use_param3;
  real sex2_use_param3 = use_param3;

  for (i in 1:P){
    for (j in 1:q[i]){
      mean_use_param1 += period_t[i,j] * use_param1[period_donor_v_idx[i,1]];
      mean_use_param2 += period_t[i,j] * use_param2[period_donor_v_idx[i,1]];
    }
      
  }

  mean_use_param1 = mean_use_param1 / sum(tot_period_t);
  mean_use_param2 = mean_use_param2 / sum(tot_period_t);
  mean_use_param3 = use_param3;

}



