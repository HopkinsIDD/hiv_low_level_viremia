functions{
  real calc_log_prob_no_intra(int trans_form, real use_param1, real use_param2, real use_param3, real t, real log10_v){
    // given a functional form of the vl x transmission rate relationship (trans_form),
    // parameters of that functional form (use_param*),
    // duration of partnership (t),
    // and log10 viral load (log10_v)
    // returns the log probability of no intra-couple transmission
    if (trans_form == 1){
      // trans_form == 1: viral establishment model
      return -108*t*use_param1*(1 - pow(1 - use_param2, use_param3*pow(10, log10_v)));
    }else if (trans_form == 2){
      // trans_form == 2: sigmoidal (hill) function
      return -t .* (use_param1 / (1 + pow(10, -use_param2 * (log10_v - use_param3))));
    }else if (trans_form == 3){
      // don't have a name for this one 
      return -108*t*(1 - pow(inv(use_param1)/(inv(use_param1) + pow(10, log10_v)), use_param2));
    }else{
      return 0;
    }
  }

  vector normalize_log_probs(vector x){
    // given a vector of mutually exclusive probabilities,
    // normalizes to sum to 1
    real tot_log_prob = log_sum_exp(x);
    return x - tot_log_prob;
  }

  array[,] real calc_trunc_norm_log_probs(real mu, real sigma, real a, int N_b, vector b, array[] int use_b){
    // given mu and sigma and a lower truncation point (a) and a vector of intermediate truncation points (b)
    // and an array of which intermediate points to use  (use_b)
    // returns the log density that is above the lower truncation point (normalization constant)
    // and the normalized log density above below each of the intermediate truncation points
    array[N_b, 3] real out;
    real log_prob_lt_a = normal_lcdf(a | mu, sigma);
    for (i in 1:N_b){
      out[i,1] = log1m_exp(log_prob_lt_a);
      if (use_b[i]==1){
        out[i,2] = log_diff_exp(normal_lcdf(b[i] | mu, sigma), log_prob_lt_a) - out[i,1];
        out[i,3] = log1m_exp(out[i,2]);
      }
    }
    return(out);
  }
}





data {
  // Individual-level viral load model
  int<lower=1> N_ind;   // number of individuals with viral load data (including pre-treatment viral load model and transmission rate model)
  int<lower=1> N; // number of viral load measurements
  array[N] int ind_idx; // index of individual assigned to each viral load measurement
  array[N] real<lower = 0> v; // viral load measurements (copies/mL), below detection = 1
  int<lower=1> L; // number of unique limits of detection (LODs)
  vector[L] lod; // unique LODs
  array[N] int lod_idx; // limit of detection index for each viral load measurement
  array[N] int<lower = 0, upper = 1> tx_avail; // boolean indicating if v is sampled during treatment availability at a population-level
  array[N,L] int prev_putative_tx; // individual level cumulative number of below detection viral loads that could have been due to unreported treatment (e.g. below-detection) for each unique LOD
  array[N] int<lower=0, upper=1> tx; // boolean indicating self-reported treatment at time that viral load v was measured
  // Transmission model
  int<lower=1> P; // number of time periods over which transmission could have occurred
  int<lower=1> Q; // maximum number of viral loads per time period (code at present only tested with Q=2)
  array[P] int<lower=1, upper=Q> q; // number of observations per period
  array[P,Q] int<lower=-1> period_donor_v_idx; // donor viral load index in the viral load arrays above
  array[P,Q] real period_t; // duration of each time period
  array[P] int<lower = 0, upper=1> period_y; // whether transmission occurred or didn't occurr in that period
  // Run settings
  int<lower = 0, upper = 1> sample_from_posterior; // whether to sample from the posterior or the prior (0 = prior, 1 = posterior)
  real v_supr_min; // minimum viral load (copies/mL) of suppressed viral loads
  real v_supr_max; // maximum viral load (copies/mL) of suppressed viral loads
  int<lower=1, upper=2> trans_vl; // whether to use all proximal viral loads (1) or just pre-treatment proximal viral loads (2)
  int<lower=0, upper=1> constrain_prob_intra_given_transmission; // whether to constrain probability of intra couple given transmission
  int<lower=0> N_constrained; // how many transmission events to constrain
  array[N_constrained] int period_constrain_idx; // period indices to use for constraint
  real N_constrained_linked; // of those couples with linkage how many are linked
  int<lower=1, upper=3> trans_form; // functional form of the VL v. transmission risk relationship. 1: Virion-establishment, 2: Hill
  int<lower=0, upper=1> allow_saturation; // whether to all saturation of transmission risk (only used when trans_form = 1)
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
  if (trans_form == 2){
    est_param1_lower = 0;
    est_param2_lower = 0;
  }

  // ---------------- //
  // VIRAL LOAD MODEL //
  // ---------------- //
  // 1. Categorize observations
  // --> 1. No self-reported treatment, above LOD (no_tx_non_bd)
  // --> 2. No self-reported treatment, below LOD (no_tx_bd)
  // --> 3. Self-reported treatment, above LOD (tx_non_bd)
  // --> 4. Self-reported treatment, below LOD (tx_bd)
  // First, count number of occurrences for each category
  int N_no_tx_non_bd = 0;
  int N_no_tx_bd = 0;
  int N_tx_non_bd = 0;
  int N_tx_bd = 0;
  for (i in 1:N){
    if (tx[i] == 0){
      if (v[i] > lod[lod_idx[i]]){
        N_no_tx_non_bd += 1;
      }
      else if (v[i] <= lod[lod_idx[i]]){
        N_no_tx_bd += 1;
      }
    }else if (tx[i] == 1){
      if (v[i] > lod[lod_idx[i]]){
        N_tx_non_bd += 1;
      }
      else if (v[i] <= lod[lod_idx[i]]){
        N_tx_bd += 1;
      }
    }
  }

  // Get a list of viral load indices belonging to each category
  array[N_no_tx_non_bd] int no_tx_non_bd_idx;
  array[N_no_tx_bd] int no_tx_bd_idx;
  array[N_tx_non_bd] int tx_non_bd_idx;
  array[N_tx_bd] int tx_bd_idx;
  int tmp1 = 1;
  int tmp2 = 1;
  int tmp3 = 1;
  int tmp4 = 1;
  for (i in 1:N){
    if (tx[i] == 0){
      if (v[i] > lod[lod_idx[i]]){
        no_tx_non_bd_idx[tmp1] = i;
        tmp1 += 1;
      }
      else if (v[i] <= lod[lod_idx[i]]){
        no_tx_bd_idx[tmp2] = i;
        tmp2 += 1;
      }
    }else if (tx[i] == 1){
      if (v[i] > lod[lod_idx[i]]){
        tx_non_bd_idx[tmp3] = i;
        tmp3 += 1;
      }
      else if (v[i] <= lod[lod_idx[i]]){
        tx_bd_idx[tmp4] = i;
        tmp4 += 1;
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

  
  // 2. Index of individuals who ever have a below-detection observation with no self-reported treatment 
  // These are the individuals who may have been on unreported treatment
  // Count number of individuals who ever had a below-detection observation with no treatment
  array[N_ind] int ever_bd_non_tx_ind_bool = rep_array(0, N_ind);
  for (i in 1:N){
    if ((v[i] < lod[lod_idx[i]]) && (tx[i] == 0)){
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

  // 3. Index all viral load observations as previous or no previous putative treatment 
  // Putative treatment indicated by a below-detection viral load measurement during population-scale treatment availability prior to self-reported treatment
  // Count viral load observations with no BDs prior to self-reported treatment and those with BDs prior to self reported treatment
  int N_no_prev_putative_tx = 0;
  int N_any_prev_putative_tx = 0;
  for (i in 1:N){
    if (sum(prev_putative_tx[i,:]) == 0){
      N_no_prev_putative_tx  += 1;
    }else if (sum(prev_putative_tx[i,:]) > 0){
      N_any_prev_putative_tx += 1;
    }
  }
  // Get indices of individuals with no previous 
  array[N] int<lower=0, upper=1> no_prev_putative_tx = rep_array(0,N);
  array[N_no_prev_putative_tx] int no_prev_putative_tx_idx;
  array[N] int<lower=0, upper=1> any_prev_putative_tx = rep_array(0,N);
  array[N_any_prev_putative_tx] int any_prev_putative_tx_idx;
  tmp1 = 1;
  tmp2 = 1;
  for (i in 1:N){
    if (sum(prev_putative_tx[i,:]) == 0){
      no_prev_putative_tx[i] = 1;
      no_prev_putative_tx_idx[tmp1] = i;
      tmp1 += 1;
    }else if (sum(prev_putative_tx[i,:]) > 0){
      any_prev_putative_tx[i] = 1;
      any_prev_putative_tx_idx[tmp2] = i;
      tmp2 += 1;
    }
  }


 
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
      if (v[period_donor_v_idx[i,j]] <= lod[lod_idx[period_donor_v_idx[i,j]]]){
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
      if (v[period_donor_v_idx[i,j]] <= lod[lod_idx[period_donor_v_idx[i,j]]]){
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


  // 2. Categorizes periods based on non-overlapping donor viral load categories
  // --> 1. above detection (uses VL as observed): period_v_ad
  // --> 3. below detection, prior to self-reported treatment (impute from intra-individual VL distribution or unreported treatment VL distribution): period_v_no_tx_bd
  // --> 4. below detection, after self-reported treatment initiation (impute from unreported treatment distribution): period_v_tx_bd
  int P_period_v_ad = 0;
  int P_period_v_no_tx_bd = 0;
  int P_period_v_tx_bd = 0;
  for (i in 1:P){
    for (j in 1:q[i]){
      if (v[period_donor_v_idx[i,j]] > lod[lod_idx[period_donor_v_idx[i,j]]]){
        P_period_v_ad += 1;
      }
      else if (v[period_donor_v_idx[i,j]] <= lod[lod_idx[period_donor_v_idx[i,j]]]){
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
  tmp3 = 1;
  for (i in 1:P){
    for (j in 1:q[i]){
      if (v[period_donor_v_idx[i,j]] > lod[lod_idx[period_donor_v_idx[i,j]]]){
        period_v_ad_idx[tmp1,1] = i;
        period_v_ad_idx[tmp1,2] = j;
        tmp1 += 1;
      }
      else if (v[period_donor_v_idx[i,j]] <= lod[lod_idx[period_donor_v_idx[i,j]]]){
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


  // 3. Categorize periods based on overlapping previous putative unreported treatment observations categories
  // --> 1. any previous putative tx
  // --> 2. treatment available
  // --> 3. self-reported treatment
  // --> 4. n_norm_log_prob_bd:
  // ----> Number of components that could give rise to below-detection observations
  // ------> intra-individual variation, unreoprte treatment, reported treatment
  array[P,Q] int<lower=0, upper=1> period_any_prev_putative_tx = rep_array(0, P, Q);
  array[P,Q] int<lower=0, upper=1> period_tx_avail = rep_array(0, P, Q);
  array[P,Q] int<lower=0, upper=1> period_tx = rep_array(0, P, Q);
  array[P,Q] int<lower=0, upper=4> n_norm_log_prob_bd = rep_array(0, P, Q);
  for (i in 1:P){
    for (j in 1:q[i]){
      if (sum(prev_putative_tx[period_donor_v_idx[i,j]]) > 0){
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
      n_norm_log_prob_bd[i,j] = 2 + tx_avail[period_donor_v_idx[i,j]] + period_any_prev_putative_tx[i,j];
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
  // whether seroconversion occurs or not
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
  // Each period is spllit into "epochs" over which viral load is considered constant
  // But partner serostatus only measured at beginning and end of period
  vector[P] tot_period_t = rep_vector(0, P);
  for (i in 1:P){
    for (j in 1:q[i]){
       tot_period_t[i] += period_t[i,j];
    }
  }

  // 6. Limit of detection transformations
  array[N] real log10_v = log10(v);
  vector[L] log10_lod = log10(lod);
  real log10_v_supr_min = log10(v_supr_min);
  real log10_v_supr_max = log10(v_supr_max);
  // which limit of detections are > upper limit of suppression distribution
  // these are the only ones that have a non-zero probabilityof generating non-suppressed below-detection values
  int L_above_v_supr_max = 0;
  for (i in 1:L){
    if (lod[i] > v_supr_max){
      L_above_v_supr_max += 1;
    }
  }
  array[L_above_v_supr_max] int lod_above_v_supr_max_idx;
  array[L] int lod_above_v_supr_max_bool = rep_array(0, L);
  tmp1 = 1;
  for (i in 1:L){
    if (lod[i] > v_supr_max){
      lod_above_v_supr_max_idx[tmp1] = i;
      tmp1 += 1;
      lod_above_v_supr_max_bool[i] = 1;
    }
  }


  // 7. Summarize linkage constraints
  real N_constrained_unlinked_plus_one = N_constrained - N_constrained_linked + 1;
  real N_constrained_linked_plus_one = N_constrained_linked + 1;


}



parameters {
  // viral load model
  real<lower=1> mu_0; // mean of individual mean pre-treatment log10 VLs
  real<lower=0> omega_mu; // individual random effects standard deviation for mean pre-treatment log10 VL
  vector[N_ind] eta_mu_i; // individual random effects for mean pre-treatment log10 VL
  real<lower=0> sigma_0; // mean of individual pre-treatment log10 VL standard deviation
  real<lower=0> omega_sigma; // individual random effects standard deviation for log10 VL standard deviation
  vector[N_ind] eta_sigma_i; // individual random effects for log10 VL standard deviation  
  real<lower=log10_v_supr_max> mu_un_tx; // mean log10 viral load of non-suppressed unreported treatment
  real<lower=0> sigma_un_tx; // variation of non-suppressed unreported treatment log10 viral loads
  real<lower=log10_v_supr_max> mu_tx; // mean log10 viral load of non-suppressed treatment-experienced
  real<lower=0> sigma_tx; // variation of non-suppressed treatment-experienced log10 viral loads
  real logit_fnr; // logit_false negative rate (only those who are pre-treatment are modeled as having false negatives)
  real logit_txr; // logit probability of unreported treatment initialization (only considered when population-level treatment availability)
  real logit_theta_un_tx; // mixing parameter between suppressed and non-suppressed viral loads for those on unreported treatment
  real logit_theta_tx; // mixing parameter between suppressed and non-suppressed viral loads for those previously on tx
  

  // transmission model parameters
  real<lower=est_param1_lower, upper=est_param1_upper> est_param1;
  real<lower=est_param2_lower, upper=est_param2_upper> est_param2;
  real<lower=est_param3_lower, upper=est_param3_upper> est_param3;
  real<lower=0> beta_extra; // rate of extra-couple transmission
  real<lower=0> beta_tx; // VL-independent transmission rate post-treatment (only if trans_vl == 2 where we only consider pre-treatment viral loads in vl x transmission rate model)
  
  // imputed viral loads
  vector[N_donor_v_no_tx_bd] 
    donor_no_tx_bd_sampled_log10_v; // A. sampled log10 viral load associated with BD observations with no self-reported treatment attributed to false negatives
  // the below are all fit on the 0,1 scale and then scaled to the appropriate limits
  vector<lower=0, upper=1>[N_donor_v_no_tx_bd] 
    raw_donor_no_tx_bd_trunc_sampled_v; // B. raw sampled viral load associated with true BD observations with no treatment (reported or unreported), follow N(mu_i, sigma_i) right-truncated at LOD
  vector<lower=0, upper=1>[N_donor_v_no_tx_bd] 
    raw_donor_un_tx_bd_imputed_v;  // C. raw imputed viral load associated with true BD observations with unreported suppressive treatment, follow Unif(v_supr_min, v_supr_max)
  vector<lower=0, upper=1>[N_donor_v_no_tx_bd] 
    raw_donor_un_tx_bd_trunc_sampled_v;  // D. raw sampled viral load associated with true BD observations with unreported unsuppressive treatment, follow N(mu_un_tx, sigma_un_tx) right-truncated at LOD
  vector<lower=0, upper=1>[N_donor_v_tx_bd] 
    raw_donor_tx_bd_imputed_v;  // E. raw imputed viral load associated with true BD observations with reported suppressive treatment, follow Unif(v_supr_min, v_supr_max)
  vector<lower=0, upper=1>[N_donor_v_tx_bd] 
    raw_donor_tx_bd_trunc_sampled_v;  // F. raw sampled viral load associated with true BD observations with reported unsuppressive treatment, follow N(mu_tx, sigma_tx) right-truncated at LOD
 }


transformed parameters{
  // ----------------------------- //
  // SHARED TRANSFORMED PARAMETERS //
  // ----------------------------- //
  // MIXTURE MODEL MIXING PARAMETERS
  real fnr = inv_logit(logit_fnr);
  real log_prob_fn = log_inv_logit(logit_fnr); // log probability of false negative 
  real log_prob_not_fn = log1m_exp(log_prob_fn); // log probability of not having a false negative 
  
  real txr = inv_logit(logit_txr);
  real log_prob_tx_init = log_inv_logit(logit_txr); // log probability of initiating treatment 
  real log_prob_not_tx_init = log1m_exp(log_prob_tx_init); // log probability of not initiating threatment

  real log_theta_un_tx = log_inv_logit(logit_theta_un_tx); // log probability of suppression given previous unreported treatment
  real log1m_theta_un_tx = log1m_exp(log_theta_un_tx); // log probability of non-suppression given previous unreported treatment

  real log_theta_tx = log_inv_logit(logit_theta_tx); // log probability of suppression given reported treatment
  real log1m_theta_tx = log1m_exp(log_theta_tx); // log probability of non-suppression given previous unreported treatment


  // COMBINED PROBABILITIES OF TREATMENT INITIATION OR FALSE NEGATIVES
  // log probability of not treatment initiation but false negative
  real log_prob_not_tx_init_fn = log_prob_not_tx_init + log_prob_fn; 
  // log probability of not treatment initiatiation and not false negative
  real log_prob_not_tx_init_not_fn = log_prob_not_tx_init + log_prob_not_fn; 
  // log probability of treatment initiation or false negative
  // these are mutually exclusive
  // either initiate treatment with probability of log_prob_tx_init
  // or DON'T initiate treatment with probability log_prob_not_tx_init and have false negative with prob log_prob_fn
  real log_prob_tx_init_or_fn = 
    log_sum_exp([
      log_prob_tx_init,
      log_prob_not_tx_init_fn]); 
  
  // INDIVIDUAL-LEVEL PRE-TREATMENT VIRAL LOAD DISTRIBUTION
  vector[N_ind] mu_i = mu_0 + omega_mu .* eta_mu_i; // individual-level log10 VL mean when pre-treatment
  vector[N_ind] sigma_i = sigma_0*exp(omega_sigma .* eta_sigma_i); // individual-level log10 VL standard deviation  when pre-treatment
  // Probability of pre-treatment below detections and probability that any given below-detection indicates unreported treatment
  // Depends on both individual-level parameters (mu_i and sigma_i) and assay LOD
  // For each individual with any BDs while not self reporting treatment, calculate probabilities for each unique assay LOD
  array[N_ind,L] real ind_log_prob_bd_given_no_tx_no_fn; // individual-level log-probability of true below-detection without treatment for each unique LOD
  array[N_ind,L] real ind_log_prob_tx_init_given_bd; // individual-level log-probability of any given BD being due to treatment for each unique LOD
  for (i in ever_bd_non_tx_ind_idx){
    for (l in 1:L){
      ind_log_prob_bd_given_no_tx_no_fn[i,l] = normal_lcdf(log10_lod[l] | mu_i[i], sigma_i[i]);
      ind_log_prob_tx_init_given_bd[i,l] = log_prob_tx_init - 
        log_sum_exp([
          log_prob_tx_init,
          log_prob_not_tx_init_fn,
          log_prob_not_tx_init_not_fn + ind_log_prob_bd_given_no_tx_no_fn[i,l]]);
    }    
  }

  // PROBABILITY OF (N)EVER HAVING BEEN ON TREATMENT, CONDITIONAL ON ABOVE PROBABILITIES 
  // probability of not being on treatment based on signal for each unique LOD,
  //  P(never_a) = (1 - P(a))^N_a, logP(never_a) = N_a * log1m_exp(P(a))
  // then sum across LODs
  //  (P(never) = P(never_a)*P(never_b), logP(never) = logP(never_a) + logP(never_b))
  vector[N] log_prob_ever_tx;
  vector[N] log_prob_never_tx;
  for (i in any_prev_putative_tx_idx){
    log_prob_never_tx[i] = 
        sum(to_vector(prev_putative_tx[i,:]) .* 
          to_vector(log1m_exp(
            ind_log_prob_tx_init_given_bd[ind_idx[i],:])));
    log_prob_ever_tx[i] = log1m_exp(log_prob_never_tx[i]);
  }
   // no prev putative tx, thus prob_never_tx = 1
  log_prob_never_tx[no_prev_putative_tx_idx] = rep_vector(0, N_no_prev_putative_tx);

  // PROBABILITIES CONDITIONAL ON PAST OR CURRENT UNREPORTED TREATMENT
  // one row for each unique LOD
  // columns:
  // log density above log10_v_supr_max (shared across LOD)
  // log density between log10_v_supr_max and LOD
  // log density between LOD and infinity
  array[L,3] real log_probs_non_supr_given_un_tx =
    calc_trunc_norm_log_probs(mu_un_tx, sigma_un_tx, log10_v_supr_max, L, log10_lod, lod_above_v_supr_max_bool);
  // For each LOD, we can also calculate the log probability of being below-detection conditional on being on unreported treatment
  array[L] real log_prob_bd_given_un_tx;
  for (l in 1:L){
      log_prob_bd_given_un_tx[l] = 
        (lod_above_v_supr_max_bool[l] == 1) ? 
        log_sum_exp([log_theta_un_tx, log1m_theta_un_tx + log_probs_non_supr_given_un_tx[l,2]]) :
        log_theta_un_tx;
  }

  // PROBABILITIES CONDITIONAL ON BEING ON REPORTED TREATMENT
  // one row for each unique LOD
  // columns:
  // log density above log10_v_supr_max (shared across LOD)
  // log density between log10_v_supr_max and LOD
  // log density between LOD and infinity
  array[L,3] real log_probs_non_supr_given_tx =
    calc_trunc_norm_log_probs(mu_tx, sigma_tx, log10_v_supr_max, L, log10_lod, lod_above_v_supr_max_bool);
  // For each LOD, we can also calculate the log probability of being below-detection conditional on being on reported treatment
  array[L] real log_prob_bd_given_tx;
  for (l in 1:L){
      log_prob_bd_given_tx[l] = 
        (lod_above_v_supr_max_bool[l] == 1) ? 
        log_sum_exp([log_theta_tx, log1m_theta_tx + log_probs_non_supr_given_tx[l,2]]) :
        log_theta_tx;
  }

  // ----------------------- //
  // TRANSMISSION RATE MODEL //
  // ----------------------- //
  // DEFINE TRANSFORMED PARAMETERS:
  // 1. PARAMETER TRANSFORMATIONS
  real use_param1 = est_param1;
  real use_param2 = est_param2;
  real use_param3 = est_param3;
  if (trans_form == 1){
    use_param1 = inv_logit(est_param1);
    use_param2 = inv_logit(est_param2);
    use_param3 = exp(est_param3);
    // if we don't allow saturation
    // then we force use_param3 to equal c
    if (allow_saturation == 0){
      use_param1 = 1;
    }
    if (c > 0){
      use_param3 = c;
    }
  }else if(trans_form == 3){
    use_param1 = exp(est_param1);
    use_param2 = exp(est_param2);
    use_param3 = 1;
  }

  // 2. TRANSFORM IMPUTED VALUES
  // B. sampled no_tx truncated log10_v when no treatment (reported or unreported)
  // sampled no_tx truncated log10_v always has a lower limit of 0 copies/mL
  // estimated on the [0,1] scale and then scaled to [0,LOD]
  // assumed uniform between [0,LOD] on the copies/mL scale
  // then we take the log
  vector[N_donor_v_no_tx_bd] scaled_donor_no_tx_bd_trunc_sampled_v =
    raw_donor_no_tx_bd_trunc_sampled_v .* lod[lod_idx[donor_v_no_tx_bd_idx]];
  vector[N_donor_v_no_tx_bd] donor_no_tx_bd_trunc_sampled_log10_v = 
    log10(scaled_donor_no_tx_bd_trunc_sampled_v);
  
  // C. imputed un_tx log10_v when unreported suppressive treatment is uniform between v_supr_min and v_supr_max
  // imputed un_tx log10_v has lower limit of log10(v_supr_min) and upper limit of log10(v_supr_max)
  vector[N_donor_v_no_tx_bd]  donor_un_tx_bd_imputed_log10_v = 
    log10(
      v_supr_min + 
        raw_donor_un_tx_bd_imputed_v .* (v_supr_max - v_supr_min));
  
  // D. sampled un_tx truncated log10_v when unreported non-suppressive treatment
  // sampled un_tx truncated log10_v has lower limit of log10(v_supr_max) and upper limit of log10(lod)
  // first we scale it and then we log10 to 
  vector[N_donor_v_no_tx_bd] scaled_donor_un_tx_bd_trunc_sampled_v = 
    v_supr_max + 
        raw_donor_un_tx_bd_trunc_sampled_v .* (lod[lod_idx[donor_v_no_tx_bd_idx]] - v_supr_max);
  vector[N_donor_v_no_tx_bd] donor_un_tx_bd_trunc_sampled_log10_v = 
    log10(scaled_donor_un_tx_bd_trunc_sampled_v);
  // normalization constant
  // values are constrained to be between v_supr_max and lod
  // for, for each lod the total density is
  // log_diff_exp(
  //    normal(log10(lod)), mu_un_tx, sigma_un_tx) -
  //    normal(log10(v_supr_max)), mu_un_tx, sigma_un_tx))
  // We use the equivalent lognromal distribution as we want to 
  // put the prior directly on the scaled vector to avoid jacobian adjustment
  vector[L] scaled_donor_un_tx_bd_trunc_sampled_v_norm;
  for (l in lod_above_v_supr_max_idx){
    scaled_donor_un_tx_bd_trunc_sampled_v_norm[l] = 
      log_diff_exp(
        lognormal_lcdf(lod[l] | mu_un_tx*log(10), sigma_un_tx*log(10)),
        lognormal_lcdf(v_supr_max | mu_un_tx*log(10), sigma_un_tx*log(10)));
  }
  
  // E. imputed tx log10_v when reported suppressive treatment is uniform between v_supr_min and v_supr_max
  // imputed tx log10_v has lower limit of log10(v_supr_max) and upper limit of log10(v_supr_max)
  // first we scale it from [0,1] scale to [v_supr_min, v_supr_max] scale
  vector[N_donor_v_tx_bd] scaled_donor_tx_bd_imputed_v = v_supr_min + 
    raw_donor_tx_bd_imputed_v .* (v_supr_max - v_supr_min);
  vector[N_donor_v_tx_bd]  donor_tx_bd_imputed_log10_v = 
    log10(scaled_donor_tx_bd_imputed_v);
  
  // F. sampled tx truncated log10_v when unreported non-suppressive treatment
  // sampled tx truncated log10_v has lower limit of log10(v_supr_max) and upper limit of log10(lod)
  // first we scale it based on LOD
  vector[N_donor_v_tx_bd] scaled_donor_tx_bd_trunc_sampled_v = v_supr_max + 
    raw_donor_tx_bd_trunc_sampled_v .* (lod[lod_idx[donor_v_tx_bd_idx] ]- v_supr_max);
  // then we log it
  vector[N_donor_v_tx_bd] donor_tx_bd_trunc_sampled_log10_v = 
    log10(scaled_donor_tx_bd_trunc_sampled_v);
  // normalization constant
  // values are constrained to be between v_supr_max and lod
  // for, for each lod the total density is
  // log_diff_exp(
  //    normal(log10(lod)), mu_tx, sigma_tx) -
  //    normal(log10(v_supr_max)), mu_tx, sigma_tx))
  // We use the equivalent lognromal distribution as we want to 
  // put the prior directly on the scaled vector to avoid jacobian adjustment
  vector[L] scaled_donor_tx_bd_trunc_sampled_v_norm;
  for (l in lod_above_v_supr_max_idx){
    scaled_donor_tx_bd_trunc_sampled_v_norm[l] = 
      log_diff_exp(
        lognormal_lcdf(lod[l] | mu_tx*log(10), sigma_tx*log(10)),
        lognormal_lcdf(v_supr_max | mu_tx*log(10), sigma_tx*log(10)));
  }

  // 3. PROBABILITY OF NO INTRA-COUPLE TRANSMISSION
  vector[P] log_prob_no_intra_transmission = rep_vector(0, P);
  vector[P] log_prob_no_transmission = rep_vector(0,P);
  vector[P] prob_intra_given_transmission = rep_vector(0,P);
  real cum_prob_intra_given_transmission = 0; // cumulative probability of intra-couple transmission given transmission occurred
  // Above detections
  // if (trans_vl == 1), then use above detection VL as observed
  if (trans_vl == 1){
    for (idx in period_v_ad_idx){
      log_prob_no_intra_transmission[idx[1]] += 
        calc_log_prob_no_intra(
          trans_form,
          use_param1, use_param2, use_param3,
          period_t[idx[1], idx[2]],
          log10_v[period_donor_v_idx[idx[1], idx[2]]]);
    }
  }
  // if (trans_vl == 2), then assign transmission to VL-independent process
  // after self-reported or inferred initiation of self-reported treatment
  if (trans_vl == 2){
    for (idx in period_v_ad_idx){
      if (tx[period_donor_v_idx[idx[1], idx[2]]] == 1){
        log_prob_no_intra_transmission[idx[1]] = -beta_tx*period_t[idx[1], idx[2]];
      }
     else if (tx[period_donor_v_idx[idx[1], idx[2]]] == 0){
        log_prob_no_intra_transmission[idx[1]] = 
          (any_prev_putative_tx[period_donor_v_idx[idx[1], idx[2]]] == 1) ? 
            log_sum_exp([
              log_prob_ever_tx[period_donor_v_idx[idx[1],idx[2]]] + (-beta_tx*period_t[idx[1], idx[2]]),
              log_prob_never_tx[period_donor_v_idx[idx[1],idx[2]]] + 
                calc_log_prob_no_intra(
                trans_form,
                use_param1, use_param2, use_param3,
                period_t[idx[1], idx[2]],
                log10_v[period_donor_v_idx[idx[1], idx[2]]])]) : 
            calc_log_prob_no_intra(
                trans_form,
                use_param1, use_param2, use_param3,
                period_t[idx[1], idx[2]],
                log10_v[period_donor_v_idx[idx[1], idx[2]]]);
      }
    }
  }
  // Below detections with no self-reported treatment
  // define some temporary data objects for each of the below indices
  vector[4] norm_log_prob_bd; // probability of each data generating process
  for (idx in period_v_no_tx_bd_idx){
    // below-detection with no reported ART
    // --> prior to unreported treatment:
      // --> 1. false negative, then impute from intra-individual VL distirbution
      // ----> ~N(mu_i, sigma_i)
      // --> 2. true negative, then impute from intra-individual VL distribution truncated at LOD (no-TX TN)
      // ----> ~N^{LOD}(mu_i, sigma_i)
    // --> after unreported treatment
      // --> 3. suppressed, then impute from tx-induced suppression distribution (TX Supr-TN) between v_supr_min and v_supr_max
      // ----> ~Uniform(v_supr_min, v_supr_max)
      // --> 4. non-suppressed, then impute tx non-suppressed distribution (TX NonSupr-TN) following N(mu_un_tx, sigma_un_tx) truncatd at v_supr_min and v_supr_max
      // ----> ~N_{v_supr_min}^LOD (mu_un_tx, sigma_un_tx)
    // 3 and 4 can only occur after treatment availability
    // normalized probabilities of four potetial data-generating processes:
    // --> 1. no-TX FN
    // when treatment not available: x prob_not_fn
    // when treatment available: x prob_not_fn x prob_not_tx_init
    norm_log_prob_bd[1] =
        period_tx_avail[idx[1], idx[2]] ? 
        log_prob_never_tx[period_donor_v_idx[idx[1],idx[2]]] + log_prob_not_tx_init_fn :
        log_prob_fn;
    // --> 2. no-TX TN
    // when treatment not available: x prob_not_fn
    // when treatment available: x prob_not_fn x prob_not_tx_init
    // + probability of BD given no treatment and no false negative
    norm_log_prob_bd[2] = 
        (period_tx_avail[idx[1], idx[2]] ? 
            log_prob_never_tx[period_donor_v_idx[idx[1],idx[2]]] + 
              log_prob_not_tx_init_not_fn :
            log_prob_not_fn) +
          ind_log_prob_bd_given_no_tx_no_fn[ind_idx[period_donor_v_idx[idx[1], idx[2]]], lod_idx[period_donor_v_idx[idx[1], idx[2]]]];
    if (period_tx_avail[idx[1], idx[2]] == 1){
      // 3. Tx Supr-TN
      // only include ever_tx term if there has been at least one previous putative TX observation
      // else, that process is not considered
      norm_log_prob_bd[3] = 
        log_sum_exp(
          [
            log_prob_never_tx[period_donor_v_idx[idx[1],idx[2]]] + 
              log_prob_tx_init,
            log_prob_ever_tx[period_donor_v_idx[idx[1],idx[2]]] + 
              log_theta_un_tx]
          [1:(1 + period_any_prev_putative_tx[idx[1], idx[2]])]);
      // 4. TX NonSupr-TN
      // only consider if at least one previous putative TX observation
      if (period_any_prev_putative_tx[idx[1], idx[2]] == 1){
        norm_log_prob_bd[4] = 
          log_prob_ever_tx[period_donor_v_idx[idx[1],idx[2]]] + 
          log1m_theta_un_tx + 
          log_probs_non_supr_given_un_tx[lod_idx[period_donor_v_idx[idx[1], idx[2]]], 2];
      }
    }    
    
    // normalize to sum 1
    // always include first two terms
    // only include third term if treatment available
    // only include fourth term if treatment availble and previous treatment use
    norm_log_prob_bd[:n_norm_log_prob_bd[idx[1], idx[2]]] = 
      normalize_log_probs(norm_log_prob_bd[:n_norm_log_prob_bd[idx[1], idx[2]]]);
    
    // calculated weighted sum of potential probabilities of no intra-couple transmission
    // drop term describing suppression on unreported treatment when that is not possible
    if (trans_vl == 1){
      log_prob_no_intra_transmission[idx[1]] += 
        log_sum_exp([
          norm_log_prob_bd[1] + 
            calc_log_prob_no_intra(
              trans_form,
              use_param1, use_param2, use_param3,
              period_t[idx[1], idx[2]],
              donor_no_tx_bd_sampled_log10_v[donor_v_no_tx_bd_rev_idx[period_donor_v_idx[idx[1], idx[2]]]]),
           norm_log_prob_bd[2] +
            calc_log_prob_no_intra(
              trans_form,
              use_param1, use_param2, use_param3,
              period_t[idx[1], idx[2]],
              donor_no_tx_bd_trunc_sampled_log10_v[donor_v_no_tx_bd_rev_idx[period_donor_v_idx[idx[1], idx[2]]]]),
            (n_norm_log_prob_bd[idx[1], idx[2]] >= 3 ? 
              norm_log_prob_bd[3] + 
                calc_log_prob_no_intra(
                  trans_form,
                  use_param1, use_param2, use_param3,
                  period_t[idx[1], idx[2]],
                   donor_un_tx_bd_imputed_log10_v[donor_v_no_tx_bd_rev_idx[period_donor_v_idx[idx[1], idx[2]]]]) : 
              0),
            (n_norm_log_prob_bd[idx[1], idx[2]] >= 4 ? 
              norm_log_prob_bd[4] + 
                calc_log_prob_no_intra(
                  trans_form,
                  use_param1, use_param2, use_param3,
                  period_t[idx[1], idx[2]],
                   donor_un_tx_bd_trunc_sampled_log10_v[donor_v_no_tx_bd_rev_idx[period_donor_v_idx[idx[1], idx[2]]]]) : 
              0)][:n_norm_log_prob_bd[idx[1], idx[2]]]);
      }else if (trans_vl == 2){
        log_prob_no_intra_transmission[idx[1]] += 
          log_sum_exp([
            norm_log_prob_bd[1] + 
              calc_log_prob_no_intra(
                trans_form,
                use_param1, use_param2, use_param3,
                period_t[idx[1], idx[2]],
                donor_no_tx_bd_sampled_log10_v[donor_v_no_tx_bd_rev_idx[period_donor_v_idx[idx[1], idx[2]]]]),
             norm_log_prob_bd[2] +
              calc_log_prob_no_intra(
                trans_form,
                use_param1, use_param2, use_param3,
                period_t[idx[1], idx[2]],
                donor_no_tx_bd_trunc_sampled_log10_v[donor_v_no_tx_bd_rev_idx[period_donor_v_idx[idx[1], idx[2]]]]),
            ((period_tx_avail[idx[1], idx[2]]>0) ? 
              log_sum_exp(norm_log_prob_bd[3:4][1:(1 + (n_norm_log_prob_bd[idx[1], idx[2]]>=4))]) + 
                -beta_tx*period_t[idx[1], idx[2]] : 
              0)][1:(2 + period_tx_avail[idx[1], idx[2]])]);
      }
  }
  // below-detection during self-reported treatment
  for (idx in period_v_tx_bd_idx){
    if (trans_vl == 1){
    // can eitehr arise from suppression
    // or arise from non-suppressed BD
    // 1. Tx Supr-TN
    norm_log_prob_bd[1] = log_theta_tx;
    // 2. Tx NonSupr-TN
    norm_log_prob_bd[2] = log1m_theta_tx + log_probs_non_supr_given_tx[lod_idx[period_donor_v_idx[idx[1], idx[2]]], 2];
    norm_log_prob_bd[1:2] = normalize_log_probs(norm_log_prob_bd[1:2]);  
    log_prob_no_intra_transmission[idx[1]] += 
      log_sum_exp([
        norm_log_prob_bd[1] + 
          calc_log_prob_no_intra(
            trans_form,
            use_param1, use_param2, use_param3,
            period_t[idx[1], idx[2]],
            donor_tx_bd_imputed_log10_v[donor_v_tx_bd_rev_idx[period_donor_v_idx[idx[1], idx[2]]]]),
        norm_log_prob_bd[2] + 
          calc_log_prob_no_intra(
            trans_form,
            use_param1, use_param2, use_param3,
            period_t[idx[1], idx[2]],
             donor_tx_bd_trunc_sampled_log10_v[donor_v_tx_bd_rev_idx[period_donor_v_idx[idx[1], idx[2]]]])]);
    }else if (trans_vl == 2){
      // if trans_vl == 2 then we chalk all transmission up to VL-independent process
      log_prob_no_intra_transmission[idx[1]] += 
        -beta_tx*period_t[idx[1], idx[2]];
    }
  }

  // 4. PROBABILITY OF NO TRANSMISSION (INTRA- OR EXTRA-)
  // incorporate extra-partner transmission rate
  log_prob_no_transmission = log_prob_no_intra_transmission + -1*tot_period_t*beta_extra;
  
  // if requested, calculate cumulative probability of intra-couple transmission for couples with linkage data
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
  mu_un_tx ~ normal(4,2); // mean log10 viral load of non-suppressed unreported treatment
  sigma_un_tx ~ normal(0,1); // variation of non-suppressed unreported treatment
  mu_tx ~ normal(4,2); // mean log10 viral load of non-suppressed treatment-experienced
  sigma_tx ~ normal(0,1); // variation of non-suppressed treatment-experienced
  logit_theta_un_tx ~ normal(0,1); // mixing parameter between suppressed and non-suppressed among unreported treatment
  logit_theta_tx ~ normal(0,1); // mixing parameter between suppressed and non-suppressed among treatment-experienced
  
  logit_fnr ~ normal(0,1); // false negative rate
  logit_txr ~ normal(0,1); // logit probability of unreported treatment initialization 
  
  // ----------------------- //
  // TRANSMISSION RATE MODEL //
  // ----------------------- //
  if (trans_form == 1){
    est_param1 ~ normal(0,5);
    est_param2 ~ normal(0,5);
    est_param3 ~ normal(0,1);
  }else if (trans_form == 2){
    est_param1 ~ normal(0, 0.5);
    est_param2 ~ lognormal(0,1);
    est_param3 ~ normal(4,2);
  }else if (trans_form == 3){
    est_param1 ~ normal(0,5);
    est_param2 ~ normal(0,5);
    est_param3 ~ normal(0,1);
  }
  

  beta_tx ~ normal(0, 0.01);
  if (constrain_prob_intra_given_transmission == 1){
    // Wiki: expected value of the posterior distribution over p, namely Beta(s+1, n−s+1)
    (cum_prob_intra_given_transmission/N_constrained) ~ beta(N_constrained_linked_plus_one, N_constrained_unlinked_plus_one);
  }else{
    beta_extra ~ normal(0, 0.01);
  }

    // prior for imputed BD values
  // A. sampled no_tx log10_v follows N(mu_i, sigma_i)
  donor_no_tx_bd_sampled_log10_v ~ 
    normal(
      mu_i[ind_idx[donor_v_no_tx_bd_idx]],
      sigma_i[ind_idx[donor_v_no_tx_bd_idx]]);

  // B. sampled no_tx truncated log10_v follows N(mu_i, sigma_i) right-truncated at LOD
  // truncation handled in transformation of raw values
  // let y = donor_no_tx_bd_trunc_sampled_log10_v
  // let x = raw_donor_no_tx_bd_trunc_sampled_v
  // let a = LOD
  // raw_donor_no_tx_bd_trunc_sampled_v*LOD = scaled_donor_no_tx_bd_trunc_sampled_v
  // y ~ Normal(mu_i, sigma_i) [this is the prior we want, but trying to avoid dealing with jacobian adjustment]
  // y = log10(a*x)
  // a*x ~ Log10normal(mu_i, sigma_i)
  // a*x ~ Lognormal(mu_i*log(10), sigma_i*log(10))
  scaled_donor_no_tx_bd_trunc_sampled_v ~ 
    lognormal(
      mu_i[ind_idx[donor_v_no_tx_bd_idx]]*log(10),
      sigma_i[ind_idx[donor_v_no_tx_bd_idx]]*log(10));
  // normalization constants because these values are truncated at LOD
  // thus, total log probability between 0 and LOD is
  // Normal_lcdf(lod | mu_i, sigma_i)
  // which is stored in ind_log_prob_bd_given_no_tx_no_fn
  for (idx in donor_v_no_tx_bd_idx){
    target += -1 * ind_log_prob_bd_given_no_tx_no_fn[ind_idx[idx], lod_idx[idx]];
  }

  // C. imputed un_tx log10_v follows Unif(v_supr_min, v_supr_max)
  // raw values are transformed to be defined on 0,1
  raw_donor_un_tx_bd_imputed_v ~ uniform(0, 1);

  // D. sampled un_tx log10_v follows N(mu_untx, sigma_untx) left-truncated at v_supr_max and right-truncated at LOD
  // truncated handled in transformation of raw values
  // log10 of the scaled values are normally distributed
  // thus, scaled values are lognormally distributed with log(10) factor on mean and sigma
  scaled_donor_un_tx_bd_trunc_sampled_v  ~ 
    lognormal(
      mu_un_tx*log(10),
      sigma_un_tx*log(10));
  // normalization constants because these values are truncated
  // between v_supr_max and  lod
  // thus, total log probability between v_supr_max and lod is 
  // log_diff_exp(normal_lcdf(lod | mu_un_tx, sigma_un_tx), normal_lcdf(v_supr_max | mu_un_tx, sigma_un_tx))
  // which is stored in scaled_donor_un_tx_bd_trunc_sampled_v_norm
  // calculate one for each LOD since only depends on LOD and not any individual-level parameters
  target += -1 * scaled_donor_un_tx_bd_trunc_sampled_v_norm[lod_idx[donor_v_no_tx_bd_idx]];


  // E. imputed tx log10_v follows Unif(v_supr_min, v_supr_max)
  // raw values are transformed to be defined on 0,1
  raw_donor_tx_bd_imputed_v ~ uniform(0, 1);

  // F. sampled tx log10_v follows N(mu_tx, sigma_tx) left-truncated at v_supr_max and right-truncated at LOD
  // truncated handled in transformation of raw values
  // log10 of the scales values are normally distributed
  // thus, scaled values are lognormally distributed wiht log(10) factor on mean and sigma
  scaled_donor_tx_bd_trunc_sampled_v  ~ 
   lognormal(
      mu_tx*log(10),
      sigma_tx*log(10));
  // normalization constants because these values are truncated
  // between v_supr_max and  lod
  // thus, total log probability between v_supr_max and lod is 
  // log_diff_exp(normal_lcdf(lod | mu_tx, sigma_tx), normal_lcdf(v_supr_max | mu_tx, sigma_tx))
  // which is stored in scaled_donor_tx_bd_trunc_sampled_v_norm
  // calculate one for each LOD since only depends on LOD and not any individual-level parameters
  target += -1 * scaled_donor_tx_bd_trunc_sampled_v_norm[lod_idx[donor_v_tx_bd_idx]];



  // LIKELIHOOD
  if (sample_from_posterior == 1){
    // ---------------- //
    // VIRAL LOAD MODEL //
    // ---------------- //
    // 1. Pre-self-reported treatment above LOD
    // -> Can arise from
    // --> 1. No false negative, no treatment initiation, no prior treatment, intra-individual VL variation
    //        log_prob_never_tx set to 0 (Prob. = 1) prior to treatment availability
    //        Treatment initiation only possible when tx is available, indicated by tx_avail
    // --> 2. Intra-individual previous unreported treatment distribution
    //        Only possible when there has been at least one previous putative treatment, indicated by any_prev_putative_tx
    //        Need to account for mixing parameter (log1m_theta_un_tx) for probability of non-suppressed given previous treatment
    //        Need to account for the fact that N(un_tx, sigma_tx) is left-truncated at v_supr_max
    //        This is given by first column of log_probs_non_supr_given_un_tx (shared across rows)
    for (idx in no_tx_non_bd_idx){
      target += log_sum_exp(
        [
          ((tx_avail[idx] == 1) ?
              log_prob_never_tx[idx] + log_prob_not_tx_init_not_fn :
              log_prob_not_fn) + 
            normal_lpdf(log10_v[idx] | mu_i[ind_idx[idx]], sigma_i[ind_idx[idx]]),
          ((any_prev_putative_tx[idx]==1) ? 
            log_prob_ever_tx[idx] + 
                log1m_theta_un_tx + 
                normal_lpdf(log10_v[idx] | mu_un_tx, sigma_un_tx) - 
                log_probs_non_supr_given_un_tx[1,1] :
            0)
        ][1:(1 + any_prev_putative_tx[idx])]);
    }

    // 2. Pre-self-reported treatment below LOD
    // -> Can arise from:
    // --> 1. No previous treatment initiation and false negative or treatment initiation at this round
    //        log_prob_never_tx set to 0 (Prob. = 1) prior to treatment availability
    //        Treatment initiation only possible when treatment available
    // --> 2. No previous treatment initiation, no false negative, no treatment initiation at this round, intra-individual VL variation
    //        log_prob_never_tx set to 0 (Prob. = 1) prior to treatment availability
    //        Treatment initiation only possible when treatment available     
    // --> 3. Intra-individual previous unreported treatment VL variation
    //        (the log probability of which occurring is given by log_prob_ever_tx,
    //        can only happen when at least 1 previous putative BD observation)
    //        Can get a below-detection form preivous unreported treatment either through
    //          1) Suppression, which occurs with log probability log_theta_un_tx
    //          2) Non-suppression below-detection,
    //            which occurs with log probability log1m_theta_un_tx + normalized density between v_supr_max and LOD
    //            normalized density between v_supr_max and LOD is second column of log_probs_non_supr_given_un_tx array
    //        Pre-calculate the sum of these two possibilities in log_prob_bd_given_un_tx array
    for (idx in no_tx_bd_idx){
      target += log_sum_exp(
        [
          (tx_avail[idx] == 1) ?
            log_prob_never_tx[idx] + log_prob_tx_init_or_fn :
            log_prob_fn,
          ((tx_avail[idx] == 1) ?
              log_prob_never_tx[idx] + log_prob_not_tx_init_not_fn :
              log_prob_not_fn) + 
            ind_log_prob_bd_given_no_tx_no_fn[ind_idx[idx], lod_idx[idx]],
          ((any_prev_putative_tx[idx]==1) ? 
            log_prob_ever_tx[idx] + 
              log_prob_bd_given_un_tx[lod_idx[idx]] :
            0)
        ][1:(2 + any_prev_putative_tx[idx])]);
    }

    // 3. reported treatment above LOD
    // -> Can arise from
    // --> 1. Intra-individual reported treatment distribution
    // need to adjust for the fact that normal distribution is truncated at log10_v_supr_max
    // normalization constant is first column of log_probs_non_supr_given_tx
    // shared across all rows so just use the first
    for (idx in tx_non_bd_idx){
      target += log1m_theta_tx  + normal_lpdf(log10_v[idx] | mu_tx, sigma_tx) - log_probs_non_supr_given_tx[1,1];
    }


    // 4. reported treatment below LOD
    // -> Can arise from
    // --> 1. Suppression
    // --> 2. Non-suppresion BD (only when LOD > v_supr_max, as indicated by lod_above_v_supr_max_bool)
    // Pre-calculate the sum of these two probabilities in log_prob_bd_given_tx array
    for (idx in tx_bd_idx){
      target += log_prob_bd_given_tx[lod_idx[idx]];
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
  real theta_un_tx = inv_logit(logit_theta_un_tx);
  real theta_tx = inv_logit(logit_theta_tx);
}



