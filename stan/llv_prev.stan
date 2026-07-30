data {
  int<lower=1> N; // number of observations
  int<lower=1> N_ind;   // number of individuals
  int<lower=1> N_strata;   // number of individuals
  array[N] int y; // whether the outcome is true
  array[N] int x1; // (individual index, 
  array[N] int x2; //strata index)
}

parameters {
  real alpha;
  sum_to_zero_vector[N_strata] beta;
  vector[N_ind] alpha_i;
  real<lower=0> alpha_i_sd;
}

transformed parameters {
  vector[N_ind] alpha_i_effect = alpha_i_sd * alpha_i;
}

model {
  vector[N] logit_mu = alpha + beta[x2] + alpha_i_effect[x1];
  alpha ~ normal(0,2);
  beta ~ normal(0,2);
  alpha_i ~ normal(0,1); // non-centered so must be standard nomral
  alpha_i_sd ~ normal(0,1);
  y ~ bernoulli_logit(logit_mu);
}

generated quantities {
  //vector[N_ind] epsilon = normal_rng(0, alpha_i_sd);
  //vector[N_ind] sim_mu_i = inv_logit(alpha+beta[x2] + normal_rng(0, alpha_i_sd));
  vector[N_strata] strata_prev=  rep_vector(0, N_strata);
  //vector[N_strata] N_per_strata = rep_vector(0, N_strata);
  int N_rep_per_strata = 100;
  {
    for (k in 1:N_strata){
      for (i in 1:N_rep_per_strata){
        strata_prev[k] += inv_logit(alpha + beta[k] + normal_rng(0, alpha_i_sd));
      }
    }
    //vector[N] mu = alpha + alpha_i_effect[x1] + beta[x2];
    
  }
  strata_prev =  strata_prev / N_rep_per_strata;
}






