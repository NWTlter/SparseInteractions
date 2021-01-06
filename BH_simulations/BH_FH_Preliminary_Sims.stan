// This script fits a Beverton-Holt generalized competition model using a Finnish (regularized) horseshoe prior (Piironen and Vehtari 2017) 
// 	following the stan implementation demonstrated on https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html

data{
  int<lower = 1> N;
  int<lower = 1> S;
  int<lower = 0> Nt[N];
  int<lower = 0> Ntp1[N];
  matrix[N,S] SpMatrix;
  vector[N] env;
  real tau0; 		// determines the scale of the global shrinkage parameter (tau)
  real slab_scale;	// scale for significant alpha_sp values
  real slab_df;		// effective degrees of freedom for significant alpha_sp values
}

transformed data{
  real slab_scale2 = square(slab_scale);
  real half_slab_df = 0.5*slab_df;
}

parameters{
  vector[3] lambdas_tilde;   // 1: lambda_max, 2: z (env. opt.), 3: sigma (niche breadth)
  vector[2] alphas_tilde;
  vector[S] alpha_hat_ij_tilde;
  vector[S] alpha_hat_eij_tilde;
  vector<lower = 0>[S] local_shrinkage_ij;
  vector<lower = 0>[S] local_shrinkage_eij;
  real<lower = 0> c2_tilde;
  real<lower = 0> tau_tilde;
}

transformed parameters{
  // Calculate the scaled parameters needed for the regularized horeshoe prior here from the normalized (and thus easier to sample)
  // 	counterparts declared in the parameters block
  real c2;
  real tau;
  vector[S] alpha_hat_ij;
  vector[S] local_shrinkage_ij_tilde;
  vector[S] alpha_hat_eij;
  vector[S] local_shrinkage_eij_tilde;
  vector[3] lambdas;
  vector[2] alphas;

  tau = tau0*tau_tilde; 	// tau ~ cauchy(0, tau0)
  c2 = slab_scale2*c2_tilde;	// c2 ~ inv_gamma(half_slab_df, half_slab_df*slab_scale2)

  // This calculation follows equation 2.8 in Piironen and Vehtari 2017
  for(s in 1:S){
    local_shrinkage_ij_tilde[s] = sqrt( c2 * square(local_shrinkage_ij[s]) / (c2 + square(tau) * square(local_shrinkage_ij[s])) );
    alpha_hat_ij[s] = tau * local_shrinkage_ij_tilde[s] * alpha_hat_ij_tilde[s];

    local_shrinkage_eij_tilde[s] = sqrt( c2 * square(local_shrinkage_eij[s]) / (c2 + square(tau) * square(local_shrinkage_eij[s])) );
    alpha_hat_eij[s] = tau * local_shrinkage_eij_tilde[s] * alpha_hat_eij_tilde[s];
  }

  // scale the lambdas and alphas values
  for(i in 1:2){
    alphas[i] = 10 * alphas_tilde[i];
    lambdas[i] = 10 * lambdas_tilde[i];
  }
  lambdas[3] = 10 * lambdas_tilde[3];
}

model{
  // Declare objects necessary for the rest of the model, including: a vector of expected fecundity values (F_hat),
  //     a matrix of the species specific alpha values for each species and plot (interaction_effects), and a matrix
  //     of the the alpha*N values for each species.
  vector[N] Ntp1_hat;
  vector[N] interaction_effects;
  matrix[N,S] alpha_eij;
  vector[N] lambda_ei;

  // set regular priors
  alphas_tilde ~ normal(0,1);
  lambdas_tilde ~ normal(0, 1);

  // set the hierarchical priors for the Finnish horseshoe (regularized horseshoe) (Piironen and Vehtari 2017)
  // Following the stan implementation from https://betanalpha.github.io/assets/case_studies/bayes_sparse_regression.html
  alpha_hat_ij_tilde ~ normal(0,1);
  local_shrinkage_ij ~ cauchy(0,1);

  alpha_hat_eij_tilde ~ normal(0,1);
  local_shrinkage_eij ~ cauchy(0,1);

  tau_tilde ~ cauchy(0,1);
  c2_tilde ~ inv_gamma(half_slab_df, half_slab_df);

  // implement the biological model
  for(i in 1:N){
    lambda_ei[i] = lambdas[1] * exp(-1*((lambdas[2] - env[i])/(2*lambdas[3]))^2);
    for(s in 1:S){
        alpha_eij[i,s] = exp(alphas[1] + alpha_hat_ij[s] + (alphas[2] + alpha_hat_eij[s]) * env[i]);
    }
    interaction_effects[i] = sum(alpha_eij[i,] .* SpMatrix[i,]);
    
    Ntp1_hat[i] = Nt[i] * lambda_ei[i] / (1 + interaction_effects[i]);
    if(Ntp1_hat[i] > 0){
      Ntp1[i] ~ poisson(Ntp1_hat[i]);
    }
  }
}