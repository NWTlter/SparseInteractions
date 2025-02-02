data{
  int<lower = 1> N;
  int<lower = 1> S;
  int Nt[N];
  int Ntp1[N];
  matrix[N,S] SpMatrix;
  vector[N] env;
  int<lower = 0> Intra[S];
  int Inclusion_ij[S];
  int Inclusion_eij[S];
}

parameters{
  vector[3] lambdas;    // 1: intercept, 2: slope
  vector[3] alpha_generic_tilde;
  vector[3] alpha_intra_tilde;
  vector[S] alpha_hat_ij;
  vector[S] alpha_hat_eij;
}

transformed parameters{
  vector[3] alpha_generic;
  vector[3] alpha_intra;

  // scale the lambdas and alphas values
  alpha_generic[1] = 3 * alpha_generic_tilde[1] - 6;
  alpha_generic[2] = alpha_generic_tilde[2] * 0.5;
  alpha_generic[3] = alpha_generic_tilde[3] * 0.5;
  alpha_intra[1] = 3 * alpha_intra_tilde[1] - 6;
  alpha_intra[2] = alpha_intra_tilde[2] * 0.5;
  alpha_intra[3] = alpha_intra_tilde[3] * 0.5;
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
  alpha_generic_tilde ~ normal(0,1);
  alpha_intra_tilde ~ normal(0,1);
  lambdas ~ normal(0,1);
  alpha_hat_ij ~ normal(0,1);
  alpha_hat_eij ~ normal(0,1);

  // implement the biological model
  for(i in 1:N){
    lambda_ei[i] = exp(lambdas[1] + lambdas[2]*env[i] + lambdas[3]*env[i]^2);
    for(s in 1:S){
      alpha_eij[i,s] = exp((1-Intra[s]) * alpha_generic[1] + Intra[s] * alpha_intra[1] + Inclusion_ij[s] * alpha_hat_ij[s] + 
      ((1-Intra[s]) * alpha_generic[2] + Inclusion_eij[s] * alpha_hat_eij[s] + Intra[s] * alpha_intra[2]) * env[i] + 
      ((1-Intra[s]) * alpha_generic[3] + Inclusion_eij[s] * alpha_hat_eij[s] + Intra[s] * alpha_intra[3]) * env[i]^2);
    }
    interaction_effects[i] = sum(alpha_eij[i,] .* SpMatrix[i,]);
    Ntp1_hat[i] = Nt[i] * lambda_ei[i] / (1 + interaction_effects[i]);
    if(Ntp1_hat[i] > 0){
      Ntp1[i] ~ poisson(Ntp1_hat[i]);
    }
  }
}

