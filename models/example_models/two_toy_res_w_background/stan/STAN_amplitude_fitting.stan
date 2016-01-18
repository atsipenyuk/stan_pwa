data {
  // Number of measured events
  int D;
  // Complex PWA amplitudes corresponding to each event
  vector[num_resonances()] amplitude_vector_data[D,2];
  // Complex normalization matrix corresponding to the model
  matrix[num_resonances(), num_resonances()] I[2];

  // Real background PWA amplitudes squared corresponding to each event
  vector[num_background()] A_v_background_abs2_data[D];
  // Real normalization background vector corresponding to the model
  vector[num_background()] I_background;
}

parameters {
  // Parameters that will be fitted
  // Total: 2 + 2 = 4 (coherent + background)
  real<lower=-5., upper=5.> theta_re_f0_1200;
  real<lower=-5., upper=5.> theta_im_f0_1200;

  real<lower=0.> theta_background_flat_abs2;
  real<lower=0., upper=5.> theta_background_rho_770_abs2;
}

transformed parameters {
  // Parameters: some fixed (reference parameters), some free
  // (these will be fitted).
  vector<lower=-5., upper=5.>[num_resonances()] theta[2]; 
  vector<lower=0.>[num_resonances()] theta_background_abs2;

  // First index denotes real/complex part, 
  // second index denotes resonance number
  theta[1,1] <- 1.0; // f0_1000 is the reference parameter
  theta[2,1] <- 0.0;
  theta[1,2] <- theta_re_f0_1200;
  theta[2,2] <- theta_im_f0_1200;

  theta_background_abs2[1] <- theta_background_flat_abs2;
  theta_background_abs2[2] <- theta_background_rho_770_abs2;
}

model {
  real logH;
  logH <- 0;

  // Sum over all events
  for (d in 1:D)
    logH <- logH + 
    log(f_genfit(amplitude_vector_data[d], theta, A_v_background_abs2_data[d], 
				      theta_background_abs2) 
	/ norm(theta, I, theta_background_abs2, I_background) );

  increment_log_prob(logH);
}
