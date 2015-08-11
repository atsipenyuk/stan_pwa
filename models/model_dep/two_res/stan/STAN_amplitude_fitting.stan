data {
  // Number of measured events
  int D;
  // Complex PWA amplitudes corresponding to each event
  vector[num_resonances()] amplitude_vector_data[D,2];
  // Complex normalization matrix corresponding to the model
  matrix[num_resonances(), num_resonances()] I[2];
}


parameters {
  // Parameters that will be fitted
  // Total: 2
  real<lower=0., upper=5.> theta_f0_1370_m;
  real<lower=-pi(), upper=pi()> theta_f0_1370_ph;
}


transformed parameters {
  // Parameters: some fixed (reference parameters), 
  // some free (these will be fitted)
  vector<lower=-5., upper=5.>[num_resonances()] theta[2];

  // First index denotes real/complex part, 
  // second index denotes resonance number
  theta[1,1] <- 1.0; // rho_770 is the reference parameter
  theta[2,1] <- 0.0;
  theta[1,2] <- theta_f0_1370_m * cos(theta_f0_1370_ph);
  theta[2,2] <- theta_f0_1370_m * sin(theta_f0_1370_ph);
}


model {
  real logH;
  logH <- 0;
  // Sum over all events
  for (d in 1:D)
    logH <- logH + log( f_genfit(amplitude_vector_data[d], theta) / norm(theta, I) );
    increment_log_prob(logH);
}






