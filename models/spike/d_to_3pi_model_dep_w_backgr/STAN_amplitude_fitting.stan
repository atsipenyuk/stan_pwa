data {
  // Number of measured events
  int D;

  // Complex PWA amplitudes corresponding to each event
  vector[num_resonances()] A_cv_data[D,2];
  // Complex normalization matrix corresponding to the model
  matrix[num_resonances(), num_resonances()] I[2];

  // Real background PWA amplitudes squared corresponding to each event
  vector[num_background()] A_v_background_abs2_data[D];
  // Real normalization background vector corresponding to the model
  vector[num_background()] I_background;
}

parameters {
  // Parameters that will be fitted
  // Total: 12 + 2 (background)
  real<lower=-5., upper=5.> theta_re_flat;
  real<lower=-5., upper=5.> theta_im_flat;
  real<lower=-5., upper=5.> theta_re_f0_980;
  real<lower=-5., upper=5.> theta_im_f0_980;
  real<lower=-5., upper=5.> theta_re_f0_600;
  real<lower=-5., upper=5.> theta_im_f0_600;
  real<lower=-5., upper=5.> theta_re_f0_1370;
  real<lower=-5., upper=5.> theta_im_f0_1370;
  real<lower=-5., upper=5.> theta_re_f0_1500;
  real<lower=-5., upper=5.> theta_im_f0_1500;
  real<lower=-5., upper=5.> theta_re_f2_1270;
  real<lower=-5., upper=5.> theta_im_f2_1270;

  real<lower=0., upper=5.> theta_background_flat_abs2;
  real<lower=0., upper=5.> theta_background_rho_770_abs2;

  // Keep coherent and total integral values for bookkeeping purposes
  real<lower=0.> I_bookkeeping_coherent;
  real<lower=0.> I_bookkeeping_total;
}


transformed parameters {
  // Parameters: some fixed (reference parameters), some free
  // (these will be fitted).
  vector<lower=-5., upper=5.>[num_resonances()] theta[2]; 
  vector<lower=0., upper=5.>[num_background()] theta_background_abs2;
  vector<lower=0.>[2] I_bookkeeping;

  theta[1,1] <- theta_re_flat;
  theta[2,1] <- theta_im_flat;
  theta[1,2] <- theta_re_f0_980;
  theta[2,2] <- theta_im_f0_980;
  theta[1,3] <- theta_re_f0_600;
  theta[2,3] <- theta_im_f0_600;
  theta[1,4] <- theta_re_f0_1370;
  theta[2,4] <- theta_im_f0_1370;
  theta[1,5] <- theta_re_f0_1500;
  theta[2,5] <- theta_im_f0_1500;
  theta[1,6] <- 1.0;
  theta[2,6] <- 0.0;
  theta[1,7] <- theta_re_f2_1270;
  theta[2,7] <- theta_im_f2_1270;

  theta_background_abs2[1] <- theta_background_flat_abs2;
  theta_background_abs2[2] <- theta_background_rho_770_abs2;

  I_bookkeeping[1] <- I_bookkeeping_coherent;
  I_bookkeeping[2] <- I_bookkeeping_total;
}

model {
  real logH;
  logH <- 0;
  // Sum over all events
  for (d in 1:D)
    I_bookkeeping <- Norm(theta, I, theta_background_abs2, I_background) );
    logH <- logH + log( f_model(A_cv_data[d], theta,
				A_v_background_abs2_data[d], theta_background_abs2) 
			/ I_bookkeeping[2];
  increment_log_prob(logH);
}
