data {
  int D; // Number of measured events
  int B1; // Number of bins
  // int B2; // Number of bins in 2nd dimension

  // Event coordinates
  vector[num_variables()] y_data[D];

  // Complex PWA amplitudes corresponding to each event
  vector[num_non_S_res()] amplitude_vector_non_S_data[D,2];

  // Dynamical factors for the binned amplitudes
  vector[D] FFZ_y1_data;
  // Only needed if binning is not symmetric in the coordinates
  //vector[D] FFZ_y2_data;

  // Real normalization vector corresponding to the (orthogonal) bins
  vector[B1] I_y1y1;

  // Complex normalization matrices corresponding to the model
  matrix[num_bins_y1(), num_bins_y2()] I_y1y2[2];
  matrix[num_bins_y1(), num_non_S_res()] I_y1res[2];
  //matrix[num_bins_y2(), num_non_S_res()] I_y2res[2];
  matrix[num_non_S_res(), num_non_S_res()] I_resres[2];

  // Real background PWA amplitudes squared corresponding to each event
  vector[num_background()] background_vector_data[D];
  // Real normalization background vector corresponding to the model
  vector[num_background()] I_background;

}

parameters {
  // Parameters that will be fitted
  // Total: 2 * B1 + 0 (rho is fixed) + 2 (f2_1270) + 2 (background)
  vector<lower=0., upper=5.>[B1] theta_bin_y1_m;
  vector<lower=-pi(), upper=pi()>[B1] theta_bin_y1_ph;

  // Only needed if the binning is not symmetric in the coordinates
  //vector<lower=0., upper=5.>[B2] theta_bin_y2_m;
  //vector<lower=-pi(), upper=pi()>[B2] theta_bin_y2_ph;

  // f2_1270
  real<lower=0., upper=5.> theta_f2_1270_m;
  real<lower=-pi(), upper=pi()> theta_f2_1270_ph;

  //real<lower=0., upper=5.> theta_background_flat_abs2;
  real<lower=0., upper=5.> theta_background_rho_770_abs2;
}


transformed parameters {
  // Parameters: some fixed (reference parameters), some free
  // (these will be fitted).

  // Binned
  vector<lower=-5., upper=5.>[num_bins_y1()] theta_bin_y1[2];
  //vector<lower=-5., upper=5.>[num_bins_y2()] theta_bin_y2[2];

  // Model-dependent
  vector<lower=-5., upper=5.>[num_non_S_res()] theta_res[2];

  // Background
  vector<lower=0., upper=5.>[num_background()] theta_background_abs2;

  // Binned content
  for (b_ in 1:B1) {
    theta_bin_y1[1,b_] <- theta_bin_y1_m[b_] * cos(theta_bin_y1_ph[b_]);
    theta_bin_y1[2,b_] <- theta_bin_y1_m[b_] * sin(theta_bin_y1_ph[b_]);
  }
  //theta_bin_y2[1] <- theta_bin_y2_m * cos(theta_bin_y2_ph);
  //theta_bin_y2[2] <- theta_bin_y2_m * sin(theta_bin_y2_ph);

  // Model-dependent content
  theta_res[1,1] <- 1.0; // The first model-dependent amplitude is rho (fixed)
  theta_res[2,1] <- 0.0;
  theta_res[1,2] <- theta_f2_1270_m * cos(theta_f2_1270_ph);
  theta_res[2,2] <- theta_f2_1270_m * sin(theta_f2_1270_ph);

  // Background content
  //theta_background_abs2[1] <- theta_background_flat_abs2;
  theta_background_abs2[1] <- theta_background_rho_770_abs2;
}


model {
  real logH; // log-likelihood

  // For each event, store its bin coordinates [#bin_y1, #bin_y2]
  int b[2];

  // Bin parameter placeholder (complex number): for each event,
  // holds the parameter of the bin containing this event.
  real theta_current_y1_bin[2]; 
  real theta_current_y2_bin[2]; 

  // Reset log-likelihood
  logH <- 0;

  // Sum over all events
  for (d in 1:D)
  {
    // Select the bin corresponding to the current event
    b <- bin(y_data[d]);

    // Select the bin amplitude corresponding to the current event
    theta_current_y1_bin[1] <- theta_bin_y1[1,b[1]];
    theta_current_y1_bin[2] <- theta_bin_y1[2,b[1]];
    // Remember that theta_bin_y2 = theta_bin_y1 in the symmetric case
    theta_current_y2_bin[1] <- theta_bin_y1[1,b[2]];
    theta_current_y2_bin[2] <- theta_bin_y1[2,b[2]];

    logH <- logH + 
            log(f_fit(FFZ_y1_data[b[1]], theta_current_y1_bin,
	              FFZ_y1_data[b[2]], theta_current_y2_bin,
		      amplitude_vector_non_S_data[d], theta_res, 
		      background_vector_data[d], theta_background_abs2) / 
                norm_bin(theta_bin_y1,
                         theta_bin_y1,
                         theta_res,
                         I_y1y1, I_y1y1,
                         I_y1y2, I_y1res, I_y1res, I_resres,
                         theta_background_abs2, I_background));
  }

  increment_log_prob(logH);
}
