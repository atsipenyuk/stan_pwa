// Generate the data according to 'f_genfit(y, theta)'
// With other words, fix 'theta' and generate a bunch
// of 'y' with the distribution 'f_genfit'.

// Warn Stan that 'theta' will be fixed.
// Tell Stan about the size of 'theta'.
// In our file, 'theta' is a complex vector
// with num_resonances() entries. 
//  - Stan does not know complex numbers, so we just tell 
// that theta is a 2-dim. array of vectors.
//  - num_resonances should be fixed in model_resonances.hpp
//  - contents of theta are fixed in STAN_data_generator.data.stan
data {
  vector[num_resonances()] theta[2];
}


// Tell Stan that we sample over 'y'.
// 'y' is a real-valued vector (y={m2_12, m2_23,...}).
// Fix the boundaries of 'y' for better sampling.
//  - num_variables() should be fixed in model_resonances.hpp
parameters {
  // You MUST adjust the boundaries to your needs!
  vector<lower=0.07, upper=2.54>[num_variables()] y;
}

// Tell Stan to sample the distribution 'f_genfit',
// fixed in lib/c_lib/model.hpp
model {

  real logH;
  logH <- 0;

  logH <- logH + log( f_genfit(amplitude_vector(y), theta) );
  increment_log_prob(logH);

}






















