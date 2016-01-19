// Generate the phase space events with coordinates y.1, y.2, etc.
// The phase_space function must be included by src/model.hpp

parameters {
  vector<lower=0., upper=3.>[num_variables()] y;
}

model {

  real logH;
  logH <- 0;

  logH <- logH + log( phase_space_gen(y) );
  increment_log_prob(logH);

}






















