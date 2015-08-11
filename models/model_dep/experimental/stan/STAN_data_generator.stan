data {
  vector[num_resonances()] theta[2];
}

parameters {
  // Adjust the boundaries of the Dalitz plot, if needed
  vector<lower=0., upper=3.>[num_variables()] y;
}

model {

  real logH;
  logH <- 0;

  //logH <- logH + log( f_genfit(amplitude_vector(y), theta) );
  logH <- logH + log(f_genfit(theta,theta));
  increment_log_prob(logH);

}






















