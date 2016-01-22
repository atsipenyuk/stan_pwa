functions{
  real my_f_gen(real x, real sigma) {
    return exp( - x*x / 2 / sigma / sigma );
  }
}
data {
  real<lower=0> sigma;
}
parameters {
  real y;
}
model {
  real log_f;
  log_f <- 0;
  log_f <- log_f + log( my_f_gen(y, sigma) );
  increment_log_prob(log_f);
}
