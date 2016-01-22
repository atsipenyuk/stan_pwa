functions{
  real my_f_gen(real x, real sigma) {
    return exp( - x*x / 2 / sigma / sigma );
  }
  real my_norm(real sigma) {
    return sqrt(2.0 * pi()) * sigma;
  }
}
data {
  int<lower=0> N; // Number of events
  real y[N];
}
parameters {
  real<lower=0.,upper=10.> sigma;
}
model {
  real log_f;
  log_f <- 0;
  for (n in 1:N)
    log_f <- log_f + 
             log( my_f_gen(y[n], sigma) / my_norm(sigma) );
  increment_log_prob(log_f);
}









