#ifndef STAN_PWA__SRC__STRUCTURES__STRUCT_PARTICLES_HPP
#define STAN_PWA__SRC__STRUCTURES__STRUCT_PARTICLES_HPP

namespace stan_pwa {
  struct particle 
  {
    const double m; // Mass in GeV
    const double m2; // Squared mass
    const double r; // Radius in GeV**-1
    const double r2; // Squared radius
    const int J; // Spin

    particle(double _m, double _r, int _J) : 
      m(_m), m2(_m*_m), r(_r), r2(_r*_r), J(_J) {};
  };
}
#endif
