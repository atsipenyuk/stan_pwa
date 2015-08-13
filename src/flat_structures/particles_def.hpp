#ifndef STAN_PWA__SRC__FLAT_STRUCTURES__PARTICLES_DEF_HPP
#define STAN_PWA__SRC__FLAT_STRUCTURES__PARTICLES_DEF_HPP

namespace stan_pwa {

  struct Particle  {   
    Particle(double _m, double _r, int _J) : 
      m(_m), m2(_m*_m), r(_r), r2(_r*_r), J(_J) {};

    const double m; // Mass in GeV
    const double m2; // Squared mass
    const double r; // Radius in GeV**-1
    const double r2; // Squared radius
    const int J; // Spin

  };


  struct Particle_w : public Particle  {   
    Particle_w(const Particle P, const double _W) : 
      Particle(P), W(_W) {};

    const double W; // Width of the resonance

  };


  struct ExternalParticles4
  {
    ExternalParticles4(const Particle _P,
		       const Particle _a,
		       const Particle _b,
		       const Particle _c) :
      P(_P), a(_a), b(_b), c(_c) {};

    const Particle P;
    const Particle a;
    const Particle b;
    const Particle c;
  };
}
#endif
