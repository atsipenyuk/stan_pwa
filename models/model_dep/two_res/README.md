TWO_TOY_RES
===========

This model shows how to generate and fit data for the following process:
'D -> Rc -> abc' -- where D is a charged D meson, a,b,c are pions, and
there are two intermediate resonances R: two fictitious pseudoscalar particles
f_0(1000), f_0(1200), with respective masses 1000, 1200 MeV, and widths of
100 MeV each. The resonances are modeled using Breit-Wigner dynamical factors.

The model function 'f_model = |theta_1 * f_0(1000) + theta_2 * f_0(1200)|^2'
is defined in 'src/model.hpp'.
