#ifndef STAN_PWA__SRC__STRUCTURES_HPP
#define STAN_PWA__SRC__STRUCTURES_HPP

/*
 * structures.hpp
 *
 * Some elements of PWA could be easier stored in structures,
 * not in functions. Basically, your PWA works as follows:
 * you have a bunch of parameters and variables that are
 * passed from STAN, but you also have some constants known
 * from the previous analyses or PDG. Depending on your model,
 * it is convenient to 'glue' all the form-factors together
 * (declaring some new structure, task #1) and insert all the
 * constants (initializing this new structure, task #2).
 * 
 * The files 'struct_*' do the task #1: they contain the information
 * about the structures - for example, about the resonances,
 * which are glued from Blatt-Weisskopf, Breit-Wigner, or other 
 * functions. 
 * 
 * The other files instantiate these structures. For example, 
 * 'particles.hpp' initializes pion and kaon particles with
 * members that describe their mass, radius, etc.
 */
#include <stan_pwa/src/structures/particles_def.hpp>
#include <stan_pwa/src/structures/particles.hpp>
#include <stan_pwa/src/structures/resonances_def.hpp>
#include <stan_pwa/src/structures/three_body_resonances.hpp>
//#include <stan_pwa/src/structures/four_body_resonances.hpp>

#endif
