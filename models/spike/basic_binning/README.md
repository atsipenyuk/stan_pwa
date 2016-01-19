This is a basic model to test the model-independent resonance fitting.
First, the data is simulated using the function f_gen containing

 - f0(1000) toy resonance (S-wave);
 - f0(1200) toy resonance (S-wave);
 - rho(770) resonance (P-wave) (reference resonance, amplitude set to 1+0j);

and the incoherent background

 - rho(770)

(that shall be first set to 0).

After the data is generated, it is processed in a slightly different
way than in the model-dependent case. 