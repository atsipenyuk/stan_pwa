src
===

You must add the link to the file 'model.hpp' to this folder. This file
is the core file provided to Stan; all functions describing the PWA
(such as amplitude_vector, f_model, norm, etc.) either  
 - must be defined in 'model.hpp' or  
 - must be included by 'model.hpp'.  


If you are storing your model file in models/MODEL_DIR/src/model.hpp
(as you should), you can create the necessary link here by calling
'relink_model.sh' from the MODEL_DIR directory. Simply call something like  
`./../../../relink_model.sh`


Reminder: if you want to introduce your own functions to Stan, you can  
 - define them in 'model.hpp';
 - define them in 'your_name.hpp'. 
 
In the latter case, you need to modify the 'makefile reload_libraries'
command. After that, you need to make the necessary additions in
'stan_pwa_function_signatures.h' and run 'make reload_libraries'.
