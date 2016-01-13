# Makefile for stan_pwa
# CAVEAT. Heavily dependent on CmdStan version.
# Current support: CmdStan 2.9.0.
# Remark. Does NOT compile c++ code whatsoever. 
# Just moves files around to the right places.

##
# Install the package.
##
install:
	# Tell the STAN makefile to link cmdstan folder 
	# when building C++ files (this allows us to use
	# '#include <stan_pwa/..>' statements in C++ files)
	# Backup original makefile
	$(shell [ ! -f install_files/makefile_vanilla_stan ] && cp ../makefile install_files/makefile_vanilla_stan)
	# Replace stan makefile by stan_pwa modified version
	$(shell cp install_files/makefile_modified ../makefile)
	# Add functions with 8 arguments to STAN ast.hpp and ast_def.cpp
	# Backup original ast.hpp
	$(shell [ ! -f install_files/ast_vanilla_stan.hpp] mv ../stan_2.9.0/src/stan/lang/ast.hpp install_files/ast_vanilla_stan.hpp)
	# Backup original ast_def.cpp
	$(shell [ ! -f install_files/ast_def_vanilla_stan.cpp] mv ../stan_2.9.0/src/stan/lang/ast_def.cpp install_files/ast_def_vanilla_stan.cpp)
	# Replace ast.hpp and ast_def.hpp
	$(shell cp install_files/ast_modified.hpp ../stan_2.9.0/src/stan/lang/ast.hpp && \
	cp install_files/ast_def_modified.cpp ../stan_2.9.0/src/stan/lang/ast_def.cpp)
	# Create a backup of STAN function_signatures.h file
	$(shell [ ! -f install_files/function_signatures_vanilla_stan.h ] mv ../stan_2.9.0/src/stan/lang/function_signatures.h install_files/function_signatures_vanilla_stan.h)
	$(shell cp install_files/function_signatures_modified.h ../stan/src/stan/lang/function_signatures.h)
	# Expose PWA math functions to STAN parser
	make reload_pwa_library

##
# Revert the changes made to STAN files
##
uninstall:
	$(shell	mv install_files/ast_vanilla_stan.hpp  ../stan/src/stan/lang/ast.hpp && \
	mv install_files/ast_def_vanilla_stan.cpp  ../stan/src/stan/lang/ast_def.cpp && \
	mv install_files/makefile_vanilla_stan ../makefile && \
	mv install_files/function_signatures_vanilla_stan.h ../stan/src/stan/lang/function_signatures.h)
	$(shell	mv install_files/math_vanilla_stan.hpp ../stan/lib/stan_math_2.7.0/stan/math.hpp)
	$(shell	mv install_files/makefile_vanilla_stan ../makefile)


##
# Reload stan_pwa libraries into STAN
##
reload_pwa_library:
	# Make the necessary changes in 'function_signatures.h'
	$(shell [ -f ! install_files/function_signatures_vanilla_stan.hpp ] mv ../stan_2.9.0/src/stan/lang/function_signatures.h install/function_signatures_vanilla_stan.h)
	$(shell cp install/function_signatures_modified.h ../stan_2.9.0/src/stan/lang/function_signatures.h)
	# Tell Stan where the files of stan_pwa are stored
	$(shell [ -f ! install_files/math_vanilla_stan.hpp ] mv ../stan_2.9.0/lib/stan_math_2.9.0/stan/math.hpp install_files/math_vanilla_stan.hpp)
	$(shell cp install_files/math_modified.hpp ../stan_2.9.0/lib/stan_math_2.9.0/stan/math.hpp)
	# STAN binaries must be rebuild
	cd ..;          \
	make clean-all; \
	cd stan_pwa







