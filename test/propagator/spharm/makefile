include ../../var.mk
include ../../../var.mk

TDSE_LIB_DIR = ../../../lib
TDSE_LIB_FLAGS = '-Wl,-rpath,$$ORIGIN/$(TDSE_LIB_DIR)' -L$(TDSE_LIB_DIR) -ltdse

bins = spharm spharm-coulomb

all: $(bins)

spharm: spharm.cc libtdse.so
	$(CXX) -o $@$(SUFFIX) $< $(OPT_FLAGS) $(WARN_FLAGS) \
		$(TDSE_LIB_FLAGS) $(PARAM_INC_FLAGS) $(PARAM_LIB_FLAGS) \
		$(DEFS)

spharm-coulomb:
	$(MAKE) spharm DEFS=-DCOULOMB SUFFIX=-coulomb

libtdse.so:
	$(MAKE) -C ../../../src

clean: 
	$(RM) $(bins)

