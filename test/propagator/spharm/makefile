include ../../var.mk
include ../../../var.mk

TDSE_LIB_DIR = ../../../lib
TDSE_LIB_FLAGS = '-Wl,-rpath,$$ORIGIN/$(TDSE_LIB_DIR)' -L$(TDSE_LIB_DIR) -ltdse

bins = spharm

all: $(bins)

spharm: spharm.cc
	$(CXX) -o $@ $< $(OPT_FLAGS) $(WARN_FLAGS) \
		$(TDSE_LIB_FLAGS) $(PARAM_INC_FLAGS) $(PARAM_LIB_FLAGS)

clean: 
	$(RM) $(bins)

