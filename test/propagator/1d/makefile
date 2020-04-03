include ../../../var.mk
include ../../var.mk

TDSE_LIB_DIR = ../../../lib
TDSE_LIB_FLAGS = '-Wl,-rpath,$$ORIGIN/$(TDSE_LIB_DIR)' -L$(TDSE_LIB_DIR) -ltdse

bins = 1d 1d-under-field 1d-with-imagpot

all: $(bins)

1d: 1d.cc libtdse.so libparam.so
	$(CXX) -o $@$(SUFFIX) $< $(OPT_FLAGS) $(WARN_FLAGS) \
		$(TDSE_LIB_FLAGS) $(PARAM_INC_FLAGS) $(PARAM_LIB_FLAGS) \
		$(DEFINE_FLAGS)

1d-under-field:
	$(MAKE) 1d SUFFIX=-under-field DEFINE_FLAGS=-DFIELD

1d-with-imagpot:
	$(MAKE) 1d SUFFIX=-with-imagpot DEFINE_FLAGS=-DIMAGPOT

#1d-under-field: 1d.cc libtdse.so libparam.so
#	$(CXX) -o $@ $< $(OPT_FLAGS) $(WARN_FLAGS) \
#		$(TDSE_LIB_FLAGS) $(PARAM_INC_FLAGS) $(PARAM_LIB_FLAGS) -DFIELD

libtdse.so:
	$(MAKE) -C ../../../src

libparam.so:
	$(MAKE) -C $(PARAMHOME)

clean:
	$(RM) $(bins) 

