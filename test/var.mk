PARAM_INC_FLAGS = -I$(PARAMHOME)/include
PARAM_LIB_DIR = $(PARAMHOME)/lib
PARAM_LIB_FLAGS = -Wl,-rpath,$(PARAM_LIB_DIR) -L$(PARAM_LIB_DIR) -lparam

