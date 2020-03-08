TDSE_LIB_DIR = ../../../lib

TDSE_LIB_FLAGS = '-Wl,-rpath,$$ORIGIN/$(TDSE_LIB_DIR)' -L$(TDSE_LIB_DIR) -ltdse

all: 1d

1d: 1d.cc libtdse.so
	$(CXX) -o $@ $< $(TDSE_LIB_FLAGS)

libtdse.so:
	$(MAKE) -C ../../../src

clean:
	$(RM) 1d
