TDSE-LDFLAGS = '-Wl,-rpath,$$ORIGIN/../../lib' -L../../lib -ltdse

all: tridiag-diag

tridiag-diag: tridiag-diag.cc libtdse.so
	$(CXX) -o $@ $< $(TDSE-LDFLAGS)

libtdse.so: ../../lib/libtdse.so
	$(MAKE) -C ../../src

clean:
	$(RM) tridiag-diag

