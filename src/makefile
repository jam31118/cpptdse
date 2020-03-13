include var.mk

LIB_DIR=../lib
WF_LIB_DIR = ./wf
PROPAGATOR_LIB_DIR = ./propagator

all: install

install: libtdse.so
	if [ ! -d "$$LIB_DIR" ]; then mkdir -p $(LIB_DIR); mv $^ $(LIB_DIR); fi

libtdse.so: tridiag.o array.o libpropagator.a libwf.a
	$(CXX) -shared *.o -o $@ \
		-Wl,--whole-archive $(PROPAGATOR_LIB_DIR)/libpropagator.a $(WF_LIB_DIR)/libwf.a \
		-Wl,--no-whole-archive $(MATRIX_LIB_FLAGS)

tridiag.o: tridiag.cc
	$(CXX) -c -o $@ $< -fPIC $(OPT_FLAGS) $(MATRIX_INC_FLAGS)

array.o: array.cc
	$(CXX) -c -o $@ $< -fPIC $(OPT_FLAGS)

libpropagator.a:
	$(MAKE) -C $(PROPAGATOR_LIB_DIR)

libwf.a:
	$(MAKE) -C $(WF_LIB_DIR)

clean:
	$(RM) *.o *.so $(LIB_DIR)/*.so
	if [ -d "$$LIB_DIR" ]; then rmdir $(LIB_DIR); fi
	$(MAKE) clean -C $(PROPAGATOR_LIB_DIR)
	$(MAKE) clean -C $(WF_LIB_DIR)

