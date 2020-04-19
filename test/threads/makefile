include ../../var.mk

bins = threads

all: $(bins)

threads: threads.cc
	$(CXX) -o $@ $< $(OPT_FLAGS) $(WARN_FLAGS) -lpthread

clean:
	$(RM) $(bins)

