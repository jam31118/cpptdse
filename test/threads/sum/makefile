include ../../../var.mk

bins = threads-sum

all: $(bins)

threads-sum: threads-sum.cc
	$(CXX) -o $@ $< $(OPT_FLAGS) $(WARN_FLAGS) -lpthread

clean:
	$(RM) $(bins)

