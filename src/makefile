include ../var.mk

LIB_DIR=../lib

OPTIMIZE_FLAGS = -O3

all: install

install: $(LIB_DIR)/libcbohm.so

$(LIB_DIR)/libcbohm.so: velocity.o log.o lapack.o lm.o hydrogen.o ode.o
	mkdir -p $(LIB_DIR)
	$(CXX) -shared -fPIC $^ -o $@ $(LAPACK_LIB_FLAGS) $(GSL_LIB_FLAGS) $(QPROP_LIB_FLAGS)

debug:
	$(MAKE) DEBUG_FLAGS=-DDEBUG

velocity.o: velocity.cc
	$(CXX) -fPIC -c $< -o $@ $(OPTIMIZE_FLAGS) $(DEBUG_FLAGS) $(GSL_INC_FLAGS)

log.o: log.cc
	$(CXX) -fPIC -c $< -o $@ $(OPTIMIZE_FLAGS) $(DEBUG_FLAGS)

lapack.o: lapack.cc
	$(CXX) -fPIC -c $< -o $@ $(OPTIMIZE_FLAGS) $(DEBUG_FLAGS)

lm.o: lm.cc
	$(CXX) -fPIC -c $< -o $@ $(OPTIMIZE_FLAGS) $(DEBUG_FLAGS)

hydrogen.o: hydrogen.cc
	$(CXX) -fPIC -c $< -o $@ $(OPTIMIZE_FLAGS) $(DEBUG_FLAGS) $(QPROP_INC_FLAGS)

ode.o: ode.cc
	$(CXX) -fPIC -c $< -o $@ $(OPTIMIZE_FLAGS) $(DEBUG_FLAGS) $(GSL_INC_FLAGS)

clean:
	$(RM) -f *.a *.o *.so $(LIB_DIR)/*.so $(LIB_DIR)/*.a

