include ../var.mk

LIB_DIR=../lib/

OPTIMIZE_FLAGS = -O3
DEBUG_FLAGS = #-DDEBUG

all: install

install: libcbohm.so
	mkdir -p $(LIB_DIR) && mv $^ $(LIB_DIR)

libcbohm.so: velocity.o log.o lapack.o
	$(CXX) -shared -fPIC $^ -o $@ $(LAPACK_LIB_FLAGS)

velocity.o: velocity.cc
	$(CXX) -fPIC -c $< -o $@ $(OPTIMIZE_FLAGS) $(DEBUG_FLAGS)

log.o: log.cc
	$(CXX) -fPIC -c $< -o $@ $(OPTIMIZE_FLAGS) $(DEBUG_FLAGS)

lapack.o: lapack.cc
	$(CXX) -fPIC -c $< -o $@ $(OPTIMIZE_FLAGS) $(DEBUG_FLAGS)

clean:
	$(RM) -f *.a *.o *.so $(LIB_DIR)/*.so $(LIB_DIR)/*.a

