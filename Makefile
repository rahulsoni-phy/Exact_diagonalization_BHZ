# Platform: linux

EXENAME= ED_exe
### ------ Personal PC compilation ------------
CXX = g++
CPPFLAGS = -std=c++11
#CPPFLAGS += -g3
CPPFLAGS += -O3
CPPFLAGS += -Isrc

LDFLAGS = -llapack -lblas
#CPPFLAGS += -I/usr/include/mkl/
CPPFLAGS += #-fopenmp
LDFLAGS += # -lmkl_intel_lp64 -lmkl_core #-lmkl_intel_thread -lpthread -lm -ldl
LDFLAGS += #-fopenmp

STRIP_COMMAND = true

$(EXENAME): clean main.o
			$(CXX) $(CPPFLAGS) -o $(EXENAME)  main.o $(LDFLAGS)
			$(STRIP_COMMAND) $(EXENAME)

main.o: main.cpp
		$(CXX) $(CPPFLAGS) -c main.cpp -o main.o

all: $(EXENAME)

clean:
			rm -f $(EXENAME) *.o
