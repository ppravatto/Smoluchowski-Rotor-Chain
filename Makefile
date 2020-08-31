CPPC=g++
CFLAGS= -std=c++11 -Wall -lgsl -llapack

serial: main.cpp
	$(CPPC) main.cpp -o main.x -lcblas $(CFLAGS)
serial-gsl: main.cpp
	$(CPPC) main.cpp -o main.x -lgslcblas -lm $(CFLAGS)
parallel: main.cpp
	$(CPPC) main.cpp -o main.x -lcblas -fopenmp $(CFLAGS)
parallel-gsl: main.cpp
	$(CPPC) main.cpp -o main.x -lgslcblas -lm -fopenmp $(CFLAGS)
debug: main.cpp
	$(CPPC) -g -O0 main.cpp -o main.x -lcblas $(CFLAGS)
debug-gsl: main.cpp
	$(CPPC) -g -O0 main.cpp -o main.x -lgslcblas -lm $(CFLAGS)
clean:
	rm main.x
