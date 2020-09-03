CPPC=g++
CFLAGS= -std=c++11 -Wall -lgsl -llapack

serial: main.cpp
	$(CPPC) main.cpp -o main.exe -lcblas $(CFLAGS)
serial-gsl: main.cpp
	$(CPPC) main.cpp -o main.exe -lgslcblas -lm $(CFLAGS)
parallel: main.cpp
	$(CPPC) main.cpp -o main.exe -lcblas -fopenmp $(CFLAGS)
parallel-gsl: main.cpp
	$(CPPC) main.cpp -o main.exe -lgslcblas -lm -fopenmp $(CFLAGS)
debug: main.cpp
	$(CPPC) -g -O0 main.cpp -o main.exe -lcblas $(CFLAGS)
debug-gsl: main.cpp
	$(CPPC) -g -O0 main.cpp -o main.exe -lgslcblas -lm $(CFLAGS)
clean:
	rm main.exe
