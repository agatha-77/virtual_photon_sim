freq:
	g++ -L/usr/local/lib/ frequency_spectrum.cpp -lgsl -lgslcblas -lm -o freq_exec.o
teste:
	g++ teste.cpp -o exec_teste.o
all:
	make freq 
	make teste

