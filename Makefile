all:
	g++ -L/usr/local/lib/ frequency_spectrum.cpp -lgsl -lgslcblas -lm -o freq_exec

