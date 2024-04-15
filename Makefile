all:
	make clean
	make freq_exec 
	make exec_teste
	make cross_exec

cross_exec:
	g++ -Wall -O -L/usr/local/lib cross_section.cpp -lgsl -lgslcblas -lm -o cross_exec

freq_exec: 
	g++ -Wall -O -L/usr/local/lib/ frequency_spectrum.cpp -lgsl -lgslcblas -lm -o freq_exec

exec_teste:
	g++ teste.cpp -o exec_teste

clean:
	rm freq_exec
	rm exec_teste
	rm cross_exec
