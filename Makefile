all:
	make frequency_spectrum
	make teste
	make cross_section

cross_section: cross_section.cpp headers/cross_section.hpp headers/phys_const.hpp headers/point_like_charge.hpp
	g++ -Wall cross_section.cpp -lgsl -lm -o cross_section

frequency_spectrum: frequency_spectrum.cpp headers/point_like_charge.hpp
	g++ -Wall -O frequency_spectrum.cpp -lgsl -lm -o frequency_spectrum

teste: teste.cpp headers/phys_const.hpp headers/point_like_charge.hpp
	g++ -lgsl -lm teste.cpp -o teste

clean:
	rm frequency_spectrum
	rm teste
	rm cross_section
