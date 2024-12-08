all:
	make frequency_spectrum
	make teste
	make cross_section
	make cross_section_monte_vegas
	make cross_section_mass
	make cross_section_EPA_qag

cross_section: cross_section.cpp headers/cross_section.hpp headers/phys_const.hpp headers/point_like_charge.hpp headers/electron_flux.hpp
	g++ -Wall cross_section.cpp -lgsl -lm -o cross_section

cross_section_monte_vegas: cross_section_MVEGAS.cpp headers/cross_section_monte_carlo.hpp headers/phys_const.hpp headers/point_like_charge.hpp 
	g++ -Wall cross_section_MVEGAS.cpp -lgsl -o cross_section_monte_vegas

cross_section_mass: cross_section_mass.cpp headers/cross_section.hpp headers/phys_const.hpp headers/point_like_charge.hpp headers/electron_flux.hpp
	g++ -Wall cross_section_mass.cpp -lgsl -lm -o cross_section_mass

cross_section_EPA_qag: cross_section_EPA_qag.cpp headers/cross_section.hpp headers/phys_const.hpp headers/point_like_charge.hpp headers/electron_flux.hpp headers/point_fraction_flux.hpp
	g++ -Wall cross_section_EPA_qag.cpp -lgsl -lm -o cross_section_EPA_qag

frequency_spectrum: frequency_spectrum.cpp headers/point_like_charge.hpp headers/electron_flux.hpp
	g++ -Wall -O frequency_spectrum.cpp -lgsl -lm -o frequency_spectrum

teste: teste.cpp headers/phys_const.hpp headers/point_like_charge.hpp
	g++ -lgsl -lm teste.cpp -o teste

clean:
	rm frequency_spectrum
	rm teste
	rm cross_section
