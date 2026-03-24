INCLUDE_FLAGS=-I/usr/local/include -L/usr/local/bin
COMMONFLAGS=-Wall -lgsl -lgslcblas -lm -O2 -v

all:
	make frequency_spectrum
	make teste
	make cross_section
	make cross_section_EPA_qag
	make photon_spectrum
	make cross_section_EPA_rapidity

cross_section: cross_section.cpp headers/cross_section.hpp headers/phys_const.hpp headers/point_like_charge.hpp headers/electron_flux.hpp
	g++ ${INCLUDE_FLAGS} cross_section.cpp ${COMMONFLAGS} -o cross_section

cross_section_EPA_qag: cross_section_EPA_qag.cpp headers/cross_section.hpp headers/phys_const.hpp headers/point_like_charge.hpp headers/electron_flux.hpp
	g++ ${INCLUDE_FLAGS} cross_section_EPA_qag.cpp ${COMMONFLAGS} -o cross_section_EPA_qag

frequency_spectrum: frequency_spectrum.cpp headers/point_like_charge.hpp headers/electron_flux.hpp headers/extended_photon_fluxes.hpp
	g++ ${INCLUDE_FLAGS} frequency_spectrum.cpp ${COMMONFLAGS} -o frequency_spectrum

teste: teste.cpp headers/phys_const.hpp headers/point_like_charge.hpp
	g++ ${INCLUDE_FLAGS} teste.cpp ${COMMONFLAGS} -o teste

photon_spectrum: photon_spectrum_fraction.cpp headers/point_fraction_flux.hpp headers/phys_const.hpp
	g++ ${INCLUDE_FLAGS} photon_spectrum_fraction.cpp ${COMMONFLAGS} -o photon_spectrum

cross_section_EPA_rapidity: cross_section_EPA_rapidity.cpp headers/point_like_charge.hpp headers/phys_const.hpp headers/cross_section.hpp
	g++ ${INCLUDE_FLAGS} cross_section_EPA_rapidity.cpp ${COMMONFLAGS} -o cross_section_EPA_rapidity

clean:
	rm frequency_spectrum
	rm teste
	rm cross_section
	rm cross_section_monte_vegas
	rm cross_section_EPA_qag
	rm photon_spectrum
	rm cross_section_EPA_rapidity
