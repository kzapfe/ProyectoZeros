libs= -larmadillo  -lgsl
incluir= -I /usr/include/gsl

4Bolla03.x: simplectic01.hpp intervalitos01.hpp \
ComparaWeylenradio01.hpp RutinasExternas02.hpp PuntosCuerdasEsfera03.hpp \
Cuerda3esfera3.cpp 
	g++ $(libs) Cuerda3esfera3.cpp -o 4Bolla03.x

PNelson01.x: ProbarPoblarNelson01.cpp RutinasNelson01.hpp
	g++ $(libs) ProbarPoblarNelson01.cpp -o PNelson01.x


NelsonS01.x: NelsonSimple01.cpp intervalitos01.hpp RutinasNelson01.hpp
	g++ $(libs) NelsonSimple01.cpp -o NelsonS01.x

NelsonPlan01.x: NelsonPlano01.cpp intervalito02.hpp RutinasNelson01.hpp
	g++ $(libs) NelsonPlano01.cpp -o NelsonPlan01.x

NelsonPlan02.x: NelsonPlano02.cpp intervalito02.hpp RutinasNelson02.hpp \
ParametrosGlobales.hpp
	g++ $(libs) NelsonPlano02.cpp -o NelsonPlan02.x

Zeros01.x: NelZeros01.cpp intervalito02.hpp RutinasNelson02.hpp \
ParametrosGlobales.hpp
	g++ $(libs) NelZeros01.cpp -o Zeros01.x

Zeros02.x: NelZeros02.cpp intervalito02.hpp RutinasNelson02.hpp \
ParametrosGlobales.hpp
	g++ $(libs) NelZeros02.cpp -o Zeros02.x

NelStat01.x: NelsonStatistics.cpp RutinasNelson02.hpp ParametrosGlobales.hpp \
ChordAprox01.hpp
	g++ $(libs) NelsonStatistics.cpp -o NelStat01.x

NelCumulantes01.x: NelsonCumulantesArbitrarios01.cpp \
	RutinasNelson02.hpp ParametrosGlobales.hpp \
	ChordAprox01.hpp BinomialCoefficient01.hpp
	g++ $(libs) $< -o $@




Scars01.x: NelsonScars01.cpp intervalito02.hpp RutinasNelson02.hpp \
ParametrosGlobales.hpp
	g++ $(libs) NelsonScars01.cpp -o Scars01.x


Scars02.x: NelsonScars02.cpp intervalito02.hpp RutinasNelson02.hpp \
ParametrosGlobales.hpp
	g++ $(libs) NelsonScars02.cpp -o Scars02.x


NelsonPlan03.x: NelsonMaxwell01.cpp RutinasNelson03.hpp ParametrosGlobales.hpp
	g++ $(libs) NelsonMaxwell01.cpp -o NelsonPlan03.x

# Este lleva Deltas Negativas. Implementacion dirty, by the moment. 
NelsonPlan04.x: NelsonMaxwell02.cpp RutinasNelson03.hpp ParametrosGlobales.hpp
	g++ $(libs) NelsonMaxwell02.cpp -o NelsonPlan04.x


# Estos son programas en fortran para un quick and dirty use of numerical data

WignerfromWeyl01.x: fftsymplback01.f90
	gfortran fftsymplback01.f90 -o WignerfromWeyl01.x

ParalelFTforward01.x: openmpfftforward02.f90 RutinaFourierForParalel.f90 constantes01.f90
	gfortran openmpfftforward02.f90 RutinaFourierForParalel.f90 constantes01.f90 -fopenmp -O3 \
	-o ParalelFTforward01.x



