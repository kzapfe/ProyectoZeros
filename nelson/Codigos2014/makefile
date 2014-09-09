##Todos los programas AHORA van a requerir usar los mismos Centros
## Si los centros NO existen, se corre el programa para crearlos

libs= -larmadillo  -lgsl
incluir= -I /usr/include/gsl
FuentesCom=ParametrosGlobales.hpp simplectic01.hpp


NelCumulantes01.x: NelsonCumulantesArbitrarios01.cpp \
	$(FuentesCom) ChordAprox01.hpp BinomialCoefficient01.hpp \
	 #(Todos los programas usan los mismos Centros!)
	g++ $(libs) $< -o $@


CentrosWigner.dat: JustPopulate01.x
	JustPopulate01.x

JustPopulate01.x: PoblarNelson01.cpp $(FuentesCom)
	g++ $(libs) $< -o $@

JustPopulateRight01.x: LlenarNelsonMonteCarlo02.cpp
	g++ $(libs) $< -o $@	

PopulateWithEnergy01.x: PopulateArgument01.cpp
	g++ $(libs) $< -o $@	


JustCuerdas01.x: CuerdasFromExteriorFile01.cpp $(FuentesCom)
	g++ $(libs) $< -o $@	

