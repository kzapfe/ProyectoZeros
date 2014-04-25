LIBS=-lgsl -larmadillo
FUENTECOM= simplectic01.hpp QuantumConstants.hpp SeccionesPointillistas01.hpp 
FUENTESESFERA=PopulateSphericalShells01.hpp

TestPassVector01.x: TestPassVector01.cpp RutinaExternaEstupida01.hpp
		g++ $(LIBS)   $<  -o $@		

#Not really testing anymore
SphereParalel01.x: EsferaParalel01.cpp $(FUENTECOM)  $(FUENTESESFERA)
	g++ $(LIBS) -fopenmp   $<  -o $@	

PointsGenParalel01.x: PuntillistaGeneralParalel01.cpp $(FUENTECOM)  $(FUENTESESFERA)
	g++ $(LIBS) -fopenmp -larmadillo  $<  -o $@	
