libs= -larmadillo  -lgsl

Puntos01.x: PuntosCuerdas01.cpp simplectic01.hpp RutinasExternas01.hpp
	g++ $(libs) PuntosCuerdas01.cpp -o Puntos01.x

Puntos02.x: PuntosCuerdas02.cpp simplectic01.hpp RutinasExternas01.hpp PuntosCuerdasRutina01.hpp  
	g++ $(libs) PuntosCuerdas02.cpp -o Puntos02.x

Puntos03.x: PuntosCuerdas03.cpp simplectic01.hpp RutinasExternas01.hpp PuntosCuerdasRutina02.hpp
	g++ $(libs) PuntosCuerdas03.cpp -o Puntos03.x

PRadial01.x: PuntosCuerdas04.cpp simplectic01.hpp RutinasExternas01.hpp PuntosCuerdasRadial01.hpp
	g++ $(libs) PuntosCuerdas04.cpp -o PRadial01.x


4Bolla01.x: simplectic01.hpp intervalitos01.hpp \
ComparaWeylenradio01.hpp RutinasExternas02.hpp PuntosCuerdasEsfera01.hpp \
Cuerda3esfera.cpp 
	g++ $(libs) Cuerda3esfera.cpp -o 4Bolla01.x

4Bolla02.x: simplectic01.hpp intervalitos01.hpp \
ComparaWeylenradio01.hpp RutinasExternas02.hpp PuntosCuerdasEsfera02.hpp \
Cuerda3esfera2.cpp 
	g++ $(libs) Cuerda3esfera2.cpp -o 4Bolla02.x

4Bolla03.x: simplectic01.hpp intervalitos01.hpp \
ComparaWeylenradio01.hpp RutinasExternas02.hpp PuntosCuerdasEsfera03.hpp \
Cuerda3esfera3.cpp 
	g++ $(libs) Cuerda3esfera3.cpp -o 4Bolla03.x
