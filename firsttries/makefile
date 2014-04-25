libs=-lm -lgsl

Jogo.x:	coherentstate01.hpp  simplectic01.hpp Jogo01.cpp
	g++ $(libs) Jogo01.cpp -o Jogo.x

JogoSuper.x:	coherentstate01.hpp  simplectic01.hpp Jogo02.cpp coherentsuperposition.hpp
	g++ $(libs) Jogo02.cpp -o JogoSuper.x

ZeroZ01.x: coherentstate01.hpp  simplectic01.hpp ZeroZambrano01.cpp
	g++ $(libs) ZeroZambrano01.cpp -o ZeroZ01.x