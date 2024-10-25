all:
	g++ -std=c++17 -O2 -o tour tour.cpp
	./tour < graph.txt
	rm tour
