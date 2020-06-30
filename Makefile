
output: main.o particle.o force.o box.o
	g++ main.o particle.o force.o box.o -o output -O3

main.o: main.cpp
	g++ -c main.cpp -O3

particle.o: particle.cpp
	g++ -c particle.cpp -O3

force.o: force.cpp
	g++ -c force.cpp -O3

box.o: box.cpp
	g++ -c box.cpp -O3

clean: 
	rm *.o output

