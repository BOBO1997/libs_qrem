test_qrem_filter: combinations.o hamming.o sgs_algorithm.o mitigation.o test_qrem_filter.o
	g++ combinations.o hamming.o sgs_algorithm.o mitigation.o test_qrem_filter.o
test_qrem_filter.o: test_qrem_filter.cpp
	g++ -c test_qrem_filter.cpp -o test_qrem_filter.o
mitigation.o: mitigation.cpp
	g++ -c mitigation.cpp -o mitigation.o
hamming.o: hamming.cpp
	g++ -c hamming.cpp -o hamming.o
sgs_algorithm.o: sgs_algorithm.cpp
	g++ -c sgs_algorithm.cpp -o sgs_algorithm.o
combinations.o: combinations.cpp
	g++ -c combinations.cpp -o combinations.o
clean:
	rm *.o *.out