test_qrem_filter: combinations.o hamming.o sgs_algorithm.o eigen_utils.o qrem_filter_base.o qrem_filter.o qrem_filter_mooney_etal.o test_qrem_filter.o
	g++ -std=c++14 combinations.o hamming.o sgs_algorithm.o eigen_utils.o qrem_filter_base.o qrem_filter.o qrem_filter_mooney_etal.o test_qrem_filter.o
test_qrem_filter.o: test_qrem_filter.cpp
	g++ -std=c++14 -c test_qrem_filter.cpp -o test_qrem_filter.o
qrem_filter.o: ../cpp/qrem_filter.cpp
	g++ -std=c++14 -c ../cpp/qrem_filter.cpp -o qrem_filter.o
qrem_filter_mooney_etal.o: ../cpp/qrem_filter_mooney_etal.cpp
	g++ -std=c++14 -c ../cpp/qrem_filter_mooney_etal.cpp -o qrem_filter_mooney_etal.o
qrem_filter_base.o: ../cpp/qrem_filter_base.cpp
	g++ -std=c++14 -c ../cpp/qrem_filter_base.cpp -o qrem_filter_base.o
eigen_utils.o: ../cpp/eigen_utils.cpp
	g++ -std=c++14 -c ../cpp/eigen_utils.cpp -o eigen_utils.o
hamming.o: ../cpp/hamming.cpp
	g++ -std=c++14 -c ../cpp/hamming.cpp -o hamming.o
sgs_algorithm.o: ../cpp/sgs_algorithm.cpp
	g++ -std=c++14 -c ../cpp/sgs_algorithm.cpp -o sgs_algorithm.o
combinations.o: ../cpp/combinations.cpp
	g++ -std=c++14 -c ../cpp/combinations.cpp -o combinations.o
clean:
	rm *.o *.out