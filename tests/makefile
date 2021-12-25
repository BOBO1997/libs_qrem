test_qrem_filter: combinations.o hamming.o sgs_algorithm.o eigen_utils.o qrem_filter.o delta_filter.o least_norm_filter.o mooney_etal_filter.o test_qrem_filter.o
	g++ -std=c++14 combinations.o hamming.o sgs_algorithm.o eigen_utils.o qrem_filter.o delta_filter.o least_norm_filter.o mooney_etal_filter.o test_qrem_filter.o
test_qrem_filter.o: test_qrem_filter.cpp
	g++ -std=c++14 -c test_qrem_filter.cpp -o test_qrem_filter.o
delta_filter.o: ../cpp/delta_filter.cpp
	g++ -std=c++14 -c ../cpp/delta_filter.cpp -o delta_filter.o
least_norm_filter.o: ../cpp/least_norm_filter.cpp
	g++ -std=c++14 -c ../cpp/least_norm_filter.cpp -o least_norm_filter.o
mooney_etal_filter.o: ../cpp/mooney_etal_filter.cpp
	g++ -std=c++14 -c ../cpp/mooney_etal_filter.cpp -o mooney_etal_filter.o
qrem_filter.o: ../cpp/qrem_filter.cpp
	g++ -std=c++14 -c ../cpp/qrem_filter.cpp -o qrem_filter.o
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