CPP = g++
CPPFLAGS = -std=c++14

SRCS = combinations.cpp \
hamming.cpp \
sgs_algorithm.cpp \
harger_higham.cpp \
eigen_utils.cpp \
qrem_filter.cpp \
delta_filter.cpp \
least_norm_filter.cpp \
mooney_etal_filter.cpp \
nation_etal_filter.cpp \
test_qrem_filter.cpp \
dummy_data.cpp

OBJS = $(SRCS:%.cpp=%.o)

test_qrem_filter.out: $(OBJS)
	$(CPP) $(CPPFLAGS) -o $@ $^

test_qrem_filter.o: test_qrem_filter.cpp
	$(CPP) $(CPPFLAGS) -c test_qrem_filter.cpp -o test_qrem_filter.o

dummy_data.o: dummy_data.cpp
	$(CPP) $(CPPFLAGS) -c dummy_data.cpp -o dummy_data.o

delta_filter.o: ../../cpp/delta_filter.cpp
	$(CPP) $(CPPFLAGS) -c ../../cpp/delta_filter.cpp -o delta_filter.o

least_norm_filter.o: ../../cpp/least_norm_filter.cpp
	$(CPP) $(CPPFLAGS) -c ../../cpp/least_norm_filter.cpp -o least_norm_filter.o

mooney_etal_filter.o: ../../cpp/mooney_etal_filter.cpp
	$(CPP) $(CPPFLAGS) -c ../../cpp/mooney_etal_filter.cpp -o mooney_etal_filter.o

nation_etal_filter.o: ../../cpp/nation_etal_filter.cpp
	$(CPP) $(CPPFLAGS) -c ../../cpp/nation_etal_filter.cpp -o nation_etal_filter.o

qrem_filter.o: ../../cpp/qrem_filter.cpp
	$(CPP) $(CPPFLAGS) -c ../../cpp/qrem_filter.cpp -o qrem_filter.o

eigen_utils.o: ../../cpp/eigen_utils.cpp
	$(CPP) $(CPPFLAGS) -c ../../cpp/eigen_utils.cpp -o eigen_utils.o

hamming.o: ../../cpp/hamming.cpp
	$(CPP) $(CPPFLAGS) -c ../../cpp/hamming.cpp -o hamming.o

sgs_algorithm.o: ../../cpp/sgs_algorithm.cpp
	$(CPP) $(CPPFLAGS) -c ../../cpp/sgs_algorithm.cpp -o sgs_algorithm.o

combinations.o: ../../cpp/combinations.cpp
	$(CPP) $(CPPFLAGS) -c ../../cpp/combinations.cpp -o combinations.o

harger_higham.o: ../../cpp/harger_higham.cpp
	$(CPP) $(CPPFLAGS) -c ../../cpp/harger_higham.cpp -o harger_higham.o

clean:
	rm -rf *.o *.out

echo_objs: $(OBJS)
	@echo $^