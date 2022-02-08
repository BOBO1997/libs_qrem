CC=clang
CXX=clang++
CXXFLAGS=-std=c++17 -O3 -Wall -pthread -I $(pwd)/../../eigen

SRCS_LIB = ../../cpp/combinations.cpp \
../../cpp/hamming.cpp \
../../cpp/sgs_algorithm.cpp \
../../cpp/harger_higham.cpp \
../../cpp/eigen_utils.cpp \
../../cpp/qrem_filter.cpp \
../../cpp/ignis_filter.cpp \
../../cpp/delta_filter.cpp \
../../cpp/least_norm_filter.cpp \
../../cpp/mooney_etal_filter.cpp \
../../cpp/nation_etal_filter.cpp

SRCS_QREM_FILTER = test_qrem_filter.cpp \
dummy_data50.cpp \
# dummy_data12.cpp

SRCS_GENERAL_CAL_MAT = test_general_cal_mat.cpp \
dummy_data12.cpp \
# dummy_data50.cpp

OBJS_LIB = $(SRCS_LIB:%.cpp=%.o)
OBJS_QREM_FILTER = $(SRCS_QREM_FILTER:%.cpp=%.o)
OBJS_GENERAL_CAL_MAT = $(SRCS_GENERAL_CAL_MAT:%.cpp=%.o)

all: test_qrem_filter.out test_general_cal_mat.out

test_qrem_filter.out: $(OBJS_QREM_FILTER) $(OBJS_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $^

test_general_cal_mat.out: $(OBJS_GENERAL_CAL_MAT) $(OBJS_LIB)
	$(CXX) $(CXXFLAGS) -o $@ $^

$(OBJS_LIB): %.o : %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

test_general_cal_mat.o: test_general_cal_mat.cpp
	$(CXX) $(CXXFLAGS) -c test_general_cal_mat.cpp -o test_general_cal_mat.o

test_qrem_filter.o: test_qrem_filter.cpp
	$(CXX) $(CXXFLAGS) -c test_qrem_filter.cpp -o test_qrem_filter.o

dummy_data12.o: dummy_data12.cpp
	$(CXX) $(CXXFLAGS) -c dummy_data12.cpp -o dummy_data12.o

dummy_data50.o: dummy_data50.cpp
	$(CXX) $(CXXFLAGS) -c dummy_data50.cpp -o dummy_data50.o

clean:
	rm -rf *.o *.out
	rm -rf $(OBJS_LIB)

clean_general:
	rm -rf test_general_cal_mat.o*

clean_qrem:
	rm -rf test_qrem_filter.o*

echo_objs: $(OBJS)
	@echo $^
