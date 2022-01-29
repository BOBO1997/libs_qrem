CPP = g++
CPPFLAGS = -std=c++17 -O3 -Wall -pthread

SRCS_LIB = ../../combinations.cpp \
../../hamming.cpp \
../../sgs_algorithm.cpp \
../../harger_higham.cpp \
../../eigen_utils.cpp \
../../qrem_filter.cpp \
../../ignis_filter.cpp \
../../delta_filter.cpp \
../../least_norm_filter.cpp \
../../mooney_etal_filter.cpp \
../../nation_etal_filter.cpp

SRCS_QREM_FILTER = test_qrem_filter.cpp \
dummy_data50.cpp \
# dummy_data12.cpp

SRCS_GENERAL_CAL_MAT = test_general_cal_mat.cpp \
dummy_data12.cpp \
# dummy_data50.cpp

OBJS_LIB = $(SRCS_LIB:%.cpp=%.o)
OBJS_QREM_FILTER = $(SRCS_QREM_FILTER:%.cpp=%.o)
OBJS_GENERAL_CAL_MAT = $(SRCS_GENERAL_CAL_MAT:%.cpp=%.o)

all: test_qrem_filter test_general_cal_mat

test_qrem_filter.out: $(OBJS_QREM_FILTER) $(OBJS_LIB:../../%=%)
	$(CPP) $(CPPFLAGS) -o $@ $^

test_general_cal_mat.out: $(OBJS_GENERAL_CAL_MAT) $(OBJS_LIB:../../%=%)
	$(CPP) $(CPPFLAGS) -o $@ $^

# %.o: $(SRCS:%.cpp=%.cpp)
# 	$(CPP) $(CPPFLAGS) -c $@ -o $^

# %.o: $(LIB_SRCS:%.cpp=%.cpp)
# 	$(CPP) $(CPPFLAGS) -c $@ -o $^

test_general_cal_mat.o: test_general_cal_mat.cpp
	$(CPP) $(CPPFLAGS) -c test_general_cal_mat.cpp -o test_general_cal_mat.o

test_qrem_filter.o: test_qrem_filter.cpp
	$(CPP) $(CPPFLAGS) -c test_qrem_filter.cpp -o test_qrem_filter.o

dummy_data12.o: dummy_data12.cpp
	$(CPP) $(CPPFLAGS) -c dummy_data12.cpp -o dummy_data12.o

dummy_data50.o: dummy_data50.cpp
	$(CPP) $(CPPFLAGS) -c dummy_data50.cpp -o dummy_data50.o

ignis_filter.o: ../../cpp/ignis_filter.cpp
	$(CPP) $(CPPFLAGS) -c ../../cpp/ignis_filter.cpp -o ignis_filter.o

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

clean_general:
	rm -rf test_general_cal_mat.o*

clean_qrem:
	rm -rf test_qrem_filter.o*

echo_objs: $(OBJS)
	@echo $^