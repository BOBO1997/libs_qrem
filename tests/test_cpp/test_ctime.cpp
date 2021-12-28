#include <ctime>
#include <iostream>

int main() {

        std::clock_t c_start = std::clock();
        // your_algorithm
        int sum = 0;
        for (int i = 0; i < 10000; i++) {
                sum += i;
        }
        std::clock_t c_end = std::clock();

        std::cout << "CPU time used: " << c_end - c_start << " ms\n";
        std::cout << "CPU time used: " << (double)(c_end - c_start) << " ms\n";

}
