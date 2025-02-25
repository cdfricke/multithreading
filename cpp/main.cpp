#include <iostream>
#include <thread>
#include <vector>
#include "Stats.h"


static bool s_Finished = false;

void DoWork() {
    using namespace std::literals::chrono_literals;
        
    while (!s_Finished) {
        std::cout << "Working...\n"; 
        std::this_thread::sleep_for(1s);
    }
    return;
}

int main() {

    std::vector<std::vector<double>> mat = {{1.0, 5.0, 9.0},
                                            {2.0, 6.0, 10.0},
                                            {4.0, 8.0, 12.0}};
    std::vector<double> means = stats::mean(mat, 1);
    for (auto val: means) {
        std::cout << val << std::endl;
    }
}
    