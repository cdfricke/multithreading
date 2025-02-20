#include <iostream>
#include <thread>
#include "Integrate.h"

#define EX1


#ifdef EX1

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

        std::thread worker(DoWork);
        std::cin.get();
        s_Finished = true;

        worker.join();
        std::cout << "Finished!\n";
        std::cin.get();
        return 0;
    }
#endif