#include <iostream>
#include <omp.h>

using namespace std;

int main() {
    omp_set_num_threads(3);
#pragma omp parallel
    {
        cout << "number: " << rand() << " " << endl;
    };
    return 0;
}