#include "Dataset.h"
#include <chrono>
using namespace std::chrono;

int main() {
    srand(time(0));
    string fname;
    int points;
    cout << "Enter filename : ";
    cin >> fname;
    cout << "Enter Total Datapoints : ";
    cin >> points;
    Dataset X("TestFiles/" + fname, points);
    
    auto start = high_resolution_clock::now();
    X.FMST();
    auto stop = high_resolution_clock::now();

    auto duration = duration_cast<microseconds>(stop - start);
    cout << "Time taken by function: " << duration.count() << " microseconds" << endl;

    return 0;
}