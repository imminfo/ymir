//#include <omp.h>
//#include <stdio.h>
//
//int main() {
//#pragma omp parallel
//    printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
//}


#include <iostream>

#include "assemblygraph.h"

using namespace std;

int main() {
    /*
    argument parsing
     */
    // parse input arguments and choose what to do


    /*
    artifical sequences assembling
     */
    // parse input folder with model
    // generate repertoire
    // optionally compute probabilities of generated sequences
    // write the repertoire to the output file


    /*
    probabilities computing
     */
    // parse input with repertoire
    // parse input folder with model
    // get layouts for all sequences which haven't been aligned


    /*
    statistical inference
     */
    // parse input with repertoire
    // parse input folder with model
    // get layouts for all sequences which haven't been aligned
    // ???


    cout << "Hello, World!" << endl;
    return 0;
}