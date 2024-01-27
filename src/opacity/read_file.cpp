#include "read_file.hpp"

pair <vector <double>, vector <double>> read_file (string file1, string file2) {
    std::ifstream file_1 {file1};  // open file
    std::ifstream file_2 {file2};  // open file
    vector <double> output_1;
    vector <double> output_2;
    // if stream OK = file readable
    if ( file_1.good() ) {
        double x;
        // as long as any 2 values readable
        while(file_1 >> x) {
        // print pairs (x,y) 
        output_1.push_back (x);
        }
    }
    if ( file_2.good() ) {
        double y;
        // as long as any 2 values readable
        while(file_2 >> y) {
        // print pairs (x,y) 
        output_2.push_back (y);
        }
    }
    pair <vector <double>, vector <double>> combined_output(output_1, output_2);
    return combined_output;
}