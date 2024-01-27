#include "read_file.hpp"

std::pair <std::vector <double>, std::vector <double>> read_file (std::string file1, std::string file2) {
    std::ifstream file_1 {file1};  // open file
    std::ifstream file_2 {file2};  // open file
    std::vector <double> output_1;
    std::vector <double> output_2;
    // if stream OK = file readable
        double x;
        // as long as next value readable
        while(file_1 >> x) {
            //put it into output1
            output_1.push_back (x);
        }
        
        double y;
        file_2.ignore(2);
        // as long as next value readable
        while(file_2 >> y) {
            //put it into output2
            output_2.push_back (y);
        }
    std::pair <std::vector <double>, std::vector <double>> combined_output = {output_1, output_2};
    return combined_output;
}