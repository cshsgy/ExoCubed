#include "cia_read.hpp"

AthenaArray<Real> reform_read (std::string filename) {
    std::ifstream file {filename};  // open file
    AthenaArray<Real> data; // create storage array
    std::string line1; // storage for the first line
    std::string line2; // storage for the second line
    if (file.good()) {
        int nx; // number of spectral points, horizontal
        int ny; // number of temperature points, vertical
        file >> ny >> nx;
        std::getline(file, line1); // Skip the first line
        std::getline(file, line2); // Skip the second line
        data.NewAthenaArray(ny, nx);
        for (int j = 0; j < ny; ++j) {
            for (int i = 0; i < nx; ++i) {
                file >> data(j, i);
            }
        }
    } else {
        throw std::runtime_error("Unable to open " + filename);
    }
    return data;
}

//we are going to read the file twice, first to count # of rows and columns and second time 
AthenaArray<Real> ff_read (std::string filename) {
    int num_of_row = 0;
    int num_of_column = 1;
    AthenaArray<Real> data; // create storage array
    std::ifstream file {filename};  // open file
    if (file.good()) {
        std::string line;
        std::getline(file, line); //get the first line
        std::getline(file, line); //get the second line
        //calculate the # of columns base on the # of space
        for (int j = 0; j < line.size(); ++j) {
            if (line[j] == ' ') {
                ++ num_of_column;
            }
        }
        while (std::getline(file, line)) {
            ++ num_of_row;
        }
    } else {
        throw std::runtime_error("Unable to open " + filename);
    }
    file.close(); // close it
    file.open(filename);  // open it again
    if (file.good()) {
        int skiped;
        std::string line1; // storage for the first line
        std::string line2; // storage for the second line
        std::getline(file, line1); // Skip the first line
        std::getline(file, line2); // Skip the second line
        int ny = num_of_row;
        int nx = num_of_column;
        data.NewAthenaArray(ny, nx);
        for (int j = 0; j < ny; ++j) {
            file >> skiped; // skip first double
            for (int i = 0; i < nx; ++i) {
                file >> data(j,i);
            }
        }
    } else {
        throw std::runtime_error("Unable to open " + filename);
    }
    file.close();
    return data;
}