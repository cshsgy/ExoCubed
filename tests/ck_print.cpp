// C/C++
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

// external
#include <gtest/gtest.h>
#include <application/application.hpp>

// athena
#include <athena/athena_arrays.hpp>

// opacity
#include <opacity/absorber_ck.hpp>

std::string data_folder = "ck_data_01242024/ck/";
std::string output_folder = "output/";

TEST(LoadCoefficient, test_case1) {
    auto app = Application::GetInstance();

    // Open a file in write mode.
    auto outfile = app->FindResource(output_folder + "output.txt");
    auto correct_outfile = app->FindResource(output_folder + "correct_output.txt");
    auto data_file = app->FindResource(data_folder + "PM_ck_HELIO5K_cond_11_nOPT_wcia.txt");
    HeliosCKPremix ck_file("PM_ck_HELIO5K_cond_11_nOPT_wcia");

    std::ofstream outFile(outfile);

    // Call the function
    ck_file.LoadCoefficient(data_file, 1, outFile);

    // Read back the contents of the output files for comparison
    std::ifstream out(outfile);
    std::ifstream correct_out(correct_outfile);
    std::string out_contents((std::istreambuf_iterator<char>(out)), std::istreambuf_iterator<char>());
    std::string correct_out_contents((std::istreambuf_iterator<char>(correct_out)), std::istreambuf_iterator<char>());
    std::cout << out_contents;

    // Close the files to flush buffers
    outFile.close();

    // Compare the contents of the two files
    EXPECT_EQ(out_contents, correct_out_contents);
}

int main(int argc, char **argv) {
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}