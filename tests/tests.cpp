#include "../src/include/new_abaqus.h"
#include "../src/include/node.h"
#include "../src/include/mid.h"
#include "../src/include/pid.h"
#include "../src/include/misc_string_functions.h"

#include <gtest/gtest.h>
#include <cstdio>
#include <string>
#include <iostream>
#include <fstream>


class unit_tests : public ::testing::Test
{
};

TEST(analysis,example2_tuning_fork)
{
    New_abaqus na;
    na.mesh.set_analysis_name("example2_tuning_fork");
    na.mesh.read_file_new_method("example_runfiles/example2_tuning_fork.inp");
    na.mesh.assemble();
    na.mesh.solve();
    // Check for correct results in logfile
    std::ifstream infile;
    infile.open("example2_tuning_fork.log");
    for (std::string line; std::getline(infile, line); ) 
    {
        if(line.find("E I G E N V A L U E    O U T P U T") != std::string::npos)
        {
            // Correct answers:
            std::array<float, 9> correct_values{170.92, 180.57, 321.48, 429.65, 996.12, 1161.13, 1208.24, 1749.29, 1807.28};
            // Advance a few lines and check values
            for (size_t i = 0; i < 6; i++)
            {
                std::getline(infile, line);
            }
            auto values = misc::split_on(line, ' ');
            for (size_t i = 0; i < correct_values.size(); i++)
            {
                auto value = std::stof(misc::split_on(line, ' ').at(2));
                
                std::getline(infile, line);
            }
            
        }
    }
}

TEST(analysis,example1_2D)
{
    New_abaqus na;
    na.mesh.set_analysis_name("example1_2D");
    na.mesh.read_file_new_method("example_runfiles/example1_2D.inp");
    na.mesh.assemble();
    na.mesh.solve();
    // Check vtk file, last line is displacement as load node
    // 0.00646062 -0.0342598 0
    std::ifstream infile;
    infile.open("example1_2D.vtk");
    // infile.seekg(0,ios_base::end)
    // använd getline en jävla massa gånger bara lol!
    size_t data_row = 3337;
    std::string line;
    for (size_t i = 0; i < data_row + 1; i++)
    {
        std::getline(infile, line);
        if (i == data_row)
        {
            auto values = misc::split_on(line, ' ');
            ASSERT_FLOAT_EQ(0.00646062f, std::stof(values[0]));
            ASSERT_FLOAT_EQ(-0.0342598f, std::stof(values[1]));
            ASSERT_FLOAT_EQ(0.f, std::stof(values[2]));
        }
    }
}

TEST(unit_tests,node)
{
    int  id = 1;
    float x = 0.0f;
    float y = 0.0f;
    float z = 0.0f;
    Node node(id,0.0f,0.0f,0.0f);
    unsigned int node_counter = 1;
    ASSERT_FLOAT_EQ(id,node.id);
    ASSERT_FLOAT_EQ(node.original_x,x);
    ASSERT_FLOAT_EQ(node.original_y,y);
    ASSERT_FLOAT_EQ(node.original_z,z);
    int id2 = 2;
    Node node_2d(id2,x,y);
    ASSERT_FLOAT_EQ(node_2d.original_z,0.0f);
}

TEST(unit_tests,mid)
{
    auto mid      = std::make_unique<Mid>("mid");
    auto E        = 210000.0f;
    auto v        = 0.3f;
    auto density  = 7.2e-9;
    mid->set_E(E);
    mid->set_v(v);
    mid->set_density(density);
    EXPECT_FLOAT_EQ(mid->get_E(),E);
    EXPECT_FLOAT_EQ(mid->get_v(),v);
    EXPECT_FLOAT_EQ(mid->get_density(),density);
}

TEST(unit_tests,pid){
    auto mid      = std::make_unique<Mid>("mid");
    mid->set_E(210000);
    mid->set_v(0.3f);
    mid->set_density(7.2e-9);
    auto pid  = std::make_unique<Pid>("pid",mid.get());
    EXPECT_EQ(pid->get_mid(),mid.get());
}


