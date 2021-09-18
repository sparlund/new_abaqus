#include "../src/new_abaqus.h"
#include "../src/node.h"
#include "../src/mid.h"
#include "../src/pid.h"
#include "../src/elements/CPS4.h"

#include <gtest/gtest.h>
#include <cstdio>

//  g++ unit_tests.cpp -o unit_tests /usr/src/gtest/src/gtest_main.cc /usr/src/gtest/src/gtest-all.cc ../build/{node,dof,mid,pid,mesh}.o -I /usr/include -I /usr/src/gtest -L /usr/local/lib -lpthread 

class unit_tests : public ::testing::Test
{
protected:
    static void SetUpTestSuite() {
        //code here
    }

    static void TearDownTestSuite() {
        //code here
    }
    // TODO: why doesnt this work?
    // void SetUp() override {
        // std::freopen("testlog.txt", "w", stdout);
        // std::cout << "aa" << std::endl;
    // };
    // void TearDown() override {
        // std::fclose(stdout);
    // };

};

TEST(elements,CPS4)
{
    New_abaqus na;
    na.mesh.read_file_new_method("resources/CPS4.inp");
    na.mesh.assemble();
    na.mesh.print_matrix_to_mtx(*na.mesh.get_K(),"test_output/CPS4.mtx");

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
    ASSERT_FLOAT_EQ(node.x,x);
    ASSERT_FLOAT_EQ(node.y,y);
    ASSERT_FLOAT_EQ(node.z,z);
    int id2 = 2;
    Node node_2d(id2,x,y);
    ASSERT_FLOAT_EQ(node_2d.z,0.0f);
}

TEST(unit_tests,reading_file){
    // Read include file and check that we found the expected entities
    New_abaqus na;
    na.mesh.set_analysis_name("test");
    na.mesh.read_file_new_method("resources/simple_inputfile.inp");
    ASSERT_EQ(6,na.mesh.get_number_of_nodes());
    auto pid = na.mesh.get_pid_by_name("pid");
    auto mid = na.mesh.get_mid_by_name("mid");
    // these are the values from the input file
    ASSERT_FLOAT_EQ(pid->get_mid()->get_E(),210000.0f);
    ASSERT_FLOAT_EQ(pid->get_mid()->get_v(),0.3);
    ASSERT_FLOAT_EQ(pid->get_mid()->get_density(),7.85E-9);
    ASSERT_EQ(na.mesh.get_number_of_elements(),2);
    ASSERT_EQ(na.mesh.get_bc()->size(),4);
    auto nsets = na.mesh.get_nsets();
    auto esets = na.mesh.get_esets();
    ASSERT_EQ(na.mesh.get_nsets()->size(),1);
    na.mesh.assemble();
    na.mesh.solve();
    // TODO: create logic for element sets
    // EXPECT_EQ(na.mesh.get_esets()->size(),1);
    // EXPECT_EQ(na.mesh.f.nonZeros(),1);
    // ASSERT_EQ(na.mesh.f.sum(),-1e5);
    ASSERT_EQ(na.mesh.get_number_of_dofs(),18);
    auto f_to_be_added = na.mesh.get_f_to_be_added();
    ASSERT_EQ(f_to_be_added.size(),1);
    ASSERT_EQ(f_to_be_added.at(0).first,16);
    ASSERT_FLOAT_EQ(f_to_be_added.at(0).second,-1e5f);
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


