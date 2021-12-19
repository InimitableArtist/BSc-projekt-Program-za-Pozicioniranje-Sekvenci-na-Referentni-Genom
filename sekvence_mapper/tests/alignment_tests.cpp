#include <gtest/gtest.h>
#include "../src/alignment.h"


TEST(align_tests, local_align_test_1) {
    std::string* cigar = new std::string;
    unsigned int* target_begin = new unsigned int;
    int result = Align("ACCTAAGG", 8, "GGCTCAATCA", 10, LOCAL, 2, -1, -2, cigar, target_begin);
    EXPECT_EQ(result, 6);
    EXPECT_STREQ((*cigar).c_str(), "2M1I2M");
    EXPECT_EQ(*target_begin, 2);
    delete cigar;
    delete target_begin;
}
TEST(align_tests, local_align_test_2) {
    std::string* cigar = new std::string;
    unsigned int* target_begin = new unsigned int;
    int result = Align("GATCATATT", 9, "TCGTAGCG", 8, LOCAL, 2, -1, -2, cigar, target_begin);
    EXPECT_EQ(result, 7);
    EXPECT_STREQ((*cigar).c_str(), "2M1X2M");
    EXPECT_EQ(*target_begin, 0);
    delete cigar;
    delete target_begin;
}
TEST(align_tests, local_align_test_3) {
    std::string* cigar = new std::string;
    unsigned int* target_begin = new unsigned int;
    int result = Align("AGGTTG", 6, "TCAGTTGCC", 9, LOCAL, 1, -2, -2, cigar, target_begin);
    EXPECT_EQ(result, 4);
    EXPECT_STREQ((*cigar).c_str(), "4M");
    EXPECT_EQ(*target_begin, 3);
    delete cigar;
    delete target_begin;
}

TEST(align_tests, semi_global_align_test_1) {
    std::string* cigar = new std::string;
    unsigned int* target_begin = new unsigned int;
    int result = Align("CACTG", 5, "GGTTA", 5, SEMI_GLOBAL, 2, -1, -3, cigar, target_begin);
    EXPECT_EQ(result, 2);
    EXPECT_STREQ((*cigar).c_str(), "1M");
    EXPECT_EQ(*target_begin, 0);
    delete cigar;
    delete target_begin;
}