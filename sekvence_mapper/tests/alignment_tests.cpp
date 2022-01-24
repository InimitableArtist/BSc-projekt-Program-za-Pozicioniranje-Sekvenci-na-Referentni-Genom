#include <gtest/gtest.h>
#include "../src/alignment.h"


TEST(align_tests, local_align_test_1) {
    std::string* cigar = new std::string;
    unsigned int* target_begin = new unsigned int;
    const char* query = "ACCTAAGG";
    const char* target = "GGCTCAATCA";
    unsigned int query_len = 8;
    unsigned int target_len = 10;
    int result = Align(query, query_len, target, target_len, LOCAL, 2, -1, -2, cigar, target_begin);
    EXPECT_EQ(result, 6);
    EXPECT_STREQ((*cigar).c_str(), "2M1D2M");
    EXPECT_EQ(*target_begin, 2);
}
TEST(align_tests, local_align_test_2) {
    std::string* cigar = new std::string;
    unsigned int* target_begin = new unsigned int;
    const char* query = "GATCATATT";
    const char* target = "TCGTAGCG";
    unsigned int query_len = 9;
    unsigned int target_len = 8;
    int result = Align(query, query_len, target, target_len, LOCAL, 2, -1, -2, cigar, target_begin);
    EXPECT_EQ(result, 7);
    EXPECT_STREQ((*cigar).c_str(), "2M1X2M");
    EXPECT_EQ(*target_begin, 0);
}
TEST(align_tests, local_align_test_3) {
    std::string* cigar = new std::string;
    unsigned int* target_begin = new unsigned int;
    const char* query = "AGGTTG";
    const char* target = "TCAGTTGCC";
    unsigned int query_len = 6;
    unsigned int target_len = 9;
    int result = Align(query, query_len, target, target_len, LOCAL, 1, -2, -2, cigar, target_begin);
    EXPECT_EQ(result, 4);
    EXPECT_STREQ((*cigar).c_str(), "4M");
    EXPECT_EQ(*target_begin, 3);
}

TEST(align_tests, semi_global_align_test_1) {
    std::string* cigar = new std::string;
    unsigned int* target_begin = new unsigned int;
    const char* query = "CACGGCTTG";
    const char* target = "GGTT";
    unsigned int query_len = 9;
    unsigned int target_len = 4;
    int result = Align(query, query_len, target, target_len, SEMI_GLOBAL, 2, -1, -2, cigar, target_begin);
    EXPECT_EQ(result, 6);
    EXPECT_STREQ((*cigar).c_str(), "2M1I2M");
    EXPECT_EQ(*target_begin, 0);
}

TEST(align_tests, semi_global_align_test_2) {
    std::string* cigar = new std::string;
    unsigned int* target_begin = new unsigned int;
    const char* query = "ACTCCGAT";
    const char* target = "TCCG";
    unsigned int query_len = 8;
    unsigned int target_len = 4;
    int result = Align(query, query_len, target, target_len, SEMI_GLOBAL, 2, -1, -2, cigar, target_begin);
    EXPECT_EQ(result, 8);
    EXPECT_STREQ((*cigar).c_str(), "4M");
    EXPECT_EQ(*target_begin, 0);
}

TEST(align_tests, global_align_test_1) {
    std::string* cigar = new std::string;
    unsigned int* target_begin = new unsigned int;
    const char* query = "CACTG";
    const char* target = "GGTTA";
    unsigned int query_len = 5;
    unsigned int target_len = 5;
    int result = Align(query, query_len, target, target_len, GLOBAL, 2, -1, -3, cigar, target_begin);
    EXPECT_EQ(result, -2);
    EXPECT_STREQ((*cigar).c_str(), "3X1M1X");
    EXPECT_EQ(*target_begin, 0);
}

TEST(align_tests, global_align_test_2) {
    std::string* cigar = new std::string;
    unsigned int* target_begin = new unsigned int;
    const char* query = "CCTAGCGCCTATATTTGTCT";
    const char* target = "TAGCCATTTATCGGAGCGCC";
    unsigned int query_len = 20;
    unsigned int target_len = 20;
    int result = Align(query, query_len, target, target_len, GLOBAL, 2, -1, -2, cigar, target_begin);
    EXPECT_EQ(result, 4);
    EXPECT_STREQ((*cigar).c_str(), "2I4M1I1M1X1M1X3M2D2X1M1D1X1M1X");
    EXPECT_EQ(*target_begin, 0);
}

TEST(align_tests, global_align_test_3) {
    std::string* cigar = new std::string;
    unsigned int* target_begin = new unsigned int;
    const char* query = "CACA";
    const char* target = "CAAA";
    unsigned int query_len = 4;
    unsigned int target_len = 4;
    int result = Align(query, query_len, target, target_len, GLOBAL, 2, -1, -2, cigar, target_begin);
    EXPECT_EQ(result, 5);
    EXPECT_STREQ((*cigar).c_str(), "2M1X1M");
    EXPECT_EQ(*target_begin, 0);
}

TEST(align_tests, global_align_test_4) {
    std::string* cigar = new std::string;
    unsigned int* target_begin = new unsigned int;
    const char* query = "CACACA";
    const char* target = "CAAA";
    unsigned int query_len = 6;
    unsigned int target_len = 4;
    int result = Align(query, query_len, target, target_len, GLOBAL, 2, -1, -2, cigar, target_begin);
    EXPECT_EQ(result, 4);
    EXPECT_STREQ((*cigar).c_str(), "2M1I1M1I1M");
    EXPECT_EQ(*target_begin, 0);
}

TEST(align_tests, hirschberg_align_test_1) {
    std::string* cigar = new std::string;
    const char* query = "TATGC";
    const char* target = "AAGTACGCA";
    unsigned int query_len = 5;
    unsigned int target_len = 9;
    int result = 0;
    result = LinearnaSloz(query, query_len, target, target_len, 2, -1, -2, cigar);
    EXPECT_EQ(result, -1);
    EXPECT_STREQ((*cigar).c_str(), "3D2M1X2M1D");
    
}

TEST(align_tests, hirschberg_align_test_2) {
    std::string* cigar = new std::string;
    const char* query = "TATGC";
    const char* target = "AGTACGCA";
    unsigned int query_len = 5;
    unsigned int target_len = 8;
    int result = 0;
    result = LinearnaSloz(query, query_len, target, target_len, 2, -1, -2, cigar);
    EXPECT_EQ(result, 1);
    EXPECT_STREQ((*cigar).c_str(), "2D2M1X2M1D");
    
}

TEST(align_tests, hirschberg_align_test_3) {
    std::string* cigar = new std::string;
    const char* query = "CCTAGCGCCTATATTTGTCT";
    const char* target = "TAGCCATTTATCGGAGCGCC";
    unsigned int query_len = 20;
    unsigned int target_len = 20;
    int result = LinearnaSloz(query, query_len, target, target_len,  2, -1, -2, cigar);
    EXPECT_EQ(result, 4);
    EXPECT_STREQ((*cigar).c_str(), "2I4M1I1M1X1M1X3M2D2X1M1D1X1M1X");
    
}