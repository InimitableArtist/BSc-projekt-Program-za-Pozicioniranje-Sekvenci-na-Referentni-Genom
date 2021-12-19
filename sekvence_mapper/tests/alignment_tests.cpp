#include <gtest/gtest.h>
#include "../src/alignment.h"


TEST(AlignTests, LocalAlign1) {
    std::string* cigar = new std::string;
    unsigned int* target_begin = new unsigned int;
    int result = Align("ACCTAAGG", 8, "GGCTCAATCA", 10, LOCAL, 2, -1, -2, cigar, target_begin);
    EXPECT_EQ(result, 6);
    EXPECT_STREQ((*cigar).c_str(), "2M1I2M");
    EXPECT_EQ(*target_begin, 2);
    delete cigar;
    delete target_begin;
}