#include <iostream>
#include "alignment.h"

int Align(const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    AlignmentType type,
    int match,
    int mismatch,
    int gap,
    std::string* cigar = nullptr,
    unsigned int* target_begin = nullptr) {   

    int V[query_len + 1][target_len + 1];

    //Smith-Waterman
    if (type == LOCAL) {
        int M = 0;
        V[0][0] = 0;

        for (int i = 0; i <= query_len; i++) {
            V[i][0] = 0;
        }
        for (int j = 0; j <= target_len; j++) {
            V[0][j] = 0;
        }
        for (int i = 1; i <= query_len; i++) {
            for (int j = 1; j <= target_len; j++) {
            
            }
        }

    }

    else if (type == GLOBAL) {

    }

    else if (type == SEMI_GLOBAL) {

    }


}