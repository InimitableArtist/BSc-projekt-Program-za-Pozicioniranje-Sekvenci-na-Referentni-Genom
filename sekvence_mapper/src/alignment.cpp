#include <iostream>
#include "alignment.h"
#include <algorithm>

int Align(const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    AlignmentType type,
    int match,
    int mismatch,
    int gap,
    std::string* cigar = nullptr,
    unsigned int* target_begin = nullptr) {   

    int V[query_len + 1][target_len + 1];
    int allign = 0;
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
    //Needleman-Wunsch
    else if (type == GLOBAL) {
        V[0][0] = 0;
        int poravnanje = 0;
        for (int i = 1; i <= query_len; i++) {
            V[i][0] = V[i-1][0] + gap;
        }
        for (int j = 1; j <= target_len; j++) {
            V[0][j] = V[0][j-1] + gap;
        }
        for (int i = 1; i <= query_len; i++) {
            for (int j = 1; j <= target_len; j++) {
                
                if (query[i] == target[j]) {
                    poravnanje = match;
                }
                else{
                    poravnanje = mismatch;
                }

                int var1 = V[i-1][j-1] + poravnanje; //M
                int var2 = V[i-1][j] + gap; //I
                int var3 = V[i][j-1] + gap; //D
                
                V[i][j] = std::max({var1, var2, var3});
            }
        }
        allign = V[query_len][target_len];
    }

    else if (type == SEMI_GLOBAL) {

    }

    return allign;
}