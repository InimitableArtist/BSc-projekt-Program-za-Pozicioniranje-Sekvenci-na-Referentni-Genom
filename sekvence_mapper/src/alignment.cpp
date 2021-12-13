#include <iostream>
#include "alignment.h"
#include <algorithm>

std::pair<int, ParentTrack> find_max(int values[], int length) {
    int max = values[0];
    int index = 0;
    ParentTrack track;
    for (int i = 0; i < length; i++) {
        if (values[i] > max) {
            max = values[i];
            index = i;
        }
    }
    if (index == 0) track = DIAGONAL;
    else if (index == 1) track = LEFT;
    else if (index == 2) track = UP;
    else track = NONE;
    return std::make_pair(max, track);
}

int Align(const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    AlignmentType type,
    int match,
    int mismatch,
    int gap,
    std::string* cigar = nullptr,
    unsigned int* target_begin = nullptr) {   
    
    //Dynamic programming table
    Cell V[query_len + 1][target_len + 1];
    int align = 0;

    //Table coords of the result
    int res_row = 0;
    int res_column = 0;

    //Smith-Waterman
    if (type == LOCAL) {
        int match_cost = 0;
        int MAX = 0;
        int left, up;
        V[0][0].cost = 0;
        V[0][0].parent = NONE;

        for (int i = 0; i <= query_len; i++) {
            V[i][0].cost = 0;
            V[i][0].parent = NONE;
        }
        for (int j = 0; j <= target_len; j++) {
            V[0][j].cost = 0;
            V[0][j].parent = NONE;
        }
        for (int i = 1; i <= query_len; i++) {
            for (int j = 1; j <= target_len; j++) {
                if (query[i - 1] == target[j - 1]) {
                    match_cost = V[i - 1][j - 1].cost + match;
                } 
                else match_cost = V[i - 1][j - 1].cost + mismatch;

                left = V[i][j - 1].cost + gap;
                up = V[i - 1][j].cost + gap;

                int values[4] = {match_cost, left, up, 0};
                std::pair<int, ParentTrack> res = find_max(values, 4);
                V[i][j].cost = res.first;
                V[i][j].parent = res.second;

                //Update MAX coords
                if (V[i][j].cost > MAX) {
                    MAX = V[i][j].cost;
                    res_row = i;
                    res_column = j;
                }


            }
        }
        align = V[res_row][res_column].cost;

    }
    //Needleman-Wunsch
    else if (type == GLOBAL) {
        V[0][0].cost = 0;
        V[0][0].parent = NONE;
        int poravnanje = 0;
        for (int i = 1; i <= query_len; i++) {
            V[i][0].cost = V[i-1][0].cost + gap;
        }
        for (int j = 1; j <= target_len; j++) {
            V[0][j].cost = V[0][j-1].cost + gap;
        }
        for (int i = 1; i <= query_len; i++) {
            for (int j = 1; j <= target_len; j++) {
                
                if (query[i] == target[j]) {
                    poravnanje = match;
                }
                else{
                    poravnanje = mismatch;
                }

                int var1 = V[i-1][j-1].cost + poravnanje; //M
                int var2 = V[i-1][j].cost + gap; //I
                int var3 = V[i][j-1].cost + gap; //D

                int vars[3] = {var1, var2, var3};
                std::pair<int, ParentTrack> res = find_max(vars, 3);
                
                V[i][j].cost = res.first;
                V[i][j].parent = res.second;
            }
        }
        align = V[query_len][target_len].cost;
    }

    else if (type == SEMI_GLOBAL) {

    }

    return align;
}