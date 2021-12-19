#include <iostream>
#include "alignment.h"
#include <algorithm>
#include <string>

using namespace std;

std::pair<int, ParentTrack> find_max(int values[], int length, bool mtch) {
    int max = values[0];
    int index = 0;
    ParentTrack track;
    for (int i = 0; i < length; i++) {
        if (values[i] > max) {
            max = values[i];
            index = i;
        }
    }
    if (index == 0) {
        if (mtch) {
            track = MATCH;
        }
        else {
            track = MISMATCH;
        }
    }
    else if (index == 1) track = INSERTION;
    else if (index == 2) track = DELETION;
    else track = NONE;
    return std::make_pair(max, track);
}

std::pair<int, int> get_target_begin(int res_row,int res_column,struct Cell **V){
    ParentTrack now = V[res_row][res_column].parent;
    while(now != NONE){
        if(now == MATCH || now == MISMATCH) {
            res_row--;
            res_column--;
        }
        else if(now == INSERTION) res_column--;    
        else res_row--;
        
        now = V[res_row][res_column].parent;
    }
    return std::make_pair(res_row, res_column);
}

std::string get_cigar(int res_row,int res_column, Cell **V){
    string str = "";
    ParentTrack now = V[res_row][res_column].parent;
    ParentTrack prev = NONE;
    int dist=1;

    while(now != NONE){
        if(now == MATCH || now == MISMATCH) {
            res_row--;
            res_column--;
        }
        else if(now == INSERTION) res_column--;    
        else res_row--;
        prev = now;
        now = V[res_row][res_column].parent;
        
        if(prev == now) dist++;
        else{
            str += std::to_string(dist);
            if (prev == MATCH) str += "M";
            else if (prev == MISMATCH) str += "X";
            else if (prev == INSERTION) str += "I";
            else str += "D";
            dist = 1;
        }
    }
    return str;
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
    Cell **V = new Cell*[query_len + 1];
    for (unsigned int i = 0; i <= query_len; i++) {
        V[i] = new Cell[target_len + 1];
    }
    int align = 0;

    //Table coords of the result
    int res_row = 0;
    int res_column = 0;
    int match_cost = 0;

    //Smith-Waterman
    if (type == LOCAL) {
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
                bool mtch = false;
                if (query[i - 1] == target[j - 1]) {
                    match_cost = V[i - 1][j - 1].cost + match;
                    mtch = true;
                } 
                else match_cost = V[i - 1][j - 1].cost + mismatch;

                left = V[i][j - 1].cost + gap;
                up = V[i - 1][j].cost + gap;

                int values[4] = {match_cost, left, up, 0};
                std::pair<int, ParentTrack> res = find_max(values, 4, mtch);
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
                bool mtch = false;
                
                if (query[i] == target[j]) {
                    poravnanje = match;
                    mtch = true;
                }
                else{
                    poravnanje = mismatch;
                }

                int var1 = V[i-1][j-1].cost + poravnanje; //M
                int var2 = V[i-1][j].cost + gap; //I
                int var3 = V[i][j-1].cost + gap; //D

                int vars[3] = {var1, var2, var3};
                std::pair<int, ParentTrack> res = find_max(vars, 3, mtch);
                
                V[i][j].cost = res.first;
                V[i][j].parent = res.second;
            }
        }
        res_row = query_len;
        res_column = target_len;
        align = V[res_row][res_column].cost;
    }

    else if (type == SEMI_GLOBAL) {
        V[0][0].cost = 0;
        V[0][0].parent = NONE;
        int MAX = INT32_MIN;
        int left, up;

        for (int i = 1; i <= query_len; i++) {
            V[i][0].cost = 0;
            V[i][0].parent = NONE;
        }
        for (int j = 1; j <= target_len; j++) {
            V[0][j].cost = 0;
            V[0][j].parent = NONE;
        }

        for (int i = 1; i <= query_len; i++) {
            for (int j = 1; j <= target_len; j++) {
                bool mtch = false;
                if(query[i - 1] == target[j - 1]) {
                    match_cost = V[i - 1][j - 1].cost + match; 
                    mtch = true;
                } else {
                    match_cost = V[i - 1][j - 1].cost + mismatch;
                }
                left = V[i][j - 1].cost + gap;
                up = V[i - 1][j].cost + gap;

                int values[4] = {match_cost, left, up}; 
                std::pair<int, ParentTrack> res = find_max(values, 3, mtch);
                V[i][j].cost = res.first;
                V[i][j].parent = res.second;

                if (V[i][j].cost > MAX && (i == query_len) || (j == target_len)) {
                    MAX = V[i][j].cost;
                    res_row = i;
                    res_column = j;
                }
            }
        }

    }

    std::pair<int, int> coords = get_target_begin(res_row, res_column, V);
    if (target_begin != nullptr) {
        if (type == GLOBAL) {
            *target_begin = 0;
        }
        else if (type == SEMI_GLOBAL) {
            *target_begin = coords.second;
        }
        else {
            *target_begin = coords.second;
        }
    }
    *cigar = get_cigar(res_row, res_column, V);

    return align;
}