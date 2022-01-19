#include <iostream>

enum AlignmentType {GLOBAL, LOCAL, SEMI_GLOBAL};
enum ParentTrack {NONE, INSERTION, DELETION, MATCH, MISMATCH};

struct Cell {
    int cost;
    ParentTrack parent;
};

std::pair<int, ParentTrack> find_max(int values[], int length, bool mtch);

std::pair<int, int> get_target_begin(int res_row,int res_column,struct Cell **V);

std::string get_cigar(int res_row,int res_column, Cell **V);


int Align(const char* query, unsigned int query_len,
          const char* target, unsigned int target_len,
          AlignmentType type, int match, int mismatch, int gap,
          std::string* cigar,
          unsigned int* target_begin);

int LinearnaSloz(const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    int match,
    int mismatch,
    int gap,
    std::string* cigar);
