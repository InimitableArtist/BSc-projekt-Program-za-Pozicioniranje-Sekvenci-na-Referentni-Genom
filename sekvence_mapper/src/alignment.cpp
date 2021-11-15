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


            if (type == LOCAL) {
                
            }

            else if (type == GLOBAL) {

            }

            else if (type == SEMI_GLOBAL) {

            }


          }