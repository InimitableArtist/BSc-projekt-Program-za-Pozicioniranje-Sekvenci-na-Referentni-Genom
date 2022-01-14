#include <iostream>
#include <getopt.h>
#include <stdlib.h>
#include <string.h>
#include <memory>
#include <unordered_map>
#include <algorithm>
#include <cmath>

#include "src/alignment.h"
#include "src/minimizers.h"

bool cigar_flag;
int kmer_len = 15;
int window_len = 5;
double minimizer_freq = 0.001;
int match_cost;
int mismatch_cost;
int gap_cost;
int thread_num = 1;
AlignmentType type;


static std::string HELP = "-h or --help for help\n"
                           "-v or --version for version\n"
                           "-m value for matching\n"
                           "-n value for mismatching\n"
                           "-g value for gap\n"
                           "-a alignment type\n"
                           "-c cigar string enabled\n"
                           "-k k-mer length\n"
                           "-w window length\n"
                           "-f top f frequent minimizers that will not be taken in account\n";

class Sequence {
    public:
        Sequence(const char* name, std::uint32_t name_len,
                const char* data, std::uint32_t data_len):
                names(name, name_len),
                all_data(data, data_len)
                {}
                std::string names, all_data;
};

void make_minimizer_index(const std::unique_ptr<Sequence>& sequence,
                        std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>>& index) {
                              std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers = Minimize(sequence->all_data.c_str(),
                              sequence->all_data.size(), kmer_len, window_len);
                        for(auto minimizer : minimizers) {
                            index[std::get<0>(minimizer)].emplace_back(std::make_pair(std::get<1>(minimizer), std::get<2>(minimizer)));
                        }
}

std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>> remove_frequent_minimizers(std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>>& index,
int skip, std::vector<std::pair<unsigned int, unsigned int>> min_occ) {
    for (int i = 0; i < skip; i++) {
        index.erase(min_occ[i].second);
    }    
    return index;
}


void make_reference_index(const std::unique_ptr<Sequence>& sequence,
std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>>& index) {
    make_minimizer_index(sequence, index);
    std::vector<std::pair<unsigned int, unsigned int>> min_occ;
    min_occ.reserve(index.size());
    int num_singl = 0;
    for (auto& i : index) {
        if (i.second.size() == 1) {
            num_singl += 1;
        }
        min_occ.emplace_back(std::make_pair(i.second.size(), i.first));
    }
    std::sort(min_occ.begin(), min_occ.end(), std::greater<std::pair<unsigned int, unsigned int>>());
    size_t skip = std::ceil(index.size() * minimizer_freq);

    int index_size = index.size();
    if (skip >= index_size) {
        skip = index_size - 1;
    }
    index = remove_frequent_minimizers(index, skip, min_occ);


}


void display_version() {
    std::cout << "v0.1.0" << "\n";
}

void display_help() {
    std::cout <<  HELP << "\n";
}

int main(int argc, char* argv[]) {
    int option;
    const char* optstring = "m:g:n:a:k:w:f:t:hvc";
    while((option = getopt(argc, argv, optstring)) != -1) {
        switch (option)
        {
        case 'v':
            display_version();
            return 0;
            break;
        case 'h':
            display_help();
            return 0;
            break;
        case 'm':
            match_cost = atoi(optarg);
            break;
        case 'n':
            mismatch_cost = atoi(optarg);
            break;
        case 'g':
            gap_cost = atoi(optarg);
            break;
        case 'a':
            type = static_cast<AlignmentType>(atoi(optarg));
            break;
        case 'k':
            kmer_len = atoi(optarg);
            break;
        case 'w':
            window_len = atoi(optarg);
            break;
        case 'f':
            minimizer_freq = atof(optarg);
            break;
        case 't':
            thread_num = atoi(optarg);
            break;
        case 'c':
            cigar_flag = true;
            break;
            

        
        default:
            exit(1);
        }
        
    }

    if (optind < argc) {
        std::cout << argv[optind++] << "\n";
        std::cout << argv[optind] << "\n";
    }
    
    
    
    
    return 0;
}