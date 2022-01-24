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
#include <bits/stdc++.h>


#include "bioparser/fasta_parser.hpp"


using namespace std;

int countDigit(unsigned long int n)
{
    int count = 0;
    while (n != 0)
    {
        n = n / 10;
        ++count;
    }
    return count;
}
void longestSubsequence(unsigned long int sequence, int n)
{
    int a[n]; // integer u int[]

    int count = n;
    while (sequence != 0)
    {
        a[count - 1] = sequence % 10;
        sequence = sequence / 10;
        count--;
    }

    unordered_map<int, int> mp;


    int dp[n];
    memset(dp, 0, sizeof(dp));

    int maximum = INT_MIN;


    
    int index = -1;
    for (int i = 0; i < n; i++)
    {

        if (mp.find(a[i] - 1) != mp.end())
        {

            int lastIndex = mp[a[i] - 1] - 1;

            dp[i] = 1 + dp[lastIndex];
        }
        else
            dp[i] = 1;

        mp[a[i]] = i + 1;

        if (maximum < dp[i])
        {
            maximum = dp[i];
            index = i;
        }
    }

    for (int curr = a[index] - maximum + 1;
         curr <= a[index]; curr++)
        cout << curr << " ";

}

bool cigar_flag;
int kmer_len = 10;
int window_len = 15;
double minimizer_freq = 0.001;
int match_cost = 1;
int mismatch_cost = - 1;
int gap_cost = - 1;
int thread_num = 1;
AlignmentType type = GLOBAL;
int algorithm;

static std::string HELP = "-h or --help for help\n"
                          "-v or --version for version\n"
                          "-m value for matching\n"
                          "-n value for mismatching\n"
                          "-g value for gap\n"
                          "-a alignment type\n"
                          "-c cigar string enabled\n"
                          "-k k-mer length\n"
                          "-w window length\n"
                          "-f top f frequent minimizers that will not be taken in account\n"
                          "query name of the file with redings in FASTQ format\n"
                          "target name of the reference file in FASTA format\n";

class Sequence
{
public:
    Sequence(const char *name, std::uint32_t name_len,
             const char *data, std::uint32_t data_len) : names(name, name_len),
                                                         all_data(data, data_len)
    {
    }
    std::string names, all_data;
};

void make_minimizer_index(const std::unique_ptr<Sequence> &sequence,
                          std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>> &index)
{
    std::vector<std::tuple<unsigned int, unsigned int, bool>> minimizers = Minimize(sequence->all_data.c_str(),
                                                                                    sequence->all_data.size(), kmer_len, window_len);
    for (auto minimizer : minimizers)
    {
        index[std::get<0>(minimizer)].emplace_back(std::make_pair(std::get<1>(minimizer), std::get<2>(minimizer)));
    }
}

std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>> remove_frequent_minimizers(std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>> &index,
                                                                                                        int skip, std::vector<std::pair<unsigned int, unsigned int>> min_occ)
{
    for (int i = 0; i < skip; i++)
    {
        index.erase(min_occ[i].second);
    }
    return index;
}

void make_reference_index(const std::unique_ptr<Sequence> &sequence,
                          std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>> &index)
{
    
    make_minimizer_index(sequence, index);
    
    std::vector<std::pair<unsigned int, unsigned int>> min_occ;
    min_occ.reserve(index.size());
    int num_singl = 0;
    for (auto &i : index)
    {
        if (i.second.size() == 1)
        {
            num_singl += 1;
        }
        min_occ.emplace_back(std::make_pair(i.second.size(), i.first));
    }
    std::sort(min_occ.begin(), min_occ.end(), std::greater<std::pair<unsigned int, unsigned int>>());
    size_t skip = std::ceil(index.size() * minimizer_freq);

    int index_size = index.size();
    if (skip >= index_size)
    {
        skip = index_size - 1;
    }
    for (int i = 0; i < skip; i++) {
        index.erase(min_occ[i].second);
    }
    //index = remove_frequent_minimizers(index, skip, min_occ);
}
bool compare_matches(const tuple<bool, int, unsigned int>& a1, const tuple<bool, int, unsigned int>& a2) {
    if (get<0>(a1) != get<0>(a2)) {
        return !get<0>(a1);
    } else {
        if (get<1>(a1) != get<1>(a2)) {
            return get<1>(a1) < get<1>(a2);
        } else {
            return get<2>(a1) < get<2>(a2);
        }
    }
}

void best_match_cluster(unordered_map<unsigned int, vector<pair<unsigned int, bool>>>& fragment_index,
                        unordered_map<unsigned int, vector<pair<unsigned int, bool>>>& reference_index,
                        vector<tuple<bool, int, unsigned int>>& match_cluster) {
        
    vector<tuple<bool, int, unsigned int>> matches;
    for (auto e : fragment_index) {
        if (reference_index.count(e.first) != 0) {
            for (auto location : e.second) {
                for (auto r_location : reference_index[e.first]) {
                    bool strand = location.second ^ r_location.second;
                    int rel_pos;
                    if (strand) {
                        rel_pos = location.first + r_location.first;
                    } else {
                        rel_pos = location.first - r_location.first;
                    }
                    matches.emplace_back(make_tuple(strand, rel_pos, r_location.first));
                }
            }
        }
    }
    if (!matches.empty()) {
        sort(matches.begin(), matches.end(), compare_matches);
        
    }
                        }


void display_version()
{
    std::cout << "v0.1.0"
              << "\n";
}

void display_help()
{
    std::cout << HELP
              << "\n";
}

int main(int argc, char *argv[])
{
    int option;
    const char *optstring = "m:g:n:a:k:w:f:t:hvc";

    while ((option = getopt(argc, argv, optstring)) != -1)
    {
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
            algorithm = atoi(optarg);
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
            std::cout<<HELP;
                     
            exit(1);
        }
    }

    if (optind < argc)
    {
        std::string path = argv[optind];   
        
    }
	
	auto ref = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[argc - 2]);
    auto reference = ref->Parse(-1);
    int ref_size = (int)reference.size();

    auto frag = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[argc - 1]);
    auto fragments = frag->Parse(-1);
    int frag_size = (int)fragments.size();
    
    int sum = 0;
    for (int i = 0; i < frag_size; i++) {
        sum += fragments[i]->all_data.size();
    } 
    float avg_size = sum / frag_size;
    std::vector<size_t> fragment_vector;
    for (int i = 0; i < frag_size; i++) {
        fragment_vector.push_back(fragments[i]->all_data.size());
    }
    std::sort(fragment_vector.begin(), fragment_vector.end());

    
    
    string cigar;
    unsigned int target_begin;

    
    
    srand(time(NULL));
    int query_index = rand() % (frag_size);
    int target_index = rand() % (frag_size);

    //query_index = 0;
    //target_index = 0;

    switch (algorithm)
    {
        case 0:
            type = GLOBAL;
            break;
        case 1:
            type = LOCAL;
            break;
        case 2:
            type = SEMI_GLOBAL;
            break;
        default:
            break;
    }
    //cout << "type: " << type << "\n";
    int align_score = Align(fragments[query_index]->all_data.c_str(), fragments[query_index]->all_data.size(),
                                                fragments[target_index]->all_data.c_str(),
                                                fragments[target_index]->all_data.size(), type, match_cost, mismatch_cost, gap_cost, &cigar, &target_begin);
    cout << "Alignment score: " << align_score << "\n";
    if (cigar_flag) {
        cout << "Cigar: " << cigar << "\n";
    }

    std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>> reference_index;
    make_reference_index(reference.front(), reference_index);
    
    string res = "";
    for (int i = 1; i < 20; i++) {
        unordered_map<unsigned int, vector<pair<unsigned int, bool>>> fragment_index;
        make_minimizer_index(fragments[i], fragment_index);
        vector<tuple<bool, int, unsigned int>> match_cluster;
        best_match_cluster(fragment_index, reference_index, match_cluster);
    }
    

    
    return 0;
}