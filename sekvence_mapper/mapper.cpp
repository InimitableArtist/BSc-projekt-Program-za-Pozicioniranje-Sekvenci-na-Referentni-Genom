#include <iostream>
#include <getopt.h>

#include <stdlib.h>
#include <string.h>
#include <memory>
#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <stdio.h>
#include <fstream>

#include "src/alignment.h"
#include "src/minimizers.h"
#include <bits/stdc++.h>


#include "bioparser/fasta_parser.hpp"
#include "thread_pool/thread_pool.hpp"

const int EPSILON = 100;

using namespace std;

int get_ceil_index(vector<tuple<bool, int, unsigned int>>& matches,
                    vector<int>& T, int l, int r, int k) {
    
    while (r - l > 1) {
        int m = l+ (r - l) / 2;
        if (get<2>(matches[T[m]]) >= k) {
            r = m;
        } else {
            l = m;
        }
    }


    return r;
                    }
vector<tuple<bool, int,unsigned int>> longest_increasing_subsequence(vector<tuple<bool, int, unsigned int>>& matches,
                                    vector<int>& tail_indices,
                                    vector<int>& prev_indices) {   
    
    int len = 1;
    
    vector<tuple<bool, int, unsigned int>> result;
    for (int i = 1; i < matches.size(); i++) {
        if (get<2>(matches[i]) < get<2>(matches[tail_indices[0]])) {
            tail_indices[0] = i;
        }
        else if (get<2>(matches[i]) > get<2>(matches[tail_indices[len - 1]])) {
            prev_indices[i] = tail_indices[len - 1];
            tail_indices[len++] = i;
        } else {
            int position = get_ceil_index(matches, tail_indices, -1, len - 1, get<2>(matches[i]));
            cout << "\nPOZIIICIJA: " << position << "\n";
            if (position == 0) {
                prev_indices[i] = -1;
            } else {
                prev_indices[i] = tail_indices[position - 1];
            }
            tail_indices[position] = i;
        }
    }
    vector<int> lis_ind;
    lis_ind.reserve(len);
    for (int i = tail_indices[len - 1]; i >= 0; i = prev_indices[i]) {
        lis_ind.emplace_back(i);
    }
    reverse(lis_ind.begin(), lis_ind.end());
    for (auto a : lis_ind) {
        result.emplace_back(matches[a]);
    }
    
    return result;

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
bool optimize = false;

static std::string HELP = "-h or --help for help\n"
                          "-v or --version for version\n"
                          "-m value for matching\n"
                          "-n value for mismatching\n"
                          "-g value for gap\n"
                          "-a alignment type (0 - GLOBAL, 1 - LOCAL, 2 - SEMI-GLOBAL)\n"
                          "-c cigar string enabled\n"
                          "-k k-mer length\n"
                          "-w window length\n"
                          "-f top f frequent minimizers that will not be taken in account\n"
                          "-o turn the optimization algorithm on (only for GLOBAL align)\n"
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
    
    //cout << "velicina minimizera: " << minimizers.size() << "\n";
    if (minimizers.size() > kmer_len * window_len) return;
    for (auto minimizer : minimizers)
    {
        //cout << "minimizer: " << std::get<0>(minimizer) << "\n";
        


        index[std::get<0>(minimizer)].emplace_back(std::make_pair(std::get<1>(minimizer), std::get<2>(minimizer)));
        
    }
    //cout << "prosao sam petlju...\n";
    
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
    cerr << "Making reference index...\n"; 
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
    cerr << "Sorting minimizers by occurence...\n";
    std::sort(min_occ.begin(), min_occ.end(), std::greater<std::pair<unsigned int, unsigned int>>());
    size_t skip = std::ceil(index.size() * minimizer_freq);

    int index_size = index.size();
    if (skip >= index_size)
    {
        skip = index_size - 1;
    }

    cerr << "Number of distinct minimizers: " << index.size() << "\n";
    cerr << "Fraction of singletons: " << ((double) num_singl / index.size()) << "\n"; 
    cerr << "NUmber of occurences of the most frequent minimizers with top " << minimizer_freq * 100 << "% most frequent ignored: " << min_occ[skip].first << "\n";

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

    //vector<tuple<bool, unsigned int, int>> result;
        
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
        vector<tuple<bool, int, unsigned int>> current_cluster;
        int max = 0;
        current_cluster.emplace_back(matches[0]);
        for (int i = 1; i <= int(matches.size()); i++) {
            if (i == int(matches.size()) || get<0>(matches[i]) != get<0>(matches[i - 1]) || get<1>(matches[i]) - get<1>(matches[i - 1]) >= EPSILON) {
                vector<int> list_ind;
                int n = int(current_cluster.size());
                vector<int> tail_indices(n, 0);
                vector<int> prev_indices(n, -1);
                current_cluster = longest_increasing_subsequence(current_cluster, tail_indices, prev_indices);
                
                

            }
        }
    }
    
                        }

string map_frags_to_ref(const vector<unique_ptr<Sequence>>& fragments,
                        const unique_ptr<Sequence>& reference,
                        unordered_map<unsigned int, vector<pair<unsigned int, bool>>>& reference_index, int frag_begin, int frag_end) {

    string res = "";
    for (int i = frag_begin; i < frag_end; i++) {
        unordered_map<unsigned int, vector<pair<unsigned int, bool>>> fragment_index;
        make_minimizer_index(fragments[i], fragment_index);
        vector<tuple<bool, int, unsigned int>> match_cluster;
        best_match_cluster(fragment_index, reference_index, match_cluster);
        string paf = "";
        res = paf;
    }             
    
    return res;           

}


string paf_format(const unique_ptr<Sequence>& reference, 
const unique_ptr<Sequence>& fragment, unsigned int target_origin ,bool origin, string* cigar){
    
    string result;
    int tr_broj = 0;
    int ukuZb = 0; 
    int matchZb = 0;
    
    
    
    for (char& i : *cigar){
        
        if (isdigit(i)){
            tr_broj = (tr_broj* 10) + (int(i)- 48); 
        }
        else{
            ukuZb+= tr_broj;
            if (i == 'M'){
                matchZb+=tr_broj;
            }
            tr_broj=0;
        }
    }
    result = "";
    
    result += reference->names.c_str();
    result += '\t';

    result += to_string(reference->all_data.size());
    result += '\t';
    result += '\n';
    
    result += fragment->names.c_str();
    result += '\t';
    result += to_string(fragment->all_data.size());
    result += '\t';
    
    
    result += to_string(target_origin);
    result += '\t';
    if (origin == true) result += "+";
    else result += "-";
    result += '\t';
    

    result += to_string(matchZb);
    result += '\t';
    result += to_string(ukuZb);
    result += '\t';
    result += '\n';
    if (cigar_flag) {
        result += "cg:Z";
        result += *cigar;
    }
    
    result += '\n'; 
    

    //result = "sdfsdf";
    
    return result;
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
    const char *optstring = "m:g:n:a:k:w:f:t:hvco";

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
        case 'o':
            optimize = true;
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

    std::cerr << "Reference genome sequences: \n";
    for (int i = 0; i < int(reference.size()); i++) {
        std::cerr << "\t" << reference[i]->names.c_str() << ", length = " << reference[i]->all_data.size() << "\n";
        cerr << "\n";
    }

    auto frag = bioparser::Parser<Sequence>::Create<bioparser::FastaParser>(argv[argc - 1]);
    vector<unique_ptr<Sequence>> fragments = frag->Parse(-1);
    int frag_size = (int)fragments.size();

    cerr << "Number of fragments: " << fragments.size() << "\n";
    uint64_t len_sum = 0;
    vector<size_t> lengths(fragments.size());
    for (int i = 0; i < int(fragments.size()); i++) {
        lengths[i] = fragments[i]->all_data.size();
        len_sum += lengths[i];
    }
    sort(lengths.begin(), lengths.end(), greater<size_t>());

    uint64_t N50 = -1, tmp_sum = 0;
    for (int i = 0; i < int(fragments.size()); i++) {
        tmp_sum += lengths[i];
        if (tmp_sum * 2 >= len_sum) {
            N50 = lengths[i];
            break;
        }
    }

    cerr << "Average length: " << len_sum * 1.0 / fragments.size() << "\n";
    cerr << "N50 length: " << N50 << "\n";
    cerr << "Minimal length: " << lengths.back() << "\n";
    cerr << "Maximal length: " << lengths.front() << "\n";
    
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
    cout << "\nAligning two random sequences: \n";
    cout << "Target sequence:\n";
    cout << fragments[target_index]->names.c_str() << "\n";
    cout << "Query sequence: \n";
    cout << fragments[query_index]->names.c_str() << "\n";
    int align_score;
    if (optimize && type == GLOBAL) {
        align_score = LinearnaSloz(fragments[query_index]->all_data.c_str(), fragments[query_index]->all_data.size(),
        fragments[target_index]->all_data.c_str(), fragments[target_index]->all_data.size(), match_cost, mismatch_cost, gap_cost, &cigar);
    } else {
        align_score = Align(fragments[query_index]->all_data.c_str(), fragments[query_index]->all_data.size(),
                                                fragments[target_index]->all_data.c_str(),
                                                fragments[target_index]->all_data.size(), type, match_cost, mismatch_cost, gap_cost, &cigar, &target_begin);
    }
    
    cout << "Alignment score: " << align_score << "\n";
    if (cigar_flag) {
        cout << "Cigar: " << cigar << "\n";
    }
    cout << "Target begin: " << target_begin << "\n";




    std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>> reference_index;
    make_reference_index(reference.front(), reference_index);
        
    
    string res = "";
    for (int i = 1; i < 20; i++) {

        std::unordered_map<unsigned int, std::vector<std::pair<unsigned int, bool>>> fragment_index;
        //cout << fragments[i]->all_data.size() << "\n";
        make_minimizer_index(fragments[i], fragment_index);
        vector<tuple<bool, int, unsigned int>> match_cluster;
        
        cout << "Mapping fragments: \n";
        
        best_match_cluster(fragment_index, reference_index, match_cluster);
        
        
        string res = paf_format(reference.front(), fragments[i], target_begin, true, &cigar);
        cout << res;
        //string paf = "";
        //res = paf;
    }

    /*
    thread_pool::ThreadPool thread_pool{};
    vector<std::future<std::string>> ftrs;
    int frag_per_thread = int(fragments.size() / (double)thread_num);
    int frag_begin = 0;
    for (int i = 0; i < (thread_num - 1); i++) {
        ftrs.emplace_back(thread_pool.Submit(map_frags_to_ref, ref(fragments), ref(reference.front()), ref(reference_index), frag_begin, frag_per_thread));
        frag_begin += frag_per_thread;
    }
    ftrs.emplace_back(thread_pool.Submit(map_frags_to_ref, ref(fragments), ref(reference.front()), ref(reference_index), frag_begin, int(fragments.size())));

    for (auto& f : ftrs) {
        string final = f.get();
        cout << final;
    }
    */
    

    
    return 0;
}