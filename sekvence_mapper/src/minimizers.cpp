#include <iostream>
#include <vector>
#include <tuple>
#include <cstring>
#include <list>
#include <iomanip>
#include <sstream>

unsigned int get_kmer_v(char* kmer, unsigned int kmer_len) {
    unsigned int value = 0;
    unsigned int i = 0;

    while (i < kmer_len) {
        switch (kmer[i]) {
            case 'A':
                value = (value << 2) + 0;
                break;
            case 'C':
                value = (value << 2) + 1;
                break;
            case 'G':
                value = (value << 2) + 2;
                break;
            case 'T':
                value = (value << 2) + 3;
                break;
        }
        i++;
    }
    return value;
}

using namespace std;

std::vector<std::tuple<unsigned long int, unsigned int, bool>>
Minimize(const char *sequence, unsigned int sequence_len, unsigned int kmer_len, unsigned int window_len)
{
    
    vector<tuple<unsigned long int, unsigned int, bool>> v1;
    vector<tuple<unsigned long int, unsigned int, bool>> v2;
    char rev[sequence_len + 1];
    rev[sequence_len] = '\0';

    vector<tuple<unsigned long int, unsigned int, bool>> minimizers;


    for (int i = sequence_len - 1; i >= 0; i--)
        rev[sequence_len - 1 - i] = sequence[i];

    for (int i = 0; i < sequence_len; i++)
    {
        if (rev[i] == '0')
            rev[i] = '3';
        else if (rev[i] == '3')
            rev[i] = '0';
        if (rev[i] == '1')
            rev[i] = '2';
        else if (rev[i] == '2')
            rev[i] = '1';
    }

    unsigned int fract_len;
    bool origin;

    char* kmer = (char*) malloc(kmer_len);
    unsigned int kmer_value;
    for (unsigned int i = 0; i < window_len - 1; i++) {
        strncpy(kmer, sequence + i, kmer_len);
        kmer_value = get_kmer_v(kmer, kmer_len);
        origin = true;

        if (i == 0 || kmer_value <= get<0>(minimizers.back())) {
            minimizers.push_back(make_tuple(kmer_value, i, origin));
        }
    }

    fract_len = window_len + kmer_len - 1;
    for (unsigned int i = 0; i <= sequence_len - fract_len; i++) {
        if(get<1>(minimizers.back()) >= 1) {
            strncpy(kmer, sequence + i + fract_len - kmer_len, kmer_len);
            kmer_value = get_kmer_v(kmer, kmer_len);
            origin = true;

            if (kmer_value <= get<0>(minimizers.back())) {
                minimizers.push_back(make_tuple(kmer_value, i + fract_len - kmer_len, origin));
            }
        } else {
            unsigned int min_value, min_position;
            bool min_origin;

            for (unsigned int j = 0; j <= fract_len - kmer_len; j++) {
                strncpy(kmer, sequence + i + j, kmer_len);
                kmer_value = get_kmer_v(kmer, kmer_len);
                origin = true;

                if (j == 0 || kmer_value < min_value) {
                    min_value = kmer_value;
                    min_position = j;
                    min_origin = origin;
                }
            }
            minimizers.push_back(make_tuple(min_value, min_position + i, min_origin));
        }
    }

    for (unsigned int i = sequence_len - window_len - kmer_len + 2; i <= sequence_len - kmer_len; i++) {
        if (get<1>(minimizers.back()) < i) {
            unsigned int min_value, min_position;
            bool min_origin;

            for (unsigned int j = i; j <= sequence_len - kmer_len; j++) {
                strncpy(kmer, sequence + j, kmer_len);
                kmer_value = get_kmer_v(kmer, kmer_len);
                origin = true;

                if (j == i || kmer_value < min_value) {
                    min_value = kmer_value;
                    min_position = j;
                    min_origin = origin;
                }
            }
            minimizers.push_back(make_tuple(min_value, min_position, min_origin));
        }


    }
    free(kmer);
    return minimizers;
}
/*
int main()
{
    unsigned int kmer_len = 15;
    unsigned int window_len = 19;
    vector<tuple<unsigned long int, unsigned int, bool>> v1 =
    Minimize("TACGCATAAGCGCCAAAAGCACGCCGGGCGACCATAATGAGAAGATGCTCACCGCCAGCGTCAACAGTAAGTAGACGATGAGAGCGCCAGGGTAAACAGTTTGCACGGGCCGCGACGAAGGCAGTCGATAATGCCGCCATTGCTGGGATGCTCGCCCCCAGACGCGCATAGGCATAACCGGAAAACATCGCCACAATACCGCCAAAGCAAAGGCGACCCCAGGTCGAGGCTTCCATTAGCAATGCAGCCTGCCCCAGCAGCGCGAAGATCGCCGCCCCCACCATTGCCCCGGATTATCGATGGAAACGACGTTCCAGACCGCGGGTTTGTTACCGTTATTACCTTCCGTGTTTCATCATGTAATAGCAGCCCTTAGTAAACACGTTATAGCCGAAAAATTGCTTAGCGACCGATGCCACCATTTTTGACCGCGAACACTGTTGCCATCTTCGTCCAGCGGTGGTGAGAGAACGCGGCAATTCCCATCACTCCAGGGGACGACCACGCCAGAATACCGCCACCTACACCCGGTTTACGCCCGGTAAACCAACACGATACGCCCAGTCACCGGGAACATACCCAAGACTACAGCCCTTCCATCATCATTTCGTAAGAATGTACGGCACGTTGTCGGCCTGAAGAACGCGTTCCTGCGTCAACGATTCGCACACCACCTGCGCCGCGCTTGTCGCGCCAAGCGTTGCCAGTTCGCTACGTGGTCGAGGAGCCCTTGGAGCACTGACGGGTATACACGTCACAGGCTTCCCATTGCATCACAATAGAGATATCCGGCGGAGTAAGCAGCCAGGCTATGGCCCGGTTATGGAGACTTTATTTGTTTGTTCCGACTGGTTAACTTCGTCAGAGAGCGCTACCTGCTCGCCAGCAGTTTGCTTTGGATATGTAAATTCGCTGCCAGCGTTGTTCAACAGTTTTCAGCGTTAATCAGGCTATCTTAGGCAATAGCGCCAGGATAGATTCACCAGTGGCGGGAAATCGTTTGCCGCGGATGCAACTCTAAGGCGATAACTGAGTTACGGGGCAAGATCCGGTCGGGTCCCCACTAGGGATCTTTGTGTCGTACCGCCTGCGGGCCGCACATCTTCCAACGCAAGGGCTAACGTACAGAGACTTTCGAGATGGATTCCAGTGCAAAGCGTAATCACTGTCACCGCACTACGACGTTGCCATCGCGCAGTACACGATAGCCACTGCCGCCAGTTGACCTGGTACATTCGCCGAAGGGAATGTAATCGGGCATTTTGTCCGCCGTTAAGTGAGTGAGAAATTGGGTGTAAGCCTGATCTGTTACGGAATTGCCTGCTGTAATTTGTTTGCATCTAACAATATTTTCTTGTTAACTCCTTTTATAAGTCTCGGGAGGTAATTCCTCACCGGCTGGTGCCGATTTCAGGCATCCTGATTTAACTTAGCACCGGAGACTTAACTACAGGAAAACACAAAGAGATAAATGTCTAATCCTGATGCAAAAACTGCATCAGCAAATTTTTAATCTTTACGGACTTTTACCGCCTGGTTTATTAATTTCTTGACCTTCCCCTTGCTGAAGGTTTAACCTTTATCACAGCCAGTCAAAACCGTGTGTAAAGGGATGTTTATGTCAAACTATCGACCTGACCCTCTGGACGGCCTGTCTGCGGCGGTCACTGCGTTAAACGCGTGAAAGAAAGTCTTGAACAGCGTCGTCCGGGATGTTGAGCAGGCGGATGTGTCTATCACTGAAGCACGTTACCGGGACTGCCGAACAGGCTCGCTAAATTGAAACCATCAAACAAGCGGGTTATGACGCATCTGTAAGCCACCCAAAGGCTAAACCGCCACAGCCGTGACTGTTATGTT", 
    1860, kmer_len, window_len);
    for (int iter = 0; iter<v1.size(); iter++){
        //for(int iter2 = 0; iter2 < 2;iter2++){
            
            cout << setw(kmer_len) << setfill('0') << get<0>(v1[iter]) << " ";
            cout << get<1>(v1[iter]) << " ";
            cout << get<2>(v1[iter]) << " ";
            cout << endl;
        //}
    }
}
*/