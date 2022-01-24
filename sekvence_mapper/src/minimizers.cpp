#include <iostream>
#include <vector>
#include <tuple>
#include <cstring>
#include <list>
#include <iomanip>
#include <sstream>

using namespace std;

std::vector<std::tuple<unsigned long int, unsigned int, bool>>
Minimize(const char *sequence, unsigned int sequence_len, unsigned int kmer_len, unsigned int window_len)
{
    vector<tuple<unsigned long int, unsigned int, bool>> v1;
    vector<tuple<unsigned long int, unsigned int, bool>> v2;
    char rev[sequence_len + 1];
    rev[sequence_len] = '\0';

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

    char sub_f[window_len + 1];
    char sub_r[window_len + 1];
    sub_f[window_len] = '\0';
    sub_r[window_len] = '\0';
    std::list<string> postojeciMinim;
    std::list<string> postojeciMinimRev;
    for (int i = 0; i < sequence_len - window_len + 1; i++)
    {

        bool isp = false;
        bool isp2 = false;
        int f_brojac = 0;
        for (int j = i; j < window_len + i; j++)
        {
            sub_f[f_brojac] = sequence[j];
            f_brojac++;
        }

        f_brojac = 0;
        for (int j = sequence_len - window_len - i; j < sequence_len - i; j++)
        {
            sub_r[f_brojac] = rev[j];
            f_brojac++;
        }
        //cout << sub_f << " " << sub_r << endl;

        char sub2[kmer_len + 1];
        char sub3[kmer_len + 1];
        char min[kmer_len + 1];
        char minrev[kmer_len + 1];

        for (int p = 0; p < kmer_len; p++)
        {
            min[p] = 'Z';
            minrev[p] = 'Z';
        }
        min[kmer_len] = '\0';
        minrev[kmer_len] = '\0';

        for (int j = 0; j < window_len - kmer_len + 1; j++)
        {
            //sub2
            f_brojac = 0;
            for (int g = j; g < j + kmer_len; g++)
            {
                sub2[f_brojac] = sub_f[g];
                f_brojac++;
            }
            //sub3
            f_brojac = 0;
            for (int g = j; g < j + kmer_len; g++)
            {
                sub3[f_brojac] = sub_r[g];
                f_brojac++;
            }
            sub2[kmer_len] = '\0';
            sub3[kmer_len] = '\0';

            if (strcmp(sub2, min) < 0)
            {
                for (int g_brojac = 0; g_brojac < kmer_len; g_brojac++)
                    min[g_brojac] = sub2[g_brojac];
                isp = true;
            }
            min[kmer_len] = '\0';

            if (strcmp(sub3, minrev) < 0)
            {
                for (int g_brojac = 0; g_brojac < kmer_len; g_brojac++)
                    minrev[g_brojac] = sub3[g_brojac];
                isp2 = true;
            }
            minrev[kmer_len] = '\0';
        }
        bool vecPostoji = false;
        bool vecPostojiRev = false;
        for (string s : postojeciMinim)
        {
            if (min == s)
                vecPostoji = true;
        }
        for (string s : postojeciMinimRev)
        {
            if (minrev == s)
                vecPostojiRev = true;
        }
        if (vecPostoji == false)
        {
            postojeciMinim.push_back(min);
            unsigned long int ispis = 0;
            for (int i = 0; i < kmer_len; i++)
            {
                if (min[i] == 'A')
                {
                    ispis = ispis * 10 + 0;
                }
                else if (min[i] == 'C')
                {
                    ispis = ispis * 10 + 1;
                }
                else if (min[i] == 'G')
                {
                    ispis = ispis * 10 + 2;
                }
                else if (min[i] == 'T')
                {
                    ispis = ispis * 10 + 3;
                }
            }

            //cout << setw(kmer_len) << setfill('0') << ispis <<endl;

            std::stringstream ss;
            ss << sub_f;
            string str1 = ss.str();
            std::stringstream ss2;
            ss2 << min;
            string str2 = ss2.str();

            int position = 0;
            int index_str;
            while ((index_str = str1.find(str2, position)) != string::npos)
            {
                position = index_str + 1;
            }
            int pozicija = i + position;
            //cout << min << " " << ispis << endl;
            if (isp == 1)
                v1.push_back(tuple<unsigned long int, unsigned int, bool>(ispis, pozicija, isp));
        }
        if (vecPostojiRev == false)
        {
            postojeciMinimRev.push_back(minrev);
            unsigned long int ispis = 0;
            for (int i = 0; i < kmer_len; i++)
            {
                if (minrev[i] == 'A')
                {
                    ispis = ispis * 10 + 0;
                }
                else if (minrev[i] == 'C')
                {
                    ispis = ispis * 10 + 1;
                }
                else if (minrev[i] == 'G')
                {
                    ispis = ispis * 10 + 2;
                }
                else if (minrev[i] == 'T')
                {
                    ispis = ispis * 10 + 3;
                }
            }
            std::stringstream ss;
            ss << sub_r;
            string str1 = ss.str();
            std::stringstream ss2;
            ss2 << minrev;
            string str2 = ss2.str();
            int position = 0;
            int index_str;
            while ((index_str = str1.find(str2, position)) != string::npos)
            {
                position = index_str + 1;
            }
            position = position + kmer_len - 1;
            position = window_len - position;
            int pozicija = i + position + 1;

            //cout << min << " " << ispis << endl;
            if (isp2 == 1)
                v2.push_back(tuple<unsigned long int, unsigned int, bool>(ispis, pozicija, !isp2));
        }
    }
    v1.insert(v1.end(), v2.begin(), v2.end());
    return v1;
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