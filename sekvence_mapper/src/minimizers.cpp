#include <iostream>
#include <vector>
#include <tuple>
#include <cstring>
#include <list>

using namespace std;

std::vector<std::tuple<unsigned int, unsigned int, bool>>
Minimize(const char *sequence, unsigned int sequence_len, unsigned int kmer_len, unsigned int window_len)
{
    vector<tuple<unsigned int, unsigned int, bool>> v1;
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
    char sub_f[kmer_len + 1];
    char sub_r[kmer_len + 1];
    sub_f[kmer_len] = '\0';
    sub_r[kmer_len] = '\0';
    std::list<string> postojeciMinim;
    for (int i = 0; i < sequence_len - kmer_len + 1; i++)
    {

        bool isp = false;
        int pozicija;
        int f_brojac = 0;
        for (int j = i; j < kmer_len + i; j++)
        {
            sub_f[f_brojac] = sequence[j];
            f_brojac++;
        }

        f_brojac = 0;
        for (int j = sequence_len - kmer_len - i; j < sequence_len - i; j++)
        {
            sub_r[f_brojac] = rev[j];
            f_brojac++;
        }

        char sub2[window_len + 1];
        char min[window_len + 1];
        for (int p = 0; p < window_len; p++)
            min[p] = 'Z';
        min[window_len] = '\0';

        for (int j = 0; j < kmer_len - window_len + 1; j++)
        {

            f_brojac = 0;
            for (int g = j; g < j + window_len; g++)
            {
                sub2[f_brojac] = sub_f[g];
                f_brojac++;
            }
            sub2[window_len] = '\0';

            if (strcmp(sub2, min) < 0)
            {
                for (int g_brojac = 0; g_brojac < window_len; g_brojac++)
                    min[g_brojac] = sub2[g_brojac];
                isp = true;
            }
            min[window_len] = '\0';

            f_brojac = 0;
            for (int g = j; g < j + window_len; g++)
            {
                sub2[f_brojac] = sub_r[g];
                f_brojac++;
            }
            sub2[window_len] = '\0';

            if (strcmp(sub2, min) < 0)
            {
                for (int g_brojac = 0; g_brojac < window_len; g_brojac++)
                    min[g_brojac] = sub2[g_brojac];
                isp = false;
            }
            min[window_len] = '\0';
            pozicija = i + j;
        }
        bool vecPostoji = false;
        for (string s : postojeciMinim)
        {
            if (min == s)
                vecPostoji = true;
        }
        if (vecPostoji == false)
        {
            postojeciMinim.push_back(min);
            unsigned int ispis = atoi(min);
            v1.push_back(tuple<unsigned int, unsigned int, bool>(ispis, pozicija, isp));
        }
    }
    return v1;
}
