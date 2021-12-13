#include <iostream>
#include <vector>
#include <tuple>

using namespace std;

std::vector<std::tuple<unsigned int, unsigned int, bool>>
Minimize(const char *sequence, unsigned int sequence_len, unsigned int kmer_len, unsigned int window_len)
{

    // od reda 15 do 124 nalazi se nedovršeni kod u c++
    // kod je preveden do kraja ali još treba pronaći neke greške
    // na dnu datoteke nalazi se python kod

    char seq[101] = "ATGCGATATCGTAGGCGTCGATGGAGAGCTAGATCGATCGATCTAAATCCCGATCGATTCCGAGCGCGATCAAAGCGCGATAGGCTAGCTAAAGCTAGCA";
    char rev[101];
    for (int i = 99; i >= 0; i--)
    {
        rev[99 - i] = seq[i];
    }
    for (int i = 0; i < 100; i++)
    {
        if (rev[i] == 'A')
            rev[i] = 'X';
    }
    for (int i = 0; i < 100; i++)
    {
        if (rev[i] == 'T')
            rev[i] = 'A';
    }
    for (int i = 0; i < 100; i++)
    {
        if (rev[i] == 'X')
            rev[i] = 'T';
    }
    for (int i = 0; i < 100; i++)
    {
        if (rev[i] == 'C')
            rev[i] = 'X';
    }
    for (int i = 0; i < 100; i++)
    {
        if (rev[i] == 'G')
            rev[i] = 'C';
    }
    for (int i = 0; i < 100; i++)
    {
        if (rev[i] == 'X')
            rev[i] = 'G';
    }

    int kmer = 31;
    int M = 7;
    int L = 100;

    for (int i = 0; i < L; i++)
    {
        cout << seq[i];
    }
    cout << endl;
    for (int i = 0; i < L; i++)
    {
        cout << rev[i];
    }
    cout << endl;

    int gg = 0;

    char sub_f[kmer + 1];
    char sub_r[kmer + 1];
    sub_f[kmer] = '\0';
    sub_r[kmer] = '\0';

    for (int i = 0; i < L - kmer + 1; i++)
    {

        int f_brojac = 0;
        for (int j = i; j < kmer + i; j++)
        {
            sub_f[f_brojac] = seq[j];
            f_brojac++;
        }
        f_brojac = 0;
        for (int j = L - kmer - i; j < L - i; j++)
        {
            sub_r[f_brojac] = rev[j];
            f_brojac++;
        }

        char min[14] = "ZZZZZZZZZZZZZ";
        char sub2[M + 1];
        sub2[M] = '\0';
        f_brojac = 0;
        for (int j = 0; j < kmer - M + 1; j++)
        {

            f_brojac = 0;
            for (int g = j; g < j + M; g++)
            {
                sub2[f_brojac] = sub_f[g];
                f_brojac++;
            }
            if (strcmp(sub2, min) < 0)
            {
                for (int fff = 0; fff < f_brojac; fff++)
                    min[fff] = sub2[fff];
                for (int fff = f_brojac; fff < 14; fff++)
                    min[fff] = '\0';
            }
            for (int g_brojac = j; g_brojac < j + M; g_brojac++)
            {
                sub2[g_brojac] = sub_r[g_brojac];
            }
            if (strcmp(sub2, min) < 0)
            {
                for (int fff = 0; fff < f_brojac; fff++)
                    min[fff] = sub2[fff];
                for (int fff = f_brojac; fff < 14; fff++)
                    min[fff] = '\0';
            }
        }

        cout << sub_f << " " << min << endl;
    }

    // kod u pythonu

    /** 
                seq="ATGCGATATCGTAGGCGTCGATGGAGAGCTAGATCGATCGATCTAAATCCCGATCGATTCCGAGCGCGATCAAAGCGCGATAGGCTAGCTAAAGCTAGCA"

                rev=seq[::-1]

                rev=rev.replace("A","X")
                rev=rev.replace("T","A")
                rev=rev.replace("X","T")
                rev=rev.replace("C","X")
                rev=rev.replace("G","C")
                rev=rev.replace("X","G")

                Kmer=31
                M=7
                L=len(seq)

                for i in range(0, L-Kmer+1):

                        sub_f=seq[i:i+Kmer]
                        sub_r=rev[L-Kmer-i:L-i]

                        min="ZZZZZZZZZZZZZ"
                        for j in range(0, Kmer-M+1):
                                sub2=sub_f[j:j+M]
                                if sub2 < min:
                                        min=sub2
                                sub2=sub_r[j:j+M]
                                if sub2 < min:
                                        min=sub2

                        print (sub_f,min)  **/
}
