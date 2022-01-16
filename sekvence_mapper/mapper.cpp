#include <iostream>
#include <getopt.h>
#include <bits/stdc++.h>

/*

    int main()
    {
        
        unsigned long int sequence = 231032101233101;
        
        int n = countDigit(sequence);
        
        cout << "Za broj " << sequence << " najveći rastući niz je: " ;
        
        longestSubsequence(sequence,n);
        
    }

*/

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

void display_version()
{
    std::cout << "v0.1.0"
              << "\n";
}

void display_help()
{
    std::cout << "Help message"
              << "\n";
}

int main(int argc, char *argv[])
{
    int option;
    const char *optstring = ":hv";
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

        default:
            exit(1);
        }
    }

    if (optind < argc)
    {
        std::cout << argv[optind++] << "\n";
        std::cout << argv[optind] << "\n";
    }

    return 0;
}