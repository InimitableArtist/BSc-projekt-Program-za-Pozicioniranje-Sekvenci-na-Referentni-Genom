#include <iostream>
#include "alignment.h"
#include <algorithm>
#include <initializer_list>
#include <string>
#include <math.h>
	
using namespace std;



std::pair<int, ParentTrack> find_max(int values[], int length, bool mtch) {
    int max = values[0];
    int index = 0;
    ParentTrack track;
    for (int i = 0; i < length; i++) {
        if (values[i] > max) {
            max = values[i];
            index = i;
        }
    }
    if (index == 0) {
        if (mtch) {
            track = MATCH;
        }
        else {
            track = MISMATCH;
        }
    }
    else if (index == 1) track = INSERTION;
    else if (index == 2) track = DELETION;
    else track = NONE;
    return std::make_pair(max, track);
}

std::pair<int, int> get_target_begin(int res_row,int res_column,struct Cell **V){
    ParentTrack now = V[res_row][res_column].parent;
    while(now != NONE){
        if(now == MATCH || now == MISMATCH) {
            res_row--;
            res_column--;
        }
        else if(now == INSERTION) res_row--;    
        else res_column--;
        
        now = V[res_row][res_column].parent;
    }
    return std::make_pair(res_row, res_column);
}

std::string get_cigar(int res_row,int res_column, Cell **V){
    string str = "";
    ParentTrack now = V[res_row][res_column].parent;
    ParentTrack prev = NONE;
    int dist=1;

    while(now != NONE){
        if(now == MATCH || now == MISMATCH) {
            res_row--;
            res_column--;
        }
        else if(now == INSERTION) res_row--;    
        else res_column--;
        prev = now;
        now = V[res_row][res_column].parent;
        
        if(prev == now) dist++;
        else{
            
            if (prev == MATCH) str.insert(0, "M");
            else if (prev == MISMATCH) str.insert(0, "X");
            else if (prev == INSERTION) str.insert(0, "I");
            else str.insert(0, "D");
            str.insert(0, std::to_string(dist));
            dist = 1;
        }
    }
    return str;
}



int Align(const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    AlignmentType type,
    int match,
    int mismatch,
    int gap,
    std::string* cigar = nullptr,
    unsigned int* target_begin = nullptr) {   
    
    //Dynamic programming table
    Cell **V = new Cell*[query_len + 1];
    for (unsigned int i = 0; i <= query_len; i++) {
        V[i] = new Cell[target_len + 1];
    }
    int align = 0;

    //Table coords of the result
    int res_row = 0;
    int res_column = 0;
    int match_cost = 0;

    //Smith-Waterman
    if (type == LOCAL) {
        int MAX = 0;
        int left, up;
        V[0][0].cost = 0;
        V[0][0].parent = NONE;

        for (int i = 0; i <= query_len; i++) {
            V[i][0].cost = 0;
            V[i][0].parent = NONE;
        }
        for (int j = 0; j <= target_len; j++) {
            V[0][j].cost = 0;
            V[0][j].parent = NONE;
        }
        for (int i = 1; i <= query_len; i++) {
            for (int j = 1; j <= target_len; j++) {
                bool mtch = false;
                if (query[i - 1] == target[j - 1]) {
                    match_cost = V[i - 1][j - 1].cost + match;
                    mtch = true;
                } 
                else match_cost = V[i - 1][j - 1].cost + mismatch;

                left = V[i][j - 1].cost + gap;
                up = V[i - 1][j].cost + gap;

                int values[4] = {match_cost, up, left, 0};
                std::pair<int, ParentTrack> res = find_max(values, 4, mtch);
                V[i][j].cost = res.first;
                V[i][j].parent = res.second;

                //Update MAX coords
                if (V[i][j].cost > MAX) {
                    MAX = V[i][j].cost;
                    res_row = i;
                    res_column = j;
                }


            }
        }
        align = V[res_row][res_column].cost;

    }
    //Needleman-Wunsch
    else if (type == GLOBAL) {
        V[0][0].cost = 0;
        V[0][0].parent = NONE;
        int poravnanje = 0;
        for (int i = 1; i <= query_len; i++) {
            V[i][0].cost = V[i-1][0].cost + gap;
            V[i][0].parent = INSERTION;
        }
        for (int j = 1; j <= target_len; j++) {
            V[0][j].cost = V[0][j-1].cost + gap;
            V[0][j].parent = DELETION;
        }
        for (int i = 1; i <= query_len; i++) {
            for (int j = 1; j <= target_len; j++) {
                bool mtch = false;
                
                if (query[i-1] == target[j-1]) {
                    poravnanje = match;
                    mtch = true;
                }
                else{
                    poravnanje = mismatch;
                }

                int var1 = V[i-1][j-1].cost + poravnanje; //M
                int var2 = V[i-1][j].cost + gap; //I
                int var3 = V[i][j-1].cost + gap; //D

                int vars[3] = {var1, var2, var3};
                std::pair<int, ParentTrack> res = find_max(vars, 3, mtch);
                
                V[i][j].cost = res.first;
                V[i][j].parent = res.second;
            }
        }
        res_row = query_len;
        res_column = target_len;
        align = V[res_row][res_column].cost;
    }

    else if (type == SEMI_GLOBAL) {
        V[0][0].cost = 0;
        V[0][0].parent = NONE;
        int MAX = INT32_MIN;
        int left, up;

        for (int i = 1; i <= query_len; i++) {
            V[i][0].cost = 0;
            V[i][0].parent = NONE;
        }
        for (int j = 1; j <= target_len; j++) {
            V[0][j].cost = 0;
            V[0][j].parent = NONE;
        }

        for (int i = 1; i <= query_len; i++) {
            for (int j = 1; j <= target_len; j++) {
                bool mtch = false;
                if(query[i - 1] == target[j - 1]) {
                    match_cost = V[i - 1][j - 1].cost + match; 
                    mtch = true;
                } else {
                    match_cost = V[i - 1][j - 1].cost + mismatch;
                }
                left = V[i][j - 1].cost + gap;
                up = V[i - 1][j].cost + gap;

                int values[4] = {match_cost, up, left}; 
                std::pair<int, ParentTrack> res = find_max(values, 3, mtch);
                V[i][j].cost = res.first;
                V[i][j].parent = res.second;

                if (V[i][j].cost > MAX && ((i == query_len) || (j == target_len))) {
                    MAX = V[i][j].cost;
                    res_row = i;
                    res_column = j;
                }
            }
        }
        
	align = V[res_row][res_column].cost;
    }

    std::pair<int, int> coords = get_target_begin(res_row, res_column, V);
    if (target_begin != nullptr) {
        if (type == GLOBAL) {
            *target_begin = 0;
        }
        else if (type == SEMI_GLOBAL) {
            *target_begin = coords.second;
        }
        else {
            *target_begin = coords.second;
        }
    }
    *cigar = get_cigar(res_row, res_column, V);

    return align;   
}




int *scoreNW(const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    int match,
    int mismatch,
    int gap){
    
    //Dynamic programming table
    int V[2][query_len+1];
    V[0][0]= 0;
    int poravnanje = 0;
    for(int j=1;j<=query_len;j++){
        V[0][j] = V[0][j-1] + gap;
    }
    for(int i=1;i<=target_len;i++){
        V[1][0] = V[0][0] + gap;
        
        for (int j = 1; j <= query_len; j++) {
                
            if (query[j-1] == target[i-1]) {
                poravnanje = match;
            }
            else{
                poravnanje = mismatch;
            }
            int var1 = V[0][j-1] + poravnanje; //M
            int var2 = V[0][j] + gap; //I
            int var3 = V[1][j-1] + gap; //D

            V[1][j] = max({var1,var2,var3});
            //cout<<V[1][j];
            
        }
        //cout<<"\n";
        for(int j=0;j<=query_len;j++){
            V[0][j]=V[1][j];
        }
    }
    int *zadnjiRed =(int*) malloc(sizeof(int) * (query_len+1));
    for(int j=0;j<=query_len;j++){
        zadnjiRed[j] = V[1][j];
        
        
    }
    
    return zadnjiRed;
}


int LinearnaSloz(const char* query, unsigned int query_len,
    const char* target, unsigned int target_len,
    int match,
    int mismatch,
    int gap,
    std::string* cigar = nullptr) {
    
    
    
    if((target_len<=1) or (query_len<=1)){
        
        string* cigarA = new string;
        int score =Align(query,query_len,
        target, target_len,
        GLOBAL,
        match,
        mismatch,
        gap, 
        cigarA);
        
        string str = *cigarA;
        string str2 = *cigar;
        
        int i= 0;
        int br=0;
        while (isdigit(str[i])){
        	br=br + (((int)(str[i])-48) * (pow(10,i)));
        	      	
        	i++;
        }
        if(!(str2.empty()) and (str[i]==(str2).back())){
        	str = str.substr(i);
        	
        	int i = 2;       		
        	int br2 =0;
        			
        	while (isdigit(str2[str2.length()-i])){
        		
        		br2=((int)(str2[str2.length()-i])-48) + br2 * pow(10,i );
        		i++;	
        	}
        		
        	
        	str2 = str2.substr(0, str2.length()-(i-1));
        	*cigar = str2;
        	*cigar = *cigar + to_string(br + br2);
        	*cigar = *cigar + str;
        }
        else{
        	*cigar = *cigar + str;
        }
        return score;
    }
    
    
    else{
        int half = target_len/2;
        string queryS = query;
        string targetS = target; 
        
        string tr = targetS.substr(0,half);
        string trR = targetS.substr(half, target_len);
        reverse(trR.begin(), trR.end());
        
        string qrR = queryS;
        reverse(qrR.begin(), qrR.end());
        
        int *zr = scoreNW(query, query_len, tr.c_str(), half, match, mismatch, gap);
        
        if (target_len%2!=0){
            half++;
            
        }
  
        int *zrR = scoreNW(qrR.c_str(),query_len, trR.c_str(), half, match, mismatch, gap);
        
        if (target_len%2!=0){
            half--;    
        }
        
        int mx;
        int poz;
        for (int i=0;i<query_len;i++){
            if (i==0){
                mx = *(zr) + *(zrR+(query_len));
                
                poz = 0;
            }
            else{
                if(mx < (*(zr+i) + *(zrR+(query_len)-i))){
                    mx = *(zr+i) + *(zrR+(query_len)-i);
                    poz = i;
                    
                } 
            }
        }
        free(zr);
        free(zrR);
        
        string qr1 = queryS.substr(0,poz);
        string qr2 = queryS.substr(poz, query_len-poz);
        
        
        reverse(trR.begin(), trR.end());
        
        int d1 = LinearnaSloz(qr1.c_str(),poz, tr.c_str(), half,match, mismatch, gap, cigar);
        
        if (target_len%2!=0){
            
            half++;
        }
        int d2 = LinearnaSloz(qr2.c_str(),query_len-poz, trR.c_str(), half,match,mismatch,gap, cigar);
           
        return (d1+d2);
    }
    
}	



