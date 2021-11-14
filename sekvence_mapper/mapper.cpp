#include <iostream>
#include <getopt.h>

void display_version() {
    std::cout << "v0.1.0" << "\n";
}

void display_help() {
    std::cout << "Help message" << "\n";
}

int main(int argc, char* argv[]) {
    int option;
    const char* optstring = "h:v";

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

        
        default:
            return 1;
        }
    }
    
    
    
    
    return 0;
}