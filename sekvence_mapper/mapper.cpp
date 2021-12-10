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
    const char* optstring = ":hv";
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
            exit(1);
        }
        
    }

    if (optind < argc) {
        std::cout << argv[optind++] << "\n";
        std::cout << argv[optind] << "\n";
    }
    
    
    
    
    return 0;
}