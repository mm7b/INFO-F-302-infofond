#include <iostream>
#include <string>
#include <sstream>
#include <utility>
#include <exception>
#include <stdexcept>
#include <cctype>
#include <stdlib.h>
#include "Solver.hpp"

template<typename T>
std::string to_string(const T& value){
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

class ParseException : public std::logic_error {
public:
    ParseException(const std::string& msg) : std::logic_error(msg) {}
};

struct Rectangle {
    int length, width, height, index;

    Rectangle() : length(-1), width(-1), height(-1), index(-1) {}
    Rectangle(int n, int m, int h, int i) : length(m), width(n), height(h), index(i) {}

    void print() { std::cout << index << " " << length << " " << width << (height != -1 ? to_string(height) : ""); }
};

struct Square : public Rectangle{
    Square() : Rectangle() {}
    Square(int n, int i) : Rectangle(n, n, n, i) {}
};

struct OrthogonalPacking{
    int k;
    int dim;
    Rectangle* container_rect;
    Rectangle** rects;

    OrthogonalPacking(int parsed_k) : k(parsed_k) {
        container_rect = NULL;
        rects = new Rectangle*[k];
        for(int i = 0; i < k; ++i){
            rects[i] = NULL;
        }
    }

    OrthogonalPacking(const OrthogonalPacking& other) : k(other.k), dim(other.dim) {
        container_rect = new Rectangle(*(other.container_rect));
        rects = new Rectangle*[k];
        for(int i = 0; i < k; ++i){
            rects[i] = (other.rects[i] == NULL) ? NULL : new Rectangle(*(other.rects[i]));
        }
    }

    OrthogonalPacking& operator=(const OrthogonalPacking& other){
        if(this != &other){
            k = other.k; dim = other.dim;

            delete container_rect;
            container_rect = new Rectangle(*(other.container_rect));
            
            for(int i = 0; i < k; ++i){
                if(rects[i] != NULL) { delete rects[i]; }
            }
            delete[] rects;

            rects = new Rectangle*[k];
            for(int i = 0; i < k; ++i){
                rects[i] = (other.rects[i] == NULL) ? NULL : new Rectangle(*(other.rects[i]));
            }
        }
        return *this;
    }

    ~OrthogonalPacking(){
        if(container_rect != NULL) { delete container_rect; }
        for(int i = 0; i < k; ++i){
            if(rects[i] != NULL) { delete rects[i]; }
        }
        delete[] rects;
    }

};

int next_int(std::string& line){
    std::size_t i = 0;
    std::string digit = "";
    while(i < line.length() && isdigit(line[i])){
        digit += line[i];
        ++i; 
    }
    if(i == 0){ throw ParseException("Could not read int from : " + line); }

    line = (++i < line.length()) ? line.substr(i, line.length() - i) : "";
    return atoi(digit.c_str());
}

OrthogonalPacking input_reader(bool three_dim = false){
    std::string input_line;
    std::getline(std::cin, input_line);
    int k = next_int(input_line);
    std::getline(std::cin, input_line);
    int n = next_int(input_line);
    std::getline(std::cin, input_line);
    int m = next_int(input_line);
    int h = -1;
    if(three_dim){
        std::getline(std::cin, input_line);
        h = next_int(input_line);
    }
    OrthogonalPacking data(k);
    data.container_rect = new Rectangle(n, m, h, 0);
    for(int i = 0; i < k; ++i){
        std::getline(std::cin, input_line);
        int index = next_int(input_line);
        if(index != i + 1) { throw ParseException("Wrong index : expected " + to_string(i + 1) + " found " + to_string(index) + " instead"); }
        int m = next_int(input_line);
        int n = next_int(input_line);   
        int h = three_dim ? next_int(input_line) : -1;
        data.rects[i] = new Rectangle(n, m, h, index);
    }
    return data;
}


int main() {
    try{
        OrthogonalPacking data(input_reader());
        data.container_rect->print();
        for(int i = 0 ; i < data.k ; ++i){
            data.rects[i]->print();
        }
        return 0;
    }catch(const std::exception& e){
        std::cout << e.what() << std::endl;
        return 1;
    }
}

/*
int main() {
    Solver s;
    vec<Lit> lits;

    int X, Y, K;
    istream & input = cin;
    next_int(input, K, "Nombre de rectangles");
    next_int(input, X, "Largeur");
    next_int(input, Y, "Longueur");

    FOR(r, 1, K){

    }


    // exemple si le R est de taille 6-7-1 (2D quoi)
    int prop[X][Y][1];

    FOR(x,1,X-1){
        FOR(y,1,Y-1){
            FOR(z, 1,1){
                prop[x][y][z] = s.newVar();
            }
        }
    }

    //ajout des contraintes


    
    // call the SAT solver
    s.solve();

    // récupération de la solution
    while (s.okay()) {
        lits.clear();
        cout << "Nouvelle solution: \n\n" ;
        s.addClause(lits);
        s.solve();
    }
    cout << "Il n'y a plus de solutions\n";
}
*/
