#include <iostream>
#include <string>
#include <sstream>
#include <exception>
#include <stdexcept>
#include <cctype>
#include <stdlib.h>
#include <algorithm>

#include "Solver.hpp"

template<typename T>
std::string to_string(const T& value){
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

struct OrthogonalPackingSolver : public Solver {

    OrthogonalPackingProblem problem;

    OrthogonalPackingSolver(const OrthogonalPackingProblem& p) : Solver(), problem(p) {
        add_constraints();

    }

    void add_constraints(){
        vec<Lit> lits;
        int mu[problem.k][problem.n][problem.m];

        /* On ne peut avoir 2 mu à vrai en même temps pour un même k, 
        un rectangle k ne pouvant avoir qu'une seule position */
        for(int k = 0; k < problem.k; ++k){
            for(int n = 0; n < problem.n; ++n){
                for(int m = 0; m < problem.m; ++m){
                    for(int i = n + 1; i < problem.n; ++i){
                        for(int j = m + 1; j < problem.m; ++j){
                            addBinary(~Lit(mu[k][n][m]), ~Lit(mu[k][i][j]));
                        }
                    }
                }
            }
        }

        /* On doit avoir au mois 1 mu à vrai pour un rectangle k donné  */
        for(int k = 0; k < problem.k; ++k){
            lits.clear();
             for(int n = 0; n < problem.n; ++n){
                for(int m = 0; m < problem.m; ++m){
                    lits.push(mu[k][n][m]);
                }
            }
            addClause(lits);
        }

        /* Pas de superposition : un rectangle est soit à gauche, soit à droite, 
        soit en haut, soit en bas, soit plus profond, soit moins profond que tous les autres */
    }

    void print_solution(std::wstream& out = std::cout){
        out << "slt sa va??";
    }
};

struct OrthogonalPackingProblem{
    int k, dim, n, m, h;
    int* lengths, widths, heights;

    OrthogonalPackingProblem(int parsed_k, int parsed_dim, int parsed_n, int parsed_m, int parsed_h) : 
        k(parsed_k), dim(parsed_dim), n(parsed_n), m(parsed_m), h(parsed_h),
        lengths(new int[k]), widths(new int[k]), heights(new int[k]) {}

    OrthogonalPackingProblem(const OrthogonalPackingProblem& other) : 
        k(other.k), dim(other.dim), n(other.n), m(other.m), h(parsed_h),
        lengths(new int[k]), widths(new int[k]), heights(new int[k]) {
            std::copy(other.lengths, other.lengths + k, lengths);
            std::copy(other.widths, other.widths + k, widths);
            std::copy(other.heights, other.heights + k, heights);
    }

    OrthogonalPackingProblem& operator=(const OrthogonalPackingProblem& other){
        if(this != &other){
            k = other.k; dim = other.dim; n = other.n; m = other.m; h = other.h;
            delete[] lengths; delete[] widths; delete[] heights;
            lengths = new int[k]; widths = new int[k]; heights = new int[k];
            std::copy(other.lengths, other.lengths + k, lengths);
            std::copy(other.widths, other.widths + k, widths);
            std::copy(other.heights, other.heights + k, heights);
        }
        return *this;
    }

    ~OrthogonalPackingProblem(){
        delete[] lengths;
        delete[] widths;
        delete[] heights;
    }


    struct Parser{
        struct ParseException : public std::runtime_error {
            ParseException(const std::string& msg) : std::runtime_error(msg) {}
        };

        inline static int next_int(std::string& line){
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

        inline static OrthogonalPackingProblem parse(std::istream& in, bool three_dim = false){
            int dim = three_dim ? 3 : 2;
            std::string input_line;
            std::getline(in, input_line);
            int k = next_int(input_line);
            std::getline(in, input_line);
            int n = next_int(input_line);
            std::getline(in, input_line);
            int m = next_int(input_line);
            int h = -1;
            if(three_dim){
                std::getline(in, input_line);
                h = next_int(input_line);
            }
            OrthogonalPackingProblem problem(k, dim, n, m, h);
            for(int i = 0; i < k; ++i){
                std::getline(in, input_line);
                int index = next_int(input_line);
                if(index != i + 1) { throw ParseException("Wrong index : expected " + to_string(i + 1) + " found " + to_string(index) + " instead"); }
                problem.lengths[i] = next_int(input_line);
                problem.widths[i] = next_int(input_line);   
                problem.heights[i] = three_dim ? next_int(input_line) : -1;
            }
            return problem;
        }
    };
};

int main() {
    try{
        OrthogonalPackingSolver solver(OrthogonalPackingProblem::Parser::parse(std::cin));
        solver.solve();
        solver.print_solution();
        return 0;
    } catch(const std::exception& e){
        std::cout << e.what() << std::endl;
        return 1;
    }
}

