#include <string>
#include <sstream>
#include <stdlib.h>
#include <exception>
#include <stdexcept>

#include "Solver.hpp"

template<typename T>
std::string to_string(const T& value){
    std::ostringstream oss;
    oss << value;
    return oss.str();
}


struct OrthogonalPackingProblem{
    int k, dim, n, m, h;
    int* lengths; int* widths; int* heights;

    OrthogonalPackingProblem(int parsed_k, int parsed_dim, int parsed_n, int parsed_m, int parsed_h) : 
        k(parsed_k), dim(parsed_dim), n(parsed_n), m(parsed_m), h(parsed_h),
        lengths(new int[k]), widths(new int[k]), heights(new int[k]) {}

    OrthogonalPackingProblem(const OrthogonalPackingProblem& other) : 
        k(other.k), dim(other.dim), n(other.n), m(other.m), h(other.h),
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

    void print(std::ostream& out = std::cout){
        out << "k(" + to_string(k) + ") " << "dim(" + to_string(dim) + ") "
            << "n(" + to_string(n) + ") " << "m(" + to_string(m) + ") "
            << (h > 0 ? "h(" + to_string(h) + ") " : "") << std::endl
            << "lengths(";
        for(int i = 0; i < k; ++i){ out << (i == k - 1 ? to_string(lengths[i]) + ")" : to_string(lengths[i]) + ", "); }
        out << std::endl << "widths(";
        for(int i = 0; i < k; ++i){ out << (i == k - 1 ? to_string(widths[i]) + ")" : to_string(widths[i]) + ", "); }
        out << std::endl;
        if(h > 0){
            out << "heights(";
            for(int i = 0; i < k; ++i){ out << (i == k - 1 ? to_string(heights[i]) + ")" : to_string(heights[i]) + ", "); }    
        }
        out << std::endl;
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
            std::cout << "PARSED : " << std::endl;
            problem.print(); 
            return problem;
        }
    };
};