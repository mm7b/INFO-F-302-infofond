#ifndef OrthogonalPackingSolver_h
#define OrthogonalPackingSolver_h

#include <string>
#include <sstream>
#include <iostream>
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

    virtual ~OrthogonalPackingProblem(){
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

    void selfGenerateNAndM(){
        for(int i = 0; i<k; i++){
            m += lengths[i];
            n += widths[i];
        }
        if (n>m){
            m = n;
        }else{
            n = m;
        }
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

        inline static OrthogonalPackingProblem parse(std::istream& in, bool three_dim = false, bool n_and_m_parsing = true){
            int dim = three_dim ? 3 : 2;
            std::string input_line;
            std::getline(in, input_line);
            int k = next_int(input_line);
            int n = 0;
            int m = 0;
            if(n_and_m_parsing){
                std::getline(in, input_line);
                n = next_int(input_line);
                std::getline(in, input_line);
                m = next_int(input_line);
            }
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
            if(!n_and_m_parsing){
                problem.selfGenerateNAndM();
            }
            return problem;
        }
    };
};

struct OrthogonalPackingSolution {
    static const std::string PYTHON_PLOTTER_FILENAME;

    OrthogonalPackingProblem problem;
    int** solution;
    bool exists;

    OrthogonalPackingSolution(const OrthogonalPackingProblem& p, bool e) : problem(p), solution(NULL), exists(e) {
        solution = new int*[problem.k];
        for(int k = 0; k < problem.k; ++k){
            solution[k] = new int[problem.dim];
            for(int d = 0; d < problem.dim; ++d){
                solution[k][d] = -1;
            }
        }
    }

    OrthogonalPackingSolution(const OrthogonalPackingSolution& other) : problem(other.problem), solution(NULL), exists(other.exists){
        solution = new int*[problem.k];
        for(int k = 0; k < problem.k; ++k){
            solution[k] = new int[problem.dim];
            std::copy(other.solution[k], other.solution[k] + problem.dim, solution[k]);
        }
    }

    OrthogonalPackingSolution& operator=(const OrthogonalPackingSolution& other){
        if(this != &other){
            if(solution != NULL) {
                for(int k = 0; k < problem.k; ++k){
                    if(solution[k] != NULL) { delete[] solution[k]; }
                }
                delete[] solution;
            }
            problem = other.problem; exists = other.exists;
            solution = new int*[problem.k];
            for(int k = 0; k < problem.k; ++k){
                solution[k] = new int[problem.dim];
                std::copy(other.solution[k], other.solution[k] + problem.dim, solution[k]);
            }
        }
        return *this;
    }

    int* operator[](int i){
        if(i < 0 || i >= problem.k) { throw std::out_of_range("Rectangle " + to_string(i) + " does not exist"); }
        return solution[i];
    }

    virtual ~OrthogonalPackingSolution(){ 
        if(solution != NULL) {
            for(int k = 0; k < problem.k; ++k){
                if(solution[k] != NULL) { delete[] solution[k]; }
            }
            delete[] solution;
        }
    }
};

class OrthogonalPackingSolver : public Solver {
private:
	OrthogonalPackingProblem problem;

    int*** mu;

    bool overlapping(int, int, int, int, int, int);

    bool out_of_bounds(int, int, int);

public:
	OrthogonalPackingSolver(const OrthogonalPackingProblem&);

	void add_constraints();

    OrthogonalPackingSolution get_solution();

	void print_solution(std::ostream&);
    void plot_solution();

	virtual ~OrthogonalPackingSolver();
};

#endif