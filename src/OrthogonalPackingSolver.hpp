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

    OrthogonalPackingProblem(int, int, int, int, int);
    OrthogonalPackingProblem(const OrthogonalPackingProblem&);

    OrthogonalPackingProblem& operator=(const OrthogonalPackingProblem&);

    virtual ~OrthogonalPackingProblem();

    void print(std::ostream& = std::cout);


    struct Parser{
        struct ParseException : public std::runtime_error {
            ParseException(const std::string& msg) : std::runtime_error(msg) {}
        };

        inline static int next_int(std::string&);
        inline static OrthogonalPackingProblem parse(std::istream&, bool = false);
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