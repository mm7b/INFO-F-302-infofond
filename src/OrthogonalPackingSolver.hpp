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

    OrthogonalPackingSolution(const OrthogonalPackingProblem&, bool);

    OrthogonalPackingSolution(const OrthogonalPackingSolution&);

    OrthogonalPackingSolution& operator=(const OrthogonalPackingSolution&);

    int* operator[](int);  

    virtual ~OrthogonalPackingSolution();
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