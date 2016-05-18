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

enum Dimension : bool {
    DIM_2 = true,
    DIM_3 = false
}

enum Solution : bool {
    SMALLEST = true,
    ANY = false
}

enum Height : bool {
    FLOAT = true,
    NO_FLOAT = false
}

enum Orientation : bool {
    PIVOT = true,
    FIX = false
}

enum EdgesUnit : bool {
    MINIMUM = true,
    FREE = false
}

struct OrthogonalPackingProblem{
    int k, dim, n, m, h;
    int* lengths; int* widths; int* heights;

    OrthogonalPackingProblem(int, int, int, int, int);
    OrthogonalPackingProblem(const OrthogonalPackingProblem&);
    OrthogonalPackingProblem& operator=(const OrthogonalPackingProblem&);
    void print(std::ostream& = std::cout);
    bool is_3d() const { return dim == 3; }
    virtual ~OrthogonalPackingProblem();

    struct Parser{
        struct ParseException : public std::runtime_error {
            ParseException(const std::string& msg) : std::runtime_error(msg) {}
        };

        static int next_int(std::string&);
        static OrthogonalPackingProblem parse(std::istream&, bool, bool, bool, bool, bool!);
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
    int**** mu;

    bool out_of_bounds(int, int, int, int);
    bool overlapping(int, int, int, int, int, int, int, int);

public:
	OrthogonalPackingSolver(const OrthogonalPackingProblem&);

	void add_constraints();

    OrthogonalPackingSolution get_solution();

	void print_solution(std::ostream&);
    void plot_solution();

	virtual ~OrthogonalPackingSolver();
};

#endif