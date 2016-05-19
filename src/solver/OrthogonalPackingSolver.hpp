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

enum Dimension {
    DIM_2 = true,
    DIM_3 = false
};

enum SolutionType {
    SMALLEST = true,
    ANY = false
};

enum HeightConstraint {
    FLOAT = true,
    NO_FLOAT = false
};

enum Orientation {
    PIVOT = true,
    FIX = false
};

enum EdgeContact {
    MINIMUM = true,
    FREE = false
};

struct OrthogonalPackingProblem{
    int k, dim, n, m, h;
    SolutionType solution_type; HeightConstraint height_constraint;
    Orientation orientation; EdgeContact edge_contact;
    int* lengths; int* widths; int* heights;

    OrthogonalPackingProblem(int, int, int, int, int, SolutionType, HeightConstraint, Orientation, EdgeContact);
    OrthogonalPackingProblem(const OrthogonalPackingProblem&);
    OrthogonalPackingProblem& operator=(const OrthogonalPackingProblem&);
    void print(std::ostream& = std::cout);
    void selfGenerateNAndM();
    bool is_3d() const { return dim == 3; }
    virtual ~OrthogonalPackingProblem();

    struct Parser{
        struct ParseException : public std::runtime_error {
            ParseException(const std::string& msg) : std::runtime_error(msg) {}
        };
        static int next_int(std::string&);
        static OrthogonalPackingProblem parse(std::istream&, Dimension, SolutionType, HeightConstraint, Orientation, EdgeContact);
    };
};

struct OrthogonalPackingSolution {
    static const std::string PYTHON_PLOTTER_FILENAME;

    OrthogonalPackingProblem problem;
    int** solution; int* pivot;
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
    int* pivot;

    bool out_of_bounds(int, int, int, int);
    bool pivot_out_of_bounds(int, int, int, int);
    bool overlapping(int, int, int, int, int, int, int, int);
    bool pivot_overlapping(int, int, int, int, int, int, int, int);
    bool carry(int, int, int, int, int, int, int, int);

public:
	OrthogonalPackingSolver(const OrthogonalPackingProblem&);

	void add_constraints();

    OrthogonalPackingSolution get_solution();

	void print_solution(std::ostream&);
    void plot_solution();

	virtual ~OrthogonalPackingSolver();
};

#endif