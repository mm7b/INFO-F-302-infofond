#ifndef OrthogonalPackingSolver_h
#define OrthogonalPackingSolver_h

#include <string>
#include <sstream>
#include <iostream>
#include <stdlib.h>
#include <exception>
#include <stdexcept>
#include <bitset>

#include "Solver.hpp"

#define BITSET_SIZE 32

template<typename T>
std::string to_string(const T& value){
    std::ostringstream oss;
    oss << value;
    return oss.str();
}

enum RectanglesSource {
    GENERATE = true,
    FROM_INPUT = false
};

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
    int k, dim, n, m, h, min_n;
    SolutionType solution_type; HeightConstraint height_constraint;
    Orientation orientation; EdgeContact edge_contact;
    int* lengths; int* widths; int* heights;

    OrthogonalPackingProblem(int, int, int, int, int, SolutionType, HeightConstraint, Orientation, EdgeContact);
    OrthogonalPackingProblem(const OrthogonalPackingProblem&);
    OrthogonalPackingProblem& operator=(const OrthogonalPackingProblem&);
    void print(std::ostream& = std::cout);
    void selfGenerateNAndM();
    void generateMinN();
    bool is_3d() const { return dim == 3; }
    virtual ~OrthogonalPackingProblem();

    struct Parser{
        struct ParseException : public std::runtime_error {
            ParseException(const std::string& msg) : std::runtime_error(msg) {}
        };
        static int next_int(std::string&);
        static OrthogonalPackingProblem parse(std::istream&, RectanglesSource, Dimension, SolutionType, HeightConstraint, Orientation, EdgeContact);
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
	void print(std::ostream&);
    void plot(int, int);
    virtual ~OrthogonalPackingSolution();
};

struct OrthogonalPackingSolver : public Solver {
    OrthogonalPackingProblem problem;
    int**** mu;
    int* pivot;
    int* dimension;
    int***** in_bounds;

    bool carry(int a, int b, int c, int d, int e, int f, int k, int l);

    OrthogonalPackingSolver(const OrthogonalPackingProblem&);

    void add_constraints();

    OrthogonalPackingSolution get_solution();


	virtual ~OrthogonalPackingSolver();
};

#endif