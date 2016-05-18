#include <iostream>
#include "OrthogonalPackingSolver.hpp"

int main() {
    try{
        OrthogonalPackingSolver solver(OrthogonalPackingProblem::Parser::parse(std::cin));
        solver.solve();
        solver.print_solution(std::cout);
        return 0;
    } catch(const std::exception& e){
        std::cout << e.what() << std::endl;
        return 1;
    }
}