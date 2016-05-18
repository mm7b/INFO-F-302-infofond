#include <iostream>
#include "OrthogonalPackingSolver.hpp"
#include <map>

bool is_number(const std::string& s){
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

enum ProblemType { Q3 = 1, Q4 = 2, Q5 = 3, Q6 = 4, Q7 = 5, Q8 = 6, Q9 = 7, Q10 = 8 };

OrthogonalPackingProblem (*parse)(std::istream&, Dimension, SolutionType, HeightConstraint, Orientation, EdgeContact) = OrthogonalPackingProblem::Parser::parse;

OrthogonalPackingProblem build_problem(int type, std::istream& in){
    Dimension dimension; SolutionType solution_type; 
    HeightConstraint height_constraint; 
    Orientation orientation; EdgeContact edge_contact;
    switch(type){
        case Q3:
            dimension = DIM_2; solution_type = ANY; height_constraint = FLOAT; orientation = FIX; edge_contact = FREE;
            break;
        case Q4:
            dimension = DIM_2; solution_type = SMALLEST; height_constraint = FLOAT; orientation = FIX; edge_contact = FREE;
            break;
        case Q5:
            dimension = DIM_2; solution_type = SMALLEST; height_constraint = FLOAT; orientation = FIX; edge_contact = FREE;
            break;
        case Q6:
            dimension = DIM_3; solution_type = ANY; height_constraint = FLOAT; orientation = FIX; edge_contact = FREE;
            break;
        case Q7:
            dimension = DIM_3; solution_type = ANY; height_constraint = NO_FLOAT; orientation = FIX; edge_contact = FREE;
            break;
        case Q8:
            dimension = DIM_2; solution_type = ANY; height_constraint = FLOAT; orientation = PIVOT; edge_contact = FREE;
            break;
        case Q9:
            dimension = DIM_2; solution_type = ANY; height_constraint = FLOAT; orientation = FIX; edge_contact = MINIMUM;
            break;
        case Q10:
            dimension = DIM_2; solution_type = SMALLEST; height_constraint = FLOAT; orientation = FIX; edge_contact = FREE;
            break;
    }
    return parse(in, dimension, solution_type, height_constraint, orientation, edge_contact);
}

void orthogonal_packing(const OrthogonalPackingProblem& problem){
	try{
        OrthogonalPackingSolver solver(problem);
        solver.solve();
        solver.print_solution(std::cout);
        solver.plot_solution();
    } catch(const std::exception& e){
        std::cout << e.what() << std::endl;
    }
}

const std::string prompt = ">> ";

bool display_menu(std::istream& in = std::cin){
    std::string input;
    std::cout   << "Choisissez :" << std::endl
                << "(1) Q3  : 2D orthogonal packing" << std::endl
                << "(2) Q4  : Plus petit carré" << std::endl
                << "(3) Q5  : Plus petit carré avec carrés de taille (1...n)" << std::endl
                << "(4) Q6  : 3D orthogonal packing" << std::endl
                << "(5) Q7  : 3D orthogonal packing sans pavés flottants" << std::endl
                << "(6) Q8  : 2D orthogonal packing avec pivot" << std::endl
                << "(7) Q9  : 2D orthogonal packing avec minimum p unités de contact" << std::endl
                << "(8) Q10 : 2D orthogonal packing avec utilisation de MAX-SAT" << std::endl
                << "(q) Quitter" << std::endl
                << prompt;
    std::getline(in, input);
    if(is_number(input))    { orthogonal_packing(build_problem(atoi(input.c_str()), in));}
    else if(input == "q")   { return true; }
    else{ std::cout << "Entrée non valide !" << std::endl; }
    return false;
}

int main() {
    bool stop = false;
    do{
        stop = display_menu();
    } while(!stop);

    return EXIT_SUCCESS;
}