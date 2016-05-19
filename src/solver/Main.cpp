#include <iostream>
#include <sstream>
#include <utility>
#include <algorithm>
#include "OrthogonalPackingSolver.hpp"
#include <map>

bool is_number(const std::string& s){
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

std::pair<bool, int> parse_number(const std::string& s){
    std::string number_str = "";
    std::string::const_iterator it = s.begin();
    while (it != s.end()){
        if(std::isdigit(*it)){ number_str += *it; }
        ++it;
    }
    if(number_str.empty()) { return std::pair<bool, int>(false, 0); }
    else{ return std::pair<bool, int>(true, atoi(number_str.c_str())); }
}

void orthogonal_packing(const OrthogonalPackingProblem& problem){
    OrthogonalPackingSolver solver(problem);
    solver.solve();
    if(problem.solution_type == SMALLEST){
        OrthogonalPackingSolution sol = solver.get_solution();
        int minimum_square_size = problem.n;
        while(sol.exists){
            solver.addUnit(~Lit(solver.dimension[minimum_square_size - problem.min_n]));
            solver.solve();
            sol = solver.get_solution();
            --minimum_square_size;
        }
        std::cout << "Tiniest square size: " << minimum_square_size << std::endl;
    }else{
        solver.print_solution(std::cout);
        solver.plot_solution();
    }

}

enum ProblemType { Q3 = 1, Q4 = 2, Q5 = 3, Q6 = 4, Q7 = 5, Q8 = 6, Q9 = 7, Q10 = 8, MIN_Q = Q3, MAX_Q = Q10 };

OrthogonalPackingProblem build_problem(int question, std::istream& in){
    Dimension dimension; SolutionType solution_type; 
    HeightConstraint height_constraint; 
    Orientation orientation; EdgeContact edge_contact;
    switch(question){
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
        default:
            throw std::runtime_error("Unsupported question");
            break;
    }
    return OrthogonalPackingProblem::Parser::parse(in, dimension, solution_type, height_constraint, orientation, edge_contact);
}

int from_arg(const std::string& arg, std::istream& in){
    std::pair<bool, int> parsed = parse_number(arg);
    if(!parsed.first){ throw std::runtime_error("Invalid program argument : please give a question (e.g. q3, Q3, 3D)"); }
    return parsed.second - 2;
}

const std::string prompt = ">> ";

int from_menu(std::istream& in){
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
    if(is_number(input)) {
        int question_input = atoi(input.c_str());
        return (MIN_Q <= question_input && question_input <= MAX_Q) ? question_input : -1; 
    }
    else if(input == "q")   { return 0; }
    else{ 
        std::cout << "Entrée non valide !" << std::endl;
        return -1;
    }
}

int main(int argc, const char** argv) {
    std::istream& in = std::cin;
    try{
        if(argc == 2){
            orthogonal_packing(build_problem(from_arg(std::string(argv[1]), in), in));
        }
        else{
            bool stop = false; bool valid = false;
            int input = -1;
            do{
                try{
                    input = from_menu(in);
                    stop = (input == 0); valid = (input > 0);
                    if(valid) { orthogonal_packing(build_problem(input, in)); }
                } 
                catch(const OrthogonalPackingProblem::Parser::ParseException& e){ std::cerr << e.what() << std::endl; }
            } while(!stop);
        }
        return EXIT_SUCCESS;
    }
    catch(const std::exception& e){
        std::cerr << e.what() << std::endl;
        return EXIT_FAILURE;
    }
}