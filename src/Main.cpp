#include <iostream>
#include "OrthogonalPackingSolver.hpp"
#include <map>

bool is_number(const std::string& s){
    std::string::const_iterator it = s.begin();
    while (it != s.end() && std::isdigit(*it)) ++it;
    return !s.empty() && it == s.end();
}

enum ProblemType : int { Q3 = 1, Q4 = 2, Q5 = 3, Q6 = 4, Q7 = 5, Q8 = 6, Q9 = 7, Q10 = 8 }

OrthogonalPackingProblem (*parse)(std::istream&, bool, bool, bool, bool, bool) = OrthogonalPackingProblem::Parser::parse;

OrthogonalPackingProblem build_problem(int problem_type, std::istream& in){
    if(problem_type == ProblemType.Q3)      { return parse(in, Dimension.DIM_2, Solution.ANY, Height.FLOAT, Orientation.FIX, EdgesUnit.FREE); }
    else if(problem_type == ProblemType.Q4) { return parse(in, Dimension.DIM_2, Solution.SMALLEST, Height.FLOAT, Orientation.FIX, EdgesUnit.FREE); }
    else if(problem_type == ProblemType.Q5) { return parse(in, Dimension.DIM_2, Solution.SMALLEST, Height.FLOAT, Orientation.FIX, EdgesUnit.FREE); }
    else if(problem_type == ProblemType.Q6) { return parse(in, Dimension.DIM_3, Solution.ANY, Height.FLOAT, Orientation.FIX, EdgesUnit.FREE); }
    else if(problem_type == ProblemType.Q7) { return parse(in, Dimension.DIM_3, Solution.ANY, Height.NO_FLOAT, Orientation.FIX, EdgesUnit.FREE); }
    else if(problem_type == ProblemType.Q8) { return parse(in, Dimension.DIM_2, Solution.ANY, Height.FLOAT, Orientation.PIVOT, EdgesUnit.FREE); }
    else if(problem_type == ProblemType.Q9) { return parse(in, Dimension.DIM_2, Solution.ANY, Height.FLOAT, Orientation.FIX, EdgesUnit.MINIMUM); }
    else if(problem_type == ProblemType.Q10){ return parse(in, Dimension.DIM_2, Solution.SMALLEST, Height.FLOAT, Orientation.FIX, EdgesUnit.FREE); }
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

bool display_menu(){
    std::istream in = std::cin;
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
    if(is_number(input))    { orthogonal_packing(build_problem(atoi(input), in));}
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