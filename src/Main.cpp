#include <iostream>
#include "OrthogonalPackingSolver.hpp"
#include <map>

void orthogonal_packing(bool is_3d){
	try{
        OrthogonalPackingSolver solver(OrthogonalPackingProblem::Parser::parse(std::cin, is_3d));
        solver.solve();
        solver.print_solution(std::cout);
        solver.plot_solution();
    } catch(const std::exception& e){
        std::cout << e.what() << std::endl;
    }
}

const std::string prompt = ">> ";

bool display_menu(){
    std::string input;
    std::cout   << "Choisissez :" << std::endl
                << "(1) 2D orthogonal packing solver" << std::endl
                << "(2) 3D orthogonal packing solver" << std::endl
                << "(3) Plus petit carré" << std::endl
                << "(q) Quitter" << std::endl
                << prompt;
    std::getline(std::cin, input);
    if(input == "1")        { orthogonal_packing(false); }
    else if(input == "2")   { orthogonal_packing(true); }
    else if(input == "3")   { orthogonal_packing(false); }
    else if(input == "q")   { return true; }
    else{ std::cout << "Entrée non valide !" << std::endl; }
    return false;
}

int main() {
    bool stop = false;
    do{
        stop = display_menu();
    } while(!stop);
    return 0;
}