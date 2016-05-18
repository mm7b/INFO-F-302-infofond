#include <iostream>
#include "OrthogonalPackingSolver.hpp"

int orthogonalPacking2D(){
	try{
        OrthogonalPackingSolver solver(OrthogonalPackingProblem::Parser::parse(std::cin));
        solver.solve();
        solver.print_solution(std::cout);
        solver.plot_solution();
        return 0;
    } catch(const std::exception& e){
        std::cout << e.what() << std::endl;
        return 1;
    }
}

int tiniestSquare(){
	OrthogonalPackingSolver solver(OrthogonalPackingProblem::Parser::parse(std::cin, false, false));
	solver.solve();
	int num = 0;
	solver.print_solution(std::cout);
	while(solver.get_solution().exists && num < 21){
		num++;
		solver.solve();
	}
	solver.print_solution(std::cout);
	solver.plot_solution();
	std::cout << num << std::endl;
	return 0;
}

int main() {

	std::cout << "Please choose the action to perform:" << std::endl;
	std::cout << "1) 2D orthogonal packing solver" << std::endl;
	std::cout << "2) Tiniest square to put all the rectangles in" << std::endl;

	int action = 0;
	std::string input = "";
 	while (true) {
	   std::getline(std::cin, input);
	   std::stringstream myStream(input);
	   if (myStream >> action)
	     break;
	   std::cout << "Invalid number, please try again" << std::endl;
	 }
	if(action == 1){
	    return orthogonalPacking2D();
	} else if (action == 2) {
		return tiniestSquare();
	} else {
		std::cout << "Invalid action, please choose between 1 & 2" << std::endl;
		return main();
	}
}