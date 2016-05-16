#include <iostream>
#include "Solver.hpp"

using namespace std;

// dimension de la grille, doit etre inferieure ou egale à 9
#define N 4
#define FOR(k,lb,ub) for (int k = (lb); (k) <= (ub); (k)++)


int main() {
	Solver s;
	vec<Lit> lits;
	
    // call the SAT solver
    s.solve();

    // récupération de la solution
    while (s.okay()) {
        lits.clear();
        cout << "Nouvelle solution: \n\n" ;
        s.addClause(lits);
        s.solve();
    }
    cout << "Il n'y a plus de solutions\n";
}
