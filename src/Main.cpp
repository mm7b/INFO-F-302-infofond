#include <iostream>
#include "Solver.hpp"

using namespace std;

#define FOR(k,lb,ub) for (int k = (lb); (k) <= (ub); (k)++)

/* Consume chars from stdin, return true if ';' or false if ',' */
static inline bool skip(istream & input)
{
    char c;
    while (1) {
        input >> c;
        if (c == ',') return false;
        else if (c == ';') return true;
    }
}  

static inline void next_int(istream & input, int & dest, const char *title)
{
    input >> dest;
    skip(input);
}


int main() {
    Solver s;
    vec<Lit> lits;

    int X, Y, K;
    istream & input = cin;
    next_int(input, K, "Nombre de rectangles");
    next_int(input, X, "Largeur");
    next_int(input, Y, "Longueur");

    FOR(r, 1, K){

    }


    // exemple si le R est de taille 6-7-1 (2D quoi)
    int prop[X][Y][1];

    FOR(x,1,X-1){
        FOR(y,1,Y-1){
            FOR(z, 1,1){
                prop[x][y][z] = s.newVar();
            }
        }
    }

    //ajout des contraintes


    
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

