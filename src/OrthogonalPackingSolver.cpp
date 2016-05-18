#include <cctype>
#include <algorithm>
#include "OrthogonalPackingSolver.hpp"

OrthogonalPackingSolver::OrthogonalPackingSolver(const OrthogonalPackingProblem& p) : Solver(), problem(p), mu(NULL) {

    add_constraints();

}

void OrthogonalPackingSolver::add_constraints(){

    /* Initialisation des prop */
    mu = new int**[problem.k];
    for(int k = 0; k < problem.k; ++k){
        mu[k] = new int*[problem.m];
        for(int a = 0; a < problem.m; ++a){
            mu[k][a] = new int[problem.n];
            for(int b = 0; b < problem.n; ++b){
               mu[k][a][b] = newVar();
            }
        }

    }

    /* On ne peut avoir 2 mu à vrai en même temps pour un même k, 
    un rectangle k ne pouvant avoir qu'une seule position */
    for(int k = 0; k < problem.k; ++k){
        for(int a = 0; a < problem.m; ++a){
            for(int b = 0; b < problem.n; ++b){
                for(int d = a + 1; d < problem.m; ++d){
                    for(int e = b + 1; e < problem.n; ++e){
                        addBinary(~Lit(mu[k][a][b]), ~Lit(mu[k][d][e]));
                    }
                }
            }
        }
    }

    /* On doit avoir au mois 1 mu à vrai pour un rectangle k donné  */
    vec<Lit> lits;
    for(int k = 0; k < problem.k; ++k){
        lits.clear();
         for(int a = 0; a < problem.m; ++a){
            for(int b = 0; b < problem.n; ++b){
                lits.push(Lit(mu[k][a][b]));
            }
        }
        addClause(lits);
    }

    /* (a, b) doit être dans les bornes du grand rectangle */
    for(int k = 0; k < problem.k; ++k){
         for(int a = 0; a < problem.m; ++a){
            for(int b = 0; b < problem.n; ++b){
                if(out_of_bounds(a, b, k)){
                    addUnit(~Lit(mu[k][a][b]));
                }
            }
        }
    }

    /* Pas de superposition : un rectangle est soit à gauche, soit à droite, 
    soit en haut, soit en bas, soit plus profond, soit moins profond que tous les autres */
    for(int k = 0; k < problem.k; ++k){
        for(int a = 0; a < problem.m; ++a){
            for(int b = 0; b < problem.n; ++b){
                for(int l = k + 1; l < problem.k; ++l){
                    for(int d = 0; d < problem.m; ++d){
                        for(int e = 0; e < problem.n; ++e){
                            if(overlapping(a, b, d, e, k, l)){
                                addBinary(~Lit(mu[k][a][b]), ~Lit(mu[l][d][e]));
                            }
                        }
                    }
                }
            }
        }
    }


}

/* Inutile de vérifier a >= 0 ni b >=0 car se sont des indices donc >=0 par définition */
bool OrthogonalPackingSolver::out_of_bounds(int a, int b, int k){
    return !(a + problem.lengths[k] <= problem.m && b + problem.widths[k] <= problem.n);
}

bool OrthogonalPackingSolver::overlapping(int a, int b, int d, int e, int k, int l){
    return !( a + problem.lengths[k] <= d
        ||  a >= d + problem.lengths[l]
        ||  b + problem.widths[k] <= e
        ||  b >= e + problem.widths[l]);
}

void OrthogonalPackingSolver::print_solution(std::ostream& out = std::cout){
    if(mu == NULL){ throw std::runtime_error("Unexpected error : prop vector is null"); }
    if(! okay()){
        out << 0 << std::endl;
    }
    else{
        for(int k = 0; k < problem.k; ++k){
            for(int a = 0; a < problem.m; ++a){
                for(int b = 0; b < problem.n; ++b){
                    if(model[mu[k][a][b]] == l_True){
                        out << k + 1 << " " << a << " " << b << std::endl;
                    }
                }
            }
        }
    }
}

virtual OrthogonalPackingSolver::~OrthogonalPackingSolver(){
    for(int k = 0; k < problem.k; ++k){
        for(int a = 0; a < problem.m; ++a){
            delete[] mu[k][a];
        }
        delete[] mu[k];
    }
    delete[] mu;
}