#include <cctype>
#include <algorithm>
#include <stdlib.h>
#include <sstream>
#include <string>
#include <math.h>
#include <bitset>
#include "OrthogonalPackingSolver.hpp"


/* OrthogonalPackingProblem definitions */

OrthogonalPackingProblem::OrthogonalPackingProblem(
    int parsed_k, int parsed_dim, int parsed_n, int parsed_m, int parsed_h,
    SolutionType config_solution, HeightConstraint config_height, 
    Orientation config_orientation, EdgeContact config_edge_contact) : 
        k(parsed_k), dim(parsed_dim), n(parsed_n), m(parsed_m), h(parsed_h), min_n(0),
        solution_type(config_solution), height_constraint(config_height), 
        orientation(config_orientation), edge_contact(config_edge_contact),
        lengths(new int[k]), widths(new int[k]), heights(new int[k]) {}

OrthogonalPackingProblem::OrthogonalPackingProblem(const OrthogonalPackingProblem& other) : 
    k(other.k), dim(other.dim), n(other.n), m(other.m), h(other.h), min_n(other.min_n),
    solution_type(other.solution_type), height_constraint(other.height_constraint), 
    orientation(other.orientation), edge_contact(other.edge_contact),
    lengths(new int[k]), widths(new int[k]), heights(new int[k]) {
        std::copy(other.lengths, other.lengths + k, lengths);
        std::copy(other.widths, other.widths + k, widths);
        std::copy(other.heights, other.heights + k, heights);
}

OrthogonalPackingProblem& OrthogonalPackingProblem::operator=(const OrthogonalPackingProblem& other){
    if(this != &other){
        k = other.k; dim = other.dim; n = other.n; m = other.m; h = other.h; min_n = other.min_n;
        solution_type = other.solution_type; height_constraint = other.height_constraint; 
        orientation = other.orientation; edge_contact = other.edge_contact;
        delete[] lengths; delete[] widths; delete[] heights;
        lengths = new int[k]; widths = new int[k]; heights = new int[k];
        std::copy(other.lengths, other.lengths + k, lengths);
        std::copy(other.widths, other.widths + k, widths);
        std::copy(other.heights, other.heights + k, heights);
    }
    return *this;
}

OrthogonalPackingProblem::~OrthogonalPackingProblem(){
    delete[] lengths;
    delete[] widths;
    delete[] heights;
}

void OrthogonalPackingProblem::print(std::ostream& out){
    out << "k(" + to_string(k) + ") " << "dim(" + to_string(dim) + ") "
        << "n(" + to_string(n) + ") " << "m(" + to_string(m) + ") "
        << (h > 0 ? "h(" + to_string(h) + ") " : "") << std::endl
        << "lengths(";
    for(int i = 0; i < k; ++i){ out << (i == k - 1 ? to_string(lengths[i]) + ")" : to_string(lengths[i]) + ", "); }
    out << std::endl << "widths(";
    for(int i = 0; i < k; ++i){ out << (i == k - 1 ? to_string(widths[i]) + ")" : to_string(widths[i]) + ", "); }
    out << std::endl;
    if(h > 0){
        out << "heights(";
        for(int i = 0; i < k; ++i){ out << (i == k - 1 ? to_string(heights[i]) + ")" : to_string(heights[i]) + ", "); }    
    }
    out << std::endl;
}

void OrthogonalPackingProblem::selfGenerateNAndM(){
    for(int i = 0; i<k; i++){
        m += lengths[i];
        n += widths[i];
    }
    if (n>m){
        m = n;
    }else{
        n = m;
    }
 }

void OrthogonalPackingProblem::generateMinN(){
	for(int i=0; i<k; i++){
		min_n += lengths[i]*widths[i];
	}
	min_n = sqrt(min_n);
}

/* Parser of OrthogonalPackingProblem definitions */

int OrthogonalPackingProblem::Parser::next_int(std::string& line){
    std::size_t i = 0;
    std::string digit = "";
    while(i < line.length() && isdigit(line[i])){
        digit += line[i];
        ++i; 
    }
    if(i == 0){ throw ParseException("Could not read int from : " + line); }

    line = (++i < line.length()) ? line.substr(i, line.length() - i) : "";
    return atoi(digit.c_str());
}

OrthogonalPackingProblem OrthogonalPackingProblem::Parser::parse(
    std::istream& in, RectanglesSource rectangles_source,
    Dimension dimension, SolutionType solution, HeightConstraint height, 
    Orientation orientation, EdgeContact edge_contact){
        int dim = dimension == DIM_3 ? 3 : 2;
        std::string input_line;
        std::getline(in, input_line);
        int k = next_int(input_line);
        int n = 0;
        int m = 0;
        if(solution == ANY){
            std::getline(in, input_line);
            n = next_int(input_line);
            std::getline(in, input_line);
            m = next_int(input_line);
        }
        int h = -1;
        if(dimension == DIM_3){
            std::getline(in, input_line);
            h = next_int(input_line);
        }
        OrthogonalPackingProblem problem(k, dim, n, m, h, solution, height, orientation, edge_contact);
    	for(int i = 0; i < k; ++i){
	        if(rectangles_source == FROM_INPUT){
		            std::getline(in, input_line);
		            int index = next_int(input_line);
		            if(index != i + 1) { throw ParseException("Wrong index : expected " + to_string(i + 1) + " found " + to_string(index) + " instead"); }
		            problem.lengths[i] = next_int(input_line);
		            problem.widths[i] = next_int(input_line);   
		            problem.heights[i] = dimension == DIM_3 ? next_int(input_line) : -1;
		    }else{
		    	problem.lengths[i] = i+1;
		    	problem.widths[i] = i+1;
		    	problem.heights[i] = -1;
	    	}
	    }
        if(solution == SMALLEST){
            problem.selfGenerateNAndM();
            problem.generateMinN();
        }
        return problem;
}

/* OrthogonalPackingSolution definitions */

const std::string OrthogonalPackingSolution::PYTHON_PLOTTER_FILENAME = "plot_opp.py";

OrthogonalPackingSolution::OrthogonalPackingSolution(const OrthogonalPackingProblem& p, bool e) : 
    problem(p), solution(NULL), pivot(NULL), exists(e) {
        solution = new int*[problem.k];
        for(int k = 0; k < problem.k; ++k){
            solution[k] = new int[problem.dim];
            for(int d = 0; d < problem.dim; ++d){
                solution[k][d] = -1;
            }
        }
        if(problem.orientation == PIVOT){
           pivot = new int[problem.k]; 
        }
}

OrthogonalPackingSolution::OrthogonalPackingSolution(const OrthogonalPackingSolution& other) : 
    problem(other.problem), solution(NULL), pivot(NULL), exists(other.exists){
        solution = new int*[problem.k];
        for(int k = 0; k < problem.k; ++k){
            solution[k] = new int[problem.dim];
            std::copy(other.solution[k], other.solution[k] + problem.dim, solution[k]);
        }
        if(problem.orientation == PIVOT){
           pivot = new int[problem.k]; 
           std::copy(other.pivot, other.pivot + problem.k, pivot);
        }
}

OrthogonalPackingSolution& OrthogonalPackingSolution::operator=(const OrthogonalPackingSolution& other){
    if(this != &other){
        if(solution != NULL) {
            for(int k = 0; k < problem.k; ++k){
                if(solution[k] != NULL) { delete[] solution[k]; }
            }
            delete[] solution;
        }
        if(pivot != NULL) { delete[] pivot; }
        problem = other.problem; exists = other.exists;
        solution = new int*[problem.k];
        for(int k = 0; k < problem.k; ++k){
            solution[k] = new int[problem.dim];
            std::copy(other.solution[k], other.solution[k] + problem.dim, solution[k]);
        }
        if(problem.orientation == PIVOT){
           pivot = new int[problem.k]; 
           std::copy(other.pivot, other.pivot + problem.k, pivot);
        }
        
    }
    return *this;
}

int* OrthogonalPackingSolution::operator[](int i){
    if(i < 0 || i >= problem.k) { throw std::out_of_range("Rectangle " + to_string(i) + " does not exist"); }
    return solution[i];
}

void OrthogonalPackingSolution::print(std::ostream& out = std::cout){
    if(! exists){
        out << 0 << std::endl;
    }
    else{
        for(int k = 0; k < problem.k; ++k){
            out << k + 1 << " ";
            for(int d = 0; d < problem.dim; ++d){
                out << solution[k][d] << (d < problem.dim - 1 ? " " : "");
            }
            if(!(problem.orientation == FIX)) { out << " " << pivot[k]; }
            out << std::endl;
        }
    }
}

void OrthogonalPackingSolution::plot(int n, int m){
    if(exists){
        std::ostringstream oss;
        oss << "python " << OrthogonalPackingSolution::PYTHON_PLOTTER_FILENAME << " " << problem.k << " " << n << " " 
            << m << " " << problem.h;

        oss << " \"[";
        for(int k = 0; k < problem.k; ++k){
            oss << "(";
            for(int d = 0; d < problem.dim; ++d){
                oss << (d == problem.dim - 1 ? to_string(solution[k][d]) + ")" : to_string(solution[k][d]) + ", ");
            }
            if(k < problem.k - 1) { oss << ", "; }
        }
        oss << "]\"";
        
        int length, width;

        oss << " \"[";
        for(int k = 0; k < problem.k; ++k){ 
            if(!(problem.orientation == FIX)){
                length = (pivot[k] ? problem.widths[k] : problem.lengths[k]);
            }
            else{ length = problem.lengths[k]; }
            oss << (k == problem.k - 1 ? to_string(length) : to_string(length) + ", ");
        }
        oss << "]\"";
        
        oss << " \"[";
        for(int k = 0; k < problem.k; ++k){
            if(!(problem.orientation == FIX)){
                width = (pivot[k] ? problem.lengths[k] : problem.widths[k]);
            }
            else{ width = problem.widths[k]; }
            oss << (k == problem.k - 1 ? to_string(width) : to_string(width) + ", ");
        }
        oss << "]\"";

        if(problem.dim == 3){
            oss << " \"[";
            for(int k = 0; k < problem.k; ++k){ oss << (k == problem.k - 1 ? to_string(problem.heights[k]) : to_string(problem.heights[k]) + ", "); }
            oss << "]\"";
        }

        oss << " --color=b" << " --alpha=0.4";
        pid_t pid = fork();
        if(pid < 0){
            throw std::runtime_error("Failed to execute plotting command");
        }
        else if (pid == 0){
            if(system(NULL)){
                system(oss.str().c_str());
                _exit(EXIT_SUCCESS);
            }
            else{
                _exit(EXIT_FAILURE);
            }
        }
    }
}

OrthogonalPackingSolution::~OrthogonalPackingSolution(){ 
    if(solution != NULL) {
        for(int k = 0; k < problem.k; ++k){
            if(solution[k] != NULL) { delete[] solution[k]; }
        }
        delete[] solution;
    }
    if(pivot != NULL){
        delete[] pivot;
    }
}

OrthogonalPackingSolver::OrthogonalPackingSolver(const OrthogonalPackingProblem& p) : 
	Solver(), problem(p), mu(NULL), pivot(NULL), dimension(NULL), in_bounds(NULL),
    b_encoding(NULL), binaries(NULL) {

    	add_constraints();

}

void OrthogonalPackingSolver::add_constraints(){

    vec<Lit> lits;

    /* Initialisation des prop */
    mu = new int***[problem.k];
    for(int k = 0; k < problem.k; ++k){
        mu[k] = new int**[problem.m];
        for(int a = 0; a < problem.m; ++a){
            mu[k][a] = new int*[problem.n];
            for(int b = 0; b < problem.n; ++b){
                if(problem.is_3d()){
                    mu[k][a][b] = new int[problem.h];
                    for(int c = 0; c < problem.h; ++c){
                        mu[k][a][b][c] = newVar();
                    }    
                }
                else{
                    mu[k][a][b] = new int[1];
                    mu[k][a][b][0] = newVar();
                }
            }
        }
    }

    if(!(problem.orientation == FIX)){
        pivot = new int[problem.k];
        for(int k = 0; k < problem.k; ++k){
            pivot[k] = newVar();
        }
    }

    if(problem.solution_type == SMALLEST){
    	dimension = new int[problem.n - problem.min_n];
    	for(int n = 0; n < problem.n - problem.min_n; ++n){
    		dimension[n] = newVar();
    	}
    	in_bounds = new int****[problem.k];
    	for(int k = 0; k < problem.k; ++k){
	    	in_bounds[k] = new int***[problem.m];
	        for(int a = 0; a < problem.m; ++a){
	            in_bounds[k][a] = new int**[problem.n];
	            for(int b = 0; b < problem.n; ++b){
                    in_bounds[k][a][b] = new int*[1];
                    in_bounds[k][a][b][0]= new int[problem.n - problem.min_n];
                    for(int n = 0; n < problem.n - problem.min_n; ++n){
                    	in_bounds[k][a][b][0][n] = newVar();
                    }

	            }
	        }
    	}
    }

    /* On doit avoir au moins 1 mu à vrai pour un rectangle k donné  */
    for(int k = 0; k < problem.k; ++k){
        lits.clear();
         for(int a = 0; a < problem.m; ++a){
            for(int b = 0; b < problem.n; ++b){
                for(int c = 0; c < (problem.is_3d() ? problem.h : 1); ++c){
                    lits.push(Lit(mu[k][a][b][c]));    
                }
            }
        }
        addClause(lits);
    }

    /* On ne peut avoir 2 mu à vrai en même temps pour un même k, 
    un rectangle k ne pouvant avoir qu'une seule position. Utilisation d'encodage binaire */
    int n_exactly_one = problem.m * problem.n * (problem.is_3d() ? problem.h : 1);
    int k_extactly_one = ceil(log2(n_exactly_one));

    b_encoding = new int[k_extactly_one];
    for(int k = 0; k < k_extactly_one; ++k){
        b_encoding[k] = newVar();
    }

    binaries = new std::bitset<BITSET_SIZE>**[problem.m];
    for(int a = 0; a < problem.m; ++a){
        binaries[a] = new std::bitset<BITSET_SIZE>*[problem.n];
        for(int b = 0; b < problem.n; ++b){
            if(problem.is_3d()){
                binaries[a][b] = new std::bitset<BITSET_SIZE>[problem.h];
                for(int c = 0; c < problem.h; ++c){
                    binaries[a][b][c] = std::bitset<BITSET_SIZE>((a * problem.n * problem.h) + (b * problem.h) + c);
                }
            }
            else{
                binaries[a][b] = new std::bitset<BITSET_SIZE>[1];
                binaries[a][b][0] = std::bitset<BITSET_SIZE>((a * problem.n) + b);
            }
        }
    }

    // for(int k = 0; k < problem.k; ++k){
    //     for(int a = 0; a < problem.m; ++a){
    //         for(int b = 0; b < problem.n; ++b){
    //             for(int c = 0; c < (problem.is_3d() ? problem.h : 1); ++c){
    //                 for(int k_b = 0; k_b < k_extactly_one; ++k_b){
    //                     if(binaries[a][b][c][BITSET_SIZE  - k_extactly_one + k_b] == 1){
    //                         addBinary(~Lit(mu[k][a][b][c]), Lit(b_encoding[k_b]));
    //                     }
    //                     else{
    //                         addBinary(~Lit(mu[k][a][b][c]), ~Lit(b_encoding[k_b]));
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    /* (a, b, c) doit être dans les bornes du grand rectangle */
    for(int k = 0; k < problem.k; ++k){
         for(int a = 0; a < problem.m; ++a){
            for(int b = 0; b < problem.n; ++b){
                for(int c = 0; c < (problem.is_3d() ? problem.h : 1); ++c){
                    if(problem.orientation == FIX){
                    	if(problem.solution_type == SMALLEST){
                    		lits.clear();
                    		lits.push(~Lit(mu[k][a][b][c]));
                    		for(int n = 0; n<problem.n-problem.min_n ; ++n){
                    			if(out_of_bounds(a, b, c, k, n + problem.min_n, n + problem.min_n)){
									addUnit(~Lit(in_bounds[k][a][b][c][n]));
                    			}else{
									addBinary(~Lit(in_bounds[k][a][b][c][n]), Lit(dimension[n]));
								}
								lits.push(Lit(in_bounds[k][a][b][c][n]));
                    		}
                    		addClause(lits);
                    	}else{
                        	if(out_of_bounds(a, b, c, k, problem.n, problem.m)){
                            	addUnit(~Lit(mu[k][a][b][c]));
                        	}
                        }
                    }
                    else{
                        bool out = out_of_bounds(a, b, c, k, problem.n, problem.m);
                        bool pivot_out = pivot_out_of_bounds(a, b, c, k);
                        if(out && pivot_out){
                            addUnit(~Lit(mu[k][a][b][c]));
                        }
                        else if(!out && pivot_out){
                            addBinary(~Lit(mu[k][a][b][c]), ~Lit(pivot[k]));
                        }
                        else if(out && !pivot_out){
                            addBinary(~Lit(mu[k][a][b][c]), Lit(pivot[k])); 
                        }
                    }
                }
            }
        }
    }

    if(problem.is_3d()){
        for(int k = 0; k < problem.k; ++k){
            for(int a = 0; a <= problem.m - problem.lengths[k]; ++a){
                for(int b = 0; b <= problem.n - problem.widths[k]; ++b){
                    for(int c = 0; c <= problem.h - problem.heights[k]; ++c){
                        for(int l = k + 1; l < problem.k; ++l){
                            for(int d = std::max(a - problem.lengths[l] + 1, 0); d < a + problem.lengths[k] && (d <= problem.m - problem.lengths[l]); ++d){ 
                                for(int e = std::max(b - problem.widths[l] + 1, 0); e < b + problem.widths[k] && (e <= problem.n - problem.widths[l]); ++e){
                                    for(int f = std::max(c - problem.heights[l] + 1, 0); f < c + problem.heights[k] && (f <= problem.h - problem.heights[l]); ++f){
                                        addBinary(~Lit(mu[k][a][b][c]), ~Lit(mu[l][d][e][f]));
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    else{
        for(int k = 0; k < problem.k; ++k){
            for(int a = 0; a <= problem.m - problem.lengths[k]; ++a){
                for(int b = 0; b <= problem.n - problem.widths[k]; ++b){
                    for(int l = k + 1; l < problem.k; ++l){
                        if(problem.orientation == PIVOT){
                            for(int d = std::max(a - problem.lengths[l] + 1, 0); d < a + problem.lengths[k] && (d <= problem.m - problem.lengths[l]); ++d){ 
                                for(int e = std::max(b - problem.widths[l] + 1, 0); e < b + problem.widths[k] && (e <= problem.n - problem.widths[l]); ++e){
                                    addTernary(~Lit(mu[k][a][b][0]), ~Lit(mu[l][d][e][0]), Lit(pivot[k]));
                                }
                            }
                            for(int d = std::max(b - problem.lengths[l] + 1, 0); d < b + problem.widths[k] && (d <= problem.m - problem.lengths[l]); ++d){ 
                                for(int e = std::max(a - problem.widths[l] + 1, 0); e < a + problem.lengths[k] && (e <= problem.n - problem.widths[l]); ++e){
                                    addTernary(~Lit(mu[k][b][a][0]), ~Lit(mu[l][d][e][0]), ~Lit(pivot[k]));
                                }
                            }
                        }
                        else{
                            for(int d = std::max(a - problem.lengths[l] + 1, 0); d < a + problem.lengths[k] && (d <= problem.m - problem.lengths[l]); ++d){ 
                                for(int e = std::max(b - problem.widths[l] + 1, 0); e < b + problem.widths[k] && (e <= problem.n - problem.widths[l]); ++e){
                                    addBinary(~Lit(mu[k][a][b][0]), ~Lit(mu[l][d][e][0]));
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    /* Pas de parallélépipède flottant. c commence à 1 car pour c = 0 le rectangle ne flotte pas */
    if(problem.is_3d() && problem.height_constraint == NO_FLOAT){
        for(int k = 0; k < problem.k; ++k){
            for(int a = 0; a < problem.m; ++a){
                for(int b = 0; b < problem.n; ++b){
                    for(int c = 1; c < problem.h; ++c){
                        lits.clear();
                        /* Pas de rectangle porteur encore ajouté donc pour l'instant mu[k][a][b][c] ne peut être valué à 1 */
                        lits.push(~Lit(mu[k][a][b][c]));
                        for(int l = 0; l < problem.k; ++l){
                            if(k == l){ continue; }
                            for(int d = 0; d < problem.m; ++d){
                                for(int e = 0; e < problem.n; ++e){
                                    for(int f = 0; f < problem.h; ++f){
                                        if(carry(a, b, c, d, e, f, k, l)){
                                            /* mu[l][d][e][f] est porteur de mu[k][a][b][c] donc si 
                                            mu[l][d][e][f] est valué à 1, mu[k][a][b][c] peut être à son tour valué à 1 */
                                            lits.push(Lit(mu[l][d][e][f]));
                                        }
                                    }
                                }
                            }
                        }
                        addClause(lits);
                    }
                }
            }
        }
    }
}

/* Inutile de vérifier a >= 0 ni b >= 0 ni c >= 0 car ce sont des indices donc >= 0 par définition */
bool OrthogonalPackingSolver::out_of_bounds(int a, int b, int c, int k, int n, int m){
    return !(  a + problem.lengths[k] <= m 
            && b + problem.widths[k] <= n
            && (problem.is_3d() ? c + problem.heights[k] <= problem.h : true));
}

bool OrthogonalPackingSolver::pivot_out_of_bounds(int a, int b, int c, int k){
    return !(  a + problem.widths[k] <= problem.m 
            && b + problem.lengths[k] <= problem.n
            && (problem.is_3d() ? c + problem.heights[k] <= problem.h : true));
}

bool OrthogonalPackingSolver::overlapping(int a, int b, int c, int d, int e, int f, int k, int l){
    return !(   a + problem.lengths[k] <= d
            ||  a >= d + problem.lengths[l]
            ||  b + problem.widths[k] <= e
            ||  b >= e + problem.widths[l]
            ||  (problem.is_3d() ? c + problem.heights[k] <= f : false)
            ||  (problem.is_3d() ? c >= f + problem.widths[l] : false));
}

bool OrthogonalPackingSolver::pivot_overlapping(int a, int b, int c, int d, int e, int f, int k, int l){
    return !(   a + problem.widths[k] <= d
            ||  a >= d + problem.lengths[l]
            ||  b + problem.lengths[k] <= e
            ||  b >= e + problem.widths[l]
            ||  (problem.is_3d() ? c + problem.heights[k] <= f : false)
            ||  (problem.is_3d() ? c >= f + problem.widths[l] : false));
}

/* Renvoie vrai si (d, e, f) porte (a, b, c) càd si (a, b) et (d, e) sont superposés (en 2D) et f + hauteur == c */
bool OrthogonalPackingSolver::carry(int a, int b, int c, int d, int e, int f, int k, int l){
    return ((c == f + problem.heights[l]) &&
            !(  a + problem.lengths[k] <= d
            ||  a >= d + problem.lengths[l]
            ||  b + problem.widths[k] <= e
            ||  b >= e + problem.widths[l])); 
}

OrthogonalPackingSolution OrthogonalPackingSolver::get_solution(){
    if(mu == NULL){ throw std::runtime_error("Unexpected error : prop vector is null"); }
    if(! okay()){
        return OrthogonalPackingSolution(problem, false);
    }
    else{
        OrthogonalPackingSolution sol(problem, true);
        for(int k = 0; k < problem.k; ++k){
            for(int a = 0; a < problem.m; ++a){
                for(int b = 0; b < problem.n; ++b){
                    for(int c = 0; c < (problem.is_3d() ? problem.h : 1); ++c){
                        if(model[mu[k][a][b][c]] == l_True){
                            sol[k][0] = a; sol[k][1] = b;
                            if(problem.is_3d()){ sol[k][2] = c; }
                        }
                    }
                }
            }
            if(!(problem.orientation == FIX)){ sol.pivot[k] = (model[pivot[k]] == l_True); }
        }
        return sol;
    }
}

OrthogonalPackingSolver::~OrthogonalPackingSolver(){
    if(mu != NULL){
        for(int k = 0; k < problem.k; ++k){
            for(int a = 0; a < problem.m; ++a){
                for(int b = 0; b < problem.n; ++b){
                    delete[] mu[k][a][b];
                }
                delete[] mu[k][a];
            }
            delete[] mu[k];
        }
        delete[] mu;
    }
    if(pivot != NULL){
        delete[] pivot;
    }
    if(dimension != NULL){
    	delete[] dimension;
    }
    if(in_bounds != NULL){
        for(int k = 0; k < problem.k; ++k){
            for(int a = 0; a < problem.m; ++a){
                for(int b = 0; b < problem.n; ++b){
                	delete[] in_bounds[k][a][b][0];
                    delete[] in_bounds[k][a][b];
                }
                delete[] in_bounds[k][a];
            }
            delete[] in_bounds[k];
        }
        delete[] in_bounds;
    }

    if(b_encoding != NULL){
        delete[] b_encoding;    
    }
    if(binaries != NULL){
        for(int a = 0; a < problem.m; ++a){
            for(int b = 0; b < problem.n; ++b){
                delete[] binaries[a][b];
            }
            delete[] binaries[a];
        }
        delete[] binaries;
    }
}