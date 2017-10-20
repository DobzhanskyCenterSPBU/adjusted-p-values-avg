//
// Created by EL on 09.05.17.
//

#include "PValue.h"
#include <iostream>
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
//#include "omp.h"
using namespace std;


vector<double> PValue::adjustPValue(vector<TableEntry> const &tests, InputData &G,
                                    vector<unsigned short> const &A, ExecutionParameters const &cont){

    vector<double> P_values((int)(tests.size())); // The results will be saved here
    vector<vector<unsigned short>> cur_G;
    vector<unsigned short> cur_A;
    double D_main, D_cur;
    int s, m, k, all_iter;

    if (cont.isAdaptive) k = cont.k;
    else k = cont.maxReplications;

    //#pragma omp parallel for
    for (int i = 0; i < (int)(tests.size()); ++i) {
        prepareData(cur_G, cur_A, G, A, tests[i]);
        D_main = calcPValue(cur_G, cur_A, tests[i].ID);
        s = 0; m = 0; all_iter = 0;
        while (s < cont.maxReplications && m < k && all_iter < cont.maxReplications){
            //cout << "ITETATION NUMBER: " << all_iter << endl;
            all_iter++;
            cur_A = phenotypeRandomPermutation(cur_A);
            //random_shuffle(cur_A.begin(), cur_A.end()); // Create a random permutation of the phenotype values
            D_cur = calcPValue(cur_G, cur_A, tests[i].ID);
            if (isnan(D_cur)) continue; // If D_cur is NaN go to the next iteration ((D_cur != D_cur))
            s++;
            if (D_cur > D_main) m++;
        }
        if (s == 0) P_values[i] = -1;
        else P_values[i] = (double)m/(double)s;
    }
    return P_values;
}

void PValue::prepareData(vector<vector<unsigned short>>& cur_G, vector<unsigned short>& cur_A, InputData & G,
                        vector<unsigned short> const & A, TableEntry cur_test){


    if (A.empty()){
        cerr << "Phenotype is empty!" << endl;
        exit(2);
    }
    cur_A = A; // Create a new copy of the phenotype
    cur_G = G.createGenotypeMatrix(cur_test.lower, cur_test.upper); // Read the appropriate part of the genotype
    if (cur_G.empty()){
        cerr << "Genotype file is empty!" << endl;
        exit(2);
    }
    if (hashIt(cur_test.ID) == eA){
        doubleSizeOfMatrices(cur_G, cur_A);
    }

    switch (hashIt(cur_test.ID)) { // Adjust the values of the genotype
        case eCD:
            break;
        case eR:
            for (vector<vector<unsigned short>>::size_type i = 0; i < cur_G.size(); i++) {
                for (vector<unsigned short>::size_type j = 0; j < cur_G[i].size(); j++ ) {
                    if (cur_G[i][j] == 1) cur_G[i][j] = 0;
                    if (cur_G[i][j] == 2) cur_G[i][j] = 1;
                }
            }
            break;
        case eD:
            for (vector<vector<unsigned short>>::size_type i = 0; i < cur_G.size(); i++) {
                for (vector<unsigned short>::size_type j = 0; j < cur_G[i].size(); j++ ) {
                    if (cur_G[i][j] == 2) cur_G[i][j] = 1;
                }
            }
            break;
        case eA:
            for (vector<vector<unsigned short>>::size_type i = 0; i < cur_G.size(); i++) {
                for (vector<unsigned short>::size_type j = 0; j < cur_G[i].size()/2; j++ ) {
                    if (cur_G[i][j] == 1) cur_G[i][j] = 0;
                    if (cur_G[i][j] == 2) cur_G[i][j] = 1;
                }
                for (vector<unsigned short>::size_type j = cur_G[i].size()/2; j < cur_G[i].size(); j++ ) {
                    if (cur_G[i][j] == 2) cur_G[i][j] = 1;
                }
            }
            break;
    }
}

double PValue::calcPValue(vector<vector<unsigned short>> const & cur_G, vector<unsigned short> const & cur_A, string ID){

    double D = 0.0; // Return value
    int D_num = 0; // Amount of valid D values
    int V_rows, V_cols; // Size of V matrix
    int G_car; // The cardinality of the sets of values of the phenotype and genotype
    int s;                // Parameters for calculating the gamma function
    double chi_sqr = 0.0;
    int row_num = (int)cur_G.size(); // The amount of rows in the current genotype matrix
    PhenotypeStatistics phen_stat;

    phen_stat = calcNumElem(cur_A); // Make it static?
    if (phen_stat.num_elements <= 1) return numeric_limits<double>::quiet_NaN();

    vector<vector<int>> V(4, vector<int>(phen_stat.max_element+1)); //Allocate memory for V matrix here, it will always be 4xphen_stat.max_element+1

    for (int i = 0; i < row_num; ++i) {

        if (!checkNumElem(cur_G[i], 3)) continue;
        fillVMatrix(cur_G[i], cur_A, V);
        G_car = calcNumElementsInGenotype(V); // Check this!

        // Calculate s and X^2
        s = (phen_stat.num_elements - 1)*(G_car - 1);
        chi_sqr = calculateChiSqr(V, 3, phen_stat.max_element+1); // Last line correspondong to '3' in genotype not considered

        // http://www.boost.org/doc/libs/1_49_0/libs/math/doc/sf_and_dist/html/math_toolkit/special/sf_gamma/igamma.html
        // http://keisan.casio.com/exec/system/1180573447
        // tgamma(a,z): Returns the full (non-normalised) upper incomplete gamma function of a and z
        // tgamma_lower(a,z): Returns the full (non-normalised) lower incomplete gamma function of a and z

        D += -log10(boost::math::tgamma(s, 2*chi_sqr)/tgamma(s)); // or 1-lower_gamma/full_gamma

        //cout << "D from function: " << D << endl;
        D_num++;
    }
    if (D_num == 0 || isnan(D)) return numeric_limits<double>::quiet_NaN();
    D = D/D_num;

    return D;
}

AlternativeHypothesisType PValue::hashIt (string const& inString) {
    if (inString == "cd") return eCD;
    if (inString == "r")  return eR;
    if (inString == "d")  return eD;
    if (inString == "a")  return eA;
}

PhenotypeStatistics PValue::calcNumElem(vector<unsigned short> const & phenotype){
    PhenotypeStatistics phen_stat = {.max_element = 0, .num_elements = 0};
    vector<unsigned short> elements;
    for(vector<unsigned short>::size_type j = 0; j < phenotype.size(); j++) {
        if (find(elements.begin(), elements.end(), phenotype[j]) == elements.end()) {
            elements.push_back(phenotype[j]);
            if (phenotype[j] > phen_stat.max_element) phen_stat.max_element = phenotype[j];
        }
    }
    phen_stat.num_elements = (int)(elements.size());
    return phen_stat;
}

// Determines if there are at least two distinct elements (not equal to 3) in the genotype vector
// In other words determines if the V matrix can be built and the computations can continue
bool PValue::checkNumElem(vector<unsigned short> const & genotype, unsigned short p) {
    vector<unsigned short> elements;
    for (vector<unsigned short>::size_type j = 0; j < genotype.size(); ++j) {
        if (find(elements.begin(), elements.end(), genotype[j]) == elements.end()) {
            elements.push_back(genotype[j]);
            if (elements.size() == 2 && find(elements.begin(), elements.end(), p) == elements.end() || elements.size() > 2){
                return true;
            }
        }
    }
    return false;
}

void PValue::doubleSizeOfMatrices(vector<vector<unsigned short>> & cur_G, vector<unsigned short> & cur_A){

    vector<unsigned short> buf_A = cur_A, G_line;
    vector<vector<unsigned short>> buf_G = cur_G;

    cur_A = {};
    cur_A.reserve(2*buf_A.size());
    cur_A.insert(cur_A.end(), buf_A.begin(), buf_A.end());
    cur_A.insert(cur_A.end(), buf_A.begin(), buf_A.end());

    // Fill cur_G
    cur_G = {};
    for (vector<vector<unsigned short>>::size_type i = 0; i < buf_G.size(); i++){
        G_line = {};
        G_line.reserve(2*buf_G[i].size());
        G_line.insert(G_line.end(), buf_G[i].begin(), buf_G[i].end());
        G_line.insert(G_line.end(), buf_G[i].begin(), buf_G[i].end());
        cur_G.push_back(G_line);
    }
}

vector<unsigned short> PValue::mapPhenotypeValuesToChar(vector<string> const &phenotype){
    vector<unsigned short> new_phenotype;
    vector<string> elements;
    ptrdiff_t pos;
    for (std::vector<string>::const_iterator it=phenotype.begin(); it!=phenotype.end(); ++it){
        if (find(elements.begin(), elements.end(), *it) == elements.end()) {
            elements.push_back(*it);
        }
        pos = find(elements.begin(), elements.end(), *it) - elements.begin();
        new_phenotype.push_back((unsigned short)pos);
    }
    return new_phenotype;
}

void PValue::fillVMatrix(vector<unsigned short> const & cur_G, vector<unsigned short> const & cur_A,
                                        vector<vector<int>> & V){

    int V_rows = (int)V.size();
    int V_cols = (int)V[0].size();
    int col_num = (int)cur_A.size(); // The amount of patients

    // Fill V with zeros
    for(vector<int>::size_type m = 0; m < V_rows; m++) {
        for(vector<int>::size_type n = 0; n < V_cols; n++){
            V[m][n] = 0;
        }
    }
    // Fill V (G_car x A_car)
    for(vector<unsigned char>::size_type j = 0; j < col_num; j++) {
        V[cur_G[j]][cur_A[j]] = V[cur_G[j]][cur_A[j]] + 1;
    }
}

double PValue::calculateChiSqr(vector<vector<int>> const & V, int V_rows, int V_cols){

    double chi_sqr = 0.0;
    int elem_num; // n in X^2
    int row_sum, col_sum;

    // Calc elem_num
    elem_num = 0;
    for (int m = 0; m < V_rows; ++m) {
        for (int n = 0; n < V_cols; ++n) {
            elem_num += V[m][n];
        }
    }

    for (int m = 0; m < V_rows; ++m) {
        row_sum = 0;
        for (int j = 0; j < V_cols; ++j) {
            row_sum = row_sum + V[m][j];
        }
        for (int n = 0; n < V_cols; ++n) {
            col_sum = 0;
            for (int j = 0; j < V_rows; ++j) {
                col_sum = col_sum + V[j][n];
            }
            if (row_sum*col_sum != 0) {
                chi_sqr = chi_sqr + pow((V[m][n] - ((double)(row_sum*col_sum))/elem_num),2)/
                                    ((double)(row_sum*col_sum)/elem_num);
            }
        }
    }
    return chi_sqr;
}

// Random generation function:
int PValue::myRandom(int i) {
    srand (time(NULL)); //Seeding the random generator
    return rand() % i;
}

// Permutate phenotype
vector<unsigned short> PValue::phenotypeRandomPermutation(vector<unsigned short>& phenotype) {
    random_shuffle(phenotype.begin(), phenotype.end(), myRandom); // Using myRandom(int i)
    return phenotype;
}

// Calculates the exact number of elements in a genotype vector based on the V matrix
int PValue::calcNumElementsInGenotype(vector<vector<int>> V){
    int num_elem = 0;
    int sum = 0;
    for (vector<vector<int>>::size_type i = 0; i < V.size() - 1; i++) {
        sum = 0;
        for (vector<int>::size_type j = 0; j < V[i].size(); j++) {
            sum += V[i][j];
            if (sum > 0) {num_elem++; break;}
        }
    }
    return num_elem;
}


/* Older functions

int PValue::calcNumElem(vector<vector<unsigned short>> const & genotype) {
    unsigned short p = 3; // Element that is not considered
    vector<unsigned short> elements;
    for (vector<vector<unsigned short>>::size_type i = 0; i < genotype.size(); ++i) {
        for (vector<unsigned short>::size_type j = 0; j < genotype[i].size(); ++j) {
            if (find(elements.begin(), elements.end(), genotype[i][j]) == elements.end()) {
                elements.push_back(genotype[i][j]);
            }
        }
    }
    if (find(elements.begin(), elements.end(), p) != elements.end()) return (int)(elements.size() - 1);
    else return (int)(elements.size());
}


 vector<vector<int>> PValue::fillVMatrix(vector<unsigned short> const & cur_G, vector<unsigned short> const & cur_A,
                                int V_rows, int V_cols){

    vector<vector<int>> V (V_rows, vector<int>(V_cols));
    int col_num = (int)cur_A.size(); // The amount of patients

    // Fill V with zeros
    for(vector<int>::size_type m = 0; m < V_rows; m++) {
        for(vector<int>::size_type n = 0; n < V_cols; n++){
            V[m][n] = 0;
        }
    }
    // Fill V (G_car x A_car)
    for(vector<unsigned char>::size_type j = 0; j < col_num; j++) {
        if (cur_G[j] < 3) {
            V[cur_G[j]][cur_A[j]] = V[cur_G[j]][cur_A[j]] + 1; // " - '0' " - transforms char to int
        }
    }
    return V;
}


 int PValue::calcNumElem(vector<unsigned short> const & genotype, unsigned short p) {
    //unsigned short p = 3; // Element that is not considered
    vector<unsigned short> elements;
        for (vector<unsigned short>::size_type j = 0; j < genotype.size(); ++j) {
            if (find(elements.begin(), elements.end(), genotype[j]) == elements.end()) {
                elements.push_back(genotype[j]);
            }
        }
    if (find(elements.begin(), elements.end(), p) != elements.end()) return (int)(elements.size() - 1);
    else return (int)(elements.size());
}

     // Define the size of V matrix
    //V_cols = A_car; // The number of columns in V corresponds to the amount of unique elements in cur_A
    //if (hashIt(ID) == eCD) V_rows = 3;
    //else V_rows = 2;
*/