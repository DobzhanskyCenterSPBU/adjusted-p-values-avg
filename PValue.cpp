//
// Created by EL on 09.05.17.
//

#include "PValue.h"
#include <iostream>
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
using namespace std;


vector<double> PValue::adjustPValue(vector<TestsData> const &tests, InputData &G,
                                    vector<unsigned short> const &A, ExecutionParameters const &cont){

    vector<double> P_values((int)(tests.size())); // The results will be saved here
    vector<vector<unsigned short>> cur_G;
    vector<unsigned short> cur_A;
    double D_main, D_cur;
    int s, m, k, all_iter;

    if (cont.isAdaptive) k = cont.k;
    else k = cont.maxReplications;

    for (int i = 0; i < (int)(tests.size()); ++i) {

        //cout << "Test type: " << tests[i].ID << endl;
        //cout << endl;

        prepareData(cur_G, cur_A, G, A, tests[i]);
        D_main = calcPValue(cur_G, cur_A, tests[i].ID);
        //cout << "D_main: " << D_main << endl;
        s = 0; m = 0; all_iter = 0;
        while (s < cont.maxReplications && m < k && all_iter < cont.maxReplications){
            //cout << endl << "ITETATION NUMBER: " << all_iter << endl << endl;
            all_iter++;
            random_shuffle(cur_A.begin(), cur_A.end()); // Create a random permutation of the phenotype values
            D_cur = calcPValue(cur_G, cur_A, tests[i].ID);
            //cout << "D_cur: " << D_cur << endl;
            if (isnan(D_cur)) continue; // If D_cur is NaN go to the next iteration ((D_cur != D_cur))
            s++;
            if (D_cur > D_main) m++;

            //break;
        }
        //cout << endl << "m: " << m << "    s: " << s << endl;
        if (s == 0) P_values[i] = -1;
        else P_values[i] = (double)m/(double)s;
    }
    return P_values;
}

void PValue::prepareData(vector<vector<unsigned short>>& cur_G, vector<unsigned short>& cur_A, InputData & G,
                        vector<unsigned short> const & A, TestsData cur_test){

    cur_A = A; // Create a new copy of the phenotype
    cur_G = G.createGenotypeMatrix(cur_test.lower, cur_test.upper); // Read the appropriate part of the genotype
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
    vector<vector<int>> V;
    int A_car, G_car; // The cardinality of the sets of values of the phenotype and genotype
    int s;                // Parameters for calculating the gamma function
    double chi_sqr = 0.0;
    int row_num = (int)cur_G.size(); // The amount of rows in the current genotype matrix

    A_car = calcNumElem(cur_A); // Make it static?
    if (A_car <= 1) return numeric_limits<double>::quiet_NaN();

    // Define the size of V matrix
    V_cols = A_car; // The number of columns in V corresponds to the amount of unique elements in cur_A
    if (hashIt(ID) == eCD) V_rows = 3;
    else V_rows = 2;


    for (int i = 0; i < row_num; ++i) {

        G_car = calcNumElem(cur_G);
        if (G_car <= 1) continue;

        V = fillVMatrix(cur_G[i], cur_A, V_rows, V_cols);

        // Calculate s and X^2
        s = (A_car - 1)*(G_car - 1);
        chi_sqr = calculateChiSqr(V, V_rows, V_cols);


        // http://www.boost.org/doc/libs/1_49_0/libs/math/doc/sf_and_dist/html/math_toolkit/special/sf_gamma/igamma.html
        // http://keisan.casio.com/exec/system/1180573447
        // tgamma(a,z): Returns the full (non-normalised) upper incomplete gamma function of a and z
        // tgamma_lower(a,z): Returns the full (non-normalised) lower incomplete gamma function of a and z

        D += -log10(1 - boost::math::tgamma(s, 2*chi_sqr)/(boost::math::tgamma(s, 2*chi_sqr) +
                boost::math::tgamma_lower(s, 2*chi_sqr)));
        //cout << "D from function: " << D << endl;
        D_num++;
    }
    //if (D == 0.0 || D_num == 0 || isnan(D)) return numeric_limits<double>::quiet_NaN(); // ! update
    D = D/D_num;

    return D;
}

AlternativeHypothesisType PValue::hashIt (string const& inString) {
    if (inString == "cd") return eCD;
    if (inString == "r")  return eR;
    if (inString == "d")  return eD;
    if (inString == "a")  return eA;
}


int PValue::calcNumElem(vector<unsigned short> const & phenotype){
    vector<unsigned short> elements;
    for(vector<unsigned short>::size_type j = 0; j < phenotype.size(); j++) {
        if (find(elements.begin(), elements.end(), phenotype[j]) == elements.end()) {
            elements.push_back(phenotype[j]);
        }
    }
    return (int)(elements.size());
}

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
        new_phenotype.push_back(pos);
    }
    return new_phenotype;
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

double PValue::calculateChiSqr(vector<vector<int>> V, int V_rows, int V_cols){

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
    //cout << "elem_num: " << elem_num << endl;

    for (int m = 0; m < V_rows; ++m) {
        row_sum = 0;
        for (int j = 0; j < V_cols; ++j) {
            row_sum = row_sum + V[m][j];
        }
        //cout << "row_sum: " << row_sum << endl;
        for (int n = 0; n < V_cols; ++n) {
            col_sum = 0;
            for (int j = 0; j < V_rows; ++j) {
                col_sum = col_sum + V[j][n];
            }
            //cout << "col_sum: " << col_sum << endl;
            if (row_sum*col_sum != 0) {
                chi_sqr = chi_sqr + pow((V[m][n] - ((double)(row_sum*col_sum))/elem_num),2)/
                                    ((double)(row_sum*col_sum)/elem_num);
                //cout << "chi_sqr: " << chi_sqr << endl << endl;
            }
        }
    }
    return chi_sqr;
}