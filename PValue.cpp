//
// Created by EL on 09.05.17.
//

#include "PValue.h"
#include <iostream>
#include <cmath>
#include <boost/math/special_functions/gamma.hpp>
using namespace std;


vector<double> PValue::adjustPValue(TestsData const & tests, InputData & G,
                                    vector<char> const & A, ExecutionParameters const & cont){

    vector<double> P_values(tests.ID.size()); // The results will be saved here
    vector<char> cur_G;
    vector<char> cur_A;
    double D_main, D_cur;
    int s, m, k, all_iter;

    if (cont.isAdaptive) k = cont.k;
    else k = cont.maxReplications;

    for (int i = 0; i < tests.ID.size(); ++i) {
        prepareData(cur_G, cur_A, G, A, tests.lower[i], tests.upper[i], tests.ID[i]);

        D_main = calcPValue(cur_G, cur_A);
        cout << "D_main: " << D_main << endl;

        s = 0; m = 0; all_iter = 0;
        while (s <= cont.maxReplications && m <= k && all_iter <= cont.maxReplications){
            all_iter++;
            random_shuffle(cur_A.begin(), cur_A.end()); // Create a random permutation of the phenotype values
            D_cur = calcPValue(cur_G, cur_A);
            cout << "D_cur: " << D_cur << endl;
            if (D_cur != D_cur) continue; // If D_cur is NaN go to the next iteration
            s++;
            if (D_cur > D_main) m++;
        }
        P_values[i] = m/s;
    }
    return P_values;
}

int PValue::prepareData(vector<char>& cur_G, vector<char>& cur_A, InputData & G,
                        vector<char> const & A, int low, int up, string id){

    vector<char> buf_G;
    if (hashIt(id) != eA){
        cur_A = A; // Create a new copy of the phenotype
        cur_G = G.createGenotypeVector(low,up); // Read the appropriate part of the genotype
    }
    else {
        cur_A = {};
        cur_A.reserve(2*A.size());
        cur_A.insert(cur_A.end(), A.begin(), A.end());
        cur_A.insert(cur_A.end(), A.begin(), A.end());
        cur_G = {};
        buf_G = G.createGenotypeVector(low,up); // Read the appropriate part of the genotype
        cur_G.reserve(2*buf_G.size());
        cur_G.insert(cur_G.end(), buf_G.begin(), buf_G.end());
        cur_G.insert(cur_G.end(), buf_G.begin(), buf_G.end());
    }

    switch (hashIt(id)) { // Adjust the values of the genotype
        case eCD:
            return 0;
        case eR:
            for (vector<char>::iterator it = cur_G.begin(); it != cur_G.end(); ++it){
                if (*it == '1') *it = '0';
                if (*it == '2') *it = '1';
            }
            break;
        case eD:
            for (vector<char>::iterator it = cur_G.begin(); it != cur_G.end(); ++it){
                if (*it == '2') *it = '1';
            }
            break;
        case eA:
            for(vector<char>::size_type i = 0; i < cur_G.size()/2; i++) {
                if (cur_G[i] == '1') cur_G[i] = '0';
                if (cur_G[i] == '2') cur_G[i] = '1';
            }
            for(vector<char>::size_type i = cur_G.size()/2; i < cur_G.size(); i++) {
                if (cur_G[i] == '2') cur_G[i] = '1';
            }
            break;
        default:
            return 5; // Signal a problem
    }
}

double PValue::calcPValue(vector<char> const & cur_G, vector<char> const & cur_A){

    double D = 0.0; // Return value
    int D_num = 0; // Amount of valid D values
    int V_rows = 3, V_cols = 4; // Size of V matrix
    int A_car, G_car; // The cardinality of the sets of values of the phenotype and genotype
    int s;           // Parameters for calculating the incomplete gamma function
    double chi_sqr;
    int row_sum, col_sum;
    int col_num = cur_A.size(); // The amount of patients
    int row_num = cur_G.size()/cur_A.size(); // The amount pf rows in the current genotype matrix
    int elem_num; // n in X^2

    A_car = calcNumElem(cur_A.begin(), cur_A.end());
    if (A_car <= 1) return numeric_limits<double>::quiet_NaN();



    for (int i = 0; i < row_num; ++i) {
        vector<vector<int>> V (V_rows, vector<int>(V_cols));

        G_car = calcNumElem(cur_G.begin()+i*col_num, cur_G.begin()+(i+1)*col_num, '3');
        if (G_car <= 1) continue;

        // Fill V with zeros
        for(vector<char>::size_type m = 0; m < V_rows; m++) {
            for(vector<char>::size_type n = 0; n < V_cols; n++){
                V[m][n] = 0;
            }
        }
        // Fill V (G_car x A_car)
        for(vector<char>::size_type j = 0; j < col_num; j++) {
            if (cur_G[i*col_num+j] != '3'){
                V[cur_G[i*col_num+j] - '0'][cur_A[j] - '0']++;
            }
        }
        // Calc elem_num
        elem_num = 0;
        for (int m = 0; m < V_rows; ++m) {
            for (int n = 0; n < V_cols; ++n) {
                elem_num += V[m][n];
            }
        }
        //cout << "elem_num: " << elem_num << endl;
        // Calculate s and X^2
        s = (A_car - 1)*(G_car - 1);
        chi_sqr = 0.0;
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

        cout << "boost::math::tgamma_lower(s, 2*chi_sqr):" << boost::math::tgamma_lower(s, 2*chi_sqr) << endl;
        D += -log10(1 - boost::math::tgamma_lower(s, 2*chi_sqr));
        cout << "D from function: " << D << endl;
        D_num++;
    }
    //if (D == 0.0 || D_num == 0 || D != D) return numeric_limits<double>::quiet_NaN(); // ! update
    D = D/D_num;

    return D;
}

AlternativeHypothesisType PValue::hashIt (string const& inString) {
    if (inString == "cd") return eCD;
    if (inString == "r")  return eR;
    if (inString == "d")  return eD;
    if (inString == "a")  return eA;
}

int PValue::calcNumElem(vector<char>::const_iterator first, vector<char>::const_iterator last, char p) {
    vector<char> elements;
    for (vector<char>::const_iterator it = first; it != last; ++it) {
        if (find(elements.begin(), elements.end(), *it) == elements.end()) {
            elements.push_back(*it);
        }
    }
    if (find(elements.begin(), elements.end(), p) != elements.end()) return elements.size() - 1;
    else return elements.size();
}