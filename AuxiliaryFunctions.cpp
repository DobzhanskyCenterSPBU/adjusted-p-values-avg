//
// Created by doby on 13/10/17.
//

#include <iostream>
#include <vector>
#include "PValue.h"
#include <fstream>

using namespace std;

template <typename T>
void printVector(vector<T> vec){
    for (typename vector<T>::size_type i = 0; i < vec.size(); i++) {
        cout << vec[i] << ' ';
    }
    cout << '\n';
};

template <typename T>
void printVectorOfVectors(vector<vector<T>> vec){
    for (typename vector<vector<T>>::size_type i = 0; i < vec.size(); i++) {
        for (typename vector<T>::size_type j = 0; j < vec[i].size(); j++ ) {
            cout << vec[i][j] << ' ';
        }
        cout << '\n';
    }
    cout << '\n';
};

vector<TableEntry> readTopHitsTable(string filename, int number_of_lines) {
    ifstream handle;
    handle.open(filename);
    if (!handle.is_open()) {
        cerr << "No such file!" << endl;
        exit(1); // What exit code should be used? How to skip to the next iteration?
    }

    TableEntry currentEntry;
    vector<TableEntry> result;
    string data, rank, test_index, chr_index, snp_ind_left, snp_ind_center, snp_ind_right, lower_pos, center_pos,
            upper_pos, count, neg_log_dens, dens_p_value, density;

    getline(handle, data); // Skip the first line
    for (int i = 0; i < number_of_lines; ++i) // What other way is there? while(handle.good()){}
    {
        getline(handle, rank, ';');
        getline(handle, test_index, ';'); //!
        getline(handle, chr_index, ';'); //!
        getline(handle, snp_ind_left, ';');
        getline(handle, snp_ind_center, ';');
        getline(handle, snp_ind_right, ';');
        getline(handle, lower_pos, ';'); //!
        getline(handle, center_pos, ';');
        getline(handle, upper_pos, ';'); //!
        getline(handle, count, ';');
        getline(handle, neg_log_dens, ';');
        getline(handle, dens_p_value, ';'); //!
        getline(handle, density);

        // double naive_p_value = boost::lexical_cast<double>(dens_p_value); // ???????????????????

        currentEntry = {.test_index = test_index, .chromosome_index = chr_index, .ID = "", .lower = stoi(lower_pos),
                .upper = stoi(upper_pos), .adjusted_p_value = -2};
        result.push_back(currentEntry); // Is int enough for index ?????
    }

    return result;
}
