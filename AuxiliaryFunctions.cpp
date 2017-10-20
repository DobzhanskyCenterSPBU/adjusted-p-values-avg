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

