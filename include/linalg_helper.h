#include<iostream>
#include<tuple>
#include "capd/capdlib.h"


#ifndef _LINALG_HELPER_H_
#define _LINALG_HELPER_H_

template<typename T>
inline capd::vectalg::Matrix<T,0,0> matrix_erase_cord(capd::vectalg::Matrix<T,0,0> M, int k) {
    auto dim = M.dimension();
    int n = std::get<0>(dim);
    assert(n == std::get<1>(dim));

    capd::vectalg::Matrix<T,0,0> M_res(n-1,n-1);

    for(int i = 0; i < n-1; i++) {
        int l = i < k ? i : i+1;
        for(int j = 0; j < n-1; j++) {
            int m = j < k ? j : j+1;
            M_res[i][j] = M[l][m];
        }
    }
    return M_res;
}

template<typename T>
inline capd::vectalg::Matrix<T,0,0> matrix_add_cord(capd::vectalg::Matrix<T,0,0> M, int k) {
    auto dim = M.dimension();
    int n = std::get<0>(dim);
    assert(n == std::get<1>(dim));

    capd::vectalg::Matrix<T,0,0> M_res(n+1,n+1);

    for(int i = 0; i < n+1; i++) {
        if(i == k) continue;
        int l = i < k ? i : i-1;
        for(int j = 0; j < n+1; j++) {
            if(j == k) continue;
            int m = j < k ? j : j-1;
            M_res[i][j] = M[l][m];
        }
    }
    return M_res;
}

inline std::tuple<long double, capd::LDVector> get_dominant_eigenvalue_and_eigenvector(capd::LDMatrix A) {
    int n = A.numberOfRows();
    assert(n == A.numberOfColumns());

    capd::LDVector rV(n), iV(n);
    capd::LDMatrix rVec(n,n), iVec(n,n);
    capd::alglib::computeEigenvaluesAndEigenvectors(A,rV,iV,rVec,iVec);

    int max_index = 0;
    for(int i = 1; i < n; i++) {
        long double module_i = rV[i] * rV[i] + iV[i] * iV[i];
        long double module_max_index = rV[max_index] * rV[max_index] + iV[max_index] * iV[max_index];

        if(module_i > module_max_index) max_index = i;
    }
    assert(iV[max_index] == 0);
    long double lambda = rV[max_index];
    capd::LDVector u = rVec.column(max_index);
    return {lambda, u};
}

inline capd::LDMatrix us_change_of_basis(capd::LDMatrix A) {
    std::pair<int,int> dim = A.dimension();
    int n = std::get<0>(dim);
    int m = std::get<1>(dim);
    assert(n == m);

    capd::LDVector rV(n), iV(n);
    capd::LDMatrix rVec(n,n), iVec(n,n);
    capd::alglib::computeEigenvaluesAndEigenvectors(A,rV,iV,rVec,iVec);

    for(int i = 0; i < n; i++) {
        if(iV[i] != 0) {
            rVec.column(i+1) = iVec.column(i);
            i++;
        }
    }
    return rVec;
}


#endif