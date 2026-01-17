#ifndef _COMPUTE_IMAGE_TASK_H_
#define _COMPUTE_IMAGE_TASK_H_

#include<iostream>
#include "capd/capdlib.h"
#include "cr3bp.h"
#include "linalg_helper.h"



struct ComputeImageTask : public capd::threading::Task {
    typedef capd::Interval Interval;
    typedef capd::LDVector LDVector; typedef capd::IVector IVector;
    typedef capd::LDMatrix LDMatrix; typedef capd::IMatrix IMatrix;

    enum Status {NEW, COMPLETED, FAILED};

    long double muSJ = 0.00095388114032796904;
    CR3BP<long double> vf;
    Interval theta, E;
    Interval E_remainder;
    IVector v_E;
    IMatrix L, L_inv;
    IVector x, result;
    capd::IMap dy;
    Status status = NEW;

    ComputeImageTask(Interval theta, Interval E, IVector v_E, IMatrix L) : theta(theta), E(E), v_E(v_E), L(L) {
        dy = capd::IMap(dy_map,6,1,2);
        dy.setParameter(0,muSJ);
        dy.setParameter(1,muSJ-1.);

        IMatrix L_temp = capd::matrixAlgorithms::gaussInverseMatrix(matrix_erase_cord(L,1));
        L_inv = matrix_add_cord(L_temp,1);
        E.split(E_remainder);

    }

    IMatrix get_energy_change_of_basis(Interval X_remainder, Interval Z_remainder) {
        IVector U = v_E;
        U[4] = E;
        U += IVector{X_remainder,0,Z_remainder,0,0,0};
        IVector grad = dy.derivative(U)[0];
        IMatrix T = IMatrix::Identity(6);
        T.row(4) = grad;
        return T;
    }

    bool prove_fixed_point(capd::Interval X_remainder, capd::Interval Z_remainder, capd::Interval DY_remainder, capd::IPoincareMap& pm_z) {
        
        capd::IVector u0 = v_E;
        capd::IVector U0 = v_E;
        capd::IVector U = u0 + capd::IVector{X_remainder,0,Z_remainder,0,DY_remainder,0};
        U0[4] += DY_remainder;

        capd::IMatrix D(6,6);
        capd::C1HORect2Set S(U);
        capd::IVector V = pm_z(S,D);
        D = pm_z.computeDP(V,D);

        capd::IMatrix M({{D[1][0],D[1][2]},{D[3][0],D[3][2]}});

        capd::C0HORect2Set S0(U0);
        V = pm_z(S0);
        capd::IVector V0{V[1],V[3]};
        capd::IVector N = - capd::matrixAlgorithms::gauss(M,V0);
        
        return subset(N, capd::IVector{X_remainder,Z_remainder});
    }


    Interval cone_coeff(int tau_div, IMatrix T_total, IMatrix T_total_inv, capd::IPoincareMap& pm_y) {
        double tau_eps = 1e-7;
        double gamma_eps = 1.4e-7;

        // tau_eps = 0.;
        // gamma_eps = 0.;

        double tau_delta = tau_eps / tau_div;

        IMatrix D_total(5,5);

        Interval gamma = gamma_eps * Interval(-1,1);
        Interval tau_i(0, tau_delta);

        for(int i = 0; i < tau_div; i++) {
            // std::cout << tau_i << std::endl;
            IVector u{tau_i,0, gamma, gamma, E, gamma};
            capd::C1HORect2Set S(capd::C1Rect2Set::C0BaseSet(v_E, T_total, u), capd::C1Rect2Set::C1BaseSet(T_total));
            
            IMatrix D(6,6);
            IVector y = pm_y(S,D,2);
            D = pm_y.computeDP(y,D);
            
            D = T_total_inv * D;
            IMatrix D_reduced = matrix_erase_cord(D,1);

            if(i == 0) D_total = D_reduced;
            else D_total = intervalHull(D_reduced,D_total);

            tau_i += tau_delta;
        }
        capd::IMaxNorm max_norm;
        capd::ISumNorm sum_norm;

        IMatrix H = matrix_erase_cord(D_total,3);

        Interval A = capd::abs(D_total[0][0]).left();
        Interval B = sum_norm(IVector{H[0][1],H[0][2],H[0][3]}).right();
        Interval C = max_norm(IVector{H[1][0],H[2][0],H[3][0]}).right();

        Interval K = max_norm(matrix_erase_cord(H,0)).right();

        Interval delta = power(-A + K,2) - 4 * B * C;
        
        Interval alpha1 = (A - K + sqrt(delta)) / (2 * B);
        Interval alpha2 = (A - K - sqrt(delta)) / (2 * B);

        std::cout << alpha1 << std::endl;
        std::cout << alpha2 << std::endl;

        return alpha2.right();
    }
    

    void run(unsigned threadID) {
        try {
            capd::IPoincareMap& pm_x = CR3BPPool::pool_X->getSolver(threadID).pm;
            capd::IPoincareMap& pm_y = CR3BPPool::pool_Y->getSolver(threadID).pm;
            capd::IPoincareMap& pm_z = CR3BPPool::pool_Z->getSolver(threadID).pm;
            
            IVector w{theta,0,0,0,E,0};

            Interval X_test = 1e-12 * Interval(-1,1);
            Interval Z_test = 1e-12 * Interval(-1,1);
            E = 1e-15 * Interval(-1,1);
            
            IMatrix T = get_energy_change_of_basis(X_test,Z_test);

            

            IMatrix T_inv = capd::matrixAlgorithms::gaussInverseMatrix(T);
            IMatrix T_total = T * L;
            IMatrix T_total_inv = L_inv * T_inv;
            // std::cout << T_total << std::endl;
            // std::cout << T_total_inv << std::endl;
            // return;

            // IMatrix D(6,6);
            // capd::C1HORect2Set S(v_E);
            // IVector V = pm_y(S,D,2);
            // D = pm_y.computeDP(V,D);
            // D = T_total_inv * D * T_total;
            // std::cout << D << std::endl;
            // return;
            

            cone_coeff(50,T_total, T_total_inv, pm_y);

            // std::cout << prove_fixed_point(X_test, Z_test, E, pm_z) << std::endl;

            // std::cout << "computing image of " << x << " in thread " << threadID << std::endl << std::flush;
            // capd::C0HORect2Set S(v_E, L, w);
            // pm_x(S);
            // result = pm_y(S);
            status = COMPLETED;
        } catch(...) {
            status = FAILED;
        }
    }

    
};

#endif