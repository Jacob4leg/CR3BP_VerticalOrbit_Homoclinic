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

    CR3BP<long double> vf;
    Interval theta, E;
    IVector v_E;
    IMatrix L, L_inv;
    IVector x, result;
    capd::IMap dy;
    Status status = NEW;

    ComputeImageTask(Interval theta, Interval E, IVector v_E, IMatrix L) : theta(theta), E(E), v_E(v_E), L(L) {
        dy = capd::IMap(dy_map,6,1,2);
        L_inv = capd::matrixAlgorithms::gaussInverseMatrix(L);
    }

    bool prove_fixed_point(capd::Interval X, capd::Interval Z, capd::Interval DY, capd::IPoincareMap& pm_z) {
        capd::IVector u0 = v_E;
        capd::IVector U0 = v_E;
        capd::IVector U = u0 + capd::IVector{X,0,Z,0,DY,0};
        U0[4] += DY;

        capd::IMatrix D(6,6);
        capd::C0HORect2Set S(U);
        capd::IVector V = pm_z(S,D);
        D = pm_z.computeDP(V,D);

        capd::IMatrix M({{D[1][0],D[1][2]},{D[3][0],D[3][2]}});

        capd::C0HORect2Set S0(U0);
        V = pm_z(S0);
        capd::IVector V0{V[1],V[3]};

        capd::IVector N = - capd::matrixAlgorithms::gauss(M,V0);

        return subset(N, capd::IVector{X,Z});
    }


    Interval cone_coeff(int tau_div, int E_div, capd::IPoincareMap& pm_y) {
        double tau_eps = 1e-7;
        double gamma_eps = 1e-8;

        double tau_delta = tau_eps / tau_div;
        double E_delta = (E.rightBound() - E.leftBound()) / E_div;

        Interval E_i(E.leftBound(), E.leftBound() + E_delta);

        IMatrix D_total(5,5);

        for(int i = 0; i < tau_div; i++) {
            Interval tau_j(0,tau_delta);
            
            for(int j = 0; j < E_div; j++) {
                IVector u{tau_j, 0, gamma_eps, gamma_eps, E_i, gamma_eps};
                capd::C1HORect2Set S(capd::C1Rect2Set::C0BaseSet(v_E, L, u), capd::C1Rect2Set::C1BaseSet(L));
                IMatrix D(6,6);
                IVector y = pm_y(S,D,2);
                D = L_inv * D;
                IMatrix D_reduced = matrix_erase_cord(D,1);

                if(i == 0 && j == 0) D_total = D_reduced;
                else D_total = intervalHull(D_reduced,D_total);
            }
        }

        return Interval(0.);
    }
    

    void run(unsigned threadID) {
        try {
            capd::IPoincareMap& pm_x = CR3BPPool::pool_X->getSolver(threadID).pm;
            capd::IPoincareMap& pm_y = CR3BPPool::pool_Y->getSolver(threadID).pm;
            capd::IPoincareMap& pm_z = CR3BPPool::pool_Z->getSolver(threadID).pm;
            
            IVector w{theta,0,0,0,E,0};

            // std::cout << "computing image of " << x << " in thread " << threadID << std::endl << std::flush;
            capd::C0HORect2Set S(v_E, L, w);
            pm_x(S);
            result = pm_y(S);
            status = COMPLETED;
        } catch(...) {
            status = FAILED;
        }
    }

    
};

#endif