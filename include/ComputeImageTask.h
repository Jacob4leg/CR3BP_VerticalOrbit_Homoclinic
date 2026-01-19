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
    Interval X, Z, DY;
    Interval E_remainder;
    IVector v0, V;
    IMatrix L, L_inv;
    IVector x, result;
    capd::IMap dy;
    capd::IMap E_map;
    Status status = NEW;

    ComputeImageTask(Interval theta, Interval E, IVector V, IMatrix L) : theta(theta), E(E), V(V), L(L) {
        result = IVector(2);
        split(v0,V);
        X = V[0]; Z = V[1]; DY = V[4];
        
        dy = capd::IMap(dy_map,6,1,2);
        dy.setParameter(0,muSJ);
        dy.setParameter(1,muSJ-1.);

        E_map = capd::IMap(energy,6,1,2);
        E_map.setParameter(0,muSJ);
        E_map.setParameter(1,muSJ-1.);

        IMatrix L_temp = capd::matrixAlgorithms::krawczykInverse(matrix_erase_cord(L,1));
        L_inv = matrix_add_cord(L_temp,1);
        E.split(E_remainder);

    }

    IMatrix get_energy_change_of_basis() {
        // IVector U = v0;
        // U += IVector{X_remainder,0,0,0,DY_remainder,0};
        IVector dE = E_map.derivative(V)[0];
        IVector z_grad = -IVector{dE[0],dE[1],-1,dE[3],dE[4],dE[5]} / dE[2];
        IMatrix T = IMatrix::Identity(6);
        T.row(2) = z_grad;
        return T;
    }

    bool prove_Z_bound() {
        // TODO
    }


    bool prove_fixed_point(capd::IPoincareMap& pm_z) {
        // TODO
        // capd::IVector u0 = v_E;
        // capd::IVector U0 = v_E;
        // capd::IVector U = u0 + capd::IVector{X_remainder,0,Z_remainder,0,DY_remainder,0};
        // U0[4] += DY_remainder;

        // capd::IMatrix D(6,6);
        // capd::C1HORect2Set S(U);
        // capd::IVector V = pm_z(S,D);
        // D = pm_z.computeDP(V,D);

        // capd::IMatrix M({{D[1][0],D[1][2]},{D[3][0],D[3][2]}});

        // capd::C0HORect2Set S0(U0);
        // V = pm_z(S0);
        // capd::IVector V0{V[1],V[3]};
        // capd::IVector N = - capd::matrixAlgorithms::gauss(M,V0);
        
        // return subset(N, capd::IVector{X_remainder,Z_remainder});
    }


    Interval cone_coeff(int tau_div, IMatrix T_total, IMatrix T_total_inv, capd::IPoincareMap& pm_y) {
        double tau_eps = 1e-7;
        double gamma_eps = 1e-8;

        double tau_delta = tau_eps / tau_div;

        IMatrix D_total(5,5);

        Interval gamma = gamma_eps * Interval(-1,1);
        Interval tau_i(0, tau_delta);

        for(int i = 0; i < tau_div; i++) {
            IVector u{tau_i,0, E, gamma, gamma, gamma};
            capd::C1HORect2Set S(capd::C1Rect2Set::C0BaseSet(v0, T_total, u), capd::C1Rect2Set::C1BaseSet(T_total));
            
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


        std::cout << D_total << std::endl;
        IMatrix H = matrix_erase_cord(D_total,1);

        // std::cout << H << std::endl;

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
            
            
            IMatrix T = get_energy_change_of_basis();

            IMatrix T_inv = capd::matrixAlgorithms::krawczykInverse(T);
            IMatrix T_total = T * L;
            IMatrix T_total_inv = L_inv * T_inv;

            Interval alpha = cone_coeff(20,T_total, T_total_inv, pm_y);

            Interval gamma = alpha * theta * Interval(-1,1);

            IVector w{theta,0,0,gamma,gamma,gamma};
            w += T_total_inv * V;
            

            capd::C1HORect2Set S(capd::C1Rect2Set::C0BaseSet(v0,T_total,w),capd::C1Rect2Set::C1BaseSet(T_total));
            pm_x(S);
            IVector res = pm_y(S);
            result = IVector{res[3],res[5]};
            

            status = COMPLETED;
        } catch(...) {
            status = FAILED;
        }
    }

    
};

#endif