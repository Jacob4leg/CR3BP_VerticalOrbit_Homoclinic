#ifndef _COMPUTE_IMAGE_TASK_H_
#define _COMPUTE_IMAGE_TASK_H_

#include<iostream>
#include "capd/capdlib.h"
#include "cr3bp.h"
#include "linalg_helper.h"
#include "InitConditionsGenerator.h"



struct ComputeImageTask : public capd::threading::Task {
    typedef capd::Interval Interval;
    typedef capd::LDVector LDVector; typedef capd::IVector IVector;
    typedef capd::LDMatrix LDMatrix; typedef capd::IMatrix IMatrix;

    enum Status {NEW, COMPLETED, FAILED};

    long double muSJ = 0.00095388114032796904;
    CR3BP<long double> vf;
    Interval tau, E;
    Interval X, Z, DY;
    Interval X0, X_remainder;
    Interval Z0, Z_remainder;
    Interval DY0, DY_remainder;
    Interval E0,E_remainder;
    IVector V;
    IMatrix L, L_inv, C;
    IVector x, result;
    capd::IMap dy;
    capd::IMap E_map;
    // InitConditionsGenerator& gen;
    Status status = NEW;

    ComputeImageTask(Interval tau, Interval E, IVector V, IMatrix C, InitConditionsGenerator& gen) : 
            tau(tau), E(E), V(V),C(C) {
        
        // *gen_ptr = &gen;

        result = IVector(2);
        X = V[0]; Z = V[2]; DY = V[4];
  
        dy = capd::IMap(dy_map,6,1,2);
        dy.setParameter(0,muSJ);
        dy.setParameter(1,muSJ-1.);

        E_map = capd::IMap(energy,6,1,2);
        E_map.setParameter(0,muSJ);
        E_map.setParameter(1,muSJ-1.);

        
        E.split(E0,E_remainder);
        X.split(X0,X_remainder); Z.split(Z0,Z_remainder); DY.split(DY0,DY_remainder);

        IMatrix L1 = gen.get_change_of_basis(E.leftBound());
        IMatrix L2 = gen.get_change_of_basis(E.rightBound());
        L = intervalHull(L1,L2);
        // C = capd::matrixAlgorithms::krawczykInverse(B);
        
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
        V[2] += 7 * Z_remainder;
        IVector v0 = V;
        IVector U(6);
        split(v0,U);
        v0 = V; v0[2] = V[2].mid();

        IVector E_grad = E_map.derivative(V)[0];
        Interval p0 = E_map(v0)[0] - E;

        Interval N = - p0 / E_grad[2];
        
        // std::cout << N << std::endl;
        // std::cout << U[2] << std::endl;
        // std::cout << "/////////////////////" << std::endl;

        return subset(N,U[2]);
        
        // std::cout << p0 << std::endl;
        


        // LDVector LD_u0{X0.leftBound(),0,Z0.leftBound(),0,DY0.leftBound(),0};
        // IVector u0{X0,0,Z0,0,DY0,0};
        // IVector U{X,0,Z0,0,DY0,0};
        
        // std::cout << (E_map(V) - E) << std::endl;

        // capd::LDMap LD_E_map(energy,6,1,2);
        // LD_E_map.setParameter(0,muSJ);
        // LD_E_map.setParameter(1,muSJ-1.);

        // LDVector LD_dE = LD_E_map.derivative(LD_u0)[0];
        
        // IMatrix A({{Interval(LD_dE[0]),-1,Interval(LD_dE[4])}});
        // A = A / Interval(LD_dE[2]);
        

        // IVector dE_point = E_map.derivative(u0)[0];
        // IVector dE = E_map.derivative(U)[0];

        // IMatrix Df_x({{dE[0],-1,dE[4]}});
        // IMatrix Df_z({{dE[2]}});

        // IMatrix Dg_x = Df_x - capd::vectalg::intersection(Df_z * A, (Df_z / dE_point[2]) * IMatrix({{dE_point[0],-1,dE_point[4]}}));
        
        // IVector g = (E_map(u0) - E0) + Dg_x * IVector{X_remainder,E_remainder,DY_remainder};
        // dE = E_map.derivative(V)[0];
        // IVector N = - g / dE[2];
        // std::cout << N << std::endl;
        // std::cout << Z_remainder + (A * IVector{X_remainder,E_remainder,DY_remainder})[0] << std::endl;
        // std::cout << "//////////////////////////////" << std::endl;

        
    }


    bool prove_fixed_point(capd::IPoincareMap& pm_z) {
        V[0] += 13 * 1.4 * X_remainder;
        V[4] += 9 * 6 * DY_remainder;
        // E += 9 * 1.3 * E_remainder;
        // E.split(E0,E_remainder);

        IVector v0 = V;
        IVector U(6);
        split(v0,U);
        v0[2] = V[2];
        
        IMatrix D(6,6);
        capd::C1HORect2Set S(V);
        capd::C1HORect2Set S0(v0);


        IVector Y = pm_z(S,D);
        IVector Y0 = pm_z(S0);

        IVector p0{Y0[1],Y0[3]};

        D = pm_z.computeDP(Y,D);
        IMatrix M{{D[1][0],D[1][4]},
                  {D[3][0],D[3][4]}};

        IVector N = - capd::matrixAlgorithms::gauss(M,p0);
        IVector delta{U[0],U[4]};
        // std::cout << N << std::endl;
        // std::cout << delta << std::endl;
        
        return subset(N,delta);
    }


    IVector cone_coeff(int tau_div, IMatrix T_total, IMatrix T_total_inv, capd::IPoincareMap& pm_y) {
        double tau_eps = 1e-7;
        double gamma_eps = 1e-9;

        double tau_delta = tau_eps / tau_div;

        Interval alpha_h = 0.001; Interval m_h = 1000 * 1000;
        Interval alpha_v = 0.001; Interval m_v = 1000 * 1000 + 1;

        IMatrix Q_h = -IMatrix::Identity(4); IMatrix Q_v = IMatrix::Identity(4);
        Q_h[0][0] = alpha_h;
        Q_v[1][1] = -alpha_v; Q_v[2][2] = -alpha_v; Q_v[3][3] = -alpha_v; 

        IMatrix D_total(5,5);

        Interval gamma = gamma_eps * Interval(-1,1);
        Interval tau_i(0, tau_delta);

        for(int i = 0; i < tau_div; i++) {
            
            IVector u{tau_i,0, 0, gamma, gamma, gamma};
            capd::C1HORect2Set S(capd::C1Rect2Set::C0BaseSet(V, T_total, u), capd::C1Rect2Set::C1BaseSet(T_total));
            
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
        

        IMatrix H = matrix_erase_cord(D_total,1);
        IMatrix H_T = H;
        H_T.transpose();
        IMatrix V_h = H_T * Q_h * H - m_h * Q_h;
        IMatrix V_v = H_T * Q_v * H - m_v * Q_v;
        
        for(int i = 0; i < 4; i++) {
            Interval R_v = 0; Interval R_h = 0;
            Interval a_v = V_v[i][i]; Interval a_h = V_h[i][i];
            for(int j = 0; j < 4 && j != i; j++) {
                R_v += capd::abs(V_v[i][j]);
                R_h += capd::abs(V_h[i][j]);
            }
            if( !(a_v - R_v > 0 && a_h - R_h > 0) ) {
                throw std::runtime_error("Some quadratic form is not positive definite");
            }
        }

        Interval A = capd::abs(D_total[0][0]).left();
        Interval B = sum_norm(IVector{H[0][1],H[0][2],H[0][3]}).right();
        Interval C = max_norm(IVector{H[1][0],H[2][0],H[3][0]}).right();

        Interval K = max_norm(matrix_erase_cord(H,0)).right();

        Interval delta = power(-A + K,2) - 4 * B * C;
        
        Interval alpha1 = (A - K + sqrt(delta)) / (2 * B);
        Interval alpha2 = (A - K - sqrt(delta)) / (2 * B);

        return IVector{alpha_h,alpha_v};
    }
    

    void run(unsigned threadID) {
        // prove_Z_bound();

        try {
            capd::IPoincareMap& pm_x = CR3BPPool::pool_X->getSolver(threadID).pm;
            capd::IPoincareMap& pm_y = CR3BPPool::pool_Y->getSolver(threadID).pm;
            capd::IPoincareMap& pm_z = CR3BPPool::pool_Z->getSolver(threadID).pm;

            if(E.left() == E.right()) {
                V[2] += 6e-16 * Interval(-1,1);
                Z = V[2];
                Z.split(Z0,Z_remainder);
            } 

            if(!prove_Z_bound()) throw std::runtime_error("Proper Z bound is not proven");
            
            if(E.left() == E.right()) {
                V[0] += 1e-15 * Interval(-1,1);
                V[4] += 1e-15 * Interval(-1,1);

                X = V[0]; DY = V[4];
                X.split(X0,X_remainder);
                DY.split(DY0,DY_remainder);
            }
            
            if(!prove_fixed_point(pm_z)) throw std::runtime_error("Proper enclosure of fixed points is not proven");
            

            IMatrix L_temp = capd::matrixAlgorithms::krawczykInverse(matrix_erase_cord(L,1));
            L_inv = matrix_add_cord(L_temp,1);
            

            IMatrix T = get_energy_change_of_basis();

            IMatrix T_inv = capd::matrixAlgorithms::krawczykInverse(T);
            IMatrix T_total = T * L;
            IMatrix T_total_inv = L_inv * T_inv;

            
            Interval alpha = cone_coeff(5,T_total, T_total_inv, pm_y)[0];
            // pm_y.getSolver().setOrder(15);
            alpha = 0.02;
            
            Interval gamma = alpha * tau * Interval(-1,1);
            
            IVector w{tau,0,0,gamma,gamma,gamma};
            // std::cout << w << std::endl;
            
            
            // w += T_total_inv * V_remainder;

            Interval returnTime;
            IVector zero_vector{0,0,0,0,0,0};

            IVector V0 = V;
            IVector V_remainder(6);
            split(V0,V_remainder);

            IVector w0 = w;
            IVector w_remainder(6);
            split(w0,w_remainder);

            // IMatrix T0 = T_total;
            // IMatrix T_remainder(6,6);
            // split(T0,T_remainder);

            std::cout << T_total * w_remainder << std::endl;
            std::cout << V_remainder << std::endl;

            capd::C0HOTripletonSet S(V0 + T_total * w0,T_total,w_remainder,V_remainder);
            S.setC0Factor(100);
            pm_x(S);
            
            IVector res = pm_y(S,zero_vector,C,returnTime);
            // IVector res = pm_y(S);
            
            result = IVector{res[3],res[5]};
            IVector z{0,0};
            std::cout << result << std::endl;
            if(subset(z,result)) throw std::runtime_error("Solution contains zero");
            
            

            status = COMPLETED;
        } catch(...) {
            status = FAILED;
        }
    }

    
};

#endif