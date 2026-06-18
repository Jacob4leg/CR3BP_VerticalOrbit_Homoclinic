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

    double muSJ = 0.00095388114032796904;
    double e0 = 3.0268216642156945825;
    double tau_0 = 9.19967402684294e-08;

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
    capd::LDMap LDE_map;
    // InitConditionsGenerator& gen;
    Status status = NEW;

    ComputeImageTask(Interval tau, Interval E, IVector V, IMatrix C, IMatrix L, InitConditionsGenerator& gen) : 
            tau(tau), E(E), V(V),C(C),L(L) {
        
        // *gen_ptr = &gen;
        result = IVector(2);
        X = V[0]; Z = V[2]; DY = V[4];
  
        dy = capd::IMap(dy_map,6,1,2);
        dy.setParameter(0,muSJ);
        dy.setParameter(1,muSJ-1.);

        E_map = capd::IMap(energy,6,1,2);
        E_map.setParameter(0,muSJ);
        E_map.setParameter(1,muSJ-1.);

        LDE_map = capd::LDMap(energy,6,1,2);
        LDE_map.setParameter(0,muSJ);
        LDE_map.setParameter(1,muSJ-1.);

        // E = E.mid();

        
        E.split(E0,E_remainder);
        X.split(X0,X_remainder); Z.split(Z0,Z_remainder); DY.split(DY0,DY_remainder);

        // IMatrix L1 = gen.get_change_of_basis(E.leftBound() - gen.E0);
        // IMatrix L2 = gen.get_change_of_basis(E.rightBound() - gen.E0);
        // this->L = intervalHull(L1,L2);
        
    }

    IMatrix get_energy_change_of_basis(IVector U) {
        IVector dE = E_map.derivative(U)[0];
        IVector z_grad = -IVector{dE[0],dE[1],-1,dE[3],dE[4],dE[5]} / dE[2];
        IMatrix T = IMatrix::Identity(6);
        T.row(2) = z_grad;
        return T;
    }

    long double find_approximate_z(LDVector u0, long double e) {
        
        long double z_delta;
        do{
            long double f = LDE_map(u0)[0] - e;
            LDVector E_grad = LDE_map.derivative(u0)[0];
            z_delta = f / E_grad[2];
            u0[2] -= z_delta;
        }while (abs(z_delta) > 1e-15);

        return u0[2];
    }

    bool prove_Z_bound(IVector U, Interval e) {
        IVector U0 = U;
        IVector U_remainder(6);
        split(U0,U_remainder);
        U0 = U; U0[2] = U[2].mid();

        LDVector u0(6);
        for(int i = 0; i < 6; i++) u0[i] = U0[i].mid().leftBound();
        U0[2] = find_approximate_z(u0,e.mid().leftBound());

        IVector E_grad = E_map.derivative(U)[0];
        Interval p0 = E_map(U0)[0] - e;
        Interval N = - p0 / E_grad[2];
        
        // std::cout << N << std::endl;
        // std::cout << U_remainder[2] << std::endl;
        return subset(N,U_remainder[2]);
    }

    std::pair<Interval,IMatrix> get_proper_bound_for_z(IVector U, Interval e) {
        
        if(!prove_Z_bound(U,e)) throw std::runtime_error("Implicit function z(...) is not proven");
        
        IMatrix T = get_energy_change_of_basis(U);
        
        IVector U0 = capd::vectalg::midVector(U);

        Interval Z = 1e-14 * Interval(-1,1);

        // U0[2] = z;
        U0[2] += Z;

        if(!prove_Z_bound(U0,e.mid())) throw std::runtime_error("Proper Z bound is not proven");
        
        return {U0[2], T};
    }


    bool prove_fixed_point(IVector U, Interval e, capd::IPoincareMap& pm_z) {
        IVector u0 = U;
        IVector U_remainder(6);
        split(u0,U_remainder);
        
        IMatrix D(6,6);
        capd::C1HORect2Set S(U);
        capd::C1HORect2Set S0(u0);

        Interval e_val = E_map(u0)[0] - e;
        IVector e_grad = E_map.derivative(U)[0];

        IVector Y = pm_z(S,D);
        IVector Y0 = pm_z(S0);

        IVector p0{Y0[1],Y0[3],e_val};

        D = pm_z.computeDP(Y,D);
        IMatrix M{{D[1][0],D[1][2],D[1][4]},
                  {D[3][0],D[3][2],D[3][4]},
                  {e_grad[0],e_grad[2],e_grad[4]}};

        IVector N = - capd::matrixAlgorithms::gauss(M,p0);
        IVector delta{U_remainder[0],U_remainder[2],U_remainder[4]};
        
        return subset(N,delta);
    }


    IVector cone_coeff(IVector U, int tau_div, IMatrix L, IMatrix L_inv, capd::IPoincareMap& pm_y) {
        IVector u0 = U;
        IVector U_remainder(6);
        split(u0,U_remainder);
        
        double tau_eps = 1e-7;
        double gamma_eps = 1e-10;

        double tau_delta = tau_eps / tau_div;

        Interval alpha_h = 0.001; Interval m_h = 1000 * 1000;
        Interval alpha_v = 0.001; Interval m_v = 1000 * 1000 + 1;

        IMatrix Q_h = -IMatrix::Identity(4); IMatrix Q_v = IMatrix::Identity(4);
        Q_h[0][0] = alpha_h;
        Q_v[1][1] = -alpha_v; Q_v[2][2] = -alpha_v; Q_v[3][3] = -alpha_v; 

        IMatrix D_total(5,5);

        Interval gamma = tau_eps * sqrt(alpha_h) * Interval(-1,1);
        Interval tau_i(0, tau_delta);
        
        for(int i = 0; i < tau_div; i++) {
            
            IVector w{tau_i,0, 0, gamma, gamma, gamma};
            
            IVector z_test{0,0,1.9e-8 * Interval(-1,1),0,0,0};

            IVector W = U + L * w + z_test;
            auto [z, T] = get_proper_bound_for_z(W,E);
            
            IMatrix T_inv = capd::matrixAlgorithms::krawczykInverse(T);
            IMatrix T_total = T * L;
            IMatrix T_total_inv = L_inv * T_inv;

            w[2] += E_remainder;
            IVector w0 = w;
            IVector w_remainder(6);
            split(w0,w_remainder);
      
            u0[2] = z + U_remainder[2];
            U_remainder[2] = 0;

            capd::C1HORect2Set S(capd::C1Rect2Set::C0BaseSet(u0 + T_total * w0, T_total, w_remainder, T * U_remainder), capd::C1Rect2Set::C1BaseSet(T_total));
            
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

        // Interval A = capd::abs(D_total[0][0]).left();
        // Interval B = sum_norm(IVector{H[0][1],H[0][2],H[0][3]}).right();
        // Interval C = max_norm(IVector{H[1][0],H[2][0],H[3][0]}).right();

        // Interval K = max_norm(matrix_erase_cord(H,0)).right();

        // Interval delta = power(-A + K,2) - 4 * B * C;
        
        // Interval alpha1 = (A - K + sqrt(delta)) / (2 * B);
        // Interval alpha2 = (A - K - sqrt(delta)) / (2 * B);

        return IVector{alpha_h,alpha_v};
    }
    

    void run(unsigned threadID) {

        try {

            capd::IPoincareMap& pm_x = CR3BPPool::pool_X->getSolver(threadID).pm;
            capd::IPoincareMap& pm_y = CR3BPPool::pool_Y->getSolver(threadID).pm;
            capd::IPoincareMap& pm_z = CR3BPPool::pool_Z->getSolver(threadID).pm;

            if(E.left() == E.right()) {
                V[2] += 1e-14 * Interval(-1,1);
                Z = V[2];
                Z.split(Z0,Z_remainder);

                V[0] += 1e-14 * Interval(-1,1);
                V[4] += 1e-14 * Interval(-1,1);

                X = V[0]; DY = V[4];
                X.split(X0,X_remainder);
                DY.split(DY0,DY_remainder);
            } 
            V[2] += 0.3 * Z_remainder;
            V[0] += 0.3 * X_remainder;
            V[4] += 0.3 * DY_remainder;
            
            IVector V0 = V;
            IVector V_remainder(6);
            split(V0,V_remainder);

            
            if(!prove_fixed_point(V,E,pm_z)) throw std::runtime_error("Proper enclosure of fixed points is not proven");
            

            IMatrix L_temp = capd::matrixAlgorithms::krawczykInverse(matrix_erase_cord(L,1));
            L_inv = matrix_add_cord(L_temp,1);

            
            Interval alpha = cone_coeff(V,5,L, L_inv, pm_y)[0];
            
            Interval gamma = sqrt(alpha) * tau * Interval(-1,1);
            // Interval gamma = alpha * tau * Interval(-1,1);
            
            IVector w{tau,0,0,gamma,gamma,gamma};
            

            IVector z_test{0,0,5e-9 * Interval(-1,1),0,0,0};

            IVector W = V + L * w + z_test;
            IVector W0 = W;
            IVector W_remainder(6);
            split(W0,W_remainder);

            // std::cout << "&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&" << std::endl;
            
            auto [z,T] = get_proper_bound_for_z(W,E);
            // IVector z_vector{0,0,z,0,0,0};
            
            IMatrix T_total = T * L;

            w[2] += E_remainder;

            Interval returnTime;
            IVector zero_vector{0,0,0,0,0,0};

            V0 = V;
            split(V0,V_remainder);
            V0[2] = z;
            // std::cout << V[2] - intervalHull(V[2],z) << std::endl;
            // V0 += V_remainder;
            
            IVector w0 = w;
            IVector w_remainder(6);
            split(w0,w_remainder);

            // IMatrix T0 = T_total;
            // IMatrix T_remainder(6,6);
            // split(T0,T_remainder);
            V_remainder = 0;

            capd::C0HOTripletonSet S(V0 + T_total * w0,T_total,w_remainder,T * V_remainder);
            // capd::C0HOTripletonSet S(V0 + T_total * w0,T_total,w_remainder);
            
            // S.setC0Factor(100);
            pm_x(S);
            
            IVector res = pm_y(S,zero_vector,C,returnTime);
            IMatrix B = {{1e10,0},{0,1e8}};
            
            Interval unit_interval(0,1);
            
            IVector res1 = B * IVector{tau - tau_0, E - e0};
            IVector res2 = IVector{res[3],res[5]};
            result = capd::vectalg::intervalHull(res1,res2);
            IVector zero{0,0};
            std::cout << result << std::endl;
            if(subset(zero,result)) throw std::runtime_error("Solution contains zero");
            
            

            status = COMPLETED;
        } catch(const std::exception& exc) {
            std::cout << exc.what() << std::endl;
            status = FAILED;
        }
    }

    
};

#endif