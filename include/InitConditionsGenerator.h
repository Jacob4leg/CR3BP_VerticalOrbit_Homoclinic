#include<iostream>
#include<fstream>
#include <iomanip>
#include "capd/capdlib.h"
#include "cr3bp.h"
#include "linalg_helper.h"

#ifndef _INIT_CONDITIONS_GENERATOR_H_
#define _INIT_CONDITIONS_GENERATOR_H_


struct InitConditionsGenerator
{
    typedef capd::LDVector LDVector;
    typedef capd::LDMatrix LDMatrix;

    LDVector w0{0.9468923401720671061132517L,-4.072102120831082499146823e-24L,0.05316795353478707980175375L,
        -5.553112274845604899656097e-08L,-0.01115319054270743243833971L,-9.199674025000000205149813e-08L};
    LDVector v0{0.94689235111997293925295213,0,0.053167950300199904636584778,0,-0.011153205408683909493953946,0};

    
    long double tau_0; long double tau_radius;
    long double E0; long double E_radius;
    int tau_div; int E_div;
    long double tau_delta; long double E_delta;

    // LDVector v_eig;

    // LDVector V_eig;
    CR3BP<long double> vf;

    InitConditionsGenerator(long double tau_radius, long double E_radius, int tau_div, int E_div) : 
        tau_radius(tau_radius), E_radius(E_radius), tau_div(tau_div), E_div(E_div) {
        

        tau_delta = tau_radius * 2 / tau_div;
        E_delta = E_radius * 2 / E_div;
        E0 = vf.E(w0)[0];

        // LDMatrix D1(6,6), D2(6,6), D(6,6);
        // LDVector v1 = vf.pm_y(v0,D1);
        // LDVector v2 = vf.pm_y(v1,D2);

        // D1 = vf.pm_y.computeDP(v1,D1);
        // D2 = vf.pm_y.computeDP(v2,D2);
        // D = D2 * D1;
        
        LDMatrix D = derivative(0.);

        LDVector rV(6), iV(6);
        LDMatrix rVec(6,6), iVec(6,6);
        capd::computeEigenvaluesAndEigenvectors(D,rV,iV,rVec,iVec);
        LDVector v_eig = rVec.column(0);

        tau_0 = (w0[0] - v0[0]) / v_eig[0];
    }

    LDVector get_fixed_point(long double E) {
        return vf.findVerticalLyapunovOrbit(v0,E + E0);
    }

    LDVector get_initial_point(long double tau, long double E) {
        LDVector v_E = get_fixed_point(E);
        LDMatrix D = derivative(E);
        LDVector rV(6),iV(6);
        LDMatrix rVec(6,6), iVec(6,6);
        capd::computeEigenvaluesAndEigenvectors(D,rV,iV,rVec,iVec);

        LDVector v_res = v_E + (tau_0 + tau) * rVec.column(0);
        return v_res;
    }

    LDMatrix get_isomorphism() {
        LDVector u1 = get_initial_point(tau_radius,0);
        LDVector u2 = get_initial_point(0,E_radius);

        LDVector Pu1_x = vf.pm_x(u1); LDVector Pu1_y = vf.pm_y(Pu1_x);
        LDVector Pu2_x = vf.pm_x(u2); LDVector Pu2_y = vf.pm_y(Pu2_x);

        LDMatrix res = LDMatrix::Identity(6);
        res[3][3] = Pu1_y[3]; res[3][5] = Pu2_y[3];
        res[5][3] = Pu1_y[5]; res[5][5] = Pu2_y[5];

        LDVector u = vf.pm_y(vf.pm_x(w0));

        res.column(1) = vf.vf(u);

        return res;
    }

    LDMatrix derivative(long double E) {
        LDVector v_E = get_fixed_point(E);
        LDMatrix D = LDMatrix::Identity(6), T(6,6);
        LDVector u = v_E;
        for(int i = 0; i < 2; i++) {
            u = vf.pm_y(u,T);
            D = vf.pm_y.computeDP(u,T) * D;
        }
        return D;
    }

    LDMatrix get_change_of_basis(long double E) {
        // E -= E0;
        
        LDVector v_E = get_fixed_point(E);
        // long double dy = v_E[4];
        // v_E[4] = E0 + E;
        LDMatrix T = LDMatrix::Identity(6);
        
        LDVector dE = vf.E.derivative(v_E)[0];
        LDVector z_grad = -LDVector{dE[0],dE[1],-1,dE[3],dE[4],dE[5]} / dE[2];

        T.row(2) = z_grad;

        LDMatrix D = derivative(E);
        D = capd::matrixAlgorithms::gaussInverseMatrix(T) * D * T;
        D = matrix_add_cord(matrix_erase_cord(D,2),2);
        D[2][2] = 1.;

        LDMatrix D_new = matrix_erase_cord(D,1);

        LDMatrix L = matrix_add_cord(matrix_add_cord(us_change_of_basis(matrix_erase_cord(D_new,1)),1),1);
        
        L[2][2] = 1.;

        
        // LDMatrix T_total = T * L;
        // LDMatrix T_total_inv = matrix_add_cord(capd::matrixAlgorithms::gaussInverseMatrix(matrix_erase_cord(L,1)),1) * capd::matrixAlgorithms::gaussInverseMatrix(T);
        
        // D = derivative(E);
        // std::cout << D << std::endl;
        // std::cout << T_total_inv * D * T_total << std::endl;

        return L;
    }

    

    std::string vectorToString(LDVector v) {
        std::ostringstream oss; 
        oss << std::scientific
            << std::setprecision(std::numeric_limits<long double>::max_digits10);
        for(int i = 0; i < v.dimension(); i++) {
            oss << v[i];
            if(i != v.dimension() - 1) oss << " ";
        }
        return oss.str();
    }

    /**
     * Saves initial data to rigorous solvers in format
     * tau_leftBound tau_rightBound E_leftBound E_rightBound x y z dx dy dz
     * where v_E = (x,y,z,dx,dy,dz) is an approximate fixed point on energy level (E_rightBound + E_leftBound) / 2
     */
    void save_init_data_to_file() {
        long double tau_i = -tau_radius + tau_0;
        long double E_i = -E_radius + E0;

        std::ofstream file("init_data.txt");
        file << std::scientific
             << std::setprecision(std::numeric_limits<long double>::max_digits10);

        long double fixed_point_error;
        long double total_fixed_point_error = 0.;
        long double energy_error;
        long double total_energy_error = 0.;
        
        capd::LDMaxNorm norm;

        LDVector v_E1 = get_fixed_point(-E_radius);
        fixed_point_error = norm(v_E1 - vf.pm_y(vf.pm_y(v_E1)));
        energy_error = abs(vf.E(v_E1)[0] + E_radius - E0);

        if(total_fixed_point_error < fixed_point_error) total_fixed_point_error = fixed_point_error;
        if(total_energy_error < energy_error) total_energy_error = energy_error;


        LDVector v_E2 = get_fixed_point(E_radius);
        fixed_point_error = norm(v_E2 - vf.pm_y(vf.pm_y(v_E2)));
        energy_error = abs(vf.E(v_E2)[0] - E_radius - E0);

        if(total_fixed_point_error < fixed_point_error) total_fixed_point_error = fixed_point_error;
        if(total_energy_error < energy_error) total_energy_error = energy_error;

        for(int i = 0; i < tau_div; i++) {
            file << tau_i << " " << -E_radius + E0 << " " << vectorToString(v_E1) << std::endl;
            tau_i += tau_delta;
        }
        
        file << std::endl;
        tau_i = -tau_radius + tau_0;

        for(int i = 0; i < tau_div; i++) {
            file << tau_i << " " <<  E_radius + E0 << " " << vectorToString(v_E2) << std::endl;
            tau_i += tau_delta;
        }

        file << std::endl;

        for(int i = 0; i < E_div; i++) {
            LDVector v_E = get_fixed_point(E_i - E0);
            fixed_point_error = norm(v_E - vf.pm_y(vf.pm_y(v_E)));
            energy_error = abs(vf.E(v_E)[0] - E_i);

            if(total_fixed_point_error < fixed_point_error) total_fixed_point_error = fixed_point_error;
            if(total_energy_error < energy_error) total_energy_error = energy_error;

            file << -tau_radius + tau_0 << " " << E_i << " " << vectorToString(v_E) << std::endl;
            E_i += E_delta;
        }

        file << std::endl;
        E_i = -E_radius + E0;

        for(int i = 0; i < E_div; i++) {
            LDVector v_E = get_fixed_point(E_i - E0);
            file <<  tau_radius + tau_0 << " " << E_i << " " << vectorToString(v_E) << std::endl;
            E_i += E_delta;
        }

        std::cout << "data saved to init_data.txt with error: " << std::endl;
        std::cout << "fixed_point_error: " << total_fixed_point_error << std::endl;
        std::cout << "energy_error: " << total_energy_error << std::endl;

        file.close();
    }

    void save_rectangle_to_file() {
        long double tau_i = -tau_radius;
        long double E_i = -E_radius;

        std::ofstream file("init_points.txt");
        file << std::scientific
             << std::setprecision(std::numeric_limits<long double>::max_digits10);
        
        for(int i = 0; i < tau_div; i++) {
            
            LDVector u1 = get_initial_point(tau_i, -E_radius);
            LDVector u2 = get_initial_point(tau_i, E_radius);

            file << tau_i << " " << -E_radius << " " << vectorToString(u1) << std::endl;
            file << tau_i << " " << E_radius << " " << vectorToString(u2) << std::endl;

            tau_i += tau_delta;
        }

        for(int i = 0; i < E_div; i++) {
            LDVector u1 = get_initial_point(-tau_radius, E_i);
            LDVector u2 = get_initial_point(tau_radius, E_i);

            file << -tau_radius << " " << E_i << " " << vectorToString(u1) << std::endl;
            file << tau_radius << " " << E_i << " " << vectorToString(u2) << std::endl;
        
            E_i += E_delta;
        }
        file.close();
    }

    void test() {

        std::ifstream file("init_data.txt");
        std::ofstream outfile("image_nonrig.txt");
        std::string line;

        auto B = get_isomorphism();
        auto C = capd::matrixAlgorithms::gaussInverseMatrix(B);

        for(int i = 0; i < 4; i++) {
            while(getline(file,line)) {
                if(line == "") break;
                long double* out = new long double[8];
                std::istringstream ss(line);
                for(int i = 0; i < 8; i++) {
                    if(!(ss >> out[i])) {
                        throw std::runtime_error("Problem with parsing data");
                    }
                }
                long double tau = out[0];
                long double E = out[1];
                LDVector v{out[2],out[3],out[4],out[5],out[6],out[7]};
                LDVector rV(6),iV(6);
                LDMatrix rVec(6,6), iVec(6,6);
                LDMatrix D = derivative(E - E0);
                capd::computeEigenvaluesAndEigenvectors(D,rV,iV,rVec,iVec);
                LDVector v_eig = rVec.column(0);
                
                LDVector v1 = v + (tau) * v_eig;
                LDVector u = vf.pm_y(vf.pm_x(v1));
                u = C * u;
                outfile << u[3] << " " << u[5] << std::endl;
                delete out;
            }
        }

        file.close();
        outfile.close();
    }

    // void test1() {
    //     long double tau_i = -tau_radius;
    //     long double E_i = -E_radius;

    //     std::ofstream file("output.txt");
    //     file << std::scientific
    //          << std::setprecision(std::numeric_limits<long double>::max_digits10);
        
    //     for(int i = 0; i < tau_div; i++) {
            
    //         LDVector u = get_initial_point(tau_i, -E_radius);
    //         LDVector res = vf.pm_y(vf.pm_x(u));
    //         file << res[3] << " " << res[5] << std::endl;
    //         u = get_initial_point(tau_i, E_radius);
    //         res = vf.pm_y(vf.pm_x(u));
    //         file << res[3] << " " << res[5] << std::endl;

    //         tau_i += tau_delta;
    //     }

    //     for(int i = 0; i < E_div; i++) {
    //         LDVector u = get_initial_point(-tau_radius, E_i);
    //         LDVector res = vf.pm_y(vf.pm_x(u));
    //         file << res[3] << " " << res[5] << std::endl;
    //         u = get_initial_point(tau_radius, E_i);
    //         res = vf.pm_y(vf.pm_x(u));
    //         file << res[3] << " " << res[5] << std::endl;

    //         E_i += E_delta;
    //     }
    //     file.close();
    // }
    

    // void test() {
    //     // save_rectangle_to_file();

    //     LDMatrix C = get_isomorphism();

    //     std::ifstream file("init_points.txt");
        
    //     if(!file) throw std::runtime_error("File does not exist");

    //     std::ofstream file1("output.txt");
    //     file1 << std::scientific
    //          << std::setprecision(std::numeric_limits<long double>::max_digits10);

    //     LDVector v(6);
    //     long double empty;
    //     while(file >> empty >> empty >> v[0] >> v[1] >> v[2] >> v[3] >> v[4] >> v[5]) {
    //         LDVector u_bufor = vf.pm_x(v);
    //         LDVector u = vf.pm_y(u_bufor);
    //         u = C * u;
    //         file1 << u[3] << " " << u[5] << std::endl;
    //     }
    //     file.close();
    //     file1.close();
    // }
};

// void test() {
    //     std::cout << E_radius << std::endl;

    //     LDVector u1 = get_fixed_point(-E_radius);
    //     LDVector u2 = get_fixed_point(E_radius);
        
    //     std::ofstream file("test.txt");
    //     file << std::scientific
    //          << std::setprecision(std::numeric_limits<long double>::max_digits10);

    //     long double E_i = -E_radius;
    //     for(int i = 0; i < E_div; i++) {
    //         LDVector u = get_fixed_point(E_i);
    //         file << u[2] << " " << u[4] << std::endl;
    //         E_i += E_delta;
    //     }
    //     file.close();

    //     return;

    //     LDMatrix D = derivative(0.);
    //     LDMatrix T = LDMatrix::Identity(6);
    //     LDVector dE = vf.E.derivative(v0)[0];
    //     LDVector E_grad{dE[0],dE[1],-1,dE[3],dE[4],dE[5]};
    //     LDVector z_grad = -E_grad / dE[2];

        
    //     // LDVector rV(6), iV(6);
    //     // capd::alglib::computeEigenvalues(D,rV,iV);
    //     // std::cout << rV << std::endl;
    //     // std::cout << iV << std::endl;
        
    //     T.row(2) = z_grad;
    //     // std::cout << T << std::endl;
    //     D = capd::matrixAlgorithms::gaussInverseMatrix(T) * D * T;
    //     D = matrix_add_cord(matrix_erase_cord(D,2),2);
    //     D[2][2] = 1.;
    //     // std::cout << D << std::endl;


    //     LDMatrix D_new = matrix_erase_cord(D,1);
    //     LDMatrix L = matrix_add_cord(matrix_add_cord(us_change_of_basis(matrix_erase_cord(D_new,1)),1),1);
    //     L[2][2] = 1.;
        
    //     LDMatrix T_total = T * L;
    //     LDMatrix T_total_inv = matrix_add_cord(capd::matrixAlgorithms::gaussInverseMatrix(matrix_erase_cord(L,1)),1) * capd::matrixAlgorithms::gaussInverseMatrix(T);
    //     D = derivative(0.);
    //     std::cout << T_total_inv * D * T_total << std::endl;

    //     LDMatrix D_temp(6,6);
    //     LDVector w1 = vf.pm_y(w0,D_temp);
    //     D = vf.pm_y.computeDP(w1,D_temp);
    //     LDVector w2 = vf.pm_y(w1,D_temp);
    //     D = vf.pm_y.computeDP(w2,D_temp) * D;

    //     std::cout << T_total_inv * D * T_total << std::endl;

        
    // }

#endif