#include <iostream>
#include <list>
#include "capd/capdlib.h"
#include "capd/threading/SolverPool.h"
#include "capd/threading/ThreadPool.h"


#ifndef _CR3BP_H_
#define _CR3BP_H_

const int threadsNo = std::thread::hardware_concurrency();

/**
 This routine implements vector field of the CR3BP.
The routine is written in the format required by the constructor
of the class DMap and IMap from the CAPD library.
*/
inline void cr3bpVectorField(capd::autodiff::Node t, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node params[], int noParams){
    typedef capd::autodiff::Node Node;

    Node xMu = in[0] + params[0]; // params[0] - relative mass of the first body
    Node xMj = in[0] + params[1]; // params[1] = params[0] - 1
    Node xMuSquare = xMu^2;       // square
    Node xMjSquare = xMj^2;
    Node ySquare = in[1]^2;
    Node zSquare = in[2]^2;
    Node yzSquare = ySquare + zSquare;
    
    Node factor1 = params[1]*((xMuSquare+yzSquare)^-1.5);
    Node factor2 = params[0]*((xMjSquare+yzSquare)^-1.5);
    Node factor = factor1 - factor2;

    out[0] = in[3];
    out[1] = in[4];
    out[2] = in[5];

    out[3] = in[0] + xMu*factor1 - xMj*factor2 + 2*in[4];
    out[4] = in[1]*(1+factor) - 2*in[3];
    out[5] = in[2]*factor;
}

/*
Implements energy constant as a map
*/
inline void energy(capd::autodiff::Node t, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node params[], int noParams){
    typedef capd::autodiff::Node Node;
    
    Node xMu = params[0]; // params[0] - relative mass of the first body
    Node xMj = params[1]; // params[1] = params[0] - 1

    Node x = in[0];
    Node y = in[1];
    Node z = in[2];

    Node r1 = sqrt( ((x + xMu)^2) + (y^2) + (z^2) );
    Node r2 = sqrt( ((x + xMj)^2) + (y^2) + (z^2) );

    Node omega = (x^2) + (y^2) + 2 * (-xMj/r1 + xMu/r2);
    
    Node dx = in[3];
    Node dy = in[4];
    Node dz = in[5];

    Node kinetic = (dx^2) + (dy^2) + (dz^2);

    out[0] = omega - kinetic;
}

inline void dy_map(capd::autodiff::Node t, capd::autodiff::Node in[], int dimIn, capd::autodiff::Node out[], int dimOut, capd::autodiff::Node params[], int noParams) {
    typedef capd::autodiff::Node Node;
    
    Node xMu = params[0]; // params[0] - relative mass of the first body
    Node xMj = params[1]; // params[1] = params[0] - 1

    Node x = in[0]; Node y = in[1]; Node z = in[2];
    Node dx = in[3]; Node E = in[4]; Node dz = in[5]; 

    Node r1 = sqrt( ((x + xMu)^2) + (y^2) + (z^2) );
    Node r2 = sqrt( ((x + xMj)^2) + (y^2) + (z^2) );

    Node omega = (x^2) + (y^2) + 2 * (-xMj/r1 + xMu/r2);

    out[0] = -sqrt( omega - (dx^2) - (dz^2) - E );
}


template<typename T>
struct CR3BP {
    long double muSJ = 0.00095388114032796904; // relative mass of Jupiter in Sun-Jupiter system

    typedef capd::vectalg::Matrix<T, 0U, 0U> TMatrix;
    typedef capd::vectalg::Vector<T,0U> TVector;
    typedef capd::map::Map<TMatrix> TMap;

    TMap vf;
    TMap E;
    TMap dy;
    capd::dynsys::BasicOdeSolver<TMap> solver;
    capd::poincare::CoordinateSection<TMatrix> section_y;
    capd::poincare::CoordinateSection<TMatrix> section_x;
    capd::poincare::CoordinateSection<TMatrix> section_z;
    capd::poincare::BasicPoincareMap<capd::dynsys::BasicOdeSolver<TMap>> pm_x;
    capd::poincare::BasicPoincareMap<capd::dynsys::BasicOdeSolver<TMap>> pm_z;
    capd::poincare::BasicPoincareMap<capd::dynsys::BasicOdeSolver<TMap>> pm_y;

    capd::vectalg::MaxNorm<TVector,TMatrix> norm;

    // #########################################################################
    /**
     The vertical Lyapunov orbits are double symmetric with respect to reflection z->-z and R-symmetry.
    In order to find them it is enough to fix an amplitude z\neq 0 and solve for
    (x,dy) -> (x,0,z,0,dy,0) -> P(x,0,z,0,dy,0) = (x1,y1,z1,dx1,dy1,dz1) = (x1,0,0,0,dy1,dz1)
    Setting Poincare section z=0 we obtain two equations 
    y1=0 and dx1=0
    with two variables (x,dy).
    This routine implements this search.
    @param [in] v - an approximate vertical Lyapunov orbit of the form (x,0,z,0,dy,0)
    @return - initial condition for vertical Lyapunov orbit (x,0,z,0,dy,0) for fixed z.
    */
    TVector findVerticalLyapunovOrbit(TVector v, T e){
        TVector u(3);
        do{
            TMatrix D(6,6);
            TVector y = pm_z(v,D);
            D = pm_z.computeDP(y,D);
            TVector E_grad = E.derivative(v)[0];

            // we want y=0 and dx=0
            u = TVector({y[1],y[3], E(v)[0] - e}); 
            
            // input variables are x,z (0 and 2)
            TMatrix M({{D[1][0],D[1][2],D[1][4]},
                       {D[3][0],D[3][2],D[3][4]},
                       {E_grad[0],E_grad[2],E_grad[4]}});
            
            u = capd::matrixAlgorithms::gauss(M,u);
            v[0] -= u[0];
            v[2] -= u[1];
            v[4] -= u[2];
        } while(norm(u) > 1e-15);
        
        // for(int i=0;i<15;++i){
        //     TMatrix D(6,6);
        //     TVector y = pm_z(v,D);
        //     D = pm_z.computeDP(y,D);
        //     TVector E_grad = E.derivative(v)[0];

        //     // we want y=0 and dx=0
        //     TVector u({y[1],y[3], E(v)[0] - e}); 
            
        //     // input variables are x,z (0 and 2)
        //     TMatrix M({{D[1][0],D[1][2],D[1][4]},
        //                {D[3][0],D[3][2],D[3][4]},
        //                {E_grad[0],E_grad[2],E_grad[4]}});
            
        //     u = capd::matrixAlgorithms::gauss(M,u);
        //     v[0] -= u[0];
        //     v[2] -= u[1];
        //     v[4] -= u[2];
        // }
        return v;
    }

    CR3BP() : 
        vf(cr3bpVectorField,6,6,2),
        E(energy,6,1,2),
        dy(dy_map,6,1,2),
        solver(vf,20),
        section_x(6,0),
        section_y(6,1),
        section_z(6,2),
        pm_x(solver,section_x),
        pm_y(solver,section_y,capd::poincare::PlusMinus),
        pm_z(solver,section_z) 
    {
        vf.setParameter(0,muSJ);
        vf.setParameter(1,muSJ-1.);

        E.setParameter(0,muSJ);
        E.setParameter(1,muSJ-1.);

        dy.setParameter(0,muSJ);
        dy.setParameter(1,muSJ-1.);
    }    
};

struct CR3BPPool{
    typedef capd::threading::PMap<capd::IPoincareMap,capd::ICoordinateSection> IPMap;
    typedef capd::threading::PMapFactory<capd::IPoincareMap,capd::ICoordinateSection> IPMapFactory;
    typedef capd::threading::SolverPool<IPMap> IPMapPool;

    static std::unique_ptr<IPMapPool> pool_X;
    static std::unique_ptr<IPMapPool> pool_Y;
    static std::unique_ptr<IPMapPool> pool_Z;

    static void init() {
        pool_X.reset(new IPMapPool(threadsNo, IPMapFactory(cr3bpVectorField,20,6,2,1,capd::ICoordinateSection(6,0),capd::poincare::Both)));
        pool_Y.reset(new IPMapPool(threadsNo, IPMapFactory(cr3bpVectorField,20,6,2,1,capd::ICoordinateSection(6,1),capd::poincare::PlusMinus)));
        pool_Z.reset(new IPMapPool(threadsNo, IPMapFactory(cr3bpVectorField,20,6,2,1,capd::ICoordinateSection(6,2),capd::poincare::Both)));
    }
};

std::unique_ptr<CR3BPPool::IPMapPool> CR3BPPool::pool_X;
std::unique_ptr<CR3BPPool::IPMapPool> CR3BPPool::pool_Y;
std::unique_ptr<CR3BPPool::IPMapPool> CR3BPPool::pool_Z;
// capd::IMap E_map(energy,6,1,2);


// struct ICR3BP {

//     capd::interval muSJ = 0.00095388114032796904;

//     capd::IMap vf;
//     capd::IMap E;
//     capd::IC2OdeSolver solver;
//     capd::ICoordinateSection section_y;
//     capd::ICoordinateSection section_x;
//     capd::ICoordinateSection section_z;
//     capd::IC2PoincareMap pm_x;
//     capd::IC2PoincareMap pm_z;
//     capd::IC2PoincareMap pm_y;

//     ICR3BP() : 
//         vf(cr3bpVectorField,6,6,2),
//         E(energy,6,1,2),
//         solver(vf,20),
//         section_x(6,0),
//         section_y(6,1),
//         section_z(6,2),
//         pm_x(solver,section_x),
//         pm_y(solver,section_y,capd::poincare::PlusMinus),
//         pm_z(solver,section_z) 
//     {
//         vf.setParameter(0,muSJ);
//         vf.setParameter(1,muSJ-1.);

//         E.setParameter(0,muSJ);
//         E.setParameter(1,muSJ-1.);

//         solver.setAbsoluteTolerance(1e-12);
//         solver.setRelativeTolerance(1e-12);
//     }    
// };

#endif