#include<iostream>
#include "capd/capdlib.h"
#include "cr3bp.h"
#include "linalg_helper.h"
#include "InitConditionsGenerator.h"
#include "ComputeImageTask.h"

using namespace std;
using namespace capd;

void init_cr3bp_pools(capd::interval mu) {
    CR3BPPool::init();

    CR3BPPool::pool_X->setParameter(0,mu);
    CR3BPPool::pool_X->setParameter(1,mu-1.);
    CR3BPPool::pool_X->setOrder(20);

    CR3BPPool::pool_Y->setParameter(0,mu);
    CR3BPPool::pool_Y->setParameter(1,mu-1.);
    CR3BPPool::pool_Y->setOrder(20);
    
    CR3BPPool::pool_Z->setParameter(0,mu);
    CR3BPPool::pool_Z->setParameter(1,mu-1.);
    CR3BPPool::pool_Z->setOrder(20);
}

void test1() {
    InitConditionsGenerator generator(1e-10,1e-7,300,30000);
    long double E_delta = 1e-10;
    Interval E = E_delta * Interval(-1,1);

    LDVector v0 = generator.get_fixed_point(0.);
    LDVector v1 = generator.get_fixed_point(-E_delta);
    LDVector v2 = generator.get_fixed_point(E_delta);
    IVector V = intervalHull(IVector(v1),IVector(v2));
    
    Interval X = Interval(capd::min(v1[0],v2[0]),capd::max(v1[0],v2[0]));
    Interval Z = Interval(capd::min(v1[2],v2[2]),capd::max(v1[2],v2[2]));
    Interval DY = Interval(capd::min(v1[4],v2[4]),capd::max(v1[4],v2[4]));

    LDMatrix L1 = generator.get_change_of_basis(-E_delta);
    LDMatrix L2 = generator.get_change_of_basis(E_delta);
    IMatrix L = intervalHull(IMatrix(L1),IMatrix(L2));
    IMatrix L_temp = matrixAlgorithms::gaussInverseMatrix(matrix_erase_cord(L,1));
    IMatrix L_inv = matrix_add_cord(L_temp,1);

    long double mu = generator.vf.muSJ;

    IMap E_map(energy,6,1,2);
    E_map.setParameter(0,mu);
    E_map.setParameter(1,mu-1.);

    IMatrix T = IMatrix::Identity(6);
    
    IVector dE = E_map.derivative(IVector(V))[0];
    // cout << E_map.derivative(IVector(V)) << endl;
    IVector z_grad = -IVector{dE[0],dE[1],-1,dE[3],dE[4],dE[5]} / dE[2];
    
    T.row(2) = z_grad;
    // cout << L << endl;
    // cout << T << endl;

    LDMap E_map_nonrig(energy,6,1,2);
    E_map.setParameter(0,mu);
    E_map.setParameter(1,mu-1.);
    LDMatrix L_nonrig = generator.get_change_of_basis(0.);
    LDMatrix T_nonrig = LDMatrix::Identity(6);
    LDVector dE_nonrig = generator.vf.E.derivative(v0)[0];
    LDVector z_grad_nonrig = -LDVector{dE_nonrig[0],dE_nonrig[1],-1,dE_nonrig[3],dE_nonrig[4],dE_nonrig[5]} / dE_nonrig[2];
    T_nonrig.row(2) = z_grad_nonrig;

    LDMatrix D = generator.derivative(0.);
    // cout << T_nonrig << endl;
    // cout << T << endl;
    

    // cout << matrix_add_cord(matrixAlgorithms::gaussInverseMatrix(matrix_erase_cord(L_nonrig,1)),1) * matrixAlgorithms::gaussInverseMatrix(T_nonrig) * D * T_nonrig * L_nonrig << endl;

    // return;

    IMatrix T_total = T * L;
    IMatrix T_total_inv = L_inv * matrixAlgorithms::gaussInverseMatrix(T);
    // D = generator.derivative(0.);
    // cout << T_total_inv * IMatrix(D) * T_total << endl;
    // return;
    // cout << T_total << endl;
    // cout << T_nonrig * L_nonrig << endl;

    
    IMap vf(cr3bpVectorField,6,6,2);
    vf.setParameter(0,mu);
    vf.setParameter(1,mu-1.);
    IOdeSolver solver(vf,20);
    ICoordinateSection section(6,1);
    IPoincareMap pm_y(solver,section,poincare::PlusMinus);

    // return;
    cout << v0 << endl;
    ComputeImageTask task(Interval(0.),E,IVector(v0),L);
    cout << task.cone_coeff(10,T_total,T_total_inv,pm_y) << endl;
    // task.run(1);
    
}

void test() {
    CR3BP<long double> vf;

    IMap f(cr3bpVectorField,6,6,2);
    f.setParameter(0,vf.muSJ);
    f.setParameter(1,vf.muSJ-1.);
    IOdeSolver solver(f,20);
    ICoordinateSection section(6,2);
    IPoincareMap pm_z(solver,section);

    LDMap E_map_nonrig(energy,6,1,2);
    E_map_nonrig.setParameter(0,vf.muSJ);
    E_map_nonrig.setParameter(1,vf.muSJ-1.);

    IMap E_map(energy,6,1,2);
    E_map.setParameter(0,vf.muSJ);
    E_map.setParameter(1,vf.muSJ-1.);

    InitConditionsGenerator generator(1e-10,1e-7,300,30000);
    long double E_delta = 1e-12;
    LDVector v0 = generator.get_fixed_point(0.);
    Interval x0 = v0[0];
    Interval z0 = v0[2];
    Interval dy0 = v0[4];
    Interval E0 = E_map(IVector(v0))[0].mid();

    LDVector E_grad_nonrig = E_map_nonrig.derivative(v0)[0];
    
    IMatrix A = IMatrix({{E_grad_nonrig[0],-1,E_grad_nonrig[4]}}) / Interval(E_grad_nonrig[2]);
    

    LDVector v1 = generator.get_fixed_point(-E_delta);
    LDVector v2 = generator.get_fixed_point(E_delta);
    cout << v1 - v2 << endl;

    long double eps_correction = 0.;
    Interval correction = eps_correction * Interval(-1,1);
    Interval E = E0 + E_delta * Interval(-1,1);
    Interval X = Interval(capd::min(v1[0],v2[0]),capd::max(v1[0],v2[0]));
    Interval Z = Interval(capd::min(v1[2],v2[2]),capd::max(v1[2],v2[2])) + correction;
    Interval DY = Interval(capd::min(v1[4],v2[4]),capd::max(v1[4],v2[4]));
    
    IVector vec{X,0,Z,0,DY,0};
    IVector remainder(6);
    split(vec,remainder);
    
    IMatrix E_grad = E_map.derivative(IVector{X,0,Z,0,DY,0});
    IMatrix E_grad_point = E_map.derivative(IVector{X,0,z0,0,DY,0});

    Interval Dz_f = E_grad_point[0][2];
    IMatrix Dx_f({{E_grad_point[0][0],-1,E_grad_point[0][4]}});

    IMatrix Dx_g = Dx_f - vectalg::intersection(Dz_f * A, Dz_f / (Interval(E_grad_nonrig[2])) * IMatrix({{E_grad_nonrig[0],-1,E_grad_nonrig[4]}}) );
    cout << Dx_g << endl;

    Interval g_val = (E_map_nonrig(v0)[0] - E0) + (Dx_g * IVector{remainder[0],Interval(-1e-12,1e-12),remainder[4]})[0];
    cout << g_val << endl;
    
    Interval N = - g_val / E_grad[0][2];
    cout << N << endl;

    Interval w = remainder[2] + (A * IVector{remainder[0],Interval(-1e-12,1e-12),remainder[4]})[0];
    cout << w << endl;

    return;


    // Interval N = - (E_map(IVector{X,0,z0,0,DY,0})[0] - E) / E_grad[2];
    // cout << N << endl;
    // cout << remainder[2] << endl;
    // cout << subset(N,remainder[2]) << endl;



    return;
    
    // Interval second_correction = 1e-5 * Interval(-1,1);
    // X += second_correction;
    // DY += second_correction;

    // IVector u0(v0);
    // IVector U0(v0);
    // U0[2] = Z;
    // IVector U = U0;
    // U[0] = X; U[4] = DY;
    
    // IMatrix D(6,6);
    // C1HORect2Set S(U);
    // IVector V = pm_z(S,D);
    // D = pm_z.computeDP(V,D);
    // IMatrix M({{D[1][0],D[1][4]},{D[3][0],D[3][4]}});
    // C0HORect2Set S0(U);
    // V = pm_z(S);
    // IVector V0{V[1],V[3]};
    // IVector N_vec = - matrixAlgorithms::gauss(M,V0);
    // cout << N_vec << endl;
    // cout << IVector{remainder[0],remainder[2]} << endl;
    // cout << subset(N_vec, IVector{remainder[0],remainder[2]}) << endl;
    
}

int main() {
    cout.precision(20);


    std::setprecision(std::numeric_limits<long double>::max_digits10);
    // test1();
    // return 0;

    CR3BP<long double> vf;

    InitConditionsGenerator generator(1e-10,1e-7,300,1000);
    
    IVector v_E(generator.v0);
    IMatrix L(generator.get_change_of_basis(0.));
    // generator.get_isomorphism();
    generator.save_init_data_to_file();

    return 0;
    // std::cout << v_E << std::endl;

    capd::threading::ThreadPool pool(threadsNo);
    interval mu = 0.00095388114032796904;
    init_cr3bp_pools(mu);
    

    ComputeImageTask task(Interval(0),Interval(generator.E0), v_E, L);
    pool.process(&task);
    task.join();
    

    // generator.test();
    return 0;
    
    

    LDVector w0{0.9468923401720671061132517L,-4.072102120831082499146823e-24L,0.05316795353478707980175375L,
        -5.553112274845604899656097e-08L,-0.01115319054270743243833971L,-9.199674025000000205149813e-08L};

    long double z = 0.05316795030019990465013731L;
    long double dy = -0.011153205408683909493953946L;
    LDVector v0 = vf.findVerticalLyapunovOrbit(LDVector({0.94690970780356629,0,z,0,dy,0}));
    
    cout << v0 << endl;
    cout << vf.E(v0) << endl;

    LDVector v_test = v0;
    v_test[4] = vf.E(v0)[0];

    cout << vf.dy(v_test) << endl;

    return 0;

    // LDMatrix D1(6,6), D2(6,6), D(6,6);
    // LDVector v1 = vf.pm_y(v0,D1);
    // LDVector v2 = vf.pm_y(v1,D2);

    // D1 = vf.pm_y.computeDP(v1,D1);
    // D2 = vf.pm_y.computeDP(v2,D2);
    // D = D2 * D1;
    // LDVector rV(6), iV(6);
    // LDMatrix rVec(6,6), iVec(6,6);
    // computeEigenvaluesAndEigenvectors(D,rV,iV,rVec,iVec);
    // cout << rV << endl;
    // cout << iV << endl;
    // cout << rVec[0] << endl;
    

    // LDVector rV1(5), iV1(5);
    // LDMatrix rVec1(5,5), iVec1(5,5);
    // D = matrix_erase_cord(D,4);
    // computeEigenvaluesAndEigenvectors(D,rV1,iV1,rVec1,iVec1);
    // cout << rV1 << endl;
    // cout << iV1 << endl;
    // cout << rVec1[0] << endl;
    
        // IVector W(w0);

    // IVector data[] = {W,W,W,W,W,W,W,W,W,W,W};
    // vector<ComputeImageTask*> tasks;
    // for(auto x: data) tasks.push_back(new ComputeImageTask(x));
    // for(auto* task: tasks) pool.process(task);
    // for(auto* task : tasks) task->join();

    // for(auto* task: tasks) {
    //     cout << task->status << ": " << task->result << endl;
    //     delete task;
    // }
}