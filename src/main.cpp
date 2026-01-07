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

int main() {
    cout.precision(26);
    std::setprecision(std::numeric_limits<long double>::max_digits10);
    CR3BP<long double> vf;

    InitConditionsGenerator generator(1e-10,1e-7,300,300);
    generator.get_change_of_basis(0);
    // generator.test();
    return 0;
    
    // capd::threading::ThreadPool pool(threadsNo);
    // interval mu = 0.00095388114032796904;
    // init_cr3bp_pools(mu);

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

    LDMatrix D1(6,6), D2(6,6), D(6,6);
    LDVector v1 = vf.pm_y(v0,D1);
    LDVector v2 = vf.pm_y(v1,D2);

    D1 = vf.pm_y.computeDP(v1,D1);
    D2 = vf.pm_y.computeDP(v2,D2);
    D = D2 * D1;
    LDVector rV(6), iV(6);
    LDMatrix rVec(6,6), iVec(6,6);
    computeEigenvaluesAndEigenvectors(D,rV,iV,rVec,iVec);
    cout << rV << endl;
    cout << iV << endl;
    cout << rVec[0] << endl;
    

    LDVector rV1(5), iV1(5);
    LDMatrix rVec1(5,5), iVec1(5,5);
    D = matrix_erase_cord(D,4);
    computeEigenvaluesAndEigenvectors(D,rV1,iV1,rVec1,iVec1);
    cout << rV1 << endl;
    cout << iV1 << endl;
    cout << rVec1[0] << endl;
    
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