#include<iostream>
#include<tuple>
#include<vector>
#include<fstream>
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

struct Record{
    long double tau;
    long double E;
    LDVector V;

    Record(string line) {
        long double* out = new long double[8];
        istringstream ss(line);
        for(int i = 0; i < 8; i++) {
            if(!(ss >> out[i])) {
                throw runtime_error("Problem with parsing data");
            }
        }
        tau = out[0];
        E = out[1];
        V = LDVector{out[2],out[3],out[4],out[5],out[6],out[7]};
        delete out;
    }
};

struct DataInstance{
    Interval tau;
    Interval E;
    IVector V;

    DataInstance(Interval tau, Interval E, IVector V) : tau(tau), E(E), V(V) {}
};

vector<DataInstance> set_data(string file_name, InitConditionsGenerator& gen) {
    ifstream file(file_name);
    if(!file) throw runtime_error("File does not exist");

    vector<DataInstance> res;

    // long double tau_left, tau_right;
    // long double E_left, E_right;
    // LDVector v_left(6), v_right(6);
    
    string line;
    

    for(int i = 0; i < 4; i++) {
        getline(file,line);
        if(line == "") break;
        Record left_record(line);
        
        while(getline(file,line)) {
            
            if(line == "") break;
            Record right_record(line);

            Interval tau(left_record.tau, right_record.tau);
            Interval E(left_record.E, right_record.E);
            IVector V = intervalHull(IVector(left_record.V),IVector(right_record.V));

            auto data = DataInstance(tau,E,V);
            res.push_back(data);

            left_record = right_record;
        }
        
    }
    return res;
}

int main() {
    cout.precision(20);

    std::setprecision(std::numeric_limits<long double>::max_digits10);
    
    capd::threading::ThreadPool pool(threadsNo);
    interval mu = 0.00095388114032796904;
    init_cr3bp_pools(mu);

    CR3BP<long double> vf;

    InitConditionsGenerator generator(1e-10,1e-10,300,3000);

    // generator.save_init_data_to_file();
    // return 0;
    auto data = set_data("init_data.txt", generator);


    vector<ComputeImageTask*> tasks;
    LDMatrix B = generator.get_isomorphism();

    fstream file("output.txt");
    file << std::scientific
             << std::setprecision(std::numeric_limits<long double>::max_digits10);

    for(int i = 1000; i < 1001; i++) {
        Interval tau = data[i].tau;
        Interval E = data[i].E;
        IVector V = data[i].V;

        tasks.push_back(new ComputeImageTask(tau,E,V,B,generator));
        // break;
    }
    for(auto* task: tasks) pool.process(task);
    for(auto* task: tasks) task->join();

    bool is_everything_ok = true;

    for(auto* task: tasks) {
        if(task->status != 1) is_everything_ok = false;

        Interval x = task->result[0];
        Interval y = task->result[1];

        file << x.leftBound() << " " << x.rightBound() << " " << y.leftBound() << " " << y.rightBound() << endl;
        // cout << x << " " << y << endl;
        delete task;
    }

    cout << is_everything_ok << endl;

    // file.close();

    return 0;
}

