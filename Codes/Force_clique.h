#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double length_cl(particle p1, particle p2){
    double l;
    l = sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
    return l;
}

pair<double, double> SpringForces_clique(polymer &p, int &i, double &ks, double &l0, int &N){
    double Fx, Fy, l;
    Fx = 0, Fy = 0;


    for(int j = 0; j<N; j++){
        if(i!=j){
            l = length_cl(p.chain[i],p.chain[j]);
            Fx += ks*(l-l0)*(p.chain[j].x - p.chain[i].x)/l;
            Fy += ks*(l-l0)*(p.chain[j].y - p.chain[i].y)/l;
        }

        //return {Fx, Fy};
    }
    return {Fx, Fy};
}