#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double length_s(particle p1, particle p2){
    double l;
    l = sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
    return l;
}


pair<double, double> SpringForces_star(polymer &p, int &i, double &ks, double &l0, int &N){
    double Fx, Fy, l;
    Fx = 0, Fy = 0;

    if(i==0){
        for(int j = 1; j<N; j++){
            l = length_s(p.chain[i],p.chain[j]);
            Fx += ks*(l-l0)*(p.chain[j].x - p.chain[i].x)/l;
            Fy += ks*(l-l0)*(p.chain[j].y - p.chain[i].y)/l;
        }

        //return {Fx, Fy};
    }
    
    if(i!=0){
        l = length_s(p.chain[i], p.chain[0]);
        Fx = ks*(l-l0)*(p.chain[0].x - p.chain[i].x)/l;
        Fy = ks*(l-l0)*(p.chain[0].y - p.chain[i].y)/l;

        //return {Fx, Fy};
    }
    return {Fx, Fy};
}