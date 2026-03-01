#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

double length_c(particle p1, particle p2){
    double l;
    l = sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
    return l;
}


pair<double, double> SpringForces_chain(polymer &p, int &i, double &ks, double &l0, int &N){
    double Fx, Fy, l;
    Fx = 0, Fy = 0;

    if(i!=0 && i!=p.N-1){
        l = length_c(p.chain[i], p.chain[i-1]);
        Fx = ks*(l-l0)*(p.chain[i-1].x - p.chain[i].x)/l;
        Fy = ks*(l-l0)*(p.chain[i-1].y - p.chain[i].y)/l;
        

        l = length_c(p.chain[i], p.chain[i+1]);
        Fx += ks*(l-l0)*(p.chain[i+1].x - p.chain[i].x)/l;
        Fy += ks*(l-l0)*(p.chain[i+1].y - p.chain[i].y)/l;
    
        //return {Fx, Fy};
    }

    else if(i==0){
        l = length_c(p.chain[i], p.chain[i+1]);
        Fx = ks*(l-l0)*(p.chain[i+1].x - p.chain[i].x)/l;
        Fy = ks*(l-l0)*(p.chain[i+1].y - p.chain[i].y)/l;

        //return {Fx, Fy};
    }

    else if(i==N-1){
        l = length_c(p.chain[i-1], p.chain[i]);
        Fx = ks*(l-l0)*(p.chain[i-1].x - p.chain[i].x)/l;
        Fy = ks*(l-l0)*(p.chain[i-1].y - p.chain[i].y)/l;

        //return {Fx, Fy};
    }
    return {Fx, Fy};
}