#include <iostream>
#include <cmath>
#include <vector>

using namespace std;

const double pi = 3.1415926;

class particle{
    public:
    double x,y, phi;
    int alpha;
    double or_x, or_y;
};

class polymer{
    public:
    unsigned int N;
    vector<particle> chain;

    polymer(int num){
        N = num;
        chain.resize(num);
        for(int i=0; i<N; i++){
            chain[i].alpha = 1;
        }
    }  
};

double CentreOfFriction(polymer &p, vector<double> &fric){
    double sum_x, sum_fric; sum_x = 0.0; sum_fric = 0.0;
    for(int i=0; i<p.N; i++){
        sum_x+=p.chain[i].x*fric[i];
        sum_fric += fric[i];
    }
    return sum_x/sum_fric;
}

double length(particle p1, particle p2){
    double l;
    //double dx, dy;
    //dx = p2.x - p1.x; dx -= Lx*nearbyint(dx/Lx);
    //dy = p2.y - p1.y; dy -= Lx*nearbyint(dy/Lx);
    //l = sqrt(pow(dx, 2) + pow(dy, 2));
    l = sqrt(pow(p2.x - p1.x, 2) + pow(p2.y - p1.y, 2));
    return l;
}

double angle(particle p1, particle p2){
    double dx, dy, theta;
    dy = p2.y - p1.y;
    dx = p2.x - p1.x;
    theta = atan(dy/dx);

    if(dy<0 && dx<0){
        theta = theta + pi;
    }

    if(dy>0 && dx<0){
        theta = theta + pi;
    }

    return theta;
}

double fs(double x, double Lx,double v0){
    //return (v0/2.0)*(1+ cos(2*pi*x/Lx));
    return v0*(1+sin(4.0*pi*x/Lx));
}

double fs_drift(double x, double v0){
    return v0*(1+sin(4.0*pi*x/10.0));
    //return (v0/2.0)*(1+ cos(2*pi*x/200.0));
}

double fs_tethered(double x, double Lx){
    if(x<50.0){
        return 0.0;
    }
    else{
        return 0.01*(pow((x-50.0),1.0));
    }
}

void AssignActivity(polymer &p, vector<int> &v){
    for(int i=0; i<v.size(); i++){
        p.chain[v[i] - 1].alpha = 1;
    }
}

/*pair<double, double> SpringForces(polymer &p, int &i, double &ks, double &l0){
    double Fx, Fy, l;

    if(i!=0 && i!=p.N-1){
        l = length(p.chain[i], p.chain[i-1]);
        if(l>1e-6){
            Fx = ks*(l-l0)*(p.chain[i-1].x - p.chain[i].x)/l;
            Fy = ks*(l-l0)*(p.chain[i-1].y - p.chain[i].y)/l;
        }
        else{
            Fx = 0.0;
            Fy = 0.0;
        }

        l = length(p.chain[i], p.chain[i+1]);
        if(l>1e-6){
            Fx += ks*(l-l0)*(p.chain[i+1].x - p.chain[i].x)/l;
            Fy += ks*(l-l0)*(p.chain[i+1].y - p.chain[i].y)/l;
        }
        else{
            Fx += 0.0;
            Fy += 0.0;
        }
    
        return {Fx, Fy};
    }
    else if(i==0){
        l = length(p.chain[i], p.chain[i+1]);
        if(l>1e-6){
            Fx = ks*(l-l0)*(p.chain[i+1].x - p.chain[i].x)/l;
            Fy = ks*(l-l0)*(p.chain[i+1].y - p.chain[i].y)/l;
        }
        else{
            Fx = 0.0;
            Fy = 0.0;
        }

        return {Fx, Fy};
        
    }
    else{
        l = length(p.chain[i-1], p.chain[i]);
        if(l>1e-6){
            Fx = ks*(l-l0)*(p.chain[i-1].x - p.chain[i].x)/l;
            Fy = ks*(l-l0)*(p.chain[i-1].y - p.chain[i].y)/l;
        }
        else{
            Fx = 0.0;
            Fy = 0.0;
        }       

        return {Fx, Fy};
    }
}*/