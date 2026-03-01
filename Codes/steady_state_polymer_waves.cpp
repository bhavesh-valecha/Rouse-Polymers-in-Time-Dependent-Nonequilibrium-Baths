// Program to generate the steady state profiles for active flexible chains with each monomer either being active or passive in traveling activity waves

#include <iostream>
#include <fstream>
#include <string>
#include "particle_polymer_waves.h"
#include "random_polymer_waves.h"
#include "Force_chain.h"
#include "Force_ring.h"
#include "Force_star.h"
#include "Force_clique.h"

using namespace std;

int main(int argc, char** argv){
    seed_drand48();

    // Parameters
    double ks = 5.0, Dt = 0.001, Dr = 10.0, tau = 1.0/Dr, l0 = 0.0;
    double Lx = 10.0; double Ly = 10.0;
    double v0 = 1.0; double vw = 0.01;
    

    //Intializing polymer chain
    int N = 8; //number of monomers
    polymer p(N);
    std::vector<int> alpha = {0,0,0,0,0,0,0,0,0,0}; //vector indicating the active/passive nature of the monomers
    std::vector<double> fric = {0,0,0,0,0,0,0,0,0,0}; //vector indicating the friction coefficients of the monomers

    for(int i=0; i<N;i++){
        alpha[i]=1;
    }

    for(int i=0; i<N;i++){
        fric[i]=1;
    }

    for(int i=0; i<N; i++){
        p.chain[i].x = 9.0 - i*1.0;
        p.chain[i].y = 9.0 - i*1.0;
        p.chain[i].phi = 2*pi*uniform_rand();
        p.chain[i].or_x = 0.0;
        p.chain[i].or_y = 0.0;
    }
    polymer p_new = p;

    // Other variables
    double t = 0.0; double tmax = 5000000.0; double dt = 0.01; 
    double eta_x, eta_y, eta_phi; double eta_or_x, eta_or_y;
    double l, Fx, Fy;
    

    // Averaging variables
    int timestep=0; double xc; int Nbins = 40; int count[Nbins]; double count_or[Nbins]; double lbin = Lx/double(Nbins); int j, j1; int timesteps_per_snapshot= 100;
    for(int k = 0; k<Nbins; k++){
        count[k] = 0;
        count_or[k] = 0.0;
    }
    
    int style; style = 0;//This variable determines the connectivity of the polymer - for linear chain - style = 0, for ring - style = 1, for star - style = 2, and for clique - style = 3

    // Opening files
    ofstream out1;
    out1.open("/home/valechbh/Downloads/Polymers/density_N_"+to_string(N)+"_AOUPol_test.dat");

    // Time loop
    
    while(t<tmax){
        timestep++;
        t = t+dt;

        //LINEAR CHAIN
        if(style==0){
            for(int i = 0; i<N; i++){
                Fx = 0.0; Fy = 0.0;
                Fx = SpringForces_chain(p, i, ks, l0, N).first;
                Fy = SpringForces_chain(p, i, ks, l0, N).second;
                
                eta_x = Gaussian();
                eta_y = Gaussian();
                eta_phi = Gaussian();
                eta_or_x = Gaussian();
                eta_or_y = Gaussian();

                p_new.chain[i].x += (dt/fric[i])*(Fx + alpha[i]*fs((p.chain[i].x-vw*t),Lx,v0)*p.chain[i].or_x) + sqrt(2*Dt*(dt/fric[i]))*eta_x;
                p_new.chain[i].y += (dt/fric[i])*(Fy + alpha[i]*fs((p.chain[i].x-vw*t),Lx,v0)*p.chain[i].or_y) + sqrt(2*Dt*(dt/fric[i]))*eta_y;

                p_new.chain[i].or_x += dt*(-p_new.chain[i].or_x/tau) + sqrt(dt/tau)*eta_or_x;
                p_new.chain[i].or_y += dt*(-p_new.chain[i].or_y/tau) + sqrt(dt/tau)*eta_or_y; 
            }
            p.chain = p_new.chain;

            if(timestep%timesteps_per_snapshot == 0){
                xc = CentreOfFriction(p,fric);
                xc = xc - vw*t; //go to the comoving frame      
                while(xc<0) {xc += Lx;}
                while(xc>Lx) {xc -= Lx;}
                j = int(xc/lbin); 
                count[j]++;
            }
        }


        //RING
        else if(style==1){
            for(int i = 0; i<N; i++){
                Fx = 0.0; Fy = 0.0;
                Fx = SpringForces_ring(p, i, ks, l0, N).first;
                Fy = SpringForces_ring(p, i, ks, l0, N).second;
                
                eta_x = Gaussian();
                eta_y = Gaussian();
                eta_phi = Gaussian();
                eta_or_x = Gaussian();
                eta_or_y = Gaussian();

                p_new.chain[i].x += (dt/fric[i])*(Fx + alpha[i]*fs((p.chain[i].x-vw*t),Lx,v0)*p.chain[i].or_x) + sqrt(2*Dt*(dt/fric[i]))*eta_x;
                p_new.chain[i].y += (dt/fric[i])*(Fy + alpha[i]*fs((p.chain[i].x-vw*t),Lx,v0)*p.chain[i].or_y) + sqrt(2*Dt*(dt/fric[i]))*eta_y;

                p_new.chain[i].or_x += dt*(-p_new.chain[i].or_x/tau) + sqrt(dt/tau)*eta_or_x;
                p_new.chain[i].or_y += dt*(-p_new.chain[i].or_y/tau) + sqrt(dt/tau)*eta_or_y; 
            }
            p.chain = p_new.chain;
            
            if(timestep%timesteps_per_snapshot == 0){
                xc = CentreOfFriction(p,fric);
                xc = xc - vw*t; //go to the comoving frame         
                while(xc<0) {xc += Lx;}
                while(xc>Lx) {xc -= Lx;}
                j = int(xc/lbin); 
                count[j]++;
            }
        }

        //STAR
        else if(style==2){
            for(int i = 0; i<N; i++){
                Fx = 0.0; Fy = 0.0;
                Fx = SpringForces_star(p, i, ks, l0, N).first;
                Fy = SpringForces_star(p, i, ks, l0, N).second;
                
                eta_x = Gaussian();
                eta_y = Gaussian();
                eta_phi = Gaussian();
                eta_or_x = Gaussian();
                eta_or_y = Gaussian();

                p_new.chain[i].x += (dt/fric[i])*(Fx + alpha[i]*fs((p.chain[i].x-vw*t),Lx,v0)*p.chain[i].or_x) + sqrt(2*Dt*(dt/fric[i]))*eta_x;
                p_new.chain[i].y += (dt/fric[i])*(Fy + alpha[i]*fs((p.chain[i].x-vw*t),Lx,v0)*p.chain[i].or_y) + sqrt(2*Dt*(dt/fric[i]))*eta_y;

                p_new.chain[i].or_x += dt*(-p_new.chain[i].or_x/tau) + sqrt(dt/tau)*eta_or_x;
                p_new.chain[i].or_y += dt*(-p_new.chain[i].or_y/tau) + sqrt(dt/tau)*eta_or_y; 
            }
            p.chain = p_new.chain;
            
            if(timestep%timesteps_per_snapshot == 0){
                xc = CentreOfFriction(p,fric);
                xc = xc - vw*t; //go to the comoving frame         
                while(xc<0) {xc += Lx;}
                while(xc>Lx) {xc -= Lx;}
                j = int(xc/lbin); 
                count[j]++;
            }
        }

        //CLIQUE
        else if(style==3){
            for(int i = 0; i<N; i++){
                Fx = 0.0; Fy = 0.0;
                Fx = SpringForces_clique(p, i, ks, l0, N).first;
                Fy = SpringForces_clique(p, i, ks, l0, N).second;
                
                eta_x = Gaussian();
                eta_y = Gaussian();
                eta_phi = Gaussian();
                eta_or_x = Gaussian();
                eta_or_y = Gaussian();

                p_new.chain[i].x += (dt/fric[i])*(Fx + alpha[i]*fs((p.chain[i].x-vw*t),Lx,v0)*p.chain[i].or_x) + sqrt(2*Dt*(dt/fric[i]))*eta_x;
                p_new.chain[i].y += (dt/fric[i])*(Fy + alpha[i]*fs((p.chain[i].x-vw*t),Lx,v0)*p.chain[i].or_y) + sqrt(2*Dt*(dt/fric[i]))*eta_y;

                p_new.chain[i].or_x += dt*(-p_new.chain[i].or_x/tau) + sqrt(dt/tau)*eta_or_x;
                p_new.chain[i].or_y += dt*(-p_new.chain[i].or_y/tau) + sqrt(dt/tau)*eta_or_y; 
            }
            p.chain = p_new.chain;
            
            if(timestep%timesteps_per_snapshot == 0){
                xc = CentreOfFriction(p,fric);
                xc = xc - vw*t; //go to the comoving frame         
                while(xc<0) {xc += Lx;}
                while(xc>Lx) {xc -= Lx;}
                j = int(xc/lbin); 
                count[j]++;
            } 
        }
    }

    double b = tmax/dt/timesteps_per_snapshot;
    for(int k=0; k<Nbins; k++){
        out1 << lbin*k + lbin/2.0 << "\t" << double((count[k]/b)*Nbins) << endl;
    }
}