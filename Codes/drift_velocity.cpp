// Program to calculate the average drift velocity of active flexible chains with each monomer either being active or passive, in a traveling activity wave

#include <iostream>
#include <fstream>
#include <string>
#include <array>
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
    double ks = 5.0, Dt = 0.01, Dr = 10.0, tau = 1.0/Dr, l0 = 0.0;
    double Lx = 10.0; double Ly = 10.0;
    double v0 = 1.0;


    //Intializing polymer chain
    int N = 2; //number of monomers
    polymer p(N);
    std::vector<int> alpha = {0,0,0,0,0,0,0,0,0,0}; //vector indicating the active/passive nature of the monomers
    std::vector<double> fric = {0,0,0,0,0,0,0,0,0,0}; //vector indicating the friction coefficients of the monomers

    for(int i=0; i<N;i++){
        alpha[i]=1;
    }
    //alpha[0]=1;  //Head active case

    for(int i=0; i<N;i++){
        fric[i]=1;
    }

    // Other variables
    double t; double tmax = 10000.0; double dt = 0.01;
    int N_run = 100;
    double eta_x, eta_y, eta_phi; double eta_or_x, eta_or_y;
    double l, Fx, Fy;

    // Averaging variables
    int timestep; double xc; double xc_0; int Nbins = 40; int count[Nbins]; double count_or[Nbins]; double lbin = Lx/double(Nbins); int timesteps_per_snapshot= 100;
    double vw, vd, vd_avg;
    const size_t steps = 40;
    double stepsize, vw_i, vw_f;
    vw_i = 0.001;
    vw_f = 0.3;

    std::array<double, steps+1> vw_array{};
    stepsize = std::log10(vw_f/vw_i)/(steps-1);
    for (int i=0; i<steps; ++i) {vw_array[i+1] = std::pow(10.0, std::log10(vw_i) + i*stepsize);}


    int style; style = 0;//This variable determines the connectivity of the polymer - for linear chain - style = 0, for ring - style = 1, for star - style = 2, and for clique - style = 3


    // Opening files
    ofstream out1;
    out1.open("/home/valechbh/Downloads/Polymers/drift_N_"+to_string(N)+"_AOUP_test.dat");


    // loop over different values of vw
    for(int j = 0; j < steps+1; j++){
    //for(double vw=0, increment = 0.001, counter = 1; vw <= 1.0; vw += increment){
	    vw = vw_array[j];
        vd_avg = 0.0;

        // loop over N_run independent runs
        for(int i = 0; i<N_run; i++){

            // initialize the polymer
            for(int i=0; i<N; i++){
                p.chain[i].x = 100.0 - i*uniform_rand();
                p.chain[i].y = 100.0 - i*uniform_rand();
                p.chain[i].phi = 2*pi*uniform_rand();
                p.chain[i].or_x = 0.0;
                p.chain[i].or_y = 0.0;
            }
            xc_0 = CentreOfFriction(p,fric);
            polymer p_new = p;

            t = 0.0;
            timestep = 0;
            // Time loop
            while(t<tmax){
                timestep++;
                t = t+dt;

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

                        p_new.chain[i].x += (dt/fric[i])*(Fx + alpha[i]*fs_drift((p.chain[i].x-vw*t),v0)*p.chain[i].or_x) + sqrt(2*Dt*(dt/fric[i]))*eta_x;
                        p_new.chain[i].y += (dt/fric[i])*(Fy + alpha[i]*fs_drift((p.chain[i].x-vw*t),v0)*p.chain[i].or_y) + sqrt(2*Dt*(dt/fric[i]))*eta_y;

                        p_new.chain[i].or_x += dt*(-p_new.chain[i].or_x/tau) + sqrt(dt/tau)*eta_or_x;
                        p_new.chain[i].or_y += dt*(-p_new.chain[i].or_y/tau) + sqrt(dt/tau)*eta_or_y;
                    }
                    p.chain = p_new.chain;
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

                        p_new.chain[i].x += (dt/fric[i])*(Fx + alpha[i]*fs_drift((p.chain[i].x-vw*t),v0)*p.chain[i].or_x) + sqrt(2*Dt*(dt/fric[i]))*eta_x;
                        p_new.chain[i].y += (dt/fric[i])*(Fy + alpha[i]*fs_drift((p.chain[i].x-vw*t),v0)*p.chain[i].or_y) + sqrt(2*Dt*(dt/fric[i]))*eta_y;

                        p_new.chain[i].or_x += dt*(-p_new.chain[i].or_x/tau) + sqrt(dt/tau)*eta_or_x;
                        p_new.chain[i].or_y += dt*(-p_new.chain[i].or_y/tau) + sqrt(dt/tau)*eta_or_y;
                    }
                    p.chain = p_new.chain;
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

                        p_new.chain[i].x += (dt/fric[i])*(Fx + alpha[i]*fs_drift((p.chain[i].x-vw*t),v0)*p.chain[i].or_x) + sqrt(2*Dt*(dt/fric[i]))*eta_x;
                        p_new.chain[i].y += (dt/fric[i])*(Fy + alpha[i]*fs_drift((p.chain[i].x-vw*t),v0)*p.chain[i].or_y) + sqrt(2*Dt*(dt/fric[i]))*eta_y;

                        p_new.chain[i].or_x += dt*(-p_new.chain[i].or_x/tau) + sqrt(dt/tau)*eta_or_x;
                        p_new.chain[i].or_y += dt*(-p_new.chain[i].or_y/tau) + sqrt(dt/tau)*eta_or_y;
                    }
                    p.chain = p_new.chain;
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

                        p_new.chain[i].x += (dt/fric[i])*(Fx + alpha[i]*fs_drift((p.chain[i].x-vw*t),v0)*p.chain[i].or_x) + sqrt(2*Dt*(dt/fric[i]))*eta_x;
                        p_new.chain[i].y += (dt/fric[i])*(Fy + alpha[i]*fs_drift((p.chain[i].x-vw*t),v0)*p.chain[i].or_y) + sqrt(2*Dt*(dt/fric[i]))*eta_y;

                        p_new.chain[i].or_x += dt*(-p_new.chain[i].or_x/tau) + sqrt(dt/tau)*eta_or_x;
                        p_new.chain[i].or_y += dt*(-p_new.chain[i].or_y/tau) + sqrt(dt/tau)*eta_or_y;
                    }
                    p.chain = p_new.chain;
                }

            }
            xc = CentreOfFriction(p,fric); //center of mass at the end of the trajectory
            vd = (xc - xc_0)/t; //drift velocity
            vd_avg += vd;
        }
        out1 << vw << "\t" << (vd_avg)/(N_run) << endl;
    }
}
