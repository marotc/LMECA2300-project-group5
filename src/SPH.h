#ifndef SPH_H
#define SPH_H
#include "particle.h"
#include "print_particules.h"
#include "kernel.h"
#include "derivatives.h"
#include <time.h>

typedef struct Setup Setup;
typedef struct Residual Residual;
typedef enum Free_surface_detection Free_surface_detection;

//CSF       : if  ||n_i|| > threshold => particle i belongs to interface
//DIVERGENCE: if div(pos) < threshold => particle i belongs to interface
//NONE      : no free surface detection
enum Free_surface_detection {CSF,DIVERGENCE,NONE};

struct Setup{
	int itermax;
	double timestep;
	double kh;
	Search* search;
	Kernel kernel;
	//int nThreads;
	Free_surface_detection free_surface_detection;  // Strategy to estimate if a particle should be considered on the free surface or not
	double interface_threshold; // Threshold for the detection of particles on the free surface
	double XSPH_epsilon; // Parameter between 0 and 1 multiplying the XSPH correction;0 if no correction wanted
};

struct Residual {
	double mass_eq;
	double momentum_x_eq;
	double momentum_y_eq;
};



Setup* Setup_new(int iter, double timestep,double kh, Search* search, Kernel kernel,Free_surface_detection free_surface_detection,double interface_threshold, double XSPH_epsilon);
void Setup_free(Setup* setup);

Residual* Residual_new();
void free_Residuals(Residual** residuals,int N);

void simulate_with_boundaries(Grid* grid, Particle** particles, int n_p,
	Setup* setup, Animation* animation, double eps);

void drift(Particle** particles, int N, double dt);
void density_pressure_update(Particle** particles, int N, double kh, Kernel kernel);
void apply_BC_Adami_all(Particle **particles, int N, Kernel kernel, double kh);
void apply_BC_Adami(Particle* pi, Kernel kernel, double kh);
void inter_particle_forces_pres_all(Particle** particles, int N, double kh, Kernel kernel, double eps);
void inter_particle_forces_pres(Particle* pi, double kh, Kernel kernel, double eps);
void kick(Particle **particles, int N, double dt, Kernel kernel, double kh);

double compute_admissible_dt(double safety_param, double h_p, double c_0, double rho_0, double mu, double sigma, xy* g, double U);
#endif
