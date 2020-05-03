#ifndef __CHECKDERIVATIES__
#define __CHECKDERIVATIES__

#include "BOV.h"
#include "neighborhood_search.h"
#include "kernel.h"
#include <time.h>
#include <math.h>
#include <stdbool.h>

#define M_PI 3.14159265358979323846

// extern int NPTS_DOMAIN;
// extern int NPTS_BOUNDARIES;

enum fieldNames {
  Density = 1,
  Velocity = 2,
  Pressure = 3,
  Temperature = 4
};

typedef enum particleType {
  ON_BOUNDARY = 1,
  IN_DOMAIN = 2,
} particleType;

typedef enum kernelType {
  CUBIC = 1,
  LUCY = 2,
  QUARTIC = 3,
  QUINTIC = 4
} kernelType;


typedef struct singleParticleDerivatives {
  double* divergence;
  double* gradient;
  double* laplacian;
} singleParticleDerivatives;

typedef struct mySingleParticle {
  double* coordinates; // x-y
  double* values; // density or velocity (vector!) or pressure or temperature or anything else
  int size_values;
  double mass;
  double density;

  neighborhood_options* particle_neighbours; // all the neighbouring particles infos related to the particle of interest
  singleParticleDerivatives* particle_derivatives; // divergence, gradient and laplacian (depending if the quantity carried by the particle is a scalar or vector)
} mySingleParticle;

typedef struct allParticles {
  mySingleParticle* array_of_particles;
  int nb_particles;
  int nb_particles_domain;
  int nb_particles_boundaries;
  int* index_part_in_domain; // array of the indices of all particles inside the domain
  int* index_part_on_boundaries; // array of the indices of all particles on the boundaries
//   int* index_part_global; // array of the indices of all particles (domain + boundaries)
} allParticles;

// Create an array of particles "mySingleParticle" by initializing the different 
mySingleParticle* create_array_of_particles(int nbParticles, int size_values, neighborhood* nh);

allParticles* create_all_particles(int* nbParticles, double* domain_lim, int size_values, neighborhood* nh,  particleType partType, double* args_init_fct, double (*f_init)(double*, double*));

allParticles* create_particles_in_domain(int* nbParticles, double* domain_lim, int size_values, double* args_init_fct, double (*f_init)(double*, double*), int starting_index);

allParticles* create_particles_on_boundaries(int* nbParticles, double* boundaries_lim, int size_values, int nbBoundaries, double* value_on_boundaries, int starting_index); 

// Initialize the elements of a singleParticleDerivatives
singleParticleDerivatives* initialize_particle_derivatives(int size_values);

// Assign coordinates, values, mass and density to one particle located on a 1-D segment whose limits are in x_lim, with index "index_art". 
void init1DSegmentWithParticles(double* x_lim, double* coord, double* values, double* mass, double* density, int nb_particles, int index_part, int size_values);

void initValuesParticles(double* coord, double *values, double* args, int size_values, double (*f)(double*, double*));

void initSquareWithParticles(double* x_lim, double* coord, double* values, double* mass, double* density, int nbPart_x, int nbPart_y, int index_x, int index_y, int size_values);

// Function assigning the values of the quantity carried by each particle based on the coordinates
double myFunctionToDerive(double* x, double* args);

double initFunction(double* x, double* args); 

double initFunction_GaussianSource(double* x, double* args);

// Compute the derivatives (divergence, gradients, laplacian) of a single particle with index "index_part"
// ---> Use for that the functions implemented in "neighborhood_search" and "kernel" 
void computeDerivativesOfParticleQuantity(mySingleParticle* myPart, int index_part, kernelType kType);

// Loop on all particles and compute the derivatives of the quantity for each of them
void computeDerivatiesAllParticles(mySingleParticle* myPart, int nbParticles, kernelType kType);

void computeDerivativesAllParticles(allParticles* myPart, kernelType kType);

void delete_all_particles(allParticles* myPart);

void neighborhood_update_new(neighborhood_options* options, neighborhood* nh, allParticles* myPart, int iterations);

allParticles* combine_two_particles_sets(allParticles* part1, allParticles* part2, particleType pType1, particleType pType2);

void clone_single_particle (mySingleParticle* part1, mySingleParticle* part2);

void associate_neighborhood_to_particles(allParticles* part, neighborhood* nh);

void integrate_equation(allParticles* part, double dt, double alpha);

double solution_Fourier_series_Gaussian_source(double *x_coord, double time, double* args_fct_init, int nb_terms);

double get_A_mn_series_Gaussian_source(int m, int n, double* args_fct_init);

double get_integration_Gaussian_source(int m, double x_c, double length_x, double sigma);

double fct_integration_Gaussian_source(double x, int m, double x_c, double length_x, double sigma);

#endif