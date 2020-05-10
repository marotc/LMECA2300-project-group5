#ifndef PARTICLE_H
#define PARTICLE_H

#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <omp.h>
#include "utils.h"

typedef struct Particle Particle;
typedef struct Particle_derivatives Particle_derivatives;
typedef struct Physical_parameters Physical_parameters;
typedef struct Cell Cell;
typedef struct Grid Grid;
typedef struct Search Search;
typedef struct Boundary Boundary;

struct Cell {
	int i,j;
	List* neighboring_cells;
	List* particles;
	bool visited; // for improved algorithm
};

struct Grid {
	int nCellx;
	int nCelly;
	Cell*** cells;	// 2D array of cell pointers
	double h; // side length of a cell
	double left;	// x-coordinate of left side
	double right;	// x-coordinate of right side
	double top;		// y-coordinate of top side
	double bottom;	// y-coordinate of bottom side
};

struct Physical_parameters {
	double nu;
	double rho_0;
	double gamma;
	double sound_speed;
	double sigma;
	double background_p;
	xy* gravity;
};

struct Particle {
	int index;
	double m;     // mass
	xy* pos;      // position
	xy* v;        // velocity
	xy *v_imp;	  // imposed velocity
	double rho;   // density
	double P;     // pressure
	double Cs;    // color field
	xy* normal;   // normal vector
	double kappa; // curvature
	double V;     // inverse of particle volume

	xy *a, *as; // momentum and drift acceleration
	xy *vs; // drift velocity (v tilde)

	xy* XSPH_correction; // Correction on the velocity field when updating the position of the particles
	bool on_free_surface; // boolean to know if particles is on the free surface (used for visualization)
	bool on_boundary; // boolean to know if particles is on the boundary
	Physical_parameters* param; // physical parameters associated to the particle

	Cell* cell;    // cell that the particle belongs to
	List* neighborhood; // list of neighbors
	List* potential_neighborhood; // list of potential neighbors (for Verlet)
};

struct Particle_derivatives {
	int index;
	double div_v;
	xy *grad_P;
	xy *lapl_v;
	xy *grad_Cs;
	double lapl_Cs;
};

struct Search {
	bool verlet;
	double kh;
	double L;
	int T;
};

struct Boundary {
    int nb_part_on_bound;
	Particle** part_on_bound;
	xy* v_imposed;
	xy* acc_imposed;
};

// Grid
void Cell_free(Cell* cell); // Cell destructor

Grid* Grid_new(double x1, double x2, double y1, double y2, double kh); // Grid constructor
void Grid_free(Grid* grid); // Grid destructor

// Particle
Particle* Particle_new(int index, double m, xy* pos, xy* v, xy *v_imp, double rho_0, double nu, double c_0, double gamma, double sigma, double background_p, xy* gravity, bool on_boundary);
void Particle_free(Particle* particle);
void free_particles(Particle** particles, int N);

double Particle_get_P(Particle *particle);
xy * Particle_get_v(Particle *particle);
xy * Particle_get_pos(Particle *particle);
double Particle_get_v_x(Particle *particle);
double Particle_get_v_y(Particle *particle);
double Particle_get_Cs(Particle *particle);
xy * Particle_get_normal(Particle *particle);

// Particle_derivatives
Particle_derivatives* Particle_derivatives_new(int index);
void Particle_derivatives_free(Particle_derivatives* particle_derivatives);
void Particle_derivatives_reset(Particle_derivatives *particle_derivatives);
void free_particles_derivatives(Particle_derivatives** particles_derivatives, int N);

Search* Search_new(bool verlet, double kh, double L, int T);
void Search_free(Search* search);

// Update of the particles' locations in cells
Cell* localize_particle(Grid* grid, Particle * p);
void update_cells(Grid* grid, Particle** particles, int N);

// Update of the neighborhood
void reset_grid(Grid* grid);
void reset_particles(Particle** particles, int N, int iter, Search* search);
void add_neighbors_from_cell(Particle* p, Cell* cell, double r);
void add_neighbors_from_cells(Grid* grid, Particle* p);
void update_neighborhoods_particles(Grid* grid, Particle** particles, int N, Search* search);
void update_from_potential_neighbors(Particle** particles, int N, Search* search);
void update_neighborhoods(Grid* grid, Particle** particles, int N, int iter, Search* search);


#endif
