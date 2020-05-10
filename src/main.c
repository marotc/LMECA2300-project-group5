#include "print_particules.h"
#include "particle.h"
#include "SPH.h"
#include "derivatives.h"
#include <math.h>
#include "kernel.h"
//#include "crtdbg.h" // for memory leak detection; comment if you're on Linux

void script_csf();
void script_csf_circle();
void script_circle_to_ellipse();
void script_csf_circle_paper();
void script_lid_driven_cavity();
void script1();
void script2();

int main() {

	//_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF); // comment if on Linux
	//script_csf_circle();
	//script_circle_to_ellipse();
	//script_csf_circle_paper();
	script_lid_driven_cavity();
	//script2();
	return EXIT_SUCCESS;
}

void fill_interior_rectangular(int n_x, int n_y, double start_x, double start_y, double h_x, double h_y, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals,
	double mass, double rho_0, double mu, double c_0, double gamma, double background_pressure, xy* gravity) {
	int k = 0;
	for (int i = 0; i < n_x; i++) {
		for (int j = 0; j < n_y; j++) {
			xy *pos = xy_new(start_x + i * h_x, start_y + j * h_y);
			xy *v = xy_new(0, 0); // initial velocity = 0
			xy* gravityP = xy_new(gravity->x, gravity->y);
			particles[k] = Particle_new(k, mass, pos, v, rho_0, mu, c_0, gamma, 0.0, background_pressure, gravityP, false);
			particles_derivatives[k] = Particle_derivatives_new(k);
			residuals[k] = Residual_new();
			k++;
		}
	}
}

/*
#include "omp.h"
#include <stdio.h>
int main()
{
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);
	omp_set_num_threads(6);
	#pragma omp parallel
	{
		int id = omp_get_thread_num();
		printf("Hello (%d)", id);
		printf("World (%d)\n", id);
	}
	return 0;
}*/
// Evolution of a 2D square submitted to the tension surface force
/*void script_csf() {

	// Parameters of the problem
	double L = 2; // size of the domain: [-L,L] x [-L,L]
	double l = 1; // size of the square: [-l,l] x [-l,l]
	double dt = 0.001; // physical time step
	double T = 100; // duration of simulation
	// double T = dt; // duration of simulation

	// Physical parameters
	double rho_0 = 998.0; // initial (physical) density of water at 20°C (in kg/m^3)
	double mu = 1.0016e-3; // dynamic viscosity of water at 20°C (in N.s/m^2)
	double gamma = 7; // typical value for liquid (dimensionless)
	double c_0 = 1.0;//1481; // sound speed in water at 20°C (in m/s)
	double sigma = 72.86e-3; // surface tension of water-air interface at 20°C (in N/m)
	double background_pressure = 0.0;
	xy* gravity = xy_new(0.0, 0.0);

	// SPH parameters
	int n_per_dim = 101; // number of particles per dimension
	double kh = sqrt(21) * 2 * l / n_per_dim; // kernel width to ensure 21 particles in the neighborhood
	int n_iter = (int)(T / dt); // number of iterations to perform
	Kernel kernel = Cubic; // kernel choice
	double interface_threshold = 1.5;//20; // If ||n_i|| > threshold => particle i belongs to interface (first detection approach)
	Verlet *verlet = NULL; // don't use Verlet (for now)
	double XSPH_epsilon = 0.5;
	Free_surface_detection surface_detection = DIVERGENCE;


	printf("n_iter = %d\n", n_iter);

	// Animation parameter
	double T_anim = 10; // duration of animation
	double dt_anim = T_anim / n_iter; // time step of animation

	// Initialize particles on a square
	int n_p = squared(n_per_dim); // total number of particles
	double h = 2 * l / (n_per_dim - 1); // step between neighboring particles
	double m = rho_0 * h*h;

	Particle** particles = (Particle**)malloc(n_p * sizeof(Particle*));
	Particle_derivatives** particles_derivatives = malloc(n_p * sizeof(Particle_derivatives*));
	Residual** residuals = malloc(n_p * sizeof(Residual*));
	for (int i = 0; i < n_per_dim; i++) {
		for (int j = 0; j < n_per_dim; j++) {
			int index = i * n_per_dim + j;
			xy *pos = xy_new(-l + i * h, -l + j * h);
			xy *v = xy_new(0, 0); // initial velocity = 0
			particles[index] = Particle_new(index, m, pos, v, rho_0, mu, c_0, gamma, sigma, background_pressure, gravity);
			particles_derivatives[index] = Particle_derivatives_new(index);
			residuals[index] = Residual_new();
		}
	}

	// Setup grid
	Grid *grid = Grid_new(-L, L, -L, L, kh);
	// Setup setup
	Setup *setup = Setup_new(n_iter, dt, kh, verlet, kernel, surface_detection, interface_threshold, XSPH_epsilon);
	// Setup animation
	Animation *animation = Animation_new(n_p, dt_anim, grid, 1);

	// Start simulation
	simulate(grid, particles, particles_derivatives, residuals, n_p, update_positions_seminar_5, setup, animation);

	// Free stuff
	free_particles(particles, n_p);
	free_particles_derivatives(particles_derivatives, n_p);
	free_Residuals(residuals, n_p);
	Grid_free(grid);
	Setup_free(setup);
	Animation_free(animation);

}


// Evolution of a 2D circle submitted to the tension surface force
void script_csf_circle() {

	// Computational Parameters
	double L = 100; // size of the domain: [-L,L] x [-L,L]
	double l = 50; // size of the square: [-l,l] x [-l,l]
	int n_per_dim = 51; // number of particles per dimension
	double kh = 0.1*L; // kernel width
	double dt = 1.0; // physical time step
	double dt_anim = 0.001; // time step for animation
	int n_iter = 2000; // number of iterations to perform
	Verlet *verlet = NULL; // don't use Verlet (for now)
	Kernel kernel = Cubic; // kernel choice
	double interface_threshold = 0.08; // If ||n_i|| > threshold => particle i belongs to interface (first detection approach)
	double XSPH_epsilon = 0.5;
	Free_surface_detection surface_detection = CSF;

	// Physical parameters
	double rho_0 = 998.0; // initial (physical) density of water at 20°C (in kg/m^3)
	double mu = 1.0016e-3; // dynamic viscosity of water at 20°C (in N.s/m^2)
	double gamma = 7.0; // typical value for liquid (dimensionless)
	double c_0 = 1.0;//1481; // sound speed in water at 20°C (in m/s)
	double sigma = 72.86e-3; // surface tension of water-air interface at 20°C (in N/m)
	double background_pressure = 0.0;
	xy* gravity = xy_new(0.0, 0.0);

	// Initialize particles in a circle
	int n_p = squared(n_per_dim); // total number of particles
	double h = 2 * l / (n_per_dim - 1); // step between neighboring particles
	double m = rho_0 * M_PI * l * l / n_p;
	double deltat_max = 0.25 * h / c_0;
	// 	printf(">>>>>>>>>>>>>> deltat_max = %2.6f \n", deltat_max);

	Particle** particles = (Particle**)malloc(n_p * sizeof(Particle*));
	Particle_derivatives** particles_derivatives = malloc(n_p * sizeof(Particle_derivatives*));
	Residual** residuals = (Residual**)malloc(n_p * sizeof(Residual*));

	int k = 0;
	int alpha = 2;
	int nb = (int)(alpha*n_per_dim);
	for (int i = 0; i < n_per_dim; i++) {
		for (int j = 0; j < n_per_dim; j++) {
			int index = i * n_per_dim + j;
			xy *pos = generate_circle(k, n_p, nb, l);
			// 			xy *pos_circle = map_to_circle(pos);
			// 			pos->x = pos_circle->x;
			// 			pos->y = pos_circle->y;
			xy *v = xy_new(0, 0); // initial velocity = 0
			particles[index] = Particle_new(index, m, pos, v, rho_0, mu, c_0, gamma, sigma, background_pressure, gravity);
			particles_derivatives[index] = Particle_derivatives_new(index);
			residuals[index] = Residual_new();
			k++;
		}
	}


	// Setup grid
	Grid *grid = Grid_new(-L, L, -L, L, kh);
	// Setup animation
	Animation *animation = Animation_new(n_p, dt_anim, grid, 1);
	// Setup setup
	Setup *setup = Setup_new(n_iter, dt, kh, verlet, kernel, surface_detection, interface_threshold, XSPH_epsilon);

	// Start simulation
	simulate(grid, particles, particles_derivatives, residuals, n_p, update_positions_seminar_5, setup, animation);

	// Free stuff
	free_particles(particles, n_p);
	free_particles_derivatives(particles_derivatives, n_p);
	free_Residuals(residuals, n_p);
	Grid_free(grid);
	Setup_free(setup);
	Animation_free(animation);

}

// Evolution of a 2D circle with non-zero initial velocities (no surface tension force)
// Test case from "Simulating Free Surface Flows with SPH", Monaghan (1994)
void script_circle_to_ellipse() {

	// Parameters of the problem
	double l = 1.0; // radius of the circle
	double L = 2.0*l; // size of the domain: [-L,L] x [-L,L]
	double dt = 1.0e-5; // physical time step
	double T = 0.0076; // duration of simulation

	// Physical parameters
	double rho_0 = 1000.0; // initial (physical) density of water at 20°C (in kg/m^3)
	double mu = 1.0016e-3; // dynamic viscosity of water at 20°C (in N.s/m^2)
	double gamma = 7.0; // typical value for liquid (dimensionless)
	double c_0 = 1400.0;//1481; // sound speed in water at 20°C (in m/s)
	double sigma = 0.0; // surface tension of water-air interface at 20°C (in N/m)
	double background_pressure = 0.0;
	xy* gravity = xy_new(0.0, 0.0);

	// SPH parameters
	Verlet *verlet = NULL; // don't use Verlet (for now)
	Kernel kernel = Cubic; // kernel choice
	double interface_threshold = 1.5;//1000.0; // If ||n_i|| > threshold => particle i belongs to interface (first detection approach)
	double XSPH_epsilon = 0.5;
	Free_surface_detection surface_detection = DIVERGENCE;
	int N_c = 30; // number of circonferences on which points are placed
	int N_p = 6; // number of points on the first circonference (doubled for every circonference)
	int N_tot = 1; // total number of points
	for (int i = 1; i < N_c; i++) {
		N_tot += i * N_p;
	}
	printf("N_tot = %d \n", N_tot);
	int n_iter = (int)(T / dt); // number of iterations to perform
	double kh = 0.2*l;// is ideal to reach t = 0.0076; // kernel width


	// Animation parameter
	double T_anim = 0.1; // duration of animation
	double dt_anim = T_anim / n_iter; // time step of animation

	// Initialize particles in a circle
	double m = rho_0 * M_PI * l * l / N_tot; // mass of each particle

	Particle** particles = (Particle**)malloc(N_tot * sizeof(Particle*));
	Particle_derivatives** particles_derivatives = malloc(N_tot * sizeof(Particle_derivatives*));
	Residual** residuals = malloc(N_tot * sizeof(Residual*));

	// parameters defining the circle
	double b, delta_s, k, theta;
	delta_s = l / ((double)N_c - 1.0);
	theta = (2 * M_PI) / ((double)N_p);

	int index = 0;
	for (int i = 0; i < N_c; i++) {
		b = i;
		if (b == 0) {
			xy *pos = xy_new(0.0, 0.0);
			xy *v = xy_new(0.0, 0.0);
			particles[index] = Particle_new(index, m, pos, v, rho_0, mu, c_0, gamma, sigma, background_pressure, gravity);
			particles_derivatives[index] = Particle_derivatives_new(index);
			residuals[index] = Residual_new();
			index++;
		}
		else {
			for (int j = 0; j < i*N_p; j++) {
				k = (double)j / b;
				xy *pos = xy_new(b*delta_s*cos(k*theta), b*delta_s*sin(k*theta));
				xy *v = xy_new(-100.0*pos->x, 100.0*pos->y);
				particles[index] = Particle_new(index, m, pos, v, rho_0, mu, c_0, gamma, sigma, background_pressure, gravity);
				particles_derivatives[index] = Particle_derivatives_new(index);
				residuals[index] = Residual_new();
				index++;
			}
		}
	}


	// Setup grid
	Grid *grid = Grid_new(-L, L, -L, L, kh);
	// Setup animation
	Animation *animation = Animation_new(N_tot, dt_anim, grid, 1);
	// Setup setup
	Setup *setup = Setup_new(n_iter, dt, kh, verlet, kernel, surface_detection, interface_threshold, XSPH_epsilon);

	// Start simulation
	simulate(grid, particles, particles_derivatives, residuals, N_tot, update_positions_ellipse, setup, animation);

	// Free stuff
	free_particles(particles, N_tot);
	free_particles_derivatives(particles_derivatives, N_tot);
	free_Residuals(residuals, N_tot);
	Grid_free(grid);
	Setup_free(setup);
	Animation_free(animation);
}

// Evolution of a 2D circle : test case based on "Modelling surface tension of two-dimensional droplet using
// Smoothed Particle Hydrodynamics", Nowoghomwenma, 2018

void script_csf_circle_paper() {

	//	// Nowoghomwenma's paper
	// 	// Parameters of the problem
	// 	double l = 1e-3; // radius of the circle
	// 	double L = 2.0*l; // size of the domain: [-L,L] x [-L,L]
	// 	double dt = 0.000001; // physical time step
	// 	double T = 0.01; // duration of simulation
	// 	
	// 	// Physical parameters
	// 	double rho_0 = 1000.0; // initial (physical) density of water at 20°C (in kg/m^3)
	// 	double mu = 0.01; // dynamic viscosity of water at 20°C (in N.s/m^2)
	// 	double gamma = 7.0; // typical value for liquid (dimensionless)
	// 	double c_0 = 5.0;//5.0;//1481; // sound speed in water at 20°C (in m/s)
	// 	double sigma = 1.0; // surface tension of water-air interface at 20°C (in N/m)

		// Brackbill's paper
		// Parameters of the problem
	double l = 2e-2; // radius of the circle
	double L = 2.0*l; // size of the domain: [-L,L] x [-L,L]
	double dt = 0.001; // physical time step
	double T = 1.0; // duration of simulation

	// Physical parameters
	double rho_0 = 1000.0; // initial (physical) density of water at 20°C (in kg/m^3)
	double mu = 0.001;//0.01; // dynamic viscosity of water at 20°C (in N.s/m^2)
	double gamma = 7.0; // typical value for liquid (dimensionless)
	double c_0 = 10.0;//5.0;//1481; // sound speed in water at 20°C (in m/s)
	double sigma = 23.61e-3; // surface tension of water-air interface at 20°C (in N/m)
	double background_pressure = 0.0;
	xy* gravity = xy_new(0.0, 0.0);

	// SPH parameters
	Verlet *verlet = NULL; // don't use Verlet (for now)
	Kernel kernel = Cubic; // kernel choice
	double interface_threshold = 0.9; // If ||n_i|| > threshold => particle i belongs to interface (first detection approach)
	double XSPH_epsilon = 1.0;
	Free_surface_detection surface_detection = DIVERGENCE;
	int N_c = 11;//21; // number of circonferences on which points are placed
	int N_p = 6; // number of points on the first circonference (doubled for every circonference)
	int N_tot = 1; // total number of points
	for (int i = 1; i < N_c; i++) {
		N_tot += i * N_p;
	}
	printf("N_tot = %d \n", N_tot);
	int n_iter = (int)(T / dt); // number of iterations to perform
	double kh = 1.0*l; // kernel width

	// Animation parameter
	double T_anim = 10; // duration of animation
	double dt_anim = T_anim / n_iter; // time step of animation


	// Initialize particles in a circle
	double m = rho_0 * M_PI * l * l / N_tot; // mass of each particle

	Particle** particles = (Particle**)malloc(N_tot * sizeof(Particle*));
	Particle_derivatives** particles_derivatives = malloc(N_tot * sizeof(Particle_derivatives*));
	Residual** residuals = malloc(N_tot * sizeof(Residual*));

	// parameters defining the circle
	double b, delta_s, k, theta;
	delta_s = l / ((double)N_c - 1.0);
	theta = (2 * M_PI) / ((double)N_p);
	double ecc = 0.6;

	int index = 0;
	for (int i = 0; i < N_c; i++) {
		b = i;
		if (b == 0) {
			xy *pos = xy_new(0.0, 0.0);
			xy *v = xy_new(0.0, 0.0);
			particles[index] = Particle_new(index, m, pos, v, rho_0, mu, c_0, gamma, sigma, background_pressure, gravity);
			particles_derivatives[index] = Particle_derivatives_new(index);
			residuals[index] = Residual_new();
			index++;
		}
		else {
			for (int j = 0; j < i*N_p; j++) {
				k = (double)j / b;
				xy *pos = xy_new(b*delta_s*cos(k*theta), b*delta_s*sin(k*theta));
				if (ecc != 0.0) {
					pos->x *= sqrt(2.0 / sin(ecc*M_PI)) * sin(0.5*ecc*M_PI);
					pos->y *= sqrt(2.0 / sin(ecc*M_PI)) * cos(0.5*ecc*M_PI);
				}
				xy *v = xy_new(0.0, 0.0);
				particles[index] = Particle_new(index, m, pos, v, rho_0, mu, c_0, gamma, sigma, background_pressure, gravity);
				particles_derivatives[index] = Particle_derivatives_new(index);
				residuals[index] = Residual_new();
				// 		if (b == N_c - 1) particles[index]->on_free_surface = true; // WARNING: to be removed
				index++;
			}
		}
	}

	//Estimate maximum admissible time step for stability
	double h_p = l / sqrt(N_tot);
	double safety_param = 0.8;
	dt = compute_admissible_dt(safety_param, h_p, c_0, rho_0, mu, sigma, gravity);
	n_iter = (int)(T / dt);

	// Setup grid
	Grid *grid = Grid_new(-L, L, -L, L, kh);
	// Setup animation
	Animation *animation = Animation_new(N_tot, dt_anim, grid, 1);
	// Setup setup
	Setup *setup = Setup_new(n_iter, dt, kh, verlet, kernel, surface_detection, interface_threshold, XSPH_epsilon);


	simulate(grid, particles, particles_derivatives, residuals, N_tot, update_positions_seminar_5, setup, animation);
	// 	simulate(grid, particles, particles_derivatives, residuals, N_tot, update_positions_test_static_bubble, setup, animation);

		// Free stuff
	free_particles(particles, N_tot);
	free_particles_derivatives(particles_derivatives, N_tot);
	free_Residuals(residuals, N_tot);
	Grid_free(grid);
	Setup_free(setup);
	Animation_free(animation);

}

*/
int nb_part_on_bound(int nb_part_first_row, int nb_rows_per_bound) {
	return nb_rows_per_bound * (nb_part_first_row + nb_rows_per_bound - 1);
}
void free_tab(xy** tab, int n)
{
	for (int i = 0;i < n;i++)
		free(tab[i]);
	free(tab);
}

void script_lid_driven_cavity() {

	// Test case taken from Adami et al., "A transport-velocity formulation for smoothed particle hydrodynamics" (2013)

	// Parameters of the problem
	double L_y = 1.0;
	double L_x = 1.0;
	double dt = 0.0; // physical time step --> set below in "compute_admissible_dt"
	double T = 1.0; // duration of simulation

	// Physical parameters
	double speed_upper_boundary = 1.0;
	double Re = 400.0;
	double rho_0 = 1.0;
	double mu = (rho_0*speed_upper_boundary*L_x) / Re;
	double gamma = 7.0; // typical value for liquid (dimensionless)
	double c_0 = 100.0*speed_upper_boundary;
	double background_pressure = 0.0;
	xy* gravity = xy_new(0.0, 0.0);

	// SPH parameters
	Search *search = Search_new(false,0, 0, 0); 
	Kernel kernel = Cubic; // kernel choice
	double XSPH_epsilon = 1.0;
	Free_surface_detection surface_detection = NONE;

	// -------- NUMBER OF PARTICLES AND SIZE --------

	// 1/ *** In the domain ***
	int N_x = 25;
	int N_y = 25;
	double h_x = L_x / ((double)N_x - 1.0);
	double h_y = L_y / ((double)N_y - 1.0);
	int N_x_corrected = N_x - 2; // remove particles on horizontal boundaries
	int N_y_corrected = N_y - 2; // remove particles on vertical boundaries
	int total_nb_part_in_domain = N_x_corrected * N_y_corrected;

	double mass = rho_0 * h_x * h_y;
	double kh = sqrt(21) * 2 * L_x / N_x;


	// 2/ *** On the boundaries ****
	int nb_boundaries = 4;
	int nb_rows_per_bound = 5; // WARNING: we should add a number of rows such that the kernel of a particle next to a boundary will not be truncated

	// Coordinates of the corners
	xy* pos_start_left = xy_new(0.0, 0.0);
	xy* pos_start_upper = xy_new(0.0, L_y);
	xy* pos_start_right = xy_new(L_x, 0.0);
	xy* pos_start_lower = xy_new(0.0, 0.0);
	xy** coord_corners = (xy**)malloc(nb_boundaries * sizeof(xy*));
	coord_corners[0] = pos_start_left, coord_corners[1] = pos_start_upper;
	coord_corners[2] = pos_start_right, coord_corners[3] = pos_start_lower;
	
	// Solid wall B.C. to be imposed on the boundaries // WARNING: only solid wall B.C. are possible for the moment. Might be good to have periodic and inflow/outflow B.C.
	xy** vel_BC = (xy**)malloc(nb_boundaries * sizeof(xy*));
	vel_BC[0] = xy_new(0.0, 0.0); // Boundary 1
	vel_BC[1] = xy_new(speed_upper_boundary, 0.0); // Boundary 2
	vel_BC[2] = xy_new(0.0, 0.0); // Boundary 3
	vel_BC[3] = xy_new(0.0, 0.0); // Boundary 4
	xy** acc_BC = (xy**)malloc(nb_boundaries * sizeof(xy*));
	acc_BC[0] = xy_new(0.0, 0.0); // Boundary 1
	acc_BC[1] = xy_new(0.0, 0.0); // Boundary 2
	acc_BC[2] = xy_new(0.0, 0.0); // Boundary 3
	acc_BC[3] = xy_new(0.0, 0.0); // Boundary 4


	// 3/ *** Combined ***
	int nb_part_per_bound_x = nb_part_on_bound(N_x, nb_rows_per_bound); //for horizontal boundary
	int nb_part_per_bound_y = nb_part_on_bound(N_y, nb_rows_per_bound); //for vertical boundary
	int total_nb_part_on_bound = 2 * (nb_part_per_bound_x+ nb_part_per_bound_y);
	int N_tot = total_nb_part_in_domain + total_nb_part_on_bound;

	Particle** particles = (Particle**)malloc(N_tot * sizeof(Particle*));
	Particle_derivatives** particles_derivatives = malloc(total_nb_part_in_domain * sizeof(Particle_derivatives*));
	Residual** residuals = malloc(total_nb_part_in_domain * sizeof(Residual*));

	// -------- PROPERTIES OF PARTICLES 
	// 1/ *** In the domain ***-
	fill_interior_rectangular(N_x_corrected, N_y_corrected, h_x, h_y, h_x, h_y, particles, particles_derivatives, residuals,
		mass, rho_0, mu, c_0, gamma, background_pressure, gravity);
	// 2/ *** On the boundaries ****
	Boundary* boundaries = build_rectangular_boundaries(coord_corners, particles, total_nb_part_in_domain, nb_part_per_bound_x, nb_part_per_bound_y,N_x,N_y, nb_rows_per_bound, vel_BC, acc_BC, h_x, h_y,
		mass, rho_0, mu, c_0, gamma, background_pressure, gravity);
	free(gravity);free(vel_BC);free(acc_BC);free_tab(coord_corners, nb_boundaries);

	//Estimate maximum admissible time step for stability
	double h_p = h_x;
	double safety_param = 1.0;
	dt = compute_admissible_dt(safety_param, h_p, c_0, rho_0, mu, 0.0, gravity);
	int n_iter = (int)(T / dt);


	// Animation parameter
	double T_anim = 1.0; // duration of animation
	double dt_anim = T_anim / n_iter; // time step of animation

	// Setup grid
	Grid *grid = Grid_new(-L_x, 2.0*L_x, -L_y, 2.0*L_y, kh);
	// Setup animation
	Animation *animation = Animation_new(N_tot, dt_anim, grid, 1);
	// Setup setup
	Setup *setup = Setup_new(n_iter, dt, kh, search, kernel, surface_detection, 0.0, XSPH_epsilon);

	printf(">>>>> dt_admissible = %2.6f <<<<<< \n", dt);
	printf(">>>>> kh = %2.6f <<<<<< \n", kh);
	printf(">>>>> N_part_domain: %d, N_part_boundaries: %d, N_part_total: %d <<<<<< \n", total_nb_part_in_domain, total_nb_part_on_bound, N_tot);
	
	omp_set_num_threads(4);//set the number of threads
	simulate_with_boundaries(grid, particles, particles_derivatives, residuals, N_tot, update_positions_project, setup, animation, total_nb_part_in_domain, boundaries, nb_boundaries);

	// Free stuff
	free_particles(particles, N_tot);
	free_particles_derivatives(particles_derivatives, total_nb_part_in_domain);
	free_Residuals(residuals, total_nb_part_in_domain);
	Grid_free(grid);
	Search_free(search);
	Setup_free(setup);
	Animation_free(animation);
	free_boundaries(boundaries, nb_boundaries);
}


// Without Verlet
void script1() {
	int n_p = 1000;
	Particle** particles = build_particles(n_p, 100);
	double kh = 30;
	double timestep = 1;
	Search* search = Search_new(false,kh,0,0);
	Grid* grid = Grid_new(-100, 100, -100, 100, kh);
	Animation* animation = Animation_new(n_p, 0.05, grid, 1);
	Kernel kernel = Cubic;
	Setup *setup = Setup_new(50, timestep, kh, search, kernel, NONE, 0, 0);

	omp_set_num_threads(4);//set the number of threads
	simulate(grid, particles, NULL, NULL, n_p, random_moves, setup, animation);
	free_particles(particles, n_p);
	Search_free(search);
	Setup_free(setup);
	Grid_free(grid);
	Animation_free(animation);
}

// With Verlet
void script2() {
	int n_p = 2000;
	Particle** particles = build_particles(n_p, 100);
	double kh = 20;
	double vmax = 2;
	int T = 4;
	double timestep = 1;
	double L = 2 * T*vmax*timestep;
	Search* search = Search_new(true,kh,L,T);

	Grid* grid = Grid_new(-100, 100, -100, 100, kh + L);
	Animation* animation = Animation_new(n_p, 0.2, grid, 1);
	Kernel kernel = Cubic;
	Setup *setup = Setup_new(50, timestep, kh, search, kernel, NONE, 0, 0);

	omp_set_num_threads(1);//set the number of threads
	simulate(grid, particles, NULL, NULL, n_p, random_moves, setup, animation);
	free_particles(particles, n_p);
	Search_free(search);
	Setup_free(setup);
	Grid_free(grid);
	Animation_free(animation);
}