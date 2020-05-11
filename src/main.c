#include "print_particules.h"
#include "particle.h"
#include "SPH.h"
#include "derivatives.h"
#include "kernel.h"

#include <math.h>
#include <omp.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
// #include "crtdbg.h" // for memory leak detection; comment if you're on Linux

void script_lid_driven_cavity();

int main() {
	//_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF); // comment if on Linux
	script_lid_driven_cavity();
	return EXIT_SUCCESS;
}

void free_tab(xy** tab, int n) {
	for (int i = 0;i < n;i++)
		free(tab[i]);
	free(tab);
}

void script_lid_driven_cavity() {

	// Parameters of the problem
	double L_y = 1.0;
	double L_x = 1.0;
	double T = 120; // duration of simulation

	// Physical parameters
	double U = 1;
	double Re = 1000;
	double rho_0 = 1.0; // reference density
	double nu = U * L_x / Re; // kinematic viscosity
	double c_0 = 10 * U; // speed of sound
	double gamma = 1;
	double sigma = 0; // surface tension
	double p_0 = gamma * squared(c_0) / rho_0; // reference pressure
	double p_b = p_0; // "We use a background pressure which is on the order of the reference pressure"
	double chi = 0.0; // see Adami2013 (background pressure in the equation of state)
	xy* gravity = xy_new(0.0, 0.0);
	double eps = 1e-6; // TODO: no idea what I should put here


	// SPH parameters
	Search *search = Search_new(false, 0, 0, 0);
	Kernel kernel = Cubic; // kernel choice
	double XSPH_epsilon = 1.0; // XSPH parameter
	Free_surface_detection surface_detection = NONE;

	// -------- NUMBER OF PARTICLES AND SIZE --------

	int N_x_fluid = 50, N_y_fluid = 50;
	double h_x = L_x / (N_x_fluid + 1), h_y = L_y / (N_y_fluid + 1);
	int dN = 5; // number of ghost particle rows (>= 1)
	int N_x = N_x_fluid + 2 * dN, N_y = N_y_fluid + 2 * dN;
	int N_tot = N_x * N_y; // total number of particles

	double m = rho_0 * h_x * h_y; // particle mass
	double kh = sqrt(21 * h_x * h_y);

	Particle** particles = (Particle**)malloc(N_tot * sizeof(Particle*));
	int k = 0;
	double x0 = -h_x * (dN - 1), y0 = -h_y * (dN - 1); // coordinate of most bottom-left particle
	for (int i = 0; i < N_x; i++) {
		for (int j = 0; j < N_y; j++) {
			xy *pos = xy_new(x0 + i * h_x, y0 + j * h_y);
			bool on_boundary = !(pos->x > 0 && pos->x < L_x && pos->y > 0 && pos->y < L_y);
			xy *v = xy_new(0, 0); // initial velocity
			xy *v_imp = NULL;
			if (on_boundary)
				v_imp = xy_new(U * (pos->y <= 0), 0); // 1 moving lid
                // v_imp = xy_new(U * (pos->y >= L_y || pos->y <= 0), 0); // 2 moving lids
			xy* gravity_p = xy_new(gravity->x, gravity->y);
			particles[k] = Particle_new(k, m, pos, v, v_imp, rho_0, nu, c_0, gamma, 0.0, chi, gravity_p, on_boundary);
			k++;
		}
	}

	//Estimate maximum admissible time step for stability
	double h_p = h_x;
	double safety_param = 1.0;
	double dt = compute_admissible_dt(safety_param, h_p, c_0, rho_0, nu, 0.0, gravity, U);
	free(gravity);
	int n_iter = (int)(T / dt);

	// printf("dt_CFL = %lf, dt_visc = %lf\n", 0.25*h_p/(c_0+U), 0.25*h_p*h_p / nu);
	//
	// printf("%lf\n", dt);
	// exit(0);

	// Animation parameter
	double T_anim = 1.0; // duration of animation
	double dt_anim = T_anim / n_iter; // time step of animation

	// Setup grid
	double x1 = -(dN - 1) * h_x;
	double x2 = L_x + (dN - 1) * h_x;
	double y1 = -(dN - 1) * h_y;
	double y2 = L_y + (dN - 1) * h_y;
	Grid *grid = Grid_new(x1, x2, y1, y2, kh);
	// Setup animation
	Animation *animation = Animation_new(N_tot, dt_anim, NULL, h_x/3); // grid == NULL to not plot grid, and only plot bulk particles
	// Setup setup
	Setup *setup = Setup_new(n_iter, dt, kh, search, kernel, surface_detection, 0.0, XSPH_epsilon);

	printf(">>>>> dt_admissible = %2.6f <<<<<< \n", dt);
	printf(">>>>> kh = %2.6f <<<<<< \n", kh);
	printf(">>>>> N_part_domain: %d, N_part_boundaries: %d, N_part_total: %d <<<<<< \n", N_x_fluid * N_y_fluid, N_tot - N_x_fluid * N_y_fluid, N_tot);

	char folder_name[100];
	sprintf(folder_name, "../data/N=%d,Re=%d", N_x_fluid, (int)Re);
	mkdir(folder_name, 0777);

	omp_set_num_threads(2);//set the number of threads
	simulate_with_boundaries(grid, particles, N_tot, setup, animation, eps, folder_name);

	// Free stuff
	free_particles(particles, N_tot);
	Grid_free(grid);
	Search_free(search);
	Setup_free(setup);
	Animation_free(animation);
}
