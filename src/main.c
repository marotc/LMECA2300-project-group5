#include "print_particules.h"
#include "particle.h"
#include "SPH.h"
#include "derivatives.h"
#include "kernel.h"

#include <math.h>
#include <omp.h>
#include <sys/stat.h>

void script_anim_from_data();
void script_lid_driven_cavity(int N_x_fluid, double Re);

int main() {
	int use_data;
	printf("Use existing data? [1,0] ");
	scanf("%d", &use_data);
	if(use_data)
		script_anim_from_data();
	else {
		printf("Number of particles per dimension (<=50 is recommended)? ");
		int N; scanf("%d", &N);
		printf("Reynolds number (<=1000 is recommended)? ");
		double Re; scanf("%lf", &Re);
		printf("Number of threads? ");
		int n_threads;scanf("%d", &n_threads);
		omp_set_num_threads(n_threads); // set the number of threads
		script_lid_driven_cavity(N, Re);
	}
	return EXIT_SUCCESS;
}

void script_anim_from_data() {
	bool take_screenshot = true;

	int N_x_fluid = 50, N_y_fluid = 50, Re = 1000;
	double h_p = 1.0/(N_x_fluid+1), U = 1, L_x = 1, L_y = 1;
	double c_0 = 10*U, nu = U*L_x / Re;
	double dt = compute_admissible_dt(-1, h_p, c_0, -1, nu, 0.0, NULL, U);
	double dt_save = 0.01; // max time between each save
	int save_period = (int) (dt_save / dt);
	double real_dt_save = save_period * dt; // true time between each save
	double T = 120;

	double h_x = L_x / (N_x_fluid + 1), h_y = L_y / (N_y_fluid + 1);
	int dN = 5; // number of ghost particle rows (>= 1)
	int N_x = N_x_fluid + 2 * dN, N_y = N_y_fluid + 2 * dN;
	int N_tot = N_x * N_y; // total number of particles

	// Create particles
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
			particles[k] = Particle_new(k, -1, pos, v, v_imp, -1, nu, c_0, -1, 0.0, -1, NULL, on_boundary);
			k++;
		}
	}

	// Animation parameter
	double T_anim = T; // duration of animation (realtime :-))
	double dt_anim = real_dt_save; // time step of animation (realtime :-))
	Animation *animation = Animation_new(N_tot, dt_anim, NULL, h_x/3); // grid == NULL to not plot grid, and only plot bulk particles

	char folder_name[100];
	sprintf(folder_name, "../data/N=%d,Re=%d", N_x_fluid, (int)Re);

	for(int i_save = 0; save_period*i_save*dt < T; i_save++) {
		double t = save_period*i_save*dt;
		char file_path[100];
		sprintf(file_path, "%s/t=%lf.txt", folder_name, t);

		FILE * fp = fopen(file_path, "r");
		if (fp == NULL) {
			printf("file does not exist :-(\n");
	        exit(EXIT_FAILURE);
		}
		char * line = NULL;
	    size_t len = 0;
	    ssize_t read;

		int k, on_boundary;
		double pos_x, pos_y, init_pos_x, init_pos_y, v_x, v_y;
		while ((read = getline(&line, &len, fp)) != -1) {
			sscanf(line, "%d %d %lf %lf %lf %lf %lf %lf", &k, &on_boundary, &pos_x, &pos_y, &init_pos_x, &init_pos_y, &v_x, &v_y);
			particles[k]->pos->x = pos_x;
			particles[k]->pos->y = pos_y;
			particles[k]->init_pos->x = init_pos_x;
			particles[k]->init_pos->y = init_pos_y;
			particles[k]->v->x = v_x;
			particles[k]->v->y = v_y;
	    }
	    fclose(fp);
	    if (line)
	        free(line);

		display_particles(particles, animation, false, take_screenshot, t);
	}
	display_particles(particles, animation, true, take_screenshot, t);
}

void script_lid_driven_cavity(int N_x_fluid, double Re) {

	// Parameters of the problem
	double L_y = 1.0;
	double L_x = 1.0;
	double T = 120; // duration of simulation

	// Physical parameters
	double U = 1;
	double rho_0 = 1.0; // reference density
	double nu = U * L_x / Re; // kinematic viscosity
	double c_0 = 10 * U; // speed of sound
	double gamma = 1;
	double sigma = 0; // surface tension
	double p_0 = gamma * squared(c_0) / rho_0; // reference pressure
	double p_b = p_0; // "We use a background pressure which is on the order of the reference pressure"
	double chi = 0.0; // see Adami2013 (background pressure in the equation of state)
	xy* gravity = xy_new(0.0, 0.0);
	double eps = 1e-6;


	// SPH parameters
	Search *search = Search_new(false, 0, 0, 0);
	Kernel kernel = Cubic; // kernel choice
	double XSPH_epsilon = 1.0; // XSPH parameter
	Free_surface_detection surface_detection = NONE;

	// -------- NUMBER OF PARTICLES AND SIZE --------

	int N_y_fluid = N_x_fluid;
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

	simulate_with_boundaries(grid, particles, N_tot, setup, animation, eps, folder_name);

	// Free stuff
	free_particles(particles, N_tot);
	Grid_free(grid);
	Search_free(search);
	Setup_free(setup);
	Animation_free(animation);
}
