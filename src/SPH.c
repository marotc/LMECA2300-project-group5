#include "SPH.h"
#include <assert.h>
#include <time.h>

Setup* Setup_new(int iter, double timestep,double kh,Search* search,Kernel kernel, Free_surface_detection free_surface_detection, double interface_threshold,double XSPH_epsilon) {
	Setup* setup = (Setup*)malloc(sizeof(Setup));
	setup->itermax = iter;
	setup->timestep = timestep;
	setup->kh = kh;
	setup->search = search;
	setup->kernel = kernel;
	setup->free_surface_detection = free_surface_detection;
	setup->interface_threshold = interface_threshold;
	setup->XSPH_epsilon = XSPH_epsilon;
	return setup;
}

void Setup_free(Setup* setup) {
	free(setup);
}

Residual* Residual_new(){
	Residual* residual = (Residual*)malloc(sizeof(Residual));
	residual->mass_eq = 0;
	residual->momentum_x_eq = 0;
	residual->momentum_y_eq = 0;
	return residual;
}

void free_Residuals(Residual** residuals, int N) {
	for (int i = 0;i < N;i++)
		free(residuals[i]);
	free(residuals);
}

void simulate_with_boundaries(Grid* grid, Particle** particles, int n_p,
	Setup* setup, Animation* animation, double eps) {

	double current_time = 0.0;
	double dt = setup->timestep;
	double kh = setup->kh;
	Kernel kernel = setup->kernel;

	int dummy;
	double START = omp_get_wtime();
	for (int iter = 0; iter < setup->itermax; iter++) {
		printf("----------------------------------------------------- \n");
		printf("iter %d / %d @ t = %lf \n", iter, setup->itermax, current_time);

		double start = omp_get_wtime();


		// Step 1 (Kick): Compute intermediate momentum velocity (Eq. 14) and shifting velocity (Eq. 15)
		// Step 2 (Drift): Shift the particles to their new positions (Eq. 16)
		drift(particles,n_p,dt);

		// Todo after particles have moved
		double start_nbh = omp_get_wtime();
		update_cells(grid, particles, n_p);
		update_neighborhoods(grid, particles, n_p, iter, setup->search); 
		double time_nbh = omp_get_wtime() - start_nbh;
		printf("%lf\n", time_nbh);

		// Step 3: Compute density (Eq. 5) and pressure (Eq. 10)
		density_pressure_update(particles,n_p,kh,kernel);

		// Assign density, velocity and pressure to boundary particles
		apply_BC_Adami_all(particles, n_p, setup->kernel, setup->kh);

		// Step 4: Compute inter-particle forces f_pres (Eq. 8), f_visc (Eq. 8), and f_bpres (Eq. 13)
		inter_particle_forces_pres_all(particles,n_p,kh,kernel,eps);

		// Step 5 (Kick): Compute full time-step velocity (Eq. 17)
		kick(particles, n_p, dt, kernel, kh);

		if (animation != NULL) display_particles(particles, animation, false, iter);

		double time_total = omp_get_wtime() - start;
		printf("Time for iteration = %lf s, frac nbh = %.2f\n", time_total, time_nbh / time_total);

		// update_positions(grid, particles, particles_derivatives, residuals, n_p, setup, boundaries, index_part_in_domain);
		current_time += setup->timestep;
	}
	update_cells(grid, particles, n_p);
	update_neighborhoods(grid, particles, n_p, 0, setup->search);
	double TOTAL = omp_get_wtime() - START;
	printf("\n TOTAL = %lf", TOTAL);

	if (animation != NULL)
		display_particles(particles, animation, true, -1);
}


void drift(Particle** particles, int N, double dt) {
	#pragma omp parallel
	{
		int start, end;
		split_thread(N, &start, &end);
		for (int i = start; i <= end; i++) {
			Particle *pi = particles[i];
			if (!pi->on_boundary) { // if bulk particle
				// Eq. 14
				pi->v->x += (dt / 2) * pi->a->x;
				pi->v->y += (dt / 2) * pi->a->y;
				// Eq. 15
				pi->vs->x = pi->v->x + (dt / 2) * pi->as->x;
				pi->vs->y = pi->v->y + (dt / 2) * pi->as->y;
				// Eq. 16
				pi->pos->x += dt * pi->vs->x;
				pi->pos->y += dt * pi->vs->y;
			}
		}
	}
}

void density_pressure_update(Particle** particles, int N, double kh, Kernel kernel) {
	#pragma omp parallel
	{
		int start, end;
		split_thread(N, &start, &end);
		for (int i = start; i <= end; i++) {
			Particle *pi = particles[i];
			if (!pi->on_boundary) { // only for fluid particles
				pi->rho = pi->V = 0;
				ListNode *node = pi->neighborhood->head;
				while (node != NULL) {
					Particle *pj = node->v;
					double Wij = eval_kernel(pi->pos, pj->pos, kh, kernel);
					// PySPH's SummationDensity
					pi->rho += pi->m * Wij;
					pi->V += Wij; // inverse of particle volume
					node = node->next;
				}
				// PySPH's StateEquation
				pi->P = squared(pi->param->sound_speed) * (pi->rho - pi->param->rho_0);
				// printf("(%lf, %lf): m = %lf, rho = %lf, V = %lf, p = %lf\n", pi->pos->x, pi->pos->y, pi->m, pi->rho, pi->V, pi->P);
			}
		}
	}
}

void apply_BC_Adami_all(Particle **particles, int N, Kernel kernel, double kh) {
	#pragma omp parallel
	{
		int start, end;
		split_thread(N, &start, &end);
		for (int i = start; i <= end; i++)
			apply_BC_Adami(particles[i], kernel, kh);
	}
}

void apply_BC_Adami(Particle* pi, Kernel kernel, double kh) {

	// loop on the particles on this boundary
	// PySPH's SetWallVelocity
	Particle *pj;
	if(pi->on_boundary) {
		xy *vf = xy_new(0,0); // filtered velocity
		double p_w = 0.0; // future pressure of particle
		double sumW = 0; // normalization factor is different from V as the particles near the boundary do not have full kernel support

		double Wij, gdotxij;

		ListNode *node = pi->neighborhood->head;

		// loop on the neighbours of the boundary particle to interpolate the velocity and pressure from the particles in the domain
		while (node != NULL) {
			pj = node->v;
			if (!pj->on_boundary) { // we interpolate only from the particles inside the domain, we don't take the contributions from the neighbouring boundary particles!
				Wij = eval_kernel(pi->pos, pj->pos, kh, kernel);
				// Normalization
				sumW += eval_kernel(pi->pos, pj->pos, kh, kernel);
				// Velocity (Adami2012, Eq. 22)
				vf->x += pj->v->x * Wij;
				vf->y += pj->v->y * Wij;
				// Pressure
				gdotxij = (pj->param->gravity->x) * (pi->pos->x - pj->pos->x)
				  + (pj->param->gravity->y) * (pi->pos->y - pj->pos->y);
				assert(gdotxij == 0);
				p_w += pj->P * Wij + pj->rho * gdotxij * Wij; // Eq. 27 in Adami 2012
			}
			node = node->next;
		}
		// calculation is done only for the relevant boundary particles.
	    // denom (and v_a_tilde, p_w) is 0 for solid particles sufficiently away from the
	    // solid-fluid interface
		if (sumW > 1e-12) {
			vf->x /= sumW;
			vf->y /= sumW;
			p_w /= sumW;
		}
		// Velocity
		pi->v->x = 2.0*pi->v_imp->x - vf->x; // Eq. 23 in Adami 2012
		pi->v->y = 2.0*pi->v_imp->y - vf->y; // Eq. 23 in Adami 2012
		// Pressure
		pi->P = p_w; // Eq. 27 in Adami 2012
		// Density
		double rho_0 = pi->param->rho_0;
		double c = pi->param->sound_speed;
		double p_0 = squared(c) * rho_0;
		pi->rho = rho_0 * (p_w/p_0 + 1); // Eq. 28 in Adami 2012

		// Cleanup
		free(vf);
	}
}

void inter_particle_forces_pres_all(Particle** particles, int N, double kh, Kernel kernel, double eps) {
	#pragma omp parallel
	{
		int start, end;
		split_thread(N, &start, &end);
		for (int i = start; i <= end; i++)
			inter_particle_forces_pres(particles[i], kh, kernel, eps);
	}
}

void inter_particle_forces_pres(Particle* pi, double kh, Kernel kernel, double eps) {
	double Vi, Vj, Vi2Vj2, r2ij, p_ij, tmp, etai, etaj, eta_ij, F_ij;
	xy *DWij;
	Particle *pj;
	if (!pi->on_boundary) {
		pi->a->x = pi->a->y = 0; // momentum acceleration
		pi->as->x = pi->as->y = 0; // drift acceleration
		ListNode *node = pi->neighborhood->head;
		int cnt = 0;
		while (node != NULL) {
			pj = node->v;
			cnt++;

			// PySPH
			Vi = 1.0 / pi->V;
			Vj = 1.0 / pj->V;
			Vi2Vj2 = Vi * Vi + Vj * Vj;
			r2ij = squared(pi->pos->x - pj->pos->x) + squared(pi->pos->y - pj->pos->y); // squared distance between pi and pj
			DWij = grad_kernel(pi->pos, pj->pos, kh, kernel);

			// PySPH's MomentumEquationPressureGradient
			p_ij = (pj->rho * pi->P + pi->rho * pj->P) / (pi->rho + pj->rho); // density-weighted pressure
			tmp = (-p_ij / pi->m) * Vi2Vj2;

			pi->a->x += tmp * DWij->x;
			pi->a->y += tmp * DWij->y;

			tmp = -pi->param->background_p / pi->m * Vi2Vj2;
			pi->as->x += tmp * DWij->x;
			pi->as->y += tmp * DWij->y;

			// PySPH's MomentumEquationViscosity AND SolidWallNoSlipBC
			etai = pi->param->nu * pi->rho;
			etaj = pj->param->nu * pj->rho;
			eta_ij = 2 * etai * etaj / (etai + etaj);
			F_ij = DWij->x * (pi->pos->x - pj->pos->x) + DWij->y * (pi->pos->y - pj->pos->y);
			tmp = (1.0 / pi->m) * Vi2Vj2 * eta_ij * F_ij / (r2ij + eps);
			// if(r2ij == 0) assert(pi == pj);
			pi->a->x += tmp * (pi->v->x - pj->v->x);
			pi->a->y += tmp * (pi->v->y - pj->v->y);

			// Good luck to proofread this...
			// PySPH's MomentumEquationArtificialStress
			if (!pj->on_boundary) { // wall particles have no drift velocity vs, so they don't contribute!!
				tmp = 1.0 / pi->m * Vi2Vj2;
				pi->a->x += tmp * (
				(pi->rho * pi->v->x * (pi->vs->x - pi->v->x) + pj->rho * pj->v->x * (pj->vs->x - pj->v->x)) / 2 * DWij->x +
				(pi->rho * pi->v->x * (pi->vs->y - pi->v->y) + pj->rho * pj->v->x * (pj->vs->y - pj->v->y)) / 2 * DWij->y
				);
				pi->a->y += tmp * (
				(pi->rho * pi->v->y * (pi->vs->x - pi->v->x) + pj->rho * pj->v->y * (pj->vs->x - pj->v->x)) / 2 * DWij->x +
				(pi->rho * pi->v->y * (pi->vs->y - pi->v->y) + pj->rho * pj->v->y * (pj->vs->y - pj->v->y)) / 2 * DWij->y
				);
			}
			free(DWij);
			node = node->next;
		}
		assert(cnt >= 21);
	}
}

void kick(Particle **particles, int N, double dt, Kernel kernel, double kh) {
	#pragma omp parallel
	{
		int start, end;
		split_thread(N, &start, &end);
		for (int i = start; i <= end; i++) {
			Particle *pi = particles[i];
			if (!pi->on_boundary) {
				pi->v->x += (dt / 2) * pi->a->x;
				pi->v->y += (dt / 2) * pi->a->y;
			}
		}
	}
}
double compute_admissible_dt(double safety_param, double h_p, double c_0, double rho_0, double nu, double sigma, xy* g, double U) {
	return fmin(0.25 * h_p / (c_0 + U), 0.25 * h_p * h_p / nu);
}
