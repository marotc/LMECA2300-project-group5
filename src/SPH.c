#include "SPH.h"

Setup* Setup_new(int iter, double timestep,double kh,Verlet* verlet,Kernel kernel, Free_surface_detection free_surface_detection, double interface_threshold,double XSPH_epsilon) {
	Setup* setup = (Setup*)malloc(sizeof(Setup));
	setup->itermax = iter;
	setup->timestep = timestep;
	setup->kh = kh;
	setup->verlet = verlet;
	setup->kernel = kernel;
	setup->free_surface_detection = free_surface_detection;
	setup->interface_threshold = interface_threshold;
	setup->XSPH_epsilon = XSPH_epsilon;
	return setup;
}

void Setup_free(Setup* setup) {
	if(setup->verlet != NULL)
		free(setup->verlet);
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

void simulate(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, update_positions update_positions, Setup* setup, Animation* animation) {
	double current_time = 0.0;
	printf("%d\n", setup->itermax);
	for (int iter = 0; iter < setup->itermax; iter++) {
		printf("----------------------------------------------------- \n");
		printf("iter %d / %d @ t = %lf \n", iter, setup->itermax, current_time);
		update_cells(grid, particles, n_p);
		update_neighborhoods(grid, particles, n_p, iter, setup->verlet);

		if (animation != NULL)
			display_particles(particles, animation, false, iter);
				 
		update_positions(grid, particles, particles_derivatives, residuals, n_p, setup);
		current_time += setup->timestep;
	}
	update_cells(grid, particles, n_p);
	update_neighborhoods(grid, particles, n_p, 0, setup->verlet);
	if (animation != NULL)
		display_particles(particles, animation, true, -1);
}


void simulate_with_boundaries(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, 
			      update_positions_with_boundaries update_positions, Setup* setup, Animation* animation, Boundary** boundaries, int* index_part_in_domain) {
	double current_time = 0.0;

	for (int iter = 0; iter < setup->itermax; iter++) {
		printf("----------------------------------------------------- \n");
		printf("iter %d / %d @ t = %lf \n", iter, setup->itermax, current_time);
		update_cells(grid, particles, n_p);
		update_neighborhoods(grid, particles, n_p, iter, setup->verlet);
			
		if (animation != NULL)
			display_particles(particles, animation, false, iter);

		update_positions(grid, particles, particles_derivatives, residuals, n_p, setup, boundaries, index_part_in_domain);
		current_time += setup->timestep;
	}
	update_cells(grid, particles, n_p);
	update_neighborhoods(grid, particles, n_p, 0, setup->verlet);
	if (animation != NULL)
		display_particles(particles, animation, true, -1);
}


//move randomly each particles
void random_moves(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup) {
	double max_speed = 2;
	for (int i = 0; i < n_p; i++) {
		double angle = rand_interval(0, 2)*M_PI;
		double speed = rand_interval(0, max_speed);
		Particle *p = particles[i];
		p->v->x = speed * cos(angle);
		p->v->y = speed * sin(angle);
		p->pos->x += p->v->x * setup->timestep;
		p->pos->y += p->v->y * setup->timestep;

		double s = 2;
		//bouncing with the wall
		if (p->pos->x < grid->left)
			p->pos->x = grid->left + s;
		if (p->pos->x > grid->right)
			p->pos->x = grid->right - s;
		if (p->pos->y < grid->bottom)
			p->pos->y = grid->bottom + s;
		if (p->pos->y > grid->top)
			p->pos->y = grid->top - s;
	}
}

void update_positions_seminar_5(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup) {

	// Compute Cs, the XSPH correction on the velocity, and the divergence of the positions
	for (int i = 0; i < n_p; i++) {
		compute_Cs(particles[i], setup->kernel, setup->kh);
		if (setup->XSPH_epsilon != 0.0) compute_XSPH_correction(particles[i], setup->kernel, setup->kh,setup->XSPH_epsilon);
	}

	// Compute derivatives and normal
	for (int i = 0; i < n_p; i++) {
		particles_derivatives[i]->div_v = compute_div(particles[i], Particle_get_v, setup->kernel, setup->kh);
		particles_derivatives[i]->lapl_v->x = compute_lapl(particles[i], Particle_get_v_x, setup->kernel, setup->kh);
		particles_derivatives[i]->lapl_v->y = compute_lapl(particles[i], Particle_get_v_y, setup->kernel, setup->kh);
		compute_grad(particles[i], Particle_get_P, setup->kernel, setup->kh, particles_derivatives[i]->grad_P);
		compute_grad(particles[i], Particle_get_Cs, setup->kernel, setup->kh, particles_derivatives[i]->grad_Cs);
		particles_derivatives[i]->lapl_Cs = compute_lapl(particles[i], Particle_get_Cs, setup->kernel, setup->kh);
		// assemble_residual_NS(particles[i], particles_derivatives[i], residuals[i], setup);
		compute_normal(particles[i], particles_derivatives[i]);
	}

	// Assemble residual and compute curvature
	for (int i = 0; i < n_p; i++) {
	    particles[i]->kappa = 2.0*compute_div(particles[i], Particle_get_normal, setup->kernel, setup->kh);
	    assemble_residual_NS(particles[i], particles_derivatives[i], residuals[i], setup);
	}

	// Integrate (obtain new values, i.e. density, velocities, pressure and positions, at time t+1)
	for (int i = 0; i < n_p; i++)
		time_integrate(particles[i], residuals[i], setup->timestep);
}

void update_positions_project(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup, Boundary** boundaries, int* index_part_in_domain) {
	
	// Assign density, velocity and pressure to boundary particles
	int nb_boundaries = boundaries[0]->nb_boundaries;
	for (int i_b = 0; i_b < nb_boundaries; i_b++) { // loop on the bounadaries
	  apply_BC_Adami(boundaries[i_b], setup->kernel, setup->kh);
	}
	
	int ind;
	int n_p_domain = boundaries[0]->nb_part_in_domain;
	
	// Compute derivatives of the bulk particles and the XSPH correction term --> that's where the effects of the B.C. particles will be feeled
	for (int i = 0; i < n_p_domain; i++) {
	        ind = index_part_in_domain[i];
		if (setup->XSPH_epsilon != 0.0) compute_XSPH_correction(particles[ind], setup->kernel, setup->kh,setup->XSPH_epsilon);
		particles_derivatives[ind]->div_v = compute_div(particles[ind], Particle_get_v, setup->kernel, setup->kh);
		particles_derivatives[ind]->lapl_v->x = compute_lapl(particles[ind], Particle_get_v_x, setup->kernel, setup->kh);
		particles_derivatives[ind]->lapl_v->y = compute_lapl(particles[ind], Particle_get_v_y, setup->kernel, setup->kh);
		compute_grad(particles[ind], Particle_get_P, setup->kernel, setup->kh, particles_derivatives[ind]->grad_P);
// 		compute_grad(particles[ind], Particle_get_Cs, setup->kernel, setup->kh, particles_derivatives[ind]->grad_Cs);
// 		particles_derivatives[ind]->lapl_Cs = compute_lapl(particles[ind], Particle_get_Cs, setup->kernel, setup->kh);
	}

	// Assemble residual based on the computed derivatives of the bulk particles
	for (int i = 0; i < n_p_domain; i++) {
	    ind = index_part_in_domain[i];
	    assemble_residual_NS(particles[ind], particles_derivatives[ind], residuals[ind], setup);
	}

	// Integrate to update the values (i.e. density, velocities, pressure and positions, at time t+1) carried by the bulk particles
	for (int i = 0; i < n_p_domain; i++) {
	    ind = index_part_in_domain[i];
	    time_integrate(particles[ind], residuals[ind], setup->timestep);
	    if (particles[ind]->pos->x < 0.0 || particles[ind]->pos->x > 1.0 || particles[ind]->pos->y < 0.0 || particles[ind]->pos->y > 1.0)
	      printf("!!! Part #%d at (x,y) = (%2.6f, %2.6f) outside the domain !!!\n", i, particles[ind]->pos->x, particles[ind]->pos->y);
// 	    if (norm(particles[ind]->v) > 0.0)
// 	      printf("part #%d: pos (%2.3f, %2.3f), vel (%2.8f, %2.8f) \n", ind, particles[ind]->pos->x, particles[ind]->pos->y, particles[ind]->v->x, particles[ind]->v->y);
	}
	
}

void apply_BC_Adami(Boundary* boundary, Kernel kernel, double kh) {
  
      // loop on the particles on this boundary
      for (int i_p = 0; i_p < boundary->nb_part_on_bound; i_p++) {
	xy* v_a_tilde = xy_new(0.0,0.0);
	double p_w = 0.0;
	double gdotxij = 0.0;
	double denom = 0.0;
	Particle *pi = boundary->part_on_bound[i_p];
	ListNode *node = pi->neighborhood->head;
	
	// loop on the neighbours of the boundary particle to interpolate the velocity and pressure from the particles in the domain
	while (node != NULL) {
		Particle *pj = node->v;
		if (!pj->on_free_surface) { // we interpolate only from the particles inside the domain, we don't take the contributions from the neighbouring boundary particles!
		  // Velocity
		  v_a_tilde->x += pj->v->x * eval_kernel(pi->pos, pj->pos, kh, kernel); // Eq. 22 in Adami 2012
		  v_a_tilde->y += pj->v->y * eval_kernel(pi->pos, pj->pos, kh, kernel); // Eq. 22 in Adami 2012
		  // Pressure
		  gdotxij = (pj->param->gravity->x - boundary->acc_imposed->x) * (pi->pos->x - pj->pos->x)
			  + (pj->param->gravity->y - boundary->acc_imposed->y) * (pi->pos->y - pj->pos->y);
		  p_w += pj->P * eval_kernel(pi->pos, pj->pos, kh, kernel) + pj->rho * gdotxij * eval_kernel(pi->pos, pj->pos, kh, kernel); // Eq. 27 in Adami 2012
		  // Normalization
		  denom += eval_kernel(pi->pos, pj->pos, kh, kernel);
		}
		
		node = node->next;
	}
	
	 // Velocity
	 boundary->part_on_bound[i_p]->v->x = 2.0*boundary->v_imposed->x - v_a_tilde->x / denom; // Eq. 23 in Adami 2012
	 boundary->part_on_bound[i_p]->v->y = 2.0*boundary->v_imposed->y - v_a_tilde->y / denom; // Eq. 23 in Adami 2012
	 // Pressure 
	 boundary->part_on_bound[i_p]->P = p_w / denom; // Eq. 27 in Adami 2012
	 // Density
	 double rho_0 = boundary->part_on_bound[i_p]->param->rho_0;
	 double gamma = boundary->part_on_bound[i_p]->param->gamma;
	 double c = boundary->part_on_bound[i_p]->param->sound_speed;
	 double background_p = boundary->part_on_bound[i_p]->param->background_p;
	 double p_0 = squared(c) * rho_0 / gamma;
	 boundary->part_on_bound[i_p]->rho = rho_0 * pow(((p_w-background_p)/p_0)+1.0, 1.0/gamma); // Eq. 28 in Adami 2012
	 
// 	 if (norm(boundary->part_on_bound[i_p]->v) > 0.0)
// 	      printf("Boundary part #%d: pos (%2.3f, %2.3f), vel (%2.8f, %2.8f) \n", i_p, boundary->part_on_bound[i_p]->pos->x, boundary->part_on_bound[i_p]->pos->y, boundary->part_on_bound[i_p]->v->x, boundary->part_on_bound[i_p]->v->y);
      }
}

void compute_Cs(Particle *particle, Kernel kernel, double kh) {
	particle->Cs = 0;
	Particle *pi = particle;
	ListNode *node = pi->neighborhood->head;
	while (node != NULL) {
		Particle *pj = node->v;
		particle->Cs += (pj->m / pj->rho) * eval_kernel(pi->pos, pj->pos, kh, kernel);
		node = node->next;
	}
	// printf("pos = (%lf, %lf), Cs = %lf\n", particle->pos->x, particle->pos->y, particle->Cs);
}

void compute_normal(Particle *particle, Particle_derivatives* particle_derivatives) {
	particle->normal = xy_new(0.0, 0.0);
	xy *n = particle_derivatives->grad_Cs; // surface normal inward
	double norm_n = norm(n); // norm of n
	particle->normal->x = n->x / norm_n;
	particle->normal->y = n->y / norm_n;
}

// Assemble the residual of the (incompressible) Navier-Stokes equations based on the derivatives available
void assemble_residual_NS(Particle* particle, Particle_derivatives* particle_derivatives, Residual* residual,Setup* setup) {
	double mu_i = particle->param->dynamic_viscosity;

	double rho_i = particle->rho;
	double div_vel_i = particle_derivatives->div_v;
	xy* grad_P = particle_derivatives->grad_P;
	xy* lapl_v = particle_derivatives->lapl_v;
	
	// Body forces to be added
	double fs_x = particle->param->gravity->x ; double fs_y = particle->param->gravity->y;
	
	// --- Surface tension effect ---
	if (setup->free_surface_detection != NONE) { 

	  // Compute UNIT normal vector
	  xy *n = particle_derivatives->grad_Cs; // surface normal inward
	  double norm_n = norm(n);
	  n->x /= norm_n, n->y /= norm_n;

	  double lapl_Cs = particle_derivatives->lapl_Cs;
	  // Choose between curvature estimated with Laplacian of colour field or with divergence of normal
	  // 	double kappa = - lapl_Cs / norm_n; // curvature with Laplacian of colour field

	  
	  // Apply surface tension only on particles in the vicinity the interface
	  bool criterion;
	  // Identification based on the norm of the normal
	  if (setup->free_surface_detection == CSF)
		  criterion = norm_n > setup->interface_threshold;
	  // Identification based on the divergence of the position vector
	  else if (setup->free_surface_detection == DIVERGENCE)
		  criterion = compute_div(particle, Particle_get_pos, setup->kernel, setup->kh) <= setup->interface_threshold;
	  else
		  criterion = false;
	  if (criterion) {
		  particle->on_free_surface = true;
  // 	      fs_x = - particle->param->sigma * lapl_Cs * n->x / norm_n;
  // 	      fs_y = - particle->param->sigma * lapl_Cs * n->y / norm_n;
		double kappa = particle->kappa;//compute_curvature(particle, setup, 0.5);
		fs_x += - particle->param->sigma * kappa * n->x ;
		fs_y += - particle->param->sigma * kappa * n->y ;
		    // double kappa = - lapl_Cs / norm_n; // curvature with Laplacian of colour field
  // 		  printf("pos = (%lf, %lf), n = (%lf, %lf), div_r = %lf, kappa = %lf\n", particle->pos->x, particle->pos->y, n->x, n->y, compute_div(particle, Particle_get_pos, setup->kernel, setup->kh), kappa);
		
	  }
	  else
		  particle->on_free_surface = false;

  // 	// Exact values of normal and curvature for a circle centered in (0,0)
  // 	xy* n_exact = xy_new(particle->pos->x, particle->pos->y);
  // 	double norm_n_exact = norm(n_exact);
  // 	double circle_radius = 1e-3;
  // 	double epsilon = 1e-5;
  // 	double kappa_exact = 1.0 / circle_radius;
  // 	// To print quantities on the surface of the circle
  // 	if (pow(particle->pos->x,2) + pow(particle->pos->y,2) <= pow(circle_radius+epsilon,2) &&  pow(particle->pos->x,2) + pow(particle->pos->y,2) >= pow(circle_radius-epsilon,2)) {
  // 	  printf("pos = (%lf, %lf), n_exact = (%lf, %lf), n = (%lf, %lf), ||n|| = %lf, fs = (%lf, %lf), kappa_exact = %2.3f, kappa = %2.6f \n", particle->pos->x, particle->pos->y,-n_exact->x / norm_n_exact, -n_exact->y / norm_n_exact, n->x / norm_n, n->y / norm_n, norm_n, fs_x, fs_y, kappa_exact, kappa);
  // 	}
	  // To print quantities on the surface of the square
  // 	double x_pos = particle->pos->x, y_pos = particle->pos->y;
  // 	if (x_pos == 50.0 || x_pos == -50.0 || y_pos == 50.0 || y_pos == -50.0) {
  // 	  printf("pos = (%lf, %lf), n = (%lf, %lf), ||n|| = %lf, fs = (%lf, %lf), kappa = %2.6f \n",   particle->pos->x, particle->pos->y, n->x / norm_n, n->y / norm_n, norm_n, fs->x, fs->y, kappa);
  // 	}
  // 	printf("Pressure = %lf, fs_norm = %lf \n", particle->P, sqrt(squared(fs_x)+squared(fs_y)));
	  
	}


	residual->mass_eq = -rho_i * div_vel_i;
	residual->momentum_x_eq = (-1.0/rho_i) * grad_P->x + (mu_i/rho_i) * lapl_v->x + fs_x;
	residual->momentum_y_eq = (-1.0/rho_i) * grad_P->y + (mu_i/rho_i) * lapl_v->y + fs_y;

}



// Time integrate the Navier-Stokes equations based on the residual already assembled
void time_integrate(Particle* particle, Residual* residual, double delta_t) {

	// Update position with an Euler explicit scheme
	particle->pos->x += delta_t * particle->v->x - delta_t * particle->XSPH_correction->x;
	particle->pos->y += delta_t * particle->v->y - delta_t * particle->XSPH_correction->y;

	// Update density and velocity with an Euler explicit scheme (TODO: implement more accurate and more stable schemes)
	particle->rho += delta_t * residual->mass_eq;
	particle->v->x += delta_t * residual->momentum_x_eq;
	particle->v->y += delta_t * residual->momentum_y_eq;
	
// 	// Update position with an Euler explicit scheme // WARNING: before (explicit) or after (implicit) the updates of the field variables?
// 	particle->pos->x += delta_t * particle->v->x - delta_t * particle->XSPH_correction->x;
// 	particle->pos->y += delta_t * particle->v->y - delta_t * particle->XSPH_correction->y;

	// Update pressure with Tait's equation of state
	double B = squared(particle->param->sound_speed) * particle->param->rho_0 / particle->param->gamma;
	particle->P = B * (pow(particle->rho / particle->param->rho_0, particle->param->gamma) - 1) + particle->param->background_p;
	

}

// Normal should be available everywhere!
double compute_curvature(Particle *particle, Setup *setup, double epsilon) {
	double num = epsilon * 2 * compute_div(particle, Particle_get_normal, setup->kernel, setup->kh);
	printf("%lf\n", compute_div(particle, Particle_get_normal, setup->kernel, setup->kh));
	double denom = 0;
	Particle *pi = particle;
	ListNode *node = pi->neighborhood->head;
	while (node != NULL) {
		Particle *pj = node->v;
		//printf("%lf %lf\n", pj->normal->x, pj->normal->y);
		xy *grad_W = grad_kernel(pi->pos, pj->pos, setup->kh, setup->kernel);
		denom += sqrt(squared(pi->pos->x - pj->pos->x) + squared(pi->pos->y - pj->pos->y)) *
			(pj->m / pj->rho) * norm(grad_W);
		free(grad_W);
		node = node->next;
	}
	return num / denom;
}

void compute_XSPH_correction(Particle *pi, Kernel kernel, double kh, double epsilon) {
	xy_reset(pi->XSPH_correction);
	ListNode *node = pi->neighborhood->head;
	while (node != NULL) {
		Particle *pj = node->v;
		pi->XSPH_correction->x += (pj->m / pj->rho) * (pi->v->x - pj->v->x) * eval_kernel(pi->pos, pj->pos, kh, kernel);
		pi->XSPH_correction->y += (pj->m / pj->rho) * (pi->v->y - pj->v->y) * eval_kernel(pi->pos, pj->pos, kh, kernel);
		//printf("%lf\n", pj->m);
		node = node->next;
	}
	pi->XSPH_correction->x *= epsilon;
	pi->XSPH_correction->y *= epsilon;
	//printf("pos = (%lf, %lf), Cs = %lf\n", particle->pos->x, particle->pos->y, particle->Cs);
}

void update_positions_ellipse(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup) {

  	// Compute Cs, the XSPH correction on the velocity, and the divergence of the positions
	for (int i = 0; i < n_p; i++) {
		compute_Cs(particles[i], setup->kernel, setup->kh);
		if (setup->XSPH_epsilon != 0.0) compute_XSPH_correction(particles[i], setup->kernel, setup->kh,setup->XSPH_epsilon);
	}
	// Compute derivatives and residuals
	for (int i = 0; i < n_p; i++) {
		particles_derivatives[i]->div_v = compute_div(particles[i], Particle_get_v, setup->kernel, setup->kh);
		particles_derivatives[i]->lapl_v->x = compute_lapl(particles[i], Particle_get_v_x, setup->kernel, setup->kh);
		particles_derivatives[i]->lapl_v->y = compute_lapl(particles[i], Particle_get_v_y, setup->kernel, setup->kh);
		compute_grad(particles[i], Particle_get_P, setup->kernel, setup->kh, particles_derivatives[i]->grad_P);
		compute_grad(particles[i], Particle_get_Cs, setup->kernel, setup->kh, particles_derivatives[i]->grad_Cs);
		particles_derivatives[i]->lapl_Cs = compute_lapl(particles[i], Particle_get_Cs, setup->kernel, setup->kh);
		assemble_residual_NS(particles[i], particles_derivatives[i], residuals[i],setup);
	}
	int index_x_max, index_x_min, index_y_max, index_y_min;
	double pos_x_max = -INFINITY, pos_x_min = INFINITY, pos_y_max = -INFINITY, pos_y_min = INFINITY;
	// Integrate (new values (i.e. density, velocities) at time t+1)
	for (int i = 0; i < n_p; i++) {
		time_integrate(particles[i], residuals[i], setup->timestep);
		if (particles[i]->pos->x > pos_x_max) pos_x_max = particles[i]->pos->x, index_x_max = i;
		if (particles[i]->pos->x < pos_x_min) pos_x_min = particles[i]->pos->x, index_x_min = i;
		if (particles[i]->pos->y > pos_y_max) pos_y_max = particles[i]->pos->y, index_y_max = i;
		if (particles[i]->pos->y < pos_y_min) pos_y_min = particles[i]->pos->y, index_y_min = i;
	}
	// Compute semi-major axis of an ellipse
	double a_ellipse = particles[index_x_max]->pos->x - particles[index_x_min]->pos->x;
	double b_ellipse = particles[index_y_max]->pos->y - particles[index_y_min]->pos->y;
	printf("a = %lf, b = %lf\n", a_ellipse * 0.5, b_ellipse * 0.5);
}

void update_positions_test_static_bubble(Grid* grid, Particle** particles, Particle_derivatives** particles_derivatives, Residual** residuals, int n_p, Setup* setup) {

	// Compute Cs, the XSPH correction on the velocity, and the divergence of the positions
	int index_x_max, index_x_min, index_y_max, index_y_min;
	double pos_x_max = -INFINITY, pos_x_min = INFINITY, pos_y_max = -INFINITY, pos_y_min = INFINITY;
	for (int i = 0; i < n_p; i++) {
		compute_Cs(particles[i], setup->kernel, setup->kh);
		if (setup->XSPH_epsilon != 0.0) compute_XSPH_correction(particles[i], setup->kernel, setup->kh,setup->XSPH_epsilon);
		// Compute radius of circle
		if (particles[i]->pos->x > pos_x_max) pos_x_max = particles[i]->pos->x, index_x_max = i;
		if (particles[i]->pos->x < pos_x_min) pos_x_min = particles[i]->pos->x, index_x_min = i;
		if (particles[i]->pos->y > pos_y_max) pos_y_max = particles[i]->pos->y, index_y_max = i;
		if (particles[i]->pos->y < pos_y_min) pos_y_min = particles[i]->pos->y, index_y_min = i;
	}
	double radius_circle = 0.5*(particles[index_x_max]->pos->x - particles[index_x_min]->pos->x);

	// Compute derivatives and normal
	for (int i = 0; i < n_p; i++) {
		particles_derivatives[i]->div_v = compute_div(particles[i], Particle_get_v, setup->kernel, setup->kh);
		particles_derivatives[i]->lapl_v->x = compute_lapl(particles[i], Particle_get_v_x, setup->kernel, setup->kh);
		particles_derivatives[i]->lapl_v->y = compute_lapl(particles[i], Particle_get_v_y, setup->kernel, setup->kh);
		compute_grad(particles[i], Particle_get_P, setup->kernel, setup->kh, particles_derivatives[i]->grad_P);
		compute_grad(particles[i], Particle_get_Cs, setup->kernel, setup->kh, particles_derivatives[i]->grad_Cs);
		particles_derivatives[i]->lapl_Cs = compute_lapl(particles[i], Particle_get_Cs, setup->kernel, setup->kh);
// 		assemble_residual_NS(particles[i], particles_derivatives[i], residuals[i], setup);
// 		assemble_residual_NS_test(particles[i], particles_derivatives[i], residuals[i], radius_circle, setup);
		compute_normal(particles[i], particles_derivatives[i]);
	}

	// Assemble residual and compute curvature
	for (int i = 0; i < n_p; i++) {
	    particles[i]->kappa = 2.0*compute_div(particles[i], Particle_get_normal, setup->kernel, setup->kh);
	    assemble_residual_NS_test(particles[i], particles_derivatives[i], residuals[i], radius_circle, setup);
	}

	// Integrate (obtain new values, i.e. density, velocities, pressure and positions, at time t+1)
	for (int i = 0; i < n_p; i++)
		time_integrate(particles[i], residuals[i], setup->timestep);
}

// Assemble the residual of the (incompressible) Navier-Stokes equations based on the derivatives available
void assemble_residual_NS_test(Particle* particle, Particle_derivatives* particle_derivatives, Residual* residual, double radius_circle, Setup* setup) {
	double mu_i = particle->param->dynamic_viscosity;

	double rho_i = particle->rho;
	double div_vel_i = particle_derivatives->div_v;
	xy* grad_P = particle_derivatives->grad_P;
	xy* lapl_v = particle_derivatives->lapl_v;


	xy *n = particle_derivatives->grad_Cs; // surface normal
	double norm_n = norm(n); // norm of n
	double lapl_Cs = particle_derivatives->lapl_Cs;
	double kappa = - lapl_Cs / norm_n; // curvature
	double kappa_2 = particle->kappa;

	// Exact values of normal and curvature for a circle centered in (0,0)
	xy* n_exact = xy_new(particle->pos->x, particle->pos->y);
	double norm_n_exact = norm(n_exact);
	double kappa_exact = 1.0 / radius_circle;

	double fs_x = 0; double fs_y = 0;
	// Apply surface tension only on particles in the vicinity the interface
	if (particle->on_free_surface) {
	    fs_x = - particle->param->sigma * kappa_exact * n->x / norm_n;
	    fs_y = - particle->param->sigma * kappa_exact * n->y / norm_n;
// 	  printf("pos = (%lf, %lf), n_exact = (%lf, %lf), n = (%lf, %lf), ||n|| = %lf, fs = (%lf, %lf), kappa_exact = %2.3f, kappa = %2.6f \n", particle->pos->x, particle->pos->y,-n_exact->x / norm_n_exact, -n_exact->y / norm_n_exact, n->x / norm_n, n->y / norm_n, norm_n, fs_x, fs_y, kappa_exact, kappa);
	  printf("kappa_exact = %2.3f, kappa = %2.6f, kappa_div_n = %2.6f \n", kappa_exact, kappa, kappa_2);


	}

	residual->mass_eq = -rho_i * div_vel_i;
	residual->momentum_x_eq = (-1.0/rho_i) * grad_P->x + (mu_i/rho_i) * lapl_v->x + fs_x;
	residual->momentum_y_eq = (-1.0/rho_i) * grad_P->y + (mu_i/rho_i) * lapl_v->y + fs_y;

}

double compute_admissible_dt(double safety_param, double h_p, double c_0, double rho_0, double mu, double sigma, xy* g) {
  // Relations from "Simulation of surface tension in 2D and 3D with smoothed particle hydrodynamics method", Zhang (2010)
  double dt_1 = 0.25 * h_p / c_0; // propagation of sound waves
  double dt_2 = INFINITY;
  if (mu > 0.0) dt_2 = 0.125 * (h_p*h_p) / (mu / rho_0); // viscous diffusion
  double dt_3 = INFINITY;
  if (sigma > 0.0) dt_3 = 0.25 * sqrt((rho_0*pow(h_p,3))/(2*M_PI*sigma)); // surface tension (capillary waves)
  double dt_4 = INFINITY;
  if (g->x > 0.0 || g->y > 0.0) {
    dt_4 = 0.25 * sqrt(h_p/norm(g)); // surface tension (capillary waves)
  }
  

  double dt_min_interm = fmin(dt_1, dt_2);
  dt_min_interm = fmin(dt_min_interm, dt_3);
  return safety_param * fmin(dt_min_interm, dt_4);

}
