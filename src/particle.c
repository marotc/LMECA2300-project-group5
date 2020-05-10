#include "particle.h"

// Private functions
bool check_distance(xy *p, xy *q, double r);

//////////////////////PARTICLES//////////////////////
Particle* Particle_new(int index, double m, xy* pos, xy* v, double rho_0, double mu, double c_0, double gamma, double sigma, double background_p, xy* gravity, bool on_boundary) {
	Particle *particle = malloc(sizeof(Particle));
	particle->index = index;
	particle->m = m;
	particle->pos = pos;
	particle->rho = rho_0;
	particle->v = v;
	particle->on_boundary = on_boundary;
	particle->P = 0.0; // assuming that the fluid is at rest (P is the dynamic pressure and not the absolute one!)

	particle->normal = xy_new(0.0, 0.0);
	particle->XSPH_correction = xy_new(0.0, 0.0);
	particle->on_free_surface = false;

	particle->param = malloc(sizeof(Physical_parameters));
	particle->param->rho_0 = rho_0;
	particle->param->dynamic_viscosity = mu;
	particle->param->gamma = gamma;
	particle->param->sound_speed = c_0;
	particle->param->sigma = sigma;
	particle->param->background_p = background_p;
	particle->param->gravity = gravity;

	particle->cell = NULL;
	particle->neighborhood = List_new();
	particle->potential_neighborhood = List_new();
	return particle;
}

void Particle_free(Particle* particle) {
	free(particle->pos);
	free(particle->v);
	free(particle->normal);
	free(particle->param->gravity);
	free(particle->param);
	free(particle->XSPH_correction);
	List_free(particle->neighborhood, NULL);
	List_free(particle->potential_neighborhood, NULL);
	free(particle);
}

void free_particles(Particle** particles, int N) {
	for (int i = 0; i < N; i++)
		Particle_free(particles[i]);
	free(particles);
}

Particle_derivatives* Particle_derivatives_new(int index) {
	Particle_derivatives *particle_derivatives = malloc(sizeof(Particle_derivatives));
	particle_derivatives->index = index;
	particle_derivatives->div_v = 0;
	particle_derivatives->grad_P = xy_new(0, 0);
	particle_derivatives->lapl_v = xy_new(0, 0);
	particle_derivatives->grad_Cs = xy_new(0, 0);
	particle_derivatives->lapl_Cs = 0;
	return particle_derivatives;
}

void Particle_derivatives_free(Particle_derivatives* particle_derivatives) {
	free(particle_derivatives->grad_P);
	free(particle_derivatives->lapl_v);
	free(particle_derivatives->grad_Cs);
	free(particle_derivatives);
}
void Particle_derivatives_reset(Particle_derivatives *particle_derivatives) {
	particle_derivatives->div_v = 0;
	xy_reset(particle_derivatives->grad_P);
	xy_reset(particle_derivatives->lapl_v);
	xy_reset(particle_derivatives->grad_Cs);
	particle_derivatives->lapl_Cs = 0;
}

void free_particles_derivatives(Particle_derivatives** particles_derivatives, int N) {
	for (int i = 0;i < N;i++)
		Particle_derivatives_free(particles_derivatives[i]);
	free(particles_derivatives);
}

double Particle_get_P(Particle *particle) { return particle->P; }
xy * Particle_get_v(Particle *particle) { return particle->v; }
xy * Particle_get_pos(Particle *particle) { return particle->pos; }
double Particle_get_v_x(Particle *particle) { return particle->v->x; }
double Particle_get_v_y(Particle *particle) { return particle->v->y; }
double Particle_get_Cs(Particle *particle) { return particle->Cs; }
xy * Particle_get_normal(Particle *particle) { return particle->normal; }

////////////////////////SEARCH///////////////////////
Grid* Grid_new(double x1, double x2, double y1, double y2, double h) {
	// Build the grid
	int nCellx = ceil((x2-x1) / h);
	int nCelly = ceil((y2-y1) / h);
	printf("Grid size: (%d,%d)\n", nCellx, nCelly);
	Cell*** cells = (Cell***) malloc(nCellx * sizeof(Cell**));

	for(int i = 0; i < nCellx; i++) {
		cells[i] = (Cell**) malloc(nCelly * sizeof(Cell*));
		for(int j = 0; j < nCelly; j++)
			cells[i][j] = (Cell*) malloc(sizeof(Cell));
	}

	// Build links between cells
	for(int i = 0; i < nCellx; i++) {
		for(int j = 0; j < nCelly; j++) {
			cells[i][j]->i = i, cells[i][j]->j = j;
			cells[i][j]->particles = List_new();
			cells[i][j]->visited = false;
			// Assign neighbor cells
			cells[i][j]->neighboring_cells = List_new();
			for (int di = -1; di <= 1; di++) for (int dj = -1; dj <= 1; dj++)
				if ((di || dj) && i+di >= 0 && i+di < nCellx && j+dj >= 0 && j+dj < nCelly)
					List_append(cells[i][j]->neighboring_cells, cells[i+di][j+dj]);
		}
	}
	Grid* grid = (Grid*)malloc(sizeof(Grid));
	grid->left = x1,	grid->right = x1 + nCellx*h; // not very elegant but ok...
	grid->bottom = y1,	grid->top = y1 + nCelly*h;
	grid->nCellx = nCellx;
	grid->nCelly = nCelly;
	grid->h = h;
	grid->cells = cells;
	return grid;
}

void Grid_free(Grid* grid) {
	for(int i = 0; i < grid->nCellx; i++) {
		for(int j = 0; j < grid->nCelly; j++)
			Cell_free(grid->cells[i][j]);
		free(grid->cells[i]);
	}
	free(grid->cells);
	free(grid);
}

void Cell_free(Cell* cell) {
	List_free(cell->neighboring_cells, NULL);
	List_free(cell->particles, NULL);
	free(cell);
}

Search* Search_new(bool verlet, double kh, double L, int T) {
	Search *s = (Search*)malloc(sizeof(Search));
	s->verlet = verlet;
	s->kh = kh;
	s->L = L;
	s->T = T;
	return s;
}

void Search_free(Search* search){
	free(search);
}

///////////////////////////update cells///////////////////////////
Cell* localize_particle(Grid *grid, Particle *p) {
	int i = floor((p->pos->x - grid->left) / grid->h);
	int j = floor((p->pos->y - grid->bottom) / grid->h);
	if(i < 0 || i >= grid->nCellx || j < 0 || j >= grid->nCelly) {
		fprintf(stderr, "ERROR: Particle is outside the grid :(\n");
		exit(0);
	}
	return grid->cells[i][j];
}

// Update links between cells and particles
void update_cells(Grid* grid, Particle** particles, int N) {
	// Clean the grid before update
	reset_grid(grid);
	for(int i = 0; i < N; i++){
		Cell* cell = localize_particle(grid, particles[i]);
		particles[i]->cell = cell; // link cell to particle
		List_append(cell->particles, particles[i]); // link particle to cell
	}
}

/////////////////////////update neighborhood////////////////////////
// Empty the list of particles inside each Cell
void reset_grid(Grid* grid) {
	for (int i = 0; i < grid->nCellx; i++) {
		for (int j = 0; j < grid->nCelly; j++) {
			Cell* cell = grid->cells[i][j];
			List_free(cell->particles, NULL);

			cell->particles = List_new();
			cell->visited = false;
		}
	}
}

// Empty neighborhood of each particle
void reset_particles(Particle** particles, int N, int iter, Search* search) {
	#pragma omp parallel
	{
		int start, end;
		split_thread(N, &start, &end);
		for (int i = start; i <= end; i++) {
			Particle* p = particles[i];
			List_free(p->neighborhood, NULL);
			p->neighborhood = List_new();
			// If in Verlet mode, empty potential nbh
			if (search->verlet && iter%search->T == 0) {
				List_free(p->potential_neighborhood, NULL);
				p->potential_neighborhood = List_new();
			}
		}
	}
}

// Add to the neighbors of particle p all particles q in cell s.t. |p-q| <= r
void add_neighbors_from_cell(Particle* p, Cell* cell , double r) {
	// Iterate over particles in cell
	ListNode *node = cell->particles->head;
	while (node != NULL) {
		Particle* q = (Particle*)node->v;
		if(check_distance(p->pos, q->pos, r))
			List_append(p->neighborhood, q);
		node = node->next;
	}
}

// Add to particle p all its neighbors (from 9 cells)
void add_neighbors_from_cells(Grid* grid, Particle* p) {
	add_neighbors_from_cell(p, p->cell, grid->h);
	ListNode *node = p->cell->neighboring_cells->head;
	while (node != NULL) {
		Cell* cell = (Cell*) node->v;
		add_neighbors_from_cell(p, cell, grid->h);
		node = node->next;
	}
}

void update_neighborhoods_particles(Grid* grid, Particle** particles, int N, Search* search) {
	#pragma omp parallel
	{
		int start, end;
		split_thread(N, &start, &end);
		for (int i = start; i <= end; i++) {
			add_neighbors_from_cells(grid, particles[i]);
			if (search->verlet) {
				List* l = particles[i]->potential_neighborhood;
				particles[i]->potential_neighborhood = particles[i]->neighborhood;
				particles[i]->neighborhood = l;
			}
		}
	}
}


// Among potential neighbors, filter the valid ones
void update_from_potential_neighbors(Particle** particles, int N, Search* search) {
	#pragma omp parallel
	{
		int start, end;
		split_thread(N, &start, &end);
		for (int i = start; i <= end; i++) {
			Particle* p = particles[i];
			ListNode *node = p->potential_neighborhood->head;
			while (node != NULL) {
				Particle* q = (Particle*)node->v;
				if (check_distance(p->pos, q->pos, search->kh))
					List_append(p->neighborhood, q);
				node = node->next;
			}
		}
	}
}

void update_neighborhoods(Grid* grid, Particle** particles, int N, int iter, Search* search) {
	// Clean the particles before update
	reset_particles(particles, N, iter, search);
	if(!search->verlet) {
		update_neighborhoods_particles(grid, particles, N, search);
	} else {
		if (iter%search->T == 0) {
			update_neighborhoods_particles(grid, particles, N, search);
		}
		update_from_potential_neighbors(particles, N, search);
	}
}

/////////////////////////BOUNDARIES///////////////////////////
Boundary* Boundary_new(int index_start_boundary, xy** pos, Particle** part, int nb_part_per_bound, xy* vel_BC, xy* acc_BC,
	double m, double rho_0, double mu, double c_0, double gamma, double background_p, xy* gravity) {
	Boundary* boundary = (Boundary*)malloc(sizeof(Boundary));
	boundary->nb_part_on_bound = nb_part_per_bound;
	boundary->part_on_bound = (Particle**)malloc(nb_part_per_bound * sizeof(Particle*));
	boundary->v_imposed = vel_BC;
	boundary->acc_imposed = acc_BC;
	int j = 0;
	for (int i = index_start_boundary; i < index_start_boundary + nb_part_per_bound; i++) {
		xy* gravityP = xy_new(gravity->x, gravity->y);
		part[i] = Particle_new(i, m, pos[j], xy_new(0.0, 0.0), rho_0, mu, c_0, gamma, 0.0, background_p, gravityP, true); // Initialize particles associated to the boundary in the whole set of particles
		boundary->part_on_bound[j] = part[i]; // Make the "part_on_boun" pointer of the boundary point to the particles associated to the boundary in the whole set of particles
		j++;
	}
	return boundary;
}

void Boundary_free(Boundary* boundary)
{
	free(boundary->part_on_bound);
	free(boundary->v_imposed);
	free(boundary->acc_imposed);
	free(boundary);
}

//upper = 1 if upper boundary, -1 if bottom boundary
xy** build_pos_on_horizontal(double x_min, double y_min, double h_x, double h_y, int n, int nb_part_first_row, int n_row, int upper) {
	xy** pos = (xy**)malloc(n * sizeof(xy*));
	int k = 0;
	for (int i = 0; i < n_row; i++) {
		for (int j = 0;j < nb_part_first_row + 2 * i;j++) {
			pos[k] = xy_new(x_min - i * h_x + j * h_x, y_min + i * h_y*upper);
			k++;
		}
	}
	return pos;
}

//upper = 1 if right boundary, -1 if left boundary
xy** build_pos_on_vertical(double x_min, double y_min, double h_x, double h_y, int n, int nb_part_first_row, int n_col, int right) {
	xy** pos = (xy**)malloc(n * sizeof(xy*));
	int k = 0;
	for (int i = 0; i < n_col; i++) {
		for (int j = 0;j < nb_part_first_row + 2 * i;j++) {
			pos[k] = xy_new(x_min + i * h_x*right, y_min - i * h_y + j * h_y);
			k++;
		}
	}
	return pos;
}
//rectangle boundaries
Boundary* build_rectangular_boundaries(xy** coord, Particle** part, int index_start_boundary, int nb_part_per_bound_x, int nb_part_per_bound_y, int nb_part_first_row_x, int nb_part_first_row_y, int nb_rows_per_bound, xy** vel_BC, xy** acc_BC, double h_x, double h_y,
	double m, double rho_0, double mu, double c_0, double gamma, double background_p, xy* gravity) {
	Boundary** boundaries = (Boundary**)malloc(4 * sizeof(Boundary*));

	//left boundary
	xy** pos = build_pos_on_vertical(coord[0]->x, coord[0]->y, h_x, h_y, nb_part_per_bound_y, nb_part_first_row_y, nb_rows_per_bound, -1);
	boundaries[0] = Boundary_new(index_start_boundary, pos, part, nb_part_per_bound_y, vel_BC[0], acc_BC[0], m, rho_0, mu, c_0, gamma, background_p, gravity);
	free(pos);
	//upper boundary
	pos = build_pos_on_horizontal(coord[1]->x, coord[1]->y, h_x, h_y, nb_part_per_bound_x, nb_part_first_row_x, nb_rows_per_bound, 1);
	boundaries[1] = Boundary_new(index_start_boundary + nb_part_per_bound_y, pos, part, nb_part_per_bound_x, vel_BC[1], acc_BC[1], m, rho_0, mu, c_0, gamma, background_p, gravity);
	free(pos);
	//right boundary
	pos = build_pos_on_vertical(coord[2]->x, coord[2]->y, h_x, h_y, nb_part_per_bound_y, nb_part_first_row_y, nb_rows_per_bound, 1);
	boundaries[2] = Boundary_new(index_start_boundary + nb_part_per_bound_y + nb_part_per_bound_x, pos, part, nb_part_per_bound_y, vel_BC[2], acc_BC[2], m, rho_0, mu, c_0, gamma, background_p, gravity);
	free(pos);
	//lower boundary
	pos = build_pos_on_horizontal(coord[3]->x, coord[3]->y, h_x, h_y, nb_part_per_bound_x, nb_part_first_row_x, nb_rows_per_bound, -1);
	boundaries[3] = Boundary_new(index_start_boundary + 2 * nb_part_per_bound_y + nb_part_per_bound_x, pos, part, nb_part_per_bound_x, vel_BC[3], acc_BC[3], m, rho_0, mu, c_0, gamma, background_p, gravity);
	free(pos);

	return boundaries;
}

void free_boundaries(Boundary** boundaries, int nb_boundaries) {
	for (int i = 0; i < nb_boundaries; i++)
		Boundary_free(boundaries[i]);
	free(boundaries);
}

////////////////////////////////////////////////
// Check if |p-q| <= r
bool check_distance(xy *p, xy *q, double r) {
	return squared(p->x - q->x) + squared(p->y - q->y) <= squared(r);
}

// Generate N particles randomly located on [-L,L] x [-L,L]
// Velocity, rho and e are zero.
Particle** build_particles(int N, double L) {
	Particle** particles = (Particle**)malloc(N * sizeof(Particle*));
	for (int i = 0; i < N; i++) {
		double x = rand_interval(-L, L);
		double y = rand_interval(-L, L);
		xy* pos = xy_new(x, y);
		xy* vel = xy_new(0, 0);
		xy* gravity = xy_new(0.0, 0.0);
		particles[i] = Particle_new(i, 0, pos, vel, 0, 0, 0, 0, 0, 0, gravity, false);
	}
	return particles;
}
