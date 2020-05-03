#include "particle.h"

// Private functions
bool check_distance(xy *p, xy *q, double r);

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

Particle* Particle_new(int index, double m, xy* pos, xy* v, double rho_0, double mu, double c_0, double gamma, double sigma, double background_p, xy* gravity) {
	Particle *particle = malloc(sizeof(Particle));
	particle->index = index;
	particle->m = m;
	particle->pos = pos;
	particle->rho = rho_0;
	particle->v = v;
	particle->P = 0.0; // assuming that the fluid is at rest (P is the dynamic pressure and not the absolute one!)

	particle->normal = xy_new(0.0,0.0);
	particle->XSPH_correction = xy_new(0.0,0.0);
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
	particle_derivatives->grad_P = xy_new(0,0);
	particle_derivatives->lapl_v = xy_new(0,0);
	particle_derivatives->grad_Cs = xy_new(0,0);
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

double Particle_get_P(Particle *particle) {	return particle->P; }
xy * Particle_get_v(Particle *particle) { return particle->v; }
xy * Particle_get_pos(Particle *particle) { return particle->pos; }
double Particle_get_v_x(Particle *particle) { return particle->v->x; }
double Particle_get_v_y(Particle *particle) { return particle->v->y; }
double Particle_get_Cs(Particle *particle) { return particle->Cs; }
xy * Particle_get_normal(Particle *particle) { return particle->normal; }


int put_particles_left(double x, xy** coord, int nb_b) {
  for(int i=0; i<nb_b; i++) {
    if (coord[i]->x < x) return 0;  
  }
  return 1;
}

int put_particles_down(double y, xy** coord, int nb_b) {
  for(int i=0; i<nb_b; i++) {
    if (coord[i]->y < y) return 0;  
  }
  return 1;
}

void fill_pos_on_vertical_row(double x_min, double y_min, double h, int nb_part ,int ind, xy** pos) {
  for (int i=0; i<nb_part; i++) {
    pos[ind+i] = xy_new(x_min,y_min+i*h);
  }
}

void fill_pos_on_horizontal_row(double x_min, double y_min, double h, int nb_part ,int ind, xy** pos) {
  for (int i=0; i<nb_part; i++) {
    pos[ind+i] = xy_new(x_min+i*h,y_min);
//     printf("pos = %2.6f \n",pos[ind+i]->x);
  }
}

Boundary* Boundary_new(int index_b, int nb_b, xy** coord, Particle** part, int index_start_boundary, int nb_part_per_bound, int nb_rows_per_bound, double m, double rho_0, 
		       double mu, double c_0, double gamma, double background_p, xy* gravity, xy* vel_BC, xy* acc_BC, int nb_p_domain) {
  Boundary* boundary = malloc(sizeof(Boundary));
  boundary->nb_part_on_bound = nb_part_per_bound;
  boundary->nb_part_in_domain = nb_p_domain;
  boundary->part_on_bound = (Particle**)malloc(nb_part_per_bound * sizeof(Particle*));
  boundary->v_imposed = xy_new(vel_BC->x, vel_BC->y);
  boundary->acc_imposed = xy_new(acc_BC->x, acc_BC->y);
  boundary->nb_boundaries = nb_b;
  
  
  int nb_part_to_add_row = 0;
  for (int j_row = 1; j_row < nb_rows_per_bound; j_row++) nb_part_to_add_row += j_row*2;
  int nb_part_first_row = (int)(nb_part_per_bound-nb_part_to_add_row)/nb_rows_per_bound;
  xy* coord_1 = coord[2*index_b];
  xy* coord_2 = coord[2*index_b+1];
//   printf("index_b = %d\n",index_b);
//   printf("coord_1_x = %2.3f\n",coord_1->x);
//   printf("coord_1_y = %2.3f\n",coord_1->y);
//   printf("coord_2_x = %2.3f\n",coord_2->x);
//   printf("coord_2_y = %2.3f\n",coord_2->y);
  xy** pos = (xy**)malloc(nb_part_per_bound*sizeof(xy*));
  int i=0;
  
  if (coord_1->x == coord_2->x){ // vertical boundary
    double L_y = abs(coord_1->y - coord_2->y);
    double h_y = L_y / ((double)nb_part_first_row - 1.0);
    double y_min = coord_1->y;
    if (y_min > coord_2->y) y_min = coord_2->y;
    if (put_particles_left(coord_1->x, coord, nb_b)) { // put the rows of particles on the left of the domain boundary
      double x_max = coord_1->x;
      for(int r=0;r<nb_rows_per_bound;r++) { // loop on the rows
	fill_pos_on_vertical_row(x_max-r*h_y, y_min-r*h_y, h_y, nb_part_first_row+2*r,i, pos);
	i += nb_part_first_row+2*r;
      }
    }
    else { // put the rows of particles on the right of the domain boundary
      double x_min = coord_1->x;
      for(int r=0;r<nb_rows_per_bound;r++) { // loop on the rows
	fill_pos_on_vertical_row(x_min+r*h_y, y_min-r*h_y, h_y, nb_part_first_row+2*r,i, pos);
	i += nb_part_first_row+2*r;
      }
    }
  }
  else if (coord_1->y == coord_2->y){ // horizontal boundary
    double L_x = abs(coord_1->x - coord_2->x);
    double h_x = L_x / ((double)nb_part_first_row - 1.0);
    double x_min = coord_1->x;
    if (x_min > coord_2->x) x_min = coord_2->x;
    if (put_particles_down(coord_1->y, coord, nb_b)) { // put the rows of particles below the domain boundary
      double y_max = coord_1->y;
      for(int r=0;r<nb_rows_per_bound;r++) { // loop on the rows
	fill_pos_on_horizontal_row(x_min-r*h_x, y_max-r*h_x, h_x, nb_part_first_row+2*r,i, pos);
	i += nb_part_first_row+2*r;
      }
    }
    else { // put the rows of particles above of the domain boundary
      double y_min = coord_1->y;
      for(int r=0;r<nb_rows_per_bound;r++) { // loop on the rows
	fill_pos_on_horizontal_row(x_min-r*h_x, y_min+r*h_x, h_x, nb_part_first_row+2*r,i, pos);
	i += nb_part_first_row+2*r;
      }
    }
  }
  else {
    printf("WARNING: only horizontal or vertical boundaries are supported for the moment");
  }
    
  int j = 0;
  for (int i=index_start_boundary; i<index_start_boundary+nb_part_per_bound; i++) {
//     printf("pos = %2.6f \n",pos[j]->x);
      part[i] = Particle_new(i, m, pos[j], xy_new(0.0, 0.0), rho_0, mu, c_0, gamma, 0.0, background_p, gravity); // Initialize particles associated to the boundary in the whole set of particles
      part[i]->on_free_surface = true;
      boundary->part_on_bound[j] = part[i]; // Make the "part_on_boun" pointer of the boundary point to the particles associated to the boundary in the whole set of particles
      j++;
  }
  

  return boundary;
  
}

void Boundary_free(Boundary* boundary)
{
  free_particles(boundary->part_on_bound,boundary->nb_part_on_bound);
  free(boundary->v_imposed);
  free(boundary->acc_imposed);
}

void free_boundaries(Boundary** boundaries, int N) {
  for (int i = 0; i < N; i++)
		Boundary_free(boundaries[i]);
	free(boundaries);
}


Verlet* Verlet_new(double kh, double L, int T) {
	Verlet *v = malloc(sizeof *v);
	v->kh = kh;
	v->L = L;
	v->T = T;
	return v;
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
// Add to the neighbors of particle p all particles q in cell s.t. |p-q| <= r
void add_neighbors_from_cell(Particle* p, Cell* cell , double r) {
	// Iterate over particles in cell
	ListNode *node = cell->particles->head;
	while (node != NULL) {
		Particle* q = (Particle*)node->v;
		// if((p->index != q->index) && check_distance(p->pos, q->pos, r))
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

// Among potential neighbors, filter the valid ones
void update_from_potential_neighbors(Particle** particles, int N, double r) {
	for (int i = 0; i < N; i++) {
		Particle* p = particles[i];
		ListNode *node = p->potential_neighborhood->head;
		while (node != NULL) {
			Particle* q = (Particle*)node->v;
			if(check_distance(p->pos, q->pos, r))
				List_append(p->neighborhood, q);
			node = node->next;
		}
	}
}


void update_neighborhoods(Grid* grid, Particle** particles, int N, int iter, Verlet* verlet) {
	// Clean the particles before update
	reset_particles(particles, N, iter, verlet);
	if(verlet==NULL) {
		//update_neighborhoods_improved(grid);
		for (int i = 0 ; i < N; i++)
			add_neighbors_from_cells(grid, particles[i]);
	} else {
		if (iter%verlet->T == 0) {
			// update_neighborhoods_improved(grid);
			for (int i = 0; i < N; i++) {
				add_neighbors_from_cells(grid, particles[i]); // TODO: shouldn't verlet->L be used here?
				// Swap nbh and potential nbh
				// TODO: this is very ugly
				List* l = particles[i]->potential_neighborhood;
				particles[i]->potential_neighborhood = particles[i]->neighborhood;
				particles[i]->neighborhood = l;
			}
		}
		update_from_potential_neighbors(particles, N, verlet->kh);
	}
}

//////////////////update neighborhood-IMPROVED algorithm///////////////

// Compute all neighbors inside a cell
void update_neighbors_1(Cell* cell,double r) {
	ListNode *node = cell->particles->head;
	while (node != NULL) {
		Particle* p = (Particle*)node->v;

		ListNode *node2 = node->next;
		while (node2 != NULL) {
			Particle* q = (Particle*)node2->v;
			if(check_distance(p->pos, q->pos, r)) {
				List_append(p->neighborhood, q);
				List_append(q->neighborhood, p);
			}
			node2 = node2->next;
		}
		node = node->next;
	}
}

//compute all neighbors between 2 cells
void neighbors_in_2_cells(Cell* cell1, Cell* cell2, double r) {
	ListNode *node1 = cell1->particles->head;
	while (node1 != NULL) {
		Particle* p1 = (Particle*)node1->v;

		ListNode *node2 = cell2->particles->head;
		while (node2 != NULL) {
			Particle* p2 = (Particle*)node2->v;
			if(check_distance(p1->pos, p2->pos, r)) {
				List_append(p1->neighborhood, p2);
				List_append(p2->neighborhood, p1);
			}
			node2 = node2->next;
		}

		node1 = node1->next;
	}
}

// Compute all neighbors from this cell's particles
//// Compute the neighbors of particules of cell with the neighboring cells
void update_neighbors_2(Cell* cell, double r) {
	ListNode *node = cell->neighboring_cells->head;
	while (node != NULL) {
		Cell* neighborCell = (Cell*)node->v;
		if (!neighborCell->visited) { // if the neighboring cell was not already visited
			neighbors_in_2_cells(cell, neighborCell, r);
		}
		node = node->next;
	}
}

void update_neighborhoods_improved(Grid* grid) {
	Cell*** cells = grid->cells;
	for(int i = 0 ;i < grid->nCellx; i++) {
		for(int j = 0; j < grid->nCelly; j++) {
			update_neighbors_1(cells[i][j], grid->h);
			update_neighbors_2(cells[i][j], grid->h);
			cells[i][j]->visited = true;
		}
	}
}


////////////////////////////////////////////////////

// Empty the list of particles inside each Cell
void reset_grid(Grid* grid) {
	for(int i = 0; i < grid->nCellx; i++) {
		for(int j = 0; j < grid->nCelly; j++) {
			Cell* cell = grid->cells[i][j];
			List_free(cell->particles, NULL);

			cell->particles = List_new();
			cell->visited = false;
		}
	}
}

// Empty neighborhood of each particle
void reset_particles(Particle** particles, int N, int iter, Verlet* verlet) {
	for (int i = 0; i < N; i++) {
		Particle* p = particles[i];
		List_free(p->neighborhood, NULL);
		p->neighborhood = List_new();
		// If in Verlet mode, empty potential nbh
		if(verlet != NULL && iter%verlet->T == 0) {
			List_free(p->potential_neighborhood, NULL);
			p->potential_neighborhood = List_new();
		}
	}
}

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
		xy* gravity = xy_new(0.0,0.0);
		particles[i] = Particle_new(i, 0, pos, vel, 0, 0, 0, 0, 0, 0, gravity);
	}
	return particles;
}
