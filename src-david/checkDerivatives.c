
#include "checkDerivatives.h"

singleParticleDerivatives* initialize_particle_derivatives(int size_values)
{
  singleParticleDerivatives* my_deriv = malloc(sizeof(singleParticleDerivatives));

  if (size_values == 1) { // scalar quantity
      my_deriv->divergence = NULL; //  divergence of a scalar has no meaning (= gradient)
      my_deriv->gradient = (double*) calloc(2,sizeof(double));
      my_deriv->laplacian = (double*) calloc(1, sizeof(double));
  }
  else if (size_values == 2) { // vector quantity
      my_deriv->divergence = (double*) calloc(1, sizeof(double));
      my_deriv->gradient = (double*) calloc(4,sizeof(double));
      my_deriv->laplacian = (double*) calloc(2, sizeof(double));
  }
  else {
    printf("---------------  The size of the quantity for which you try to compute the derivatives is not correct ------------------------");
    return 0;
  }

  return my_deriv;
}


mySingleParticle* create_array_of_particles(int nbParticles, int size_values, neighborhood* nh)
{

  mySingleParticle* my_particles = malloc(nbParticles*sizeof(mySingleParticle));


//   double x_lim[2] = {-1.0, 1.0};
//   for (int i = 0; i < nbParticles; i++) {
//     my_particles[i].coordinates = calloc(2,sizeof(double));
//     my_particles[i].values = calloc(size_values,sizeof(double));
//     my_particles[i].mass = 0.0;
//     my_particles[i].density = 0.0;
//     my_particles[i].size_values = size_values;
//     my_particles[i].particle_neighbours = neighborhood_options_init(0.0, 0.0);
//     my_particles[i].particle_neighbours->nh = &(nh[i]);
//     my_particles[i].particle_derivatives = initialize_particle_derivatives(size_values);
//
//     init1DSegmentWithParticles(x_lim, my_particles[i].coordinates, my_particles[i].values, &(my_particles[i].mass), &(my_particles[i].density), nbParticles, i, size_values);
//   }
  double x_lim[4] = {-1.0, 1.0, -1.0, 1.0};
  int nbPart_x = sqrt(nbParticles);
  int nbPart_y = nbPart_x;
  int i = 0;
  for (int ind_y = 0; ind_y < nbPart_y; ind_y++) {
      for (int ind_x = 0; ind_x < nbPart_x; ind_x++) {
	my_particles[i].coordinates = calloc(2,sizeof(double));
	my_particles[i].values = calloc(size_values,sizeof(double));
	my_particles[i].mass = 0.0;
	my_particles[i].density = 0.0;
	my_particles[i].size_values = size_values;
	my_particles[i].particle_neighbours = neighborhood_options_init(0.0, 0.0);
	my_particles[i].particle_neighbours->nh = &(nh[i]);
	my_particles[i].particle_derivatives = initialize_particle_derivatives(size_values);

	initSquareWithParticles(x_lim, my_particles[i].coordinates, my_particles[i].values, &(my_particles[i].mass), &(my_particles[i].density), nbPart_x, nbPart_y, ind_x, ind_y, size_values);
	i++;
      }
   }

//    allParticles* my_particles = malloc(sizeof(allParticles));
  return my_particles;
}


allParticles* create_particles_in_domain(int* nbParticles, double* domain_lim, int size_values, double* args_init_fct, double (*f_init)(double*, double*), int starting_index)
{

  allParticles* my_particles = malloc(sizeof(allParticles));

  int nbPart_x = nbParticles[0];
  int nbPart_y = nbParticles[1];
  my_particles->array_of_particles = malloc(nbPart_x*nbPart_y*sizeof(mySingleParticle));
  int* index_domain = calloc(nbPart_x*nbPart_y, sizeof(int));

  int i = 0;
  for (int ind_y = 0; ind_y < nbPart_y; ind_y++) {
      for (int ind_x = 0; ind_x < nbPart_x; ind_x++) {
	my_particles->array_of_particles[i].coordinates = calloc(2,sizeof(double));
	my_particles->array_of_particles[i].values = calloc(size_values,sizeof(double));
	my_particles->array_of_particles[i].mass = 0.0;
	my_particles->array_of_particles[i].density = 0.0;
	my_particles->array_of_particles[i].size_values = size_values;
	my_particles->array_of_particles[i].particle_neighbours = neighborhood_options_init(0.0, 0.0);
	my_particles->array_of_particles[i].particle_derivatives = initialize_particle_derivatives(size_values);

	initSquareWithParticles(domain_lim, my_particles->array_of_particles[i].coordinates, my_particles->array_of_particles[i].values, &(my_particles->array_of_particles[i].mass), &(my_particles->array_of_particles[i].density),
				nbPart_x, nbPart_y, ind_x, ind_y, size_values);

	initValuesParticles(my_particles->array_of_particles[i].coordinates, my_particles->array_of_particles[i].values, args_init_fct, size_values, f_init);
	index_domain[i] = i+starting_index;
	i++;
      }
  }
  my_particles->nb_particles = i;
  my_particles->nb_particles_domain = i;
  my_particles->nb_particles_boundaries = 0;
  my_particles->index_part_in_domain = index_domain;
  my_particles->index_part_on_boundaries = NULL;

  return my_particles;
}

allParticles* create_particles_on_boundaries(int* nbParticles, double* boundaries_lim, int size_values, int nbBoundaries, double* value_on_boundaries, int starting_index)
{

  allParticles* my_particles = malloc(sizeof(allParticles));

  int nbPart_tot = 0;
  for (int i=0; i<nbBoundaries; i++) nbPart_tot += nbParticles[i];
  my_particles->array_of_particles = malloc(nbPart_tot*sizeof(mySingleParticle));
  int* index_domain = calloc(nbPart_tot, sizeof(int));

  int i = 0;
  for (int ind = 0; ind < nbBoundaries; ind++) {

    double* bound_lim_local = calloc(4, sizeof(double));
    bound_lim_local[0] = *((boundaries_lim+ind*nbBoundaries) + 0);
    bound_lim_local[1] = *((boundaries_lim+ind*nbBoundaries) + 1);
    bound_lim_local[2] = *((boundaries_lim+ind*nbBoundaries) + 2);
    bound_lim_local[3] = *((boundaries_lim+ind*nbBoundaries) + 3);

    double* value_bound_local = calloc(size_values, sizeof(double));
    value_bound_local[0] = value_on_boundaries[size_values*ind];
    if (size_values > 1) value_bound_local[1] = value_on_boundaries[size_values*ind+1];

    for (int j = 0; j < nbParticles[ind]; j++) {
	my_particles->array_of_particles[i].coordinates = calloc(2,sizeof(double));
	my_particles->array_of_particles[i].values = calloc(size_values,sizeof(double));
	my_particles->array_of_particles[i].mass = 0.0;
	my_particles->array_of_particles[i].density = 0.0;
	my_particles->array_of_particles[i].size_values = size_values;
	my_particles->array_of_particles[i].particle_neighbours = neighborhood_options_init(0.0, 0.0);
	my_particles->array_of_particles[i].particle_derivatives = initialize_particle_derivatives(size_values);

	init1DSegmentWithParticles(bound_lim_local, my_particles->array_of_particles[i].coordinates, my_particles->array_of_particles[i].values, &(my_particles->array_of_particles[i].mass), &(my_particles->array_of_particles[i].density),
				nbParticles[ind], j, size_values);

	my_particles->array_of_particles[i].values[0] = value_bound_local[0];
	if (size_values > 1) my_particles->array_of_particles[i].values[1] = value_bound_local[1];

	index_domain[i] = i+starting_index;
	i++;
    }
  }
  my_particles->nb_particles = i;
  my_particles->nb_particles_domain = 0;
  my_particles->nb_particles_boundaries = i;
  my_particles->index_part_in_domain = NULL;
  my_particles->index_part_on_boundaries = index_domain;


  return my_particles;
}

allParticles* create_all_particles(int* nbParticles, double* domain_lim, int size_values, neighborhood* nh, particleType partType, double* args_init_fct, double (*f_init)(double*, double*))
{

  allParticles* my_particles = malloc(sizeof(allParticles));
  if (partType == IN_DOMAIN) {
    int nbPart_x = nbParticles[0];
    int nbPart_y = nbParticles[1];
    my_particles->array_of_particles = malloc(nbPart_x*nbPart_y*sizeof(mySingleParticle));

    int i = 0;
    for (int ind_y = 0; ind_y < nbPart_y; ind_y++) {
	for (int ind_x = 0; ind_x < nbPart_x; ind_x++) {
	  my_particles->array_of_particles[i].coordinates = calloc(2,sizeof(double));
	  my_particles->array_of_particles[i].values = calloc(size_values,sizeof(double));
	  my_particles->array_of_particles[i].mass = 0.0;
	  my_particles->array_of_particles[i].density = 0.0;
	  my_particles->array_of_particles[i].size_values = size_values;
	  my_particles->array_of_particles[i].particle_neighbours = neighborhood_options_init(0.0, 0.0);
	  if (nh) {
	    my_particles->array_of_particles[i].particle_neighbours->nh = &(nh[i]);
	  }
	  my_particles->array_of_particles[i].particle_derivatives = initialize_particle_derivatives(size_values);

	  initSquareWithParticles(domain_lim, my_particles->array_of_particles[i].coordinates, my_particles->array_of_particles[i].values, &(my_particles->array_of_particles[i].mass), &(my_particles->array_of_particles[i].density),
				  nbPart_x, nbPart_y, ind_x, ind_y, size_values);

	  initValuesParticles(my_particles->array_of_particles[i].coordinates, my_particles->array_of_particles[i].values, args_init_fct, size_values, f_init);

	  i++;
	}
    }
    my_particles->nb_particles = i;
 }
 else if (partType == ON_BOUNDARY) {

 }
 else {
   printf("Unknown type of particles. Particles should be in the domain or on the boundaries.");
 }

  return my_particles;
}

void initSquareWithParticles(double* x_lim, double* coord, double* values, double* mass, double* density, int nbPart_x, int nbPart_y, int index_x, int index_y, int size_values)// double (*myFun)(double*))
{
  double x_min = x_lim[0];
  double x_max = x_lim[1];
  double y_min = x_lim[2];
  double y_max = x_lim[3];
  if (x_min > x_max) {
      double x_temp = x_min;
      x_min = x_max;
      x_max = x_temp;
  }
  if (y_min > y_max) {
      double y_temp = y_min;
      y_min = y_max;
      y_max = y_temp;
  }
  double length_x = fabs(x_max-x_min);
  double delta_x = length_x / ((double)nbPart_x - 1.0);
  double length_y = fabs(y_max-y_min);
  double delta_y = length_y / ((double)nbPart_y - 1.0);

  coord[0] = x_min + index_x*delta_x;
  coord[1] = y_min + index_y*delta_y;

//   // function values
//   values[0] = myFunctionToDerive(coord);//(myFun)(&x_coord);
//   if (size_values > 1) values[1] = 0.0; // 1-D segment so don't care about the second dimension
  /**mass = 1.0;
  *density = *mass / (delta_x*delta_y);//WARNING: only valid for cartesian grid of particle!*/
  *density = 1.0;
  *mass = *density * (delta_x*delta_y);
}

void initValuesParticles(double* coord, double *values, double* args, int size_values, double (*f)(double*, double*)) {
  values[0] = (f)(coord, args);//(myFun)(&x_coord);
  if (size_values > 1) values[1] = 0.0;
}

void init1DSegmentWithParticles(double* x_lim, double* coord, double* values, double* mass, double* density, int nb_particles, int index_part, int size_values)// double (*myFun)(double*))
{
  double x_min = x_lim[0];
  double x_max = x_lim[1];
  double y_min = x_lim[2];
  double y_max = x_lim[3];
  if (x_min > x_max) {
      double x_temp = x_min;
      x_min = x_max;
      x_max = x_temp;
  }
  if (y_min > y_max) {
      double y_temp = y_min;
      y_min = y_max;
      y_max = y_temp;
  }
  double delta;
  if (y_min == y_max) { // horizontal segment
    double length = fabs(x_max-x_min);
    delta = length / ((double)nb_particles - 1.0);
    double x_coord = x_min + index_part*delta; // uniformly distributed points on the segment
    // coordinates
    coord[0] = x_coord ;
    coord[1] = y_min;
  //   // function values
  //   values[0] = myFunctionToDerive(coord);//(myFun)(&x_coord);
  //   if (size_values > 1) values[1] = 0.0; // 1-D segment so don't care about the second dimension
  }
  else if (x_min == x_max) { // vertical segment
    double length = fabs(y_max-y_min);
    delta = length / ((double)nb_particles - 1.0);
    double y_coord = y_min + index_part*delta; // uniformly distributed points on the segment
    // coordinates
    coord[0] = x_min ;
    coord[1] = y_coord;
  }
  else {
    printf("Only able to treat horizontal or vertical boundaries for the moment");
  }

  /**mass = 1.0;
  *density = *mass / (delta*delta);*/
  *density = 1.0;
  *mass = *density * (delta*delta);
}

double myFunctionToDerive(double* x, double* args)
{
   return x[0]*x[0];
}

double initFunction(double* x, double* args)
{
   double value_to_be_imposed = args[0];
   return value_to_be_imposed;
//    return x[0]*x[0];
}

double initFunction_GaussianSource(double* x_coord, double* args)
{
   double intensity = args[0];
   double sigma = args[1];
   double x_0 = args[2];
   double y_0 = args[3];
   double x = x_coord[0];
   double y = x_coord[1];
   double arg_exp = (1/(2*sigma*sigma))*((x-x_0)*(x-x_0) + (y-y_0)*(y-y_0));
   return intensity* exp(-arg_exp);
//    return x[0]*x[0];
}

void computeDerivativesOfParticleQuantity(mySingleParticle* myPart, int index_part, kernelType kType) {

  mySingleParticle* local_part = &(myPart[index_part]);
  double * divergence_part = local_part->particle_derivatives->divergence;
  double * grad_part = local_part->particle_derivatives->gradient;
  double * laplacian_part = local_part->particle_derivatives->laplacian;

  // Reset all the derivatives to zero before computing their new values
  // --- Divergence ---
  if (myPart->size_values == 2) divergence_part[0] = 0.0;
  // --- Gradient & Laplacian ---
  for (int k=0; k < myPart->size_values; k++) {
    grad_part[2*k] = 0.0;
    grad_part[2*k+1] = 0.0;
    laplacian_part[k] = 0.0;
  }

  int nNeigh = local_part->particle_neighbours->nh->nNeighbours;
  neighbours* List = local_part->particle_neighbours->nh->list;
  double kh = local_part->particle_neighbours->kh;
  for (int j = 0; j < nNeigh; j++) {
      int index_j = List->index;
      double d_ij = List->distance; // WARNING: distance = 0 in 1-D if computed this way, why???
      double d_x_ij =  myPart[index_j].coordinates[0] - local_part->coordinates[0];// WARNING: absolute value or not here?
//       double d_ij = fabs(d_x_ij); // only valid in 1-D, should use List->distance instead but = 0 for the moment (see line above)
      double d_y_ij =  myPart[index_j].coordinates[1] - local_part->coordinates[1];// WARNING: absolute value or not here?

      double weight_x;
      double weight_y;
      if (kType == CUBIC) {
	weight_x = grad_w_cubic(d_ij, kh, d_x_ij);
	weight_y = grad_w_cubic(d_ij, kh, d_y_ij);
      }
      else if (kType == LUCY) {
	weight_x = grad_w_lucy(d_ij, kh, d_x_ij);
	weight_y = grad_w_lucy(d_ij, kh, d_y_ij);
      }
      else if (kType == QUARTIC) {
	weight_x = grad_w_newquartic(d_ij, kh, d_x_ij);
	weight_y = grad_w_newquartic(d_ij, kh, d_y_ij);
      }
      else if (kType == QUINTIC) {
        weight_x = grad_w_quinticspline(d_ij, kh, d_x_ij);
        weight_y = grad_w_quinticspline(d_ij, kh, d_y_ij);
      }

      // --- Divergence ---
      if (myPart->size_values == 2) divergence_part[0] += (1.0 / local_part->density) * ( (myPart[index_j].values[0] - local_part->values[0]) * weight_x
							  + (myPart[index_j].values[1] - local_part->values[1]) * weight_y ) * myPart[index_j].mass;
      // --- Gradient & Laplacian ---
      for (int k=0; k < myPart->size_values; k++) {
	// d/dx value_k
	grad_part[2*k] += -local_part->density * myPart[index_j].mass * ((local_part->values[k] / (local_part->density*local_part->density))
								      + (myPart[index_j].values[k] / (myPart[index_j].density*myPart[index_j].density))) * weight_x;
	// d/dy value_k
	grad_part[2*k+1] += -local_part->density * myPart[index_j].mass * ((local_part->values[k] / (local_part->density*local_part->density))
								      + (myPart[index_j].values[k] / (myPart[index_j].density*myPart[index_j].density))) * weight_y;
	// d^2/dx^2 value_k + d^2/dy^2 value_k
	laplacian_part[k] += 2.0 * (myPart[index_j].mass / myPart[index_j].density) * (local_part->values[k] - myPart[index_j].values[k]) * (d_x_ij * weight_x + d_y_ij * weight_y) /(d_ij*d_ij);
      }
//       if (index_part == 1225) printf("index_j = %2.6f \n", laplacian_part[0]);
//       double test = 2.0 * (myPart[index_j].mass / myPart[index_j].density) * (local_part->values[0] - myPart[index_j].values[0]) * (d_x_ij * weight_x + d_y_ij * weight_y) /(d_ij*d_ij);
//       if (index_part == 5) printf("laplacian_part[k] = %2.6f \n",  d_ij);//laplacian_part[0]);

      List = List->next;
  }
// Don't free otherwize you loose the pointers to local_part->particle_derivatives
//   free(divergence_part);
//   free(grad_part);
//   free(laplacian_part);

}


void computeDerivatiesAllParticles(mySingleParticle* myPart, int nbParticles, kernelType kType) {
  // WARNING: this routine only valid for scalar quantity for the moment
  printf("\ni  (X,Y)			f(x)		Grad_x			Grad_y			Laplacian\n\n");
  for (int i=0; i<nbParticles; i++) {

      computeDerivativesOfParticleQuantity(myPart, i, kType);
      double* x = myPart[i].coordinates;
      double* grad = myPart[i].particle_derivatives->gradient;
      double* lapl = myPart[i].particle_derivatives->laplacian;
      double value = myPart[i].values[0];
      // Just print one line of particle along the x-direction
      if (x[1] > 0.0 && x[1] < 0.1) {
	printf("%d  (%2.3f,%2.3f)		%2.6f        %2.6f		%2.6f		%2.6f\n",i,x[0],x[1],value,grad[0],grad[1],lapl[0]);
      }

     }
}


void computeDerivativesAllParticles(allParticles* myPart, kernelType kType) {
  // WARNING: this routine only valid for scalar quantity for the moment
//   printf("\ni  (X,Y)			f(x)		Grad_x			Grad_y			Laplacian\n\n");
  int* index_in_domain = myPart->index_part_in_domain;
  int nbPart = myPart->nb_particles_domain;//NPTS_DOMAIN*NPTS_DOMAIN;//sizeof(index_in_domain)/sizeof(index_in_domain[0]);
  for (int i=0; i<nbPart; i++) {
      int local_index = index_in_domain[i];
      computeDerivativesOfParticleQuantity(myPart->array_of_particles, local_index, kType);
//       double* x = myPart->array_of_particles[i].coordinates;
//       double* grad = myPart->array_of_particles[i].particle_derivatives->gradient;
//       double* lapl = myPart->array_of_particles[local_index].particle_derivatives->laplacian;
//       double value = myPart->array_of_particles[i].values[0];
     // Just print one line of particle along the x-direction
//       if (x[1] > 0.5 && x[1] < 0.55) {
// 	printf("%d  (%2.3f,%2.3f)		%2.6f        %2.6f		%2.6f		%2.6f\n",i,x[0],x[1],value,grad[0],grad[1],lapl[0]);
//       }
//       printf(">>>>>>>>>>>> i = %d \n", i);
     }
}


void delete_all_particles(allParticles* myPart) {
  free(myPart->array_of_particles);
  free(myPart);
}

allParticles* combine_two_particles_sets(allParticles* part1, allParticles* part2, particleType pType1, particleType pType2) {
  allParticles* part_combined = malloc(sizeof(allParticles));
  int size_part1 = part1->nb_particles;
  int size_part2 = part2->nb_particles;
  int size_part_combined = size_part1 + size_part2;

  int* index_part1;
  if (pType1 == IN_DOMAIN) index_part1 = part1->index_part_in_domain;
  else if (pType1 == ON_BOUNDARY) index_part1 = part1->index_part_on_boundaries;

  int* index_part2;
  if (pType2 == IN_DOMAIN) index_part2 = part2->index_part_in_domain;
  else if (pType2 == ON_BOUNDARY) index_part2 = part2->index_part_on_boundaries;

  part_combined->array_of_particles = malloc(size_part_combined*sizeof(mySingleParticle));

  for (int i=0; i < size_part1; i++) {
    int local_index = index_part1[i];
    clone_single_particle(&(part_combined->array_of_particles[local_index]), &(part1->array_of_particles[i]));
  }

  for (int i=0; i < size_part2; i++) {
    int local_index = index_part2[i];
    clone_single_particle(&(part_combined->array_of_particles[local_index]), &(part2->array_of_particles[i]));
  }

  part_combined->nb_particles = size_part_combined;
  // WARNING: only work for one set of particle in the domain and the other one on the boundaries (no other combinations possible)
  if (pType1 == IN_DOMAIN) {
    part_combined->index_part_in_domain = index_part1;
    part_combined->index_part_on_boundaries = index_part2;
    part_combined->nb_particles_domain = part1->nb_particles_domain;
  }
  else{
    part_combined->index_part_in_domain = index_part2;
    part_combined->index_part_on_boundaries = index_part1;
    part_combined->nb_particles_boundaries = part2->nb_particles_boundaries;
  }
}


void clone_single_particle (mySingleParticle* part1, mySingleParticle* part2) {
  part1->coordinates = part2->coordinates;
  part1->values = part2->values;
  part1->size_values = part2->size_values;
  part1->density = part2->density;
  part1->mass = part2->mass;
  if (part2->particle_neighbours) part1->particle_neighbours = part2->particle_neighbours;
  if (part2->particle_derivatives) part1->particle_derivatives = part2->particle_derivatives;

}

void associate_neighborhood_to_particles(allParticles* part, neighborhood* nh) {
  int* index_part_domain = NULL;
  int* index_part_boundary = NULL;



  int nbPart_domain = part->nb_particles_domain;//NPTS_DOMAIN;
  int nbPart_boundary = part->nb_particles_boundaries;//NPTS_BOUNDARIES;

  if(part->index_part_in_domain) {
    index_part_domain = part->index_part_in_domain;
    for (int i=0; i<nbPart_domain; i++) {
      int local_index = index_part_domain[i];
      part->array_of_particles[local_index].particle_neighbours->nh = &(nh[local_index]);
    }
  }

  if(part->index_part_on_boundaries) {
    index_part_boundary = part->index_part_on_boundaries;
    for (int i=0; i<nbPart_boundary; i++) {
      int local_index = index_part_boundary[i];
      part->array_of_particles[local_index].particle_neighbours->nh = &(nh[local_index]);
    }
  }

}


void integrate_equation(allParticles* part, double dt, double alpha) {
  int* index_in_domain = part->index_part_in_domain;
  int nbPart = part->nb_particles_domain;//NPTS_DOMAIN*NPTS_DOMAIN;//sizeof(index_in_domain)/sizeof(index_in_domain[0]);
  for (int i = 0; i<nbPart; i++) {
      int local_index = index_in_domain[i];
      part->array_of_particles[local_index].values[0] += dt*alpha*part->array_of_particles[local_index].particle_derivatives->laplacian[0];

      if (part->array_of_particles[local_index].size_values > 1) {
	part->array_of_particles[local_index].values[1] += dt*alpha*part->array_of_particles[local_index].particle_derivatives->laplacian[1];
      }
  }
}


double solution_Fourier_series_Gaussian_source(double *x_coord, double time, double* args_fct_init, int nb_terms) {
    double x = x_coord[0];
    double y = x_coord[1];
    double length_x = 1.0;
    double length_y = 1.0;
    double T = 0.0;
    for (int m=0; m<nb_terms; m++) {
	for (int n=0; n<nb_terms; n++) {
	    double A_mn = get_A_mn_series_Gaussian_source(m,n, args_fct_init);
	    double lambda_mn = M_PI*M_PI*((m*m/(length_x*length_x)) + (n*n/(length_y*length_y)));
	    T += A_mn * sin((m*M_PI*x)/length_x) * sin((n*M_PI*y)/length_y) * exp(-lambda_mn * time);
	}
    }
    return T;
}

double get_A_mn_series_Gaussian_source(int m, int n, double* args_fct_init) {
  double length_x = 1.0;
  double length_y = 1.0;
  double intensity = args_fct_init[0];
  double sigma = args_fct_init[1];
  double x_c = args_fct_init[2];
  double y_c = args_fct_init[3];
  double first_integr = get_integration_Gaussian_source(m, x_c, length_x, sigma);
  double second_integr = get_integration_Gaussian_source(n, y_c, length_y, sigma);

  double A_mn = ((4*intensity) / (length_x*length_y)) * first_integr * second_integr;

  return A_mn;
}

double get_integration_Gaussian_source(int m, double x_c, double length_x, double sigma) {
  int n = 50;
  double a = 0.0;
  double b = length_x;
  double h=fabs(length_x-0.0)/n;
  double sum = 0;
  for(int i=1;i<n;i++){
    double x = a + i*h;
    sum=sum + fct_integration_Gaussian_source(x, m, x_c, length_x, sigma);
  }
  double integral = (h/2)*(fct_integration_Gaussian_source(a, m, x_c, length_x, sigma)+fct_integration_Gaussian_source(b, m, x_c, length_x, sigma)+2*sum);

  return integral;

}

double fct_integration_Gaussian_source(double x, int m, double x_c, double length_x, double sigma) {
  double term1 = exp((-1.0/(2*sigma*sigma))*(x-x_c)*(x-x_c));
  double term2 = sin((m*M_PI*x)/length_x);
  return term1*term2;
}
