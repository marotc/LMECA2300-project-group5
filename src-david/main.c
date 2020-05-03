#include "neighborhood_search.h"
#include "kernel.h"
#include "checkDerivatives.h"
#include "visualization.h"
#include <math.h>

int NPTS;

  // ***********************************************************************************************
  // ***************************** DERIVATIVES OF ANALYTICAL FUNCTION ******************************
  // ***********************************************************************************************

void deriveAnalyticalFunction() {

    // Number of particles
    int nb_particles_domain = 50*50;
    int nb_particles_boundaries = 0;
    int nbBoundaries = 0;
    int total_nb_particles = nb_particles_domain + nbBoundaries*nb_particles_boundaries;
    NPTS = total_nb_particles;

    // *** 1 *** ASSIGN POSITIONS, DENSITY, AND MASS TO EVERY PARTICLES (neighbourhoods not yet defined)

    // Creation of particles in the domain
    int nbParticles_domain[2] = {(int) sqrt(nb_particles_domain), (int) sqrt(nb_particles_domain)}; // number of particles in each direction (WARNING: same number of particles should be chosen for the moment)
    double dx =  1.0 / ((double)sqrt(nb_particles_domain) - 1.0);
    double domain_lim[4] = {0.0,1.0,0.0,1.0}; // limits of the computational domain
    int starting_index_part_in_domain = 0;
    double args_init_function[1] = {0.0}; // value of the temperature to be imposed everywhere in the domain. This argument is passed to the "initFunction" routine which will specify the values of each particle quantity
    int size_values = 1; // scalar temperature field with a single component

    allParticles* particles_everywhere = create_particles_in_domain(nbParticles_domain, domain_lim, size_values, args_init_function, myFunctionToDerive, starting_index_part_in_domain);

    // *** 2 *** CREATION OF THE NEIGHBOURHOODS OF EACH PARTICLE (in the domain and on the boundaries)

    // Creation of neighborhoods
    double timestep = 0.0;
    double maxspeed = 0.0;
    neighborhood_options* options = neighborhood_options_init(timestep, maxspeed);
    neighborhood* nh = options->nh;
    neighborhood_update_new(options, nh, particles_everywhere, 0);

    // Associate the computed neighbourhoods to each existing particle
    associate_neighborhood_to_particles(particles_everywhere, nh);


    // *** 3 *** COMPUTATION OF THE DERIVATIVES OF EVERY PARTICLES

    computeDerivativesAllParticles(particles_everywhere, CUBIC);

    // Print the values of the derivatives along an horizontal line of particules
    int* index_in_domain = particles_everywhere->index_part_in_domain;
    int nbPart = particles_everywhere->nb_particles_domain;//NPTS_DOMAIN*NPTS_DOMAIN;//sizeof(index_in_domain)/sizeof(index_in_domain[0]);
    printf("\ni  (X,Y)			f(x)		Grad_x			Grad_y			Laplacian\n\n");
    for (int i=0; i<nbPart; i++) {
      int local_index = index_in_domain[i];
      double* x = particles_everywhere->array_of_particles[i].coordinates;
      double* grad = particles_everywhere->array_of_particles[i].particle_derivatives->gradient;
      double* lapl = particles_everywhere->array_of_particles[local_index].particle_derivatives->laplacian;
      double value = particles_everywhere->array_of_particles[i].values[0];
      // Just print one line of particle along the x-direction
      if (x[1] > 0.5*domain_lim[1] && x[1] < 0.5*domain_lim[1]+dx) {
	printf("%d  (%2.3f,%2.3f)		%2.6f        %2.6f		%2.6f		%2.6f\n",i,x[0],x[1],value,grad[0],grad[1],lapl[0]);
      }
    }

    // Draw the particles with their values equal to the values of the analytical function we derived
    GLfloat(*data)[8] = malloc(sizeof(data[0]) * NPTS);
    CHECK_MALLOC(data);
    double extrema[2] = {0.0, 2.0};
    myFillData(data, particles_everywhere, extrema);
    draw_particles(domain_lim, data, total_nb_particles);

    // Free the different structures
    neighborhood_options_delete(options,nh);
    free(data);
    delete_all_particles(particles_everywhere);
//     return EXIT_SUCCESS;

}


  // ***********************************************************************************************
  // ************************** HEAT EQUATION IN A SQUARE WITH DIRICHLET B.C. **********************
  // ***********************************************************************************************

void solveHeatEquationInASquare(int problemChoice) {

    // Number of particles
    int nb_particles_domain = 50*50;
    int nb_particles_boundaries = 100;
    int nbBoundaries = 4;
    int total_nb_particles = nb_particles_domain + nbBoundaries*nb_particles_boundaries;
    NPTS = total_nb_particles;

    int problem_choice = problemChoice; // 1: heat equation with a hot boundary, 2: heat equation with a heat source in the center


    // *** 1 *** ASSIGN POSITIONS, DENSITY, AND MASS TO EVERY PARTICLES (neighbourhoods not yet defined) + BOUNDARY CONDITIONS

    // Creation of particles in the domain
    allParticles* particles_in_domain = NULL;

    int nbParticles_domain[2] = {(int) sqrt(nb_particles_domain), (int) sqrt(nb_particles_domain)}; // number of particles in each direction (WARNING: same number of particles should be chosen for the moment)
    double dx = 0.5 * 1.0 / ((double)sqrt(nb_particles_domain) - 1.0); // particles in the domain a bit shifted towards the inside of the domain to avoid overlapping with the particles on the boundaries defined after
    double domain_lim[4] = {0.0+dx,1.0-dx,0.0+dx,1.0-dx}; // limits of the computational domain
    int starting_index_part_in_domain = 0; // particles in the domain are the first ones in the global set of particles
    int size_values = 1; // scalar temperature field with a single component
    double alpha; // Thermal diffusivity
    double* args_init_function = NULL;
    if (problem_choice == 1) { // heat equation with a hot boundary and no heat source
      args_init_function = calloc(1, sizeof(double));
      *args_init_function = 0.0; // value of the temperature to be imposed everywhere in the domain. This argument is passed to the "initFunction" routine which will specify the values of each particle quantity
      particles_in_domain = create_particles_in_domain(nbParticles_domain, domain_lim, size_values, args_init_function, initFunction, starting_index_part_in_domain);
      alpha = 1.0;
    }
    else if (problem_choice == 2) { // heat equation with homogeneous Dirichlet B.C. and a heat source
      args_init_function = calloc(4, sizeof(double)); // arguments for the Gaussian function initializing the temperature profile in the domain: intensity, variance, x & y coordinates of its center
      *args_init_function = 1.0, *(args_init_function+1) = 0.1, *(args_init_function+2) = 0.5, *(args_init_function+3) = 0.5;
      particles_in_domain = create_particles_in_domain(nbParticles_domain, domain_lim, size_values, args_init_function, initFunction_GaussianSource, starting_index_part_in_domain);
      alpha = 1.0;
    }

    // Creation of particles on the boundaries
    allParticles* particles_on_boundaries = NULL;

    int nbParticles_boundaries[4] = {nb_particles_boundaries, nb_particles_boundaries, nb_particles_boundaries, nb_particles_boundaries}; // number of particles on each boundary (WARNING: same number of particles should be chosen for the moment)
    double boundaries_lim[4][4] = {{0.0, 0.0, 0.0, 1.0},
				   {0.0, 1.0, 1.0, 1.0},
				   {1.0, 1.0, 0.0, 1.0},
				   {0.0, 1.0, 0.0, 0.0}}; // (x_min, x_max, y_min, y_max) of each boundary
    int starting_index_part_on_bound = nb_particles_domain; // particles on the boundaries come after the ones in the domain in the global set of particles
    if (problem_choice == 1) {
      double values_Dirichlet[4] = {0.0, 0.0, 0.0, 100.0}; // Hot boundary initiliazed with a temperature = 100
      particles_on_boundaries = create_particles_on_boundaries(nbParticles_boundaries, (double *)boundaries_lim, size_values, nbBoundaries, values_Dirichlet, starting_index_part_on_bound);
    }
    else if (problem_choice == 2) {
      double values_Dirichlet[4] = {0.0, 0.0, 0.0, 0.0}; // Homogeneous Dirichlet B.C.
      particles_on_boundaries = create_particles_on_boundaries(nbParticles_boundaries, (double *)boundaries_lim, size_values, nbBoundaries, values_Dirichlet, starting_index_part_on_bound);
    }

    // Assemble the particles in the domain and on the boundaries in a single structure, for the computation of the neighbourhoods just after
    allParticles* particles_everywhere = combine_two_particles_sets(particles_in_domain, particles_on_boundaries, IN_DOMAIN, ON_BOUNDARY);


    // *** 2 *** CREATION OF THE NEIGHBOURHOODS OF EACH PARTICLE (in the domain and on the boundaries)

    // Creation of neighborhoods
    double timestep = 0.0;
    double maxspeed = 0.0;
    neighborhood_options* options = neighborhood_options_init(timestep, maxspeed);
    neighborhood* nh = options->nh;

    neighborhood_update_new(options, nh, particles_everywhere, 0);

    // Associate the computed neighbourhoods to each existing particle
    associate_neighborhood_to_particles(particles_everywhere, nh);

    int index_particle_to_check = nb_particles_domain * 0.5 - (int)(nbParticles_domain[0]*0.5);
//     display_neighbourhood_one_particle(particles_everywhere, index_particle_to_check, total_nb_particles); // Verification that the neighbourhood of a single chosen particle is well computed


    // *** 3 *** TIME LOOP TO RESOLVE THE EQUATIONS

    // Time stepping scheme parameters
    double time = 0.0;
    int nb_time_step = 0;
    double dt;
    double time_max;
    int nb_time_step_max;
    int print_every_time_step;
    double extrema[2] = {0.0, 1.0};
    if (problem_choice == 1) {
      dt = 1.0E-5;
      time_max = 1000*dt;
      nb_time_step_max = time_max / dt;
      print_every_time_step = 20;
      extrema[0] = 0.0, extrema[1] = 100.0;
    }
    else if (problem_choice == 2) {
      dt = 1.0E-5;
      time_max = 200*dt;
      nb_time_step_max = time_max / dt;
      print_every_time_step = 20;
      extrema[0] = 0.0, extrema[1] = 1.0;
    }

//     int size_solution_vector = (int)(nb_time_step_max / print_every_time_step) + 1;
//     double* solution_vector = malloc(size_values*size_solution_vector*sizeof(double));

    print_error_heat_equation(problemChoice, index_particle_to_check, nb_time_step, nb_time_step_max, time, particles_everywhere, args_init_function);
    // Time loop
    while (time < time_max) {
      // Compute derivatives of all the particles inside the domain, in particular the Laplacian here
      computeDerivativesAllParticles(particles_everywhere, CUBIC);
      // Update the temperature of every particles with an Euler explicit time integration scheme
      integrate_equation(particles_everywhere, dt, alpha);

      time += dt;
      nb_time_step++;
      // Print infos
      if (nb_time_step%print_every_time_step == 0)
	print_error_heat_equation(problemChoice, index_particle_to_check, nb_time_step, nb_time_step_max, time, particles_everywhere, args_init_function);
      // NOTE: No need to update the positions of the particles or to update the neighbourhoods since the particles are fixed
    }

    // Draw the particles with their values equal to the values of the analytical function we derived
    GLfloat(*data)[8] = malloc(sizeof(data[0]) * total_nb_particles);
    CHECK_MALLOC(data);
    myFillData(data, particles_everywhere, extrema);
    draw_particles(domain_lim, data, total_nb_particles);

    // Free the different structures
    neighborhood_options_delete(options,nh);
    free(data);
    delete_all_particles(particles_in_domain);
    delete_all_particles(particles_on_boundaries);
    delete_all_particles(particles_everywhere);
//     return EXIT_SUCCESS;

}

  // ***********************************************************************************************
  // ************************** MAIN FUNCTION - Choose the simu you want to run ********************
  // ***********************************************************************************************

int main()
{

//   deriveAnalyticalFunction();

  solveHeatEquationInASquare(2);

  return EXIT_SUCCESS;

}
