#ifndef VISUALIZATION
#define VISUALIZATION

#include "BOV.h"
#include "checkDerivatives.h"

void colormap(float v, float color[3]);

void myColormap(float v, float color[3], float v_max, float v_min);

void fillData(GLfloat(* data)[8], int nbParticles);

void myFillData(GLfloat(* data)[8], allParticles* my_array_of_particles, double* extrema);

void draw_particles(double* x_lim, GLfloat(* data)[8], int nbParticles);

void create_window_animation(GLfloat(* data)[8], bov_window_t* window, bov_points_t *particles, int nbParticles);

void display_particles(GLfloat(* data)[8], bov_window_t* window, bov_points_t *particles, int end, int nbParticles);

void display_neighbourhood_one_particle(allParticles* allPart, int index_part, int nbParticles);

void print_error_heat_equation(int problemChoice, int index_particle_to_check, int nb_time_step, int nb_time_step_max, double time, allParticles* particles_everywhere, double* args_init_function);

#endif