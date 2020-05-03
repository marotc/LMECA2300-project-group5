#include "visualization.h"

// for v, a floating point value between 0 and 1, this function fills color with
// the improved jet colormap color corresponding to v
void colormap(float v, float color[3])
{
	float v1 = 3.5 * (v - 0.7);
	float v2 = 1.25 * v;
	float v3 = fminf(0.5, v) * 2.0;

	color[0] = -v1 * v1 + 1.0f;
	color[1] = 6.0f * v2 * v2 * (1.0f - v2);
	color[2] = 5.5f * v3 * (1.0f - v3) * (1.0f - v3);

	// alternative: classical jet colormap
	// color[0] = 1.5 - 4.0 * fabs(v - 0.75);
	// color[1] = 1.5 - 4.0 * fabs(v - 0.5 );
	// color[2] = 1.5 - 4.0 * fabs(v - 0.25);
}

void myColormap(float v, float color[3], float v_max, float v_min)
{
	v -= v_min;
	v /= (v_max - v_min);

	// alternative: classical jet colormap
	color[0] = 1.5 - 4.0 * fabs(v - 0.75);
	color[1] = 1.5 - 4.0 * fabs(v - 0.5 );
	color[2] = 1.5 - 4.0 * fabs(v - 0.25);
}

// function to fill the data table of the nPoints particles positions, speeds, colors and transparency and the coord table with the nPoints particles positions used to draw;
// data[i][0] == coord[i][0] && data[i][1] == coord[i][1]
void fillData(GLfloat(* data)[8], int nbParticles)
{
	double x_lim[4] = {0.0, 1.0, 0.0, 1.0};
	int i=0;
	int nbPart_x = sqrt(nbParticles);
	int nbPart_y = nbPart_x;
	for (int ind_y = 0; ind_y < nbPart_y; ind_y++) {
	  for (int ind_x = 0; ind_x < nbPart_x; ind_x++) {
// 	  double coord[2] = {data[i][0], data[i][1]};
	  double* coord = calloc(2,sizeof(double));
// 	  double values[2] = {0.0,0.0};
	  double* values = calloc(2,sizeof(double));
	    double mass = 1.0;
	    double density = 1.0;

// 	    init1DSegmentWithParticles(x_lim, coord, values, &mass, &density, nbParticles, i, 1);
	    initSquareWithParticles(x_lim, coord, values, &mass, &density, nbPart_x, nbPart_y, ind_x, ind_y, 1);
	    data[i][0] = coord[0];
	    data[i][1] = coord[1];

	    data[i][2] = 0.0;
	    data[i][3] = 0.0;
	    double r = values[0];//sqrt(values[0] * values[0]+ values[1] * values[1]);
	    myColormap(r, &data[i][4], 2.0, -2.0); // fill color
	    data[i][7] = 0.8f; // transparency

	    i++;
// 	    printf("(values[0], values[1]) = (%2.6f,%2.6f) \n", values[0], values[1]);
	}
      }
}

void myFillData(GLfloat(* data)[8], allParticles* my_array_of_particles, double* extrema)
{
	int nbPart = my_array_of_particles->nb_particles;
	for (int i=0; i<nbPart; i++) {
	    data[i][0] = my_array_of_particles->array_of_particles[i].coordinates[0];
	    data[i][1] = my_array_of_particles->array_of_particles[i].coordinates[1];
	    double r = my_array_of_particles->array_of_particles[i].values[0];
	    myColormap(r, &data[i][4], extrema[1], extrema[0]);
	    data[i][7] = 0.8f; // transparency
	}
}

void draw_particles(double* x_lim, GLfloat(* data)[8], int nbParticles) {
      bov_window_t* window = bov_window_new(1024, 780, "ANM Project: SPH");
      bov_window_set_color(window, (GLfloat[]){0.9f, 0.85f, 0.8f, 0.0f});

      /* send data to GPU, and receive reference to those data in a points object */
      bov_points_t *particles = bov_particles_new(data, nbParticles, GL_STATIC_DRAW);

      /* setting particles appearance */
      bov_points_set_width(particles, 0.01);
      bov_points_set_outline_width(particles, 0.0025);
//       bov_points_set_outline_color(particles, ());


      /* fit particles in a [-0.8 0.8]x[-0.8 0.8] square */
      bov_points_scale(particles, (GLfloat[2]) {1.0, 1.0});
      bov_points_set_pos(particles, (GLfloat[2]) {0.0, -0.25});

      /* we got 0.2 at the top to write something. The screen goes from -1 to 1 */
      bov_text_t* msg =  bov_text_new(
	      (GLubyte[]) {"Temperature field"},
	      GL_STATIC_DRAW);
      bov_text_set_pos(msg, (GLfloat[2]){-0.95, 0.82});
      bov_text_set_fontsize(msg, 0.1);

      while(!bov_window_should_close(window)){
// 	      fillData(data);
// 	      bov_particles_update(particles, data, NPTS);

	      bov_particles_draw(window, particles, 0, BOV_TILL_END);
	      // bov_points_draw(window, particles, 0, BOV_TILL_END);
	      // bov_lines_draw(window, particles, 0, BOV_TILL_END);
	      // bov_triangles_draw(window, particles, 0, BOV_TILL_END);

	      bov_text_draw(window, msg);

	      // In your actual project, don't wait for events => bov_window_update(window)
	      bov_window_update_and_wait_events(window);
// 	      bov_window_update(window);
      }


      bov_text_delete(msg);
      bov_points_delete(particles);
      bov_window_delete(window);
}

void create_window_animation(GLfloat(* data)[8], bov_window_t* window, bov_points_t *particles, int nbParticles) {
      window = bov_window_new(1024, 780, "ANM Project: SPH");
      bov_window_set_color(window, (GLfloat[]){0.9f, 0.85f, 0.8f, 0.0f});

      /* send data to GPU, and receive reference to those data in a points object */
      particles = bov_particles_new(data, nbParticles, GL_STATIC_DRAW);

      /* setting particles appearance */
      bov_points_set_width(particles, 0.01);
      bov_points_set_outline_width(particles, 0.0025);
//       bov_points_set_outline_color(particles, ());


      /* fit particles in a [-0.8 0.8]x[-0.8 0.8] square */
      bov_points_scale(particles, (GLfloat[2]) {1.0, 1.0});
      bov_points_set_pos(particles, (GLfloat[2]) {0.0, -0.25});

      /* we got 0.2 at the top to write something. The screen goes from -1 to 1 */
      bov_text_t* msg =  bov_text_new(
	      (GLubyte[]) {"Temperature field"},
	      GL_STATIC_DRAW);
      bov_text_set_pos(msg, (GLfloat[2]){-0.95, 0.82});
      bov_text_set_fontsize(msg, 0.1);
//       while(!bov_window_should_close(window)){
      double tbegin = bov_window_get_time(window);
      while (bov_window_get_time(window) - tbegin < 2.0) {
	  bov_particles_draw(window, particles, 0, BOV_TILL_END);
	  bov_text_draw(window, msg);
	  // In your actual project, don't wait for events => bov_window_update(window)
// 	  bov_window_update_and_wait_events(window);
	      bov_window_update(window);
      }
}

void display_particles(GLfloat(* data)[8], bov_window_t* window, bov_points_t *particles, int end, int nbParticles) {
  float transition_time = 1.0;
  bov_points_t* new_particles = bov_particles_update(particles,data,nbParticles);
  bov_window_t* new_window = window;
  double tbegin = bov_window_get_time(new_window);
  if (!end) {
    while (bov_window_get_time(new_window) - tbegin < transition_time) {
      bov_particles_draw(new_window, new_particles, 0, BOV_TILL_END);
      bov_window_update(new_window);
    }
  }
  else {
    while (!bov_window_should_close(new_window)) {
      bov_particles_draw(new_window, new_particles, 0, BOV_TILL_END);
      bov_window_update_and_wait_events(new_window);
    }
  }
}


void display_neighbourhood_one_particle(allParticles* allPart, int index_part, int nbParticles) {
  mySingleParticle* local_part = &(allPart->array_of_particles[index_part]);
  int nNeigh = local_part->particle_neighbours->nh->nNeighbours;
  neighbours* List = local_part->particle_neighbours->nh->list;
  double kh = local_part->particle_neighbours->kh;

  local_part->values[0] = 1.0;
  for (int j = 0; j < nNeigh; j++) {
      int index_j = List->index;
      printf("index = %d \n", index_j);
      allPart->array_of_particles[index_j].values[0] = 0.5;
      List = List->next;
  }
  double extrema[2] = {0.0, 1.0};
  GLfloat(*data)[8] = malloc(sizeof(data[0]) * nbParticles);
  CHECK_MALLOC(data);
  myFillData(data, allPart, extrema);
  double domain_lim[4] = {0.0, 1.0, 0.0, 1.0};
  draw_particles(domain_lim, data, nbParticles);
}


void print_error_heat_equation(int problemChoice, int index_particle_to_check, int nb_time_step, int nb_time_step_max, double time, allParticles* particles_everywhere, double* args_init_function) {
  if (problemChoice == 1) {
    printf("------ Time step %d over %d -------\n", nb_time_step, nb_time_step_max);
  }
  else if (problemChoice == 2) {
    printf("------ Time step %d over %d -------\n", nb_time_step, nb_time_step_max);
    double T_SPH = particles_everywhere->array_of_particles[index_particle_to_check].values[0];
    double* x_SPH = particles_everywhere->array_of_particles[index_particle_to_check].coordinates;
    printf("T_SPH = %2.10f @ (x,y) = (%2.3f, %2.3f) \n", T_SPH, x_SPH[0], x_SPH[1]);
    double T_exact = solution_Fourier_series_Gaussian_source(x_SPH, time, args_init_function, 10);
    printf("T_exact = %2.10f @ (x,y) = (%2.3f, %2.3f) \n", T_exact, x_SPH[0], x_SPH[1]);
    double relative_error = (fabs(T_SPH-T_exact)/T_exact) * 100;
    printf("||error|| = %2.10f percent \n", relative_error);
  }

}
