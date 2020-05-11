#include "print_particules.h"

void colormap_Cs(Particle *p, float color[3]);

// Fills data with particle data
void fillData(GLfloat(*data)[8], Particle** particles, int N) {
// 	float rmax = 100.0*sqrtf(2.0f);
	double max_norm_vel = -1;
	int i_max = -1;
	double vel_norm_local;
	for (int i = 0; i < N; i++) {
		vel_norm_local = norm(particles[i]->v);
		if (vel_norm_local > max_norm_vel) {
			max_norm_vel = vel_norm_local;
			i_max = i;
		}
	}
	// printf("max velocity norm = %lf, pos = (%lf, %lf)\n", max_norm_vel, particles[i_max]->pos->x, particles[i_max]->pos->y);
	for (int i = 0; i < N; i++) {
		Particle* p = particles[i];
		if(!p->on_boundary) {
			data[i][0] = p->pos->x;
			data[i][1] = p->pos->y;
			data[i][2] = p->v->x;
			data[i][3] = p->v->y;
			// colormap_velocity(p, &data[i][4], max_norm_vel);
			colormap_init_pos(p, &data[i][4]);
			data[i][7] = 1; // transparency
		}
	}
}

bov_points_t * load_box() {
	GLfloat(*data)[2] = malloc(sizeof(data[0]) * 2 * 4);
	data[0][0] = 0, data[0][1] = 0;
	data[1][0] = 1, data[1][1] = 0;
	data[2][0] = 1, data[2][1] = 0;
	data[3][0] = 1, data[3][1] = 1;
	data[4][0] = 1, data[4][1] = 1;
	data[5][0] = 0, data[5][1] = 1;
	data[6][0] = 0, data[6][1] = 1;
	data[7][0] = 0, data[7][1] = 0;

	bov_points_t *points = bov_points_new(data, 2 * 4, GL_STATIC_DRAW);
	bov_points_set_width(points, 0.004);
	free(data);
	return points;
}


bov_points_t * load_Grid(Grid* grid,double scale)
{
	int nLines = (grid->nCellx + 1) + (grid->nCelly + 1);
	GLfloat(*data)[2] = malloc(sizeof(data[0]) * 2 * nLines);
	for (int i = 0;i < (grid->nCellx + 1);i++)
	{
		data[2 * i][0] = grid->left + i * grid->h;
		data[2 * i][1] = grid->bottom;
		data[2 * i + 1][0] = grid->left + i * grid->h;
		data[2 * i + 1][1] = grid->top;
	}
	int shift = 2 * (grid->nCellx + 1);
	for (int j = 0;j < (grid->nCelly + 1);j++)
	{
		data[shift + 2 * j][0] = grid->left;
		data[shift + 2 * j][1] = grid->bottom + j * grid->h;
		data[shift + 2 * j + 1][0] = grid->right;
		data[shift + 2 * j + 1][1] = grid->bottom + j * grid->h;
	}
	bov_points_t *points = bov_points_new(data, 2 * nLines, GL_STATIC_DRAW);
	bov_points_set_width(points, 0.005);
	double L = grid->h * grid->nCellx;
	// bov_points_scale(points, (GLfloat[2]){0.8/L*scale, 0.8/L*scale});
	// bov_points_scale(points, (GLfloat[2]) { 0.008, 0.008 });
	free(data);
	return points;
}

Animation* Animation_new(int N, double timeout,Grid* grid, double particle_width)
{
	Animation* animation = (Animation*)malloc(sizeof(Animation));
	animation->window = bov_window_new(1024, 780, "ANM Project: SPH");
	bov_window_translate(animation->window, (GLfloat[]) { -0.5f, -0.5f });
	bov_window_set_zoom(animation->window, 1.5);
	// 	bov_window_set_color(animation->window, (GLfloat[]) { 0.9f, 0.85f, 0.8f, 0.0f });
	bov_window_set_color(animation->window, (GLfloat[]) { 1.0f, 1.0f, 1.0f, 0.0f });
	bov_window_enable_help(animation->window);
	animation->N = N;
	animation->timeout = timeout;
	// double L = grid->h*grid->nCellx;
	////set-up particles////
	GLfloat(*data)[8] = malloc(sizeof(data[0])*N);
	bov_points_t *particles = bov_particles_new(data, N, GL_STATIC_DRAW);
	free(data);
	// setting particles appearance
	bov_points_set_width(particles, particle_width);
	// bov_points_set_outline_width(particles, 0.0025);

	double c = 4;
	// bov_points_scale(particles, (GLfloat[2]){0.4*c/L*scale, 0.4*c/L*scale});//0.8
	//bov_points_scale(particles, (GLfloat[2]){ 0.008, 0.008 });
	animation->particles = particles;
	////set-up grid////
	if (grid != NULL)
		animation->grid = load_Grid(grid, c);
	else
		animation->grid = NULL;
	animation->box = load_box();
	return animation;
}
void Animation_free(Animation* animation)
{
	bov_points_delete(animation->particles);
	if(animation->grid != NULL)
		bov_points_delete(animation->grid);
	bov_window_delete(animation->window);
	free(animation);
}

void display_particles(Particle** particles, Animation* animation, bool end, bool take_screenshot, double t) {
	int N = animation->N;
	GLfloat(*data)[8] = malloc(sizeof(data[0])*N);
	fillData(data, particles, N);
	animation->particles = bov_particles_update(animation->particles,data,N);
	free(data);

	// path for screenshot
	char screenshot_path[128];
	sprintf(screenshot_path, "../screenshots/t=%010lf", t);

	// prepare annotations
	GLubyte message_time[128], message_upper[128], message_lower[128];
	sprintf(message_time, "t = %.2lf", t);
	if(t < 60) sprintf(message_lower, "Lower lid is on");
	else sprintf(message_lower, "Lower lid is off");
	if(t >= 30 && t < 90) sprintf(message_upper, "Upper lid is on");
	else sprintf(message_upper, "Upper lid is off");

	bov_text_t* text_time = bov_text_new(message_time, GL_STATIC_DRAW);
	bov_text_t* text_upper = bov_text_new(message_upper, GL_STATIC_DRAW);
	bov_text_t* text_lower = bov_text_new(message_lower, GL_STATIC_DRAW);
	bov_text_set_pos(text_time, (GLfloat[]){0,1.07});
	bov_text_set_pos(text_upper, (GLfloat[]){0,1.01});
	bov_text_set_pos(text_lower, (GLfloat[]){0,-0.04});


	bov_window_t* window = animation->window;
	double tbegin = bov_window_get_time(window);
	if (!end){
		while (bov_window_get_time(window) - tbegin < animation->timeout) {
			bov_lines_draw(window,animation->box,0, BOV_TILL_END);
			if(animation->grid != NULL)
				bov_lines_draw(window,animation->grid,0, BOV_TILL_END);
			bov_particles_draw(window, animation->particles, 0, BOV_TILL_END);
			bov_text_draw(window, text_time);
			bov_text_draw(window, text_upper);
			bov_text_draw(window, text_lower);
			if(take_screenshot) bov_window_screenshot(window, screenshot_path);

			bov_window_update(window);
		}
	}
	else {
		// we want to keep the window open with everything displayed
		while (!bov_window_should_close(window)) {
			bov_lines_draw(window,animation->box,0, BOV_TILL_END);
			if (animation->grid != NULL)
				bov_lines_draw(window, animation->grid, 0, BOV_TILL_END);
			bov_particles_draw(window, animation->particles, 0, BOV_TILL_END);
			bov_text_draw(window, text_time);
			bov_text_draw(window, text_upper);
			bov_text_draw(window, text_lower);
			if(take_screenshot) bov_window_screenshot(window, screenshot_path);
			bov_window_update_and_wait_events(window);
		}
	}
}

// Fills color with the color to use for particle p
void colormap_cell(Particle* p, float color[3]) {
	if (p->cell == NULL) {
		color[0] = 0;color[1] = 0;color[2] = 0;
	}
	else {
		if (p->cell->i % 2 == 0) {
			color[0] = 20;
			if (p->cell->j % 2 == 0)
				color[1] = 20;
			else
				color[1] = 0;
		}
		else {
			color[0] = 0;
			if (p->cell->j % 2 == 0)
				color[1] = 20;
			else
				color[1] = 0;
		}
		color[2] = 0;
	}
}

void colormap_uni_color(float color[3])
{
	color[0] = 0;color[1] = 10;color[2] = 20;

}

void colormap_uni_color_2(float color[3]) {
	color[0] = 0;color[1] = 0;color[2] = 20;
}

void colormap_Cs(Particle *p, float color[3]) {
	color[0] = 20*squared(p->Cs);
	color[1] = 0;
	color[2] = 20*squared(1.0-p->Cs);
}

void colormap_fs(Particle *p, float color[3], double max_norm) {
	xy* fs = xy_new(-p->param->sigma * p->normal->x * p->kappa, -p->param->sigma * p->normal->y * p->kappa);
// 	double fs_norm = 1.0;//norm(fs);
	color[0] = 10*squared(fs->x/max_norm) + 10*squared(fs->y/max_norm);
	color[1] = 0.0;
	color[2] = 0.0;//20*squared(fs->y/max_norm);
}

void colormap_jet(double x, float color[3]) {
	color[0] = fmin(4*x - 1.5, -4*x + 4.5);
	color[1] = fmin(4*x - 0.5, -4*x + 3.5);
	color[2] = fmin(4*x + 0.5, -4*x + 2.5);
}

void colormap_velocity(Particle *p, float color[3], double max_norm) {
	double x = norm(p->v);
	colormap_jet(x, color);
}

void colormap_pressure(Particle *p, float color[3], double max_P) {
	color[0] = 20*p->v->x / max_P;
	if (p->on_free_surface)
	  color[1] = 1.0;//10*(abs(p->rho - max_P)) / max_P;
	else
	  color[1] = 0.0;
	color[2] = 0;
}

void colours_neighbors(GLfloat(*data)[8], Particle** particles, int index)
{
	Particle* p = particles[index];
	ListNode *node = p->neighborhood->head;
	data[index][4] = 20;data[index][5] = 20;data[index][6] = 20;
	while (node != NULL) {
		Particle* q = (Particle*)node->v;
		int i = q->index;
		data[i][4] = 0;data[i][5] = 20;data[i][6] = 20;
		node = node->next;
	}
}

void colormap_init_pos(Particle *p, float color[3]) {
	colormap_jet(p->init_pos->y, color);
}
