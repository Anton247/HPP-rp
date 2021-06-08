#include <assert.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define FIELD(i, j) field[(i)*(WIDTH) + j]
#define NEW_FIELD(i, j) new_field[(i)*(WIDTH) + j]

#define RIGHT 2
#define DOWN 4
#define LEFT 8
#define UP 1
#define RP1 16
#define RP2 32
#define RP3 64
#define RP4 128

char *field, *new_field;
int ensemble;

void fill(int i0, int j0, int i1, int j1,
	double r, double d, double l, double u,
	double r1, double r2, double r3, double r4);

#include "config.c"
#include "collide.c"

inline int get_move(char cell)
{
	int mass=0;
	if (cell&RIGHT) mass+=1;
	if (cell&DOWN) mass+=1;
	if (cell&LEFT) mass+=1;
	if (cell&UP) mass+=1;
	return mass;
}

inline int get_rest(char cell)
{
	int mass=0;
	if (cell&RP1) mass+=2;
	if (cell&RP2) mass+=4;
	if (cell&RP3) mass+=8;
	if (cell&RP4) mass+=16;
	return mass;
}

void fill(int i0, int j0, int i1, int j1,
	double r, double d, double l, double u,
	double r1, double r2, double r3, double r4)
{
	int i, j;
	for (i=i0; i<i1; i++) {
		for (j=j0; j<j1; j++) {
			char cell=0;
			if (drand48() < r) { cell += RIGHT; }
			if (drand48() < d) { cell += DOWN; }
			if (drand48() < l) { cell += LEFT; }
			if (drand48() < u) { cell += UP; }
			if (drand48() < r1) { cell += RP1; }
			if (drand48() < r2) { cell += RP2; }
			if (drand48() < r3) { cell += RP3; }
			if (drand48() < r4) { cell += RP4; }
			FIELD(i, j)=cell;
		}
	}
}

void sum_mass(int row, int col, int* rest, int* move)
{
	int i, j;
	for (i=-AVERAGING_RADIUS; i<=AVERAGING_RADIUS; i++) {
		for (j=-AVERAGING_RADIUS; j<=AVERAGING_RADIUS; j++) {
			*rest+=get_rest(FIELD(i+row, j+col));
			*move+=get_move(FIELD(i+row, j+col));
		}
	}
}

double* load(char const* path)
{
	int j;
	FILE *fl;
	double *data=(double*)malloc(
		(WIDTH-2*AVERAGING_RADIUS-1)*sizeof(double)
	);

	fl=fopen(path, "r");																																																																			

	for(j=AVERAGING_RADIUS; j<WIDTH-AVERAGING_RADIUS-1; j++) {
		int row;
		double val;
		if (fscanf(fl, "%d\t%lf\n", &row, &val)!=2 || row!=j) {
			fprintf(stderr, "ERROR: File '%s' has illegal format",
				path);
			abort();
		} else {
			data[j-AVERAGING_RADIUS]=val;
		}
	}

	fclose(fl);
	return data;
}

void write(FILE *fl, int j, int ensemble, double *old, double val)
{
	assert(ensemble? !!old: !old);

	fprintf(fl, "%d\t%lf\n", j, 
		ensemble?
			(val+ensemble*old[j-AVERAGING_RADIUS])/(ensemble+1) 
			: val
	);
}

void save(int iter)
{
	int i, j;
	char density_path[100], move_path[100], rest_path[100];
	FILE *fl_density, *fl_move, *fl_rest;
	double square=(AVERAGING_RADIUS*2+1)*(AVERAGING_RADIUS*2+1);
	double *old_density=0, *old_move=0, *old_rest=0;

	sprintf(density_path, "density/%06d.xls", iter);
	sprintf(move_path, "move/%06d.xls", iter);
	sprintf(rest_path, "rest/%06d.xls", iter);


	if (ensemble>0) {
		old_density=load(density_path);
		old_move=load(move_path);
		old_rest=load(rest_path);
	}

	fl_density=fopen(density_path, "w");
	fl_move=fopen(move_path, "w");
	fl_rest=fopen(rest_path, "w");

	for (j=AVERAGING_RADIUS; j<WIDTH-AVERAGING_RADIUS-1; j++) {
		int rest=0, move=0;
		double new_density, new_move, new_rest;
		for (i=AVERAGING_RADIUS; i<HEIGHT-AVERAGING_RADIUS-1; i++) {
			sum_mass(i, j, &rest, &move);
		}

		
		new_density=1.0*(rest+move)/(HEIGHT-AVERAGING_RADIUS*2-1)/square;
		new_move=1.0*move/(HEIGHT-AVERAGING_RADIUS*2-1)/square;
		new_rest=1.0*rest/(HEIGHT-AVERAGING_RADIUS*2-1)/square;

		write(fl_density, j, ensemble, old_density, new_density);
		write(fl_move, j, ensemble, old_move, new_move);
		write(fl_rest, j, ensemble, old_rest, new_rest);

		
	}

	fclose(fl_density);
	fclose(fl_rest);
	fclose(fl_move);

	free(old_density);
	free(old_move);
	free(old_rest);
}

void swap_buffers()
{
	char* tmp;
	tmp=field;
	field=new_field;
	new_field=tmp;
}

void print_field(int e, int iter){
	char field_path[100];
	FILE *fl_field;
	sprintf(field_path, "Field/field%d_%d.txt", e, iter);
	fl_field=fopen(field_path, "w");
	int i0= 0,j0= 0, i1= HEIGHT, j1= WIDTH;
	for (int i=i0; i<i1; i++) {
		for (int j=j0; j<j1; j++) {
			int d = FIELD(i, j);
			fprintf(fl_field, "%d\n", d);
		}
	}
	fclose(fl_field);
}

int main()
{
	struct timespec t0, t1;

	field=(char*)malloc(WIDTH*HEIGHT);
	new_field=(char*)malloc(WIDTH*HEIGHT);

	clock_gettime(CLOCK_MONOTONIC_RAW, &t0);

	#pragma omp parallel num_threads(THREADS)
	{
	int iter;
	#pragma omp single
	printf("Running %d threads\n", omp_get_num_threads());

	for (ensemble=0; ensemble<ENSEMBLE_SIZE;)
	{
		#pragma omp single
		printf("ENSEMBLE: %d\n", ensemble);

		init();


		for (iter=0; iter<ITERS; iter++) {
			int i;

			#pragma omp single
			sources(iter);
			#pragma omp single
			print_field(ensemble, iter);


			#pragma omp for private(i) schedule(static)
			for(i=0; i<HEIGHT; i++) {
				int j;
				for (j=0; j<WIDTH; j++) {
					FIELD(i, j)=collide(FIELD(i, j));
				}
			}

			#pragma omp barrier

			#pragma omp for private(i) schedule(static)
			for (i=0; i<HEIGHT; i++) {
				int j;
				for (j=0; j<WIDTH; j++) {
					char cell=FIELD(i, j)&0xf0; // keep rest mass
					
					// propogate particles
					if (FIELD(i, (j+WIDTH-1)%WIDTH)&RIGHT) cell+=RIGHT; //field[(i)*(WIDTH) + j]
					if (FIELD(i, (j+1)%WIDTH)&LEFT) cell+=LEFT;
					if (FIELD((i+HEIGHT-1)%HEIGHT, j)&UP) cell+=UP;
					if (FIELD((i+1)%HEIGHT, j)&DOWN) cell+=DOWN;

					NEW_FIELD(i, j)=cell;
				}
			}

			#pragma omp barrier
			#pragma omp single
			{
			swap_buffers();

			if (iter%SAVE_PERIOD==0) {
				save(iter);
				printf("%d/%d\r", iter, ITERS);
				fflush(stdout);
			}
			}
			#pragma omp barrier
		}

		#pragma omp single
		ensemble++;
	}

	}

	clock_gettime(CLOCK_MONOTONIC_RAW, &t1);

	printf("Done %d iters, %lf MCell/sec\n", ITERS,
			1.0*ITERS*HEIGHT*WIDTH/
			((t1.tv_sec-t0.tv_sec)+0.000000001*(t1.tv_nsec-t0.tv_nsec))
			/1000000);

	return 0;
}
