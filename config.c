// Generated by system/gen_config_flatwave.py
#define THREADS 8
#define WIDTH 100
#define HEIGHT 100
#define ITERS 10
#define SAVE_PERIOD 5
#define AVERAGING_RADIUS 1
#define ENSEMBLE_SIZE 2
#define SOURCE_WIDTH 4

// Field initialization
void init()
{
	fill(0, 0, HEIGHT, WIDTH, 0.7, 0.7, 0.7, 0.7, 0.25, 0, 0, 0);
}

// Sources rules
void sources(int iter)
{
	if (iter==0) {
		fill(0, WIDTH/2-SOURCE_WIDTH/2, HEIGHT, WIDTH/2+SOURCE_WIDTH/2, 1, 1, 1, 1, 0.75, 0, 0, 0);
	}
}