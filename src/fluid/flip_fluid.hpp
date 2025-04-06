#ifndef __FLIP_FLUID_H__
#define __FLIP_FLUID_H__

#include <stdint.h>
#include <vector>

#define U_FIELD 0
#define V_FIELD 1

#define FLUID_CELL 0
#define AIR_CELL 1
#define SOLID_CELL 2

class FlipFluid
{
private:
	/* data */
public:
	float density;
	uint32_t fluid_num_x;
	uint32_t fluid_num_y;
	float h;
	float h1;
	float h2;
	uint32_t fluid_num_cells;

	std::vector<float> u;
	std::vector<float> v;
	std::vector<float> du;
	std::vector<float> dv;
	std::vector<float> prevU;
	std::vector<float> prevV;
	std::vector<float> p;
	std::vector<float> s;
	std::vector<uint32_t> cell_type;
	std::vector<float> cell_color;

	// Particle
	uint32_t max_particles;

	std::vector<float> particle_pos;
	std::vector<float> particle_color;
	std::vector<float> particle_vel;
	std::vector<float> particle_density;
	std::vector<float> particle_accel;
	float particle_rest_density;

	float particle_radius;
	float particle_inv_spacing;
	uint32_t particle_num_x;
	uint32_t particle_num_y;
	uint32_t particle_num_cells;

	std::vector<uint32_t> numCellParticles;
	std::vector<uint32_t> firstCellParticle;
	std::vector<uint32_t> cellParticleIds;

	uint32_t particle_nums;

public:
	FlipFluid(float density, float width, float height, float spacing, float particleRaius, uint32_t maxParticles);
	void integrate_particles(float dt, float gravity_x, float gravity_y);
	void push_particle_apart(uint32_t numIters);
	void handle_particle_collisions();
	void update_particle_density();
	void transfer_velocities(uint8_t toGrid, float flipRatio);
	void solve_incompressibility(uint32_t numIters, float dt, float overRelaxation, uint8_t compensateDrift);
	void updateParticleColor();
	void setSciColor(uint32_t cellNr, float val, float minVal, float maxVal);
	void updateCellColors();
	void simulate(float dt, float gravity_x, float gravity_y, float flipRatio, uint32_t numPressureIters, uint32_t numParticlesIters, float overRelaxation, uint8_t compensateDrift, uint8_t seperateParticles, float obstacleX, float obstacleY, float obstacleRadius);
};

#endif // __FLIP_FLUID_H__