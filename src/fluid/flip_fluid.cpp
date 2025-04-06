#include "flip_fluid.hpp"

static inline float clamp(float x, float a, float b)
{
	return x < a ? a : (x > b ? b : x);
}

static inline float max(float a, float b)
{
	return a > b ? a : b;
}
FlipFluid::FlipFluid(float density, float width, float height, float spacing, float particleRaius, uint32_t maxParticles)
{
	// Fluid
	this->density = density;
	this->fluid_num_x = static_cast<uint32_t>(floor(width / spacing)) + 1;
	this->fluid_num_y = static_cast<uint32_t>(floor(height / spacing)) + 1;
	this->h = max(width / this->fluid_num_x, height / this->fluid_num_y);
	this->h1 = 1.0f / this->h;
	this->h2 = 0.5f * this->h;
	this->fluid_num_cells = this->fluid_num_x * this->fluid_num_y;

	this->u.resize(this->fluid_num_cells, 0.0f);
	this->v.resize(this->fluid_num_cells, 0.0f);
	this->du.resize(this->fluid_num_cells, 0.0f);
	this->dv.resize(this->fluid_num_cells, 0.0f);
	this->prevU.resize(this->fluid_num_cells, 0.0f);
	this->prevV.resize(this->fluid_num_cells, 0.0f);
	this->p.resize(this->fluid_num_cells, 0.0f);
	this->s.resize(this->fluid_num_cells, 0.0f);
	this->cell_type.resize(this->fluid_num_cells, FLUID_CELL);
	this->cell_color.resize(3 * this->fluid_num_cells, 0.0f);

	// Particle
	this->max_particles = maxParticles;

	this->particle_pos.resize(2 * this->max_particles, 0.0f);
	this->particle_color.resize(3 * this->max_particles, 0.0f);
	for (uint32_t i = 0; i < this->max_particles; i++)
	{
		this->particle_color[3 * i + 2] = 1.0f;
	}
	this->particle_vel.resize(2 * this->max_particles, 0.0f);
	this->particle_accel.resize(2 * this->max_particles, 0.0f);
	this->particle_density.resize(this->fluid_num_cells, 0.0f);
	this->particle_rest_density = 0.0f;

	this->particle_radius = particleRaius;
	this->particle_inv_spacing = 1.0f / (2.2f * particle_radius);
	this->particle_num_x = static_cast<uint32_t>(floor(width * this->particle_inv_spacing) + 1);
	this->particle_num_y = static_cast<uint32_t>(floor(height * this->particle_inv_spacing) + 1);
	this->particle_num_cells = this->particle_num_x * this->particle_num_y;

	this->numCellParticles.resize(this->particle_num_cells, 0);
	this->firstCellParticle.resize(this->particle_num_cells + 1, 0);
	this->cellParticleIds.resize(maxParticles, 0);

	this->particle_nums = 0;
}

void FlipFluid::integrate_particles(float dt, float gravity_x, float gravity_y)
{
	for (uint32_t i = 0; i < this->particle_nums; i++)
	{
		/* Use full Verlet velocity */
		this->particle_pos[2 * i] += this->particle_vel[2 * i] * dt + 0.5f * this->particle_accel[2 * i] * dt * dt;
		this->particle_pos[2 * i + 1] += this->particle_vel[2 * i + 1] * dt + 0.5f * this->particle_accel[2 * i + 1] * dt * dt;
		this->particle_vel[2 * i] += 0.5f * (this->particle_accel[2 * i] + gravity_x) * dt;
		this->particle_vel[2 * i + 1] += 0.5f * (this->particle_accel[2 * i + 1] + gravity_y) * dt;
		this->particle_accel[2 * i] = gravity_x;
		this->particle_accel[2 * i + 1] = gravity_y;
	}
}

void FlipFluid::push_particle_apart(uint32_t numIters)
{
	const float colorDiffusionCoeff = 0.001f;
	std::fill(this->numCellParticles.begin(), this->numCellParticles.end(), 0);
	
	/* Count particles in each cell */
	for (uint32_t i = 0; i < this->particle_nums; i++)
	{
		const float x = this->particle_pos[2 * i];
		const float y = this->particle_pos[2 * i + 1];
		const uint32_t xi = (uint32_t)clamp(floor(x * this->particle_inv_spacing), 0, this->particle_num_x - 1);
		const uint32_t yi = (uint32_t)clamp(floor(y * this->particle_inv_spacing), 0, this->particle_num_y - 1);
		const uint32_t cellNr = xi * this->particle_num_y + yi;
		this->numCellParticles[cellNr]++;
	}

	// Calculate partial sums
	uint32_t partial_sums = 0;
	for (uint32_t i = 0; i < this->particle_num_cells; i++)
	{
		partial_sums += this->numCellParticles[i];
		this->firstCellParticle[i] = partial_sums;
	}
	this->firstCellParticle[this->particle_num_cells] = partial_sums; // guard

	// Fill Particle in cell
	for (uint32_t i = 0; i < this->particle_nums; i++)
	{
		const float x = this->particle_pos[2 * i];
		const float y = this->particle_pos[2 * i + 1];
		const uint32_t xi = (uint32_t)clamp(floor(x * this->particle_inv_spacing), 0, this->particle_num_x - 1);
		const uint32_t yi = (uint32_t)clamp(floor(y * this->particle_inv_spacing), 0, this->particle_num_y - 1);
		const uint32_t cellNr = xi * this->particle_num_y + yi;
		this->firstCellParticle[cellNr]--;
		this->cellParticleIds[this->firstCellParticle[cellNr]] = i;
	}

	// Push particles apart
	const float minDist = 2.0f * this->particle_radius;
	const float minDist2 = minDist * minDist;
	for (uint32_t iter = 0; iter < numIters; iter++)
	{
		for (uint32_t i = 0; i < this->particle_nums; i++)
		{
			const float px = this->particle_pos[2 * i];
			const float py = this->particle_pos[2 * i + 1];
			const uint32_t pxi = (uint32_t)floor(px * this->particle_inv_spacing);
			const uint32_t pyi = (uint32_t)floor(py * this->particle_inv_spacing);
			const uint32_t x0 = (uint32_t)max(pxi - 1, 0);
			const uint32_t y0 = (uint32_t)max(pyi - 1, 0);
			const uint32_t x1 = (uint32_t)std::min(pxi + 1, this->particle_num_x - 1);
			const uint32_t y1 = (uint32_t)std::min(pyi + 1, this->particle_num_y - 1);

			for (uint32_t xi = x0; xi <= x1; xi++)
			{
				for (uint32_t yi = y0; yi <= y1; yi++)
				{
					const uint32_t cellNr = xi * this->particle_num_y + yi;
					const uint32_t first = this->firstCellParticle[cellNr];
					const uint32_t last = this->firstCellParticle[cellNr + 1];
					for (uint32_t j = first; j < last; j++)
					{
						const uint32_t id = this->cellParticleIds[j];
						if (id == i)
							continue;
						const float qx = this->particle_pos[2 * id];
						const float qy = this->particle_pos[2 * id + 1];

						float dx = qx - px;
						float dy = qy - py;
						const float d2 = dx * dx + dy * dy;
						if (d2 > minDist2 || d2 == 0.0f)
							continue;
						const float d = sqrt(d2);
						const float s = 0.5f * (minDist - d) / d;
						dx *= s;
						dy *= s;
						this->particle_pos[2 * i] -= dx;
						this->particle_pos[2 * i + 1] -= dy;
						this->particle_pos[2 * id] += dx;
						this->particle_pos[2 * id + 1] += dy;

						// Diffuse color
						for (uint32_t k = 0; k < 3; k++)
						{
							const float color0 = this->particle_color[3 * i + k];
							const float color1 = this->particle_color[3 * id + k];
							const float color = (color0 + color1) * 0.5f;
							this->particle_color[3 * i + k] = color0 + (color - color0) * colorDiffusionCoeff;
							this->particle_color[3 * id + k] = color1 + (color - color1) * colorDiffusionCoeff;
						}
					}
				}
			}
		}
	}
}

void FlipFluid::handle_particle_collisions()
{
	const float h = this->h;
	const float r = this->particle_radius;

	const float minX = h + r;
	const float maxX = (this->fluid_num_x - 1) * h - r;
	const float minY = h + r;
	const float maxY = (this->fluid_num_y - 1) * h - r;

	for (uint32_t i = 0; i < this->particle_nums; i++)
	{
		float x = this->particle_pos[2 * i];
		float y = this->particle_pos[2 * i + 1];

		// wall collisions
		if (x < minX)
		{
			x = minX;
			this->particle_vel[2 * i] = 0.0f;
		}
		if (x > maxX)
		{
			x = maxX;
			this->particle_vel[2 * i] = 0.0f;
		}
		if (y < minY)
		{
			y = minY;
			this->particle_vel[2 * i + 1] = 0.0f;
		}
		if (y > maxY)
		{
			y = maxY;
			this->particle_vel[2 * i + 1] = 0.0f;
		}
		this->particle_pos[2 * i] = x;
		this->particle_pos[2 * i + 1] = y;
	}
}

void FlipFluid::update_particle_density()
{
	const uint32_t n = this->fluid_num_y;
	const float h = this->h;
	const float h1 = this->h1;
	const float h2 = this->h2;
	std::vector<float> &d = particle_density;
	std::fill(d.begin(), d.end(), 0.0f);

	for (uint32_t i = 0; i < this->particle_nums; i++)
	{
		float x = this->particle_pos[2 * i];
		float y = this->particle_pos[2 * i + 1];
		x = clamp(x, h, (this->fluid_num_x - 1) * h);
		y = clamp(y, h, (this->fluid_num_y - 1) * h);

		const uint32_t x0 = (uint32_t)floor((x - h2) * h1);
		const uint32_t x1 = std::min(x0 + 1, this->fluid_num_x - 2);
		const uint32_t y0 = (uint32_t)floor((y - h2) * h1);
		const uint32_t y1 = std::min(y0 + 1, this->fluid_num_y - 2);

		const float tx = ((x - h2) - x0 * h) * h1;
		const float ty = ((y - h2) - y0 * h) * h1;
		const float sx = 1.0f - tx;
		const float sy = 1.0f - ty;

		if (x0 < this->fluid_num_x && y0 < this->fluid_num_y)
			d[x0 * n + y0] += sx * sy;
		if (x1 < this->fluid_num_x && y0 < this->fluid_num_y)
			d[x1 * n + y0] += tx * sy;
		if (x1 < this->fluid_num_x && y1 < this->fluid_num_y)
			d[x1 * n + y1] += tx * ty;
		if (x0 < this->fluid_num_x && y1 < this->fluid_num_y)
			d[x0 * n + y1] += sx * ty;
	}

	if (this->particle_rest_density == 0.0f)
	{
		float sum = 0.0f;
		uint32_t numFluidCells = 0;
		for (uint32_t i = 0; i < this->fluid_num_cells; i++)
		{
			if (this->cell_type[i] == FLUID_CELL)
			{
				sum += d[i];
				numFluidCells++;
			}
		}

		if (numFluidCells > 0)
			this->particle_rest_density = sum / numFluidCells;
	}
}

void FlipFluid::transfer_velocities(uint8_t toGrid, float flipRatio)
{
	const uint32_t n = this->fluid_num_y;
	const float h = this->h;
	const float h1 = this->h1;
	const float h2 = this->h2;
	if (toGrid)
	{
		this->prevU = this->u; // Copy contents of u into prevU
		this->prevV = this->v; // Copy contents of v into prevV

		std::fill(this->du.begin(), this->du.end(), 0.0f);
		std::fill(this->dv.begin(), this->dv.end(), 0.0f);
		std::fill(this->u.begin(), this->u.end(), 0.0f);
		std::fill(this->v.begin(), this->v.end(), 0.0f);

		for (uint32_t i = 0; i < this->fluid_num_cells; i++)
			this->cell_type[i] = this->s[i] == 0.0f ? SOLID_CELL : AIR_CELL;

		for (uint32_t i = 0; i < this->particle_nums; i++)
		{
			const float x = this->particle_pos[2 * i];
			const float y = this->particle_pos[2 * i + 1];
			const uint32_t xi = (uint32_t)clamp(floor(x * h1), 0, this->fluid_num_x - 1);
			const uint32_t yi = (uint32_t)clamp(floor(y * h1), 0, this->fluid_num_y - 1);
			const uint32_t cellNr = xi * n + yi;
			if (this->cell_type[cellNr] == AIR_CELL)
				this->cell_type[cellNr] = FLUID_CELL;
		}
	}

	for (uint32_t component = 0; component < 2; component++)
	{

		float dx = component == 0 ? 0.0f : h2;
		float dy = component == 0 ? h2 : 0.0f;

		std::vector<float> &f = component == 0 ? this->u : this->v;
		std::vector<float> &prevF = component == 0 ? this->prevU : this->prevV;
		std::vector<float> &d = component == 0 ? this->du : this->dv;

		for (uint32_t i = 0; i < this->particle_nums; i++)
		{
			float x = this->particle_pos[2 * i];
			float y = this->particle_pos[2 * i + 1];

			x = clamp(x, h, (this->fluid_num_x - 1) * h);
			y = clamp(y, h, (this->fluid_num_y - 1) * h);

			const uint32_t x0 = std::min((uint32_t)floor((x - dx) * h1), this->fluid_num_x - 2);
			const uint32_t x1 = std::min(x0 + 1, this->fluid_num_x - 2);
			const uint32_t y0 = std::min((uint32_t)floor((y - dy) * h1), this->fluid_num_y - 2);
			const uint32_t y1 = std::min(y0 + 1, this->fluid_num_y - 2);

			const float tx = ((x - dx) - x0 * h) * h1;
			const float ty = ((y - dy) - y0 * h) * h1;
			const float sx = 1.0f - tx;
			const float sy = 1.0f - ty;

			const float d0 = sx * sy;
			const float d1 = tx * sy;
			const float d2 = tx * ty;
			const float d3 = sx * ty;

			const uint32_t nr0 = x0 * n + y0;
			const uint32_t nr1 = x1 * n + y0;
			const uint32_t nr2 = x1 * n + y1;
			const uint32_t nr3 = x0 * n + y1;

			if (toGrid)
			{
				const float pv = this->particle_vel[2 * i + component];
				f[nr0] += pv * d0;
				d[nr0] += d0;
				f[nr1] += pv * d1;
				d[nr1] += d1;
				f[nr2] += pv * d2;
				d[nr2] += d2;
				f[nr3] += pv * d3;
				d[nr3] += d3;
			}
			else
			{
				uint32_t offset = component == 0 ? n : 1;
				const float valid0 = this->cell_type[nr0] != AIR_CELL || this->cell_type[nr0 - offset] != AIR_CELL ? 1.0f : 0.0f;
				const float valid1 = this->cell_type[nr1] != AIR_CELL || this->cell_type[nr1 - offset] != AIR_CELL ? 1.0f : 0.0f;
				const float valid2 = this->cell_type[nr2] != AIR_CELL || this->cell_type[nr2 - offset] != AIR_CELL ? 1.0f : 0.0f;
				const float valid3 = this->cell_type[nr3] != AIR_CELL || this->cell_type[nr3 - offset] != AIR_CELL ? 1.0f : 0.0f;

				const float v = this->particle_vel[2 * i + component];
				const float d = valid0 * d0 + valid1 * d1 + valid2 * d2 + valid3 * d3;

				if (d > 0.0f)
				{
					const float picV = (valid0 * d0 * f[nr0] + valid1 * d1 * f[nr1] + valid2 * d2 * f[nr2] + valid3 * d3 * f[nr3]) / d;
					const float corr = (valid0 * d0 * (f[nr0] - prevF[nr0]) + valid1 * d1 * (f[nr1] - prevF[nr1]) + valid2 * d2 * (f[nr2] - prevF[nr2]) + valid3 * d3 * (f[nr3] - prevF[nr3])) / d;
					const float flipV = v + corr;
					this->particle_vel[2 * i + component] = (1.0f - flipRatio) * picV + flipRatio * flipV;
				}
			}
		}

		if (toGrid)
		{
			for (uint32_t i = 0; i < f.size(); i++)
			{
				if (d[i] > 0.0f)
					f[i] /= d[i];
			}

			// restore solid cells

			for (uint32_t i = 0; i < this->fluid_num_x; i++)
			{
				for (uint32_t j = 0; j < this->fluid_num_y; j++)
				{
					const uint32_t solid = this->cell_type[i * n + j] == SOLID_CELL;
					if (solid || (i > 0 && this->cell_type[(i - 1) * n + j] == SOLID_CELL))
						this->u[i * n + j] = this->prevU[i * n + j];
					if (solid || (j > 0 && this->cell_type[i * n + j - 1] == SOLID_CELL))
						this->v[i * n + j] = this->prevV[i * n + j];
				}
			}
		}
	}
}

void FlipFluid::solve_incompressibility(uint32_t numIters, float dt, float overRelaxation, uint8_t compensateDrift)
{
	std::fill(this->p.begin(), this->p.end(), 0.0f); // Fill vector with 0.0f
	this->prevU = this->u;							 // Copy contents of `u` into `prevU`
	this->prevV = this->v;							 // Copy contents of `v` into `prevV`

	const uint32_t n = this->fluid_num_y;
	const float cp = this->density * this->h / dt;
	for (uint32_t iter = 0; iter < numIters; iter++)
	{
		for (uint32_t i = 1; i < this->fluid_num_x - 1; i++)
		{
			for (uint32_t j = 1; j < this->fluid_num_y - 1; j++)
			{
				if (this->cell_type[i * n + j] != FLUID_CELL)
					continue;

				const uint32_t center = i * n + j;
				const uint32_t left = (i - 1) * n + j;
				const uint32_t right = (i + 1) * n + j;
				const uint32_t bottom = i * n + j - 1;
				const uint32_t top = i * n + j + 1;

				const float sx0 = this->s[left];
				const float sx1 = this->s[right];
				const float sy0 = this->s[bottom];
				const float sy1 = this->s[top];
				const float s = sx0 + sx1 + sy0 + sy1;
				if (s == 0.0f)
					continue;

				float div = this->u[right] - this->u[center] +
							this->v[top] - this->v[center];

				if (this->particle_rest_density > 0.0f && compensateDrift)
				{
					const float k = 1.0f;
					const float compression = this->particle_density[i * n + j] - this->particle_rest_density;
					if (compression > 0.0f)
						div = div - k * compression;
				}

				const float p = - overRelaxation * div / s;
				this->p[center] += cp * p;

				this->u[center] -= sx0 * p;
				this->u[right] += sx1 * p;
				this->v[center] -= sy0 * p;
				this->v[top] += sy1 * p;
			}
		}
	}
}

void FlipFluid::updateParticleColor()
{
	const float h1 = this->h1;

	for (uint32_t i = 0; i < this->particle_nums; i++)
	{
		float s = 0.01f;

		this->particle_color[3 * i] = clamp(this->particle_color[3 * i] - s, 0.0f, 1.0f);
		this->particle_color[3 * i + 1] = clamp(this->particle_color[3 * i + 1] - s, 0.0f, 1.0f);
		this->particle_color[3 * i + 2] = clamp(this->particle_color[3 * i + 2] + s, 0.0f, 1.0f);

		float x = this->particle_pos[2 * i];
		float y = this->particle_pos[2 * i + 1];
		uint32_t xi = (uint32_t)clamp(floor(x * h1), 1, this->fluid_num_x - 1);
		uint32_t yi = (uint32_t)clamp(floor(y * h1), 1, this->fluid_num_y - 1);
		uint32_t cellNr = xi * this->fluid_num_y + yi;

		float d0 = this->particle_rest_density;

		if (d0 > 0.0f)
		{
			float relDensity = this->particle_density[cellNr] / d0;
			if (relDensity < 0.7f)
			{
				float s = 0.8f;
				this->particle_color[3 * i] = s;
				this->particle_color[3 * i + 1] = s;
				this->particle_color[3 * i + 2] = 1.0f;
			}
		}
	}
}

void FlipFluid::setSciColor(uint32_t cellNr, float val, float minVal, float maxVal)
{
	val = std::min(std::max(val, minVal), maxVal - 0.0001f);
	float d = maxVal - minVal;
	val = d == 0.0f ? 0.5f : (val - minVal) / d;
	float m = 0.25f;
	uint32_t num = (uint32_t)floor(val / m);
	float s = (val - num * m) / m;
	float r = 0.0f;
	float g = 0.0f;
	float b = 0.0f;

	switch (num)
	{
	case 0:
		r = 0.0f;
		g = s;
		b = 1.0f;
		break;
	case 1:
		r = 0.0f;
		g = 1.0f;
		b = 1.0f - s;
		break;
	case 2:
		r = s;
		g = 1.0f;
		b = 0.0f;
		break;
	case 3:
		r = 1.0f;
		g = 1.0f - s;
		b = 0.0f;
		break;
	}

	this->cell_color[3 * cellNr] = r;
	this->cell_color[3 * cellNr + 1] = g;
	this->cell_color[3 * cellNr + 2] = b;
}

void FlipFluid::updateCellColors()
{

	std::fill(this->cell_color.begin(), this->cell_color.end(), 0.0f);

	for (uint32_t i = 0; i < this->fluid_num_cells; i++)
	{
		if (this->cell_type[i] == SOLID_CELL)
		{
			this->cell_color[3 * i] = 0.5;
			this->cell_color[3 * i + 1] = 0.5;
			this->cell_color[3 * i + 2] = 0.5;
		}
		else if (this->cell_type[i] == FLUID_CELL)
		{
			float d = this->particle_density[i];
			if (this->particle_rest_density > 0.0)
				d /= this->particle_rest_density;
			this->setSciColor(i, d, 0.0, 2.0);
		}
	}
}

void FlipFluid::simulate(float dt, float gravity_x, float gravity_y, float flipRatio, uint32_t numPressureIters, uint32_t numParticleIters, float overRelaxation, uint8_t compensateDrift, uint8_t separateParticles, float obstacleX, float obstacleY, float obstacleRadius)
{
	uint32_t numSubSteps = 1;
	float sdt = dt / numSubSteps;

	for (uint32_t step = 0; step < numSubSteps; step++)
	{
		integrate_particles(sdt, gravity_x, gravity_y);
		if (separateParticles)
			push_particle_apart(numParticleIters);
		handle_particle_collisions();
		transfer_velocities(true, flipRatio);
		update_particle_density();
		solve_incompressibility(numPressureIters, sdt, overRelaxation, compensateDrift);
		transfer_velocities(false, flipRatio);
	}

	updateParticleColor();
	updateCellColors();
}
