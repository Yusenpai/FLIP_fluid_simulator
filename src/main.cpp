#include "main.hpp"

// Global Variables
bool g_is_applying_force = false;
uint32_t g_collision_count = 0;

FlipFluid *fluid;

struct AppContext
{
	sf::RenderWindow window;
	sf::Font font;
	sf::Text fps_text;
	sf::Clock clock;
};

struct Scene
{
	float gravity_x = 0.0f;
	float gravity_y = -9.81;
	float dt = 1.0f / 60.0f;
	float flipRatio = 0.9f;
	uint32_t num_pressure_iters = 100;
	uint32_t num_particle_iters = 2;
	float over_relaxation = 1.9f;
	uint8_t compensateDrift = 1;
	uint8_t separateParticles = 1;
	float obstacleX = 0.0f;
	float obstacleY = 0.0f;
	float obstacleRadius = 0.15f;
};

static Scene scene;
static AppContext app;

// Global Variables
float g_gravity_x = scene.gravity_x;
float g_gravity_y = scene.gravity_y;
bool g_draw_particles = true;
bool g_draw_cells = true;

// Function Declarations
static void Initialize();
static void Run_Simulation();
static void updatePhysics();
static void Render();
static void setupFluid();

static sf::Clock speed_clock;

int main()
{
	Initialize();
	while (app.window.isOpen())
	{
		Window_Process_Events(app.window);
		Run_Simulation();
		Render();
	}

	return 0;
}

static void Initialize()
{
	auto &window = app.window;
	auto &font = app.font;
	auto &fps_text = app.fps_text;

	window.create({conf::window_size.x, conf::window_size.y}, "Physics Simulator", sf::Style::Default);
	window.setFramerateLimit(conf::max_frame_rate);

	if (!font.loadFromFile("res/roboto.ttf"))
	{
		printf("Failed to load font!\n");
		exit(-1);
	}

	// FPS Text
	fps_text.setFont(font);
	fps_text.setCharacterSize(20);
	fps_text.setFillColor(sf::Color::White);
	fps_text.setPosition(10, 10);

	/* Setup Fluid */
	setupFluid();
}

static void Run_Simulation()
{
	auto &clock = app.clock;
	auto &window = app.window;
	auto &fps_text = app.fps_text;

	float elapsed_time = clock.restart().asSeconds();
	fps_text.setString("FPS: " + std::to_string(static_cast<int>(1.0f / elapsed_time)));
	updatePhysics();
}

static void updatePhysics()
{
	// speed_clock.restart();
	fluid->simulate(scene.dt, scene.gravity_x, scene.gravity_y, scene.flipRatio, scene.num_pressure_iters, scene.num_particle_iters, scene.over_relaxation, scene.compensateDrift, scene.separateParticles, scene.obstacleX, scene.obstacleY, scene.obstacleRadius);
	// float elapsed_time = speed_clock.restart().asMicroseconds();
	// printf("Elapsed time: %f\n", elapsed_time);
}

static void Render()
{
	auto &window = app.window;
	auto &font = app.font;
	auto &fps_text = app.fps_text;
	window.clear();
	drawCanvas(window);
	// if (g_draw_cells)
	// 	drawCell(window, fluid);
	// if (g_draw_particles)
	drawParticle(window, fluid);
	window.draw(fps_text);
	window.display();
}

static void setupFluid()
{
	scene.obstacleRadius = 0.15;
	scene.over_relaxation = 1.9;
	scene.dt = 1.0 / 60.0;
	scene.dt *= 10;
	scene.num_pressure_iters = 10;
	scene.num_particle_iters = 2;

	uint32_t res = 100;
	float tankHeight = 900.0;
	float tankWidth = 1600.0;
	float h = tankHeight / res;
	float density = 1000.0f;

	float relWaterHeight = 0.8f;
	float relWaterWidth = 0.6f;

	// dam break

	// compute number of particles

	// r = 0.3 * h = 0.3 * 62.5 = 18.75
	float r = 0.3f * h;
	float dx = 2.0f * r;
	float dy = std::sqrt(3.0f) / 2.0f * dx;

	uint32_t numX = (uint32_t)floor((relWaterWidth * tankWidth - 2.0f * h - 2.0 * r) / dx);
	uint32_t numY = (uint32_t)floor((relWaterHeight * tankHeight - 2.0f * h - 2.0 * r) / dy);
	uint32_t maxParticles = numX * numY;

	// create fluid
	fluid = new FlipFluid(density, tankWidth, tankHeight, h, r, maxParticles);

	// create particles
	fluid->particle_nums = numX * numY;
	uint32_t p = 0;
	for (uint32_t i = 0; i < numX; i++)
	{
		for (uint32_t j = 0; j < numY; j++)
		{
			uint32_t index = i * numY + j;
			fluid->particle_pos[2 * index] = h + r + dx * i + (j % 2 == 0 ? 0.0f : r);
			fluid->particle_pos[2 * index + 1] = h + r + dy * j;
			p++;
		}
	}

	// setup grid cells for tank

	uint32_t n = fluid->fluid_num_y;

	for (uint32_t i = 0; i < fluid->fluid_num_x; i++)
	{
		for (uint32_t j = 0; j < fluid->fluid_num_y; j++)
		{
			float s = 1.0f; // fluid
			if (i == 0 || i == fluid->fluid_num_x - 1 || j == 0 || j == fluid->fluid_num_y - 1)
			{
				s = 0.0f; // solid
			}
			fluid->s[i * n + j] = s;
		}
	}
}