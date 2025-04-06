#pragma once
#include <SFML/Window.hpp>
#include "conf.hpp"

extern bool g_is_applying_force;
extern float g_gravity_x;
extern float g_gravity_y;
extern bool g_draw_particles;
extern bool g_draw_cells;

void Window_Process_Events(sf::Window &window);