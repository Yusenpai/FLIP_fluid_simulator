#pragma once

#include <SFML/Graphics.hpp>
#include "conf.hpp"
#include <algorithm>
#include "fluid/flip_fluid.hpp"

/// @brief Draw a canvas to run simulation in.
/// @param window Reference to RenderWindow object.
void drawCanvas(sf::RenderWindow &window);

/// @brief Draw particle
/// @param window Reference to RenderWindow object.
/// @param fluid Pointer to FlipFluid object
void drawParticle(sf::RenderWindow &window, FlipFluid *fluid);

void drawCell(sf::RenderWindow &window, FlipFluid *fluid);