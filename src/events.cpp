#include "events.hpp"

void Window_Process_Events(sf::Window &window)
{
	for (auto event = sf::Event{}; window.pollEvent(event);)
	{
		if (event.type == sf::Event::Closed)
		{
			window.close();
		}
		else if (event.type == sf::Event::KeyPressed)
		{
			switch (event.key.code)
			{
			case sf::Keyboard::Escape:
				window.close();
				break;

			case sf::Keyboard::F:
				g_is_applying_force = !g_is_applying_force;
				break;

			case sf::Keyboard::G:
				static uint8_t state = 0;
				switch (state)
				{
				case 0:
					g_gravity_x = 0.0f;
					g_gravity_y = -1000.0f;
					state++;
					break;
				case 1:
					g_gravity_x = 1000.0f;
					g_gravity_y = 0.0f;
					state++;
					break;
				case 2:
					g_gravity_x = 0.0f;
					g_gravity_y = 1000.0f;
					state++;
					break;
				case 3:
					g_gravity_x = -1000.0f;
					g_gravity_y = 0.0f;
					state = 0;
					break;
				default:
					break;
				}
				break;
			case sf::Keyboard::P:
				g_draw_particles = !g_draw_particles;
				break;
			case sf::Keyboard::C:
				g_draw_cells = !g_draw_cells;
				break;
			default:
				break;
			}
		}
	}
}
