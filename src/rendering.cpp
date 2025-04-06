#include "rendering.hpp"

static inline uint32_t circleIndex(uint32_t j, uint32_t edge_count)
{
    return 3 * edge_count * j;
}

static sf::Vector2f convertPosition(sf::Vector2f point)
{
    sf::Vector2f converted_point;

    converted_point.x = conf::canvas_position_f.x + point.x;
    converted_point.y = conf::canvas_position_f.y +
                        conf::canvas_size_f.y - point.y;

    return converted_point;
}

sf::Color getSciColor(float p, float minP, float maxP)
{
    float t = (p - minP) / (maxP - minP);
    t = std::clamp(t, 0.0f, 1.0f);

    sf::Color color;

    if (t < 0.25f)
    {
        // Red → Yellow (0.0 to 0.25)
        float ratio = t / 0.25f;
        color.r = 255;
        color.g = static_cast<sf::Uint8>(255 * ratio);
        color.b = 0;
    }
    else if (t < 0.50f)
    {
        // Yellow → Green (0.25 to 0.50)
        float ratio = (t - 0.25f) / 0.25f;
        color.r = static_cast<sf::Uint8>(255 * (1.0f - ratio));
        color.g = 255;
        color.b = 0;
    }
    else if (t < 0.75f)
    {
        // Green → Cyan (0.50 to 0.75)
        float ratio = (t - 0.50f) / 0.25f;
        color.r = 0;
        color.g = 255;
        color.b = static_cast<sf::Uint8>(255 * ratio);
    }
    else
    {
        // Cyan → Blue (0.75 to 1.0)
        float ratio = (t - 0.75f) / 0.25f;
        color.r = 0;
        color.g = static_cast<sf::Uint8>(255 * (1.0f - ratio));
        color.b = 255;
    }

    return color;
}

void drawCanvas(sf::RenderWindow &window)
{
    sf::RectangleShape canvas{conf::canvas_size_f};
    canvas.setPosition(conf::canvas_position_f);
    canvas.setFillColor(sf::Color(36, 40, 59));
    window.draw(canvas);
}

void drawParticle(sf::RenderWindow &window, FlipFluid *fluid)
{
    const int numTriangles = 12; // More = smoother circles
    float radius = fluid->particle_radius;

    // Each particle is made up of numTriangles triangles, 3 vertices per triangle
    sf::VertexArray vertices(sf::Triangles, fluid->particle_nums * numTriangles * 3);

    for (uint32_t i = 0; i < fluid->particle_nums; i++)
    {
        sf::Vector2f center = convertPosition({fluid->particle_pos[2 * i],
                                               fluid->particle_pos[2 * i + 1]});

        uint8_t r = static_cast<uint8_t>(255.0f * fluid->particle_color[3 * i]);
        uint8_t g = static_cast<uint8_t>(255.0f * fluid->particle_color[3 * i + 1]);
        uint8_t b = static_cast<uint8_t>(255.0f * fluid->particle_color[3 * i + 2]);

        sf::Color color(r, g, b);

        for (int j = 0; j < numTriangles; j++)
        {
            float angle1 = j * (2.0f * M_PI / numTriangles);
            float angle2 = (j + 1) * (2.0f * M_PI / numTriangles);

            sf::Vector2f p1 = center + sf::Vector2f(radius * std::cos(angle1), radius * std::sin(angle1));
            sf::Vector2f p2 = center + sf::Vector2f(radius * std::cos(angle2), radius * std::sin(angle2));

            int base = (i * numTriangles + j) * 3;

            vertices[base + 0].position = center;
            vertices[base + 0].color = color;

            vertices[base + 1].position = p1;
            vertices[base + 1].color = color;

            vertices[base + 2].position = p2;
            vertices[base + 2].color = color;
        }
    }

    window.draw(vertices);
}

void drawCell(sf::RenderWindow &window, FlipFluid *fluid)
{
    sf::VertexArray cells(sf::Quads, fluid->fluid_num_cells * 4);
    const uint32_t n = fluid->fluid_num_y;
    for (uint32_t i = 0; i < fluid->fluid_num_x; i++)
    {
        for (uint32_t j = 0; j < fluid->fluid_num_y; j++)
        {
            uint32_t index = (i * n + j);
            sf::Color color;
            switch (fluid->cell_type[index])
            {
            case FLUID_CELL:
                color = sf::Color(0, 255, 0);
                break;
            case SOLID_CELL:
                color = sf::Color(100, 100, 100);
                break;
            case AIR_CELL:
                color = sf::Color(0, 0, 0);
                break;
            default:
                break;
            }
            cells[4 * index]
                .position = convertPosition({i * fluid->h, j * fluid->h});
            cells[4 * index + 1].position = convertPosition({(i + 1) * fluid->h, j * fluid->h});
            cells[4 * index + 2].position = convertPosition({(i + 1) * fluid->h, (j + 1) * fluid->h});
            cells[4 * index + 3].position = convertPosition({i * fluid->h, (j + 1) * fluid->h});
            for (int k = 0; k < 4; k++)
                cells[4 * index + k].color = color;
        }
    }
    window.draw(cells);
}