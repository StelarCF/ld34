#version 150 core

uniform vec2 u_pos;
uniform vec2 u_window;
uniform float u_radius;
uniform float u_length;

in vec2 a_tex;
out vec2 v_coord;

void main() {
    v_coord = a_tex;
    float spread = 3.0 * a_tex.y + 1.0;
    float rad = -spread * (a_tex.x * 2.0 - 1.0);
    vec2 coord = rad * vec2(-1.0, 1.0) * vec2(0.5 * sqrt(2.0));
    coord *= 2.0 * (u_radius + 1.0);
    coord += 2.0 * u_pos * vec2(1.0, -1.0);
    coord += a_tex.y * u_length;
    gl_Position = vec4((coord.xy) / u_window, 0.0, 1.0);
}
