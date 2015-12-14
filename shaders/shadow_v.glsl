#version 150 core

uniform vec2 u_pos;
uniform vec2 u_window;
uniform float u_radius;
uniform float u_length;

in vec2 a_tex;
out vec2 v_coord;

void main() {
    float spread = 0.08 * u_length * a_tex.y;
    float rad = -(a_tex.x * 2.0 - 1.0);
    rad += spread * sign(rad);
    vec2 coord = rad * vec2(-1.0, 1.0) * vec2(0.5 * sqrt(2.0));
    coord *= 2.0 * (u_radius + 1.0);
    coord += a_tex.y * u_length;
    v_coord = coord;
    coord += 2.0 * u_pos * vec2(1.0, -1.0);
    gl_Position = vec4(coord / u_window, 0.0, 1.0);
}
