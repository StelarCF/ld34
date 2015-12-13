#version 150 core

uniform sampler2D u_penumbra;

in vec2 v_coord;
out vec4 o_color;

void main() {
    //softening
    float len = 1.0 - exp(-v_coord.y);
    float spread = abs(2.0 * v_coord.x - 1.0);
    float value = texture(u_penumbra, vec2(1.0 - len, 1.0 - spread)).r + 0.1;
    value += smoothstep(0.0, 1.0, v_coord.y) + 0.05;
    o_color = vec4(vec3(value), 1.0);
}
