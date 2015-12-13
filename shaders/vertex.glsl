#version 150 core

in vec3 a_pos;
out vec2 v_coord;

void main() {
    v_coord = a_pos.xy;
    gl_Position = vec4(a_pos, 1.0);
}
