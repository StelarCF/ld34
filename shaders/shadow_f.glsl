#version 150 core

uniform float u_length;
uniform float u_radius;
uniform float u_radius_abs;

in vec2 v_coord;
out vec4 o_color;

const float s = 0.5 * sqrt(2.0);
const mat2 rotmat = mat2(s, s, -s, s);

void main() {
    vec2 coord = v_coord * rotmat;
    float len = coord.x / u_length;
    float spread = (abs(coord.y) - 2.0 * u_radius) / u_radius;
    /*vec2 pos = vec2(1.0, 1.0) * len;
    pos += vec2(1.0, -1.0) * spread;
    pos = vec2(pos.x, 1.0 - pos.y);
    float value = texture(u_penumbra, pos).r + 0.1;*/
    float value = smoothstep(-len, len, 0.01 * spread);// * u_radius_abs);
    value += smoothstep(0.3, 1.0, len) + 0.3;
    o_color = vec4(vec3(value), 1.0);
    //o_color = vec4(len, spread, 0.0, 1.0);
}
