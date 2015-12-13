#version 150 core

uniform sampler2D u_tex;

in vec2 v_coord;
out vec4 o_color;

const float ABB_SCALE[3] = float[](0.99, 0.98, 0.97);
const float ABB_INTENSITY[3] = float[](0.3, 0.2, 0.1);
const float BLUR_RADIUS = 0.01;

vec4 sample_tex(vec2 pos) {
    return texture(u_tex, pos * 0.5 + vec2(0.5, 0.5));
}

void main() {
    float dist = length(v_coord);
    //blur
    float blur_rad = BLUR_RADIUS * smoothstep(0.8 - 0.5, 0.8, dist) + 0.002;
    vec4 color = vec4(0.0);
    float total = 0.0;
    for(int dx = -5; dx <= 5; dx++) {
        for(int dy = -5; dy < 5; dy++) {
            float dx2 = blur_rad * float(dx) / 5.0;
            float dy2 = blur_rad * float(dy) / 5.0;
            vec2 coord = v_coord + vec2(dx2, dy2);
            vec4 tex = sample_tex(coord);
            //chromatic abberation
            tex.r += sample_tex(coord / ABB_SCALE[0]).r * ABB_INTENSITY[0];
            tex.r += sample_tex(coord / ABB_SCALE[1]).r * ABB_INTENSITY[1];
            tex.r += sample_tex(coord / ABB_SCALE[2]).r * ABB_INTENSITY[2];
            tex.b += sample_tex(coord * ABB_SCALE[0]).b * ABB_INTENSITY[0];
            tex.b += sample_tex(coord * ABB_SCALE[1]).b * ABB_INTENSITY[1];
            tex.b += sample_tex(coord * ABB_SCALE[2]).b * ABB_INTENSITY[2];
            float scale = smoothstep(0.0, 3.0, float(abs(dx) + abs(dy)));
                                     //length(vec2(float(dx), float(dy))));
            color += scale * tex;
            total += scale;
        }
    }
    color /= total;
    float vignette = smoothstep(1.2, 0.5, dist);
    o_color = color * vignette;
}
