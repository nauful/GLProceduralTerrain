uniform mat4 wvp;
uniform float dof_focal_distance;
uniform vec3 camera_pos;
uniform float elapsed_time;

uniform sampler2D tex_random;
uniform sampler2D tex_scheme;
uniform sampler2D tex_grass;
uniform sampler2D tex_lightmap;

uniform int use_composite_texturing;
uniform sampler2D tex_detail;
uniform sampler2D tex_moss1;
uniform sampler2D tex_moss2;
uniform sampler2D tex_rock;
uniform sampler2D tex_cliff;

uniform vec3 atmosphere_color_dense;
uniform vec3 atmosphere_color_apex;
uniform vec3 sun_dir;
uniform vec2 fog_distance;

uniform float gen_horz_scale;
uniform float gen_horz_minor_scale;
uniform float map_horz_scale;
uniform float map_vert_scale;
uniform vec2 node_base_xz;
uniform vec3 patch_color;

#if defined(VERTEX_SHADER)
attribute vec3 position;
attribute vec2 texcoord;
#endif

varying vec3 frag_pos;
varying vec2 frag_texcoord;

#if defined(PASS_RENDER)
varying float frag_fog;
#endif

vec3 cnoise_tex(vec2 pt)
{
	float i = floor(pt.x);
	float j = floor(pt.y);

	float u = pt.x - i;
	float v = pt.y - j;

    float du = 30.0 * u * u * (u * (u - 2.0) + 1.0);
    float dv = 30.0 * v * v * (v * (v - 2.0) + 1.0);

    u = u * u * u * (u * (u * 6.0 - 15.0) + 10.0);
    v = v * v * v * (v * (v * 6.0 - 15.0) + 10.0);

	float a = texture2D(tex_random, (1.0 / 128.0) * vec2(i, j)).r;
	float b = texture2D(tex_random, (1.0 / 128.0) * vec2(i + 1.0, j)).r;
	float c = texture2D(tex_random, (1.0 / 128.0) * vec2(i, j + 1.0)).r;
	float d = texture2D(tex_random, (1.0 / 128.0) * vec2(i + 1.0, j + 1.0)).r;

	float k0 = a;
	float k1 = b - a;
	float k2 = c - a;
	float k3 = a - b - c + d;
	return vec3(
		k0 + k1 * u + k2 * v + k3 * u * v,
		du * (k1 + k3 * v),
		dv * (k2 + k3 * u));
}

float fbm(vec2 c, float n)
{
	float f = 0.0;
	float w = 0.5;
	float persistence = 2.0;
	float gain = 0.5;

	vec2 d = vec2(0.0, 0.0);
	for (float i = 0.0; i < n; i += 1.001)
	{
		vec3 nr = cnoise_tex(c);
		f += w * nr.x;
		c *= persistence;
		w *= gain;
	}

	return f;
}

float terrain_feature_height(vec2 c, float n)
{
	float m = cnoise_tex(0.5 * c).x * 0.8 + 0.2;
	//m = max(0.2, m * m);

	float f = 0.0;
	float w = 0.5;
	float persistence = 2.0;
	float gain = 0.5;
	float erosion = 0.5;

	w *= m;

	vec2 d = vec2(0.0, 0.0);
	for (float i = 0.0; i < n; i += 1.001)
	{
		vec3 nr = cnoise_tex(c);
		d += nr.yz * erosion;
		f += w * nr.x / (1.0 + dot(d, d));
		c *= persistence;
		w *= gain;
	}

	return f;
}

float terrain_height(vec2 p)
{
	float h = terrain_feature_height(p, 5.0);
	h += (0.02 + max(h, 0.0) * 0.05) * (fbm(p * gen_horz_minor_scale, 5.0) - 0.5);
	return h;
}

float terrain_height_approx(vec2 p)
{
	return terrain_feature_height(p, 4.0);
}

float contrast(float value, float c)
{
	return clamp((value - 0.5) * c + 0.5, 0.0, 1.0);
}

vec3 contrast3(vec3 value, float c)
{
	return clamp((value - vec3(0.5, 0.5, 0.5)) * c + 0.5, 0.0, 1.0);
}

#if defined(VERTEX_SHADER)
void main()
{
	frag_pos = position;

#if defined(PASS_LIGHTMAP)
	frag_pos.xz += node_base_xz;
	frag_pos.y = terrain_height(frag_pos.xz * gen_horz_scale) * map_vert_scale;
	gl_Position = vec4(2.0 * (texcoord - vec2(0.5, 0.5)), 0.0, 1.0);
#elif defined(PASS_RENDER)
	frag_pos.xz += node_base_xz;
	frag_pos.y = terrain_height(frag_pos.xz * gen_horz_scale) * map_vert_scale;
	gl_Position = wvp * vec4(frag_pos, 1.0);
#else
	gl_Position = wvp * vec4(frag_pos, 1.0);
#endif

#if defined(PASS_RENDER)
	float fog_t = gl_Position.z;
	frag_fog = clamp((fog_t - fog_distance.x) / (fog_distance.y - fog_distance.x), 0.0, 1.0);
#endif

	frag_texcoord = texcoord;
}
#endif

#if defined(FRAGMENT_SHADER)
void main()
{
#if defined(PASS_LIGHTMAP)
	vec3 pos0 = vec3(frag_pos.x, 0.0, frag_pos.z);
	vec3 posx = pos0 + vec3(0.01, 0.0, 0.0);
	vec3 posz = pos0 + vec3(0.0, 0.0, 0.01);
	pos0.y = terrain_height(pos0.xz * gen_horz_scale) * map_vert_scale;
	posx.y = terrain_height(posx.xz * gen_horz_scale) * map_vert_scale;
	posz.y = terrain_height(posz.xz * gen_horz_scale) * map_vert_scale;
	vec3 n = cross(pos0 - posx, posz - pos0);
	// Prevent a face from ever facing perfectly away from the sun (fake ambient)
	n.y *= 1.7;
	n = normalize(n);

	pos0.y = terrain_height_approx(pos0.xz * gen_horz_scale) * map_vert_scale;
	posx.y = terrain_height_approx(posx.xz * gen_horz_scale) * map_vert_scale;
	posz.y = terrain_height_approx(posz.xz * gen_horz_scale) * map_vert_scale;
	vec3 n2 = cross(pos0 - posx, posz - pos0);
	n2 = normalize(n2);

	vec2 tex_param;
	tex_param.x = contrast(n2.y, 1.5);
	tex_param.y = clamp((1.25 + 0.5 * n2.y) * frag_pos.y / map_vert_scale, 0.0, 1.0);

	float d1 = max(0.0, dot(n, sun_dir));
	d1 = 0.3 + 0.7 * d1;

	gl_FragData[0].rgb = vec3(tex_param.x, tex_param.y, 0.0);
	gl_FragData[0].a = d1;
#elif defined(PASS_RENDER)
	vec4 lightmap = texture2D(tex_lightmap, frag_texcoord);
	gl_FragData[0] = vec4(1.0, 1.0, 1.0, 1.0);

	gl_FragData[0].rgb = texture2D(tex_scheme, lightmap.rg).rgb;
	gl_FragData[0].rgb *= lightmap.a;

	vec2 tex_pos = frag_pos.xz / 32.0;
	if (use_composite_texturing == 1)
	{
		gl_FragData[0].rgb *= texture2D(tex_detail, frag_pos.xz / 5.0).rgb + vec3(0.5, 0.5, 0.5);
		vec3 moss1 = texture2D(tex_moss1, tex_pos).rgb;
		vec3 rock = texture2D(tex_rock, tex_pos).rgb;
		vec3 cliff = texture2D(tex_cliff, tex_pos).rgb;

		//vec3 texbase = mix(moss1, moss2, contrast(tex_noise.z, 8.0));
		vec3 texbase = mix(moss1, cliff, contrast(lightmap.g, 4.0));
		texbase = mix(rock, texbase, contrast(lightmap.r, 4.0));
		gl_FragData[0].rgb *= texbase + vec3(0.5, 0.5, 0.5);
		gl_FragData[0].rgb *= 1.3;
	}

	gl_FragData[0].rgb = mix(gl_FragData[0].rgb, atmosphere_color_dense, min(1.0, frag_fog));
	gl_FragData[1] = vec4(clamp(dof_focal_distance / length(camera_pos - frag_pos), 0.0, 1.0), 0.0, 0.0, 1.0);
#elif defined(PASS_SKY)
	vec3 dir = normalize(frag_pos.xyz);
	vec3 res = mix(atmosphere_color_dense, atmosphere_color_apex, max(0.0, dir.y));

	float cp_t = 500.0 / dir.y;
	vec2 st = (frag_pos.xz + dir.xz * cp_t) * 0.003;
	st.x += elapsed_time * 0.05;
	float gain = 0.5;
	float r = 0.0;
	//float norm = 0.5;
	float clamp1 = -0.15;
	float clamp2 = 0.6;
	float a = 0.42;
	mat2 m = mat2(cos(a), sin(a), -sin(a), cos(a));
	for (float i = 0.0; i < 8.0; i += 1.0)
	{
		r -= gain * (2.0 * cnoise_tex(st).r - 0.85);
		st = (m * st) * 2.2;
		gain *= 0.6;
	}
	float v = clamp((r - clamp1) / (clamp2 - clamp1), 0.0, 1.0);
	float horizon_fade = clamp(dir.y * 10.0, 0.0, 1.0);
	v *= pow(horizon_fade, 3.0);

	gl_FragData[0] = vec4(mix(res, vec3(1.0, 1.0, 1.0), v), 1.0);

	float sun_t = max((dot(dir, sun_dir) - 0.997) / (1.0 - 0.997), 0.0);
	gl_FragData[0].rgb += vec3(1.0, 0.75, 0.5) * pow(sun_t, 8) * 16.0;

	gl_FragData[1] = vec4(0.0, 0.0, 0.0, 1.0);
#else
	gl_FragData[0] = vec4(1.0, 1.0, 1.0, 1.0);
#endif

#if defined(PASS_SKY) || defined(PASS_RENDER)
	// If gamma correction is needed
	// gl_FragData[0].r = pow(gl_FragData[0].r, 1.0);
	// gl_FragData[0].g = pow(gl_FragData[0].g, 1.0);
	// gl_FragData[0].b = pow(gl_FragData[0].b, 1.0);
#endif
}
#endif
