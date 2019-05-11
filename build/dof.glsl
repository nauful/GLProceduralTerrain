uniform sampler2D basetex;
uniform sampler2D blurredtex;
uniform sampler2D z_tex;

uniform float tex_w;
uniform float tex_h;

#ifdef VERTEX_SHADER
attribute vec2 position;
#endif

varying vec2 frag_texcoord;

#ifdef VERTEX_SHADER
void main()
{
	gl_Position = vec4(position.xy, 0.0, 1.0);
	frag_texcoord = 0.5 * position.xy + 0.5;
}
#endif

#ifdef FRAGMENT_SHADER
float GaussianValue(float x, float dev)
{
	return exp(-(x * x) / (2.0 * dev * dev)) / sqrt(2.0 * 3.14159 * dev * dev);
}

void main()
{
	vec4 sum = vec4(0.0, 0.0, 0.0, 0.0);

	float center_depth = clamp(texture2D(z_tex, frag_texcoord).r, 0.0, 1.0);
	float blur_radius = 1.0 - center_depth;
	const float dev = 5.0;

#if defined(PASS_HORIZONTAL)
	float weight_neg_4 = GaussianValue(-4.0, dev) * (1.0 - clamp(texture2D(z_tex, frag_texcoord + vec2(-4.0 * blur_radius, 0.0) / tex_w).r, 0.0, 1.0));
	float weight_neg_3 = GaussianValue(-3.0, dev) * (1.0 - clamp(texture2D(z_tex, frag_texcoord + vec2(-3.0 * blur_radius, 0.0) / tex_w).r, 0.0, 1.0));
	float weight_neg_2 = GaussianValue(-2.0, dev) * (1.0 - clamp(texture2D(z_tex, frag_texcoord + vec2(-2.0 * blur_radius, 0.0) / tex_w).r, 0.0, 1.0));
	float weight_neg_1 = GaussianValue(-1.0, dev) * (1.0 - clamp(texture2D(z_tex, frag_texcoord + vec2(-1.0 * blur_radius, 0.0) / tex_w).r, 0.0, 1.0));
	float weight_center = GaussianValue(0.0, dev);
	float weight_pos_1 = GaussianValue(1.0, dev) * (1.0 - clamp(texture2D(z_tex, frag_texcoord + vec2(1.0 * blur_radius, 0.0) / tex_w).r, 0.0, 1.0));
	float weight_pos_2 = GaussianValue(2.0, dev) * (1.0 - clamp(texture2D(z_tex, frag_texcoord + vec2(2.0 * blur_radius, 0.0) / tex_w).r, 0.0, 1.0));
	float weight_pos_3 = GaussianValue(3.0, dev) * (1.0 - clamp(texture2D(z_tex, frag_texcoord + vec2(3.0 * blur_radius, 0.0) / tex_w).r, 0.0, 1.0));
	float weight_pos_4 = GaussianValue(4.0, dev) * (1.0 - clamp(texture2D(z_tex, frag_texcoord + vec2(4.0 * blur_radius, 0.0) / tex_w).r, 0.0, 1.0));

	sum.rgb += texture2D(basetex, frag_texcoord + vec2(-4.0 * blur_radius, 0.0) / tex_w).rgb * weight_neg_4;
	sum.rgb += texture2D(basetex, frag_texcoord + vec2(-3.0 * blur_radius, 0.0) / tex_w).rgb * weight_neg_3;
	sum.rgb += texture2D(basetex, frag_texcoord + vec2(-2.0 * blur_radius, 0.0) / tex_w).rgb * weight_neg_2;
	sum.rgb += texture2D(basetex, frag_texcoord + vec2(-1.0 * blur_radius, 0.0) / tex_w).rgb * weight_neg_1;
	sum.rgb += texture2D(basetex, frag_texcoord + vec2(0.0, 0.0) / tex_w).rgb * weight_center;
	sum.rgb += texture2D(basetex, frag_texcoord + vec2(1.0 * blur_radius, 0.0) / tex_w).rgb * weight_pos_1;
	sum.rgb += texture2D(basetex, frag_texcoord + vec2(2.0 * blur_radius, 0.0) / tex_w).rgb * weight_pos_2;
	sum.rgb += texture2D(basetex, frag_texcoord + vec2(3.0 * blur_radius, 0.0) / tex_w).rgb * weight_pos_3;
	sum.rgb += texture2D(basetex, frag_texcoord + vec2(4.0 * blur_radius, 0.0) / tex_w).rgb * weight_pos_4;

	sum.a += weight_neg_4;
	sum.a += weight_neg_3;
	sum.a += weight_neg_2;
	sum.a += weight_neg_1;
	sum.a += weight_center;
	sum.a += weight_pos_1;
	sum.a += weight_pos_2;
	sum.a += weight_pos_3;
	sum.a += weight_pos_4;
#endif

#if defined(PASS_VERTICAL)
	float weight_neg_4 = GaussianValue(-4.0, dev) * (1.0 - clamp(texture2D(z_tex, frag_texcoord + vec2(0.0, -4.0 * blur_radius) / tex_w).r, 0.0, 1.0));
	float weight_neg_3 = GaussianValue(-3.0, dev) * (1.0 - clamp(texture2D(z_tex, frag_texcoord + vec2(0.0, -3.0 * blur_radius) / tex_w).r, 0.0, 1.0));
	float weight_neg_2 = GaussianValue(-2.0, dev) * (1.0 - clamp(texture2D(z_tex, frag_texcoord + vec2(0.0, -2.0 * blur_radius) / tex_w).r, 0.0, 1.0));
	float weight_neg_1 = GaussianValue(-1.0, dev) * (1.0 - clamp(texture2D(z_tex, frag_texcoord + vec2(0.0, -1.0 * blur_radius) / tex_w).r, 0.0, 1.0));
	float weight_center = GaussianValue(0.0, dev);
	float weight_pos_1 = GaussianValue(1.0, dev) * (1.0 - clamp(texture2D(z_tex, frag_texcoord + vec2(0.0, 1.0 * blur_radius) / tex_w).r, 0.0, 1.0));
	float weight_pos_2 = GaussianValue(2.0, dev) * (1.0 - clamp(texture2D(z_tex, frag_texcoord + vec2(0.0, 2.0 * blur_radius) / tex_w).r, 0.0, 1.0));
	float weight_pos_3 = GaussianValue(3.0, dev) * (1.0 - clamp(texture2D(z_tex, frag_texcoord + vec2(0.0, 3.0 * blur_radius) / tex_w).r, 0.0, 1.0));
	float weight_pos_4 = GaussianValue(4.0, dev) * (1.0 - clamp(texture2D(z_tex, frag_texcoord + vec2(0.0, 4.0 * blur_radius) / tex_w).r, 0.0, 1.0));

	sum.rgb += texture2D(basetex, frag_texcoord + vec2(0.0, -4.0 * blur_radius) / tex_w).rgb * weight_neg_4;
	sum.rgb += texture2D(basetex, frag_texcoord + vec2(0.0, -3.0 * blur_radius) / tex_w).rgb * weight_neg_3;
	sum.rgb += texture2D(basetex, frag_texcoord + vec2(0.0, -2.0 * blur_radius) / tex_w).rgb * weight_neg_2;
	sum.rgb += texture2D(basetex, frag_texcoord + vec2(0.0, -1.0 * blur_radius) / tex_w).rgb * weight_neg_1;
	sum.rgb += texture2D(basetex, frag_texcoord + vec2(0.0, 0.0) / tex_w).rgb * weight_center;
	sum.rgb += texture2D(basetex, frag_texcoord + vec2(0.0, 1.0 * blur_radius) / tex_w).rgb * weight_pos_1;
	sum.rgb += texture2D(basetex, frag_texcoord + vec2(0.0, 2.0 * blur_radius) / tex_w).rgb * weight_pos_2;
	sum.rgb += texture2D(basetex, frag_texcoord + vec2(0.0, 3.0 * blur_radius) / tex_w).rgb * weight_pos_3;
	sum.rgb += texture2D(basetex, frag_texcoord + vec2(0.0, 4.0 * blur_radius) / tex_w).rgb * weight_pos_4;

	sum.a += weight_neg_4;
	sum.a += weight_neg_3;
	sum.a += weight_neg_2;
	sum.a += weight_neg_1;
	sum.a += weight_center;
	sum.a += weight_pos_1;
	sum.a += weight_pos_2;
	sum.a += weight_pos_3;
	sum.a += weight_pos_4;
#endif

#if defined(PASS_COMBINE)
	float combine_weight = pow(center_depth, 3.0);
	sum.rgb += texture2D(basetex, frag_texcoord).rgb * combine_weight;
	sum.rgb += texture2D(blurredtex, frag_texcoord).rgb * (1.0 - combine_weight);
	sum.a = 1.0;
	//sum.rgb = texture2D(blurredtex, frag_texcoord).rgb;
	//sum.rgb = combine_weight;
#endif

    gl_FragColor = vec4(sum.rgb / sum.a, 1.0);
}
#endif
