uniform sampler2D basetex;
uniform mat4 wvp;

uniform float sdf_min_t;
uniform float antialias_t;
uniform float sdf_boost_t;

#if defined(VERTEX_SHADER)
attribute vec2 position;
attribute vec2 texcoord;
attribute vec3 vertex_color;
#endif

varying vec2 frag_texcoord;
varying vec3 frag_color;

#if defined(VERTEX_SHADER)
void main()
{
	gl_Position = wvp * vec4(position.xy, 0.0, 1.0);
	frag_texcoord = texcoord;
	frag_color = vertex_color;
}
#endif

#if defined(FRAGMENT_SHADER)
void main()
{
	gl_FragColor = vec4(frag_color, 1.0);
	float d = texture2D(basetex, frag_texcoord).r;
	d = smoothstep(sdf_min_t - antialias_t, 0.5 + 0.5 * antialias_t, d);
	d *= sdf_boost_t;
	gl_FragColor.a *= d;
}
#endif
