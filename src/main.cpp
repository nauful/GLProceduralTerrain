#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <assert.h>

#include "../lib/math_util.h"
#include "../lib/gl_util.h"

int screen_w = 1280;
int screen_h = 720;

const float camera_speed = 0.08f;
const float camera_speed_fast = 0.4f;

Vector3 sun_dir(3, 5, 2);

const int dof_downsample_factor = 2;
const int dof_blur_passes = 5;
const float dof_focal_distance = 200;

bool use_composite_texturing = true;
bool use_dof = true;
bool use_fast_camera = false;
bool render_wireframe = false;

const int map_size = 4 << 20;
const float map_horz_scale = 4.0f;
const float map_horz_gen_scale = 0.01f;
const float map_horz_gen_minor_scale = 10.0f;
const float map_vert_scale = 600.0f;
const float map_visible_dist = 3000.0f;
const float patch_tex_scale = 4;
const float patch_tex_max = 1024;
const float actual_visible_dist = map_visible_dist * map_horz_scale;

const Vector3 atmosphere_color_apex(0.15f, 0.24f, 0.44f);
const Vector3 atmosphere_color_dense(0.47f, 0.64f, 0.95f);

const float frustum_far_plane = actual_visible_dist * 1.25f;
const float frustum_near_plane = 3.0f;

const float fog_start = 500.0f;
const float fog_end = 8000.0f;

// Used when rendering to multiple rendertargets simultaneously
GLenum draw_buffers[] = { GL_COLOR_ATTACHMENT0, GL_COLOR_ATTACHMENT1, GL_COLOR_ATTACHMENT2, GL_COLOR_ATTACHMENT3 };

GLFont font;

GLuint shader_terrain = 0;
GLuint shader_terrain_generate_lightmap = 0;
GLuint shader_sky = 0;

GLuint shader_dof_blur_horz = 0;
GLuint shader_dof_blur_vert = 0;
GLuint shader_dof_combine = 0;

GLuint tex_random = 0;
GLuint tex_scheme = 0;
GLuint tex_terrain[5] = { 0 };

GLuint dof_fbo[3] = { 0 };
// Two RTs per buffer (image, z from camera)
// Depth is a simple depth buffer
GLuint dof_rt[3] = { 0 };
GLuint dof_rt_z[3] = { 0 };
GLuint dof_rt_depth[3] = { 0 };
GLuint cur_fbo = 0;

int fps = 0;
int frame = 0;
int time_fps = 0;
int prev_time = 0;

Vector3 camera_pos(180, 147, -81);
// Forward, strafe
Vector3 camera_motion(0, 0, 0);

float camera_vangle = 75;
float camera_hangle = 0;
int glut_left_mouse_state = GLUT_UP;
int glut_mouse_x = 0;
int glut_mouse_y = 0;

struct Vertex
{
	float x, y, z;
	float u, v;
};

// Used to find neighbouring nodes
enum NodeLocation { NL_TL = 0, NL_TR = 1, NL_BL = 2, NL_BR = 3 };
enum NodeNeighbour { NN_T = 0, NN_R = 1, NN_B = 2, NN_L = 3 };

struct Node
{
	int level;
	// Center LOD + 4 adjacent
	int lod[5];

	int pos_x, pos_z;
	int size_x, size_z;
	bool distance_cull;

	bool has_patch;
	int patch_lod[5];
	int num_patch_indices;
	GLuint patch_vb, patch_ib, patch_lightmap;

	Node* parent;
	Node* children[4];
	bool has_children;
};

Node* terrain_root = 0;

// Frustum for frustum culling
class Frustum
{
private:
	Plane planes[6];
	Vector3 points[8];
	Vector3 nearpoint;

public:
	Frustum(const Matrix& viewproj)
	{
		points[0] = Vector3(-1, -1, 0);
		points[1] = Vector3(1, -1, 0);
		points[2] = Vector3(-1, 1, 0);
		points[3] = Vector3(1, 1, 0);
		points[4] = Vector3(-1, -1, 1);
		points[5] = Vector3(1, -1, 1);
		points[6] = Vector3(-1, 1, 1);
		points[7] = Vector3(1, 1, 1);
		nearpoint = Vector3(0.5f, 0.5f, 0);

		Matrix inv = viewproj;
		inv.Invert();
		for (int i = 0; i < 8; i++) { points[i] = inv.TransformCoord(points[i]); }
		nearpoint = inv.TransformCoord(nearpoint);

		planes[0] = Plane(points[0], points[1], points[2]);
		planes[1] = Plane(points[6], points[7], points[5]);
		planes[2] = Plane(points[2], points[6], points[4]);
		planes[3] = Plane(points[7], points[3], points[5]);
		planes[4] = Plane(points[2], points[3], points[6]);
		planes[5] = Plane(points[1], points[0], points[5]);
	}

	bool TestExtent(const Vector3& pt, const Vector3& extents) const
	{
		for (int i = 0; i < 6; i++)
		{
			float d = planes[i].x * pt.x + planes[i].y * pt.y + planes[i].z * pt.z;
			float r = extents.x * fabsf(planes[i].x) + extents.y * fabsf(planes[i].y) + extents.z * fabsf(planes[i].z);
			if (d + r < -planes[i].d) { return false; }
		}
	
		return true;
	}
	
	bool TestAABB(const Vector3& bbmin, const Vector3& bbmax) const
	{
		return TestExtent(0.5f * (bbmin + bbmax), 0.5f * (bbmax - bbmin));
	}
};

Node* node_create(Node* parent, int pos_x, int pos_z, int size_x, int size_z)
{
	Node* n = new Node();
	n->parent = parent;

	memset(n->children, 0, sizeof(Node*) * 4);
	memset(n->lod, 0, sizeof(int) * 5);

	n->pos_x = pos_x;
	n->pos_z = pos_z;
	n->size_x = size_x;
	n->size_z = size_z;
	n->distance_cull = false;

	n->level = parent ? parent->level + 1 : 0;

	n->has_patch = false;
	memset(n->patch_lod, 0, sizeof(int) * 5);
	n->num_patch_indices = 0;
	n->patch_vb = 0;
	n->patch_ib = 0;
	n->patch_lightmap = 0;
	n->has_children = false;

	return n;
}

void release_patch(Node* n)
{
	if (n->has_patch)
	{
		n->has_patch = false;
		n->num_patch_indices = 0;
		glDeleteBuffers(1, &n->patch_vb);
		glDeleteBuffers(1, &n->patch_ib);
		glDeleteTextures(1, &n->patch_lightmap);
	}
}

void node_free(Node* n)
{
	if (n->has_children)
	{
		for (int i = 0; i < 4; i++) { node_free(n->children[i]); }
	}

	if (n->has_patch)
	{
		release_patch(n);
	}
	delete n;
}

Vector3 project_to_line(const Vector3& v1, const Vector3& v2, const Vector3& pt)
{
	Vector3 pa = pt - v1;
	Vector3 pb = v2 - v1;
	float t = pa.Dot(pb) / pb.Dot(pb);
	if (t < 0.0f) { return v1; }
	else if (t > 1.0f) { return v2; }

	return Vector3::Lerp(v1, v2, t);
}

// LOD level is calculated as a combination of logarithmic and linear distance
// Coefficients were determined through trial and error
// LOD levels are in the range [0, 15]
int lod_level(float dist)
{
	assert(dist >= 0);

	float lod_log = logf(max(1, 0.05f * dist - 5.0f));
	lod_log /= logf(3);
	float lod_linear = 0.005f * dist;

	float lod = lod_log + 0.05f * (lod_linear - lod_log);
	return max(0, (int)min(15, floorf(lod)));
}

// Closest point on a patch from a source point
Vector3 closest_point(const Vector3& source, int minx, int minz, int maxx, int maxz)
{
	if (source.x >= minx && source.z >= minz &&
		source.x <= maxx && source.z <= maxz) { return source; }

	Vector3 line_seg[] =
	{
		Vector3(minx, 0, minz), Vector3(maxx, 0, minz),
		Vector3(minx, 0, maxz), Vector3(maxx, 0, maxz),
		Vector3(minx, 0, minz), Vector3(minx, 0, maxz),
		Vector3(maxx, 0, minz), Vector3(maxx, 0, maxz),
	};

	Vector3 closest_pt;
	float best_dist_sqr = 0;
	for (int i = 0; i < 4; i++)
	{
		Vector3 proj_pt = project_to_line(line_seg[i * 2], line_seg[i * 2 + 1], source);
		if (!i || (proj_pt - source).LengthSquared() < best_dist_sqr)
		{
			closest_pt = proj_pt;
			best_dist_sqr = (proj_pt - source).LengthSquared();
		}
	}
	return closest_pt;
}

// Furthest point on a patch from a point
Vector3 furthest_point(const Vector3& source, int minx, int minz, int maxx, int maxz)
{
	Vector3 pts[] =
	{
		Vector3(minx, 0, minz),
		Vector3(minx, 0, maxz),
		Vector3(maxx, 0, minz),
		Vector3(maxx, 0, maxz),
	};

	Vector3 furthest;
	float best_dist_sqr = 0;
	for (int i = 0; i < 4; i++)
	{
		if (!i || (pts[i] - source).LengthSquared() > best_dist_sqr)
		{
			furthest = pts[i];
			best_dist_sqr = (pts[i] - source).LengthSquared();
		}
	}
	return furthest;
}

void split(Node* n)
{
	assert(!n->has_children);

	n->children[NL_TL] = node_create(n, n->pos_x, n->pos_z, n->size_x / 2, n->size_z / 2);
	n->children[NL_TR] = node_create(n, n->pos_x + n->size_x / 2, n->pos_z, n->size_x / 2, n->size_z / 2);
	n->children[NL_BL] = node_create(n, n->pos_x, n->pos_z + n->size_z / 2, n->size_x / 2, n->size_z / 2);
	n->children[NL_BR] = node_create(n, n->pos_x + n->size_x / 2, n->pos_z + n->size_z / 2, n->size_x / 2, n->size_z / 2);
	n->has_children = true;
}

void merge(Node* n)
{
	if (n->has_children)
	{
		for (int i = 0; i < 4; i++) { node_free(n->children[i]); }
		n->has_children = false;
	}

	if (n->has_patch) { release_patch(n); }
}

int calculate_edge_lod(Node* n, Node* adj)
{
	// Adjacent has children which will connect to this
	if (adj->has_children) { return n->lod[0]; }

	// Smaller edges stitch to larger patches
	if (adj->level < n->level) { return adj->lod[0]; }

	// Edge patches must be from finer lod to coarser/equal lod
	return max(adj->lod[0], n->lod[0]);
}

Node* find_adj(Node* n, NodeNeighbour nn)
{
	if (n->parent)
	{
		int this_loc = -1;
		for (int i = 0; i < 4; i++)
		{
			if (n->parent->children[i] == n) { this_loc = i; }
		}
		assert(this_loc > -1);

		if (this_loc == NL_TL && nn == NN_R) { return n->parent->children[NL_TR]; }
		if (this_loc == NL_BL && nn == NN_R) { return n->parent->children[NL_BR]; }
		if (this_loc == NL_TR && nn == NN_L) { return n->parent->children[NL_TL]; }
		if (this_loc == NL_BR && nn == NN_L) { return n->parent->children[NL_BL]; }

		if (this_loc == NL_TL && nn == NN_B) { return n->parent->children[NL_BL]; }
		if (this_loc == NL_TR && nn == NN_B) { return n->parent->children[NL_BR]; }
		if (this_loc == NL_BL && nn == NN_T) { return n->parent->children[NL_TL]; }
		if (this_loc == NL_BR && nn == NN_T) { return n->parent->children[NL_TR]; }

		Node* candidate = find_adj(n->parent, nn);
		if (!candidate) { return 0; }

		if (candidate->level + 1 == n->level && candidate->has_children)
		{
			if (nn == NN_T)
			{
				if (this_loc == NL_TL) { candidate = candidate->children[NL_BL]; }
				else if (this_loc == NL_TR) { candidate = candidate->children[NL_BR]; }
				else { assert(false); }
			}
			else if (nn == NN_B)
			{
				if (this_loc == NL_BL) { candidate = candidate->children[NL_TL]; }
				else if (this_loc == NL_BR) { candidate = candidate->children[NL_TR]; }
				else { assert(false); }
			}
			else if (nn == NN_R)
			{
				if (this_loc == NL_TR) { candidate = candidate->children[NL_TL]; }
				else if (this_loc == NL_BR) { candidate = candidate->children[NL_BL]; }
				else { assert(false); }
			}
			else if (nn == NN_L)
			{
				if (this_loc == NL_TL) { candidate = candidate->children[NL_TR]; }
				else if (this_loc == NL_BL) { candidate = candidate->children[NL_BR]; }
				else { assert(false); }
			}
		}
		return candidate;
	}

	return 0;
}

void update_edge_lod(Node* n)
{
	Node* adj_t = find_adj(n, NN_T);
	Node* adj_b = find_adj(n, NN_B);
	Node* adj_r = find_adj(n, NN_R);
	Node* adj_l = find_adj(n, NN_L);
	n->lod[1 + NN_T] = adj_t && !adj_t->distance_cull ? calculate_edge_lod(n, adj_t) : n->lod[0];
	n->lod[1 + NN_R] = adj_r && !adj_r->distance_cull ? calculate_edge_lod(n, adj_r) : n->lod[0];
	n->lod[1 + NN_B] = adj_b && !adj_b->distance_cull ? calculate_edge_lod(n, adj_b) : n->lod[0];
	n->lod[1 + NN_L] = adj_l && !adj_l->distance_cull ? calculate_edge_lod(n, adj_l) : n->lod[0];

	if (n->has_children)
	{
		for (int i = 0; i < 4; i++) { update_edge_lod(n->children[i]); }
	}
}

int update_lod(Vector3 cam_pos, Node* n)
{
	Vector3 planar_pos(cam_pos.x, 0, cam_pos.z);
	float dnear = (planar_pos - closest_point(planar_pos, n->pos_x, n->pos_z, n->pos_x + n->size_x, n->pos_z + n->size_z)).Length();
	float dfar = (planar_pos - furthest_point(planar_pos, n->pos_x, n->pos_z, n->pos_x + n->size_x, n->pos_z + n->size_z)).Length();	
	int lodnear = lod_level(dnear);
	int lodfar = lod_level(dfar);

	n->lod[0] = lodnear;
	if (dnear > map_visible_dist)
	{
		n->distance_cull = true;
		merge(n);
		return 1;
	}
	else { n->distance_cull = false; }

	int dim = n->size_x / max(1.0f, powf(2, n->lod[0]));
	int tex_res = dim * patch_tex_scale;

	// Patch generation rules:
	// Surface texture too large -> split
	// Patch too small to be worth a draw call -> merge
	// LOD split
	// LOD merge
	if (!n->has_children && tex_res > patch_tex_max)
	{
		split(n);
	}
	else if (dim < 8 && n->parent)
	{
		merge(n->parent);
		return 0;
	}
	else if ((lodfar - lodnear > 0 && dim >= 128) || dim > 1024)
	{
		if (!n->has_children) { split(n); }
	}
	else if (lodfar == lodnear && n->has_children) // && tex_res <= patch_max_tex_res)
	{
		merge(n);
	}

	int num_nodes = 1;
	for (int i = 0; i < 4; i++)
	{
		// Merged
		if (!n->has_children) { break; }
		num_nodes += update_lod(cam_pos, n->children[i]);
	}
	return num_nodes;
}

void generate_mesh(Node* n)
{
	assert(!n->has_children);
	assert(!n->distance_cull);
	if (n->has_patch) { release_patch(n); }

	int lod = n->lod[0];
	int lod_top = n->lod[1 + NN_T];
	int lod_right = n->lod[1 + NN_R];
	int lod_bottom = n->lod[1 + NN_B];
	int lod_left = n->lod[1 + NN_L];

	// Edge patches must be from finer lod to coarser/equal lod
	if (lod > lod_top) { lod = lod_top; }
	if (lod > lod_right) { lod = lod_right; }
	if (lod > lod_bottom) { lod = lod_bottom; }
	if (lod > lod_left) { lod = lod_left; }

	int patch_size = n->size_x;
	// Only square patches are permitted
	assert(n->size_x == n->size_z);

	int patch_lod_step = max(1, powf(2, lod));
	int patch_top_step = max(1, powf(2, lod_top)) / patch_lod_step;
	int patch_left_step = max(1, powf(2, lod_left)) / patch_lod_step;
	int patch_right_step = max(1, powf(2, lod_right)) / patch_lod_step;
	int patch_bottom_step = max(1, powf(2, lod_bottom)) / patch_lod_step;
	
	// Can't skip the entire patch
	if (patch_lod_step >= patch_size) { patch_lod_step = patch_size; }

	int patch_lod_size = patch_size / patch_lod_step;
	assert(patch_top_step <= patch_lod_size);
	assert(patch_left_step <= patch_lod_size);
	assert(patch_right_step <= patch_lod_size);
	assert(patch_bottom_step <= patch_lod_size);

	for (int i = 0; i < 5; i++) { n->patch_lod[i] = n->lod[i]; }

	int vertex_pitch = 1 + patch_lod_size;

	Vertex* vertex_data = new Vertex[vertex_pitch * vertex_pitch * 4];
	GLuint* index_data = new GLuint[(vertex_pitch * vertex_pitch * 6 + patch_size * 64 * 4)];
	for (int z = 0; z <= patch_lod_size; z++)
	{
		for (int x = 0; x <= patch_lod_size; x++)
		{
			Vertex v;
			v.x = x * patch_lod_step;
			v.z = z * patch_lod_step;
			v.y = 0;
			v.u = (float)x / (float)patch_lod_size;
			v.v = (float)z / (float)patch_lod_size;
			vertex_data[z * vertex_pitch + x] = v;
		}
	}

#define I(x, z) ((z) * vertex_pitch + (x))
	int num_indices = 0;
	for (int z = 1; z < patch_lod_size - 1; z++)
	{
		for (int x = 1; x < patch_lod_size - 1; x++)
		{
			index_data[num_indices++] = I(x, z);
			index_data[num_indices++] = I(x + 1, z);
			index_data[num_indices++] = I(x, z + 1);

			index_data[num_indices++] = I(x + 1, z);
			index_data[num_indices++] = I(x + 1, z + 1);
			index_data[num_indices++] = I(x, z + 1);
		}
	}

	if (lod_top == lod)
	{
		index_data[num_indices++] = I(0, 0);
		index_data[num_indices++] = I(1, 0);
		index_data[num_indices++] = I(1, 1);

		index_data[num_indices++] = I(patch_lod_size - 1, 0);
		index_data[num_indices++] = I(patch_lod_size - 1, 1);
		index_data[num_indices++] = I(patch_lod_size, 0);
			
		for (int x = 1; x < patch_lod_size - 1; x++)
		{
			index_data[num_indices++] = I(x, 0);
			index_data[num_indices++] = I(x + 1, 0);
			index_data[num_indices++] = I(x, 1);

			index_data[num_indices++] = I(x + 1, 0);
			index_data[num_indices++] = I(x + 1, 1);
			index_data[num_indices++] = I(x, 1);
		}
	}
	else
	{
		for (int x = 0; x < patch_lod_size; x += patch_top_step)
		{
			index_data[num_indices++] = I(x, 0);
			index_data[num_indices++] = I(x + patch_top_step / 2, 1);
			index_data[num_indices++] = I(x + patch_top_step, 0);

			for (int j = 0; j < patch_top_step / 2; j++)
			{
				if (x == 0 && j == 0) { continue; }
				index_data[num_indices++] = I(x, 0);
				index_data[num_indices++] = I(x + j + 1, 1);
				index_data[num_indices++] = I(x + j, 1);
			}
			for (int j = patch_top_step / 2; j < patch_top_step; j++)
			{
				if (x == patch_lod_size - patch_top_step && j == patch_top_step - 1) { continue; }
				index_data[num_indices++] = I(x + patch_top_step, 0);
				index_data[num_indices++] = I(x + j + 1, 1);
				index_data[num_indices++] = I(x + j, 1);
			}
		}
	}

	if (lod_bottom == lod)
	{
		index_data[num_indices++] = I(0, patch_lod_size);
		index_data[num_indices++] = I(1, patch_lod_size - 1);
		index_data[num_indices++] = I(1, patch_lod_size);

		index_data[num_indices++] = I(patch_lod_size - 1, patch_lod_size);
		index_data[num_indices++] = I(patch_lod_size - 1, patch_lod_size - 1);
		index_data[num_indices++] = I(patch_lod_size, patch_lod_size);
			
		for (int x = 1; x < patch_lod_size - 1; x++)
		{
			index_data[num_indices++] = I(x, patch_lod_size - 1);
			index_data[num_indices++] = I(x + 1, patch_lod_size - 1);
			index_data[num_indices++] = I(x, patch_lod_size);

			index_data[num_indices++] = I(x + 1, patch_lod_size - 1);
			index_data[num_indices++] = I(x + 1, patch_lod_size);
			index_data[num_indices++] = I(x, patch_lod_size);
		}
	}
	else
	{
		for (int x = 0; x < patch_lod_size; x += patch_bottom_step)
		{
			index_data[num_indices++] = I(x, patch_lod_size);
			index_data[num_indices++] = I(x + patch_bottom_step / 2, patch_lod_size - 1);
			index_data[num_indices++] = I(x + patch_bottom_step, patch_lod_size);

			for (int j = 0; j < patch_bottom_step / 2; j++)
			{
				if (x == 0 && j == 0) { continue; }
				index_data[num_indices++] = I(x, patch_lod_size);
				index_data[num_indices++] = I(x + j + 1, patch_lod_size - 1);
				index_data[num_indices++] = I(x + j, patch_lod_size - 1);
			}
			for (int j = patch_bottom_step / 2; j < patch_bottom_step; j++)
			{
				if (x == patch_lod_size - patch_bottom_step && j == patch_bottom_step - 1) { continue; }
				index_data[num_indices++] = I(x + patch_bottom_step, patch_lod_size);
				index_data[num_indices++] = I(x + j + 1, patch_lod_size- 1);
				index_data[num_indices++] = I(x + j, patch_lod_size - 1);
			}
		}
	}

	if (lod_left == lod)
	{
		index_data[num_indices++] = I(0, 0);
		index_data[num_indices++] = I(1, 1);
		index_data[num_indices++] = I(0, 1);

		index_data[num_indices++] = I(0, patch_lod_size - 1);
		index_data[num_indices++] = I(1, patch_lod_size - 1);
		index_data[num_indices++] = I(0, patch_lod_size);
			
		for (int z = 1; z < patch_lod_size - 1; z++)
		{
			index_data[num_indices++] = I(0, z);
			index_data[num_indices++] = I(0, z + 1);
			index_data[num_indices++] = I(1, z);

			index_data[num_indices++] = I(0, z + 1);
			index_data[num_indices++] = I(1, z + 1);
			index_data[num_indices++] = I(1, z);
		}
	}
	else
	{
		for (int z = 0; z < patch_lod_size; z += patch_left_step)
		{
			index_data[num_indices++] = I(0, z);
			index_data[num_indices++] = I(1, z + patch_left_step / 2);
			index_data[num_indices++] = I(0, z + patch_left_step);

			for (int j = 0; j < patch_left_step / 2; j++)
			{
				if (z == 0 && j == 0) { continue; }
				index_data[num_indices++] = I(0, z);
				index_data[num_indices++] = I(1, z + j);
				index_data[num_indices++] = I(1, z + j + 1);
			}
			for (int j = patch_left_step / 2; j < patch_left_step; j++)
			{
				if (z == patch_lod_size - patch_left_step && j == patch_left_step - 1) { continue; }
				index_data[num_indices++] = I(0, z + patch_left_step);
				index_data[num_indices++] = I(1, z + j + 1);
				index_data[num_indices++] = I(1, z + j);
			}
		}
	}

	if (lod_right == lod)
	{
		index_data[num_indices++] = I(patch_lod_size, 0);
		index_data[num_indices++] = I(patch_lod_size - 1, 1);
		index_data[num_indices++] = I(patch_lod_size, 1);

		index_data[num_indices++] = I(patch_lod_size, patch_lod_size - 1);
		index_data[num_indices++] = I(patch_lod_size - 1, patch_lod_size - 1);
		index_data[num_indices++] = I(patch_lod_size, patch_lod_size);
			
		for (int z = 1; z < patch_lod_size - 1; z++)
		{
			index_data[num_indices++] = I(patch_lod_size - 1, z);
			index_data[num_indices++] = I(patch_lod_size - 1, z + 1);
			index_data[num_indices++] = I(patch_lod_size, z);

			index_data[num_indices++] = I(patch_lod_size - 1, z + 1);
			index_data[num_indices++] = I(patch_lod_size, z + 1);
			index_data[num_indices++] = I(patch_lod_size, z);
		}
	}
	else
	{
		for (int z = 0; z < patch_lod_size; z += patch_right_step)
		{
			index_data[num_indices++] = I(patch_lod_size, z);
			index_data[num_indices++] = I(patch_lod_size - 1, z + patch_right_step / 2);
			index_data[num_indices++] = I(patch_lod_size, z + patch_right_step);
			
			for (int j = 0; j < patch_right_step / 2; j++)
			{
				if (z == 0 && j == 0) { continue; }
				index_data[num_indices++] = I(patch_lod_size, z);
				index_data[num_indices++] = I(patch_lod_size - 1, z + j);
				index_data[num_indices++] = I(patch_lod_size - 1, z + j + 1);
			}

			for (int j = patch_right_step / 2; j < patch_right_step; j++)
			{
				if (z == patch_lod_size - patch_right_step && j == patch_right_step - 1) { continue; }
				index_data[num_indices++] = I(patch_lod_size, z + patch_right_step);
				index_data[num_indices++] = I(patch_lod_size - 1, z + j + 1);
				index_data[num_indices++] = I(patch_lod_size - 1, z + j);
			}
		}
	}
#undef I

	for (int i = 0; i < num_indices; i++) { assert(index_data[i] < vertex_pitch * vertex_pitch); }

	n->has_patch = true;
	n->num_patch_indices = num_indices;
	n->patch_vb = gl_buffer_create(GL_ARRAY_BUFFER, vertex_data, sizeof(vertex_data[0]) * vertex_pitch * vertex_pitch, GL_STATIC_DRAW);
	n->patch_ib = gl_buffer_create(GL_ELEMENT_ARRAY_BUFFER, index_data, sizeof(index_data[0]) * num_indices, GL_STATIC_DRAW);

	int tex_res = patch_lod_size * patch_tex_scale;

	glGenTextures(1, &n->patch_lightmap);
	glBindTexture(GL_TEXTURE_2D, n->patch_lightmap);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, tex_res, tex_res, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glBindTexture(GL_TEXTURE_2D, 0);

	//const char* gl_ext = (const char*)glGetString(GL_EXTENSIONS);
	//typedef GLenum (GLAPIENTRY* glCheckFramebufferStatusProc)(GLenum target);
	//glCheckFramebufferStatusProc _glCheckFramebufferStatus = (glCheckFramebufferStatusProc)wglGetProcAddress("glCheckFramebufferStatus");

	GLuint fbo = 0;
	glGenFramebuffersEXT(1, &fbo);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, fbo);
	glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, n->patch_lightmap, 0);
	if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
	{
		printf("Error: FBO for %dx%d lightmap is \"incomplete\", most likely unsupported by the hardware\n", tex_res, tex_res);
		exit(EXIT_FAILURE);
	}

	int viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	glViewport(0, 0, tex_res, tex_res);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);

	int polygon_mode[2];
	glGetIntegerv(GL_POLYGON_MODE, polygon_mode);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glClear(GL_COLOR_BUFFER_BIT);

	GLuint shader = shader_terrain_generate_lightmap;
	glUseProgram(shader);

	GLuint uniform_tex_random = glGetUniformLocation(shader, "tex_random");
	GLuint uniform_gen_horz_scale = glGetUniformLocation(shader, "gen_horz_scale");
	GLuint uniform_gen_horz_minor_scale = glGetUniformLocation(shader, "gen_horz_minor_scale");
	GLuint uniform_map_horz_scale = glGetUniformLocation(shader, "map_horz_scale");
	GLuint uniform_map_vert_scale = glGetUniformLocation(shader, "map_vert_scale");
	GLuint uniform_node_base_xz = glGetUniformLocation(shader, "node_base_xz");
	GLuint uniform_sun_dir = glGetUniformLocation(shader, "sun_dir");
	
	glUniform1f(uniform_gen_horz_scale, map_horz_gen_scale);
	glUniform1f(uniform_gen_horz_minor_scale, map_horz_gen_minor_scale);
	glUniform1f(uniform_map_horz_scale, map_horz_scale);
	glUniform1f(uniform_map_vert_scale, map_vert_scale);
	glUniform2f(uniform_node_base_xz, n->pos_x, n->pos_z);
	glUniform3f(uniform_sun_dir, sun_dir.x, sun_dir.y, sun_dir.z);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, tex_random);
	glUniform1i(uniform_tex_random, 0);

	glBindBuffer(GL_ARRAY_BUFFER, n->patch_vb);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, n->patch_ib);

	GLuint attrib_position = glGetAttribLocation(shader, "position");
	GLuint attrib_texcoord = glGetAttribLocation(shader, "texcoord");
	glEnableVertexAttribArray(attrib_position);
	glEnableVertexAttribArray(attrib_texcoord);

	glVertexAttribPointer(attrib_position, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
	glVertexAttribPointer(attrib_texcoord, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(sizeof(float) * 3));

	glDrawElements(GL_TRIANGLES, n->num_patch_indices, GL_UNSIGNED_INT, 0);

	glDisableVertexAttribArray(attrib_position);
	glDisableVertexAttribArray(attrib_texcoord);

	if (polygon_mode[0] == GL_LINE) { glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); }
	glEnable(GL_DEPTH_TEST);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, cur_fbo);
	glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
	glDeleteFramebuffers(1, &fbo);

	delete[] vertex_data;
	delete[] index_data;
	printf("Generate patch lod -> TRBL: %d -> %d %d %d %d\t%d tri, %d tex res\n", lod, lod_top, lod_right, lod_bottom, lod_left, num_indices / 3, tex_res);
}

int render_sky(const Vector3& camera_pos, const Matrix& wvp)
{
	Vertex sky_vertices[] = 
	{
		{ -100, -100, -100, 0, 0 },
		{ -100, 100, -100, 0, 0 },
		{ 100, -100, -100, 0, 0 },
		{ 100, 100, -100, 0, 0 },
		{ -100, -100, 100, 0, 0 },
		{ -100, 100, 100, 0, 0 },
		{ 100, -100, 100, 0, 0 },
		{ 100, 100, 100, 0, 0 },
	};
	
	GLushort sky_indices[] =
	{
		4, 5, 1, 0, 4, 1,
		3, 7, 6, 3, 6, 2,

		0, 1, 2, 1, 3, 2,
		6, 5, 4, 6, 7, 5,

		0, 2, 4, 2, 6, 4,
		5, 3, 1, 5, 7, 3,
	};

	GLuint vertex_buffer = gl_buffer_create(GL_ARRAY_BUFFER, sky_vertices, sizeof(sky_vertices), GL_STATIC_DRAW);
	GLuint index_buffer = gl_buffer_create(GL_ELEMENT_ARRAY_BUFFER, sky_indices, sizeof(sky_indices), GL_STATIC_DRAW);

	GLuint shader = shader_sky;
	glUseProgram(shader);
	glDisable(GL_DEPTH_TEST);

	GLuint attrib_position = glGetAttribLocation(shader, "position");
	GLuint attrib_texcoord = glGetAttribLocation(shader, "texcoord");

	GLuint uniform_camera_pos = glGetUniformLocation(shader, "camera_pos");
	GLuint uniform_elapsed_time = glGetUniformLocation(shader, "elapsed_time");
	GLuint uniform_wvp = glGetUniformLocation(shader, "wvp");

	GLuint uniform_tex_random = glGetUniformLocation(shader, "tex_random");

	GLuint uniform_fog_distance = glGetUniformLocation(shader, "fog_distance");
	GLuint uniform_dof_focal_distance = glGetUniformLocation(shader, "dof_focal_distance");
	GLuint uniform_atmosphere_color_apex = glGetUniformLocation(shader, "atmosphere_color_apex");
	GLuint uniform_atmosphere_color_dense = glGetUniformLocation(shader, "atmosphere_color_dense");
	
	GLuint uniform_sun_dir = glGetUniformLocation(shader, "sun_dir");

	Matrix m_transform;
	m_transform.Translate(camera_pos);
	m_transform = m_transform * wvp;
	
	glUniformMatrix4fv(uniform_wvp, 1, true, (const GLfloat*)m_transform.cell);
	glUniform1f(uniform_elapsed_time, glutGet(GLUT_ELAPSED_TIME) * 0.001f);
	glUniform3f(uniform_camera_pos, camera_pos.x / map_horz_scale, camera_pos.y / map_horz_scale, camera_pos.z / map_horz_scale);
	glUniform2f(uniform_fog_distance, fog_start, fog_end);
	glUniform1f(uniform_dof_focal_distance, dof_focal_distance);
	glUniform3f(uniform_atmosphere_color_apex, atmosphere_color_apex.x, atmosphere_color_apex.y, atmosphere_color_apex.z);
	glUniform3f(uniform_atmosphere_color_dense, atmosphere_color_dense.x, atmosphere_color_dense.y, atmosphere_color_dense.z);
	glUniform3f(uniform_sun_dir, sun_dir.x, sun_dir.y, sun_dir.z);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, tex_random);
	glUniform1i(uniform_tex_random, 0);
	
	glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, index_buffer);

	glEnableVertexAttribArray(attrib_position);
	glEnableVertexAttribArray(attrib_texcoord);

	glVertexAttribPointer(attrib_position, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
	glVertexAttribPointer(attrib_texcoord, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(sizeof(float) * 3));

	glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_SHORT, 0);

	glDisableVertexAttribArray(attrib_position);
	glDisableVertexAttribArray(attrib_texcoord);

	glDeleteBuffers(1, &vertex_buffer);
	glDeleteBuffers(1, &index_buffer);
	glEnable(GL_DEPTH_TEST);

	return 12;
}

int render(const Vector3& camera_pos, const Matrix& wvp, const Frustum& view_frustum, Node* n)
{
		if (n->distance_cull) { return 0; }
		if (n->has_children)
		{
			int num_tri = 0;
			for (int i = 0; i < 4; i++) { num_tri += render(camera_pos, wvp, view_frustum, n->children[i]); }
			return num_tri;
		}

	Vector3 bounds_min(n->pos_x * map_horz_scale, -map_vert_scale, n->pos_z * map_horz_scale);
	Vector3 bounds_max((n->pos_x + n->size_x) * map_horz_scale, map_vert_scale, (n->pos_z + n->size_z) * map_horz_scale);
	if (!view_frustum.TestAABB(bounds_min, bounds_max)) { return 0; }

	bool needs_update = false;
	for (int i = 0; i < 5; i++) { needs_update |= n->lod[i] != n->patch_lod[i]; }
	
	int tri_drawn = 0;
	if (!n->has_patch || needs_update) { generate_mesh(n); }

	GLuint shader = shader_terrain;
	glUseProgram(shader);
	GLuint attrib_position = glGetAttribLocation(shader, "position");
	GLuint attrib_texcoord = glGetAttribLocation(shader, "texcoord");

	GLuint uniform_camera_pos = glGetUniformLocation(shader, "camera_pos");
	GLuint uniform_wvp = glGetUniformLocation(shader, "wvp");
	
	GLuint uniform_use_composite_texturing = glGetUniformLocation(shader, "use_composite_texturing");

	GLuint uniform_tex_random = glGetUniformLocation(shader, "tex_random");
	GLuint uniform_tex_lightmap = glGetUniformLocation(shader, "tex_lightmap");
	GLuint uniform_tex_scheme = glGetUniformLocation(shader, "tex_scheme");
	GLuint uniform_tex_grass = glGetUniformLocation(shader, "tex_grass");
	GLuint uniform_tex_detail = glGetUniformLocation(shader, "tex_detail");

	GLuint uniform_tex_moss1 = glGetUniformLocation(shader, "tex_moss1");
	GLuint uniform_tex_moss2 = glGetUniformLocation(shader, "tex_moss2");
	GLuint uniform_tex_rock = glGetUniformLocation(shader, "tex_rock");
	GLuint uniform_tex_cliff = glGetUniformLocation(shader, "tex_cliff");

	GLuint uniform_fog_distance = glGetUniformLocation(shader, "fog_distance");
	GLuint uniform_dof_focal_distance = glGetUniformLocation(shader, "dof_focal_distance");
	GLuint uniform_atmosphere_color_apex = glGetUniformLocation(shader, "atmosphere_color_apex");
	GLuint uniform_atmosphere_color_dense = glGetUniformLocation(shader, "atmosphere_color_dense");
	GLuint uniform_gen_horz_scale = glGetUniformLocation(shader, "gen_horz_scale");
	GLuint uniform_gen_horz_minor_scale = glGetUniformLocation(shader, "gen_horz_minor_scale");
	GLuint uniform_map_horz_scale = glGetUniformLocation(shader, "map_horz_scale");
	GLuint uniform_map_vert_scale = glGetUniformLocation(shader, "map_vert_scale");
	GLuint uniform_node_base_xz = glGetUniformLocation(shader, "node_base_xz");
	
	glUniformMatrix4fv(uniform_wvp, 1, true, (const GLfloat*)wvp.cell);
	glUniform3f(uniform_camera_pos, camera_pos.x / map_horz_scale, camera_pos.y / map_horz_scale, camera_pos.z / map_horz_scale);
	glUniform2f(uniform_fog_distance, fog_start, fog_end);
	glUniform1f(uniform_dof_focal_distance, dof_focal_distance);
	glUniform3f(uniform_atmosphere_color_apex, atmosphere_color_apex.x, atmosphere_color_apex.y, atmosphere_color_apex.z);
	glUniform3f(uniform_atmosphere_color_dense, atmosphere_color_dense.x, atmosphere_color_dense.y, atmosphere_color_dense.z);
	glUniform1f(uniform_gen_horz_scale, map_horz_gen_scale);
	glUniform1f(uniform_gen_horz_minor_scale, map_horz_gen_minor_scale);
	glUniform1f(uniform_map_horz_scale, map_horz_scale);
	glUniform1f(uniform_map_vert_scale, map_vert_scale);
	glUniform2f(uniform_node_base_xz, n->pos_x, n->pos_z);
	glUniform1i(uniform_use_composite_texturing, use_composite_texturing && !render_wireframe ? 1 : 0);

	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, n->patch_lightmap);
	glUniform1i(uniform_tex_lightmap, 1);

	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, tex_scheme);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	glUniform1i(uniform_tex_scheme, 2);

	for (int i = 0; i < 5; i++)
	{
		glActiveTexture(GL_TEXTURE3 + i);
		glBindTexture(GL_TEXTURE_2D, tex_terrain[i]);
	}

	glUniform1i(uniform_tex_detail, 3);
	glUniform1i(uniform_tex_moss1, 4);
	glUniform1i(uniform_tex_moss2, 5);
	glUniform1i(uniform_tex_rock, 6);
	glUniform1i(uniform_tex_cliff, 7);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, tex_random);
	glUniform1i(uniform_tex_random, 0);
	
	glBindBuffer(GL_ARRAY_BUFFER, n->patch_vb);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, n->patch_ib);

	glEnableVertexAttribArray(attrib_position);
	glEnableVertexAttribArray(attrib_texcoord);

	glVertexAttribPointer(attrib_position, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)0);
	glVertexAttribPointer(attrib_texcoord, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*)(sizeof(float) * 3));

	glDrawElements(GL_TRIANGLES, n->num_patch_indices, GL_UNSIGNED_INT, 0);

	glDisableVertexAttribArray(attrib_position);
	glDisableVertexAttribArray(attrib_texcoord);
	
	tri_drawn += n->num_patch_indices / 3;
	return tri_drawn;
}

void render_targets_free()
{
	glDeleteFramebuffers(3, dof_fbo);
	glDeleteTextures(3, dof_rt);
	glDeleteTextures(3, dof_rt_z);
	glDeleteTextures(3, dof_rt_depth);
}

void render_targets_create()
{
	for (int i = 0; i < 3; i++)
	{
		int tex_w = i == 0 ? screen_w : screen_w / dof_downsample_factor;
		int tex_h = i == 0 ? screen_h : screen_h / dof_downsample_factor;

		glGenTextures(1, &dof_rt[i]);
		glBindTexture(GL_TEXTURE_2D, dof_rt[i]);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, tex_w, tex_h, 0, GL_BGRA, GL_UNSIGNED_BYTE, 0);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glBindTexture(GL_TEXTURE_2D, 0);

		glGenTextures(1, &dof_rt_z[i]);
		glBindTexture(GL_TEXTURE_2D, dof_rt_z[i]);
		//glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, tex_w, tex_h, 0, GL_RED, GL_FLOAT, 0);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA8, tex_w, tex_h, 0, GL_BGRA, GL_UNSIGNED_BYTE, 0);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
		glBindTexture(GL_TEXTURE_2D, 0);

		glGenTextures(1, &dof_rt_depth[i]);
		glBindTexture(GL_TEXTURE_2D, dof_rt_depth[i]);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_STENCIL_EXT, tex_w, tex_h, 0, GL_DEPTH_STENCIL_EXT, GL_UNSIGNED_INT_24_8_EXT, 0);
		glBindTexture(GL_TEXTURE_2D, 0);

		glGenFramebuffersEXT(1, &dof_fbo[i]);
		glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, dof_fbo[i]);
		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT0_EXT, GL_TEXTURE_2D, dof_rt[i], 0);
		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_COLOR_ATTACHMENT1_EXT, GL_TEXTURE_2D, dof_rt_z[i], 0);
		glFramebufferTexture2DEXT(GL_FRAMEBUFFER_EXT, GL_DEPTH_ATTACHMENT_EXT, GL_TEXTURE_2D, dof_rt_depth[i], 0);
		if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
		{
			printf("Error: Failed to create FBO for depth of field\n");
			exit(EXIT_FAILURE);
		}
		glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
	}
}

void textures_create()
{
	//srand(0);
	float random_data[128 * 128 * 3];
	//for (int i = 0; i < 128 * 128 * 3; i++) { random_data[i] = (rand() % 10000) / 10000.0f; }

	memset(random_data, 0, sizeof(float) * 128 * 128 * 3);
	FILE* noise_file = fopen("noise_seed", "rb");
	if (!noise_file)
	{
		printf("Error: Noise seed file missing\n");
		exit(EXIT_FAILURE);
	}
	int noise_w = 0;
	int noise_h = 0;
	fread(&noise_w, sizeof(int), 1, noise_file);
	fread(&noise_h, sizeof(int), 1, noise_file);
	if (noise_w != 128 || noise_h != 128)
	{
		printf("Error: Seed file must be 128x128x3\n");
		exit(EXIT_FAILURE);
	}

	for (int i = 0; i < 128 * 128; i++)
	{
		float v[4];
		fread(v, sizeof(float), 4, noise_file);
		random_data[i * 3] = v[0];
		random_data[i * 3 + 1] = v[1];
		random_data[i * 3 + 2] = v[2];
	}
	fclose(noise_file);

	glEnable(GL_TEXTURE_2D);
	glGenTextures(1, &tex_random);
	glBindTexture(GL_TEXTURE_2D, tex_random);
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, 128, 128, 0, GL_RGB, GL_FLOAT, random_data);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
	glDisable(GL_TEXTURE_2D);

	tex_scheme = gl_image_load("TerrainScheme.ftex");
	const char* terrain_textures[] =
	{
		"tex_detail.ftex",
		"tex_moss.ftex",
		"tex_moss2.ftex",
		"tex_rock.ftex",
		"tex_cliff.ftex",
	};

	for (int i = 0; i < 5; i++) { tex_terrain[i] = gl_image_load(terrain_textures[i]); }
}

void dof_combine_pass(GLuint src_tex, GLuint src_z, GLuint blurred_tex)
{
	float vertex_data[] =
	{
		-1, -1,
		1, -1,
		1, 1,
		-1, 1,
	};

	GLushort index_data[] = { 0, 1, 2, 0, 2, 3 };

	glUseProgram(shader_dof_combine);
	GLuint attrib_pos_loc = glGetAttribLocation(shader_dof_combine, "position");
	GLuint uniform_basetex_loc = glGetUniformLocation(shader_dof_combine, "basetex");
	GLuint uniform_blurredtex_loc = glGetUniformLocation(shader_dof_combine, "blurredtex");
	GLuint uniform_z_tex_loc = glGetUniformLocation(shader_dof_combine, "z_tex");

	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, src_z);
	glUniform1i(uniform_z_tex_loc, 1);

	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, blurred_tex);
	glUniform1i(uniform_blurredtex_loc, 2);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, src_tex);
	glUniform1i(uniform_basetex_loc, 0);

	GLuint vb = gl_buffer_create(GL_ARRAY_BUFFER, vertex_data, sizeof(float) * 8, GL_STATIC_DRAW);
	GLuint ib = gl_buffer_create(GL_ELEMENT_ARRAY_BUFFER, index_data, sizeof(GLushort) * 6, GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, vb);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ib);

	glEnableVertexAttribArray(attrib_pos_loc);
	glVertexAttribPointer(attrib_pos_loc, 2, GL_FLOAT, GL_FALSE, sizeof(float) * 2, (void*)0);

	glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_SHORT, 0);

	glDisableVertexAttribArray(attrib_pos_loc);

	glDeleteBuffers(1, &vb);
	glDeleteBuffers(1, &ib);
}

void dof_blur_pass(GLuint src_tex, GLuint dest_fbo, GLuint z_tex, int src_w, int src_h, int dest_w, int dest_h, bool is_horz)
{
	GLuint shader = is_horz ? shader_dof_blur_horz : shader_dof_blur_vert;

	glBindFramebuffer(GL_FRAMEBUFFER_EXT, dest_fbo);
	glViewport(0, 0, dest_w, dest_h);

	float vertex_data[] =
	{
		-1, -1,
		1, -1,
		1, 1,
		-1, 1,
	};

	GLushort index_data[] = { 0, 1, 2, 0, 2, 3 };

	glUseProgram(shader);
	GLuint attrib_pos_loc = glGetAttribLocation(shader, "position");
	GLuint uniform_tex_w_loc = glGetUniformLocation(shader, "tex_w");
	GLuint uniform_tex_h_loc = glGetUniformLocation(shader, "tex_h");
	GLuint uniform_tex_loc = glGetUniformLocation(shader, "basetex");
	GLuint uniform_z_tex_loc = glGetUniformLocation(shader, "z_tex");

	glUniform1f(uniform_tex_w_loc, src_w);
	glUniform1f(uniform_tex_h_loc, src_h);

	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, z_tex);
	glUniform1i(uniform_z_tex_loc, 1);

	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, src_tex);
	glUniform1i(uniform_tex_loc, 0);

	GLuint vb = gl_buffer_create(GL_ARRAY_BUFFER, vertex_data, sizeof(float) * 8, GL_STATIC_DRAW);
	GLuint ib = gl_buffer_create(GL_ELEMENT_ARRAY_BUFFER, index_data, sizeof(GLushort) * 6, GL_STATIC_DRAW);

	glBindBuffer(GL_ARRAY_BUFFER, vb);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ib);

	glEnableVertexAttribArray(attrib_pos_loc);
	glVertexAttribPointer(attrib_pos_loc, 2, GL_FLOAT, GL_FALSE, sizeof(float) * 2, (void*)0);

	glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_SHORT, 0);

	glDisableVertexAttribArray(attrib_pos_loc);

	glDeleteBuffers(1, &vb);
	glDeleteBuffers(1, &ib);
}

void display()
{
	int tris_drawn = 0;
	frame++;
	glClearColor(0, 0, 0, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	glViewport(0, 0, screen_w, screen_h);

	int num_nodes = update_lod(camera_pos * (1.0f / map_horz_scale), terrain_root);
	update_edge_lod(terrain_root);

	Matrix m_camera1, m_camera2;
	m_camera1.RotateZ(-camera_hangle);
	m_camera2.RotateY(-camera_vangle);
	Vector3 camera_up = (m_camera1 * m_camera2).TransformCoord(Vector3(0, 1, 0));
	camera_up.Normalize();

	Vector3 camera_dir = (m_camera1 * m_camera2).TransformCoord(Vector3(1, 0, 0));
	camera_dir.Normalize();

	Matrix m_view, m_proj, m_viewproj, m_world_scale;
	m_proj.Perspective(PI / 4.0f, (float)screen_w / (float)screen_h, frustum_near_plane, frustum_far_plane);
	m_view.LookAt(camera_pos, camera_pos + camera_dir, camera_up);

	m_viewproj = m_view * m_proj;
	m_world_scale.Scale(map_horz_scale, 1, map_horz_scale);

	Frustum view_frustum(m_viewproj);

	glEnable(GL_DEPTH_TEST);

	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
	cur_fbo = 0;

	int viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	//glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	//tris_drawn += render(camera_pos, m_world_scale * m_viewproj, view_frustum, terrain_root);
	//glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	if (use_dof && !render_wireframe)
	{
		glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, dof_fbo[0]);
		cur_fbo = dof_fbo[0];
		glDrawBuffers(2, draw_buffers);
		glViewport(0, 0, screen_w, screen_h);
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	}

	if (render_wireframe) { glPolygonMode(GL_FRONT_AND_BACK, GL_LINE); }
	else { tris_drawn += render_sky(camera_pos, m_viewproj); }
	tris_drawn += render(camera_pos, m_world_scale * m_viewproj, view_frustum, terrain_root);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	glDrawBuffers(1, draw_buffers);
	glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
	cur_fbo = 0;

	glDisable(GL_DEPTH_TEST);

	if (use_dof && !render_wireframe)
	{
		for (int i = 0; i < dof_blur_passes; i++)
		{
			dof_blur_pass(dof_rt[i == 0 ? 0 : 1], dof_fbo[2], dof_rt_z[0], screen_w / dof_downsample_factor, screen_h / dof_downsample_factor, screen_w / dof_downsample_factor, screen_h / dof_downsample_factor, true);
			dof_blur_pass(dof_rt[2], dof_fbo[1], dof_rt_z[0], screen_w / dof_downsample_factor, screen_h / dof_downsample_factor, screen_w / dof_downsample_factor, screen_h / dof_downsample_factor, false);
		}

		glBindFramebufferEXT(GL_FRAMEBUFFER_EXT, 0);
		cur_fbo = 0;
		glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
		dof_combine_pass(dof_rt[0], dof_rt_z[0], dof_rt[1]);
	}

	char printbuf[2048];
	snprintf(printbuf, sizeof(printbuf), 
		"$cFFFFPS: $c8F8%d$cFFF\n"
		"%d triangles from %d nodes\n"
		"Render mode: $c8F8%s$cFFF $cFF0[1/2]$cFFF\n"
		"Fast camera:  $c8F8%s$cFFF $cFF0[q]$cFFF\n"
		"Composite texturing: $c8F8%s$cFFF $cFF0[o]$cFFF\n"
		"Depth of field: $c8F8%s$cFFF $cFF0[p]$cFFF\n",
		fps,
		tris_drawn, num_nodes,
		render_wireframe ? "Wireframe" : "Solid",
		use_fast_camera ? "Yes" : "No",
		use_composite_texturing ? "Yes" : "No",
		use_dof ? "Yes" : "No");
	gl_font_draw(&font, printbuf, 4, 4, 0.3f, screen_w, screen_h);
	glutSwapBuffers();
}

void window_resized(int w, int h)
{
	screen_w = w;
	screen_h = h;

	// Render targets were set up for the previous screen resolution
	render_targets_free();
	render_targets_create();
}

void mouse_input(int button, int state, int x, int y)
{
	if (button == GLUT_LEFT_BUTTON)
	{
		int cur_tick = glutGet(GLUT_ELAPSED_TIME);
		glut_mouse_x = x;
		glut_mouse_y = y;

		glut_left_mouse_state = state;
	}
}

void mouse_motion_input(int x, int y)
{
	if (glut_left_mouse_state == GLUT_DOWN)
	{
		camera_vangle += (x - glut_mouse_x) / 2.0f;
		camera_hangle -= (y - glut_mouse_y) / 2.0f;
		if (camera_hangle < -89) { camera_hangle = -89; }
		if (camera_hangle > 89) { camera_hangle = 89; }
	}
	glut_mouse_x = x;
	glut_mouse_y = y;
}

void key_input(unsigned char key, int x, int y)
{
	const unsigned char WIN_Escape = 27;
	if (key >= 'a' && key <= 'z') { key -= 'a' - 'A'; }
	if (key == WIN_Escape) { exit(EXIT_SUCCESS); }

	if (key == 'W') { camera_motion.x = 1; }
	if (key == 'S') { camera_motion.x = -1; }
	if (key == 'A') { camera_motion.y = 1; }
	if (key == 'D') { camera_motion.y = -1; }
	if (key == 'O') { use_dof = !use_dof; }
	if (key == 'P') { use_composite_texturing = !use_composite_texturing; }
	if (key == 'Q') { use_fast_camera = !use_fast_camera; }
	if (key == '1') { render_wireframe = false; }
	if (key == '2') { render_wireframe = true; }
}

void key_up_input(unsigned char key, int x, int y)
{
	if (key >= 'a' && key <= 'z') { key -= 'a' - 'A'; }
	if (key == 'W' || key == 'S') { camera_motion.x = 0; }
	if (key == 'A' || key == 'D') { camera_motion.y = 0; }
}

void idle()
{
	int time = glutGet(GLUT_ELAPSED_TIME);
	int dt = time - prev_time;

	Matrix m_camera1, m_camera2;
	m_camera1.RotateZ(-camera_hangle);
	m_camera2.RotateY(-camera_vangle);
	Vector3 camera_up = (m_camera1 * m_camera2).TransformCoord(Vector3(0, 1, 0));
	camera_up.Normalize();

	Vector3 camera_dir = (m_camera1 * m_camera2).TransformCoord(Vector3(1, 0, 0));
	camera_dir.Normalize();

	Vector3 camera_strafe = camera_dir.Cross(camera_up);
	camera_strafe.Normalize();

	camera_pos += dt * (use_fast_camera ? camera_speed_fast : camera_speed)
		* (camera_motion.x * camera_dir + camera_motion.y * camera_strafe);

	prev_time = time;
	if (time_fps == 0) { time_fps = time; }
	if (time - time_fps > 1000)
	{
		fps = (int)(frame * 1000.0f / (time - time_fps));
		frame = 0;
		time_fps = time;
	}

	glutPostRedisplay();
}

int main(int argc, char** argv)
{
	srand(0);

	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitWindowSize(screen_w, screen_h);
	glutCreateWindow("Nauful's final project");

	glewInit();

	font = gl_font_load("font.df", "font_shader.glsl");

	shader_terrain_generate_lightmap = gl_program_create(
		gl_shader_create("terrain.glsl", GL_VERTEX_SHADER, "#define VERTEX_SHADER\r\n#define PASS_LIGHTMAP\r\n"),
		gl_shader_create("terrain.glsl", GL_FRAGMENT_SHADER, "#define FRAGMENT_SHADER\r\n#define PASS_LIGHTMAP\r\n"));

	shader_terrain = gl_program_create(
		gl_shader_create("terrain.glsl", GL_VERTEX_SHADER, "#define VERTEX_SHADER\r\n#define PASS_RENDER\r\n"),
		gl_shader_create("terrain.glsl", GL_FRAGMENT_SHADER, "#define FRAGMENT_SHADER\r\n#define PASS_RENDER\r\n"));

	shader_sky = gl_program_create(
		gl_shader_create("terrain.glsl", GL_VERTEX_SHADER, "#define VERTEX_SHADER\r\n#define PASS_SKY\r\n"),
		gl_shader_create("terrain.glsl", GL_FRAGMENT_SHADER, "#define FRAGMENT_SHADER\r\n#define PASS_SKY\r\n"));

	shader_dof_combine = gl_program_create(
		gl_shader_create("dof.glsl", GL_VERTEX_SHADER, "#define VERTEX_SHADER\r\n#define PASS_COMBINE\r\n"),
		gl_shader_create("dof.glsl", GL_FRAGMENT_SHADER, "#define FRAGMENT_SHADER\r\n#define PASS_COMBINE\r\n"));

	shader_dof_blur_horz = gl_program_create(
		gl_shader_create("dof.glsl", GL_VERTEX_SHADER, "#define VERTEX_SHADER\r\n#define PASS_HORIZONTAL\r\n"),
		gl_shader_create("dof.glsl", GL_FRAGMENT_SHADER, "#define FRAGMENT_SHADER\r\n#define PASS_HORIZONTAL\r\n"));

	shader_dof_blur_vert = gl_program_create(
		gl_shader_create("dof.glsl", GL_VERTEX_SHADER, "#define VERTEX_SHADER\r\n#define PASS_VERTICAL\r\n"),
		gl_shader_create("dof.glsl", GL_FRAGMENT_SHADER, "#define FRAGMENT_SHADER\r\n#define PASS_VERTICAL\r\n"));

	sun_dir.Normalize();

	render_targets_create();
	textures_create();

	terrain_root = node_create(0, -map_size / 2, -map_size / 2, map_size, map_size);

	glutIdleFunc(idle);
	glutDisplayFunc(display);
	glutKeyboardFunc(key_input);
	glutKeyboardUpFunc(key_up_input);
	glutMouseFunc(mouse_input);
	glutMotionFunc(mouse_motion_input);
	glutReshapeFunc(window_resized);
	glutMainLoop();

	return 0;
}
