#include "math_util.h"

Vector3::Vector3()
	: x(0), y(0), z(0)
{
}

Vector3::Vector3(float x, float y, float z)
{
	this->x = x;
	this->y = y;
	this->z = z;
}

float Vector3::Length() const
{
	return sqrtf(x * x + y * y + z * z);
}

float Vector3::LengthSquared() const
{
	return x * x + y * y + z * z;
}

void Vector3::Normalize()
{
	float len = Length();
	if (len > 0)
	{
		x /= len;
		y /= len;
		z /= len;
	}
}

Vector3 Vector3::Cross(const Vector3& b) const
{
	return Vector3(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x);
}

float Vector3::Dot(const Vector3& b) const
{
	return x * b.x + y * b.y + z * b.z;
}

Vector3 Vector3::Minimize(const Vector3& a, const Vector3& b)
{
	return Vector3(min(a.x, b.x), min(a.y, b.y), min(a.z, b.z));
}

Vector3 Vector3::Maximize(const Vector3& a, const Vector3& b)
{
	return Vector3(max(a.x, b.x), max(a.y, b.y), max(a.z, b.z));
}

Vector3 Vector3::Lerp(const Vector3& a, const Vector3& b, float t)
{
	return a + (b - a) * t;
}

Vector3 Vector3::SplineCatmullRom(const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3, float t)
{
	return 0.5f * (
			(2.0f * p1) +
			(-p0 + p2) * t +
			(2.0f * p0 - 5.0f * p1 + 4.0f * p2 - p3) * t * t +
			(-p0 + 3.0f * p1 - 3.0f * p2 + p3) * t * t * t);
}

float Vector3::MaxComponent() const
{
	return max(max(x, y), z);
}

float Vector3::MinComponent() const
{
	return min(min(x, y), z);
}

Vector3 Vector3::ProductMultiply(const Vector3& a) const
{
	return Vector3(x * a.x, y * a.y, z * a.z);
}

Matrix::Matrix()
{
	Identity();
}

Matrix::Matrix(float* m)
{
	Set(m);
}

void Matrix::Identity()
{
	for (int i = 0; i < 16; i++) { cell[i] = 0; }
	cell[0] = cell[5] = cell[10] = cell[15] = 1;
}

void Matrix::Set(float* m)
{
	for (int i = 0; i < 16; i++) { cell[i] = m[i]; }
}

void Matrix::Transpose()
{
	float mt[16];
	for (int x = 0; x < 4; x++)
	{
		for (int y = 0; y < 4; y++)
		{
			mt[x * 4 + y] = m[y][x];
		}
	}
	Set(mt);
}

void Matrix::Perspective(float fov, float aspect_ratio, float znear, float zfar)
{
	float s = cosf(0.5f * fov * aspect_ratio) / sinf(0.5f * fov);
	float zdiff = zfar - znear;

	float m[] = 
	{
		s / aspect_ratio, 0, 0, 0,
		0, s, 0, 0,
		0, 0, zfar / zdiff, (-zfar * znear) / zdiff,
		0, 0, 1, 0
	};

	Set(m);
}

void Matrix::LookAt(const Vector3& position, const Vector3& target, const Vector3& up)
{
	Vector3 zaxis = target - position;
	zaxis.Normalize();

	Vector3 xaxis = up.Cross(zaxis);
	Vector3 yaxis = zaxis.Cross(xaxis);

	float m[] =
	{
		xaxis.x, xaxis.y, xaxis.z, -position.Dot(xaxis),
		yaxis.x, yaxis.y, yaxis.z, -position.Dot(yaxis),
		zaxis.x, zaxis.y, zaxis.z, -position.Dot(zaxis),
		0, 0, 0, 1,
	};

	Set(m);
}

void Matrix::Concatenate(const Matrix& m2)
{
	Matrix res;
	for (int c = 0; c < 4; c++)
	{
		for (int r = 0; r < 4; r++)
		{
			res.cell[r * 4 + c] = cell[r * 4] * m2.cell[c] +
								  cell[r * 4 + 1] * m2.cell[c + 4] +
								  cell[r * 4 + 2] * m2.cell[c + 8] +
								  cell[r * 4 + 3] * m2.cell[c + 12];
		}
	}
	for (int c = 0; c < 16; c++) { cell[c] = res.cell[c]; }
}

void Matrix::RotateX(float r)
{
	float sx = (float)sinf(r * PI / 180);
	float cx = (float)cosf(r * PI / 180);
	Identity();
	cell[5] = cx;
	cell[6] = sx;
	cell[9] = -sx;
	cell[10] = cx;
}

void Matrix::RotateY(float r)
{
	float sy = (float)sinf(r * PI / 180);
	float cy = (float)cosf(r * PI / 180);
	Identity();
	cell[0] = cy;
	cell[2] = -sy;
	cell[8] = sy;
	cell[10] = cy;
}

void Matrix::RotateZ(float r)
{
	float sz = (float)sinf(r * PI / 180);
	float cz = (float)cosf(r * PI / 180);
	Identity();
	cell[0] = cz;
	cell[1] = sz;
	cell[4] = -sz;
	cell[5] = cz;
}

void Matrix::RotateAxis(float angle, Vector3 dir)
{
	dir.Normalize();
	Identity();

	float s = (float)sinf(0.5f * angle * PI / 180);
	float c = (float)cosf(0.5f * angle * PI / 180);
	Identity();
	dir.Normalize();

	float x = dir.x * s;
	float y = dir.y * s;
	float z = dir.z * s;

	float x2 = x + x;
	float y2 = y + y;
	float z2 = z + z;

	float xx = x * x2;
	float xy = x * y2;
	float xz = x * z2;

	float yy = y * y2;
	float yz = y * z2;
	float zz = z * z2;

	float wx = c * x2;
	float wy = c * y2;
	float wz = c * z2;

	m[0][0] = 1.0f - (yy + zz);
	m[0][1] = xy - wz;
	m[0][2] = xz + wy;

	m[1][0] = xy + wz;
	m[1][1] = 1.0f - (xx + zz);
	m[1][2] = yz - wx;

	m[2][0] = xz - wy;
	m[2][1] = yz + wx;
	m[2][2] = 1.0f - (xx + yy);

	m[3][3] = 1.0f;
}

void Matrix::Invert()
{
	float det, inv_det;

	float det2_01_01 = m[0][0] * m[1][1] - m[0][1] * m[1][0];
	float det2_01_02 = m[0][0] * m[1][2] - m[0][2] * m[1][0];
	float det2_01_03 = m[0][0] * m[1][3] - m[0][3] * m[1][0];
	float det2_01_12 = m[0][1] * m[1][2] - m[0][2] * m[1][1];
	float det2_01_13 = m[0][1] * m[1][3] - m[0][3] * m[1][1];
	float det2_01_23 = m[0][2] * m[1][3] - m[0][3] * m[1][2];

	float det3_201_012 = m[2][0] * det2_01_12 - m[2][1] * det2_01_02 + m[2][2] * det2_01_01;
	float det3_201_013 = m[2][0] * det2_01_13 - m[2][1] * det2_01_03 + m[2][3] * det2_01_01;
	float det3_201_023 = m[2][0] * det2_01_23 - m[2][2] * det2_01_03 + m[2][3] * det2_01_02;
	float det3_201_123 = m[2][1] * det2_01_23 - m[2][2] * det2_01_13 + m[2][3] * det2_01_12;

	det = -det3_201_123 * m[3][0] + det3_201_023 * m[3][1] - det3_201_013 * m[3][2] + det3_201_012 * m[3][3];

	if (det == 0.0f) { return; }

	inv_det = 1.0f / det;

	float det2_03_01 = m[0][0] * m[3][1] - m[0][1] * m[3][0];
	float det2_03_02 = m[0][0] * m[3][2] - m[0][2] * m[3][0];
	float det2_03_03 = m[0][0] * m[3][3] - m[0][3] * m[3][0];
	float det2_03_12 = m[0][1] * m[3][2] - m[0][2] * m[3][1];
	float det2_03_13 = m[0][1] * m[3][3] - m[0][3] * m[3][1];
	float det2_03_23 = m[0][2] * m[3][3] - m[0][3] * m[3][2];

	float det2_13_01 = m[1][0] * m[3][1] - m[1][1] * m[3][0];
	float det2_13_02 = m[1][0] * m[3][2] - m[1][2] * m[3][0];
	float det2_13_03 = m[1][0] * m[3][3] - m[1][3] * m[3][0];
	float det2_13_12 = m[1][1] * m[3][2] - m[1][2] * m[3][1];
	float det2_13_13 = m[1][1] * m[3][3] - m[1][3] * m[3][1];
	float det2_13_23 = m[1][2] * m[3][3] - m[1][3] * m[3][2];

	float det3_203_012 = m[2][0] * det2_03_12 - m[2][1] * det2_03_02 + m[2][2] * det2_03_01;
	float det3_203_013 = m[2][0] * det2_03_13 - m[2][1] * det2_03_03 + m[2][3] * det2_03_01;
	float det3_203_023 = m[2][0] * det2_03_23 - m[2][2] * det2_03_03 + m[2][3] * det2_03_02;
	float det3_203_123 = m[2][1] * det2_03_23 - m[2][2] * det2_03_13 + m[2][3] * det2_03_12;

	float det3_213_012 = m[2][0] * det2_13_12 - m[2][1] * det2_13_02 + m[2][2] * det2_13_01;
	float det3_213_013 = m[2][0] * det2_13_13 - m[2][1] * det2_13_03 + m[2][3] * det2_13_01;
	float det3_213_023 = m[2][0] * det2_13_23 - m[2][2] * det2_13_03 + m[2][3] * det2_13_02;
	float det3_213_123 = m[2][1] * det2_13_23 - m[2][2] * det2_13_13 + m[2][3] * det2_13_12;

	float det3_301_012 = m[3][0] * det2_01_12 - m[3][1] * det2_01_02 + m[3][2] * det2_01_01;
	float det3_301_013 = m[3][0] * det2_01_13 - m[3][1] * det2_01_03 + m[3][3] * det2_01_01;
	float det3_301_023 = m[3][0] * det2_01_23 - m[3][2] * det2_01_03 + m[3][3] * det2_01_02;
	float det3_301_123 = m[3][1] * det2_01_23 - m[3][2] * det2_01_13 + m[3][3] * det2_01_12;

	m[0][0] = -det3_213_123 * inv_det; m[1][0] = det3_213_023 * inv_det; m[2][0] = -det3_213_013 * inv_det; m[3][0] = det3_213_012 * inv_det;
	m[0][1] = det3_203_123 * inv_det; m[1][1] = -det3_203_023 * inv_det; m[2][1] = det3_203_013 * inv_det; m[3][1] = -det3_203_012 * inv_det;
	m[0][2] = det3_301_123 * inv_det; m[1][2] = -det3_301_023 * inv_det; m[2][2] = det3_301_013 * inv_det; m[3][2] = -det3_301_012 * inv_det;
	m[0][3] = -det3_201_123 * inv_det; m[1][3] = det3_201_023 * inv_det; m[2][3] = -det3_201_013 * inv_det; m[3][3] = det3_201_012 * inv_det;
}

Matrix operator *(const Matrix& a, const Matrix& b)
{
	Matrix r = b;
	r.Concatenate(a);
	return r;
}

Vector3 Matrix::FirstRow3() const
{
	return Vector3(m[0][0], m[1][0], m[2][0]);
}

// Transform with w-divide
Vector3 Matrix::TransformCoord(const Vector3& v) const
{
	float w = m[3][0] * v.x + m[3][1] * v.y + m[3][2] * v.z + m[3][3];
	if (w == 0.0f) { return Vector3(0, 0, 0); }

	float x = m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z + m[0][3]; x /= w;
	float y = m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z + m[1][3]; y /= w;
	float z = m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z + m[2][3]; z /= w;
	return Vector3(x, y, z);
}

Vector3 Matrix::TransformVector(const Vector3& v) const
{
	float x = m[0][0] * v.x + m[0][1] * v.y + m[0][2] * v.z;
	float y = m[1][0] * v.x + m[1][1] * v.y + m[1][2] * v.z;
	float z = m[2][0] * v.x + m[2][1] * v.y + m[2][2] * v.z;
	return Vector3(x, y, z);
}

void Matrix::SetRows(const Vector3& v0, const Vector3& v1, const Vector3& v2)
{
	// Transposed
	float m[] = 
	{
		v0.x, v1.x, v2.x, 0,
		v0.y, v1.y, v2.y, 0,
		v0.z, v1.z, v2.z, 0,
		0, 0, 0, 1,
	};
	Set(m);
}

void Matrix::Scale(float x, float y, float z)
{
	cell[0] *= x;
	cell[5] *= y;
	cell[10] *= z;
}

void Matrix::Translate(const Vector3& p)
{
	cell[3] += p.x;
	cell[7] += p.y;
	cell[11] += p.z;
}

Plane::Plane()
	: x(0), y(0), z(0), d(0)
{
}

Plane::Plane(Vector3 v, float d)
	: x(v.x), y(v.y), z(v.z), d(d)
{
}

Plane::Plane(Vector3 p1, Vector3 p2, Vector3 p3)
{
	Vector3 edge1 = p2 - p1;
	Vector3 edge2 = p3 - p1;
	Vector3 normal = edge1.Cross(edge2);
	normal.Normalize();

	x = normal.x; y = normal.y; z = normal.z;
	d = -normal.Dot(p1);
}

Plane::Plane(Vector3 point, Vector3 normal)
{
	normal.Normalize();
	this->x = normal.x;
	this->y = normal.y;
	this->z = normal.z;
	d = -normal.Dot(point);
}

Vector3 Plane::Normal() const
{
	return Vector3(x, y, z);
}

float Plane::Dot(const Vector3& v) const
{
	return v.x * x + v.y * y + v.z * z + d;
}

//origin + t * dir = p
//origin.n + t * dir.n = p.n
//origin.n + t * dir.n = -d
//t * dir.n = (-d - origin.n)
//t = -(d + origin.n)/(dir.n)

Vector3 Plane::IntersectRay(const Vector3& p1, const Vector3& p2) const
{
	Vector3 normal(x, y, z);
	Vector3 direction = p2 - p1;
	float t = normal.Dot(direction);
	if (t != 0.0f) { t = -(normal.Dot(p1) + d) / t; }

	return p1 + direction * t;
}

bool Plane::IntersectLine(const Vector3& p1, const Vector3& p2, Vector3& hitpt, float epsilon) const
{
	Vector3 normal(x, y, z);
	Vector3 direction = p2 - p1;
	float t = normal.Dot(direction);
	if (fabsf(t) < epsilon) 
	{
		return false; 
	}
	t = -(normal.Dot(p1) + d) / t;

	if (t < -epsilon || t > 1 + epsilon) { return false; }
	hitpt = p1 + direction * t;

	return true;
}