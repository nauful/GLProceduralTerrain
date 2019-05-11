#ifndef MATH_UTIL_H
#define MATH_UTIL_H

#include <math.h>

#ifndef PI
#define PI 3.14159265
#endif

#ifndef max
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef min
#define min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

struct Vector3
{
	union
	{
		struct { float x, y, z; };
		struct { float r, g, b; };
		struct { float cell[3]; };
	};

	Vector3();
	Vector3(float x, float y, float z);

	void Normalize();
	float Length() const;
	float LengthSquared() const;

	Vector3 Cross(const Vector3& b) const;
	float Dot(const Vector3& b) const;

	static Vector3 Minimize(const Vector3& a, const Vector3& b);
	static Vector3 Maximize(const Vector3& a, const Vector3& b);
	static Vector3 Lerp(const Vector3& a, const Vector3& b, float t);
	static Vector3 SplineCatmullRom(const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3, float t);

	Vector3 ProductMultiply(const Vector3& a) const;

	float MaxComponent() const;
	float MinComponent() const;

	void operator +=(const Vector3& V) { x += V.x; y += V.y; z += V.z; }
	void operator +=(const Vector3* V) { x += V->x; y += V->y; z += V->z; }
	void operator -=(const Vector3& V) { x -= V.x; y -= V.y; z -= V.z; }
	void operator -=(const Vector3* V) { x -= V->x; y -= V->y; z -= V->z; }
	void operator *=(float f) { x *= f; y *= f; z *= f; }
	void operator *=(const Vector3& V) { x *= V.x; y *= V.y; z *= V.z; }
	void operator *=(const Vector3* V) { x *= V->x; y *= V->y; z *= V->z; }
	Vector3 operator-() const { return Vector3(-x, -y, -z); }
	friend Vector3 operator +(const Vector3& v1, const Vector3& v2) { return Vector3(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z); }
	friend Vector3 operator -(const Vector3& v1, const Vector3& v2) { return Vector3(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z); }
	friend Vector3 operator +(const Vector3& v1, Vector3* v2) { return Vector3(v1.x + v2->x, v1.y + v2->y, v1.z + v2->z); }
	friend Vector3 operator -(const Vector3& v1, Vector3* v2) { return Vector3(v1.x - v2->x, v1.y - v2->y, v1.z - v2->z); }
	friend Vector3 operator *(const Vector3& v, float f) { return Vector3(v.x * f, v.y * f, v.z * f); }
	friend Vector3 operator *(const Vector3& v1, Vector3& v2) { return Vector3(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z); }
	friend Vector3 operator *(float f, const Vector3& v) { return Vector3(v.x * f, v.y * f, v.z * f); }

	operator const float*() const { return cell; }
	operator float*() { return cell; }
};

struct Matrix
{
	union
	{
		// cells stored in column-major
		struct { float cell[16]; };
		// m[col][row] or m[x][y]
		struct { float m[4][4]; };
	};

	Matrix();
	Matrix(float* m);

	void Identity();
	void Set(float* m);
	void SetRows(const Vector3& v0, const Vector3& v1, const Vector3& v2);
	void Invert();
	void Transpose();
	void Concatenate(const Matrix& m2);
	void Scale(float x, float y, float z);
	void Translate(const Vector3& p);

	void Perspective(float fov, float aspect_ratio, float znear, float zfar);
	void LookAt(const Vector3& position, const Vector3& target, const Vector3& up);
	void RotateX(float r); // degrees
	void RotateY(float r); // degrees
	void RotateZ(float r); // degrees
	void RotateAxis(float angle, Vector3 dir); // degrees around dir

	Vector3 FirstRow3() const;
	Vector3 TransformCoord(const Vector3& v) const;
	Vector3 TransformVector(const Vector3& v) const;

	friend Matrix operator *(const Matrix& a, const Matrix& b);
};

struct Plane
{
	float x, y, z, d;

	Plane();
	Plane(Vector3 v, float d);
	Plane(Vector3 point, Vector3 normal);
	Plane(Vector3 v0, Vector3 v1, Vector3 v2);
	Vector3 Normal() const;
	float Dot(const Vector3& v) const;

	Vector3 IntersectRay(const Vector3& p1, const Vector3& p2) const;
	bool IntersectLine(const Vector3& p1, const Vector3& p2, Vector3& hitpt, float epsilon = 0) const;
};

#endif