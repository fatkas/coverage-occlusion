#pragma once

#include <math.h>
#include <string.h>

struct vec2_t
{
	union
	{
		struct
		{
			float x;
			float y;
		};
		float v[2];
	};	
};

struct vec3_t
{
    union
    {
        struct
        {
            float x;
            float y;
            float z;
        };
        float v[3];
    };
};

#define __forceinline __attribute__((always_inline))
#define atomic_add(a,b) __atomic_add_fetch(&a, b, __ATOMIC_RELAXED) - (b);

#define USE_SSE !defined(__arm64__)
#define USE_NEON defined(__arm64__)

#if USE_SSE
#include "rasterizer_math_sse.inl"
#endif

#if USE_NEON
#include "rasterizer_math_neon.inl"
#endif

inline vec4_t Vector3Dotv(const vec4_t& a, const vec4_t& b)
{
    vec4_t ab = VecMul(a, b);
    vec4_t x = VecShuffle(ab, ab, VecShuffleMask(0, 0, 0, 0));
    vec4_t y = VecShuffle(ab, ab, VecShuffleMask(1, 1, 1, 1));
    vec4_t z = VecShuffle(ab, ab, VecShuffleMask(2, 2, 2, 2));
	return VecAdd(x, VecAdd(y, z));
}

inline float Vector3Dot(vec4_t a, vec4_t b)
{
	float f;
    VecStoreS(&f, Vector3Dotv(a, b));
	return f;
}

inline vec4_t Vector3Lengthv(vec4_t v)
{
	return VecSqrt(Vector3Dotv(v, v));
}

inline float Vector3Length(vec4_t v)
{
	float f;
    VecStoreS(&f, Vector3Lengthv(v));
	return f;
}

inline vec4_t Vector3Cross(vec4_t a, vec4_t b)
{
    vec4_t v1, v2;
	v1 = VecShuffle(a, a, VecShuffleMask(1, 2, 0, 3));
	v2 = VecMul(VecShuffle(b, b, VecShuffleMask(2, 0, 1, 3)), v1);
	v1 = VecMul(v1, b);
	v1 = VecShuffle(v1, v1, VecShuffleMask(1, 2, 0, 3));
	return VecSub(v2, v1);
}

inline vec4_t Vector3Normalize(vec4_t a)
{
	return VecDiv(a, Vector3Lengthv(a));
}

__forceinline inline void Vector3TransformCoord4(const vec4_t* m, const vec4_t* src, vec2_t* dst)
{
    vec4_t src0 = src[0];
    vec4_t src1 = src[1];
    vec4_t src2 = src[2];
    vec4_t src3 = src[3];

    vec4_t tmp0 = VecUnpackLo(src0, src1); // x0 x1 y0 y1
    vec4_t tmp1 = VecUnpackHi(src0, src1); // z0 z1 w0 w1
    vec4_t tmp2 = VecUnpackLo(src2, src3); // x2 x3 y2 y3
    vec4_t tmp3 = VecUnpackHi(src2, src3); // z2 z3 w2 w3

    vec4_t xxxx = VecMoveLH(tmp0, tmp2); // x0 x1 x2 x3
    vec4_t yyyy = VecMoveHL(tmp0, tmp2); // y0 y1 y2 y3
    vec4_t zzzz = VecMoveLH(tmp1, tmp3); // z0 z1 z2 z3

    vec4_t tx = VecMad(m[ 8], zzzz, VecMad(m[4], yyyy, VecMad(m[0], xxxx, m[12])));
    vec4_t ty = VecMad(m[ 9], zzzz, VecMad(m[5], yyyy, VecMad(m[1], xxxx, m[13])));
    vec4_t tw = VecMad(m[11], zzzz, VecMad(m[7], yyyy, VecMad(m[3], xxxx, m[15])));

	tmp0 = VecDiv(VecUnpackLo(tx, ty), VecShuffle(tw, tw, VecShuffleMask(0, 0, 1, 1))); // x0 y0 x1 y1
	tmp1 = VecDiv(VecUnpackHi(tx, ty), VecShuffle(tw, tw, VecShuffleMask(2, 2, 3, 3))); // x2 y2 x3 y3
	VecStore(dst, tmp0);
	VecStore(dst + 2, tmp1);
}

__forceinline inline void Vector3TransformCoord4Homogeneous(const vec4_t* m, const vec4_t* src, vec4_t* dst)
{
    vec4_t src0 = src[0];
    vec4_t src1 = src[1];
    vec4_t src2 = src[2];
    vec4_t src3 = src[3];

    vec4_t tmp0 = VecUnpackLo(src0, src1); // x0 x1 y0 y1
    vec4_t tmp1 = VecUnpackHi(src0, src1); // z0 z1 w0 w1
    vec4_t tmp2 = VecUnpackLo(src2, src3); // x2 x3 y2 y3
    vec4_t tmp3 = VecUnpackHi(src2, src3); // z2 z3 w2 w3

    vec4_t xxxx = VecMoveLH(tmp0, tmp2); // x0 x1 x2 x3
    vec4_t yyyy = VecMoveHL(tmp0, tmp2); // y0 y1 y2 y3
    vec4_t zzzz = VecMoveLH(tmp1, tmp3); // z0 z1 z2 z3

    vec4_t tx = VecMad(m[ 8], zzzz, VecMad(m[4], yyyy, VecMad(m[0], xxxx, m[12])));
    vec4_t ty = VecMad(m[ 9], zzzz, VecMad(m[5], yyyy, VecMad(m[1], xxxx, m[13])));
    vec4_t tz = VecMad(m[10], zzzz, VecMad(m[6], yyyy, VecMad(m[2], xxxx, m[14])));
    vec4_t tw = VecMad(m[11], zzzz, VecMad(m[7], yyyy, VecMad(m[3], xxxx, m[15])));

	tmp0 = VecUnpackLo(tx, ty); // x0 y0 x1 y1
	tmp1 = VecUnpackHi(tx, ty); // x2 y2 x3 y3
	tmp2 = VecUnpackLo(tz, tw); // z0 w0 z1 w1
	tmp3 = VecUnpackHi(tz, tw); // z2 w2 z3 w3

	dst[0] = VecMoveLH(tmp0, tmp2);
	dst[1] = VecMoveHL(tmp0, tmp2);
	dst[2] = VecMoveLH(tmp1, tmp3);
	dst[3] = VecMoveHL(tmp1, tmp3);
}

__forceinline inline vec4_t Vector3TransformCoord(const Matrix& m, vec4_t vec)
{
    vec4_t xxxx = VecShuffle( vec, vec, VecShuffleMask(0,0,0,0) );
    vec4_t yyyy = VecShuffle( vec, vec, VecShuffleMask(1,1,1,1) );
    vec4_t zzzz = VecShuffle( vec, vec, VecShuffleMask(2,2,2,2) );
    vec4_t vvvv = VecAdd(VecAdd(m.r[3], VecMul(zzzz, m.r[2])), VecAdd(VecMul(yyyy, m.r[1]), VecMul(xxxx, m.r[0])));
    vec4_t transformed_vvvv = VecDiv(vvvv, VecShuffle(vvvv, vvvv, VecShuffleMask(3, 3, 3, 3)));
	return VecShuffle(transformed_vvvv, vvvv, VecShuffleMask(0, 1, 2, 3));
}

inline Matrix MatrixIdentity()
{
	Matrix m;
	m.r[0] = Vector4(1, 0, 0, 0);
	m.r[1] = Vector4(0, 1, 0, 0);
	m.r[2] = Vector4(0, 0, 1, 0);
	m.r[3] = Vector4(0, 0, 0, 1);
	return m;
}

inline Matrix MatrixSet(const float* mat)
{
    Matrix m;
    memcpy(&m, mat, sizeof(m));
    return m;
}

inline Matrix operator * ( const Matrix& a, const Matrix& b )
{
	Matrix m;
    vec4_t xxxx, yyyy, zzzz, wwww;
#define MUL_LINE(d,s) \
	xxxx = VecShuffle(s, s, VecShuffleMask(0,0,0,0)); \
	yyyy = VecShuffle(s, s, VecShuffleMask(1,1,1,1)); \
	zzzz = VecShuffle(s, s, VecShuffleMask(2,2,2,2)); \
	wwww = VecShuffle(s, s, VecShuffleMask(3,3,3,3)); \
	d = VecAdd(VecAdd(VecMul(b.r[3], wwww), VecMul(b.r[2], zzzz)), VecAdd(VecMul(b.r[1], yyyy), VecMul(b.r[0], xxxx)));

	MUL_LINE(m.r[0],a.r[0]);
	MUL_LINE(m.r[1],a.r[1]);
	MUL_LINE(m.r[2],a.r[2]);
	MUL_LINE(m.r[3],a.r[3]);
#undef MUL_LINE
	return m;
}

inline Matrix MatrixOrtho(float w, float h, float znear, float zfar)
{
	Matrix m;
	m.r[0] = Vector4(2 / w, 0, 0, 0);
	m.r[1] = Vector4(0, - 2 / h, 0, 0);
	m.r[2] = Vector4(0, 0, 1 / (zfar - znear ), 0);
	m.r[3] = Vector4(-1, 1, -znear / (zfar - znear), 1);
	return m;
}

inline Matrix MatrixTranslation(vec4_t v)
{
	Matrix m = MatrixIdentity();
	m.r[3] = v;
	return m;
}

inline Matrix MatrixPerspective(float fov, float aspect, float znear, float zfar)
{
	float h = 1 / tanf(fov / 2);
	float w = h / aspect;

	Matrix m;
	m.r[0] = Vector4(w, 0, 0, 0);
	m.r[1] = Vector4(0, h, 0, 0);
	m.r[2] = Vector4(0, 0, zfar / (zfar - znear ), 1);
	m.r[3] = Vector4(0, 0, znear * zfar / (znear - zfar), 0);

	return m;
}

inline Matrix MatrixLookAt(vec4_t from, vec4_t at, vec4_t up )
{
    vec4_t zaxis, yaxis, xaxis;

	zaxis = Vector3Normalize(VecSub(at, from));
	xaxis = Vector3Normalize(Vector3Cross(up, zaxis));
	yaxis = Vector3Normalize(Vector3Cross(zaxis, xaxis));

    vec4_t v0 = VecUnpackLo(xaxis, yaxis); // x0 x1 y0 y1
    vec4_t v1 = VecUnpackLo(zaxis, VecZero()); // x2 0 y2 0
    vec4_t v2 = VecUnpackHi(xaxis, yaxis); // z0 z1 ? ?
    vec4_t v3 = VecUnpackHi(zaxis, VecZero()); // z2 0 ? 0

	Matrix m;
	m.r[0] = VecMoveLH(v0, v1);
	m.r[1] = VecMoveHL(v0, v1);
	m.r[2] = VecMoveLH(v2, v3);
	m.r[3] = Vector4(-Vector3Dot(xaxis, from), -Vector3Dot(yaxis, from), -Vector3Dot(zaxis, from), 1);

	return m;
}

inline Matrix MatrixScaling(float scale)
{
	Matrix m = MatrixIdentity();
	m.r[0] = Vector4(scale, 0, 0, 0);
	m.r[1] = Vector4(0, scale, 0, 0);
	m.r[2] = Vector4(0, 0, scale, 0);
	return m;
}

inline Matrix MatrixScaling(float scalex, float scaley, float scalez)
{
	Matrix m = MatrixIdentity();
	m.r[0] = Vector4(scalex, 0, 0, 0);
	m.r[1] = Vector4(0, scaley, 0, 0);
	m.r[2] = Vector4(0, 0, scalez, 0);
	return m;
}
