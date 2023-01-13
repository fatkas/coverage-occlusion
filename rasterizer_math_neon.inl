#pragma once

#include <arm_neon.h>

#define vec4_t float32x4_t
#define vec4i_t int32x4_t

#define ALIGN16 alignas(16)

struct Matrix
{
    vec4_t r[4];
};

inline vec4_t Vector4(float a, float b, float c, float d)
{
    ALIGN16 float data[4] = {a, b, c, d};
    return vld1q_f32(data);
}

inline vec4_t Vector4(float a)
{
	return vdupq_n_f32(a);
}

#define VecShuffleMask(a0,b0,c0,d0) (((d0) << 6) | ((c0) << 4) | ((b0) << 2) | ((a0)))
#define VecShuffle(a,b,mask) __builtin_shufflevector(a, b, (mask) & (0x3), ((mask) >> 2) & 0x3, (((mask) >> 4) & 0x3) + 4, (((mask) >> 6) & 0x3) + 4)

inline vec4_t VecZero()
{
    float32x2_t z = vcreate_f32( 0 );
    return vcombine_f32( z, z );
}

inline vec4_t VecMad(vec4_t a, vec4_t b, vec4_t c)
{
    return vfmaq_f32(c, a, b);
}
inline vec4_t VecMul(vec4_t a, vec4_t b)
{
    return vmulq_f32(a, b);
}
inline vec4_t VecDiv(vec4_t a, vec4_t b)
{
    return vdivq_f32(a, b);
}
inline vec4_t VecRcp(vec4_t a)
{
    return vrecpeq_f32(a);
}
inline vec4_t VecAdd(vec4_t a, vec4_t b)
{
    return vaddq_f32(a, b);
}
inline vec4_t VecSub(vec4_t a, vec4_t b)
{
    return vsubq_f32(a, b);
}
inline vec4_t VecSqrt(vec4_t a)
{
    return vsqrtq_f32(a);
}
inline vec4_t VecMax(vec4_t a, vec4_t b)
{
    return vmaxq_f32(a, b);
}
inline vec4_t VecMin(vec4_t a, vec4_t b)
{
    return vminq_f32(a, b);
}
inline vec4_t VecOr(vec4_t a, vec4_t b)
{
    return vorrq_s32(a, b);
}
inline vec4_t VecAnd(vec4_t a, vec4_t b)
{
    return vandq_s32(a, b);
}
inline vec4_t VecAndNot(vec4_t a, vec4_t b)
{
    return vbicq_s32(b, a);
}
inline int VecMask(vec4_t a)
{
    const int32x4_t shift = {0, 1, 2, 3};
    uint32x4_t tmp = vshrq_n_u32(a, 31);
    return vaddvq_u32(vshlq_u32(tmp, shift));
}
inline vec4_t VecCmpGt(vec4_t a, vec4_t b)
{
    return vcgtq_f32(a, b);
}
inline vec4_t VecCmpLt(vec4_t a, vec4_t b)
{
    return vcltq_f32(a, b);
}
inline vec4_t VecCmpLe(vec4_t a, vec4_t b)
{
    return vcleq_f32(a, b);
}
inline vec4_t VecUnpackLo(vec4_t a, vec4_t b)
{
    return vzip1q_f32(a, b);
}
inline vec4_t VecUnpackHi(vec4_t a, vec4_t b)
{
    return vzip2q_f32(a, b);
}
inline vec4_t VecMoveLH(vec4_t a, vec4_t b)
{
    return vcombine_f32(vget_low_f32(a), vget_low_f32(b));
}
inline vec4_t VecMoveHL(vec4_t a, vec4_t b)
{
    return vzip2q_u64(a, b);
}
inline void VecStore(void* where, vec4_t what)
{
    vst1q_f32((float*)where, what);
}
inline void VecStoreU(void* where, vec4_t what)
{
    vst1q_f32((float*)where, what);
}
inline void VecStoreS(void* where, vec4_t what)
{
    vst1q_lane_f32((float*)where, what, 0);
}
inline vec4_t VecLoad(const void* what)
{
    return vld1q_f32((const float*)what);
}
inline vec4_t VecLoadU(const void* what)
{
    return vld1q_f32((const float*)what);
}

inline vec4i_t Vector4Int(int a, int b, int c, int d)
{
    ALIGN16 int32_t data[4] = {a, b, c, d};
    return vld1q_s32(data);
}
inline vec4i_t VecIntZero()
{
    return vdupq_n_s32(0);
}
inline vec4i_t VecIntAnd(vec4i_t a, vec4i_t b)
{
    return vandq_s32(a, b);
}
inline vec4i_t VecIntOr(vec4i_t a, vec4i_t b)
{
    return vorrq_s32(a, b);
}
inline vec4i_t VecIntXor(vec4i_t a, vec4i_t b)
{
    return veorq_s32(a, b);
}
inline vec4i_t VecFloat2Int(vec4_t a)
{
    return vcvtq_s32_f32(a);
}
inline vec4i_t VecIntLoad(const void* p)
{
    return vld1q_s32((const int32_t *) p);
}
inline void VecIntStore(void* where, vec4i_t a)
{
    vst1q_s32((int32_t *)where, a);
}
inline vec4i_t VecIntCmpEqual(vec4i_t a, vec4i_t b)
{
    return vceqq_s32(a, b);
}
inline int VecIntMask(vec4i_t a)
{
    uint8x16_t input = a;
    uint16x8_t high_bits = vshrq_n_u8(input, 7);
    uint32x4_t paired16 = vsraq_n_u16(high_bits, high_bits, 7);
    uint64x2_t paired32 = vsraq_n_u32(paired16, paired16, 14);
    uint8x16_t paired64 = vsraq_n_u64(paired32, paired32, 28);
    return vgetq_lane_u8(paired64, 0) | ((int) vgetq_lane_u8(paired64, 8) << 8);
}
