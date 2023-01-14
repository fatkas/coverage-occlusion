#pragma once

#include <xmmintrin.h>
#include <emmintrin.h>
#include <smmintrin.h>

#define vec4_t __m128
#define vec4i_t __m128i

#define ALIGN16 alignas(16)

struct Matrix
{
    vec4_t r[4];
};

inline vec4_t Vector4( float a, float b, float c, float d )
{
	return _mm_setr_ps(a,b,c,d);
}

inline vec4_t Vector4( float a )
{
	return _mm_set1_ps(a);
}

#define VecShuffleMask(a0,b0,c0,d0) _MM_SHUFFLE(d0,c0,b0,a0)
#define VecShuffle(a,b,mask) _mm_shuffle_ps((a),(b),mask)
#define VecZero() _mm_setzero_ps()
#define VecMad(a,b,c) _mm_add_ps(c,_mm_mul_ps(a,b))
#define VecMul(a,b) _mm_mul_ps(a,b)
#define VecDiv(a,b) _mm_div_ps(a,b)
#define VecRcp(a) _mm_rcp_ps(a)
#define VecAdd(a,b) _mm_add_ps(a,b)
#define VecSub(a,b) _mm_sub_ps(a,b)
#define VecSqrt(a) _mm_sqrt_ps(a)
#define VecMax(a,b) _mm_max_ps(a,b)
#define VecMin(a,b) _mm_min_ps(a,b)
#define VecOr(a,b) _mm_or_ps(a,b)
#define VecAnd(a,b) _mm_and_ps(a,b)
#define VecAndNot(a,b) _mm_andnot_ps(a,b)
#define VecMask(a) _mm_movemask_ps(a)
#define VecCmpGt(a,b) _mm_cmpgt_ps(a,b)
#define VecCmpLt(a,b) _mm_cmplt_ps(a,b)
#define VecCmpLe(a,b) _mm_cmple_ps(a,b)
#define VecUnpackLo(a,b) _mm_unpacklo_ps(a,b)
#define VecUnpackHi(a,b) _mm_unpackhi_ps(a,b)
#define VecMoveLH(a,b) _mm_movelh_ps(a,b)
#define VecMoveHL(a,b) _mm_movehl_ps(b,a)
#define VecStore(where,what) _mm_store_ps((float*)(where),(what))
#define VecStoreU(where,what) _mm_storeu_ps((float*)(where),(what))
#define VecStoreS(where,what) _mm_store_ss((float*)(where),(what))
#define VecLoad(what) _mm_load_ps((const float*)(what))
#define VecLoadU(what) _mm_loadu_ps((const float*)(what))

#define Vector4Int(a,b,c,d) _mm_setr_epi32(a,b,c,d)
#define VecIntShuffle(a,mask) _mm_shuffle_epi32((a),mask)
#define VecIntZero() _mm_setzero_si128()
#define VecIntAnd(a,b) _mm_and_si128(a,b)
#define VecIntOr(a,b) _mm_or_si128(a,b)
#define VecIntAndNot(a,b) _mm_andnot_si128(a,b)
#define VecIntXor(a,b) _mm_xor_si128(a,b)
#define VecIntMul(a,b) _mm_mul_epu32(a,b)
#define VecIntAdd(a,b) _mm_add_epi32(a,b)
#define VecIntSub(a,b) _mm_sub_epi32(a,b)
#define VecIntMax(a,b) _mm_max_epi16(a,b)
#define VecIntMin(a,b) _mm_min_epi16(a,b)
#define VecFloat2Int(a) _mm_cvttps_epi32(a)
#define VecIntLoad(where) _mm_load_si128( (const __m128i*)(where))
#define VecIntStore(where,what) _mm_store_si128( (__m128i*)(where),(what))
#define VecIntSRL64(a,b) _mm_srl_epi64((a),(b))
#define VecIntSRL32(a,b) _mm_srl_epi32((a),(b))
#define VecIntSLL64(a,b) _mm_sll_epi64((a),(b))
#define VecIntSLL32(a,b) _mm_sll_epi32((a),(b))
#define VecIntCmpEqual(a,b) _mm_cmpeq_epi32((a),(b))
#define VecIntUnpackLo(a,b) _mm_unpacklo_epi16((a),(b))
#define VecIntUnpackHi(a,b) _mm_unpackhi_epi16((a),(b))
#define VecIntMask(a) _mm_movemask_epi8(a)
#define VecInt2Float(a) _mm_cvtepi32_ps(a)
#define VecIntPack16(a,b) _mm_packus_epi32(a,b)
