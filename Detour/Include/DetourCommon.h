//
// Copyright (c) 2009-2010 Mikko Mononen memon@inside.org
//
// This software is provided 'as-is', without any express or implied
// warranty.  In no event will the authors be held liable for any damages
// arising from the use of this software.
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
// 1. The origin of this software must not be misrepresented; you must not
//    claim that you wrote the original software. If you use this software
//    in a product, an acknowledgment in the product documentation would be
//    appreciated but is not required.
// 2. Altered source versions must be plainly marked as such, and must not be
//    misrepresented as being the original software.
// 3. This notice may not be removed or altered from any source distribution.
//

#ifndef DETOURCOMMON_H
#define DETOURCOMMON_H

#include "DetourMath.h"
#include <stddef.h>

/**
@defgroup detour Detour

Members in this module are used to create, manipulate, and query navigation 
meshes.

@note This is a summary list of members.  Use the index or search 
feature to find minor members.
*/

/// @name General helper functions
/// @{

/// Used to ignore a function parameter.  VS complains about unused parameters
/// and this silences the warning.
/// 用于忽略函数参数。VS 抱怨未使用的参数，这会消除警告。
template<class T> void dtIgnoreUnused(const T&) { }

/// Swaps the values of the two parameters.
///  @param[in,out]	a	Value A
///  @param[in,out]	b	Value B
/// 交换两个参数的值。
template<class T> inline void dtSwap(T& a, T& b) { T t = a; a = b; b = t; }

/// Returns the minimum of two values.
///  @param[in]		a	Value A
///  @param[in]		b	Value B
///  @return The minimum of the two values.
/// 返回两个值中的最小值。
template<class T> inline T dtMin(T a, T b) { return a < b ? a : b; }

/// Returns the maximum of two values.
///  @param[in]		a	Value A
///  @param[in]		b	Value B
///  @return The maximum of the two values.
/// 返回两个值中的最大值。
template<class T> inline T dtMax(T a, T b) { return a > b ? a : b; }

/// Returns the absolute value.
///  @param[in]		a	The value.
///  @return The absolute value of the specified value.
/// 返回绝对值。
template<class T> inline T dtAbs(T a) { return a < 0 ? -a : a; }

/// Returns the square of the value.
///  @param[in]		a	The value.
///  @return The square of the value.
/// 返回值的平方。
template<class T> inline T dtSqr(T a) { return a*a; }

/// Clamps the value to the specified range.
///  @param[in]		v	The value to clamp.
///  @param[in]		mn	The minimum permitted return value.
///  @param[in]		mx	The maximum permitted return value.
///  @return The value, clamped to the specified range.
/// 将值限制在指定范围内。
///  @param[in]		v	要限制的值。
///  @param[in]		mn	允许的最小返回值。
///  @param[in]		mx	允许的最大返回值。
template<class T> inline T dtClamp(T v, T mn, T mx) { return v < mn ? mn : (v > mx ? mx : v); }

/// @}
/// @name Vector helper functions.
/// @{

/// Derives the cross product of two vectors. (@p v1 x @p v2)
///  @param[out]	dest	The cross product. [(x, y, z)]
///  @param[in]		v1		A Vector [(x, y, z)]
///  @param[in]		v2		A vector [(x, y, z)]
/// 导出两个向量的叉积。 （@p v1 x @p v2）
inline void dtVcross(float* dest, const float* v1, const float* v2)
{
	dest[0] = v1[1]*v2[2] - v1[2]*v2[1];
	dest[1] = v1[2]*v2[0] - v1[0]*v2[2];
	dest[2] = v1[0]*v2[1] - v1[1]*v2[0]; 
}

/// Derives the dot product of two vectors. (@p v1 . @p v2)
///  @param[in]		v1	A Vector [(x, y, z)]
///  @param[in]		v2	A vector [(x, y, z)]
/// @return The dot product.
/// 导出两个向量的点积。 （@p v1 。@p v2）
inline float dtVdot(const float* v1, const float* v2)
{
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

/// Performs a scaled vector addition. (@p v1 + (@p v2 * @p s))
///  @param[out]	dest	The result vector. [(x, y, z)]
///  @param[in]		v1		The base vector. [(x, y, z)]
///  @param[in]		v2		The vector to scale and add to @p v1. [(x, y, z)]
///  @param[in]		s		The amount to scale @p v2 by before adding to @p v1.
/// 执行缩放向量加法。 (@p v1 + (@p v2 * @p s))
inline void dtVmad(float* dest, const float* v1, const float* v2, const float s)
{
	dest[0] = v1[0]+v2[0]*s;
	dest[1] = v1[1]+v2[1]*s;
	dest[2] = v1[2]+v2[2]*s;
}

/// Performs a linear interpolation between two vectors. (@p v1 toward @p v2)
///  @param[out]	dest	The result vector. [(x, y, x)]
///  @param[in]		v1		The starting vector.
///  @param[in]		v2		The destination vector.
///	 @param[in]		t		The interpolation factor. [Limits: 0 <= value <= 1.0]
/// 在两个向量之间执行线性插值。 （@p v1 朝向@p v2）
inline void dtVlerp(float* dest, const float* v1, const float* v2, const float t)
{
	dest[0] = v1[0]+(v2[0]-v1[0])*t;
	dest[1] = v1[1]+(v2[1]-v1[1])*t;
	dest[2] = v1[2]+(v2[2]-v1[2])*t;
}

/// Performs a vector addition. (@p v1 + @p v2)
///  @param[out]	dest	The result vector. [(x, y, z)]
///  @param[in]		v1		The base vector. [(x, y, z)]
///  @param[in]		v2		The vector to add to @p v1. [(x, y, z)]
/// 执行向量加法。 （@p v1 + @p v2）
inline void dtVadd(float* dest, const float* v1, const float* v2)
{
	dest[0] = v1[0]+v2[0];
	dest[1] = v1[1]+v2[1];
	dest[2] = v1[2]+v2[2];
}

/// Performs a vector subtraction. (@p v1 - @p v2)
///  @param[out]	dest	The result vector. [(x, y, z)]
///  @param[in]		v1		The base vector. [(x, y, z)]
///  @param[in]		v2		The vector to subtract from @p v1. [(x, y, z)]
/// 执行矢量减法。 （@p v1 - @p v2）
inline void dtVsub(float* dest, const float* v1, const float* v2)
{
	dest[0] = v1[0]-v2[0];
	dest[1] = v1[1]-v2[1];
	dest[2] = v1[2]-v2[2];
}

/// Scales the vector by the specified value. (@p v * @p t)
///  @param[out]	dest	The result vector. [(x, y, z)]
///  @param[in]		v		The vector to scale. [(x, y, z)]
///  @param[in]		t		The scaling factor.
/// 按指定值缩放向量。 (@p v * @p t)
inline void dtVscale(float* dest, const float* v, const float t)
{
	dest[0] = v[0]*t;
	dest[1] = v[1]*t;
	dest[2] = v[2]*t;
}

/// Selects the minimum value of each element from the specified vectors.
///  @param[in,out]	mn	A vector.  (Will be updated with the result.) [(x, y, z)]
///  @param[in]	v	A vector. [(x, y, z)]
/// 从指定向量中选择每个分量的最小值。
inline void dtVmin(float* mn, const float* v)
{
	mn[0] = dtMin(mn[0], v[0]);
	mn[1] = dtMin(mn[1], v[1]);
	mn[2] = dtMin(mn[2], v[2]);
}

/// Selects the maximum value of each element from the specified vectors.
///  @param[in,out]	mx	A vector.  (Will be updated with the result.) [(x, y, z)]
///  @param[in]		v	A vector. [(x, y, z)]
/// 从指定向量中选择每个分量的最大值。
inline void dtVmax(float* mx, const float* v)
{
	mx[0] = dtMax(mx[0], v[0]);
	mx[1] = dtMax(mx[1], v[1]);
	mx[2] = dtMax(mx[2], v[2]);
}

/// Sets the vector elements to the specified values.
///  @param[out]	dest	The result vector. [(x, y, z)]
///  @param[in]		x		The x-value of the vector.
///  @param[in]		y		The y-value of the vector.
///  @param[in]		z		The z-value of the vector.
/// 将向量元素设置为指定值。
inline void dtVset(float* dest, const float x, const float y, const float z)
{
	dest[0] = x; dest[1] = y; dest[2] = z;
}

/// Performs a vector copy.
///  @param[out]	dest	The result. [(x, y, z)]
///  @param[in]		a		The vector to copy. [(x, y, z)]
/// 执行矢量复制。
inline void dtVcopy(float* dest, const float* a)
{
	dest[0] = a[0];
	dest[1] = a[1];
	dest[2] = a[2];
}

/// Derives the scalar length of the vector.
///  @param[in]		v The vector. [(x, y, z)]
/// @return The scalar length of the vector.
/// 返回向量的标量长度。
inline float dtVlen(const float* v)
{
	return dtMathSqrtf(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
}

/// Derives the square of the scalar length of the vector. (len * len)
///  @param[in]		v The vector. [(x, y, z)]
/// @return The square of the scalar length of the vector.
/// 导出向量标量长度的平方。 （长度*长度）
inline float dtVlenSqr(const float* v)
{
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
}

/// Returns the distance between two points.
///  @param[in]		v1	A point. [(x, y, z)]
///  @param[in]		v2	A point. [(x, y, z)]
/// @return The distance between the two points.
/// 返回两点之间的距离。
inline float dtVdist(const float* v1, const float* v2)
{
	const float dx = v2[0] - v1[0];
	const float dy = v2[1] - v1[1];
	const float dz = v2[2] - v1[2];
	return dtMathSqrtf(dx*dx + dy*dy + dz*dz);
}

/// Returns the square of the distance between two points.
///  @param[in]		v1	A point. [(x, y, z)]
///  @param[in]		v2	A point. [(x, y, z)]
/// @return The square of the distance between the two points.
/// 返回两点之间距离的平方。
inline float dtVdistSqr(const float* v1, const float* v2)
{
	const float dx = v2[0] - v1[0];
	const float dy = v2[1] - v1[1];
	const float dz = v2[2] - v1[2];
	return dx*dx + dy*dy + dz*dz;
}

/// Derives the distance between the specified points on the xz-plane.
///  @param[in]		v1	A point. [(x, y, z)]
///  @param[in]		v2	A point. [(x, y, z)]
/// @return The distance between the point on the xz-plane.
///
/// The vectors are projected onto the xz-plane, so the y-values are ignored.
/// 导出 xz 平面上指定点之间的距离。
/// 矢量投影到 xz 平面上，因此忽略 y 值。
inline float dtVdist2D(const float* v1, const float* v2)
{
	const float dx = v2[0] - v1[0];
	const float dz = v2[2] - v1[2];
	return dtMathSqrtf(dx*dx + dz*dz);
}

/// Derives the square of the distance between the specified points on the xz-plane.
///  @param[in]		v1	A point. [(x, y, z)]
///  @param[in]		v2	A point. [(x, y, z)]
/// @return The square of the distance between the point on the xz-plane.
/// 返回 xz 平面上指定点之间距离的平方。
inline float dtVdist2DSqr(const float* v1, const float* v2)
{
	const float dx = v2[0] - v1[0];
	const float dz = v2[2] - v1[2];
	return dx*dx + dz*dz;
}

/// Normalizes the vector.
///  @param[in,out]	v	The vector to normalize. [(x, y, z)]
/// 标准化向量。
inline void dtVnormalize(float* v)
{
	float d = 1.0f / dtMathSqrtf(dtSqr(v[0]) + dtSqr(v[1]) + dtSqr(v[2]));
	v[0] *= d;
	v[1] *= d;
	v[2] *= d;
}

/// Performs a 'sloppy' colocation check of the specified points.
///  @param[in]		p0	A point. [(x, y, z)]
///  @param[in]		p1	A point. [(x, y, z)]
/// @return True if the points are considered to be at the same location.
///
/// Basically, this function will return true if the specified points are
/// close enough to eachother to be considered colocated.
/// 基本上，如果指定点彼此足够接近以被视为共置，则此函数将返回 true。
inline bool dtVequal(const float* p0, const float* p1)
{
	static const float thr = dtSqr(1.0f/16384.0f);
	const float d = dtVdistSqr(p0, p1);
	return d < thr;
}

/// Checks that the specified vector's components are all finite.
///  @param[in]		v	A point. [(x, y, z)]
/// @return True if all of the point's components are finite, i.e. not NaN
/// or any of the infinities.
/// 检查指定向量的分量是否都是有限的。
/// 如果该点的所有分量都是有限的，即不是 NaN 或任何无穷大，则为真。
inline bool dtVisfinite(const float* v)
{
	bool result =
		dtMathIsfinite(v[0]) &&
		dtMathIsfinite(v[1]) &&
		dtMathIsfinite(v[2]);

	return result;
}

/// Checks that the specified vector's 2D components are finite.
///  @param[in]		v	A point. [(x, y, z)]
/// 检查指定向量的 2D 分量是否有限。
inline bool dtVisfinite2D(const float* v)
{
	bool result = dtMathIsfinite(v[0]) && dtMathIsfinite(v[2]);
	return result;
}

/// Derives the dot product of two vectors on the xz-plane. (@p u . @p v)
///  @param[in]		u		A vector [(x, y, z)]
///  @param[in]		v		A vector [(x, y, z)]
/// @return The dot product on the xz-plane.
///
/// The vectors are projected onto the xz-plane, so the y-values are ignored.
/// 导出 xz 平面上两个向量的点积。 （@p u 。@p v）
/// 矢量投影到 xz 平面上，因此忽略 y 值。
inline float dtVdot2D(const float* u, const float* v)
{
	return u[0]*v[0] + u[2]*v[2];
}

/// Derives the xz-plane 2D perp product of the two vectors. (uz*vx - ux*vz)
///  @param[in]		u		The LHV vector [(x, y, z)]
///  @param[in]		v		The RHV vector [(x, y, z)]
/// @return The perp dot product on the xz-plane.
///
/// The vectors are projected onto the xz-plane, so the y-values are ignored.
/// 导出两个向量的 xz 平面 2D perp 乘积(叉积，外积，向量积)。 (uz*vx - ux*vz)
/// 矢量投影到 xz 平面上，因此忽略 y 值。
inline float dtVperp2D(const float* u, const float* v)
{
	return u[2]*v[0] - u[0]*v[2];
}

/// @}
/// @name Computational geometry helper functions.
/// @{

/// Derives the signed xz-plane area of the triangle ABC, or the relationship of line AB to point C.
///  @param[in]		a		Vertex A. [(x, y, z)]
///  @param[in]		b		Vertex B. [(x, y, z)]
///  @param[in]		c		Vertex C. [(x, y, z)]
/// @return The signed xz-plane area of the triangle.
/// 导出三角形 ABC 的带符号 xz 平面面积，或直线 AB 与点 C 的关系。
/// 三角形的有符号 xz 平面面积。(这里是叉积所以是平行四边形形面积，三角形需要除以2)
inline float dtTriArea2D(const float* a, const float* b, const float* c)
{
	const float abx = b[0] - a[0];
	const float abz = b[2] - a[2];
	const float acx = c[0] - a[0];
	const float acz = c[2] - a[2];
	return acx*abz - abx*acz;
}

/// Determines if two axis-aligned bounding boxes overlap.
///  @param[in]		amin	Minimum bounds of box A. [(x, y, z)]
///  @param[in]		amax	Maximum bounds of box A. [(x, y, z)]
///  @param[in]		bmin	Minimum bounds of box B. [(x, y, z)]
///  @param[in]		bmax	Maximum bounds of box B. [(x, y, z)]
/// @return True if the two AABB's overlap.
/// @see dtOverlapBounds
/// 确定两个轴对齐的边界框是否重叠。(判断AABB盒重叠)
inline bool dtOverlapQuantBounds(const unsigned short amin[3], const unsigned short amax[3],
								 const unsigned short bmin[3], const unsigned short bmax[3])
{
	bool overlap = true;
	overlap = (amin[0] > bmax[0] || amax[0] < bmin[0]) ? false : overlap;
	overlap = (amin[1] > bmax[1] || amax[1] < bmin[1]) ? false : overlap;
	overlap = (amin[2] > bmax[2] || amax[2] < bmin[2]) ? false : overlap;
	return overlap;
}

/// Determines if two axis-aligned bounding boxes overlap.
///  @param[in]		amin	Minimum bounds of box A. [(x, y, z)]
///  @param[in]		amax	Maximum bounds of box A. [(x, y, z)]
///  @param[in]		bmin	Minimum bounds of box B. [(x, y, z)]
///  @param[in]		bmax	Maximum bounds of box B. [(x, y, z)]
/// @return True if the two AABB's overlap.
/// @see dtOverlapQuantBounds
/// 确定两个轴对齐的边界框是否重叠。(判断AABB盒重叠)
inline bool dtOverlapBounds(const float* amin, const float* amax,
							const float* bmin, const float* bmax)
{
	bool overlap = true;
	overlap = (amin[0] > bmax[0] || amax[0] < bmin[0]) ? false : overlap;
	overlap = (amin[1] > bmax[1] || amax[1] < bmin[1]) ? false : overlap;
	overlap = (amin[2] > bmax[2] || amax[2] < bmin[2]) ? false : overlap;
	return overlap;
}

/// Derives the closest point on a triangle from the specified reference point.
///  @param[out]	closest	The closest point on the triangle.	
///  @param[in]		p		The reference point from which to test. [(x, y, z)]
///  @param[in]		a		Vertex A of triangle ABC. [(x, y, z)]
///  @param[in]		b		Vertex B of triangle ABC. [(x, y, z)]
///  @param[in]		c		Vertex C of triangle ABC. [(x, y, z)]
/// 从指定参考点得出三角形上最近的点。
void dtClosestPtPointTriangle(float* closest, const float* p,
							  const float* a, const float* b, const float* c);

/// Derives the y-axis height of the closest point on the triangle from the specified reference point.
///  @param[in]		p		The reference point from which to test. [(x, y, z)]
///  @param[in]		a		Vertex A of triangle ABC. [(x, y, z)]
///  @param[in]		b		Vertex B of triangle ABC. [(x, y, z)]
///  @param[in]		c		Vertex C of triangle ABC. [(x, y, z)]
///  @param[out]	h		The resulting height.
/// 导出三角形上距离指定参考点最近的点的 y 轴高度。
bool dtClosestHeightPointTriangle(const float* p, const float* a, const float* b, const float* c, float& h);

// 线段和多边形在xz平面上的投影是否相交
// 并且返回
bool dtIntersectSegmentPoly2D(const float* p0, const float* p1,
							  const float* verts, int nverts,
							  float& tmin, float& tmax,
							  int& segMin, int& segMax);

/// 判断线段是否相交
///  @param[in] ap 线段AB的A点坐标[(x, y, z)]
///  @param[in] aq 线段AB的B点坐标[(x, y, z)]
///  @param[in] bp 线段CD的C点坐标[(x, y, z)]
///  @param[in] bq 线段CD的D点坐标[(x, y, z)]
///  @param[in] s  如果相交赋值焦点到A点的距离
///  @param[in] t  如果相交赋值焦点到C点的距离
/// @return 返回是否相交
bool dtIntersectSegSeg2D(const float* ap, const float* aq,
						 const float* bp, const float* bq,
						 float& s, float& t);

/// Determines if the specified point is inside the convex polygon on the xz-plane.
///  @param[in]		pt		The point to check. [(x, y, z)]
///  @param[in]		verts	The polygon vertices. [(x, y, z) * @p nverts]
///  @param[in]		nverts	The number of vertices. [Limit: >= 3]
/// @return True if the point is inside the polygon.
/// 确定指定点是否位于投影在 xz 平面上的凸多边形内部。
bool dtPointInPolygon(const float* pt, const float* verts, const int nverts);

/// 1. 返回点(pt)是否在多边形内部
/// 2. 输出点(pt)到多边形每条边的最短距离(参数 ed)
/// 3. 输出点(pt)到多边形每条边最短距离点与边起点(逆时针)的距离[0,1]
///  @param pt[in]     点
///  @param verts[in]  多边形顶点列表
///  @param nverts[in] 多边形顶点数
///  @param ed[out]    点(pt)到多边形每条边最短距离的平方
///  @param et[out]    点(pt)到多边形每条边最短距离点到边起点的距离[0,1]
bool dtDistancePtPolyEdgesSqr(const float* pt, const float* verts, const int nverts,
							  float* ed, float* et);

/// 点到线段距离的平方
///  @param pt[in] 需要求距离的点
///  @param p[in]  需要求距离的线段PQ中的p点
///  @param q[in]  需要求距离的线段PQ中的q点
///  @param t[out] 返回点pt与线段PQ距离最短的点到p点的距离[0,1]
/// @return 返回点pt与线段PQ距离最短的平方
float dtDistancePtSegSqr2D(const float* pt, const float* p, const float* q, float& t);

/// Derives the centroid of a convex polygon.
///  @param[out]	tc		The centroid of the polgyon. [(x, y, z)]
///  @param[in]		idx		The polygon indices. [(vertIndex) * @p nidx]
///  @param[in]		nidx	The number of indices in the polygon. [Limit: >= 3]
///  @param[in]		verts	The polygon vertices. [(x, y, z) * vertCount]
/// 导出凸多边形的质心。
///  @param[out]	tc		多边形的质心 [(x, y, z)]
///  @param[in]		idx		多边形索引. [(vertIndex) * @p nidx]
///  @param[in]		nidx	多边形中的索引数. [Limit: >= 3]
///  @param[in]		verts	多边形顶点. [(x, y, z) * vertCount]
void dtCalcPolyCenter(float* tc, const unsigned short* idx, int nidx, const float* verts);

/// Determines if the two convex polygons overlap on the xz-plane.
///  @param[in]		polya		Polygon A vertices.	[(x, y, z) * @p npolya]
///  @param[in]		npolya		The number of vertices in polygon A.
///  @param[in]		polyb		Polygon B vertices.	[(x, y, z) * @p npolyb]
///  @param[in]		npolyb		The number of vertices in polygon B.
/// @return True if the two polygons overlap.
/// 确定两个凸多边形是否在 xz 平面上重叠。
/// 所有顶点都投影到 xz 平面上，因此 y 值将被忽略。
///  @param[in]		polya		多边形A的顶点 [(x, y, z) * @p npolya]
///  @param[in]		npolya		多边形A的顶点数
///  @param[in]		polyb		多边形B的顶点 [(x, y, z) * @p npolyb]
///  @param[in]		npolyb		多边形B的顶点数
/// @return 如果两个多边形重叠，则为 True。
bool dtOverlapPolyPoly2D(const float* polya, const int npolya,
						 const float* polyb, const int npolyb);

/// @}
/// @name Miscellanious functions.
/// @{

//求一个数的临近的较大的2的整数次幂
inline unsigned int dtNextPow2(unsigned int v)
{
	v--;
	v |= v >> 1;
	v |= v >> 2;
	v |= v >> 4;
	v |= v >> 8;
	v |= v >> 16;
	v++;
	return v;
}

// 求一个数取log2的较小的整数
inline unsigned int dtIlog2(unsigned int v)
{
	unsigned int r;
	unsigned int shift;
	r = (v > 0xffff) << 4; v >>= r;
	shift = (v > 0xff) << 3; v >>= shift; r |= shift;
	shift = (v > 0xf) << 2; v >>= shift; r |= shift;
	shift = (v > 0x3) << 1; v >>= shift; r |= shift;
	r |= (v >> 1);
	return r;
}

inline int dtAlign4(int x) { return (x+3) & ~3; }

inline int dtOppositeTile(int side) { return (side+4) & 0x7; }

inline void dtSwapByte(unsigned char* a, unsigned char* b)
{
	unsigned char tmp = *a;
	*a = *b;
	*b = tmp;
}

inline void dtSwapEndian(unsigned short* v)
{
	unsigned char* x = (unsigned char*)v;
	dtSwapByte(x+0, x+1);
}

inline void dtSwapEndian(short* v)
{
	unsigned char* x = (unsigned char*)v;
	dtSwapByte(x+0, x+1);
}

inline void dtSwapEndian(unsigned int* v)
{
	unsigned char* x = (unsigned char*)v;
	dtSwapByte(x+0, x+3); dtSwapByte(x+1, x+2);
}

inline void dtSwapEndian(int* v)
{
	unsigned char* x = (unsigned char*)v;
	dtSwapByte(x+0, x+3); dtSwapByte(x+1, x+2);
}

inline void dtSwapEndian(float* v)
{
	unsigned char* x = (unsigned char*)v;
	dtSwapByte(x+0, x+3); dtSwapByte(x+1, x+2);
}

/// 返回凸多边形中的随机点。
///  @param pts[in]   多边形点列表
///  @param npts[in]  多边形点数量
///  @param areas[in] 三角形区域列表
///  @param s[in]     区域的加权
///  @param t[in]     区域
///  @param out[out]  区域
void dtRandomPointInConvexPoly(const float* pts, const int npts, float* areas,
							   const float s, const float t, float* out);

template<typename TypeToRetrieveAs>
TypeToRetrieveAs* dtGetThenAdvanceBufferPointer(const unsigned char*& buffer, const size_t distanceToAdvance)
{
	TypeToRetrieveAs* returnPointer = reinterpret_cast<TypeToRetrieveAs*>(buffer);
	buffer += distanceToAdvance;
	return returnPointer;
}

template<typename TypeToRetrieveAs>
TypeToRetrieveAs* dtGetThenAdvanceBufferPointer(unsigned char*& buffer, const size_t distanceToAdvance)
{
	TypeToRetrieveAs* returnPointer = reinterpret_cast<TypeToRetrieveAs*>(buffer);
	buffer += distanceToAdvance;
	return returnPointer;
}


/// @}

#endif // DETOURCOMMON_H

///////////////////////////////////////////////////////////////////////////

// This section contains detailed documentation for members that don't have
// a source file. It reduces clutter in the main section of the header.

/**

@fn float dtTriArea2D(const float* a, const float* b, const float* c)
@par

The vertices are projected onto the xz-plane, so the y-values are ignored.
顶点投影到 xz 平面上，因此忽略 y 值。

This is a low cost function than can be used for various purposes.  Its main purpose
is for point/line relationship testing.
这是一个低成本函数，可用于各种目的。 其主要目的是进行点/线关系测试。

In all cases: A value of zero indicates that all vertices are collinear or represent the same point.
(On the xz-plane.)
在所有情况下：值为零表示所有顶点共线或表示同一点。 （在 xz 平面上。）

When used for point/line relationship tests, AB usually represents a line against which
the C point is to be tested.  In this case:
当用于点/线关系测试时，AB通常表示要测试C点的一条线。 在这种情况下：

A positive value indicates that point C is to the left of line AB, looking from A toward B.<br/>
A negative value indicates that point C is to the right of lineAB, looking from A toward B.
正值表示从 A 向 B 方向看，C 点位于 AB 线的左侧。<br/>
负值表示从 A 向 B 方向看，C 点位于 AB 线的右侧。

When used for evaluating a triangle:
当用于评估三角形时：

The absolute value of the return value is two times the area of the triangle when it is
projected onto the xz-plane.
返回值的绝对值是三角形投影到 xz 平面上时面积的两倍。

A positive return value indicates:
正返回值表示：

<ul>
<li>The vertices are wrapped in the normal Detour wrap direction.</li>
<li>The triangle's 3D face normal is in the general up direction.</li>
</ul>
<ul>
<li>顶点以正常的 Detour 包裹方向包裹。</li>
<li>三角形的 3D 面法线总体向上。</li>
</ul>

A negative return value indicates:
负返回值表示：

<ul>
<li>The vertices are reverse wrapped. (Wrapped opposite the normal Detour wrap direction.)</li>
<li>The triangle's 3D face normal is in the general down direction.</li>
</ul>
<ul>
<li>顶点是反向包裹的。 （与正常绕行缠绕方向相反的缠绕方向。）</li>
<li>三角形的 3D 面法线大致向下。</li>
</ul>

*/
