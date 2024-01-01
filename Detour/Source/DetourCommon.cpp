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

#include "DetourCommon.h"
#include "DetourMath.h"

//////////////////////////////////////////////////////////////////////////////////////////

/// 从指定参考点得出三角形上最近的点。
void dtClosestPtPointTriangle(float* closest, const float* p,
							  const float* a, const float* b, const float* c)
{
	// Check if P in vertex region outside A
	// 检查 P 是否位于 A 之外的顶点区域
	float ab[3], ac[3], ap[3];
	dtVsub(ab, b, a);
	dtVsub(ac, c, a);
	dtVsub(ap, p, a);
	float d1 = dtVdot(ab, ap);
	float d2 = dtVdot(ac, ap);
	if (d1 <= 0.0f && d2 <= 0.0f)
	{
		// barycentric coordinates (1,0,0)
		// 重心坐标
		dtVcopy(closest, a);
		return;
	}

	// Check if P in vertex region outside B
	// 检查 P 是否位于 B 之外的顶点区域
	float bp[3];
	dtVsub(bp, p, b);
	float d3 = dtVdot(ab, bp);
	float d4 = dtVdot(ac, bp);
	// CB=(AB-AC) => CB*PB>=0 => (AB-AC)*PB>=0 => AB*PB>=AC*PB => d3>d4
	if (d3 >= 0.0f && d4 <= d3)
	{
		// barycentric coordinates (0,1,0)
		// 重心坐标
		dtVcopy(closest, b);
		return;
	}

	// Check if P in edge region of AB, if so return projection of P onto AB
	// 检查 P 是否在 AB 的边缘区域，如果是，则返回 P 到 AB 的投影
	float vc = d1*d4 - d3*d2;
	if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f)
	{
		// barycentric coordinates (1-v,v,0)
		// 重心坐标
		float v = d1 / (d1 - d3);
		closest[0] = a[0] + v * ab[0];
		closest[1] = a[1] + v * ab[1];
		closest[2] = a[2] + v * ab[2];
		return;
	}

	// Check if P in vertex region outside C
	// 检查 P 是否位于 C 之外的顶点区域
	float cp[3];
	dtVsub(cp, p, c);
	float d5 = dtVdot(ab, cp);
	float d6 = dtVdot(ac, cp);
	if (d6 >= 0.0f && d5 <= d6)
	{
		// barycentric coordinates (0,0,1)
		// 重心坐标
		dtVcopy(closest, c);
		return;
	}

	// Check if P in edge region of AC, if so return projection of P onto AC
	// 检查 P 是否在 AC 的边缘区域，如果是，则返回 P 到 AC 上的投影
	float vb = d5*d2 - d1*d6;
	if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f)
	{
		// barycentric coordinates (1-w,0,w)
		// 重心坐标
		float w = d2 / (d2 - d6);
		closest[0] = a[0] + w * ac[0];
		closest[1] = a[1] + w * ac[1];
		closest[2] = a[2] + w * ac[2];
		return;
	}

	// Check if P in edge region of BC, if so return projection of P onto BC
	// 检查 P 是否在 BC 的边缘区域，如果是，则返回 P 在 BC 上的投影
	float va = d3*d6 - d5*d4;
	if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f)
	{
		// barycentric coordinates (0,1-w,w)
		// 重心坐标
		float w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
		closest[0] = b[0] + w * (c[0] - b[0]);
		closest[1] = b[1] + w * (c[1] - b[1]);
		closest[2] = b[2] + w * (c[2] - b[2]);
		return;
	}

	// P inside face region. Compute Q through its barycentric coordinates (u,v,w)
	// P 内面区域。 通过其重心坐标 (u,v,w) 计算 Q
	float denom = 1.0f / (va + vb + vc);
	float v = vb * denom;
	float w = vc * denom;
	closest[0] = a[0] + ab[0] * v + ac[0] * w;
	closest[1] = a[1] + ab[1] * v + ac[1] * w;
	closest[2] = a[2] + ab[2] * v + ac[2] * w;
}

/// 线段和多边形在xz平面上的投影是否相交
///  @param p0[in]      线段AB的A点坐标
///  @param p1[in]      线段AB的B点坐标
///  @param verts[in]   多边形顶点坐标(凸多边形)
///  @param nverts[in]  多边形顶点数
///  @param tmin[out]   线段与多边形的第一个交点到p0的距离(范围[0, 1])
///  @param tmax[out]   线段与多边形的第二个交点到p0的距离(范围[0, 1])
///  @param segMin[out] 线段与多边形的第一个交点相交的边(范围[0, nverts-1])
///  @param segMax[out] 线段与多边形的第二个交点相交的边(范围[0, nverts-])
/// @return 是否相交
bool dtIntersectSegmentPoly2D(const float* p0, const float* p1,
							  const float* verts, int nverts,
							  float& tmin, float& tmax,
							  int& segMin, int& segMax)
{
	static const float EPS = 0.000001f;

	tmin = 0;
	tmax = 1;
	segMin = -1;
	segMax = -1;

	float dir[3];
	dtVsub(dir, p1, p0); //线段AB向量

	for (int i = 0, j = nverts-1; i < nverts; j=i++)
	{
		float edge[3], diff[3];
		dtVsub(edge, &verts[i*3], &verts[j*3]); //边CD向量
		dtVsub(diff, p0, &verts[j*3]); //CA向量
		const float n = dtVperp2D(edge, diff); //CDxCA
		const float d = dtVperp2D(dir, edge); //-det => ABxCD
		if (fabsf(d) < EPS)
		{
			// S is nearly parallel to this edge
			// S 几乎平行于这条边
			if (n < 0)
				// 线段在多边形外
				return false;
			else
				// 线段在边下方, 或者接近和这条边重合
				continue;
		}
		// 面积之比 即 高之比 即 交点到端点距离与线段之比
		// t 表示的是 （s与edge的交点 到 s的起点）的线段长度/ s的长度
		const float t = n / d; //交点到A(p0)的距离
		if (d < 0)
		{
			// segment S is entering across this edge
			// 线段 S 穿过该边缘进入(AB->CD夹角大于180度)
			// 更新的是 tmin 应该取最大值，因为是靠近，所以交点肯定是越远离起点越好 越远离起点 t值越大
			// 最大的tmin 应该是 edge 与 s的交点在edge上 而非是edge的延长线上。
			if (t > tmin)
			{
				tmin = t;
				segMin = j;
				// S enters after leaving polygon
				// 线段S离开多边形后进入
				if (tmin > tmax)
					return false;
			}
		}
		else
		{
			// segment S is leaving across this edge
			// 线段 S 穿过该边缘离开
			// 更新的是 tmin 应该取最小值 因为是远离，所以交点肯定是越靠近起点越好 越靠近起点 t值越小
			// 最小的tmax 应该是 edge 与 s的交点在edge上 而非是edge的延长线上。
			// 如果 t 一直是 大于1 说明 终点在多边形之内 那么 segMax = -1 tmax = 1
			// 如果 segMax 不等于 -1 ,那么 tmax 肯定是小于1
			if (t < tmax)
			{
				tmax = t;
				segMax = j;
				// S leaves before entering polygon
				// 线段S在进入多边形之前离开
				if (tmax < tmin)
					return false;
			}
		}
	}

	return true;
}

/// 点到线段距离的平方
float dtDistancePtSegSqr2D(const float* pt, const float* p, const float* q, float& t)
{
	float pqx = q[0] - p[0];
	float pqz = q[2] - p[2];
	float dx = pt[0] - p[0];
	float dz = pt[2] - p[2];
	float d = pqx*pqx + pqz*pqz;
	t = pqx*dx + pqz*dz;
	if (d > 0) t /= d; //t: 点pt到直线pq垂线焦点到p点的距离(范围[0,1])
	if (t < 0) t = 0; //距离最短的点在p点
	else if (t > 1) t = 1; //距离最短的点在q点
	//计算pt与线段pq的距离最短的点
	dx = p[0] + t*pqx - pt[0];
	dz = p[2] + t*pqz - pt[2];
	return dx*dx + dz*dz;
}

/// 导出凸多边形的质心。
void dtCalcPolyCenter(float* tc, const unsigned short* idx, int nidx, const float* verts)
{
	tc[0] = 0.0f;
	tc[1] = 0.0f;
	tc[2] = 0.0f;
	for (int j = 0; j < nidx; ++j)
	{
		const float* v = &verts[idx[j]*3];
		tc[0] += v[0];
		tc[1] += v[1];
		tc[2] += v[2];
	}
	const float s = 1.0f / nidx;
	tc[0] *= s;
	tc[1] *= s;
	tc[2] *= s;
}

/// 导出三角形上距离指定参考点最近的点的 y 轴高度。
/*
 * 一个极坐标关于点是否在三角形内的判断方法
The advantage of the method above is that it's very simple to understand so that once you read it you should be able to remember it forever and code it up at any time without having to refer back to anything. It's just - hey the point has to be on the same side of each line as the triangle point that's not in the line. Cake.

Well, there's another method that is also as easy conceptually but executes faster. The downside is there's a little more math involved, but once you see it worked out it should be no problem.

So remember that the three points of the triangle define a plane in space. Pick one of the points and we can consider all other locations on the plane as relative to that point. Let's go with A -- it'll be our origin on the plane. Now what we need are basis vectors so we can give coordinate values to all the locations on the plane. We'll pick the two edges of the triangle that touch A, (C - A) and (B - A). Now we can get to any point on the plane just by starting at A and walking some distance along (C - A) and then from there walking some more in the direction (B - A).

With that in mind we can now describe any point on the plane as

    P = A + u * (C - A) + v * (B - A)
Notice now that if u or v < 0 then we've walked in the wrong direction and must be outside the triangle. Also if u or v > 1 then we've walked too far in a direction and are outside the triangle. Finally if u + v > 1 then we've crossed the edge BC again leaving the triangle.

Given u and v we can easily calculate the point P with the above equation, but how can we go in the reverse direction and calculate u and v from a given point P? Time for some math!

    P = A + u * (C - A) + v * (B - A)       // Original equation
    (P - A) = u * (C - A) + v * (B - A)     // Subtract A from both sides
    v2 = u * v0 + v * v1                    // Substitute v0, v1, v2 for less writing

    // We have two unknowns (u and v) so we need two equations to solve
    // for them.  Dot both sides by v0 to get one and dot both sides by
    // v1 to get a second.
    (v2) . v0 = (u * v0 + v * v1) . v0
    (v2) . v1 = (u * v0 + v * v1) . v1

    // Distribute v0 and v1
    v2 . v0 = u * (v0 . v0) + v * (v1 . v0)
    v2 . v1 = u * (v0 . v1) + v * (v1 . v1)

    // Now we have two equations and two unknowns and can solve one 
    // equation for one variable and substitute into the other.  Or
    // if you're lazy like me, fire up Mathematica and save yourself
    // some handwriting.
    Solve[v2.v0 == {u(v0.v0) + v(v1.v0), v2.v1 == u(v0.v1) + v(v1.v1)}, {u, v}]
    u = ((v1.v1)(v2.v0)-(v1.v0)(v2.v1)) / ((v0.v0)(v1.v1) - (v0.v1)(v1.v0))
    v = ((v0.v0)(v2.v1)-(v0.v1)(v2.v0)) / ((v0.v0)(v1.v1) - (v0.v1)(v1.v0))
Here's an implementation in Flash that you can play with. :)
// Compute vectors
v0 = C - A
v1 = B - A
v2 = P - A

// Compute dot products
dot00 = dot(v0, v0)
dot01 = dot(v0, v1)
dot02 = dot(v0, v2)
dot11 = dot(v1, v1)
dot12 = dot(v1, v2)

// Compute barycentric coordinates
invDenom = 1 / (dot00 * dot11 - dot01 * dot01)
u = (dot11 * dot02 - dot01 * dot12) * invDenom
v = (dot00 * dot12 - dot01 * dot02) * invDenom

// Check if point is in triangle
return (u >= 0) && (v >= 0) && (u + v < 1)
 */
bool dtClosestHeightPointTriangle(const float* p, const float* a, const float* b, const float* c, float& h)
{
	const float EPS = 1e-6f;
	float v0[3], v1[3], v2[3];

	dtVsub(v0, c, a); //AC
	dtVsub(v1, b, a); //AB
	dtVsub(v2, p, a); //AP

	// Compute scaled barycentric coordinates
	// AC和AB在xz平面上的投影的叉乘
	float denom = v0[0] * v1[2] - v0[2] * v1[0]; //ACxAB
	if (fabsf(denom) < EPS)
		return false; //三点共线，不是一个三角形

	float u = v1[2] * v2[0] - v1[0] * v2[2]; //APxAB
	float v = v0[0] * v2[2] - v0[2] * v2[0]; //ACxAB

	if (denom < 0) {
		denom = -denom;
		u = -u;
		v = -v;
	}

	// If point lies inside the triangle, return interpolated ycoord.
	// 如果点位于三角形内，则返回插值 y 坐标。
	if (u >= 0.0f && v >= 0.0f && (u + v) <= denom) {
		h = a[1] + (v0[1] * u + v1[1] * v) / denom;
		return true;
	}
	return false;
}

/// @par
///
/// All points are projected onto the xz-plane, so the y-values are ignored.
/// 所有点都是投影到xz平面，所以y坐标被忽略
/// 1. 遍历多边形每一条边
/// 2. 对每一条边判断点是否在边和z轴形成的梯形中
/// 3. 计算点在边与z轴形成的梯形中出现的次数
/// 4. 如果为次数为单次则在多边形内
bool dtPointInPolygon(const float* pt, const float* verts, const int nverts)
{
	// TODO: Replace pnpoly with triArea2D tests?
	int i, j;
	bool c = false;
	for (i = 0, j = nverts-1; i < nverts; j = i++)
	{
		const float* vi = &verts[i*3];
		const float* vj = &verts[j*3];
		if (((vi[2] > pt[2]) != (vj[2] > pt[2])) &&
			(pt[0] < (vj[0]-vi[0]) * (pt[2]-vi[2]) / (vj[2]-vi[2]) + vi[0]) )
			c = !c;
	}
	return c;
}

/// 1. 返回点(pt)是否在多边形内部
/// 2. 输出点(pt)到多边形每条边的最短距离(参数 ed)
/// 3. 输出点(pt)到多边形每条边最短距离点与边起点(逆时针)的距离[0,1]
bool dtDistancePtPolyEdgesSqr(const float* pt, const float* verts, const int nverts,
							  float* ed, float* et)
{
	// TODO: Replace pnpoly with triArea2D tests?
	int i, j;
	bool c = false;
	for (i = 0, j = nverts-1; i < nverts; j = i++)
	{
		const float* vi = &verts[i*3];
		const float* vj = &verts[j*3];
		if (((vi[2] > pt[2]) != (vj[2] > pt[2])) &&
			(pt[0] < (vj[0]-vi[0]) * (pt[2]-vi[2]) / (vj[2]-vi[2]) + vi[0]) )
			c = !c;
		ed[j] = dtDistancePtSegSqr2D(pt, vj, vi, et[j]);
	}
	return c;
}

//计算多边形在某条轴上的相对投影范围
static void projectPoly(const float* axis, const float* poly, const int npoly,
						float& rmin, float& rmax)
{
	// 求最大和最小的点积值 相当于 多边形在 轴上的投影范围（真正的这个范围需要除以|axis|，因为都乘了无所谓）
	rmin = rmax = dtVdot2D(axis, &poly[0]);
	for (int i = 1; i < npoly; ++i)
	{
		const float d = dtVdot2D(axis, &poly[i*3]);
		rmin = dtMin(rmin, d);
		rmax = dtMax(rmax, d);
	}
}

// 计算是否重合
inline bool overlapRange(const float amin, const float amax,
						 const float bmin, const float bmax,
						 const float eps)
{
	return ((amin+eps) > bmax || (amax-eps) < bmin) ? false : true;
}

/// @par
///
/// All vertices are projected onto the xz-plane, so the y-values are ignored.
/// 所有顶点都投影到 xz 平面上，因此 y 值将被忽略。
bool dtOverlapPolyPoly2D(const float* polya, const int npolya,
						 const float* polyb, const int npolyb)
{
	const float eps = 1e-4f;
	
	for (int i = 0, j = npolya-1; i < npolya; j=i++)
	{
		const float* va = &polya[j*3];
		const float* vb = &polya[i*3];
		// 边AB在xz平面上投影的法线(与边垂直的向量)，作为分离轴
		const float n[3] = { vb[2]-va[2], 0, -(vb[0]-va[0]) };
		float amin,amax,bmin,bmax;
		projectPoly(n, polya, npolya, amin,amax); //多边形A在n上的投影*|n|
		projectPoly(n, polyb, npolyb, bmin,bmax); //多边形B在n上的投影*|n|
		if (!overlapRange(amin,amax, bmin,bmax, eps))
		{
			// Found separating axis
			return false;
		}
	}
	for (int i = 0, j = npolyb-1; i < npolyb; j=i++)
	{
		const float* va = &polyb[j*3];
		const float* vb = &polyb[i*3];
		// 边AB在xz平面上投影的法线(与边垂直的向量)，作为分离轴
		const float n[3] = { vb[2]-va[2], 0, -(vb[0]-va[0]) };
		float amin,amax,bmin,bmax;
		projectPoly(n, polya, npolya, amin,amax); //多边形A在n上的投影*|n|
		projectPoly(n, polyb, npolyb, bmin,bmax); //多边形B在n上的投影*|n|
		if (!overlapRange(amin,amax, bmin,bmax, eps))
		{
			// Found separating axis
			return false;
		}
	}
	//所有的边的xz平面法线上的投影都重合，代表两个多边形重合
	return true;
}

// Returns a random point in a convex polygon.
// Adapted from Graphics Gems article.
// 返回凸多边形中的随机点。
void dtRandomPointInConvexPoly(const float* pts, const int npts, float* areas,
							   const float s, const float t, float* out)
{
	// Calc triangle araes
	float areasum = 0.0f;
	for (int i = 2; i < npts; i++) {
		//计算三角形
		areas[i] = dtTriArea2D(&pts[0], &pts[(i-1)*3], &pts[i*3]);
		areasum += dtMax(0.001f, areas[i]);
	}
	// Find sub triangle weighted by area.
	// 找到按面积加权的子三角形。
	const float thr = s*areasum;
	float acc = 0.0f;
	float u = 1.0f;
	int tri = npts - 1;
	for (int i = 2; i < npts; i++) {
		const float dacc = areas[i];
		if (thr >= acc && thr < (acc+dacc))
		{
			u = (thr - acc) / dacc;
			tri = i;
			break;
		}
		acc += dacc;
	}

	float v = dtMathSqrtf(t);

	const float a = 1 - v;
	const float b = (1 - u) * v;
	const float c = u * v;
	const float* pa = &pts[0];
	const float* pb = &pts[(tri-1)*3];
	const float* pc = &pts[tri*3];

	out[0] = a*pa[0] + b*pb[0] + c*pc[0];
	out[1] = a*pa[1] + b*pb[1] + c*pc[1];
	out[2] = a*pa[2] + b*pb[2] + c*pc[2];
}

inline float vperpXZ(const float* a, const float* b) { return a[0]*b[2] - a[2]*b[0]; }

/// 判断线段是否相交
bool dtIntersectSegSeg2D(const float* ap, const float* aq,
						 const float* bp, const float* bq,
						 float& s, float& t)
{
	float u[3], v[3], w[3];
	dtVsub(u,aq,ap); //AB向量
	dtVsub(v,bq,bp); //CD向量
	dtVsub(w,ap,bp); //CA向量
	float d = vperpXZ(u,v); //-det => ABxCD
	if (fabsf(d) < 1e-6f) return false;
	s = vperpXZ(v,w) / d; //(CDxCA)/(ABxCD)
	t = vperpXZ(u,w) / d; //(ABxCA)/(ABxCD)
	return true;
}

