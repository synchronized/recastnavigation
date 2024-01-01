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

#include <math.h>
#include "Recast.h"
#include "RecastAlloc.h"
#include "RecastAssert.h"

/// Check whether two bounding boxes overlap
///
/// @param[in]	aMin	Min axis extents of bounding box A
/// @param[in]	aMax	Max axis extents of bounding box A
/// @param[in]	bMin	Min axis extents of bounding box B
/// @param[in]	bMax	Max axis extents of bounding box B
/// @returns true if the two bounding boxes overlap.  False otherwise.
static bool overlapBounds(const float* aMin, const float* aMax, const float* bMin, const float* bMax)
{
	return
		aMin[0] <= bMax[0] && aMax[0] >= bMin[0] &&
		aMin[1] <= bMax[1] && aMax[1] >= bMin[1] &&
		aMin[2] <= bMax[2] && aMax[2] >= bMin[2];
}

/// Allocates a new span in the heightfield.
/// Use a memory pool and free list to minimize actual allocations.
///
/// @param[in]	heightfield		The heightfield
/// @returns A pointer to the allocated or re-used span memory.
static rcSpan* allocSpan(rcHeightfield& heightfield)
{
	// If necessary, allocate new page and update the freelist.
	if (heightfield.freelist == NULL || heightfield.freelist->next == NULL)
	{
		// Create new page.
		// Allocate memory for the new pool.
		rcSpanPool* spanPool = (rcSpanPool*)rcAlloc(sizeof(rcSpanPool), RC_ALLOC_PERM);
		if (spanPool == NULL)
		{
			return NULL;
		}

		// Add the pool into the list of pools.
		spanPool->next = heightfield.pools;
		heightfield.pools = spanPool;

		// Add new spans to the free list.
		rcSpan* freeList = heightfield.freelist;
		rcSpan* head = &spanPool->items[0];
		rcSpan* it = &spanPool->items[RC_SPANS_PER_POOL];
		do
		{
			--it;
			it->next = freeList;
			freeList = it;
		}
		while (it != head);
		heightfield.freelist = it;
	}

	// Pop item from the front of the free list.
	rcSpan* newSpan = heightfield.freelist;
	heightfield.freelist = heightfield.freelist->next;
	return newSpan;
}

/// Releases the memory used by the span back to the heightfield, so it can be re-used for new spans.
/// @param[in]	heightfield		The heightfield.
/// @param[in]	span	A pointer to the span to free
static void freeSpan(rcHeightfield& heightfield, rcSpan* span)
{
	if (span == NULL)
	{
		return;
	}
	// Add the span to the front of the free list.
	span->next = heightfield.freelist;
	heightfield.freelist = span;
}

/// Adds a span to the heightfield.  If the new span overlaps existing spans,
/// it will merge the new span with the existing ones.
///
/// @param[in]	heightfield					Heightfield to add spans to
/// @param[in]	x					The new span's column cell x index
/// @param[in]	z					The new span's column cell z index
/// @param[in]	min					The new span's minimum cell index
/// @param[in]	max					The new span's maximum cell index
/// @param[in]	areaID				The new span's area type ID
/// @param[in]	flagMergeThreshold	How close two spans maximum extents need to be to merge area type IDs
/// 向高度场添加跨度。 如果新的 Span 与现有的 Span 重叠，则会将新的 Span 与现有的 Span 合并。
static bool addSpan(rcHeightfield& heightfield,
                    const int x, const int z,
                    const unsigned short min, const unsigned short max,
                    const unsigned char areaID, const int flagMergeThreshold)
{
	// Create the new span.
	rcSpan* newSpan = allocSpan(heightfield);
	if (newSpan == NULL)
	{
		return false;
	}
	newSpan->smin = min;
	newSpan->smax = max;
	newSpan->area = areaID;
	newSpan->next = NULL;

	const int columnIndex = x + z * heightfield.width;
	rcSpan* previousSpan = NULL;
	rcSpan* currentSpan = heightfield.spans[columnIndex];

	// Insert the new span, possibly merging it with existing spans.
	// 插入新的span，可能会将其与现有的span合并。
	while (currentSpan != NULL)
	{
		if (currentSpan->smin > newSpan->smax)
		{
			// Current span is completely after the new span, break.
			// 当前span在新span之上。
			break;
		}

		if (currentSpan->smax < newSpan->smin)
		{
			// Current span is completely before the new span.  Keep going.
			// 当前span在新span之下, 查找下一个
			previousSpan = currentSpan;
			currentSpan = currentSpan->next;
		}
		else
		{
			// The new span overlaps with an existing span.  Merge them.
			// 新span与现有span重叠。 合并它们。
			if (currentSpan->smin < newSpan->smin)
			{
				newSpan->smin = currentSpan->smin;
			}
			if (currentSpan->smax > newSpan->smax)
			{
				newSpan->smax = currentSpan->smax;
			}

			// Merge flags.
			if (rcAbs((int)newSpan->smax - (int)currentSpan->smax) <= flagMergeThreshold)
			{
				// Higher area ID numbers indicate higher resolution priority.
				// area ID 编号越高表示分辨率优先级越高。
				newSpan->area = rcMax(newSpan->area, currentSpan->area);
			}

			// Remove the current span since it's now merged with newSpan.
			// Keep going because there might be other overlapping spans that also need to be merged.
			// 删除当前跨度，因为它现已与 newSpan 合并。
			// 继续，因为可能还有其他重叠的跨度也需要合并。
			rcSpan* next = currentSpan->next;
			freeSpan(heightfield, currentSpan);
			if (previousSpan)
			{
				previousSpan->next = next;
			}
			else
			{
				heightfield.spans[columnIndex] = next;
			}
			currentSpan = next;
		}
	}

	// Insert new span after prev
	if (previousSpan != NULL)
	{
		newSpan->next = previousSpan->next;
		previousSpan->next = newSpan;
	}
	else
	{
		// This span should go before the others in the list
		// 该span应位于列表中其他span之前
		newSpan->next = heightfield.spans[columnIndex];
		heightfield.spans[columnIndex] = newSpan;
	}

	return true;
}

bool rcAddSpan(rcContext* context, rcHeightfield& heightfield,
               const int x, const int z,
               const unsigned short spanMin, const unsigned short spanMax,
               const unsigned char areaID, const int flagMergeThreshold)
{
	rcAssert(context);

	if (!addSpan(heightfield, x, z, spanMin, spanMax, areaID, flagMergeThreshold))
	{
		context->log(RC_LOG_ERROR, "rcAddSpan: Out of memory.");
		return false;
	}

	return true;
}

enum rcAxis
{
	RC_AXIS_X = 0,
	RC_AXIS_Y = 1,
	RC_AXIS_Z = 2
};

/// Divides a convex polygon of max 12 vertices into two convex polygons
/// across a separating axis.
///
/// @param[in]	inVerts			The input polygon vertices
/// @param[in]	inVertsCount	The number of input polygon vertices
/// @param[out]	outVerts1		Resulting polygon 1's vertices
/// @param[out]	outVerts1Count	The number of resulting polygon 1 vertices
/// @param[out]	outVerts2		Resulting polygon 2's vertices
/// @param[out]	outVerts2Count	The number of resulting polygon 2 vertices
/// @param[in]	axisOffset		THe offset along the specified axis
/// @param[in]	axis			The separating axis
/// 将最多 12 个顶点的凸多边形沿分离轴划分为两个凸多边形。
static void dividePoly(const float* inVerts, int inVertsCount,
                       float* outVerts1, int* outVerts1Count,
                       float* outVerts2, int* outVerts2Count,
                       float axisOffset, rcAxis axis)
{
	rcAssert(inVertsCount <= 12);

	// How far positive or negative away from the separating axis is each vertex.
	// 取出多边形每个点相对于分割线的距离
	float inVertAxisDelta[12];
	for (int inVert = 0; inVert < inVertsCount; ++inVert)
	{
		inVertAxisDelta[inVert] = axisOffset - inVerts[inVert * 3 + axis];
	}

	int poly1Vert = 0;
	int poly2Vert = 0;
	// 遍历所有的边与切割线的相交关系，inVertA和inVertB两个点即为一条边
	for (int inVertA = 0, inVertB = inVertsCount - 1; inVertA < inVertsCount; inVertB = inVertA, ++inVertA)
	{
		// If the two vertices are on the same side of the separating axis
		// 如果两个顶点位于分离轴的同一侧
		bool sameSide = (inVertAxisDelta[inVertA] >= 0) == (inVertAxisDelta[inVertB] >= 0);

		// 如果边上的2个点在切割线两边，则新产生一个顶点，即切割点，切割点属于被切开的两个多边形，所以outVerts1和outVerts2中都要加入此切割点
		if (!sameSide)
		{
			// 计算切割比（由于inVerts和inVertb两个点在分割线两侧，所以inVertAxisDelta[inVertB] - inVertAxisDelta[inVertA]为到分割线距离的和）
			float s = inVertAxisDelta[inVertB] / (inVertAxisDelta[inVertB] - inVertAxisDelta[inVertA]);
			// 根据切割比，计算切割点的三维坐标
			outVerts1[poly1Vert * 3 + 0] = inVerts[inVertB * 3 + 0] + (inVerts[inVertA * 3 + 0] - inVerts[inVertB * 3 + 0]) * s;
			outVerts1[poly1Vert * 3 + 1] = inVerts[inVertB * 3 + 1] + (inVerts[inVertA * 3 + 1] - inVerts[inVertB * 3 + 1]) * s;
			outVerts1[poly1Vert * 3 + 2] = inVerts[inVertB * 3 + 2] + (inVerts[inVertA * 3 + 2] - inVerts[inVertB * 3 + 2]) * s;
			rcVcopy(&outVerts2[poly2Vert * 3], &outVerts1[poly1Vert * 3]);
			poly1Vert++;
			poly2Vert++;

			// add the inVertA point to the right polygon. Do NOT add points that are on the dividing line
			// since these were already added above
			// 将 inVertA 点添加到正确的多边形。 不要添加分界线上的点，因为这些点已经在上面添加了
			if (inVertAxisDelta[inVertA] > 0)
			{
				rcVcopy(&outVerts1[poly1Vert * 3], &inVerts[inVertA * 3]);
				poly1Vert++;
			}
			else if (inVertAxisDelta[inVertA] < 0)
			{
				rcVcopy(&outVerts2[poly2Vert * 3], &inVerts[inVertA * 3]);
				poly2Vert++;
			}
		}
		else
		{
			// add the inVertA point to the right polygon. Addition is done even for points on the dividing line
			// 将 inVertA 点添加到正确的多边形。 即使分界线上的点也进行添加
			if (inVertAxisDelta[inVertA] >= 0)
			{
				rcVcopy(&outVerts1[poly1Vert * 3], &inVerts[inVertA * 3]);
				poly1Vert++;
				if (inVertAxisDelta[inVertA] != 0)
				{
					continue;
				}
			}
			rcVcopy(&outVerts2[poly2Vert * 3], &inVerts[inVertA * 3]);
			poly2Vert++;
		}
	}

	*outVerts1Count = poly1Vert;
	*outVerts2Count = poly2Vert;
}

///	Rasterize a single triangle to the heightfield.
///
///	This code is extremely hot, so much care should be given to maintaining maximum perf here.
///
/// @param[in] 	v0					Triangle vertex 0
/// @param[in] 	v1					Triangle vertex 1
/// @param[in] 	v2					Triangle vertex 2
/// @param[in] 	areaID				The area ID to assign to the rasterized spans
/// @param[in] 	heightfield			Heightfield to rasterize into
/// @param[in] 	heightfieldBBMin	The min extents of the heightfield bounding box
/// @param[in] 	heightfieldBBMax	The max extents of the heightfield bounding box
/// @param[in] 	cellSize			The x and z axis size of a voxel in the heightfield
/// @param[in] 	inverseCellSize		1 / cellSize
/// @param[in] 	inverseCellHeight	1 / cellHeight
/// @param[in] 	flagMergeThreshold	The threshold in which area flags will be merged
/// @returns true if the operation completes successfully.  false if there was an error adding spans to the heightfield.
///	将单个三角形光栅化到高度场。
///
///	这段代码非常热门，因此应该非常小心地保持这里的最大性能。
static bool rasterizeTri(const float* v0, const float* v1, const float* v2,
                         const unsigned char areaID, rcHeightfield& heightfield,
                         const float* heightfieldBBMin, const float* heightfieldBBMax,
                         const float cellSize, const float inverseCellSize, const float inverseCellHeight,
                         const int flagMergeThreshold)
{
	// Calculate the bounding box of the triangle.
	// 计算三角形碰撞盒
	float triBBMin[3];
	rcVcopy(triBBMin, v0);
	rcVmin(triBBMin, v1);
	rcVmin(triBBMin, v2);

	float triBBMax[3];
	rcVcopy(triBBMax, v0);
	rcVmax(triBBMax, v1);
	rcVmax(triBBMax, v2);

	// If the triangle does not touch the bounding box of the heightfield, skip the triangle.
	// 如果三角形与高度场没有重叠，则跳过该三角形
	if (!overlapBounds(triBBMin, triBBMax, heightfieldBBMin, heightfieldBBMax))
	{
		return true;
	}

	const int w = heightfield.width;
	const int h = heightfield.height;
	const float by = heightfieldBBMax[1] - heightfieldBBMin[1];

	// Calculate the footprint of the triangle on the grid's z-axis
	// 计算三角形z轴格子的覆盖区
	int z0 = (int)((triBBMin[2] - heightfieldBBMin[2]) * inverseCellSize);
	int z1 = (int)((triBBMax[2] - heightfieldBBMin[2]) * inverseCellSize);

	// use -1 rather than 0 to cut the polygon properly at the start of the tile
	// 使用 -1 而不是 0 在图块的开头正确切割多边形
	z0 = rcClamp(z0, -1, h - 1);
	z1 = rcClamp(z1, 0, h - 1);

	// Clip the triangle into all grid cells it touches.
	// 将三角形夹入其接触的所有网格单元中。
	// 定义一块连续内存 buf[7*3*4]
	// 7的意义:对三角形切割时，切出来的一个格子，最多可变成一个7边形
	// 3的意义:存的是顶点，所以3个float， x y z
    // 4的意义:下面的in inrow p1 p2四份等长数据，这些数据里存的是多边形，即顶点坐标
	float buf[7 * 3 * 4];
	float* in = buf;
	float* inRow = buf + 7 * 3;
	float* p1 = inRow + 7 * 3;
	float* p2 = p1 + 7 * 3;

	rcVcopy(&in[0], v0);
	rcVcopy(&in[1 * 3], v1);
	rcVcopy(&in[2 * 3], v2);
	int nvRow;
	int nvIn = 3;

	for (int z = z0; z <= z1; ++z)
	{
		// Clip polygon to row. Store the remaining polygon as well
		// 将多边形剪辑到行(z轴)。 也存储剩余的多边形
		const float cellZ = heightfieldBBMin[2] + (float)z * cellSize;
		// 多边形切割，将in以cellZ+cellSize为切割线,划分为inRow和p1两个多边形（z轴递进横向切割）
		dividePoly(in, nvIn, inRow, &nvRow, p1, &nvIn, cellZ + cellSize, RC_AXIS_Z);
		// inRow在本次循环内按x0到x1切割完成，另一半p1下次循环在进行切分，所以p1和in交换下 下次再切分in即p1
		rcSwap(in, p1);

		//这次要处理的这一半多边形没切到东西，可能切到了in（多边形）的一个顶点或者与y轴平行的一条边，顶点和边没必要体素化，所以continue
		if (nvRow < 3)
		{
			continue;
		}
		if (z < 0)
		{
			continue;
		}

		// find X-axis bounds of the row
		// 查找x轴在行中的碰撞
		float minX = inRow[0];
		float maxX = inRow[0];
		// 遍历多边形inRow的顶点，找到最小x和最大x
		for (int vert = 1; vert < nvRow; ++vert)
		{
			if (minX > inRow[vert * 3])
			{
				minX = inRow[vert * 3];
			}
			if (maxX < inRow[vert * 3])
			{
				maxX = inRow[vert * 3];
			}
		}
		int x0 = (int)((minX - heightfieldBBMin[0]) * inverseCellSize);
		int x1 = (int)((maxX - heightfieldBBMin[0]) * inverseCellSize);
		if (x1 < 0 || x0 >= w)
		{
			continue;
		}
		x0 = rcClamp(x0, -1, w - 1);
		x1 = rcClamp(x1, 0, w - 1);

		int nv;
		int nv2 = nvRow;

		for (int x = x0; x <= x1; ++x)
		{
			// Clip polygon to column. store the remaining polygon as well
			// 将多边形剪辑到列(x轴,应该是单元格)。 也存储剩余的多边形
			// 再竖着切割多边形inRow
			const float cx = heightfieldBBMin[0] + (float)x * cellSize;
			dividePoly(inRow, nv2, p1, &nv, p2, &nv2, cx + cellSize, RC_AXIS_X);
			// p1为切割线左边的多边形，p2为切割线右边的多边形，此时p1已经被切割完毕（已经是一个格子内的多边形了）
            // p2还没切割完成，与inrow交换，下次遍历切割
			rcSwap(inRow, p2);

			if (nv < 3)
			{
				continue;
			}
			if (x < 0)
			{
				continue;
			}

			// Calculate min and max of the span.
			// smin：多边形p1的最低高度
            // smax：多边形p1的最高高度
			float spanMin = p1[1];
			float spanMax = p1[1];
			for (int vert = 1; vert < nv; ++vert)
			{
				spanMin = rcMin(spanMin, p1[vert * 3 + 1]);
				spanMax = rcMax(spanMax, p1[vert * 3 + 1]);
			}
			spanMin -= heightfieldBBMin[1];
			spanMax -= heightfieldBBMin[1];

			// Skip the span if it's completely outside the heightfield bounding box
			// 如果span在高度场 bbox 之外，则跳过span
			if (spanMax < 0.0f)
			{
				continue;
			}
			if (spanMin > by)
			{
				continue;
			}

			// Clamp the span to the heightfield bounding box.
			// 将span限制在高度域 bbox 上
			if (spanMin < 0.0f)
			{
				spanMin = 0;
			}
			if (spanMax > by)
			{
				spanMax = by;
			}

			// Snap the span to the heightfield height grid.
			// 包围盒最低体素ismin，最高体素ismax
			unsigned short spanMinCellIndex = (unsigned short)rcClamp((int)floorf(spanMin * inverseCellHeight), 0, RC_SPAN_MAX_HEIGHT);
			unsigned short spanMaxCellIndex = (unsigned short)rcClamp((int)ceilf(spanMax * inverseCellHeight), (int)spanMinCellIndex + 1, RC_SPAN_MAX_HEIGHT);

			//把切分下来的span合并到到高度场中
			if (!addSpan(heightfield, x, z, spanMinCellIndex, spanMaxCellIndex, areaID, flagMergeThreshold))
			{
				return false;
			}
		}
	}

	return true;
}

bool rcRasterizeTriangle(rcContext* context,
                         const float* v0, const float* v1, const float* v2,
                         const unsigned char areaID, rcHeightfield& heightfield, const int flagMergeThreshold)
{
	rcAssert(context != NULL);

	rcScopedTimer timer(context, RC_TIMER_RASTERIZE_TRIANGLES);

	// Rasterize the single triangle.
	const float inverseCellSize = 1.0f / heightfield.cs;
	const float inverseCellHeight = 1.0f / heightfield.ch;
	if (!rasterizeTri(v0, v1, v2, areaID, heightfield, heightfield.bmin, heightfield.bmax, heightfield.cs, inverseCellSize, inverseCellHeight, flagMergeThreshold))
	{
		context->log(RC_LOG_ERROR, "rcRasterizeTriangle: Out of memory.");
		return false;
	}

	return true;
}

bool rcRasterizeTriangles(rcContext* context,
                          const float* verts, const int /*nv*/,
                          const int* tris, const unsigned char* triAreaIDs, const int numTris,
                          rcHeightfield& heightfield, const int flagMergeThreshold)
{
	rcAssert(context != NULL);

	rcScopedTimer timer(context, RC_TIMER_RASTERIZE_TRIANGLES);

	// Rasterize the triangles.
	const float inverseCellSize = 1.0f / heightfield.cs;
	const float inverseCellHeight = 1.0f / heightfield.ch;
	for (int triIndex = 0; triIndex < numTris; ++triIndex)
	{
		const float* v0 = &verts[tris[triIndex * 3 + 0] * 3];
		const float* v1 = &verts[tris[triIndex * 3 + 1] * 3];
		const float* v2 = &verts[tris[triIndex * 3 + 2] * 3];
		if (!rasterizeTri(v0, v1, v2, triAreaIDs[triIndex], heightfield, heightfield.bmin, heightfield.bmax, heightfield.cs, inverseCellSize, inverseCellHeight, flagMergeThreshold))
		{
			context->log(RC_LOG_ERROR, "rcRasterizeTriangles: Out of memory.");
			return false;
		}
	}

	return true;
}

bool rcRasterizeTriangles(rcContext* context,
                          const float* verts, const int /*nv*/,
                          const unsigned short* tris, const unsigned char* triAreaIDs, const int numTris,
                          rcHeightfield& heightfield, const int flagMergeThreshold)
{
	rcAssert(context != NULL);

	rcScopedTimer timer(context, RC_TIMER_RASTERIZE_TRIANGLES);

	// Rasterize the triangles.
	const float inverseCellSize = 1.0f / heightfield.cs;
	const float inverseCellHeight = 1.0f / heightfield.ch;
	for (int triIndex = 0; triIndex < numTris; ++triIndex)
	{
		const float* v0 = &verts[tris[triIndex * 3 + 0] * 3];
		const float* v1 = &verts[tris[triIndex * 3 + 1] * 3];
		const float* v2 = &verts[tris[triIndex * 3 + 2] * 3];
		if (!rasterizeTri(v0, v1, v2, triAreaIDs[triIndex], heightfield, heightfield.bmin, heightfield.bmax, heightfield.cs, inverseCellSize, inverseCellHeight, flagMergeThreshold))
		{
			context->log(RC_LOG_ERROR, "rcRasterizeTriangles: Out of memory.");
			return false;
		}
	}

	return true;
}

bool rcRasterizeTriangles(rcContext* context,
                          const float* verts, const unsigned char* triAreaIDs, const int numTris,
                          rcHeightfield& heightfield, const int flagMergeThreshold)
{
	rcAssert(context != NULL);

	rcScopedTimer timer(context, RC_TIMER_RASTERIZE_TRIANGLES);

	// Rasterize the triangles.
	const float inverseCellSize = 1.0f / heightfield.cs;
	const float inverseCellHeight = 1.0f / heightfield.ch;
	for (int triIndex = 0; triIndex < numTris; ++triIndex)
	{
		const float* v0 = &verts[(triIndex * 3 + 0) * 3];
		const float* v1 = &verts[(triIndex * 3 + 1) * 3];
		const float* v2 = &verts[(triIndex * 3 + 2) * 3];
		if (!rasterizeTri(v0, v1, v2, triAreaIDs[triIndex], heightfield, heightfield.bmin, heightfield.bmax, heightfield.cs, inverseCellSize, inverseCellHeight, flagMergeThreshold))
		{
			context->log(RC_LOG_ERROR, "rcRasterizeTriangles: Out of memory.");
			return false;
		}
	}

	return true;
}
