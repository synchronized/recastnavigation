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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include "DetourNavMesh.h"
#include "DetourCommon.h"
#include "DetourMath.h"
#include "DetourNavMeshBuilder.h"
#include "DetourAlloc.h"
#include "DetourAssert.h"

static unsigned short MESH_NULL_IDX = 0xffff;

//包围层次盒
struct BVItem
{
	unsigned short bmin[3]; //最小边界
	unsigned short bmax[3]; //最大边界
	int i; //
};

//比较两个BVItem的x分量
static int compareItemX(const void* va, const void* vb)
{
	const BVItem* a = (const BVItem*)va;
	const BVItem* b = (const BVItem*)vb;
	if (a->bmin[0] < b->bmin[0])
		return -1;
	if (a->bmin[0] > b->bmin[0])
		return 1;
	return 0;
}

//比较两个BVItem的y分量
static int compareItemY(const void* va, const void* vb)
{
	const BVItem* a = (const BVItem*)va;
	const BVItem* b = (const BVItem*)vb;
	if (a->bmin[1] < b->bmin[1])
		return -1;
	if (a->bmin[1] > b->bmin[1])
		return 1;
	return 0;
}

//比较两个BVItem的z分量
static int compareItemZ(const void* va, const void* vb)
{
	const BVItem* a = (const BVItem*)va;
	const BVItem* b = (const BVItem*)vb;
	if (a->bmin[2] < b->bmin[2])
		return -1;
	if (a->bmin[2] > b->bmin[2])
		return 1;
	return 0;
}

//求items的最小(bmin)/最大(bmax)边界
static void calcExtends(BVItem* items, const int /*nitems*/, const int imin, const int imax,
						unsigned short* bmin, unsigned short* bmax)
{
	bmin[0] = items[imin].bmin[0];
	bmin[1] = items[imin].bmin[1];
	bmin[2] = items[imin].bmin[2];

	bmax[0] = items[imin].bmax[0];
	bmax[1] = items[imin].bmax[1];
	bmax[2] = items[imin].bmax[2];

	for (int i = imin+1; i < imax; ++i)
	{
		const BVItem& it = items[i];
		if (it.bmin[0] < bmin[0]) bmin[0] = it.bmin[0];
		if (it.bmin[1] < bmin[1]) bmin[1] = it.bmin[1];
		if (it.bmin[2] < bmin[2]) bmin[2] = it.bmin[2];

		if (it.bmax[0] > bmax[0]) bmax[0] = it.bmax[0];
		if (it.bmax[1] > bmax[1]) bmax[1] = it.bmax[1];
		if (it.bmax[2] > bmax[2]) bmax[2] = it.bmax[2];
	}
}

//求数值最大的坐标轴(0:x, 1:y, 2:z)
inline int longestAxis(unsigned short x, unsigned short y, unsigned short z)
{
	int	axis = 0;
	unsigned short maxVal = x;
	if (y > maxVal)
	{
		axis = 1;
		maxVal = y;
	}
	if (z > maxVal)
	{
		axis = 2;
	}
	return axis;
}

//细分bvtree
static void subdivide(BVItem* items, int nitems, int imin, int imax, int& curNode, dtBVNode* nodes)
{
	int inum = imax - imin;
	int icur = curNode;

	dtBVNode& node = nodes[curNode++];

	if (inum == 1)
	{
		// Leaf 叶子
		node.bmin[0] = items[imin].bmin[0];
		node.bmin[1] = items[imin].bmin[1];
		node.bmin[2] = items[imin].bmin[2];

		node.bmax[0] = items[imin].bmax[0];
		node.bmax[1] = items[imin].bmax[1];
		node.bmax[2] = items[imin].bmax[2];

		node.i = items[imin].i;
	}
	else
	{
		// Split 切割
		calcExtends(items, nitems, imin, imax, node.bmin, node.bmax);

		int	axis = longestAxis(node.bmax[0] - node.bmin[0],
							   node.bmax[1] - node.bmin[1],
							   node.bmax[2] - node.bmin[2]);

		if (axis == 0)
		{
			// Sort along x-axis
			qsort(items+imin, inum, sizeof(BVItem), compareItemX);
		}
		else if (axis == 1)
		{
			// Sort along y-axis
			qsort(items+imin, inum, sizeof(BVItem), compareItemY);
		}
		else
		{
			// Sort along z-axis
			qsort(items+imin, inum, sizeof(BVItem), compareItemZ);
		}

		int isplit = imin+inum/2;

		// Left
		subdivide(items, nitems, imin, isplit, curNode, nodes);
		// Right
		subdivide(items, nitems, isplit, imax, curNode, nodes);

		int iescape = curNode - icur;
		// Negative index means escape.
		node.i = -iescape;
	}
}

//创建BVTree
static int createBVTree(dtNavMeshCreateParams* params, dtBVNode* nodes, int /*nnodes*/)
{
	// Build tree
	float quantFactor = 1 / params->cs;
	BVItem* items = (BVItem*)dtAlloc(sizeof(BVItem)*params->polyCount, DT_ALLOC_TEMP);
	for (int i = 0; i < params->polyCount; i++)
	{
		BVItem& it = items[i];
		it.i = i;
		// Calc polygon bounds. Use detail meshes if available.
		// 计算多边形边界。 如果可用，请使用细节网格。
		if (params->detailMeshes)
		{
			int vb = (int)params->detailMeshes[i*4+0]; //顶点下标开始
			int ndv = (int)params->detailMeshes[i*4+1]; //顶点数
			float bmin[3];
			float bmax[3];

			const float* dv = &params->detailVerts[vb*3];
			dtVcopy(bmin, dv);
			dtVcopy(bmax, dv);

			for (int j = 1; j < ndv; j++)
			{
				dtVmin(bmin, &dv[j * 3]);
				dtVmax(bmax, &dv[j * 3]);
			}

			// BV-tree uses cs for all dimensions
			// BV-tree 对所有维度都使用 cs(cell size)
			it.bmin[0] = (unsigned short)dtClamp((int)((bmin[0] - params->bmin[0])*quantFactor), 0, 0xffff);
			it.bmin[1] = (unsigned short)dtClamp((int)((bmin[1] - params->bmin[1])*quantFactor), 0, 0xffff);
			it.bmin[2] = (unsigned short)dtClamp((int)((bmin[2] - params->bmin[2])*quantFactor), 0, 0xffff);

			it.bmax[0] = (unsigned short)dtClamp((int)((bmax[0] - params->bmin[0])*quantFactor), 0, 0xffff);
			it.bmax[1] = (unsigned short)dtClamp((int)((bmax[1] - params->bmin[1])*quantFactor), 0, 0xffff);
			it.bmax[2] = (unsigned short)dtClamp((int)((bmax[2] - params->bmin[2])*quantFactor), 0, 0xffff);
		}
		else
		{
			//多边形数据(多边形顶点索引以MESH_NULL_IDX结尾)
			const unsigned short* p = &params->polys[i*params->nvp * 2];
			it.bmin[0] = it.bmax[0] = params->verts[p[0] * 3 + 0];
			it.bmin[1] = it.bmax[1] = params->verts[p[0] * 3 + 1];
			it.bmin[2] = it.bmax[2] = params->verts[p[0] * 3 + 2];

			//计算最小/最大边界
			for (int j = 1; j < params->nvp; ++j)
			{
				if (p[j] == MESH_NULL_IDX) break;
				unsigned short x = params->verts[p[j] * 3 + 0];
				unsigned short y = params->verts[p[j] * 3 + 1];
				unsigned short z = params->verts[p[j] * 3 + 2];

				if (x < it.bmin[0]) it.bmin[0] = x;
				if (y < it.bmin[1]) it.bmin[1] = y;
				if (z < it.bmin[2]) it.bmin[2] = z;

				if (x > it.bmax[0]) it.bmax[0] = x;
				if (y > it.bmax[1]) it.bmax[1] = y;
				if (z > it.bmax[2]) it.bmax[2] = z;
			}
			// Remap y, 重新计算y
			it.bmin[1] = (unsigned short)dtMathFloorf((float)it.bmin[1] * params->ch / params->cs);
			it.bmax[1] = (unsigned short)dtMathCeilf((float)it.bmax[1] * params->ch / params->cs);
		}
	}

	int curNode = 0;
	subdivide(items, params->polyCount, 0, params->polyCount, curNode, nodes);

	dtFree(items);

	return curNode;
}

//在xz平面上的投影, 计算点(pt)位于AABB盒(bmin, bmax)方位(0-7,0xff)
static unsigned char classifyOffMeshPoint(const float* pt, const float* bmin, const float* bmax)
{
	static const unsigned char XP = 1<<0;
	static const unsigned char ZP = 1<<1;
	static const unsigned char XM = 1<<2;
	static const unsigned char ZM = 1<<3;

	unsigned char outcode = 0;
	outcode |= (pt[0] >= bmax[0]) ? XP : 0;
	outcode |= (pt[2] >= bmax[2]) ? ZP : 0;
	outcode |= (pt[0] < bmin[0])  ? XM : 0;
	outcode |= (pt[2] < bmin[2])  ? ZM : 0;

	switch (outcode)
	{
	case XP: return 0;
	case XP|ZP: return 1;
	case ZP: return 2;
	case XM|ZP: return 3;
	case XM: return 4;
	case XM|ZM: return 5;
	case ZM: return 6;
	case XP|ZM: return 7;
	};

	return 0xff;
}

// TODO: Better error handling.

/// @par
/// 
/// The output data array is allocated using the detour allocator (dtAlloc()).  The method
/// used to free the memory will be determined by how the tile is added to the navigation
/// mesh.
///
/// @see dtNavMesh, dtNavMesh::addTile()
/// 输出数据数组使用detour allocator(dtAlloc())进行分配。用于释放内存的方法将取决于如何将tile添加到导航网格中。
bool dtCreateNavMeshData(dtNavMeshCreateParams* params, unsigned char** outData, int* outDataSize)
{
	if (params->nvp > DT_VERTS_PER_POLYGON)
		return false; //每个多边形顶点数大于6
	if (params->vertCount >= 0xffff)
		return false; //顶点数超过2字节
	if (!params->vertCount || !params->verts)
		return false; //没有顶点数据
	if (!params->polyCount || !params->polys)
		return false; //没有多边形数据

	const int nvp = params->nvp;

	// Classify off-mesh connection points. We store only the connections
	// whose start point is inside the tile.
	// 对off-mesh connections点进行分类。 我们只存储起点位于tile内部的连接。
	unsigned char* offMeshConClass = 0;
	int storedOffMeshConCount = 0; //保存数量
	int offMeshConLinkCount = 0;

	if (params->offMeshConCount > 0)
	{
		offMeshConClass = (unsigned char*)dtAlloc(sizeof(unsigned char)*params->offMeshConCount*2, DT_ALLOC_TEMP);
		if (!offMeshConClass)
			return false;

		// Find tight heigh bounds, used for culling out off-mesh start locations.
		// 找到严格的高度边界，用于剔除网络外连接的起始位置。
		float hmin = FLT_MAX;
		float hmax = -FLT_MAX;

		if (params->detailVerts && params->detailVertsCount)
		{
			// 如果有详细网格（三角形）的数据就用详细网格的顶点
			for (int i = 0; i < params->detailVertsCount; ++i)
			{
				const float h = params->detailVerts[i*3+1];
				hmin = dtMin(hmin,h);
				hmax = dtMax(hmax,h);
			}
		}
		else
		{
			// 没有生成详细网格，就只能用 poly 网格的顶点
			for (int i = 0; i < params->vertCount; ++i)
			{
				const unsigned short* iv = &params->verts[i*3];
				const float h = params->bmin[1] + iv[1] * params->ch;
				hmin = dtMin(hmin,h);
				hmax = dtMax(hmax,h);
			}
		}
		hmin -= params->walkableClimb;
		hmax += params->walkableClimb;
		float bmin[3], bmax[3];
		dtVcopy(bmin, params->bmin);
		dtVcopy(bmax, params->bmax);
		bmin[1] = hmin;
		bmax[1] = hmax;

		for (int i = 0; i < params->offMeshConCount; ++i)
		{
			const float* p0 = &params->offMeshConVerts[(i*2+0)*3];
			const float* p1 = &params->offMeshConVerts[(i*2+1)*3];
			// 计算点相对于该 Tile 的方位
			offMeshConClass[i*2+0] = classifyOffMeshPoint(p0, bmin, bmax);
			offMeshConClass[i*2+1] = classifyOffMeshPoint(p1, bmin, bmax);

			// Zero out off-mesh start positions which are not even potentially touching the mesh.
			// 将甚至不会接触网格的off-mesh起始位置归零。
			if (offMeshConClass[i*2+0] == 0xff)
			{
				// 如果起始点位于该 Tile 的 AABB 上方或下方，则重新修改标记为 0
				// 有个疑问，如果这种情况被标记为 0，那和 XP 情况的标记为 0 不是重复了吗
				if (p0[1] < bmin[1] || p0[1] > bmax[1])
					offMeshConClass[i*2+0] = 0;
			}

			// Cound how many links should be allocated for off-mesh connections.
			// 计算起点或终点可能需要连接到该 Tile 的数量
			if (offMeshConClass[i*2+0] == 0xff)
				offMeshConLinkCount++;
			if (offMeshConClass[i*2+1] == 0xff)
				offMeshConLinkCount++;

			// 计算起点可能需要连接到该 Tile 的数量
			if (offMeshConClass[i*2+0] == 0xff)
				storedOffMeshConCount++;
		}
	}

	// Off-mesh connections are stored as polygons, adjust values.
	// Off-mesh connection结构被看成是一个连接节点单元，这里只需要计算起点位于 Tile 内的即可(存储为多边形)
	const int totPolyCount = params->polyCount + storedOffMeshConCount;
	const int totVertCount = params->vertCount + storedOffMeshConCount*2;

	// Find portal edges which are at tile borders.
	// 查找位于tile边界处的portal edges。
	int edgeCount = 0;
	int portalCount = 0;
	for (int i = 0; i < params->polyCount; ++i)
	{
		const unsigned short* p = &params->polys[i*2*nvp];
		for (int j = 0; j < nvp; ++j)
		{
			if (p[j] == MESH_NULL_IDX) { break; }
			//正常边的数量
			edgeCount++;

			// 在前面构造 PolyMesh 并计算 MeshAdjancy 的时候对于 portal edge 其首位(short)设置为 1
			if (p[nvp+j] & 0x8000)
			{
				unsigned short dir = p[nvp+j] & 0xf;
				if (dir != 0xf) {
					// portal 边的数量
					portalCount++;
				}
			}
		}
	}

	const int maxLinkCount = edgeCount + portalCount*2 + offMeshConLinkCount*2;

	// Find unique detail vertices.
	int uniqueDetailVertCount = 0;
	// 三角形数量 如果有详细网格，则用详细网格的三角形数 如果没有，则从 nav poly 中切出三角形
	int detailTriCount = 0;
	if (params->detailMeshes)
	{
		// Has detail mesh, count unique detail vertex count and use input detail tri count.
		detailTriCount = params->detailTriCount;
		for (int i = 0; i < params->polyCount; ++i)
		{
			const unsigned short* p = &params->polys[i*nvp*2];
			int ndv = params->detailMeshes[i*4+1];
			int nv = 0;
			for (int j = 0; j < nvp; ++j)
			{
				if (p[j] == MESH_NULL_IDX) break;
				nv++;
			}
			ndv -= nv;
			uniqueDetailVertCount += ndv;
		}
	}
	else
	{
		// No input detail mesh, build detail mesh from nav polys.
		// 没有输入细节网格，从导航多边形构建细节网格。
		uniqueDetailVertCount = 0; // No extra detail verts.
		detailTriCount = 0; //三角形个数
		for (int i = 0; i < params->polyCount; ++i)
		{
			const unsigned short* p = &params->polys[i*nvp*2];
			int nv = 0;
			for (int j = 0; j < nvp; ++j)
			{
				if (p[j] == MESH_NULL_IDX) break;
				nv++;
			}
			detailTriCount += nv-2;
		}
	}

	// Calculate data size
	const int headerSize = dtAlign4(sizeof(dtMeshHeader));
	const int vertsSize = dtAlign4(sizeof(float)*3*totVertCount);
	const int polysSize = dtAlign4(sizeof(dtPoly)*totPolyCount);
	const int linksSize = dtAlign4(sizeof(dtLink)*maxLinkCount);
	const int detailMeshesSize = dtAlign4(sizeof(dtPolyDetail)*params->polyCount);
	const int detailVertsSize = dtAlign4(sizeof(float)*3*uniqueDetailVertCount);
	const int detailTrisSize = dtAlign4(sizeof(unsigned char)*4*detailTriCount);
	const int bvTreeSize = params->buildBvTree ? dtAlign4(sizeof(dtBVNode)*params->polyCount*2) : 0;
	const int offMeshConsSize = dtAlign4(sizeof(dtOffMeshConnection)*storedOffMeshConCount);

	const int dataSize = headerSize + vertsSize + polysSize + linksSize +
						 detailMeshesSize + detailVertsSize + detailTrisSize +
						 bvTreeSize + offMeshConsSize;

	unsigned char* data = (unsigned char*)dtAlloc(sizeof(unsigned char)*dataSize, DT_ALLOC_PERM);
	if (!data)
	{
		dtFree(offMeshConClass);
		return false;
	}
	memset(data, 0, dataSize);

	unsigned char* d = data;

	dtMeshHeader* header = dtGetThenAdvanceBufferPointer<dtMeshHeader>(d, headerSize);
	float* navVerts = dtGetThenAdvanceBufferPointer<float>(d, vertsSize);
	dtPoly* navPolys = dtGetThenAdvanceBufferPointer<dtPoly>(d, polysSize);
	d += linksSize; // Ignore links; just leave enough space for them. They'll be created on load.
	dtPolyDetail* navDMeshes = dtGetThenAdvanceBufferPointer<dtPolyDetail>(d, detailMeshesSize);
	float* navDVerts = dtGetThenAdvanceBufferPointer<float>(d, detailVertsSize);
	unsigned char* navDTris = dtGetThenAdvanceBufferPointer<unsigned char>(d, detailTrisSize);
	dtBVNode* navBvtree = dtGetThenAdvanceBufferPointer<dtBVNode>(d, bvTreeSize);
	dtOffMeshConnection* offMeshCons = dtGetThenAdvanceBufferPointer<dtOffMeshConnection>(d, offMeshConsSize);


	// Store header
	header->magic = DT_NAVMESH_MAGIC;
	header->version = DT_NAVMESH_VERSION;
	header->x = params->tileX;
	header->y = params->tileY;
	header->layer = params->tileLayer;
	header->userId = params->userId;
	header->polyCount = totPolyCount;
	header->vertCount = totVertCount;
	header->maxLinkCount = maxLinkCount;
	dtVcopy(header->bmin, params->bmin);
	dtVcopy(header->bmax, params->bmax);
	header->detailMeshCount = params->polyCount;
	header->detailVertCount = uniqueDetailVertCount;
	header->detailTriCount = detailTriCount;
	header->bvQuantFactor = 1.0f / params->cs;
	header->offMeshBase = params->polyCount;
	header->walkableHeight = params->walkableHeight;
	header->walkableRadius = params->walkableRadius;
	header->walkableClimb = params->walkableClimb;
	header->offMeshConCount = storedOffMeshConCount;
	header->bvNodeCount = params->buildBvTree ? params->polyCount*2 : 0;

	const int offMeshVertsBase = params->vertCount;
	const int offMeshPolyBase = params->polyCount;

	// Store vertices
	// Mesh vertices
	for (int i = 0; i < params->vertCount; ++i)
	{
		const unsigned short* iv = &params->verts[i*3];
		float* v = &navVerts[i*3];
		v[0] = params->bmin[0] + iv[0] * params->cs;
		v[1] = params->bmin[1] + iv[1] * params->ch;
		v[2] = params->bmin[2] + iv[2] * params->cs;
	}
	// Off-mesh link vertices.
	// Off-mesh connection顶点
	int n = 0;
	for (int i = 0; i < params->offMeshConCount; ++i)
	{
		// Only store connections which start from this tile.
		// 仅存储从此图块开始的连接。
		if (offMeshConClass[i*2+0] == 0xff)
		{
			const float* linkv = &params->offMeshConVerts[i*2*3];
			float* v = &navVerts[(offMeshVertsBase + n*2)*3];
			dtVcopy(&v[0], &linkv[0]);
			dtVcopy(&v[3], &linkv[3]);
			n++;
		}
	}

	// Store polygons
	// Mesh polys
	const unsigned short* src = params->polys;
	for (int i = 0; i < params->polyCount; ++i)
	{
		// 这里填充后面要使用的 Poly 对象
		dtPoly* p = &navPolys[i];
		p->vertCount = 0;
		p->flags = params->polyFlags[i];
		p->setArea(params->polyAreas[i]);
		p->setType(DT_POLYTYPE_GROUND);
		for (int j = 0; j < nvp; ++j)
		{
			if (src[j] == MESH_NULL_IDX) break;
			p->verts[j] = src[j];
			if (src[nvp+j] & 0x8000)
			{
				// Border or portal edge.
				unsigned short dir = src[nvp+j] & 0xf;
				if (dir == 0xf) { // Border
					// 纯边界边了
					p->neis[j] = 0;
				}
				else if (dir == 0) // Portal x-
					p->neis[j] = DT_EXT_LINK | 4;
				else if (dir == 1) // Portal z+
					p->neis[j] = DT_EXT_LINK | 2;
				else if (dir == 2) // Portal x+
					p->neis[j] = DT_EXT_LINK | 0;
				else if (dir == 3) // Portal z-
					p->neis[j] = DT_EXT_LINK | 6;
			}
			else
			{
				// Normal connection
				// Tile 内部的多边形之间的邻接边
				// 之前在 poly 里是记录的 poly 序号
				// polys[p0 + maxVertsPerPoly + e.polyEdge[0]] = e.poly[1];
				// polys[p1 + maxVertsPerPoly + e.polyEdge[1]] = e.poly[0];
				// 注意这里序号 +1 了，后面使用的话应该要再次减 1 吧
				p->neis[j] = src[nvp+j]+1;
			}

			p->vertCount++;
		}
		src += nvp*2;
	}
	// Off-mesh connection vertices.
	n = 0;
	for (int i = 0; i < params->offMeshConCount; ++i)
	{
		// Only store connections which start from this tile.
		// 仅存储从此图块开始的连接。
		if (offMeshConClass[i*2+0] == 0xff)
		{
			dtPoly* p = &navPolys[offMeshPolyBase+n];
			p->vertCount = 2;
			p->verts[0] = (unsigned short)(offMeshVertsBase + n*2+0);
			p->verts[1] = (unsigned short)(offMeshVertsBase + n*2+1);
			p->flags = params->offMeshConFlags[i];
			p->setArea(params->offMeshConAreas[i]);
			p->setType(DT_POLYTYPE_OFFMESH_CONNECTION);
			n++;
		}
	}

	// Store detail meshes and vertices.
	// The nav polygon vertices are stored as the first vertices on each mesh.
	// We compress the mesh data by skipping them and using the navmesh coordinates.
	//
	// 存储细节网格和顶点。
	// 导航多边形顶点存储为每个网格上的第一个顶点。
	// 我们通过跳过网格数据并使用导航网格坐标来压缩网格数据。
	if (params->detailMeshes)
	{
		unsigned short vbase = 0;
		for (int i = 0; i < params->polyCount; ++i)
		{
			dtPolyDetail& dtl = navDMeshes[i];
			const int vb = (int)params->detailMeshes[i*4+0];
			const int ndv = (int)params->detailMeshes[i*4+1];
			const int nv = navPolys[i].vertCount;
			dtl.vertBase = (unsigned int)vbase;
			dtl.vertCount = (unsigned char)(ndv-nv);
			dtl.triBase = (unsigned int)params->detailMeshes[i*4+2];
			dtl.triCount = (unsigned char)params->detailMeshes[i*4+3];
			// Copy vertices except the first 'nv' verts which are equal to nav poly verts.
			// 复制除前面“nv”顶点之外的顶点，这些顶点等于导航多边形顶点。
			if (ndv-nv)
			{
				memcpy(&navDVerts[vbase*3], &params->detailVerts[(vb+nv)*3], sizeof(float)*3*(ndv-nv));
				vbase += (unsigned short)(ndv-nv);
			}
		}
		// Store triangles.
		memcpy(navDTris, params->detailTris, sizeof(unsigned char)*4*params->detailTriCount);
	}
	else
	{
		// Create dummy detail mesh by triangulating polys.
		// 通过对多边形进行三角剖分来创建虚拟细节网格。
		int tbase = 0;
		for (int i = 0; i < params->polyCount; ++i)
		{
			dtPolyDetail& dtl = navDMeshes[i];
			const int nv = navPolys[i].vertCount;
			dtl.vertBase = 0;
			dtl.vertCount = 0;
			dtl.triBase = (unsigned int)tbase;
			dtl.triCount = (unsigned char)(nv-2);
			// Triangulate polygon (local indices).
			// 对多边形进行三角剖分（局部索引）。
			for (int j = 2; j < nv; ++j)
			{
				unsigned char* t = &navDTris[tbase*4];
				t[0] = 0;
				t[1] = (unsigned char)(j-1);
				t[2] = (unsigned char)j;
				// Bit for each edge that belongs to poly boundary.
				// 属于多边形边界的每条边的位。
				t[3] = (1<<2);
				if (j == 2) t[3] |= (1<<0); //第一个三角形
				if (j == nv-1) t[3] |= (1<<4); //最后一个三角形
				tbase++;
			}
		}
	}

	// Store and create BVtree.
	// 存储并创建 BVtree。
	if (params->buildBvTree)
	{
		createBVTree(params, navBvtree, 2*params->polyCount);
	}

	// Store Off-Mesh connections.
	// 存储离网连接。
	n = 0;
	for (int i = 0; i < params->offMeshConCount; ++i)
	{
		// Only store connections which start from this tile.
		// 仅存储从此图块开始的连接。
		if (offMeshConClass[i*2+0] == 0xff)
		{
			dtOffMeshConnection* con = &offMeshCons[n];
			con->poly = (unsigned short)(offMeshPolyBase + n);
			// Copy connection end-points.
			const float* endPts = &params->offMeshConVerts[i*2*3];
			dtVcopy(&con->pos[0], &endPts[0]);
			dtVcopy(&con->pos[3], &endPts[3]);
			con->rad = params->offMeshConRad[i];
			con->flags = params->offMeshConDir[i] ? DT_OFFMESH_CON_BIDIR : 0;
			// 因为起点已在 tile 内，所以这里只要记录终点的相对位置即可
			con->side = offMeshConClass[i*2+1];
			if (params->offMeshConUserID)
				con->userId = params->offMeshConUserID[i];
			n++;
		}
	}

	dtFree(offMeshConClass);

	*outData = data;
	*outDataSize = dataSize;

	return true;
}

bool dtNavMeshHeaderSwapEndian(unsigned char* data, const int /*dataSize*/)
{
	dtMeshHeader* header = (dtMeshHeader*)data;

	int swappedMagic = DT_NAVMESH_MAGIC;
	int swappedVersion = DT_NAVMESH_VERSION;
	dtSwapEndian(&swappedMagic);
	dtSwapEndian(&swappedVersion);

	if ((header->magic != DT_NAVMESH_MAGIC || header->version != DT_NAVMESH_VERSION) &&
		(header->magic != swappedMagic || header->version != swappedVersion))
	{
		return false;
	}

	dtSwapEndian(&header->magic);
	dtSwapEndian(&header->version);
	dtSwapEndian(&header->x);
	dtSwapEndian(&header->y);
	dtSwapEndian(&header->layer);
	dtSwapEndian(&header->userId);
	dtSwapEndian(&header->polyCount);
	dtSwapEndian(&header->vertCount);
	dtSwapEndian(&header->maxLinkCount);
	dtSwapEndian(&header->detailMeshCount);
	dtSwapEndian(&header->detailVertCount);
	dtSwapEndian(&header->detailTriCount);
	dtSwapEndian(&header->bvNodeCount);
	dtSwapEndian(&header->offMeshConCount);
	dtSwapEndian(&header->offMeshBase);
	dtSwapEndian(&header->walkableHeight);
	dtSwapEndian(&header->walkableRadius);
	dtSwapEndian(&header->walkableClimb);
	dtSwapEndian(&header->bmin[0]);
	dtSwapEndian(&header->bmin[1]);
	dtSwapEndian(&header->bmin[2]);
	dtSwapEndian(&header->bmax[0]);
	dtSwapEndian(&header->bmax[1]);
	dtSwapEndian(&header->bmax[2]);
	dtSwapEndian(&header->bvQuantFactor);

	// Freelist index and pointers are updated when tile is added, no need to swap.
	// 添加tile时，空闲列表索引和指针会更新，无需交换。

	return true;
}

/// @par
///
/// @warning This function assumes that the header is in the correct endianess already.
/// Call #dtNavMeshHeaderSwapEndian() first on the data if the data is expected to be in wrong endianess
/// to start with. Call #dtNavMeshHeaderSwapEndian() after the data has been swapped if converting from
/// native to foreign endianess.
/// @warning 该函数假设标头已经采用正确的字节顺序。
/// 如果数据预计以错误的字节顺序开始，请首先对数据调用 #dtNavMeshHeaderSwapEndian()。
/// 如果从本机字节序转换为外部字节序，则在交换数据后调用 #dtNavMeshHeaderSwapEndian()。
bool dtNavMeshDataSwapEndian(unsigned char* data, const int /*dataSize*/)
{
	// Make sure the data is in right format.
	// 确保数据格式正确。
	dtMeshHeader* header = (dtMeshHeader*)data;
	if (header->magic != DT_NAVMESH_MAGIC)
		return false;
	if (header->version != DT_NAVMESH_VERSION)
		return false;

	// Patch header pointers.
	const int headerSize = dtAlign4(sizeof(dtMeshHeader));
	const int vertsSize = dtAlign4(sizeof(float)*3*header->vertCount);
	const int polysSize = dtAlign4(sizeof(dtPoly)*header->polyCount);
	const int linksSize = dtAlign4(sizeof(dtLink)*(header->maxLinkCount));
	const int detailMeshesSize = dtAlign4(sizeof(dtPolyDetail)*header->detailMeshCount);
	const int detailVertsSize = dtAlign4(sizeof(float)*3*header->detailVertCount);
	const int detailTrisSize = dtAlign4(sizeof(unsigned char)*4*header->detailTriCount);
	const int bvtreeSize = dtAlign4(sizeof(dtBVNode)*header->bvNodeCount);
	const int offMeshLinksSize = dtAlign4(sizeof(dtOffMeshConnection)*header->offMeshConCount);

	unsigned char* d = data + headerSize;
	float* verts = dtGetThenAdvanceBufferPointer<float>(d, vertsSize);
	dtPoly* polys = dtGetThenAdvanceBufferPointer<dtPoly>(d, polysSize);
	d += linksSize; // Ignore links; they technically should be endian-swapped but all their data is overwritten on load anyway.
	// 忽略链接； 从技术上讲，它们应该是字节序交换的，但无论如何，它们的所有数据都会在加载时被覆盖。
	//dtLink* links = dtGetThenAdvanceBufferPointer<dtLink>(d, linksSize);
	dtPolyDetail* detailMeshes = dtGetThenAdvanceBufferPointer<dtPolyDetail>(d, detailMeshesSize);
	float* detailVerts = dtGetThenAdvanceBufferPointer<float>(d, detailVertsSize);
	d += detailTrisSize; // Ignore detail tris; single bytes can't be endian-swapped.
	// 忽略细节三角形； 单个字节不能进行字节序交换。
	//unsigned char* detailTris = dtGetThenAdvanceBufferPointer<unsigned char>(d, detailTrisSize);
	dtBVNode* bvTree = dtGetThenAdvanceBufferPointer<dtBVNode>(d, bvtreeSize);
	dtOffMeshConnection* offMeshCons = dtGetThenAdvanceBufferPointer<dtOffMeshConnection>(d, offMeshLinksSize);

	// Vertices
	for (int i = 0; i < header->vertCount*3; ++i)
	{
		dtSwapEndian(&verts[i]);
	}

	// Polys
	for (int i = 0; i < header->polyCount; ++i)
	{
		dtPoly* p = &polys[i];
		// poly->firstLink is update when tile is added, no need to swap.
		// poly->firstLink 在添加图块时更新，无需交换。
		for (int j = 0; j < DT_VERTS_PER_POLYGON; ++j)
		{
			dtSwapEndian(&p->verts[j]);
			dtSwapEndian(&p->neis[j]);
		}
		dtSwapEndian(&p->flags);
	}

	// Links are rebuild when tile is added, no need to swap.
	// 添加图块时会重建链接，无需交换。

	// Detail meshes
	for (int i = 0; i < header->detailMeshCount; ++i)
	{
		dtPolyDetail* pd = &detailMeshes[i];
		dtSwapEndian(&pd->vertBase);
		dtSwapEndian(&pd->triBase);
	}

	// Detail verts
	for (int i = 0; i < header->detailVertCount*3; ++i)
	{
		dtSwapEndian(&detailVerts[i]);
	}

	// BV-tree
	for (int i = 0; i < header->bvNodeCount; ++i)
	{
		dtBVNode* node = &bvTree[i];
		for (int j = 0; j < 3; ++j)
		{
			dtSwapEndian(&node->bmin[j]);
			dtSwapEndian(&node->bmax[j]);
		}
		dtSwapEndian(&node->i);
	}

	// Off-mesh Connections.
	for (int i = 0; i < header->offMeshConCount; ++i)
	{
		dtOffMeshConnection* con = &offMeshCons[i];
		for (int j = 0; j < 6; ++j)
			dtSwapEndian(&con->pos[j]);
		dtSwapEndian(&con->rad);
		dtSwapEndian(&con->poly);
	}

	return true;
}
