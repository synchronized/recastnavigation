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

#ifndef DETOURNAVMESHBUILDER_H
#define DETOURNAVMESHBUILDER_H

#include "DetourAlloc.h"

/// Represents the source data used to build an navigation mesh tile.
/// @ingroup detour
/// 表示用于构建导航网格图块的源数据。
struct dtNavMeshCreateParams
{

	/// @name Polygon Mesh Attributes
	/// Used to create the base navigation graph.
	/// See #rcPolyMesh for details related to these attributes.
	/// @name 多边形网格属性
	/// 用于创建基本导航图。
	/// 有关这些属性的详细信息，请参阅#rcPolyMesh。
	/// @{

	///< The polygon mesh vertices. [(x, y, z) * #vertCount] [Unit: vx]
	///< 多边形网格顶点。 [(x, y, z) * #vertCount] [单位：vx(体素)]
	const unsigned short* verts;
	///< The number vertices in the polygon mesh. [Limit: >= 3]
	///< 多边形网格中的顶点数。 [限制：>= 3]
	int vertCount;
	///< The polygon data. [Size: #polyCount * 2 * #nvp]
	///< 多边形数据(多边形顶点索引以MESH_NULL_IDX结尾)。 [大小：#polyCount * 2 * #nvp(多边形最多顶点数)]
	const unsigned short* polys;
	///< The user defined flags assigned to each polygon. [Size: #polyCount]
	///< 分配给每个多边形的用户定义标志。 [尺寸：#polyCount]
	const unsigned short* polyFlags;
	///< The user defined area ids assigned to each polygon. [Size: #polyCount]
	///< 分配给每个多边形的用户定义区域 ID。 [尺寸：#polyCount]
	const unsigned char* polyAreas;
	///< Number of polygons in the mesh. [Limit: >= 1]
	///< 网格中的多边形数量。 [限制：>= 1]
	int polyCount;
	///< Number maximum number of vertices per polygon. [Limit: >= 3]
	///< Number 每个多边形的最大顶点数。 [限制：>= 3]
	int nvp;

	/// @}
	/// @name Height Detail Attributes (Optional)
	/// See #rcPolyMeshDetail for details related to these attributes.
	/// @name 高度详细属性（可选）
	/// 有关这些属性的详细信息，请参阅#rcPolyMeshDetail。
	/// @{

	///< The height detail sub-mesh data. [Size: 4 * #polyCount]
	///< 高度细节子网格数据。 [尺寸：4 * #polyCount][顶点下标开始, 顶点数, ...]
	const unsigned int* detailMeshes;
	///< The detail mesh vertices. [Size: 3 * #detailVertsCount] [Unit: wu]
	///< 详细网格顶点。 [大小：3 * #detailVertsCount] [单位：wu(世界单位)]
	const float* detailVerts;
	///< The number of vertices in the detail mesh.
	///< 细节网格中的顶点数。
	int detailVertsCount;
	///< The detail mesh triangles. [Size: 4 * #detailTriCount]
	///< 细节网格三角形。 [尺寸：4 * #detailTriCount]
	const unsigned char* detailTris;
	///< The number of triangles in the detail mesh.
	///< 细节网格中三角形的数量。
	int detailTriCount;

	/// @}
	/// @name Off-Mesh Connections Attributes (Optional)
	/// Used to define a custom point-to-point edge within the navigation graph, an
	/// off-mesh connection is a user defined traversable connection made up to two vertices,
	/// at least one of which resides within a navigation mesh polygon.
	/// @name 离网连接属性（可选）
	/// 用于在导航图中定义自定义点对点边，网络外连接是用户定义的由两个顶点组成的可遍历连接，
	/// 其中至少一个顶点位于导航网格多边形内。，
	/// @{

	/// Off-mesh connection vertices. [(ax, ay, az, bx, by, bz) * #offMeshConCount] [Unit: wu]
	/// 离网格连接顶点。 [(ax, ay, az, bx, by, bz) * #offMeshConCount] [单位：wu]
	const float* offMeshConVerts;
	/// Off-mesh connection radii. [Size: #offMeshConCount] [Unit: wu]
	/// 网络外连接半径。 [大小：#offMeshConCount] [单位：wu(世界单位)]
	const float* offMeshConRad;
	/// User defined flags assigned to the off-mesh connections. [Size: #offMeshConCount]
	/// 分配给网络外连接的用户定义标志。 [大小：#offMeshConCount]
	const unsigned short* offMeshConFlags;
	/// User defined area ids assigned to the off-mesh connections. [Size: #offMeshConCount]
	/// 分配给网络外连接的用户定义区域 ID。 [大小：#offMeshConCount]
	const unsigned char* offMeshConAreas;
	/// The permitted travel direction of the off-mesh connections. [Size: #offMeshConCount]
	/// 网络外连接允许的行进方向。 [大小：#offMeshConCount]
	///
	/// 0 = Travel only from endpoint A to endpoint B.<br/>
	/// #DT_OFFMESH_CON_BIDIR = Bidirectional travel.
	/// 0 = 仅从端点 A 行驶至端点 B。<br/>
	/// #DT_OFFMESH_CON_BIDIR = 双向行驶
	const unsigned char* offMeshConDir;	
	/// The user defined ids of the off-mesh connection. [Size: #offMeshConCount]
	/// 用户定义的网络外连接 ID。 [大小：#offMeshConCount]
	const unsigned int* offMeshConUserID;
	/// The number of off-mesh connections. [Limit: >= 0]
	/// 网络外连接的数量。 [限制：>= 0]
	int offMeshConCount;

	/// @}
	/// @name Tile Attributes
	/// @note The tile grid/layer data can be left at zero if the destination is a single tile mesh.
	/// @name Tile 属性
	/// @note 如果目标是单个tile网格，则tile网格/层数据可以保留为零。
	/// @{

	///< The user defined id of the tile.
	///< 用户定义的图块 ID。
	unsigned int userId;
	///< The tile's x-grid location within the multi-tile destination mesh. (Along the x-axis.)
	///< 多tile目标网格内tile的 x 网格位置。 （沿 x 轴。）
	int tileX;
	///< The tile's y-grid location within the multi-tile destination mesh. (Along the z-axis.)
	///< 多tile目标网格内tile的 y 网格位置。 （沿 z 轴。）
	int tileY;
	///< The tile's layer within the layered destination mesh. [Limit: >= 0] (Along the y-axis.)
	///< 分层目标网格内tile的层。 [限制：>= 0]（沿 y 轴。）
	int tileLayer;
	///< The minimum bounds of the tile. [(x, y, z)] [Unit: wu]
	///< tile的最小边界。 [(x,y,z)][单位：wu]
	float bmin[3];
	///< The maximum bounds of the tile. [(x, y, z)] [Unit: wu]
	///< tile的最大边界。 [(x,y,z)][单位：wu]
	float bmax[3];

	/// @}
	/// @name General Configuration Attributes
	/// @name 常规配置属性
	/// @{

	///< The agent height. [Unit: wu]
	///< agent高度。 [单位：wu]
	float walkableHeight;
	///< The agent radius. [Unit: wu]
	///< agent半径。 [单位：wu]
	float walkableRadius;
	///< The agent maximum traversable ledge. (Up/Down) [Unit: wu]
	///< agent最大可穿越高度。 (上/下) [单位: wu]
	float walkableClimb;
	///< The xz-plane cell size of the polygon mesh. [Limit: > 0] [Unit: wu]
	///< 多边形网格的 xz 平面单元尺寸。 [限制：> 0] [单位：wu]
	float cs;
	///< The y-axis cell height of the polygon mesh. [Limit: > 0] [Unit: wu]
	///< 多边形网格的 y 轴单元高度。 [限制：> 0] [单位：wu]
	float ch;

	/// True if a bounding volume tree should be built for the tile.
	/// @note The BVTree is not normally needed for layered navigation meshes.
	/// 如果应为tile构建包围体树，则为 true。
	/// @note 分层导航网格通常不需要 BVTree。
	bool buildBvTree;

	/// @}
};

/// Builds navigation mesh tile data from the provided tile creation data.
/// @ingroup detour
///  @param[in]		params		Tile creation data.
///  @param[out]	outData		The resulting tile data.
///  @param[out]	outDataSize	The size of the tile data array.
/// @return True if the tile data was successfully created.
/// 根据提供的tile创建数据构建导航网格tile数据。
///  @param[in]		params		Tile 创建数据
///  @param[out]	outData		生成的图块数据。
///  @param[out]	outDataSize	图块数据数组的大小。
/// @return 如果已成功创建图块数据，则为 True。
bool dtCreateNavMeshData(dtNavMeshCreateParams* params, unsigned char** outData, int* outDataSize);

/// Swaps the endianess of the tile data's header (#dtMeshHeader).
///  @param[in,out]	data		The tile data array.
///  @param[in]		dataSize	The size of the data array.
/// 交换图块数据标头的字节顺序 (#dtMeshHeader)。
///  @param[in,out]	data		tile数据数组。
///  @param[in]		dataSize	数据数组的大小。
bool dtNavMeshHeaderSwapEndian(unsigned char* data, const int dataSize);

/// Swaps endianess of the tile data.
///  @param[in,out]	data		The tile data array.
///  @param[in]		dataSize	The size of the data array.
/// 交换图块数据的字节顺序。
///  @param[in,out]	data		tile数据数组。
///  @param[in]		dataSize	数据数组的大小。
bool dtNavMeshDataSwapEndian(unsigned char* data, const int dataSize);

#endif // DETOURNAVMESHBUILDER_H

// This section contains detailed documentation for members that don't have
// a source file. It reduces clutter in the main section of the header.
// 本节包含针对没有源文件的成员的详细文档。 它减少了标题主要部分的混乱。

/**

@struct dtNavMeshCreateParams
@par

This structure is used to marshal data between the Recast mesh generation pipeline and Detour navigation components.
此结构用于在 Recast 网格生成管道和 Detour 导航组件之间编组数据。

See the rcPolyMesh and rcPolyMeshDetail documentation for detailed information related to mesh structure.
有关网格结构的详细信息，请参阅 rcPolyMesh 和 rcPolyMeshDetail 文档。

Units are usually in voxels (vx) or world units (wu). The units for voxels, grid size, and cell size
are all based on the values of #cs and #ch.
单位通常采用体素 (vx) 或世界单位 (wu)。 体素、网格大小和像元大小的单位均基于 #cs 和 #ch 的值。

The standard navigation mesh build process is to create tile data using dtCreateNavMeshData, then add the tile
to a navigation mesh using either the dtNavMesh single tile <tt>init()</tt> function or the dtNavMesh::addTile()
function.
标准导航网格构建过程是使用 dtCreateNavMeshData 创建图块数据，然后使用 dtNavMesh 单图块 <tt>init()</tt> 函数或
dtNavMesh::addTile() 函数将图块添加到导航网格。

@see dtCreateNavMeshData

*/

