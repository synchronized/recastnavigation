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

#ifndef DETOURNAVMESH_H
#define DETOURNAVMESH_H

#include "DetourAlloc.h"
#include "DetourStatus.h"

// Undefine (or define in a build config) the following line to use 64bit polyref.
// Generally not needed, useful for very large worlds.
// Note: tiles build using 32bit refs are not compatible with 64bit refs!
//#define DT_POLYREF64 1

#ifdef DT_POLYREF64
// TODO: figure out a multiplatform version of uint64_t
// - maybe: https://code.google.com/p/msinttypes/
// - or: http://www.azillionmonkeys.com/qed/pstdint.h
#include <stdint.h>
#endif

// Note: If you want to use 64-bit refs, change the types of both dtPolyRef & dtTileRef.
// It is also recommended that you change dtHashRef() to a proper 64-bit hash.

/// A handle to a polygon within a navigation mesh tile.
/// @ingroup detour
/// 导航网格tile中多边形的句柄。
#ifdef DT_POLYREF64
static const unsigned int DT_SALT_BITS = 16;
static const unsigned int DT_TILE_BITS = 28;
static const unsigned int DT_POLY_BITS = 20;
typedef uint64_t dtPolyRef;
#else
typedef unsigned int dtPolyRef;
#endif

/// A handle to a tile within a navigation mesh.
/// @ingroup detour
/// 导航网格中tile的句柄。
#ifdef DT_POLYREF64
typedef uint64_t dtTileRef;
#else
typedef unsigned int dtTileRef;
#endif

/// The maximum number of vertices per navigation polygon.
/// @ingroup detour
/// 每个导航多边形的最大顶点数。
static const int DT_VERTS_PER_POLYGON = 6;

/// @{
/// @name Tile Serialization Constants
/// These constants are used to detect whether a navigation tile's data
/// and state format is compatible with the current build.
///

/// tile序列化常量
/// 这些常量用于检测导航tile的数据和状态格式是否与当前版本兼容。

/// A magic number used to detect compatibility of navigation tile data.
/// 用于检测导航tile数据兼容性的魔数。
static const int DT_NAVMESH_MAGIC = 'D'<<24 | 'N'<<16 | 'A'<<8 | 'V';

/// A version number used to detect compatibility of navigation tile data.
/// 用于检测导航tile数据兼容性的版本号。
static const int DT_NAVMESH_VERSION = 7;

/// A magic number used to detect the compatibility of navigation tile states.
/// 用于检测导航tile状态的兼容性的魔数。
static const int DT_NAVMESH_STATE_MAGIC = 'D'<<24 | 'N'<<16 | 'M'<<8 | 'S';

/// A version number used to detect compatibility of navigation tile states.
/// 用于检测导航tile状态兼容性的版本号。
static const int DT_NAVMESH_STATE_VERSION = 1;

/// @}

/// A flag that indicates that an entity links to an external entity.
/// (E.g. A polygon edge is a portal that links to another polygon.)
/// 指示实体链接到外部实体的标志。(例如，多边形边缘是连接到另一个多边形的入口。)
static const unsigned short DT_EXT_LINK = 0x8000;

/// A value that indicates the entity does not link to anything.
/// 表示实体不链接到任何东西的值。
static const unsigned int DT_NULL_LINK = 0xffffffff;

/// A flag that indicates that an off-mesh connection can be traversed in both directions. (Is bidirectional.)
/// 一个标志，表明一个网格外链接可以在两个方向上遍历。(双向)。
static const unsigned int DT_OFFMESH_CON_BIDIR = 1;

/// The maximum number of user defined area ids.
/// @ingroup detour
/// 用户自定义区域id的最大数量。
static const int DT_MAX_AREAS = 64;

/// Tile flags used for various functions and fields.
/// For an example, see dtNavMesh::addTile().
/// 用于各种函数和字段的平铺标志。示例请参见dtNavMesh::addTile()。
enum dtTileFlags
{
	/// The navigation mesh owns the tile memory and is responsible for freeing it.
	/// 导航网格拥有tile内存并负责释放它。
	DT_TILE_FREE_DATA = 0x01
};

/// Vertex flags returned by dtNavMeshQuery::findStraightPath.
/// 由dtNavMeshQuery::findStraightPath返回的顶点标志。
enum dtStraightPathFlags
{
	///< The vertex is the start position in the path.
	///< 顶点是路径的起始位置。
	DT_STRAIGHTPATH_START = 0x01,
	///< The vertex is the end position in the path.
	///< 顶点是路径的结束位置。
	DT_STRAIGHTPATH_END = 0x02,
	///< The vertex is the start of an off-mesh connection.
	///< 顶点是一个网格外链接的开始。
	DT_STRAIGHTPATH_OFFMESH_CONNECTION = 0x04
};

/// Options for dtNavMeshQuery::findStraightPath.
/// dtNavMeshQuery::findStraightPath的选项
enum dtStraightPathOptions
{
	///< Add a vertex at every polygon edge crossing where area changes.
	///< 在面积变化的每个多边形边缘交叉处添加一个顶点。
	DT_STRAIGHTPATH_AREA_CROSSINGS = 0x01,
	///< Add a vertex at every polygon edge crossing.
	///< 在每个多边形边缘交叉处添加一个顶点。
	DT_STRAIGHTPATH_ALL_CROSSINGS = 0x02
};


/// Options for dtNavMeshQuery::initSlicedFindPath and updateSlicedFindPath
/// dtNavMeshQuery::initSlicedFindPath和updateSlicedFindPath的选项
enum dtFindPathOptions
{
	///< use raycasts during pathfind to "shortcut" (raycast still consider costs)
	///< 在寻径过程中使用光线投射来“快捷”(光线投射仍然考虑成本)
	DT_FINDPATH_ANY_ANGLE	= 0x02
};

/// Options for dtNavMeshQuery::raycast
/// dtNavMeshQuery::raycast 的选项
enum dtRaycastOptions
{
	///< Raycast should calculate movement cost along the ray and fill RaycastHit::cost
	///< Raycast 应计算沿光线的移动成本并填充 RaycastHit::cost
	DT_RAYCAST_USE_COSTS = 0x01
};

enum dtDetailTriEdgeFlags
{
	///< Detail triangle edge is part of the poly boundary
	///< 细节三角形边是多边形边界的一部分
	DT_DETAIL_EDGE_BOUNDARY = 0x01
};


/// Limit raycasting during any angle pahfinding
/// The limit is given as a multiple of the character radius
/// 在任何角度寻路期间限制光线投射,限制以字符半径的倍数给出
static const float DT_RAY_CAST_LIMIT_PROPORTIONS = 50.0f;

/// Flags representing the type of a navigation mesh polygon.
/// 表示导航网格多边形类型的标志。
enum dtPolyTypes
{
	/// The polygon is a standard convex polygon that is part of the surface of the mesh.
	/// 该多边形是标准凸多边形，是网格表面的一部分。
	DT_POLYTYPE_GROUND = 0,
	/// The polygon is an off-mesh connection consisting of two vertices.
	/// 多边形是由两个顶点组成的网格外连接。
	DT_POLYTYPE_OFFMESH_CONNECTION = 1
};


/// Defines a polygon within a dtMeshTile object.
/// @ingroup detour
/// 在 dtMeshTile 对象内定义多边形。
struct dtPoly
{
	/// Index to first link in linked list. (Or #DT_NULL_LINK if there is no link.)
	/// 链表中第一个链接的索引。（如果没有链接，则为#DT_NULL_LINK。）
	/// 对应 dtMeshTile.links
	unsigned int firstLink;

	/// The indices of the polygon's vertices.
	/// The actual vertices are located in dtMeshTile::verts.
	/// 多边形顶点的索引。 实际的顶点位于 dtMeshTile::verts 中。
	unsigned short verts[DT_VERTS_PER_POLYGON];

	/// Packed data representing neighbor polygons references and flags for each edge.
	/// 表示每条边和相邻多边形references和flags。
	unsigned short neis[DT_VERTS_PER_POLYGON];

	/// The user defined polygon flags.
	/// 用户定义的多边形flags。
	unsigned short flags;

	/// The number of vertices in the polygon.
	/// 多边形中的顶点数。
	unsigned char vertCount;

	/// The bit packed area id and polygon type.
	/// @note Use the structure's set and get methods to access this value.
	/// 位封装区域 ID 和多边形类型。
	/// 注意 使用结构体的 set 和 get 方法来访问该值。
	unsigned char areaAndtype;

	/// Sets the user defined area id. [Limit: < #DT_MAX_AREAS]
	/// 设置用户定义的区域 ID。 [限制：< #DT_MAX_AREAS]
	inline void setArea(unsigned char a) { areaAndtype = (areaAndtype & 0xc0) | (a & 0x3f); }

	/// Sets the polygon type. (See: #dtPolyTypes.)
	/// 设置多边形类型。 （参见：#dtPolyTypes。）
	inline void setType(unsigned char t) { areaAndtype = (areaAndtype & 0x3f) | (t << 6); }

	/// Gets the user defined area id.
	/// 获取用户定义的区域 ID。
	inline unsigned char getArea() const { return areaAndtype & 0x3f; }

	/// Gets the polygon type. (See: #dtPolyTypes)
	/// 获取多边形类型。 （参见：#dtPolyTypes）
	inline unsigned char getType() const { return areaAndtype >> 6; }
};

/// Defines the location of detail sub-mesh data within a dtMeshTile.
/// 定义 dtMeshTile 中细节子网格数据的位置。
struct dtPolyDetail
{
	///< The offset of the vertices in the dtMeshTile::detailVerts array.
	///< dtMeshTile::detailVerts 数组中顶点的偏移量。
	unsigned int vertBase;
	///< The offset of the triangles in the dtMeshTile::detailTris array.
	///< dtMeshTile::detailTris 数组中三角形的偏移量。
	unsigned int triBase;
	///< The number of vertices in the sub-mesh.
	///< 子网格中的顶点数。
	unsigned char vertCount;
	///< The number of triangles in the sub-mesh.
	///< 子网格中三角形的数量。
	unsigned char triCount;
};

/// Defines a link between polygons.
/// @note This structure is rarely if ever used by the end user.
/// @see dtMeshTile
/// 定义多边形之间的链接。
/// 注意 此结构很少被最终用户使用。
struct dtLink
{
	///< Neighbour reference. (The neighbor that is linked to.)
	///< 邻居引用。（链接到(对方)的邻居。）
	dtPolyRef ref;
	///< Index of the next link.
	///< 下一个链接的索引。(dtMeshTile.links索引)
	unsigned int next;
	///< Index of the polygon edge that owns this link.
	///< 拥有此链接的多边形边的索引。(第几条边,从0开始)
	unsigned char edge;
	///< If a boundary link, defines on which side the link is.
	///< 如果是边界链接，则定义该链接位于哪一侧。(0-7,-1 右侧为0顺时针)
	///< 如果相邻的 poly 属于和当前 poly 位于属于同一个 MeshTile，那么 side的值则为 0xff
	///< 如果相邻的 poly 位于当前 poly 所在的 MeshTile 的周围 XZ 轴 4 个方向，则 side分别为 0(+x) 2(+z) 4(-x) 6(-z)
	unsigned char side;
	///< If a boundary link, defines the minimum sub-edge area.
	///< 如果是边界链接，则定义最小子边缘区域。
	///< 如果相邻的 poly 属于和当前 poly 位于属于同一个 MeshTile，那么 bmin和 bmax都为 0
	unsigned char bmin;
	///< If a boundary link, defines the maximum sub-edge area.
	///< 如果是边界链接，则定义最大子边缘区域。
	///< 如果相邻的 poly 位于当前 poly 所在的 MeshTile 的周围 XZ 轴 4 个方向，bmin和 bmax是以边的起始顶点为相对原点的比例（长度一样，则比值为 255）。
	unsigned char bmax;
};

/// Bounding volume node.
/// @note This structure is rarely if ever used by the end user.
/// @see dtMeshTile
/// 包围体节点
/// 注意: 最终用户很少使用这种结构。
struct dtBVNode
{
	///< Minimum bounds of the node's AABB. [(x, y, z)]
	///< 节点 AABB 的最小边界。 [（x，y，z）]
	unsigned short bmin[3];
	///< Maximum bounds of the node's AABB. [(x, y, z)]
	///< 节点 AABB 的最大边界。 [（x，y，z）]
	unsigned short bmax[3];
	///< The node's index. (Negative for escape sequence.)
	///< 节点的索引。（转义序列为负。）
	int i;
};

/// Defines an navigation mesh off-mesh connection within a dtMeshTile object.
/// An off-mesh connection is a user defined traversable connection made up to two vertices.
/// 在 dtMeshTile 对象内定义导航网格体离网格连接。
/// 网格外连接是用户定义的由两个顶点组成的可遍历连接。
struct dtOffMeshConnection
{
	/// The endpoints of the connection. [(ax, ay, az, bx, by, bz)]
	/// 连接的端点。 [(ax, ay, az, bx, by, bz)]
	float pos[6];

	/// The radius of the endpoints. [Limit: >= 0]
	/// 端点的半径。 [限制：>= 0]
	float rad;

	/// The polygon reference of the connection within the tile.
	// tile内连接的多边形引用(索引)。
	unsigned short poly;

	/// Link flags.
	/// @note These are not the connection's user defined flags. Those are assigned via the
	/// connection's dtPoly definition. These are link flags used for internal purposes.
	/// 链接标志。
	/// 注意: 这些不是连接的用户定义标志。 这些是通过分配,连接的 dtPoly 定义。 这些是用于内部目的的链接标志。
	unsigned char flags;

	/// End point side.
	/// 终点侧。
	unsigned char side;

	/// The id of the offmesh connection. (User assigned when the navigation mesh is built.)
	/// 网络外连接的 ID。 （构建navigation mesh时由用户分配。）
	unsigned int userId;
};

/// Provides high level information related to a dtMeshTile object.
/// @ingroup detour
/// 提供与 dtMeshTile 对象相关的高级信息。
struct dtMeshHeader
{
	///< Tile magic number. (Used to identify the data format.)
	///< tile魔术数字。 （用于识别数据格式。）
	int magic;
	///< Tile data format version number.
	///< tile数据格式版本号。
	int version;
	///< The x-position of the tile within the dtNavMesh tile grid. (x, y, layer)
	///< dtNavMesh tile网格内图块的 x 位置。（x，y，层）
	int x;
	///< The y-position of the tile within the dtNavMesh tile grid. (x, y, layer)
	///< dtNavMesh tile网格内图块的 y 位置。（x，y，层）
	int y;
	///< The layer of the tile within the dtNavMesh tile grid. (x, y, layer)
	///< dtNavMesh tile网格内图块的平铺层。（x，y，层）
	int layer;
	///< The user defined id of the tile.
	///< 用户定义的图块 ID。
	unsigned int userId;
	///< The number of polygons in the tile.
	///< tile中多边形的数量。
	int polyCount;
	///< The number of vertices in the tile.
	///< tile中的顶点数量
	int vertCount;
	///< The number of allocated links.
	///< 分配的链接数。
	int maxLinkCount;
	///< The number of sub-meshes in the detail mesh.
	///< 细节网格中子网格的数量。
	int detailMeshCount;

	/// The number of unique vertices in the detail mesh. (In addition to the polygon vertices.)
	/// 细节网格中唯一顶点的数量。（除了多边形顶点。）
	int detailVertCount;

	///< The number of triangles in the detail mesh.
	///< 细节网格中三角形的数量。
	int detailTriCount;
	///< The number of bounding volume nodes. (Zero if bounding volumes are disabled.)
	///< 包围体节点的数量。（如果边界体积被禁用，则为零。）
	int bvNodeCount;
	///< The number of off-mesh connections.
	///< 网络外连接的数量。
	int offMeshConCount;
	///< The index of the first polygon which is an off-mesh connection.
	///< 作为网格外连接的第一个多边形的索引。
	int offMeshBase;
	///< The height of the agents using the tile.
	///< 使用tile的agents的高度。
	float walkableHeight;
	///< The radius of the agents using the tile.
	///< 使用tile的agents的半径。
	float walkableRadius;
	///< The maximum climb height of the agents using the tile.
	///< 使用该tile的agents的最大攀爬高度。
	float walkableClimb;
	///< The minimum bounds of the tile's AABB. [(x, y, z)]
	///< tile的 AABB 的最小边界。 [（x，y，z）]
	float bmin[3];
	///< The maximum bounds of the tile's AABB. [(x, y, z)]
	///< tile的AABB 的最大边界。[(x, y, z)]
	float bmax[3];

	/// The bounding volume quantization factor.
	/// 包围体量化因子。
	float bvQuantFactor;
};

/// Defines a navigation mesh tile.
/// @ingroup detour
/// 定义导航网格tile。
struct dtMeshTile
{
	///< Counter describing modifications to the tile.
	///< 描述对tile的修改的计数器。
	///< 记录 tile 的修改次数，主要是动态更新使用(默认从1开始)
	unsigned int salt;

	///< Index to the next free link.
	///< free link链表的头节点的索引。
	unsigned int linksFreeList;
	///< The tile header.
	///< tile标题。
	dtMeshHeader* header;
	///< The tile polygons. [Size: dtMeshHeader::polyCount]
	///< tile多边形。 [size：dtMeshHeader::polyCount]
	dtPoly* polys;
	///< The tile vertices. [(x, y, z) * dtMeshHeader::vertCount]
	///< tile顶点。 [(x, y, z) * dtMeshHeader::vertCount]
	float* verts;
	///< The tile links. [Size: dtMeshHeader::maxLinkCount]
	///< tile链接。 [size：dtMeshHeader::maxLinkCount]
	dtLink* links;
	///< The tile's detail sub-meshes. [Size: dtMeshHeader::detailMeshCount]
	///< title的细节子网。 [size ：dtMeshheader::detailMeshCount]
	dtPolyDetail* detailMeshes;

	/// The detail mesh's unique vertices. [(x, y, z) * dtMeshHeader::detailVertCount]
	/// 细节网格的唯一顶点。 [(x, y, z) * dtMeshHeader::detailVertCount]
	float* detailVerts;

	/// The detail mesh's triangles. [(vertA, vertB, vertC, triFlags) * dtMeshHeader::detailTriCount].
	/// See dtDetailTriEdgeFlags and dtGetDetailTriEdgeFlags.
	/// 细节网格三角形定点索引。 [(vertA、vertB、vertC、triFlags) * dtMes dtMeshHeader::detailTriCount]。
	/// 请参见 dtDetailTriEdgeFlags 和 dtGetDetailTriEdgeFlags。
	unsigned char* detailTris;

	/// The tile bounding volume nodes. [Size: dtMeshHeader::bvNodeCount]
	/// (Will be null if bounding volumes are disabled.)
	/// 平铺包围体节点。 [大小：dtMeshHeader::bvNodeCount]
	/// （如果禁用边界体积，则将为空。）
	dtBVNode* bvTree;

	///< The tile off-mesh connections. [Size: dtMeshHeader::offMeshConCount]
	///< 网格外连接
	dtOffMeshConnection* offMeshCons;

	///< The tile data. (Not directly accessed under normal situations.)
	///< tile数据。 （正常情况下不能直接访问。）
	unsigned char* data;
	///< Size of the tile data.
	///< tile数据大小
	int dataSize;
	///< Tile flags. (See: #dtTileFlags)
	int flags;
	///< The next free tile, or the next tile in the spatial grid.
	///< 下一个空闲tile，或空间网格中的下一个图块。
	dtMeshTile* next;
private:
	dtMeshTile(const dtMeshTile&);
	dtMeshTile& operator=(const dtMeshTile&);
};

/// Get flags for edge in detail triangle.
/// @param[in]	triFlags		The flags for the triangle (last component of detail vertices above).
/// @param[in]	edgeIndex		The index of the first vertex of the edge. For instance, if 0,
///								returns flags for edge AB.
/// 获取详细三角形中边的标志。
/// @param[in]	triFlags		三角形的标志（上面细节顶点的最后一个组件）。
/// @param[in]	edgeIndex		边的第一个顶点的索引。 例如，如果为 0，则返回边 AB 的标志。
inline int dtGetDetailTriEdgeFlags(unsigned char triFlags, int edgeIndex)
{
	return (triFlags >> (edgeIndex * 2)) & 0x3;
}

/// Configuration parameters used to define multi-tile navigation meshes.
/// The values are used to allocate space during the initialization of a navigation mesh.
/// @see dtNavMesh::init()
/// @ingroup detour
/// 用于定义多tile导航网格的配置参数。
/// 这些值用于在导航网格初始化期间分配空间。
struct dtNavMeshParams
{
	///< The world space origin of the navigation mesh's tile space. [(x, y, z)]
	///< 导航网格体tile空间的世界空间原点。 [（x，y，z）]
	float orig[3];
	///< The width of each tile. (Along the x-axis.)
	///< 每块tile的宽度。（沿 x 轴。）
	float tileWidth;
	///< The height of each tile. (Along the z-axis.)
	///< 每块tile的高度。（沿 z 轴。）
	float tileHeight;
	///< The maximum number of tiles the navigation mesh can contain. This and maxPolys are used to calculate how many bits are needed to identify tiles and polygons uniquely.
	///< 导航网格可以包含的最大tile数量。 它和 maxPolys 用于计算需要多少位来唯一地识别tile和多边形。
	int maxTiles;
	///< The maximum number of polygons each tile can contain. This and maxTiles are used to calculate how many bits are needed to identify tiles and polygons uniquely.
	///< 每个tile可以包含的最大多边形数。 它和 maxTiles 用于计算需要多少位来唯一地识别tile和多边形。
	int maxPolys;
};

/// A navigation mesh based on tiles of convex polygons.
/// @ingroup detour
/// 基于凸多边形图块的导航网格。
class dtNavMesh
{
public:
	dtNavMesh();
	~dtNavMesh();

	/// @{
	/// @name Initialization and Tile Management

	/// Initializes the navigation mesh for tiled use.
	///  @param[in]	params		Initialization parameters.
	/// @return The status flags for the operation.
	/// 初始化导航网格以供tile使用。
	///  @param[in]	params		初始化参数。
	/// @return 操作的状态标志。
	dtStatus init(const dtNavMeshParams* params);

	/// Initializes the navigation mesh for single tile use.
	///  @param[in]	data		Data of the new tile. (See: #dtCreateNavMeshData)
	///  @param[in]	dataSize	The data size of the new tile.
	///  @param[in]	flags		The tile flags. (See: #dtTileFlags)
	/// @return The status flags for the operation.
	///  @see dtCreateNavMeshData
	/// 初始化导航网格以供单个tile使用。
	///  @param[in]	data		新tile的数据。 （参见：#dtCreateNavMeshData）
	///  @param[in]	dataSize	新图块的数据大小。
	///  @param[in]	flags		tile标志。（参见：#dtTileFlags）
	/// @return 操作的状态标志。
	dtStatus init(unsigned char* data, const int dataSize, const int flags);

	/// The navigation mesh initialization params.
	/// 导航网格初始化参数。
	const dtNavMeshParams* getParams() const;

	/// Adds a tile to the navigation mesh.
	///  @param[in]		data		Data for the new tile mesh. (See: #dtCreateNavMeshData)
	///  @param[in]		dataSize	Data size of the new tile mesh.
	///  @param[in]		flags		Tile flags. (See: #dtTileFlags)
	///  @param[in]		lastRef		The desired reference for the tile. (When reloading a tile.) [opt] [Default: 0]
	///  @param[out]	result		The tile reference. (If the tile was succesfully added.) [opt]
	/// @return The status flags for the operation.
	/// 将tile添加到导航网格。
	///  @param[in]		data		新图块网格的数据。 （参见：#dtCreateNavMeshData）
	///  @param[in]		dataSize	新图块网格的数据大小。
	///  @param[in]		flags		平铺标志。 （参见：#dtTileFlags）
	///  @param[in]		lastRef     所需的tile参考。 （重新加载tile时。） [可选] [默认值：0]
	///  @param[out]	result		tile引用。（如果图块已成功添加。）
	/// @return 操作的状态标志
	dtStatus addTile(unsigned char* data, int dataSize, int flags, dtTileRef lastRef, dtTileRef* result);

	/// Removes the specified tile from the navigation mesh.
	///  @param[in]		ref			The reference of the tile to remove.
	///  @param[out]	data		Data associated with deleted tile.
	///  @param[out]	dataSize	Size of the data associated with deleted tile.
	/// @return The status flags for the operation.
	/// 操作的状态标志。
	///  @param[in]		ref			要删除的tile的参考。
	///  @param[out]	data		与已删除tile关联的数据。
	///  @param[out]	dataSize	与已删除tile关联的数据大小。
	/// @return 操作的状态标志。
	dtStatus removeTile(dtTileRef ref, unsigned char** data, int* dataSize);

	/// @}

	/// @{
	/// @name Query Functions

	/// Calculates the tile grid location for the specified world position.
	///  @param[in]	pos  The world position for the query. [(x, y, z)]
	///  @param[out]	tx		The tile's x-location. (x, y)
	///  @param[out]	ty		The tile's y-location. (x, y)
	/// 计算指定世界位置的tile网格位置。
	///  @param[in]  pos 查询的世界位置。 [(x，y，z)]
	///  @param[out] tx  tile的 x 位置。 (x，y)
	///  @param[out] ty  tile的 y 位置。 (x，y)
	void calcTileLoc(const float* pos, int* tx, int* ty) const;

	/// Gets the tile at the specified grid location.
	///  @param[in]	x		The tile's x-location. (x, y, layer)
	///  @param[in]	y		The tile's y-location. (x, y, layer)
	///  @param[in]	layer	The tile's layer. (x, y, layer)
	/// @return The tile, or null if the tile does not exist.
	/// 获取指定网格位置处的tile。
	///  @param[in]	x		tile的 x 位置。 （x，y，层）
	///  @param[in]	y		tile的 y 位置。 （x，y，层）
	///  @param[in]	layer	tile的层。 （x，y，层）
	/// @return 返回tile，如果查找的tile不存在则为 null。
	const dtMeshTile* getTileAt(const int x, const int y, const int layer) const;

	/// Gets all tiles at the specified grid location. (All layers.)
	///  @param[in]		x			The tile's x-location. (x, y)
	///  @param[in]		y			The tile's y-location. (x, y)
	///  @param[out]	tiles		A pointer to an array of tiles that will hold the result.
	///  @param[in]		maxTiles	The maximum tiles the tiles parameter can hold.
	/// @return The number of tiles returned in the tiles array.
	/// 获取指定网格位置的所有图块。 （所有层。）
	///  @param[in]		x			tile的 x 位置。 （x，y）
	///  @param[in]		y			tile的 y 位置。 （x，y）
	///  @param[out]	tiles		指向将保存结果的tile数组的指针。
	///  @param[in]		maxTiles	tiles 参数可以容纳的最大tile数。
	/// @return 在tiles数组中返回的tiles的数量。
	int getTilesAt(const int x, const int y,
				   dtMeshTile const** tiles, const int maxTiles) const;

	/// Gets the tile reference for the tile at specified grid location.
	///  @param[in]	x		The tile's x-location. (x, y, layer)
	///  @param[in]	y		The tile's y-location. (x, y, layer)
	///  @param[in]	layer	The tile's layer. (x, y, layer)
	/// @return The tile reference of the tile, or 0 if there is none.
	/// 获取指定网格位置处图块的tile引用。
	///  @param[in]	x		tile的x位置。(x，y，层)
	///  @param[in]	y		tile的y位置。 (x，y，层)
	///  @param[in]	layer	tile的层。(x，y，层)
	/// @return tile的图块引用，如果没有则为 0。
	dtTileRef getTileRefAt(int x, int y, int layer) const;

	/// Gets the tile reference for the specified tile.
	///  @param[in]	tile	The tile.
	/// @return The tile reference of the tile.
	/// 获取指定图块的图块引用。
	///  @param[in]	tile	tile.
	/// @return tile引用
	dtTileRef getTileRef(const dtMeshTile* tile) const;

	/// Gets the tile for the specified tile reference.
	///  @param[in]	ref		The tile reference of the tile to retrieve.
	/// @return The tile for the specified reference, or null if the
	///		reference is invalid.
	/// 获取指定tile引用的图块。
	///  @param[in]	ref		要检索的tile的图块引用。
	/// @return 指定引用的tile，如果引用无效，则为 null。
	const dtMeshTile* getTileByRef(dtTileRef ref) const;

	/// The maximum number of tiles supported by the navigation mesh.
	/// @return The maximum number of tiles supported by the navigation mesh.
	/// 导航网格支持的最大tile数量。
	/// @return 导航网格支持的最大tile数量。
	int getMaxTiles() const;

	/// Gets the tile at the specified index.
	///  @param[in]	i		The tile index. [Limit: 0 >= index < #getMaxTiles()]
	/// @return The tile at the specified index.
	/// 获取指定索引处的tile。
	///  @param[in]	i		tile索引。 [限制：0 >= 索引 < #getMaxTiles()]
	/// @return 指定索引处的tile。
	const dtMeshTile* getTile(int i) const;

	/// Gets the tile and polygon for the specified polygon reference.
	///  @param[in]		ref		The reference for the a polygon.
	///  @param[out]	tile	The tile containing the polygon.
	///  @param[out]	poly	The polygon.
	/// @return The status flags for the operation.
	/// 获取指定多边形参考的tile和多边形。
	///  @param[in]		ref		多边形的引用。
	///  @param[out]	tile	包含多边形的tile。
	///  @param[out]	poly	多边形.
	/// @return 操作的状态标志。
	dtStatus getTileAndPolyByRef(const dtPolyRef ref, const dtMeshTile** tile, const dtPoly** poly) const;

	/// Returns the tile and polygon for the specified polygon reference.
	///  @param[in]		ref		A known valid reference for a polygon.
	///  @param[out]	tile	The tile containing the polygon.
	///  @param[out]	poly	The polygon.
	/// 返回指定多边形参考的图块和多边形。
	///  @param[in]		ref		多边形的已知有效引用。
	///  @param[out]	tile	包含多边形的tile.
	///  @param[out]	poly	多边形.
	void getTileAndPolyByRefUnsafe(const dtPolyRef ref, const dtMeshTile** tile, const dtPoly** poly) const;

	/// Checks the validity of a polygon reference.
	///  @param[in]	ref		The polygon reference to check.
	/// @return True if polygon reference is valid for the navigation mesh.
	/// 检查多边形引用的有效性。
	///  @param[in]	ref		要检查的多边形的引用.
	/// @return 如果多边形引用对于导航网格有效，则为 True.
	bool isValidPolyRef(dtPolyRef ref) const;

	/// Gets the polygon reference for the tile's base polygon.
	///  @param[in]	tile		The tile.
	/// @return The polygon reference for the base polygon in the specified tile.
	/// 获取tile的多边形的base polygon reference。(第0个多边形引用)
	///  @param[in]	tile		tile.
	/// @return 指定tile中基本多边形的多边形引用。
	dtPolyRef getPolyRefBase(const dtMeshTile* tile) const;

	/// Gets the endpoints for an off-mesh connection, ordered by "direction of travel".
	///  @param[in]		prevRef		The reference of the polygon before the connection.
	///  @param[in]		polyRef		The reference of the off-mesh connection polygon.
	///  @param[out]	startPos	The start position of the off-mesh connection. [(x, y, z)]
	///  @param[out]	endPos		The end position of the off-mesh connection. [(x, y, z)]
	/// @return The status flags for the operation.
	/// 获取网络外连接的端点，按“行进方向”排序。
	///  @param[in]		prevRef		连接前多边形的引用。
	///  @param[in]		polyRef		网格外连接多边形的引用。
	///  @param[out]	startPos	网络外连接的起始位置。[(x，y，z)]
	///  @param[out]	endPos		网络外连接的结束位置。[(x，y，z)]
	/// @return 操作的状态标志。
	dtStatus getOffMeshConnectionPolyEndPoints(dtPolyRef prevRef, dtPolyRef polyRef, float* startPos, float* endPos) const;

	/// Gets the specified off-mesh connection.
	///  @param[in]	ref		The polygon reference of the off-mesh connection.
	/// @return The specified off-mesh connection, or null if the polygon reference is not valid.
	/// 获取指定的网络外连接。
	///  @param[in]	ref		网络外连接的多边形引用。
	/// @return 指定的网络外连接，如果多边形引用无效，则为 null。
	const dtOffMeshConnection* getOffMeshConnectionByRef(dtPolyRef ref) const;

	/// @}

	/// @{
	/// @name State Management
	/// These functions do not effect #dtTileRef or #dtPolyRef's.
	/// 这些函数不会影响#dtTileRef 或#dtPolyRef。

	/// Sets the user defined flags for the specified polygon.
	///  @param[in]	ref		The polygon reference.
	///  @param[in]	flags	The new flags for the polygon.
	/// @return The status flags for the operation.
	/// 为指定的多边形设置用户定义的标志。
	///  @param[in]	ref		多边形引用
	///  @param[in]	flags	多边形的新标志。
	/// @return 操作的状态标志。
	dtStatus setPolyFlags(dtPolyRef ref, unsigned short flags);

	/// Gets the user defined flags for the specified polygon.
	///  @param[in]		ref				The polygon reference.
	///  @param[out]	resultFlags		The polygon flags.
	/// @return The status flags for the operation.
	/// 获取指定多边形的用户定义标志。
	///  @param[in]		ref				多边形应用。
	///  @param[out]	resultFlags		多边形标志。
	/// @return 操作的状态标志。
	dtStatus getPolyFlags(dtPolyRef ref, unsigned short* resultFlags) const;

	/// Sets the user defined area for the specified polygon.
	///  @param[in]	ref		The polygon reference.
	///  @param[in]	area	The new area id for the polygon. [Limit: < #DT_MAX_AREAS]
	/// @return The status flags for the operation.
	/// 设置指定多边形的用户定义区域。
	///  @param[in]	ref		多边形引用。
	///  @param[in]	area	多边形的新区域 ID。[限制：< #DT_MAX_AREAS]
	/// @return 操作的状态标志。
	dtStatus setPolyArea(dtPolyRef ref, unsigned char area);

	/// Gets the user defined area for the specified polygon.
	///  @param[in]		ref			The polygon reference.
	///  @param[out]	resultArea	The area id for the polygon.
	/// @return The status flags for the operation.
	/// 获取指定多边形的用户定义区域。
	///  @param[in]		ref			多边形引用
	///  @param[out]	resultArea	多边形的区域 ID。
	/// @return 操作的状态标志。
	dtStatus getPolyArea(dtPolyRef ref, unsigned char* resultArea) const;

	/// Gets the size of the buffer required by #storeTileState to store the specified tile's state.
	///  @param[in]	tile	The tile.
	/// @return The size of the buffer required to store the state.
	/// 获取 #storeTileState 存储指定tile状态所需的缓冲区大小。
	///  @param[in]	tile	tile.
	/// @return 存储状态所需的缓冲区大小。
	int getTileStateSize(const dtMeshTile* tile) const;

	/// Stores the non-structural state of the tile in the specified buffer. (Flags, area ids, etc.)
	///  @param[in]		tile			The tile.
	///  @param[out]	data			The buffer to store the tile's state in.
	///  @param[in]		maxDataSize		The size of the data buffer. [Limit: >= #getTileStateSize]
	/// @return The status flags for the operation.
	/// 将tile的非结构状态存储在指定的缓冲区中。（旗帜、区域 ID 等）
	///  @param[in]		tile			tile.
	///  @param[out]	data			用于存储tile状态的缓冲区。
	///  @param[in]		maxDataSize		数据缓冲区的大小。 [限制：>= #getTileStateSize]
	/// @return 操作的状态标志。
	dtStatus storeTileState(const dtMeshTile* tile, unsigned char* data, const int maxDataSize) const;

	/// Restores the state of the tile.
	///  @param[in]	tile			The tile.
	///  @param[in]	data			The new state. (Obtained from #storeTileState.)
	///  @param[in]	maxDataSize		The size of the state within the data buffer.
	/// @return The status flags for the operation.
	/// 恢复tile的状态。
	///  @param[in]	tile			tile.
	///  @param[in]	data			新的状态。 （从#storeTileState 获取。）
	///  @param[in]	maxDataSize		数据缓冲区内的状态大小。
	/// @return 操作的状态标志。
	dtStatus restoreTileState(dtMeshTile* tile, const unsigned char* data, const int maxDataSize);

	/// @}

	/// @{
	/// @name Encoding and Decoding
	/// These functions are generally meant for internal use only.
	/// 这些功能通常仅供内部使用。

	/// Derives a standard polygon reference.
	///  @note This function is generally meant for internal use only.
	///  @param[in]	salt	The tile's salt value.
	///  @param[in]	it		The index of the tile.
	///  @param[in]	ip		The index of the polygon within the tile.
	/// 导出标准多边形引用。
	///  @note 此功能通常仅供内部使用。
	///  @param[in]	salt	tile的盐值。
	///  @param[in]	it		tile的索引。
	///  @param[in]	ip		tile内多边形的索引。
	inline dtPolyRef encodePolyId(unsigned int salt, unsigned int it, unsigned int ip) const
	{
#ifdef DT_POLYREF64
		return ((dtPolyRef)salt << (DT_POLY_BITS+DT_TILE_BITS)) | ((dtPolyRef)it << DT_POLY_BITS) | (dtPolyRef)ip;
#else
		return ((dtPolyRef)salt << (m_polyBits+m_tileBits)) | ((dtPolyRef)it << m_polyBits) | (dtPolyRef)ip;
#endif
	}

	/// Decodes a standard polygon reference.
	///  @note This function is generally meant for internal use only.
	///  @param[in]	ref   The polygon reference to decode.
	///  @param[out]	salt	The tile's salt value.
	///  @param[out]	it		The index of the tile.
	///  @param[out]	ip		The index of the polygon within the tile.
	///  @see #encodePolyId
	/// 解码标准多边形引用。
	///  @note 此功能通常仅供内部使用。
	///  @param[in]	ref   要解码的多边形引用。
	///  @param[out]	salt	tile的盐值。
	///  @param[out]	it		tile的索引。
	///  @param[out]	ip		tile内多边形的索引。
	inline void decodePolyId(dtPolyRef ref, unsigned int& salt, unsigned int& it, unsigned int& ip) const
	{
#ifdef DT_POLYREF64
		const dtPolyRef saltMask = ((dtPolyRef)1<<DT_SALT_BITS)-1;
		const dtPolyRef tileMask = ((dtPolyRef)1<<DT_TILE_BITS)-1;
		const dtPolyRef polyMask = ((dtPolyRef)1<<DT_POLY_BITS)-1;
		salt = (unsigned int)((ref >> (DT_POLY_BITS+DT_TILE_BITS)) & saltMask);
		it = (unsigned int)((ref >> DT_POLY_BITS) & tileMask);
		ip = (unsigned int)(ref & polyMask);
#else
		const dtPolyRef saltMask = ((dtPolyRef)1<<m_saltBits)-1;
		const dtPolyRef tileMask = ((dtPolyRef)1<<m_tileBits)-1;
		const dtPolyRef polyMask = ((dtPolyRef)1<<m_polyBits)-1;
		salt = (unsigned int)((ref >> (m_polyBits+m_tileBits)) & saltMask);
		it = (unsigned int)((ref >> m_polyBits) & tileMask);
		ip = (unsigned int)(ref & polyMask);
#endif
	}

	/// Extracts a tile's salt value from the specified polygon reference.
	///  @note This function is generally meant for internal use only.
	///  @param[in]	ref		The polygon reference.
	///  @see #encodePolyId
	/// 从指定的多边形引用中提取tile的盐值。
	///  @note 此功能通常仅供内部使用。
	///  @param[in]	ref		多边形引用。
	inline unsigned int decodePolyIdSalt(dtPolyRef ref) const
	{
#ifdef DT_POLYREF64
		const dtPolyRef saltMask = ((dtPolyRef)1<<DT_SALT_BITS)-1;
		return (unsigned int)((ref >> (DT_POLY_BITS+DT_TILE_BITS)) & saltMask);
#else
		const dtPolyRef saltMask = ((dtPolyRef)1<<m_saltBits)-1;
		return (unsigned int)((ref >> (m_polyBits+m_tileBits)) & saltMask);
#endif
	}

	/// Extracts the tile's index from the specified polygon reference.
	///  @note This function is generally meant for internal use only.
	///  @param[in]	ref		The polygon reference.
	///  @see #encodePolyId
	/// 从指定的多边形引用中提取tile的索引。
	///  @note 此功能通常仅供内部使用。
	///  @param[in]	ref		多边形引用。
	inline unsigned int decodePolyIdTile(dtPolyRef ref) const
	{
#ifdef DT_POLYREF64
		const dtPolyRef tileMask = ((dtPolyRef)1<<DT_TILE_BITS)-1;
		return (unsigned int)((ref >> DT_POLY_BITS) & tileMask);
#else
		const dtPolyRef tileMask = ((dtPolyRef)1<<m_tileBits)-1;
		return (unsigned int)((ref >> m_polyBits) & tileMask);
#endif
	}

	/// Extracts the polygon's index (within its tile) from the specified polygon reference.
	///  @note This function is generally meant for internal use only.
	///  @param[in]	ref		The polygon reference.
	///  @see #encodePolyId
	/// 从指定的多边形引用中提取多边形的索引（在其tile内）。
	///  @note 此功能通常仅供内部使用。
	///  @param[in]	ref		多边形引用。
	inline unsigned int decodePolyIdPoly(dtPolyRef ref) const
	{
#ifdef DT_POLYREF64
		const dtPolyRef polyMask = ((dtPolyRef)1<<DT_POLY_BITS)-1;
		return (unsigned int)(ref & polyMask);
#else
		const dtPolyRef polyMask = ((dtPolyRef)1<<m_polyBits)-1;
		return (unsigned int)(ref & polyMask);
#endif
	}

	/// @}

private:
	// Explicitly disabled copy constructor and copy assignment operator.
	// 显式禁用复制构造函数和复制赋值运算符。
	dtNavMesh(const dtNavMesh&);
	dtNavMesh& operator=(const dtNavMesh&);

	/// Returns pointer to tile in the tile array.
	/// 返回指向tile 数组中tile 的指针。
	dtMeshTile* getTile(int i);

	/// Returns neighbour tile based on side.
	/// 根据边返回相邻tile。
	int getTilesAt(const int x, const int y,
				   dtMeshTile** tiles, const int maxTiles) const;

	/// Returns neighbour tile based on side.
	/// 根据边返回相邻tile。
	int getNeighbourTilesAt(const int x, const int y, const int side,
							dtMeshTile** tiles, const int maxTiles) const;

	/// Returns all polygons in neighbour tile based on portal defined by the segment.
	/// 根据段定义的入口返回相邻tile中的所有多边形。
	int findConnectingPolys(const float* va, const float* vb,
							const dtMeshTile* tile, int side,
							dtPolyRef* con, float* conarea, int maxcon) const;

	/// Builds internal polygons links for a tile.
	/// 为tile构建内部多边形链接。
	void connectIntLinks(dtMeshTile* tile);
	/// Builds internal polygons links for a tile.
	/// 为tile构建内部多边形链接。
	void baseOffMeshLinks(dtMeshTile* tile);

	/// Builds external polygon links for a tile.
	/// 为tile构建外部多边形链接。
	void connectExtLinks(dtMeshTile* tile, dtMeshTile* target, int side);
	/// Builds external polygon links for a tile.
	/// 为tile构建外部多边形链接。
	void connectExtOffMeshLinks(dtMeshTile* tile, dtMeshTile* target, int side);

	/// Removes external links at specified side.
	/// 删除指定一侧的外部链接。
	void unconnectLinks(dtMeshTile* tile, dtMeshTile* target);


	// TODO: These methods are duplicates from dtNavMeshQuery, but are needed for off-mesh connection finding.
	// 这些方法与 dtNavMeshQuery 重复，但需要用于离网连接查找。

	/// Queries polygons within a tile.
	/// 查询图块内的多边形。
	int queryPolygonsInTile(const dtMeshTile* tile, const float* qmin, const float* qmax,
							dtPolyRef* polys, const int maxPolys) const;
	/// Find nearest polygon within a tile.
	/// 查找tile内最近的多边形。
	dtPolyRef findNearestPolyInTile(const dtMeshTile* tile, const float* center,
									const float* halfExtents, float* nearestPt) const;
	/// Returns whether position is over the poly and the height at the position if so.
	/// 返回位置是否在多边形上方，如果是，则返回该位置的高度。
	bool getPolyHeight(const dtMeshTile* tile, const dtPoly* poly, const float* pos, float* height) const;
	/// Returns closest point on polygon.
	/// 返回多边形上最近的点。
	void closestPointOnPoly(dtPolyRef ref, const float* pos, float* closest, bool* posOverPoly) const;

	///< Current initialization params. TODO: do not store this info twice.
	///< 当前初始化参数。 TODO：不要将此信息存储两次。
	dtNavMeshParams m_params;
	///< Origin of the tile (0,0)
	///< tile的原点 (0,0)
	float m_orig[3];
	///< Dimensions of each tile.
	///< 每块tile的尺寸。
	float m_tileWidth, m_tileHeight;
	///< Max number of tiles.
	///< 最大tile数量。
	int m_maxTiles;
	///< Tile hash lookup size (must be pot).
	///< 平铺哈希查找大小（必须是 pot）。
	int m_tileLutSize;
	///< Tile hash lookup mask.
	///< tile哈希查找掩码。
	int m_tileLutMask;

	///< Tile hash lookup.
	///< tile哈希表链表头节点
	dtMeshTile** m_posLookup;
	///< Freelist of tiles.
	///< tile的空闲列表链表头节点。
	dtMeshTile* m_nextFree;
	///< List of tiles.
	///< tile列表。
	dtMeshTile* m_tiles;

#ifndef DT_POLYREF64
	///< tile ID 中的salt位数。
	unsigned int m_saltBits;
	///< tile ID 中的tile位数。
	unsigned int m_tileBits;
	///< tile ID 中的多边形位数。
	unsigned int m_polyBits;
#endif

	friend class dtNavMeshQuery;
};

/// Allocates a navigation mesh object using the Detour allocator.
/// @return A navigation mesh that is ready for initialization, or null on failure.
///  @ingroup detour
/// 使用 Detour 分配器分配导航网格对象。
/// @return 准备初始化的导航网格，失败时返回 null。
dtNavMesh* dtAllocNavMesh();

/// Frees the specified navigation mesh object using the Detour allocator.
///  @param[in]	navmesh		A navigation mesh allocated using #dtAllocNavMesh
///  @ingroup detour
/// 使用 Detour 分配器释放指定的导航网格对象。
///  @param[in]	navmesh		使用 #dtAllocNavMesh 分配的导航网格
void dtFreeNavMesh(dtNavMesh* navmesh);

#endif // DETOURNAVMESH_H

///////////////////////////////////////////////////////////////////////////

// This section contains detailed documentation for members that don't have
// a source file. It reduces clutter in the main section of the header.
// 本节包含针对没有源文件的成员的详细文档。 它减少了标题主要部分的混乱。

/**

@typedef dtPolyRef
@par

Polygon references are subject to the same invalidate/preserve/restore
rules that apply to #dtTileRef's.  If the #dtTileRef for the polygon's
tile changes, the polygon reference becomes invalid.
多边形引用遵循适用于#dtTileRef 的相同无效/保留/恢复规则。
如果多边形tile的 #dtTileRef 发生更改，则多边形参考将变为无效。

Changing a polygon's flags, area id, etc. does not impact its polygon
reference.
更改多边形的标志、区域 ID 等不会影响其多边形引用。

@typedef dtTileRef
@par

The following changes will invalidate a tile reference:
以下更改将使tile引用无效：

- The referenced tile has been removed from the navigation mesh.
- 引用的图块已从导航网格中删除。
- The navigation mesh has been initialized using a different set
  of #dtNavMeshParams.
- 导航网格已使用一组不同的#dtNavMeshParams 进行初始化。

A tile reference is preserved/restored if the tile is added to a navigation
mesh initialized with the original #dtNavMeshParams and is added at the
original reference location. (E.g. The lastRef parameter is used with
dtNavMesh::addTile.)
如果将tile添加到使用原始 #dtNavMeshParams 初始化的导航网格并添加到原始引用位置，
则tile引用将被保留/恢复。 （例如，lastRef 参数与 dtNavMesh::addTile 一起使用。）

Basically, if the storage structure of a tile changes, its associated
tile reference changes.
基本上，如果tile的存储结构发生变化，则其关联的tile引用也会发生变化。

@var unsigned short dtPoly::neis[DT_VERTS_PER_POLYGON]
@par

Each entry represents data for the edge starting at the vertex of the same index.
E.g. The entry at index n represents the edge data for vertex[n] to vertex[n+1].
每个条目表示从同一索引的顶点开始的边的数据。 例如。 索引n处的条目表示vertex[n]到vertex[n+1]的边数据。

A value of zero indicates the edge has no polygon connection. (It makes up the
border of the navigation mesh.)
值为零表示边没有多边形连接。（它构成了导航网格的边框。）

The information can be extracted as follows:
信息可以提取如下：
@code
// Get the neighbor polygon reference.
// 获取相邻多边形引用。
neighborRef = neis[n] & 0xff;

if (neis[n] & #DT_EX_LINK)
{
    // The edge is an external (portal) edge.
	// 该边缘是外部（入口）边缘。
}
@endcode

@var float dtMeshHeader::bvQuantFactor
@par

This value is used for converting between world and bounding volume coordinates.
该值用于在世界坐标和包围体坐标之间进行转换。
For example:
@code
const float cs = 1.0f / tile->header->bvQuantFactor;
const dtBVNode* n = &tile->bvTree[i];
if (n->i >= 0)
{
    // This is a leaf node. 叶子节点
    float worldMinX = tile->header->bmin[0] + n->bmin[0]*cs;
    float worldMinY = tile->header->bmin[0] + n->bmin[1]*cs;
    // Etc...
}
@endcode

@struct dtMeshTile
@par

Tiles generally only exist within the context of a dtNavMesh object.
tile通常仅存在于 dtNavMesh 对象的上下文中。

Some tile content is optional.  For example, a tile may not contain any
off-mesh connections.  In this case the associated pointer will be null.
某些tile内容是可选的。 例如，tile可能不包含任何网络外连接。 在这种情况下，关联的指针将为空。

If a detail mesh exists it will share vertices with the base polygon mesh.
Only the vertices unique to the detail mesh will be stored in #detailVerts.
如果存在细节网格，它将与基础多边形网格共享顶点。
只有细节网格特有的顶点才会存储在#detailVerts 中。

@warning Tiles returned by a dtNavMesh object are not guarenteed to be populated.
For example: The tile at a location might not have been loaded yet, or may have been removed.
In this case, pointers will be null.  So if in doubt, check the polygon count in the
tile's header to determine if a tile has polygons defined.
@警告 dtNavMesh 对象返回的tile不保证被填充。 例如：某个位置的tile可能尚未加载，或者可能已被删除。
在这种情况下，指针将为空。 因此，如果有疑问，请检查tile头信息中的多边形计数，以确定tile是否已定义多边形。

@var float dtOffMeshConnection::pos[6]
@par

For a properly built navigation mesh, vertex A will always be within the bounds of the mesh.
Vertex B is not required to be within the bounds of the mesh.
对于正确构建的导航网格，顶点 A 将始终位于网格的边界内。顶点 B 不需要位于网格的边界内。

*/
