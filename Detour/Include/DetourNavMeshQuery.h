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

#ifndef DETOURNAVMESHQUERY_H
#define DETOURNAVMESHQUERY_H

#include "DetourNavMesh.h"
#include "DetourStatus.h"


// Define DT_VIRTUAL_QUERYFILTER if you wish to derive a custom filter from dtQueryFilter.
// On certain platforms indirect or virtual function call is expensive. The default
// setting is to use non-virtual functions, the actual implementations of the functions
// are declared as inline for maximum speed.
// 如果您希望从 dtQueryFilter 派生自定义过滤器，请定义 DT_VIRTUAL_QUERYFILTER。
// 在某些平台上，间接或虚拟函数调用的成本很高。 默认设置是使用非虚函数，函数的实际实现被声明为内联以获得最大速度。

//#define DT_VIRTUAL_QUERYFILTER 1

/// Defines polygon filtering and traversal costs for navigation mesh query operations.
/// 定义导航网格查询操作的多边形过滤和遍历成本。
/// @ingroup detour
class dtQueryFilter
{
	///< Cost per area type. (Used by default implementation.)
	///< 每个区域类型的成本。 （默认实现使用。）
	float m_areaCost[DT_MAX_AREAS];
	///< Flags for polygons that can be visited. (Used by default implementation.)
	///< 可以访问的多边形的标志。 （默认实现使用。）
	unsigned short m_includeFlags;
	///< Flags for polygons that should not be visited. (Used by default implementation.)
	///< 不应访问的多边形的标志。 （默认实现使用。）
	unsigned short m_excludeFlags;

public:
	dtQueryFilter();

#ifdef DT_VIRTUAL_QUERYFILTER
	virtual ~dtQueryFilter() { }
#endif

	/// Returns true if the polygon can be visited.  (I.e. Is traversable.)
	///  @param[in]		ref		The reference id of the polygon test.
	///  @param[in]		tile	The tile containing the polygon.
	///  @param[in]		poly  The polygon to test.
	/// 如果可以访问多边形则返回 true。 （即是可遍历的。）
	///  @param[in]		ref		多边形测试的参考ID。
	///  @param[in]		tile	包含多边形的tile。
	///  @param[in]		poly  要测试的多边形。
#ifdef DT_VIRTUAL_QUERYFILTER
	virtual bool passFilter(const dtPolyRef ref,
							const dtMeshTile* tile,
							const dtPoly* poly) const;
#else
	bool passFilter(const dtPolyRef ref,
					const dtMeshTile* tile,
					const dtPoly* poly) const;
#endif

	/// Returns cost to move from the beginning to the end of a line segment
	/// that is fully contained within a polygon.
	///  @param[in]		pa			The start position on the edge of the previous and current polygon. [(x, y, z)]
	///  @param[in]		pb			The end position on the edge of the current and next polygon. [(x, y, z)]
	///  @param[in]		prevRef		The reference id of the previous polygon. [opt]
	///  @param[in]		prevTile	The tile containing the previous polygon. [opt]
	///  @param[in]		prevPoly	The previous polygon. [opt]
	///  @param[in]		curRef		The reference id of the current polygon.
	///  @param[in]		curTile		The tile containing the current polygon.
	///  @param[in]		curPoly		The current polygon.
	///  @param[in]		nextRef		The refernece id of the next polygon. [opt]
	///  @param[in]		nextTile	The tile containing the next polygon. [opt]
	///  @param[in]		nextPoly	The next polygon. [opt]
	/// 返回从完全包含在多边形内的线段的起点移动到终点的成本。
	///  @param[in]		pa			前一个和当前多边形边缘的起始位置。 [（x，y，z）]
	///  @param[in]		pb			当前多边形和下一个多边形边缘上的结束位置。 [（x，y，z）]
	///  @param[in]		prevRef		前一个多边形的参考 ID。[可选]
	///  @param[in]		prevTile	包含前一个多边形的tile。[可选]
	///  @param[in]		prevPoly	前一个多边形。[可选]
	///  @param[in]		curRef		当前多边形的引用ID。
	///  @param[in]		curTile		包含当前多边形的tile。
	///  @param[in]		curPoly		当前多边形.
	///  @param[in]		nextRef		下一个多边形的参考 ID。[可选]
	///  @param[in]		nextTile	包含下一个多边形的tile。[可选]
	///  @param[in]		nextPoly	下一个多边形。[可选]
#ifdef DT_VIRTUAL_QUERYFILTER
	virtual float getCost(const float* pa, const float* pb,
						  const dtPolyRef prevRef, const dtMeshTile* prevTile, const dtPoly* prevPoly,
						  const dtPolyRef curRef, const dtMeshTile* curTile, const dtPoly* curPoly,
						  const dtPolyRef nextRef, const dtMeshTile* nextTile, const dtPoly* nextPoly) const;
#else
	float getCost(const float* pa, const float* pb,
				  const dtPolyRef prevRef, const dtMeshTile* prevTile, const dtPoly* prevPoly,
				  const dtPolyRef curRef, const dtMeshTile* curTile, const dtPoly* curPoly,
				  const dtPolyRef nextRef, const dtMeshTile* nextTile, const dtPoly* nextPoly) const;
#endif

	/// @name Getters and setters for the default implementation data.
	/// 默认实现数据的 getter 和 setter。
	///@{

	/// Returns the traversal cost of the area.
	///  @param[in]		i		The id of the area.
	/// @returns The traversal cost of the area.
	/// 返回区域的遍历成本。
	///  @param[in]		i		区域的 ID。
	/// @returns 该区域的遍历成本。
	inline float getAreaCost(const int i) const { return m_areaCost[i]; }

	/// Sets the traversal cost of the area.
	///  @param[in]		i		The id of the area.
	///  @param[in]		cost	The new cost of traversing the area.
	/// 设置区域的遍历成本。
	///  @param[in]		i		区域的 ID。
	///  @param[in]		cost	穿越该区域的新成本。
	inline void setAreaCost(const int i, const float cost) { m_areaCost[i] = cost; }

	/// Returns the include flags for the filter.
	/// Any polygons that include one or more of these flags will be
	/// included in the operation.
	/// 返回过滤器的包含标志。
	/// 包含一个或多个这些标志的任何多边形都将包含在操作中。
	inline unsigned short getIncludeFlags() const { return m_includeFlags; }

	/// Sets the include flags for the filter.
	/// @param[in]		flags	The new flags.
	/// 设置过滤器的包含标志。
	/// @param[in]		flags	新的flags.
	inline void setIncludeFlags(const unsigned short flags) { m_includeFlags = flags; }

	/// Returns the exclude flags for the filter.
	/// Any polygons that include one ore more of these flags will be
	/// excluded from the operation.
	/// 返回过滤器的排除标志。
	/// 任何包含一个或多个这些标志的多边形都将从操作中排除。
	inline unsigned short getExcludeFlags() const { return m_excludeFlags; }

	/// Sets the exclude flags for the filter.
	/// @param[in]		flags		The new flags.
	/// 设置过滤器的排除标志。
	/// @param[in]		flags		新的flags
	inline void setExcludeFlags(const unsigned short flags) { m_excludeFlags = flags; }

	///@}

};

/// Provides information about raycast hit
/// filled by dtNavMeshQuery::raycast
/// @ingroup detour
/// 提供有关由 dtNavMeshQuery::raycast 填充的光线投射命中的信息
struct dtRaycastHit
{
	/// The hit parameter. (FLT_MAX if no wall hit.)
	/// 命中参数。 （如果没有撞到墙壁，则为 FLT_MAX。）
	float t;

	/// hitNormal	The normal of the nearest wall hit. [(x, y, z)]
	/// hitNormal 最近墙壁撞击的法线。 [（x，y，z）]
	float hitNormal[3];

	/// The index of the edge on the final polygon where the wall was hit.
	/// 撞到墙的最终多边形上的边的索引。
	int hitEdgeIndex;

	/// Pointer to an array of reference ids of the visited polygons. [opt]
	/// 指向已访问多边形的引用 ID 数组的指针。 [可选]
	dtPolyRef* path;

	/// The number of visited polygons. [opt]
	/// 访问过的多边形数量。 [可选]
	int pathCount;

	/// The maximum number of polygons the @p path array can hold.
	/// @p 路径数组可以容纳的最大多边形数。
	int maxPath;

	///  The cost of the path until hit.
	///  到达之前的路径成本。
	float pathCost;
};

/// Provides custom polygon query behavior.
/// Used by dtNavMeshQuery::queryPolygons.
/// @ingroup detour
/// 提供自定义多边形查询行为。
/// 被 dtNavMeshQuery::queryPolygons. 使用
class dtPolyQuery
{
public:
	virtual ~dtPolyQuery();

	/// Called for each batch of unique polygons touched by the search area in dtNavMeshQuery::queryPolygons.
	/// This can be called multiple times for a single query.
	/// 为 dtNavMeshQuery::queryPolygons 中的搜索区域触及的每批唯一多边形调用。
	/// 对于单个查询，可以多次调用此方法。
	virtual void process(const dtMeshTile* tile, dtPoly** polys, dtPolyRef* refs, int count) = 0;
};

/// Provides the ability to perform pathfinding related queries against
/// a navigation mesh.
/// @ingroup detour
/// 提供针对导航网格执行寻路相关查询的能力。
class dtNavMeshQuery
{
public:
	dtNavMeshQuery();
	~dtNavMeshQuery();

	/// Initializes the query object.
	///  @param[in]		nav			Pointer to the dtNavMesh object to use for all queries.
	///  @param[in]		maxNodes	Maximum number of search nodes. [Limits: 0 < value <= 65535]
	/// @returns The status flags for the query.
	/// 初始化查询对象。
	///  @param[in]		nav			指向用于所有查询的 dtNavMesh 对象的指针。
	///  @param[in]		maxNodes	搜索节点的最大数量。[限制：0 < 值 <= 65535]
	/// @returns 查询的状态标志。
	dtStatus init(const dtNavMesh* nav, const int maxNodes);

	/// @name Standard Pathfinding Functions
	/// 标准寻路功能
	/// @{

	/// Finds a path from the start polygon to the end polygon.
	///  @param[in]		startRef	The reference id of the start polygon.
	///  @param[in]		endRef		The reference id of the end polygon.
	///  @param[in]		startPos	A position within the start polygon. [(x, y, z)]
	///  @param[in]		endPos		A position within the end polygon. [(x, y, z)]
	///  @param[in]		filter		The polygon filter to apply to the query.
	///  @param[out]	path		An ordered list of polygon references representing the path. (Start to end.)
	///  							[(polyRef) * @p pathCount]
	///  @param[out]	pathCount	The number of polygons returned in the @p path array.
	///  @param[in]		maxPath		The maximum number of polygons the @p path array can hold. [Limit: >= 1]
	/// 查找从起始多边形到终止多边形的路径。
	///  @param[in]		startRef	起始多边形的引用 ID。
	///  @param[in]		endRef		结束多边形的引用 ID。
	///  @param[in]		startPos	起始多边形内的位置。 [（x，y，z）]
	///  @param[in]		endPos		结束多边形内的位置。 [（x，y，z）]
	///  @param[in]		filter		要应用于查询的多边形过滤器。
	///  @param[out]	path		表示路径的多边形参考的有序列表。 (从开始到结束。) [(polyRef) * @p pathCount]
	///  @param[out]	pathCount	@p 路径数组中返回的多边形数量。
	///  @param[in]		maxPath		@p 路径数组可以容纳的最大多边形数。[限制：>= 1]
	dtStatus findPath(dtPolyRef startRef, dtPolyRef endRef,
					  const float* startPos, const float* endPos,
					  const dtQueryFilter* filter,
					  dtPolyRef* path, int* pathCount, const int maxPath) const;

	/// Finds the straight path from the start to the end position within the polygon corridor.
	///  @param[in]		startPos			Path start position. [(x, y, z)]
	///  @param[in]		endPos				Path end position. [(x, y, z)]
	///  @param[in]		path				An array of polygon references that represent the path corridor.
	///  @param[in]		pathSize			The number of polygons in the @p path array.
	///  @param[out]	straightPath		Points describing the straight path. [(x, y, z) * @p straightPathCount].
	///  @param[out]	straightPathFlags	Flags describing each point. (See: #dtStraightPathFlags) [opt]
	///  @param[out]	straightPathRefs	The reference id of the polygon that is being entered at each point. [opt]
	///  @param[out]	straightPathCount	The number of points in the straight path.
	///  @param[in]		maxStraightPath		The maximum number of points the straight path arrays can hold.  [Limit: > 0]
	///  @param[in]		options				Query options. (see: #dtStraightPathOptions)
	/// @returns The status flags for the query.
	/// 查找多边形走廊内从起点到终点位置的直线路径。
	///  @param[in]		startPos			路径起始位置。 [（x，y，z）]
	///  @param[in]		endPos				路径结束位置。 [（x，y，z）]
	///  @param[in]		path				表示路径走廊的多边形引用数组。
	///  @param[in]		pathSize			@p 路径数组中多边形的数量。
	///  @param[out]	straightPath		描述直线路径的点。 [(x, y, z) * @p StraightPathCount]。
	///  @param[out]	straightPathFlags	描述每个点的标志。 （参见：#dtStraightPathFlags）[可选]
	///  @param[out]	straightPathRefs	在每个点输入的多边形的引用 ID。 [可选]
	///  @param[out]	straightPathCount	直线路径上的点数。
	///  @param[in]		maxStraightPath		直线路径数组可以容纳的最大点数。 [限制：> 0]
	///  @param[in]		options				查询可选项. (见: #dtStraightPathOptions)
	/// @returns The status flags for the query.
	dtStatus findStraightPath(const float* startPos, const float* endPos,
							  const dtPolyRef* path, const int pathSize,
							  float* straightPath, unsigned char* straightPathFlags, dtPolyRef* straightPathRefs,
							  int* straightPathCount, const int maxStraightPath, const int options = 0) const;

	///@}
	/// @name Sliced Pathfinding Functions
	/// Common use case:
	///	-# Call initSlicedFindPath() to initialize the sliced path query.
	///	-# Call updateSlicedFindPath() until it returns complete.
	///	-# Call finalizeSlicedFindPath() to get the path.
	/// 切片寻路函数
	/// 常见用例:
	///	-# 调用 initSlicedFindPath() 初始化切片路径查询。
	///	-# 调用 updateSlicedFindPath() 直到它返回完成。
	///	-# 调用finalizeSlicedFindPath()来获取路径。
	///@{

	/// Initializes a sliced path query.
	///  @param[in]		startRef	The reference id of the start polygon.
	///  @param[in]		endRef		The reference id of the end polygon.
	///  @param[in]		startPos	A position within the start polygon. [(x, y, z)]
	///  @param[in]		endPos		A position within the end polygon. [(x, y, z)]
	///  @param[in]		filter		The polygon filter to apply to the query.
	///  @param[in]		options		query options (see: #dtFindPathOptions)
	/// @returns The status flags for the query.
	/// 初始化切片路径查询。
	///  @param[in]		startRef	起始多边形的引用 ID。
	///  @param[in]		endRef		结束多边形的引用 ID。
	///  @param[in]		startPos	起始多边形内的位置。 [（x，y，z）]
	///  @param[in]		endPos		结束多边形内的位置。 [（x，y，z）]
	///  @param[in]		filter		要应用于查询的多边形过滤器。
	///  @param[in]		options		查询选项 (见: #dtFindPathOptions)
	/// @returns 查询的状态标志。
	dtStatus initSlicedFindPath(dtPolyRef startRef, dtPolyRef endRef,
								const float* startPos, const float* endPos,
								const dtQueryFilter* filter, const unsigned int options = 0);

	/// Updates an in-progress sliced path query.
	///  @param[in]		maxIter		The maximum number of iterations to perform.
	///  @param[out]	doneIters	The actual number of iterations completed. [opt]
	/// @returns The status flags for the query.
	/// 更新正在进行的切片路径查询。
	///  @param[in]		maxIter		要执行的最大迭代次数。
	///  @param[out]	doneIters	实际完成的迭代次数。 [可选]
	/// @returns 查询的状态标志。
	dtStatus updateSlicedFindPath(const int maxIter, int* doneIters);

	/// Finalizes and returns the results of a sliced path query.
	///  @param[out]	path		An ordered list of polygon references representing the path. (Start to end.)
	///  							[(polyRef) * @p pathCount]
	///  @param[out]	pathCount	The number of polygons returned in the @p path array.
	///  @param[in]		maxPath		The max number of polygons the path array can hold. [Limit: >= 1]
	/// @returns The status flags for the query.

	/// 最终确定并返回切片路径查询的结果。
	///  @param[out]	path		表示路径的多边形引用的有序列表。 （从开始到结束。）[(polyRef) * @p pathCount]
	///  @param[out]	pathCount	@p 路径数组中返回的多边形数量。
	///  @param[in]		maxPath		路径数组可以容纳的最大多边形数。 [限制：>= 1]
	/// @returns 查询的状态标志。
	dtStatus finalizeSlicedFindPath(dtPolyRef* path, int* pathCount, const int maxPath);

	/// Finalizes and returns the results of an incomplete sliced path query, returning the path to the furthest
	/// polygon on the existing path that was visited during the search.
	///  @param[in]		existing		An array of polygon references for the existing path.
	///  @param[in]		existingSize	The number of polygon in the @p existing array.
	///  @param[out]	path			An ordered list of polygon references representing the path. (Start to end.)
	///  								[(polyRef) * @p pathCount]
	///  @param[out]	pathCount		The number of polygons returned in the @p path array.
	///  @param[in]		maxPath			The max number of polygons the @p path array can hold. [Limit: >= 1]
	/// @returns The status flags for the query.

	/// 最终确定并返回不完整切片路径查询的结果，返回搜索期间访问的现有路径上最远多边形的路径。
	///  @param[in]		existing		现有路径的多边形参考数组。
	///  @param[in]		existingSize	@p 现有数组中多边形的数量。
	///  @param[out]	path			表示路径的多边形引用的有序列表。 （从开始到结束。）[(polyRef) * @p pathCount]
	///  @param[out]	pathCount		@p 路径数组中返回的多边形数量。
	///  @param[in]		maxPath			@p 路径数组可以容纳的最大多边形数。 [限制：>= 1]
	/// @returns 查询的状态标志。
	dtStatus finalizeSlicedFindPathPartial(const dtPolyRef* existing, const int existingSize,
										   dtPolyRef* path, int* pathCount, const int maxPath);

	///@}
	/// @name Dijkstra Search Functions
	/// Dijkstra 搜索功能
	/// @{

	/// Finds the polygons along the navigation graph that touch the specified circle.
	///  @param[in]		startRef		The reference id of the polygon where the search starts.
	///  @param[in]		centerPos		The center of the search circle. [(x, y, z)]
	///  @param[in]		radius			The radius of the search circle.
	///  @param[in]		filter			The polygon filter to apply to the query.
	///  @param[out]	resultRef		The reference ids of the polygons touched by the circle. [opt]
	///  @param[out]	resultParent	The reference ids of the parent polygons for each result.
	///  								Zero if a result polygon has no parent. [opt]
	///  @param[out]	resultCost		The search cost from @p centerPos to the polygon. [opt]
	///  @param[out]	resultCount		The number of polygons found. [opt]
	///  @param[in]		maxResult		The maximum number of polygons the result arrays can hold.
	/// @returns The status flags for the query.

	/// 沿着导航图查找与指定圆接触的多边形。
	///  @param[in]		startRef		搜索开始的多边形的引用 ID。
	///  @param[in]		centerPos		搜索圈的中心。 [（x，y，z）]
	///  @param[in]		radius			搜索圆的半径。
	///  @param[in]		filter			要应用于查询的多边形过滤器。
	///  @param[out]	resultRef		圆所接触的多边形的引用 ID。 [可选]
	///  @param[out]	resultParent	每个结果的父多边形的引用 ID。
	///  								如果结果多边形没有父多边形，则为零。 [可选]
	///  @param[out]	resultCost		从@p centerPos 到多边形的搜索成本。 [可选]
	///  @param[out]	resultCount		找到的多边形数量。 [可选]
	///  @param[in]		maxResult		结果数组可以容纳的最大多边形数。
	/// @returns 查询的状态标志。
	dtStatus findPolysAroundCircle(dtPolyRef startRef, const float* centerPos, const float radius,
								   const dtQueryFilter* filter,
								   dtPolyRef* resultRef, dtPolyRef* resultParent, float* resultCost,
								   int* resultCount, const int maxResult) const;

	/// Finds the polygons along the naviation graph that touch the specified convex polygon.
	///  @param[in]		startRef		The reference id of the polygon where the search starts.
	///  @param[in]		verts			The vertices describing the convex polygon. (CCW)
	///  								[(x, y, z) * @p nverts]
	///  @param[in]		nverts			The number of vertices in the polygon.
	///  @param[in]		filter			The polygon filter to apply to the query.
	///  @param[out]	resultRef		The reference ids of the polygons touched by the search polygon. [opt]
	///  @param[out]	resultParent	The reference ids of the parent polygons for each result. Zero if a
	///  								result polygon has no parent. [opt]
	///  @param[out]	resultCost		The search cost from the centroid point to the polygon. [opt]
	///  @param[out]	resultCount		The number of polygons found.
	///  @param[in]		maxResult		The maximum number of polygons the result arrays can hold.
	/// @returns The status flags for the query.

	/// 沿着导航图查找与指定凸多边形接触的多边形。
	///  @param[in]		startRef		搜索开始的多边形的引用 ID。
	///  @param[in]		verts			描述凸多边形的顶点。 （CCW）[(x, y, z) * @p nverts]
	///  @param[in]		nverts			多边形中的顶点数。
	///  @param[in]		filter			要应用于查询的多边形过滤器。
	///  @param[out]	resultRef		搜索多边形所触及的多边形的引用 ID。 [可选]
	///  @param[out]	resultParent	每个结果的父多边形的参考 ID。 如果结果多边形没有父多边形，则为零。 [可选]
	///  @param[out]	resultCost		从质心点到多边形的搜索成本。 [可选]
	///  @param[out]	resultCount		找到的多边形数量。
	///  @param[in]		maxResult		结果数组可以容纳的最大多边形数。
	/// @returns 查询的状态标志。
	dtStatus findPolysAroundShape(dtPolyRef startRef, const float* verts, const int nverts,
								  const dtQueryFilter* filter,
								  dtPolyRef* resultRef, dtPolyRef* resultParent, float* resultCost,
								  int* resultCount, const int maxResult) const;

	/// Gets a path from the explored nodes in the previous search.
	///  @param[in]		endRef		The reference id of the end polygon.
	///  @param[out]	path		An ordered list of polygon references representing the path. (Start to end.)
	///  							[(polyRef) * @p pathCount]
	///  @param[out]	pathCount	The number of polygons returned in the @p path array.
	///  @param[in]		maxPath		The maximum number of polygons the @p path array can hold. [Limit: >= 0]
	///  @returns		The status flags. Returns DT_FAILURE | DT_INVALID_PARAM if any parameter is wrong, or if
	///  				@p endRef was not explored in the previous search. Returns DT_SUCCESS | DT_BUFFER_TOO_SMALL
	///  				if @p path cannot contain the entire path. In this case it is filled to capacity with a partial path.
	///  				Otherwise returns DT_SUCCESS.
	///  @remarks		The result of this function depends on the state of the query object. For that reason it should only
	///  				be used immediately after one of the two Dijkstra searches, findPolysAroundCircle or findPolysAroundShape.

	/// 从先前搜索中探索过的节点获取路径。
	///  @param[in]		endRef		结束多边形的引用 ID。
	///  @param[out]	path		表示路径的多边形引用的有序列表。 （从开始到结束。）[(polyRef) * @p pathCount]
	///  @param[out]	pathCount	@p 路径数组中返回的多边形数量。
	///  @param[in]		maxPath		@p 路径数组可以容纳的最大多边形数。 [限制：>= 0]
	///  @returns		状态标志。 返回 DT_FAILURE | DT_INVALID_PARAM 如果任何参数错误，或者如果在之前的搜索中未探索@p endRef。
	///                 如果 @p 路径不能包含整个路径，则返回 DT_SUCCESS|DT_BUFFER_TOO_SMALL。 在这种情况下，它已被部分路径填满。 否则返回 DT_SUCCESS。
	///  @remarks		该函数的结果取决于查询对象的状态。 因此，只能在两个 Dijkstra 搜索 findPolysAroundCircle 或 findPolysAroundShape 之一之后立即使用它。
	dtStatus getPathFromDijkstraSearch(dtPolyRef endRef, dtPolyRef* path, int* pathCount, int maxPath) const;

	/// @}
	/// @name Local Query Functions
	/// 本地查询功能
	///@{

	/// Finds the polygon nearest to the specified center point.
	/// [opt] means the specified parameter can be a null pointer, in that case the output parameter will not be set.
	///
	///  @param[in]		center		The center of the search box. [(x, y, z)]
	///  @param[in]		halfExtents	The search distance along each axis. [(x, y, z)]
	///  @param[in]		filter		The polygon filter to apply to the query.
	///  @param[out]	nearestRef	The reference id of the nearest polygon. Will be set to 0 if no polygon is found.
	///  @param[out]	nearestPt	The nearest point on the polygon. Unchanged if no polygon is found. [opt] [(x, y, z)]
	/// @returns The status flags for the query.

	/// 查找最接近指定中心点的多边形。
	/// [opt] 表示指定的参数可以是空指针，在这种情况下，不会设置输出参数。
	///
	///  @param[in]		center		搜索框的中心。 [（x，y，z）]
	///  @param[in]		halfExtents	沿每个轴的搜索距离。 [（x，y，z）]
	///  @param[in]		filter		要应用于查询的多边形过滤器。
	///  @param[out]	nearestRef	最近多边形的引用 ID。 如果未找到多边形，则设置为 0。
	///  @param[out]	nearestPt	多边形上最近的点。如果未找到多边形，则保持不变。[可选] [(x, y, z)]
	/// @returns 查询的状态标志。
	dtStatus findNearestPoly(const float* center, const float* halfExtents,
							 const dtQueryFilter* filter,
							 dtPolyRef* nearestRef, float* nearestPt) const;

	/// Finds the polygon nearest to the specified center point.
	/// [opt] means the specified parameter can be a null pointer, in that case the output parameter will not be set.
	///
	///  @param[in]		center		The center of the search box. [(x, y, z)]
	///  @param[in]		halfExtents	The search distance along each axis. [(x, y, z)]
	///  @param[in]		filter		The polygon filter to apply to the query.
	///  @param[out]	nearestRef	The reference id of the nearest polygon. Will be set to 0 if no polygon is found.
	///  @param[out]	nearestPt	The nearest point on the polygon. Unchanged if no polygon is found. [opt] [(x, y, z)]
	///  @param[out]	isOverPoly 	Set to true if the point's X/Z coordinate lies inside the polygon, false otherwise. Unchanged if no polygon is found. [opt]
	/// @returns The status flags for the query.

	/// 查找最接近指定中心点的多边形。
	/// [opt] 表示指定的参数可以是空指针，在这种情况下，不会设置输出参数。
	///
	///  @param[in]		center		搜索框的中心。 [（x，y，z）]
	///  @param[in]		halfExtents	沿每个轴的搜索距离。 [（x，y，z）]
	///  @param[in]		filter		要应用于查询的多边形过滤器。
	///  @param[out]	nearestRef	最近多边形的引用 ID。 如果未找到多边形，则设置为 0。
	///  @param[out]	nearestPt	多边形上最近的点。 如果未找到多边形，则保持不变。 [选择] [(x, y, z)]
	///  @param[out]	isOverPoly 	如果点的 X/Z 坐标位于多边形内部，则设置为 true，否则设置为 false。如果未找到多边形，则保持不变。[opt]
	/// @returns 查询的状态标志。
	dtStatus findNearestPoly(const float* center, const float* halfExtents,
							 const dtQueryFilter* filter,
							 dtPolyRef* nearestRef, float* nearestPt, bool* isOverPoly) const;

	/// Finds polygons that overlap the search box.
	///  @param[in]		center		The center of the search box. [(x, y, z)]
	///  @param[in]		halfExtents		The search distance along each axis. [(x, y, z)]
	///  @param[in]		filter		The polygon filter to apply to the query.
	///  @param[out]	polys		The reference ids of the polygons that overlap the query box.
	///  @param[out]	polyCount	The number of polygons in the search result.
	///  @param[in]		maxPolys	The maximum number of polygons the search result can hold.
	/// @returns The status flags for the query.

	/// 查找与搜索框重叠的多边形。
	///  @param[in]		center		搜索框的中心。 [（x，y，z）]
	///  @param[in]		halfExtents	沿每个轴的搜索距离。 [（x，y，z）]
	///  @param[in]		filter		要应用于查询的多边形过滤器。
	///  @param[out]	polys		与查询框重叠的多边形的引用 ID。
	///  @param[out]	polyCount	搜索结果中的多边形数量。
	///  @param[in]		maxPolys	搜索结果可以容纳的最大多边形数。
	/// @returns 查询的状态标志。
	dtStatus queryPolygons(const float* center, const float* halfExtents,
						   const dtQueryFilter* filter,
						   dtPolyRef* polys, int* polyCount, const int maxPolys) const;

	/// Finds polygons that overlap the search box.
	///  @param[in]		center		The center of the search box. [(x, y, z)]
	///  @param[in]		halfExtents		The search distance along each axis. [(x, y, z)]
	///  @param[in]		filter		The polygon filter to apply to the query.
	///  @param[in]		query		The query. Polygons found will be batched together and passed to this query.

	/// 查找与搜索框重叠的多边形。
	///  @param[in]		center		搜索框的中心。 [（x，y，z）]
	///  @param[in]		halfExtents	沿每个轴的搜索距离。 [（x，y，z）]
	///  @param[in]		filter		要应用于查询的多边形过滤器。
	///  @param[in]		query		查询。 找到的多边形将被分批在一起并传递给此查询。
	dtStatus queryPolygons(const float* center, const float* halfExtents,
						   const dtQueryFilter* filter, dtPolyQuery* query) const;

	/// Finds the non-overlapping navigation polygons in the local neighbourhood around the center position.
	///  @param[in]		startRef		The reference id of the polygon where the search starts.
	///  @param[in]		centerPos		The center of the query circle. [(x, y, z)]
	///  @param[in]		radius			The radius of the query circle.
	///  @param[in]		filter			The polygon filter to apply to the query.
	///  @param[out]	resultRef		The reference ids of the polygons touched by the circle.
	///  @param[out]	resultParent	The reference ids of the parent polygons for each result.
	///  								Zero if a result polygon has no parent. [opt]
	///  @param[out]	resultCount		The number of polygons found.
	///  @param[in]		maxResult		The maximum number of polygons the result arrays can hold.
	/// @returns The status flags for the query.

	/// 查找中心位置周围局部邻域中不重叠的导航多边形。
	///  @param[in]		startRef		搜索开始的多边形的引用 ID。
	///  @param[in]		centerPos		查询圆的中心。 [（x，y，z）]
	///  @param[in]		radius			查询圆的半径。
	///  @param[in]		filter			要应用于查询的多边形过滤器。
	///  @param[out]	resultRef		圆所接触的多边形的引用 ID。
	///  @param[out]	resultParent	每个结果的父多边形的引用 ID。如果结果多边形没有父多边形，则为零。 [选择]
	///  @param[out]	resultCount		找到的多边形数量。
	///  @param[in]		maxResult		结果数组可以容纳的最大多边形数。
	/// @returns 查询的状态标志。
	dtStatus findLocalNeighbourhood(dtPolyRef startRef, const float* centerPos, const float radius,
									const dtQueryFilter* filter,
									dtPolyRef* resultRef, dtPolyRef* resultParent,
									int* resultCount, const int maxResult) const;

	/// Moves from the start to the end position constrained to the navigation mesh.
	///  @param[in]		startRef		The reference id of the start polygon.
	///  @param[in]		startPos		A position of the mover within the start polygon. [(x, y, x)]
	///  @param[in]		endPos			The desired end position of the mover. [(x, y, z)]
	///  @param[in]		filter			The polygon filter to apply to the query.
	///  @param[out]	resultPos		The result position of the mover. [(x, y, z)]
	///  @param[out]	visited			The reference ids of the polygons visited during the move.
	///  @param[out]	visitedCount	The number of polygons visited during the move.
	///  @param[in]		maxVisitedSize	The maximum number of polygons the @p visited array can hold.
	/// @returns The status flags for the query.

	/// 从开始位置移动到受限于导航网格的结束位置。
	///  @param[in]		startRef		起始多边形的引用 ID。
	///  @param[in]		startPos		移动器在起始多边形内的位置。 [（x，y，x）]
	///  @param[in]		endPos			所需的动子最终位置。 [（x，y，z）]
	///  @param[in]		filter			要应用于查询的多边形过滤器。
	///  @param[out]	resultPos		动子的结果位置。 [（x，y，z）]
	///  @param[out]	visited			移动期间访问的多边形的引用 ID。
	///  @param[out]	visitedCount	移动过程中访问的多边形数量。
	///  @param[in]		maxVisitedSize	@p 访问过的数组可以容纳的最大多边形数。
	/// @returns 查询的状态标志。
	dtStatus moveAlongSurface(dtPolyRef startRef, const float* startPos, const float* endPos,
							  const dtQueryFilter* filter,
							  float* resultPos, dtPolyRef* visited, int* visitedCount, const int maxVisitedSize) const;

	/// Casts a 'walkability' ray along the surface of the navigation mesh from
	/// the start position toward the end position.
	/// @note A wrapper around raycast(..., RaycastHit*). Retained for backward compatibility.
	///  @param[in]		startRef	The reference id of the start polygon.
	///  @param[in]		startPos	A position within the start polygon representing
	///  							the start of the ray. [(x, y, z)]
	///  @param[in]		endPos		The position to cast the ray toward. [(x, y, z)]
	///  @param[in]		filter		The polygon filter to apply to the query.
	///  @param[out]	t			The hit parameter. (FLT_MAX if no wall hit.)
	///  @param[out]	hitNormal	The normal of the nearest wall hit. [(x, y, z)]
	///  @param[out]	path		The reference ids of the visited polygons. [opt]
	///  @param[out]	pathCount	The number of visited polygons. [opt]
	///  @param[in]		maxPath		The maximum number of polygons the @p path array can hold.
	/// @returns The status flags for the query.

	/// 沿着导航网格的表面从起始位置向结束位置投射一条“walkability”光线。
	/// @note raycast(..., RaycastHit*) 的包装。 保留是为了向后兼容。
	///  @param[in]		startRef	起始多边形的引用 ID。
	///  @param[in]		startPos	起始多边形内表示射线起点的位置。 [（x，y，z）]
	///  @param[in]		endPos		光线投射的位置。 [（x，y，z）]
	///  @param[in]		filter		要应用于查询的多边形过滤器。
	///  @param[out]	t			命中参数。 （如果没有撞到墙壁，则为 FLT_MAX。）
	///  @param[out]	hitNormal	最近的墙壁撞击的法线。 [（x，y，z）]
	///  @param[out]	path		访问过的多边形的引用 ID。 [opt]
	///  @param[out]	pathCount	访问过的多边形数量。[opt]
	///  @param[in]		maxPath		@p 路径数组可以容纳的最大多边形数。
	/// @returns 查询的状态标志。
	dtStatus raycast(dtPolyRef startRef, const float* startPos, const float* endPos,
					 const dtQueryFilter* filter,
					 float* t, float* hitNormal, dtPolyRef* path, int* pathCount, const int maxPath) const;

	/// Casts a 'walkability' ray along the surface of the navigation mesh from
	/// the start position toward the end position.
	///  @param[in]		startRef	The reference id of the start polygon.
	///  @param[in]		startPos	A position within the start polygon representing
	///  							the start of the ray. [(x, y, z)]
	///  @param[in]		endPos		The position to cast the ray toward. [(x, y, z)]
	///  @param[in]		filter		The polygon filter to apply to the query.
	///  @param[in]		options		govern how the raycast behaves. See dtRaycastOptions
	///  @param[out]	hit			Pointer to a raycast hit structure which will be filled by the results.
	///  @param[in]		prevRef		parent of start ref. Used during for cost calculation [opt]
	/// @returns The status flags for the query.

	/// 沿着导航网格的表面从起始位置向结束位置投射一条“walkability”光线。
	///  @param[in]		startRef	起始多边形的引用 ID。
	///  @param[in]		startPos	起始多边形内表示射线起点的位置。 [（x，y，z）]
	///  @param[in]		endPos		光线投射的位置。 [（x，y，z）]
	///  @param[in]		filter		要应用于查询的多边形过滤器。
	///  @param[in]		options		控制光线投射的行为方式。 请参阅 dtRaycastOptions
	///  @param[out]	hit			指向将由结果填充的光线投射命中结构的指针。
	///  @param[in]		prevRef		起始引用的父级。 用于成本计算 [opt]
	/// @returns 查询的状态标志。
	dtStatus raycast(dtPolyRef startRef, const float* startPos, const float* endPos,
					 const dtQueryFilter* filter, const unsigned int options,
					 dtRaycastHit* hit, dtPolyRef prevRef = 0) const;


	/// Finds the distance from the specified position to the nearest polygon wall.
	///  @param[in]		startRef		The reference id of the polygon containing @p centerPos.
	///  @param[in]		centerPos		The center of the search circle. [(x, y, z)]
	///  @param[in]		maxRadius		The radius of the search circle.
	///  @param[in]		filter			The polygon filter to apply to the query.
	///  @param[out]	hitDist			The distance to the nearest wall from @p centerPos.
	///  @param[out]	hitPos			The nearest position on the wall that was hit. [(x, y, z)]
	///  @param[out]	hitNormal		The normalized ray formed from the wall point to the
	///  								source point. [(x, y, z)]
	/// @returns The status flags for the query.

	/// 查找从指定位置到最近的多边形墙的距离。
	///  @param[in]		startRef		包含@p centerPos的多边形的参考id。
	///  @param[in]		centerPos		搜索圈的中心。 [（x，y，z）]
	///  @param[in]		maxRadius		搜索圆的半径。
	///  @param[in]		filter			要应用于查询的多边形过滤器。
	///  @param[out]	hitDist			从 @p centerPos 到最近的墙壁的距离。
	///  @param[out]	hitPos			墙上距离被击中最近的位置。 [（x，y，z）]
	///  @param[out]	hitNormal		从壁点到源点形成的归一化射线。 [（x，y，z）]
	/// @returns 查询的状态标志。
	dtStatus findDistanceToWall(dtPolyRef startRef, const float* centerPos, const float maxRadius,
								const dtQueryFilter* filter,
								float* hitDist, float* hitPos, float* hitNormal) const;

	/// Returns the segments for the specified polygon, optionally including portals.
	///  @param[in]		ref				The reference id of the polygon.
	///  @param[in]		filter			The polygon filter to apply to the query.
	///  @param[out]	segmentVerts	The segments. [(ax, ay, az, bx, by, bz) * segmentCount]
	///  @param[out]	segmentRefs		The reference ids of each segment's neighbor polygon.
	///  								Or zero if the segment is a wall. [opt] [(parentRef) * @p segmentCount]
	///  @param[out]	segmentCount	The number of segments returned.
	///  @param[in]		maxSegments		The maximum number of segments the result arrays can hold.
	/// @returns The status flags for the query.

	/// 返回指定多边形的线段，可选包含portals。
	///  @param[in]		ref				多边形的引用 ID。
	///  @param[in]		filter			要应用于查询的多边形过滤器。
	///  @param[out]	segmentVerts	片段。 [(ax, ay, az, bx, by, bz) * segmentCount]
	///  @param[out]	segmentRefs		每个线段的相邻多边形的引用ID。如果该段是墙，则为零。[opt] [(parentRef) * @p segmentCount]
	///  @param[out]	segmentCount	返回的段数。
	///  @param[in]		maxSegments		结果数组可以容纳的最大段数。
	/// @returns 查询的状态标志。
	dtStatus getPolyWallSegments(dtPolyRef ref, const dtQueryFilter* filter,
								 float* segmentVerts, dtPolyRef* segmentRefs, int* segmentCount,
								 const int maxSegments) const;

	/// Returns random location on navmesh.
	/// Polygons are chosen weighted by area. The search runs in linear related to number of polygon.
	///  @param[in]		filter			The polygon filter to apply to the query.
	///  @param[in]		frand			Function returning a random number [0..1).
	///  @param[out]	randomRef		The reference id of the random location.
	///  @param[out]	randomPt		The random location.
	/// @returns The status flags for the query.

	/// 返回导航网格上的随机位置。
	/// 多边形是按面积加权选择的。 搜索以与多边形数量线性相关的方式运行。
	///  @param[in]		filter			要应用于查询的多边形过滤器。
	///  @param[in]		frand			函数返回一个随机数 [0..1)。
	///  @param[out]	randomRef		随机位置的引用 ID。
	///  @param[out]	randomPt		随机位置。
	/// @returns 查询的状态标志。
	dtStatus findRandomPoint(const dtQueryFilter* filter, float (*frand)(),
							 dtPolyRef* randomRef, float* randomPt) const;

	/// Returns random location on navmesh within the reach of specified location.
	/// Polygons are chosen weighted by area. The search runs in linear related to number of polygon.
	/// The location is not exactly constrained by the circle, but it limits the visited polygons.
	///  @param[in]		startRef		The reference id of the polygon where the search starts.
	///  @param[in]		centerPos		The center of the search circle. [(x, y, z)]
	///  @param[in]		maxRadius		The radius of the search circle. [Units: wu]
	///  @param[in]		filter			The polygon filter to apply to the query.
	///  @param[in]		frand			Function returning a random number [0..1).
	///  @param[out]	randomRef		The reference id of the random location.
	///  @param[out]	randomPt		The random location. [(x, y, z)]
	/// @returns The status flags for the query.

	/// 返回导航网格上指定位置范围内的随机位置。
	/// 多边形是按面积加权选择的。 搜索以与多边形数量线性相关的方式运行。
	/// 该位置并不完全受圆约束，但它限制了访问的多边形。
	///  @param[in]		startRef		搜索开始的多边形的引用 ID。
	///  @param[in]		centerPos		搜索圈的中心。 [（x，y，z）]
	///  @param[in]		maxRadius		搜索圆的半径。 [单位：wu]
	///  @param[in]		filter			要应用于查询的多边形过滤器。
	///  @param[in]		frand			函数返回一个随机数 [0..1)。
	///  @param[out]	randomRef		随机位置的引用 ID。
	///  @param[out]	randomPt		随机位置。 [（x，y，z）]
	/// @returns 查询的状态标志。
	dtStatus findRandomPointAroundCircle(dtPolyRef startRef, const float* centerPos, const float maxRadius,
										 const dtQueryFilter* filter, float (*frand)(),
										 dtPolyRef* randomRef, float* randomPt) const;

	/// Finds the closest point on the specified polygon.
	///  @param[in]		ref			The reference id of the polygon.
	///  @param[in]		pos			The position to check. [(x, y, z)]
	///  @param[out]	closest		The closest point on the polygon. [(x, y, z)]
	///  @param[out]	posOverPoly	True of the position is over the polygon.
	/// @returns The status flags for the query.

	/// 查找指定多边形上最近的点。
	///  @param[in]		ref			多边形的参考 ID。
	///  @param[in]		pos			要检查的位置。 [（x，y，z）]
	///  @param[out]	closest		多边形上最近的点。 [（x，y，z）]
	///  @param[out]	posOverPoly	正确的位置是在多边形上方。
	/// @returns 查询的状态标志。
	dtStatus closestPointOnPoly(dtPolyRef ref, const float* pos, float* closest, bool* posOverPoly) const;

	/// Returns a point on the boundary closest to the source point if the source point is outside the
	/// polygon's xz-bounds.
	///  @param[in]		ref			The reference id to the polygon.
	///  @param[in]		pos			The position to check. [(x, y, z)]
	///  @param[out]	closest		The closest point. [(x, y, z)]
	/// @returns The status flags for the query.

	/// 如果源点位于多边形的 xz 边界之外，则返回边界上最接近源点的点。
	///  @param[in]		ref			多边形的引用 ID。
	///  @param[in]		pos			要检查的位置。 [（x，y，z）]
	///  @param[out]	closest		最近的点。 [（x，y，z）]
	/// @returns 查询的状态标志。
	dtStatus closestPointOnPolyBoundary(dtPolyRef ref, const float* pos, float* closest) const;

	/// Gets the height of the polygon at the provided position using the height detail. (Most accurate.)
	///  @param[in]		ref			The reference id of the polygon.
	///  @param[in]		pos			A position within the xz-bounds of the polygon. [(x, y, z)]
	///  @param[out]	height		The height at the surface of the polygon.
	/// @returns The status flags for the query.

	/// 使用高度详细信息获取指定位置处多边形的高度。 （最准确。）
	///  @param[in]		ref			多边形的引用 ID。
	///  @param[in]		pos			多边形 xz 边界内的位置。 [（x，y，z）]
	///  @param[out]	height		多边形表面的高度。
	/// @returns 查询的状态标志。
	dtStatus getPolyHeight(dtPolyRef ref, const float* pos, float* height) const;

	/// @}
	/// @name Miscellaneous Functions
	/// 杂项功能
	/// @{

	/// Returns true if the polygon reference is valid and passes the filter restrictions.
	///  @param[in]		ref			The polygon reference to check.
	///  @param[in]		filter		The filter to apply.
	/// 如果多边形引用有效并通过过滤器限制，则返回 true。
	///  @param[in]		ref			要检查的多边形引用。
	///  @param[in]		filter		要应用的过滤器。
	bool isValidPolyRef(dtPolyRef ref, const dtQueryFilter* filter) const;

	/// Returns true if the polygon reference is in the closed list.
	///  @param[in]		ref		The reference id of the polygon to check.
	/// @returns True if the polygon is in closed list.
	/// 如果多边形引用位于关闭列表中，则返回 true。
	///  @param[in]		ref		要检查的多边形的引用 ID。
	/// @returns 如果多边形位于封闭列表中，则为 True。
	bool isInClosedList(dtPolyRef ref) const;

	/// Gets the node pool.
	/// @returns The node pool.
	/// 获取节点池。
	/// @returns 节点池。
	class dtNodePool* getNodePool() const { return m_nodePool; }

	/// Gets the navigation mesh the query object is using.
	/// @return The navigation mesh the query object is using.
	/// 获取查询对象正在使用的导航网格。
	/// @return 查询对象正在使用的导航网格。
	const dtNavMesh* getAttachedNavMesh() const { return m_nav; }

	/// @}

  private:

	// Explicitly disabled copy constructor and copy assignment operator
	// 显式禁用复制构造函数和复制赋值运算符
	dtNavMeshQuery(const dtNavMeshQuery&);
	dtNavMeshQuery& operator=(const dtNavMeshQuery&);

	/// Queries polygons within a tile.
	/// 查询图块内的多边形。
	void queryPolygonsInTile(const dtMeshTile* tile, const float* qmin, const float* qmax,
							 const dtQueryFilter* filter, dtPolyQuery* query) const;

	/// Returns portal points between two polygons.
	/// 返回两个多边形之间的入口点。
	dtStatus getPortalPoints(dtPolyRef from, dtPolyRef to, float* left, float* right,
							 unsigned char& fromType, unsigned char& toType) const;
	dtStatus getPortalPoints(dtPolyRef from, const dtPoly* fromPoly, const dtMeshTile* fromTile,
							 dtPolyRef to, const dtPoly* toPoly, const dtMeshTile* toTile,
							 float* left, float* right) const;

	/// Returns edge mid point between two polygons.
	/// 返回两个多边形之间的边中点。
	dtStatus getEdgeMidPoint(dtPolyRef from, dtPolyRef to, float* mid) const;
	dtStatus getEdgeMidPoint(dtPolyRef from, const dtPoly* fromPoly, const dtMeshTile* fromTile,
							 dtPolyRef to, const dtPoly* toPoly, const dtMeshTile* toTile,
							 float* mid) const;

	// Appends vertex to a straight path
	// 将顶点附加到直线路径
	dtStatus appendVertex(const float* pos, const unsigned char flags, const dtPolyRef ref,
						  float* straightPath, unsigned char* straightPathFlags, dtPolyRef* straightPathRefs,
						  int* straightPathCount, const int maxStraightPath) const;

	// Appends intermediate portal points to a straight path.
	// 将中间入口点附加到直线路径。
	dtStatus appendPortals(const int startIdx, const int endIdx, const float* endPos, const dtPolyRef* path,
						   float* straightPath, unsigned char* straightPathFlags, dtPolyRef* straightPathRefs,
						   int* straightPathCount, const int maxStraightPath, const int options) const;

	// Gets the path leading to the specified end node.
	// 获取通往指定结束节点的路径。
	dtStatus getPathToNode(struct dtNode* endNode, dtPolyRef* path, int* pathCount, int maxPath) const;

	///< Pointer to navmesh data.
	///< 指向导航网格数据的指针。
	const dtNavMesh* m_nav;

	struct dtQueryData
	{
		dtStatus status;
		struct dtNode* lastBestNode;
		float lastBestNodeCost;
		dtPolyRef startRef, endRef;
		float startPos[3], endPos[3];
		const dtQueryFilter* filter;
		unsigned int options;
		float raycastLimitSqr;
	};
	///< 切片查询状态。
	dtQueryData m_query;

	///< Pointer to small node pool.
	///< 指向小节点池的指针。
	class dtNodePool* m_tinyNodePool;
	///< Pointer to node pool.
	///< 指向节点池的指针。
	class dtNodePool* m_nodePool;
	///< Pointer to open list queue.
	///< 指向打开列表队列的指针。
	class dtNodeQueue* m_openList;
};

/// Allocates a query object using the Detour allocator.
/// @return An allocated query object, or null on failure.
/// @ingroup detour
/// 使用 Detour 分配器分配查询对象。
/// @return 分配的查询对象，或失败时为 null。
dtNavMeshQuery* dtAllocNavMeshQuery();

/// Frees the specified query object using the Detour allocator.
///  @param[in]		query		A query object allocated using #dtAllocNavMeshQuery
/// @ingroup detour
/// 使用 Detour 分配器释放指定的查询对象。
///  @param[in]		query		使用 #dtAllocNavMeshQuery 分配的查询对象
void dtFreeNavMeshQuery(dtNavMeshQuery* query);

#endif // DETOURNAVMESHQUERY_H
