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

#ifndef RECAST_H
#define RECAST_H

/// The value of PI used by Recast.
static const float RC_PI = 3.14159265f;

/// Used to ignore unused function parameters and silence any compiler warnings.
/// 用于消除未使用的函数参数产生的编译器警告。
template<class T> void rcIgnoreUnused(const T&) { }

/// Recast log categories.
/// @see rcContext
enum rcLogCategory
{
	RC_LOG_PROGRESS = 1,	///< A progress log entry.
	RC_LOG_WARNING,			///< A warning log entry.
	RC_LOG_ERROR			///< An error log entry.
};

/// Recast performance timer categories.
/// @see rcContext
enum rcTimerLabel
{
	/// The user defined total time of the build.
	/// 用户定义的构建总时间。
	RC_TIMER_TOTAL,
	/// A user defined build time.
	/// 用户定义的构建时间。
	RC_TIMER_TEMP,
	/// The time to rasterize the triangles. (See: #rcRasterizeTriangle)
	/// 光栅化三角形的时间。 （参见：#rcRasterizeTriangle）
	RC_TIMER_RASTERIZE_TRIANGLES,
	/// The time to build the compact heightfield. (See: #rcBuildCompactHeightfield)
	/// 构建紧凑高度场的时间。 （参见：#rcBuildCompactHeightfield）
	RC_TIMER_BUILD_COMPACTHEIGHTFIELD,
	/// The total time to build the contours(轮廓). (See: #rcBuildContours)
	/// 构建轮廓的总时间。 （参见：#rcBuildContours）
	RC_TIMER_BUILD_CONTOURS,
	/// The time to trace the boundaries of the contours. (See: #rcBuildContours)
	/// 追踪轮廓边界的时间。 （参见：#rcBuildContours）
	RC_TIMER_BUILD_CONTOURS_TRACE,
	/// The time to simplify the contours. (See: #rcBuildContours)
	/// 是时候简化轮廓了。 （参见：#rcBuildContours）
	RC_TIMER_BUILD_CONTOURS_SIMPLIFY,
	/// The time to filter ledge spans. (See: #rcFilterLedgeSpans)
	/// 过滤壁架跨度的时间。 （参见：#rcFilterLedgeSpans）
	RC_TIMER_FILTER_BORDER,
	/// The time to filter low height spans. (See: #rcFilterWalkableLowHeightSpans)
	/// 过滤低高度跨度的时间。 （参见：#rcFilterWalkableLowHeightSpans）
	RC_TIMER_FILTER_WALKABLE,
	/// The time to apply the median filter. (See: #rcMedianFilterWalkableArea)
	/// 应用中值过滤器的时间。 （参见：#rcMedianFilterWalkableArea）
	RC_TIMER_MEDIAN_AREA,
	/// The time to filter low obstacles. (See: #rcFilterLowHangingWalkableObstacles)
	/// 过滤低障碍物的时间。 （参见：#rcFilterLowHangingWalkableObstacles）
	RC_TIMER_FILTER_LOW_OBSTACLES,
	/// The time to build the polygon mesh. (See: #rcBuildPolyMesh)
	/// 构建多边形网格的时间。 （参见：#rcBuildPolyMesh）
	RC_TIMER_BUILD_POLYMESH,
	/// The time to merge polygon meshes. (See: #rcMergePolyMeshes)
	/// 合并多边形网格的时间。 （参见：#rcMergePolyMeshes）
	RC_TIMER_MERGE_POLYMESH,
	/// The time to erode the walkable area. (See: #rcErodeWalkableArea)
	/// 裁剪可行走区域时间。 （参见：#rcErodeWalkableArea）
	RC_TIMER_ERODE_AREA,
	/// The time to mark a box area. (See: #rcMarkBoxArea)
	/// 标记体素掩码区域的时间。 （参见：#rcMarkBoxArea）
	RC_TIMER_MARK_BOX_AREA,
	/// The time to mark a cylinder area. (See: #rcMarkCylinderArea)
	/// 标记圆柱区域的时间。 （参见：#rcMarkCylinderArea）
	RC_TIMER_MARK_CYLINDER_AREA,
	/// 标记凸多边形区域的时间。 （参见：#rcMarkConvexPolyArea）
	RC_TIMER_MARK_CONVEXPOLY_AREA,
	/// The total time to build the distance field. (See: #rcBuildDistanceField)
	/// 构建距离场的总时间。 （参见：#rcBuildDistanceField）
	RC_TIMER_BUILD_DISTANCEFIELD,
	/// The time to build the distances of the distance field. (See: #rcBuildDistanceField)
	/// 建立距离场距离的时间。 （参见：#rcBuildDistanceField）
	RC_TIMER_BUILD_DISTANCEFIELD_DIST,
	/// The time to blur the distance field. (See: #rcBuildDistanceField)
	/// 模糊距离场的时间。 （参见：#rcBuildDistanceField）
	RC_TIMER_BUILD_DISTANCEFIELD_BLUR,
	/// The total time to build the regions. (See: #rcBuildRegions, #rcBuildRegionsMonotone)
	/// 构建区域的总时间。 （参见：#rcBuildRegions、#rcBuildRegionsMonotone）
	RC_TIMER_BUILD_REGIONS,
	/// The total time to apply the watershed algorithm. (See: #rcBuildRegions)
	/// 应用分水岭算法的总时间。 （参见：#rcBuildRegions）
	RC_TIMER_BUILD_REGIONS_WATERSHED,
	/// The time to expand regions while applying the watershed algorithm. (See: #rcBuildRegions)
	/// 应用分水岭算法时扩展区域的时间。 （参见：#rcBuildRegions）
	RC_TIMER_BUILD_REGIONS_EXPAND,
	/// The time to flood regions while applying the watershed algorithm. (See: #rcBuildRegions)
	/// 应用分水岭算法时淹没区域的时间。 （参见：#rcBuildRegions）
	RC_TIMER_BUILD_REGIONS_FLOOD,
	/// The time to filter out small regions. (See: #rcBuildRegions, #rcBuildRegionsMonotone)
	/// 是时候过滤掉小区域了。 （参见：#rcBuildRegions、#rcBuildRegionsMonotone）
	RC_TIMER_BUILD_REGIONS_FILTER,
	/// The time to build heightfield layers. (See: #rcBuildHeightfieldLayers)
	/// 构建高度场图层的时间。 （参见：#rcBuildHeightfieldLayers）
	RC_TIMER_BUILD_LAYERS,
	/// The time to build the polygon mesh detail. (See: #rcBuildPolyMeshDetail)
	/// 构建多边形网格细节的时间。 （参见：#rcBuildPolyMeshDetail）
	RC_TIMER_BUILD_POLYMESHDETAIL,
	/// The time to merge polygon mesh details. (See: #rcMergePolyMeshDetails)
	/// 合并多边形网格细节的时间。 （参见：#rcMergePolyMeshDetails）
	RC_TIMER_MERGE_POLYMESHDETAIL,
	/// The maximum number of timers.  (Used for iterating timers.)
	RC_MAX_TIMERS
};

/// Provides an interface for optional logging and performance tracking of the Recast
/// build process.
///
/// This class does not provide logging or timer functionality on its
/// own.  Both must be provided by a concrete implementation
/// by overriding the protected member functions.  Also, this class does not
/// provide an interface for extracting log messages. (Only adding them.)
/// So concrete implementations must provide one.
///
/// If no logging or timers are required, just pass an instance of this
/// class through the Recast build process.
///
/// @ingroup recast
///
/// 提供用于 Recast 构建过程的可选日志记录和性能跟踪的接口。
///
/// 此类本身不提供日志记录或计时器功能。 两者都必须通过覆盖受保护的成员函数的具体实现来提供。
/// 此外，此类不提供用于提取日志消息的接口。 （仅添加它们。）因此具体实现必须提供一个。
///
/// 如果不需要日志记录或计时器，只需通过 Recast 构建过程传递此类的实例即可。
///
class rcContext
{
public:
	/// Constructor.
	///  @param[in]		state	TRUE if the logging and performance timers should be enabled.  [Default: true]
	inline rcContext(bool state = true) : m_logEnabled(state), m_timerEnabled(state) {}
	virtual ~rcContext() {}

	/// Enables or disables logging.
	///  @param[in]		state	TRUE if logging should be enabled.
	inline void enableLog(bool state) { m_logEnabled = state; }

	/// Clears all log entries.
	inline void resetLog() { if (m_logEnabled) doResetLog(); }

	/// Logs a message.
	///
	/// Example:
	/// @code
	/// // Where ctx is an instance of rcContext and filepath is a char array.
	/// ctx->log(RC_LOG_ERROR, "buildTiledNavigation: Could not load '%s'", filepath);
	/// @endcode
	///
	/// @param[in]		category	The category of the message.
	/// @param[in]		format		The message.
	void log(const rcLogCategory category, const char* format, ...);

	/// Enables or disables the performance timers.
	///  @param[in]		state	TRUE if timers should be enabled.
	inline void enableTimer(bool state) { m_timerEnabled = state; }

	/// Clears all performance timers. (Resets all to unused.)
	inline void resetTimers() { if (m_timerEnabled) doResetTimers(); }

	/// Starts the specified performance timer.
	/// @param	label	The category of the timer.
	inline void startTimer(const rcTimerLabel label) { if (m_timerEnabled) doStartTimer(label); }

	/// Stops the specified performance timer.
	/// @param	label	The category of the timer.
	inline void stopTimer(const rcTimerLabel label) { if (m_timerEnabled) doStopTimer(label); }

	/// Returns the total accumulated time of the specified performance timer.
	/// @param	label	The category of the timer.
	/// @return The accumulated time of the timer, or -1 if timers are disabled or the timer has never been started.
	/// @return 定时器的累计时间，如果定时器被禁用或定时器从未启动过，则为 -1
	inline int getAccumulatedTime(const rcTimerLabel label) const { return m_timerEnabled ? doGetAccumulatedTime(label) : -1; }

protected:
	/// Clears all log entries.
	virtual void doResetLog();

	/// Logs a message.
	/// @param[in]		category	The category of the message.
	/// @param[in]		msg			The formatted message.
	/// @param[in]		len			The length of the formatted message.
	virtual void doLog(const rcLogCategory category, const char* msg, const int len) { rcIgnoreUnused(category); rcIgnoreUnused(msg); rcIgnoreUnused(len); }

	/// Clears all timers. (Resets all to unused.)
	virtual void doResetTimers() {}

	/// Starts the specified performance timer.
	/// @param[in]		label	The category of timer.
	virtual void doStartTimer(const rcTimerLabel label) { rcIgnoreUnused(label); }

	/// Stops the specified performance timer.
	/// @param[in]		label	The category of the timer.
	virtual void doStopTimer(const rcTimerLabel label) { rcIgnoreUnused(label); }

	/// Returns the total accumulated time of the specified performance timer.
	/// @param[in]		label	The category of the timer.
	/// @return The accumulated time of the timer, or -1 if timers are disabled or the timer has never been started.
	virtual int doGetAccumulatedTime(const rcTimerLabel label) const { rcIgnoreUnused(label); return -1; }

	/// True if logging is enabled.
	bool m_logEnabled;

	/// True if the performance timers are enabled.
	bool m_timerEnabled;
};

/// A helper to first start a timer and then stop it when this helper goes out of scope.
/// @see rcContext
/// 一个助手首先启动一个计时器，然后在该助手超出范围时停止它。
class rcScopedTimer
{
  public:
	/// Constructs an instance and starts the timer.
	///  @param[in]		ctx		The context to use.
	///  @param[in]		label	The category of the timer.
	inline rcScopedTimer(rcContext* ctx, const rcTimerLabel label) : m_ctx(ctx), m_label(label) { m_ctx->startTimer(m_label); }
	inline ~rcScopedTimer() { m_ctx->stopTimer(m_label); }

  private:
	// Explicitly disabled copy constructor and copy assignment operator.
	rcScopedTimer(const rcScopedTimer&);
	rcScopedTimer& operator=(const rcScopedTimer&);

	rcContext* const m_ctx;
	const rcTimerLabel m_label;
};

/// Specifies a configuration to use when performing Recast builds.
/// 指定执行 Recast 构建时要使用的配置。
/// @ingroup recast
struct rcConfig
{
	/// The width of the field along the x-axis. [Limit: >= 0] [Units: vx]
	/// 沿 x 轴的宽度。 [限制：>= 0] [单位：vx]
	int width;

	/// The height of the field along the z-axis. [Limit: >= 0] [Units: vx]
	/// 沿 z 轴的高度。 [限制：>= 0] [单位：vx]
	int height;

	/// The width/height size of tile's on the xz-plane. [Limit: >= 0] [Units: vx]
	/// xz 平面上tile的宽度/高度尺寸。 [限制：>= 0] [单位：vx]
	int tileSize;

	/// The size of the non-navigable border around the heightfield. [Limit: >=0] [Units: vx]
	/// 高度场周围不可导航边框的大小。 [限制：>=0] [单位：vx]
	int borderSize;

	/// The xz-plane cell size to use for fields. [Limit: > 0] [Units: wu]
	/// xz 平面的cell size。 [限制：> 0] [单位：wu]
	float cs;

	/// The y-axis cell size to use for fields. [Limit: > 0] [Units: wu]
	/// y轴的 cell size . [Limit: > 0] [Units: wu]
	float ch;

	/// The minimum bounds of the field's AABB. [(x, y, z)] [Units: wu]
	/// AABB盒 的最小界限。 [(x, y, z)] [单位：wu]
	float bmin[3];

	/// The maximum bounds of the field's AABB. [(x, y, z)] [Units: wu]
	/// AABB盒 的最大界限。 [(x, y, z)] [单位：wu]
	float bmax[3];

	/// The maximum slope that is considered walkable. [Limits: 0 <= value < 90] [Units: Degrees]
	/// 可行走平面的最大坡度。 [限制：0 <= 值 < 90] [单位：Degrees]
	float walkableSlopeAngle;

	/// Minimum floor to 'ceiling' height that will still allow the floor area to
	/// be considered walkable. [Limit: >= 3] [Units: vx]
	/// 可行走的地板到“天花板”的最小高度。 [限制：>= 3] [单位：vx]
	int walkableHeight;

	/// Maximum ledge height that is considered to still be traversable. [Limit: >=0] [Units: vx]
	/// 可行走的相邻地区平面最大落差高度。 [限制：>=0] [单位：vx]
	int walkableClimb;

	/// The distance to erode/shrink the walkable area of the heightfield away from
	/// obstructions.  [Limit: >=0] [Units: vx]
	/// 可行走半径，缩小可行走区域轮廓。 [限制：>=0] [单位：vx]
	int walkableRadius;

	/// The maximum allowed length for contour edges along the border of the mesh. [Limit: >=0] [Units: vx]
	/// 沿网格边界的轮廓边缘的最大允许长度。 [限制：>=0] [单位：vx]
	int maxEdgeLen;

	/// The maximum distance a simplified contour's border edges should deviate
	/// the original raw contour. [Limit: >=0] [Units: vx]
	/// 简化轮廓的边界边缘应偏离原始原始轮廓的最大距离。 [限制：>=0] [单位：vx]
	float maxSimplificationError;

	/// The minimum number of cells allowed to form isolated island areas. [Limit: >=0] [Units: vx]
	/// 允许形成孤岛区域的最小单元数。 [限制：>=0] [单位：vx]
	int minRegionArea;

	/// Any regions with a span count smaller than this value will, if possible,
	/// be merged with larger regions. [Limit: >=0] [Units: vx]
	/// 如果可能，任何span数量小于此值的区域都将与更大的区域合并。 [限制：>=0] [单位：vx]
	int mergeRegionArea;

	/// The maximum number of vertices allowed for polygons generated during the
	/// contour to polygon conversion process. [Limit: >= 3]
	/// 轮廓到多边形转换过程中生成的多边形允许的最大顶点数。 [限制：>= 3]
	int maxVertsPerPoly;

	/// Sets the sampling distance to use when generating the detail mesh.
	/// (For height detail only.) [Limits: 0 or >= 0.9] [Units: wu]
	/// 设置生成细节网格时使用的采样距离。 （仅适用于高度详细信息。） [限制：0 或 >= 0.9] [单位：wu]
	float detailSampleDist;

	/// The maximum distance the detail mesh surface should deviate from heightfield
	/// data. (For height detail only.) [Limit: >=0] [Units: wu]
	/// 细节网格表面应偏离高度场数据的最大距离。 （仅适用于高度详细信息。） [限制：>=0] [单位：wu]
	float detailSampleMaxError;
};

/// Defines the number of bits allocated to rcSpan::smin and rcSpan::smax.
/// 定义分配给 rcSpan::smin 和 rcSpan::smax 的二进制位数。
static const int RC_SPAN_HEIGHT_BITS = 13;
/// Defines the maximum value for rcSpan::smin and rcSpan::smax.
/// 定义 rcSpan::smin 和 rcSpan::smax 的最大值
static const int RC_SPAN_MAX_HEIGHT = (1 << RC_SPAN_HEIGHT_BITS) - 1;

/// The number of spans allocated per span spool.
/// @see rcSpanPool
/// 每个span池分配的span数。
static const int RC_SPANS_PER_POOL = 2048;

/// Represents a span in a heightfield.
/// @see rcHeightfield
/// 表示高度场中的跨度。
struct rcSpan
{
	unsigned int smin : RC_SPAN_HEIGHT_BITS; ///< The lower limit of the span. [Limit: < #smax]
	unsigned int smax : RC_SPAN_HEIGHT_BITS; ///< The upper limit of the span. [Limit: <= #RC_SPAN_MAX_HEIGHT]
	unsigned int area : 6;                   ///< The area id assigned to the span.(分配给跨度的区域 ID)
	rcSpan* next;                            ///< The next span higher up in column.(列中较高位置的下一个跨度)
};

/// A memory pool used for quick allocation of spans within a heightfield.
/// @see rcHeightfield
/// 用于在高度场(heightfield)内快速分配跨度的内存池。
struct rcSpanPool
{
	rcSpanPool* next;					///< The next span pool.
	rcSpan items[RC_SPANS_PER_POOL];	///< Array of spans in the pool.
};

/// A dynamic heightfield representing obstructed space.
/// @ingroup recast
/// 代表受阻空间的动态高度场(heightfield)。
struct rcHeightfield
{
	rcHeightfield();
	~rcHeightfield();

	///< The width of the heightfield. (Along the x-axis in cell units.)
	///< 高度场(heightfield)的宽度。 （沿 x 轴，以单元格为单位。）
	int width;
	///< The height of the heightfield. (Along the z-axis in cell units.)
	///< 高度场(heightfield)的高度。 （沿 z 轴，以单元格为单位。）
	int height;
	///< The minimum bounds in world space. [(x, y, z)]
	///< 世界空间的最小边界。 [（x，y，z）]
	float bmin[3];
	///< The maximum bounds in world space. [(x, y, z)]
	///< 世界空间的最大边界。 [（x，y，z）]
	float bmax[3];
	///< The size of each cell. (On the xz-plane.)
	///< 每个单元格的大小。 （在 xz 平面上。）
	float cs;
	///< The height of each cell. (The minimum increment along the y-axis.)
	///< 每个单元格的高度。 （沿 y 轴的最小增量。）
	float ch;
	///< Heightfield of spans (width*height).
	///< spans 链表头节点指针相当于指针数组
	rcSpan** spans;

	// memory pool for rcSpan instances.
	// rcSpan 实例的内存池。
	rcSpanPool* pools;	///< Linked list of span pools.(span 池的链表)
	rcSpan* freelist;	///< The next free span.(free span的链表)

  private:
	// Explicitly-disabled copy constructor and copy assignment operator.
	rcHeightfield(const rcHeightfield&);
	rcHeightfield& operator=(const rcHeightfield&);
};

/// Provides information on the content of a cell column in a compact heightfield.
/// 提供有关紧凑高度场(heightfield)中单元格列内容的信息。
struct rcCompactCell
{
	///< Index to the first span in the column.
	///< 列中第一个span的索引。
	unsigned int index : 24;
	///< Number of spans in the column.
	///< 列中spans 的数量
	unsigned int count : 8;
};

/// Represents a span of unobstructed space within a compact heightfield.
/// 表示紧凑高度场(compact heightfield)内的一段无障碍空间(空心的部分)。
struct rcCompactSpan
{
	///< The lower extent of the span. (Measured from the heightfield's base.)
	///< 跨度的下限。 （从高度场(heightfield)的底部开始测量。）
	unsigned short y;
	///< The id of the region the span belongs to. (Or zero if not in a region.)
	///< 该跨度所属区域的 ID。 （如果不在某个区域，则为零。）
	unsigned short reg;
	///< Packed neighbor connection data.
	///< 打包的邻居连接数据。
	unsigned int con : 24;
	///< The height of the span.  (Measured from #y.)
	///< 跨度的高度。 （从 #y 开始测量。）
	unsigned int h : 8;
};

/// A compact, static heightfield representing unobstructed space.
/// @ingroup recast
/// 代表无障碍空间的紧凑、静态高度场。
struct rcCompactHeightfield
{
	rcCompactHeightfield();
	~rcCompactHeightfield();

	///< The width of the heightfield. (Along the x-axis in cell units.)
	///< 高度场(heightfield)的宽度。 （沿 x 轴，以单元格为单位。）
	int width;
	///< The height of the heightfield. (Along the z-axis in cell units.)
	///< 高度场(heightfield)的高度。 （沿 z 轴，以单元格为单位。）
	int height;
	///< The number of spans in the heightfield.
	///< 高度场(heightfield)中spans的数量
	int spanCount;
	///< The walkable height used during the build of the field.  (See: rcConfig::walkableHeight)
	///< 可通行高度。 （参见：rcConfig::walkableHeight）
	int walkableHeight;
	///< The walkable climb used during the build of the field. (See: rcConfig::walkableClimb)
	///< 可通行体素高度落差(See: rcConfig::walkableClimb)
	int walkableClimb;
	///< The AABB border size used during the build of the field. (See: rcConfig::borderSize)
	///< 构建数据期间使用的 AABB 边框大小。 （参见：rcConfig::borderSize）
	int borderSize;
	///< The maximum distance value of any span within the field.
	///< 高度场(heightfield)内任意span的最大距离值.
	unsigned short maxDistance;
	///< The maximum region id of any span within the field.
	///< 高度场(heightfield)内任意 span 的最大region id。
	unsigned short maxRegions;
	///< The minimum bounds in world space. [(x, y, z)]
	///< 高度场(heightfield)AABB盒最小边界(世界空间)
	float bmin[3];
	///< The maximum bounds in world space. [(x, y, z)]
	///< 高度场(heightfield)AABB盒最大边界(世界空间)
	float bmax[3];
	///< The size of each cell. (On the xz-plane.)
	///< 每个cell的宽度. (在 xz-plane.)
	float cs;
	///< The height of each cell. (The minimum increment along the y-axis.)
	///< 每个cell的高度。 （沿 y 轴的最小增量。）
	float ch;
	///< Array of cells. [Size: #width*#height]
	///< cell的数组. [Size: #width*#height]
	rcCompactCell* cells;
	///< span的数组. [Size: #spanCount]
	rcCompactSpan* spans;
	///< Array containing border distance data. [Size: #spanCount]
	///< 包含边界距离数据的数组 [Size: #spanCount]
	unsigned short* dist;
	///< Array containing area id data. [Size: #spanCount]
	///< 包含area id 数据的数组 [Size: #spanCount]
	unsigned char* areas;

  private:
	// Explicitly-disabled copy constructor and copy assignment operator.
	rcCompactHeightfield(const rcCompactHeightfield&);
	rcCompactHeightfield& operator=(const rcCompactHeightfield&);
};

/// Represents a heightfield layer within a layer set.
/// @see rcHeightfieldLayerSet
/// 表示图层集(layer set)中的高度场(heightfield)的layer。
struct rcHeightfieldLayer
{
	float bmin[3];				///< The minimum bounds in world space. [(x, y, z)]
	float bmax[3];				///< The maximum bounds in world space. [(x, y, z)]
	float cs;					///< The size of each cell. (On the xz-plane.)
	float ch;					///< The height of each cell. (The minimum increment along the y-axis.)
	int width;					///< The width of the heightfield. (Along the x-axis in cell units.)
	int height;					///< The height of the heightfield. (Along the z-axis in cell units.)
	int minx;					///< The minimum x-bounds of usable data.
	int maxx;					///< The maximum x-bounds of usable data.
	int miny;					///< The minimum y-bounds of usable data. (Along the z-axis.)
	int maxy;					///< The maximum y-bounds of usable data. (Along the z-axis.)
	int hmin;					///< The minimum height bounds of usable data. (Along the y-axis.)
	int hmax;					///< The maximum height bounds of usable data. (Along the y-axis.)
	unsigned char* heights;		///< The heightfield. [Size: width * height]
	unsigned char* areas;		///< Area ids. [Size: Same as #heights]
	unsigned char* cons;		///< Packed neighbor connection information. [Size: Same as #heights]
};

/// Represents a set of heightfield layers.
/// @ingroup recast
/// @see rcAllocHeightfieldLayerSet, rcFreeHeightfieldLayerSet
struct rcHeightfieldLayerSet
{
	rcHeightfieldLayerSet();
	~rcHeightfieldLayerSet();

	rcHeightfieldLayer* layers;			///< The layers in the set. [Size: #nlayers]
	int nlayers;						///< The number of layers in the set.

private:
	// Explicitly-disabled copy constructor and copy assignment operator.
	rcHeightfieldLayerSet(const rcHeightfieldLayerSet&);
	rcHeightfieldLayerSet& operator=(const rcHeightfieldLayerSet&);
};

/// Represents a simple, non-overlapping contour in field space.
/// 表示field space中的简单、不重叠的轮廓。
struct rcContour
{
	///< Simplified contour vertex and connection data. [Size: 4 * #nverts]
	///< 简化的轮廓顶点和连接数据。[Size: 4 * #nverts]
	int* verts;
	///< The number of vertices in the simplified contour.
	///< 简化轮廓中的顶点数。
	int nverts;
	///< Raw contour vertex and connection data. [Size: 4 * #nrverts]
	///< 原始(Raw)轮廓顶点和连接数据。 [大小：4 * #nrverts]
	int* rverts;
	///< 原始轮廓中的顶点数。
	int nrverts;
	///< The region id of the contour.
	///< 轮廓的region ID。
	unsigned short reg;
	///< The area id of the contour.
	///< 轮廓的area ID。
	unsigned char area;
};

/// Represents a group of related contours.
/// @ingroup recast
/// 表示一组相关的轮廓。
struct rcContourSet
{
	rcContourSet();
	~rcContourSet();

	rcContour* conts;	///< An array of the contours in the set. [Size: #nconts]
	int nconts;			///< The number of contours in the set.
	float bmin[3];  	///< The minimum bounds in world space. [(x, y, z)]
	float bmax[3];		///< The maximum bounds in world space. [(x, y, z)]
	float cs;			///< The size of each cell. (On the xz-plane.)
	float ch;			///< The height of each cell. (The minimum increment along the y-axis.)
	int width;			///< The width of the set. (Along the x-axis in cell units.)
	int height;			///< The height of the set. (Along the z-axis in cell units.)
	int borderSize;		///< The AABB border size used to generate the source data from which the contours were derived.
	///< The max edge error that this contour set was simplified with.
	///< 该轮廓集简化后的最大边缘误差。
	float maxError;

  private:
	// Explicitly-disabled copy constructor and copy assignment operator.
	rcContourSet(const rcContourSet&);
	rcContourSet& operator=(const rcContourSet&);
};

/// Represents a polygon mesh suitable for use in building a navigation mesh.
/// @ingroup recast
/// 表示适合用于构建导航网格的多边形网格。
struct rcPolyMesh
{
	rcPolyMesh();
	~rcPolyMesh();

	///< The mesh vertices. [Form: (x, y, z) * #nverts]
	unsigned short* verts;
	///< Polygon and neighbor data. [Length: #maxpolys * 2 * #nvp]
	///< 多边形和邻居数据 [Length: #maxpolys * 2 * #nvp]
	unsigned short* polys;
	///< The region id assigned to each polygon. [Length: #maxpolys]
	///< 分配给每个多边形的区域 ID。 [Length: #maxpolys]
	unsigned short* regs;
	///< The user defined flags for each polygon. [Length: #maxpolys]
	///< 用户为每个多边形定义的标志。[Length: #maxpolys]
	unsigned short* flags;
	///< The area id assigned to each polygon. [Length: #maxpolys]
	///< 分配给每个多边形的区域 ID。[Length: #maxpolys]
	unsigned char* areas;
	///< The number of vertices.
	int nverts;
	///< The number of polygons.
	int npolys;
	///< The number of allocated polygons.
	int maxpolys;
	///< The maximum number of vertices per polygon.
	///< 每个多边形的最大顶点数。
	int nvp;
	///< The minimum bounds in world space. [(x, y, z)]
	///< 世界空间中AABB盒最小边界. [(x, y, z)]
	float bmin[3];
	///< The maximum bounds in world space. [(x, y, z)]
	///< 世界空间中AABB盒最大边界. [(x, y, z)]
	float bmax[3];
	///< The size of each cell. (On the xz-plane.)
	float cs;
	///< The height of each cell. (The minimum increment along the y-axis.)
	float ch;
	///< The AABB border size used to generate the source data from which the mesh was derived.
	///< 用于生成从中导出网格的源数据的 AABB 边界大小。
	int borderSize;
	///< The max error of the polygon edges in the mesh.
	///< 网格中多边形边的最大误差。
	float maxEdgeError;

private:
	// Explicitly-disabled copy constructor and copy assignment operator.
	rcPolyMesh(const rcPolyMesh&);
	rcPolyMesh& operator=(const rcPolyMesh&);
};

/// Contains triangle meshes that represent detailed height data associated
/// with the polygons in its associated polygon mesh object.
/// @ingroup recast
/// 包含三角形网格，表示与其关联的多边形网格对象中的多边形关联的详细高度数据。
struct rcPolyMeshDetail
{
	rcPolyMeshDetail();

	///< The sub-mesh data. [Size: 4*#nmeshes]
	unsigned int* meshes;
	///< The mesh vertices. [Size: 3*#nverts]
	float* verts;
	///< The mesh triangles. [Size: 4*#ntris]
	unsigned char* tris;
	///< The number of sub-meshes defined by #meshes.
	int nmeshes;
	///< The number of vertices in #verts.
	int nverts;
	///< The number of triangles in #tris.
	int ntris;

private:
	// Explicitly-disabled copy constructor and copy assignment operator.
	rcPolyMeshDetail(const rcPolyMeshDetail&);
	rcPolyMeshDetail& operator=(const rcPolyMeshDetail&);
};

/// @name Allocation Functions
/// Functions used to allocate and de-allocate Recast objects.
/// 用于分配和取消分配 Recast 对象的函数。
/// @see rcAllocSetCustom
/// @{

/// Allocates a heightfield object using the Recast allocator.
/// @return A heightfield that is ready for initialization, or null on failure.
/// @ingroup recast
/// @see rcCreateHeightfield, rcFreeHeightField
rcHeightfield* rcAllocHeightfield();

/// Frees the specified heightfield object using the Recast allocator.
/// @param[in]		heightfield	A heightfield allocated using #rcAllocHeightfield
/// @ingroup recast
/// @see rcAllocHeightfield
void rcFreeHeightField(rcHeightfield* heightfield);

/// Allocates a compact heightfield object using the Recast allocator.
/// @return A compact heightfield that is ready for initialization, or null on failure.
/// @ingroup recast
/// @see rcBuildCompactHeightfield, rcFreeCompactHeightfield
rcCompactHeightfield* rcAllocCompactHeightfield();

/// Frees the specified compact heightfield object using the Recast allocator.
/// @param[in]		compactHeightfield		A compact heightfield allocated using #rcAllocCompactHeightfield
/// @ingroup recast
/// @see rcAllocCompactHeightfield
void rcFreeCompactHeightfield(rcCompactHeightfield* compactHeightfield);

/// Allocates a heightfield layer set using the Recast allocator.
/// @return A heightfield layer set that is ready for initialization, or null on failure.
/// @ingroup recast
/// @see rcBuildHeightfieldLayers, rcFreeHeightfieldLayerSet
rcHeightfieldLayerSet* rcAllocHeightfieldLayerSet();

/// Frees the specified heightfield layer set using the Recast allocator.
/// @param[in]		layerSet	A heightfield layer set allocated using #rcAllocHeightfieldLayerSet
/// @ingroup recast
/// @see rcAllocHeightfieldLayerSet
void rcFreeHeightfieldLayerSet(rcHeightfieldLayerSet* layerSet);

/// Allocates a contour set object using the Recast allocator.
/// @return A contour set that is ready for initialization, or null on failure.
/// @ingroup recast
/// @see rcBuildContours, rcFreeContourSet
rcContourSet* rcAllocContourSet();

/// Frees the specified contour set using the Recast allocator.
/// @param[in]		contourSet	A contour set allocated using #rcAllocContourSet
/// @ingroup recast
/// @see rcAllocContourSet
void rcFreeContourSet(rcContourSet* contourSet);

/// Allocates a polygon mesh object using the Recast allocator.
/// @return A polygon mesh that is ready for initialization, or null on failure.
/// @ingroup recast
/// @see rcBuildPolyMesh, rcFreePolyMesh
rcPolyMesh* rcAllocPolyMesh();

/// Frees the specified polygon mesh using the Recast allocator.
/// @param[in]		polyMesh	A polygon mesh allocated using #rcAllocPolyMesh
/// @ingroup recast
/// @see rcAllocPolyMesh
void rcFreePolyMesh(rcPolyMesh* polyMesh);

/// Allocates a detail mesh object using the Recast allocator.
/// @return A detail mesh that is ready for initialization, or null on failure.
/// @ingroup recast
/// @see rcBuildPolyMeshDetail, rcFreePolyMeshDetail
rcPolyMeshDetail* rcAllocPolyMeshDetail();

/// Frees the specified detail mesh using the Recast allocator.
/// @param[in]		detailMesh	A detail mesh allocated using #rcAllocPolyMeshDetail
/// @ingroup recast
/// @see rcAllocPolyMeshDetail
void rcFreePolyMeshDetail(rcPolyMeshDetail* detailMesh);

/// @}

/// Heightfield border flag.
/// If a heightfield region ID has this bit set, then the region is a border
/// region and its spans are considered un-walkable.
/// (Used during the region and contour build process.)
/// @see rcCompactSpan::reg
/// 海特菲尔德边界标志。 如果高度场region ID 设置了此位，则该区域是边界区域，并且其span被认为是不可行走的。 （在区域和轮廓构建过程中使用。）
static const unsigned short RC_BORDER_REG = 0x8000;

/// Polygon touches multiple regions.
/// If a polygon has this region ID it was merged with or created
/// from polygons of different regions during the polymesh
/// build step that removes redundant border vertices.
/// (Used during the polymesh and detail polymesh build processes)
/// @see rcPolyMesh::regs
/// 多边形接触多个区域。 如果多边形具有此区域 ID，则在删除冗余边界顶点的多边形网格构建步骤期间，它将与不同区域的多边形合并或从不同区域的多边形创建。
/// （在多边形网格和细节多边形网格构建过程中使用）
static const unsigned short RC_MULTIPLE_REGS = 0;

/// Border vertex flag.
/// If a region ID has this bit set, then the associated element lies on
/// a tile border. If a contour vertex's region ID has this bit set, the
/// vertex will later be removed in order to match the segments and vertices
/// at tile boundaries.
/// (Used during the build process.)
/// @see rcCompactSpan::reg, #rcContour::verts, #rcContour::rverts
/// 边界顶点标志。 如果区域 ID 设置了此位，则关联元素位于tile边框上。 如果轮廓顶点的区域 ID 设置了此位，则稍后将删除该顶点，
/// 以便匹配图块边界处的线段和顶点。 （在构建过程中使用。）
static const int RC_BORDER_VERTEX = 0x10000;

/// Area border flag.
/// If a region ID has this bit set, then the associated element lies on
/// the border of an area.
/// (Used during the region and contour build process.)
/// @see rcCompactSpan::reg, #rcContour::verts, #rcContour::rverts
/// area边界标志。 如果region ID 设置了该位，则关联元素位于area的边界上。 （在region和轮廓构建过程中使用。）
static const int RC_AREA_BORDER = 0x20000;

/// Contour build flags.
/// @see rcBuildContours
enum rcBuildContoursFlags
{
	///< Tessellate solid (impassable) edges during contour simplification.
	///< 在轮廓简化过程中对实体（不可通过的）边缘进行细分。
	RC_CONTOUR_TESS_WALL_EDGES = 0x01,
	///< Tessellate edges between areas during contour simplification.
	///< 在轮廓简化期间对area之间的边缘进行细分。
	RC_CONTOUR_TESS_AREA_EDGES = 0x02
};

/// Applied to the region id field of contour vertices in order to extract the region id.
/// The region id field of a vertex may have several flags applied to it.  So the
/// fields value can't be used directly.
/// @see rcContour::verts, rcContour::rverts
/// 应用于轮廓顶点的region id 字段以提取region id。 顶点的region id 字段可能应用了多个标志。 所以不能直接使用字段值。
static const int RC_CONTOUR_REG_MASK = 0xffff;

/// An value which indicates an invalid index within a mesh.
/// @note This does not necessarily indicate an error.
/// @see rcPolyMesh::polys
/// 指示网格内无效索引的值。
/// @note 这并不一定表示错误。
static const unsigned short RC_MESH_NULL_IDX = 0xffff;

/// Represents the null area.
/// When a data element is given this value it is considered to no longer be
/// assigned to a usable area.  (E.g. It is un-walkable.)
/// 代表空区。 当一个数据元素被赋予这个值时，它被认为不再被分配到可用区域。 （例如，这里不适合步行。）
static const unsigned char RC_NULL_AREA = 0;

/// The default area id used to indicate a walkable polygon.
/// This is also the maximum allowed area id, and the only non-null area id
/// recognized by some steps in the build process.
/// 用于指示可步行多边形的默认区域 ID。 这也是允许的最大区域 ID，也是构建过程中某些步骤识别的唯一非空区域 ID。
static const unsigned char RC_WALKABLE_AREA = 63;

/// The value returned by #rcGetCon if the specified direction is not connected
/// to another span. (Has no neighbor.)
/// 如果指定方向未连接到另一个span，则返回 #rcGetCon 的值。 （没有邻居。）
static const int RC_NOT_CONNECTED = 0x3f;

/// @name General helper functions
/// @{

/// Swaps the values of the two parameters.
/// @param[in,out]	a	Value A
/// @param[in,out]	b	Value B
template<class T> inline void rcSwap(T& a, T& b) { T t = a; a = b; b = t; }

/// Returns the minimum of two values.
/// @param[in]		a	Value A
/// @param[in]		b	Value B
/// @return The minimum of the two values.
template<class T> inline T rcMin(T a, T b) { return a < b ? a : b; }

/// Returns the maximum of two values.
/// @param[in]		a	Value A
/// @param[in]		b	Value B
/// @return The maximum of the two values.
template<class T> inline T rcMax(T a, T b) { return a > b ? a : b; }

/// Returns the absolute value.
/// @param[in]		a	The value.
/// @return The absolute value of the specified value.
template<class T> inline T rcAbs(T a) { return a < 0 ? -a : a; }

/// Returns the square of the value.
/// @param[in]		a	The value.
/// @return The square of the value.
template<class T> inline T rcSqr(T a) { return a * a; }

/// Clamps the value to the specified range.
/// @param[in]		value			The value to clamp.
/// @param[in]		minInclusive	The minimum permitted return value.
/// @param[in]		maxInclusive	The maximum permitted return value.
/// @return The value, clamped to the specified range.
template<class T> inline T rcClamp(T value, T minInclusive, T maxInclusive)
{
	return value < minInclusive ? minInclusive: (value > maxInclusive ? maxInclusive : value);
}

/// Returns the square root of the value.
///  @param[in]		x	The value.
///  @return The square root of the vlaue.
/// 返回值的平方根。
float rcSqrt(float x);

/// @}
/// @name Vector helper functions.
/// @{

/// Derives the cross product of two vectors. (@p v1 x @p v2)
/// @param[out]		dest	The cross product. [(x, y, z)]
/// @param[in]		v1		A Vector [(x, y, z)]
/// @param[in]		v2		A vector [(x, y, z)]
/// 导出两个向量的叉积。 （@p v1 x @p v2）
inline void rcVcross(float* dest, const float* v1, const float* v2)
{
	dest[0] = v1[1]*v2[2] - v1[2]*v2[1];
	dest[1] = v1[2]*v2[0] - v1[0]*v2[2];
	dest[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

/// Derives the dot product of two vectors. (@p v1 . @p v2)
/// @param[in]		v1	A Vector [(x, y, z)]
/// @param[in]		v2	A vector [(x, y, z)]
/// @return The dot product.
/// 导出两个向量的点积。 （@p v1 。@p v2）
inline float rcVdot(const float* v1, const float* v2)
{
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2];
}

/// Performs a scaled vector addition. (@p v1 + (@p v2 * @p s))
/// @param[out]		dest	The result vector. [(x, y, z)]
/// @param[in]		v1		The base vector. [(x, y, z)]
/// @param[in]		v2		The vector to scale and add to @p v1. [(x, y, z)]
/// @param[in]		s		The amount to scale @p v2 by before adding to @p v1.
/// 执行缩放向量加法。 (@p v1 + (@p v2 * @p s))
inline void rcVmad(float* dest, const float* v1, const float* v2, const float s)
{
	dest[0] = v1[0]+v2[0]*s;
	dest[1] = v1[1]+v2[1]*s;
	dest[2] = v1[2]+v2[2]*s;
}

/// Performs a vector addition. (@p v1 + @p v2)
/// @param[out]		dest	The result vector. [(x, y, z)]
/// @param[in]		v1		The base vector. [(x, y, z)]
/// @param[in]		v2		The vector to add to @p v1. [(x, y, z)]
/// 执行向量加法。 （@p v1 + @p v2）
inline void rcVadd(float* dest, const float* v1, const float* v2)
{
	dest[0] = v1[0]+v2[0];
	dest[1] = v1[1]+v2[1];
	dest[2] = v1[2]+v2[2];
}

/// Performs a vector subtraction. (@p v1 - @p v2)
/// @param[out]		dest	The result vector. [(x, y, z)]
/// @param[in]		v1		The base vector. [(x, y, z)]
/// @param[in]		v2		The vector to subtract from @p v1. [(x, y, z)]
/// 执行矢量减法。 （@p v1 - @p v2）
inline void rcVsub(float* dest, const float* v1, const float* v2)
{
	dest[0] = v1[0]-v2[0];
	dest[1] = v1[1]-v2[1];
	dest[2] = v1[2]-v2[2];
}

/// Selects the minimum value of each element from the specified vectors.
/// @param[in,out]	mn	A vector.  (Will be updated with the result.) [(x, y, z)]
/// @param[in]		v	A vector. [(x, y, z)]
/// 从指定向量中选择每个分量的最小值。
inline void rcVmin(float* mn, const float* v)
{
	mn[0] = rcMin(mn[0], v[0]);
	mn[1] = rcMin(mn[1], v[1]);
	mn[2] = rcMin(mn[2], v[2]);
}

/// Selects the maximum value of each element from the specified vectors.
/// @param[in,out]	mx	A vector.  (Will be updated with the result.) [(x, y, z)]
/// @param[in]		v	A vector. [(x, y, z)]
/// 从指定向量中选择每个分量的最大值。
inline void rcVmax(float* mx, const float* v)
{
	mx[0] = rcMax(mx[0], v[0]);
	mx[1] = rcMax(mx[1], v[1]);
	mx[2] = rcMax(mx[2], v[2]);
}

/// Performs a vector copy.
/// @param[out]		dest	The result. [(x, y, z)]
/// @param[in]		v		The vector to copy. [(x, y, z)]
inline void rcVcopy(float* dest, const float* v)
{
	dest[0] = v[0];
	dest[1] = v[1];
	dest[2] = v[2];
}

/// Returns the distance between two points.
/// @param[in]		v1	A point. [(x, y, z)]
/// @param[in]		v2	A point. [(x, y, z)]
/// @return The distance between the two points.
/// 返回两点之间的距离。
inline float rcVdist(const float* v1, const float* v2)
{
	float dx = v2[0] - v1[0];
	float dy = v2[1] - v1[1];
	float dz = v2[2] - v1[2];
	return rcSqrt(dx*dx + dy*dy + dz*dz);
}

/// Returns the square of the distance between two points.
/// @param[in]		v1	A point. [(x, y, z)]
/// @param[in]		v2	A point. [(x, y, z)]
/// @return The square of the distance between the two points.
/// 返回两点之间距离的平方。
inline float rcVdistSqr(const float* v1, const float* v2)
{
	float dx = v2[0] - v1[0];
	float dy = v2[1] - v1[1];
	float dz = v2[2] - v1[2];
	return dx*dx + dy*dy + dz*dz;
}

/// Normalizes the vector.
/// @param[in,out]	v	The vector to normalize. [(x, y, z)]
/// 归一化向量
inline void rcVnormalize(float* v)
{
	float d = 1.0f / rcSqrt(rcSqr(v[0]) + rcSqr(v[1]) + rcSqr(v[2]));
	v[0] *= d;
	v[1] *= d;
	v[2] *= d;
}

/// @}
/// @name Heightfield Functions
/// @see rcHeightfield
/// @{

/// Calculates the bounding box of an array of vertices.
/// @ingroup recast
/// @param[in]		verts		An array of vertices. [(x, y, z) * @p nv]
/// @param[in]		numVerts	The number of vertices in the @p verts array.
/// @param[out]		minBounds	The minimum bounds of the AABB. [(x, y, z)] [Units: wu]
/// @param[out]		maxBounds	The maximum bounds of the AABB. [(x, y, z)] [Units: wu]
/// 计算顶点数组的边界框。
void rcCalcBounds(const float* verts, int numVerts, float* minBounds, float* maxBounds);

/// Calculates the grid size based on the bounding box and grid cell size.
/// @ingroup recast
/// @param[in]		minBounds	The minimum bounds of the AABB. [(x, y, z)] [Units: wu]
/// @param[in]		maxBounds	The maximum bounds of the AABB. [(x, y, z)] [Units: wu]
/// @param[in]		cellSize	The xz-plane cell size. [Limit: > 0] [Units: wu]
/// @param[out]		sizeX		The width along the x-axis. [Limit: >= 0] [Units: vx]
/// @param[out]		sizeZ		The height along the z-axis. [Limit: >= 0] [Units: vx]
/// 根据边界框和网格单元大小计算网格大小。
void rcCalcGridSize(const float* minBounds, const float* maxBounds, float cellSize, int* sizeX, int* sizeZ);

/// Initializes a new heightfield.
/// See the #rcConfig documentation for more information on the configuration parameters.
///
/// @see rcAllocHeightfield, rcHeightfield
/// @ingroup recast
///
/// @param[in,out]	context		The build context to use during the operation.
/// @param[in,out]	heightfield	The allocated heightfield to initialize.
/// @param[in]		sizeX		The width of the field along the x-axis. [Limit: >= 0] [Units: vx]
/// @param[in]		sizeZ		The height of the field along the z-axis. [Limit: >= 0] [Units: vx]
/// @param[in]		minBounds	The minimum bounds of the field's AABB. [(x, y, z)] [Units: wu]
/// @param[in]		maxBounds	The maximum bounds of the field's AABB. [(x, y, z)] [Units: wu]
/// @param[in]		cellSize	The xz-plane cell size to use for the field. [Limit: > 0] [Units: wu]
/// @param[in]		cellHeight	The y-axis cell size to use for field. [Limit: > 0] [Units: wu]
/// @returns True if the operation completed successfully.
bool rcCreateHeightfield(rcContext* context, rcHeightfield& heightfield, int sizeX, int sizeZ,
						 const float* minBounds, const float* maxBounds,
						 float cellSize, float cellHeight);

/// Sets the area id of all triangles with a slope below the specified value
/// to #RC_WALKABLE_AREA.
///
/// Only sets the area id's for the walkable triangles.  Does not alter the
/// area id's for un-walkable triangles.
///
/// See the #rcConfig documentation for more information on the configuration parameters.
///
/// @see rcHeightfield, rcClearUnwalkableTriangles, rcRasterizeTriangles
///
/// @ingroup recast
/// @param[in,out]	context				The build context to use during the operation.
/// @param[in]		walkableSlopeAngle	The maximum slope that is considered walkable.
/// 									[Limits: 0 <= value < 90] [Units: Degrees]
/// @param[in]		verts				The vertices. [(x, y, z) * @p nv]
/// @param[in]		numVerts			The number of vertices.
/// @param[in]		tris				The triangle vertex indices. [(vertA, vertB, vertC) * @p nt]
/// @param[in]		numTris				The number of triangles.
/// @param[out]		triAreaIDs			The triangle area ids. [Length: >= @p nt]
///
/// 将所有斜率低于指定值的三角形的区域 ID 设置为 #RC_WALKABLE_AREA。
/// 仅设置可步行三角形的area ID。 不改变不可行走三角形的area ID。
///
/// 有关配置参数的更多信息，请参阅#rcConfig 文档。
void rcMarkWalkableTriangles(rcContext* context, float walkableSlopeAngle, const float* verts, int numVerts,
							 const int* tris, int numTris, unsigned char* triAreaIDs);

/// Sets the area id of all triangles with a slope greater than or equal to the specified value to #RC_NULL_AREA.
///
/// Only sets the area id's for the un-walkable triangles.  Does not alter the
/// area id's for walkable triangles.
///
/// See the #rcConfig documentation for more information on the configuration parameters.
///
/// @see rcHeightfield, rcClearUnwalkableTriangles, rcRasterizeTriangles
///
/// @ingroup recast
/// @param[in,out]	context				The build context to use during the operation.
/// @param[in]		walkableSlopeAngle	The maximum slope that is considered walkable.
/// 									[Limits: 0 <= value < 90] [Units: Degrees]
/// @param[in]		verts				The vertices. [(x, y, z) * @p nv]
/// @param[in]		numVerts			The number of vertices.
/// @param[in]		tris				The triangle vertex indices. [(vertA, vertB, vertC) * @p nt]
/// @param[in]		numTris				The number of triangles.
/// @param[out]		triAreaIDs			The triangle area ids. [Length: >= @p nt]
///
/// 将所有斜率大于或等于指定值的三角形的面积id设置为#RC_NULL_AREA。
/// 只为不可行走的三角形设置area ID。 不改变可步行三角形的area ID。
/// 有关配置参数的更多信息，请参阅#rcConfig 文档。
void rcClearUnwalkableTriangles(rcContext* context, float walkableSlopeAngle, const float* verts, int numVerts,
								const int* tris, int numTris, unsigned char* triAreaIDs);

/// Adds a span to the specified heightfield.
///
/// The span addition can be set to favor flags. If the span is merged to
/// another span and the new @p spanMax is within @p flagMergeThreshold units
/// from the existing span, the span flags are merged.
///
/// @ingroup recast
/// @param[in,out]	context				The build context to use during the operation.
/// @param[in,out]	heightfield			An initialized heightfield.
/// @param[in]		x					The column x index where the span is to be added.
/// 									[Limits: 0 <= value < rcHeightfield::width]
/// @param[in]		z					The column z index where the span is to be added.
/// 									[Limits: 0 <= value < rcHeightfield::height]
/// @param[in]		spanMin				The minimum height of the span. [Limit: < @p spanMax] [Units: vx]
/// @param[in]		spanMax				The maximum height of the span. [Limit: <= #RC_SPAN_MAX_HEIGHT] [Units: vx]
/// @param[in]		areaID				The area id of the span. [Limit: <= #RC_WALKABLE_AREA)
/// @param[in]		flagMergeThreshold	The merge threshold(合并阈值). [Limit: >= 0] [Units: vx]
/// @returns True if the operation completed successfully.
///
/// 将apan添加到指定的高度场。
///
/// 可以设置span添加以支持标志。 如果该span合并到另一个span，并且新的 @p spanMax
/// 位于现有span的 @p flagMergeThreshold 单位内，则span标志将被合并。
bool rcAddSpan(rcContext* context, rcHeightfield& heightfield,
	           int x, int z,
               unsigned short spanMin, unsigned short spanMax,
               unsigned char areaID, int flagMergeThreshold);

/// Rasterizes a single triangle into the specified heightfield.
///
/// Calling this for each triangle in a mesh is less efficient than calling rcRasterizeTriangles
///
/// No spans will be added if the triangle does not overlap the heightfield grid.
///
/// @see rcHeightfield
/// @ingroup recast
/// @param[in,out]	context				The build context to use during the operation.
/// @param[in]		v0					Triangle vertex 0 [(x, y, z)]
/// @param[in]		v1					Triangle vertex 1 [(x, y, z)]
/// @param[in]		v2					Triangle vertex 2 [(x, y, z)]
/// @param[in]		areaID				The area id of the triangle. [Limit: <= #RC_WALKABLE_AREA]
/// @param[in,out]	heightfield			An initialized heightfield.
/// @param[in]		flagMergeThreshold	The distance where the walkable flag is favored over the non-walkable flag.
/// 									[Limit: >= 0] [Units: vx]
/// @returns True if the operation completed successfully.
///
/// 将单个三角形光栅化到指定的高度场中。
/// 对网格中的每个三角形调用此方法的效率低于调用 rcRasterizeTriangles
/// 如果三角形不与高度场网格重叠，则不会添加span。
bool rcRasterizeTriangle(rcContext* context,
                         const float* v0, const float* v1, const float* v2,
                         unsigned char areaID, rcHeightfield& heightfield, int flagMergeThreshold = 1);

/// Rasterizes an indexed triangle mesh into the specified heightfield.
///
/// Spans will only be added for triangles that overlap the heightfield grid.
///
/// @see rcHeightfield
/// @ingroup recast
/// @param[in,out]	context				The build context to use during the operation.
/// @param[in]		verts				The vertices. [(x, y, z) * @p nv]
/// @param[in]		numVerts			The number of vertices. (unused) TODO (graham): Remove in next major release
/// @param[in]		tris				The triangle indices. [(vertA, vertB, vertC) * @p nt]
/// @param[in]		triAreaIDs			The area id's of the triangles. [Limit: <= #RC_WALKABLE_AREA] [Size: @p nt]
/// @param[in]		numTris				The number of triangles.
/// @param[in,out]	heightfield			An initialized heightfield.
/// @param[in]		flagMergeThreshold	The distance where the walkable flag is favored over the non-walkable flag.
///										[Limit: >= 0] [Units: vx]
/// @returns True if the operation completed successfully.
///
/// 将索引三角形网格栅格化到指定的高度场中。
/// 仅对与高度场网格重叠的三角形添加span。
bool rcRasterizeTriangles(rcContext* context,
                          const float* verts, int numVerts,
                          const int* tris, const unsigned char* triAreaIDs, int numTris,
                          rcHeightfield& heightfield, int flagMergeThreshold = 1);

/// Rasterizes an indexed triangle mesh into the specified heightfield.
///
/// Spans will only be added for triangles that overlap the heightfield grid.
///
/// @see rcHeightfield
/// @ingroup recast
/// @param[in,out]	context				The build context to use during the operation.
/// @param[in]		verts				The vertices. [(x, y, z) * @p nv]
/// @param[in]		numVerts			The number of vertices. (unused) TODO (graham): Remove in next major release
/// @param[in]		tris				The triangle indices. [(vertA, vertB, vertC) * @p nt]
/// @param[in]		triAreaIDs			The area id's of the triangles. [Limit: <= #RC_WALKABLE_AREA] [Size: @p nt]
/// @param[in]		numTris				The number of triangles.
/// @param[in,out]	heightfield			An initialized heightfield.
/// @param[in]		flagMergeThreshold	The distance where the walkable flag is favored over the non-walkable flag.
/// 									[Limit: >= 0] [Units: vx]
/// @returns True if the operation completed successfully.
///
/// 将索引三角形网格栅格化到指定的高度场中。
/// 仅对与高度场网格重叠的三角形添加span。
bool rcRasterizeTriangles(rcContext* context,
                          const float* verts, int numVerts,
                          const unsigned short* tris, const unsigned char* triAreaIDs, int numTris,
                          rcHeightfield& heightfield, int flagMergeThreshold = 1);

/// Rasterizes a triangle list into the specified heightfield.
///
/// Expects each triangle to be specified as three sequential vertices of 3 floats.
///
/// Spans will only be added for triangles that overlap the heightfield grid.
///
/// @see rcHeightfield
/// @ingroup recast
/// @param[in,out]	context				The build context to use during the operation.
/// @param[in]		verts				The triangle vertices. [(ax, ay, az, bx, by, bz, cx, by, cx) * @p nt]
/// @param[in]		triAreaIDs			The area id's of the triangles. [Limit: <= #RC_WALKABLE_AREA] [Size: @p nt]
/// @param[in]		numTris				The number of triangles.
/// @param[in,out]	heightfield			An initialized heightfield.
/// @param[in]		flagMergeThreshold	The distance where the walkable flag is favored over the non-walkable flag.
/// 									[Limit: >= 0] [Units: vx]
/// @returns True if the operation completed successfully.
///
/// 将三角形列表光栅化到指定的高度场中。
/// 期望将每个三角形指定为 3 个浮点数的三个连续顶点。
/// 仅对与高度场网格重叠的三角形添加span。
bool rcRasterizeTriangles(rcContext* context,
                          const float* verts, const unsigned char* triAreaIDs, int numTris,
                          rcHeightfield& heightfield, int flagMergeThreshold = 1);

/// Marks non-walkable spans as walkable if their maximum is within @p walkableClimb of the span below them.
///
/// This removes small obstacles and rasterization artifacts that the agent would be able to walk over
/// such as curbs.  It also allows agents to move up terraced structures like stairs.
///
/// Obstacle spans are marked walkable if: <tt>obstacleSpan.smax - walkableSpan.smax < walkableClimb</tt>
///
/// @warning Will override the effect of #rcFilterLedgeSpans.  If both filters are used, call #rcFilterLedgeSpans only after applying this filter.
///
/// @see rcHeightfield, rcConfig
///
/// @ingroup recast
/// @param[in,out]	context			The build context to use during the operation.
/// @param[in]		walkableClimb	Maximum ledge height that is considered to still be traversable.
/// 								[Limit: >=0] [Units: vx]
/// @param[in,out]	heightfield		A fully built heightfield.  (All spans have been added.)
///
/// 如果不可通行span的最大值在其下方span的 @p walkableClimb 范围内，则将不可通行span标记为可通行。
/// 这消除了Agent能够走过的小障碍物和光栅化伪像，例如路缘石。 它还允许Agent向上移动楼梯等阶梯结构。
/// 如果满足以下条件，则障碍跨度被标记为可步行：<tt>obstacleSpan.smax - walkableSpan.smax < walkableClimb</tt>
/// @warning 将覆盖 #rcFilterLedgeSpans 的效果。 如果使用两个过滤器，则仅在应用此过滤器后调用#rcFilterLedgeSpans。
///
/// @see rcHeightfield, rcConfig
void rcFilterLowHangingWalkableObstacles(rcContext* context, int walkableClimb, rcHeightfield& heightfield);

/// Marks spans that are ledges as not-walkable.
///
/// A ledge is a span with one or more neighbors whose maximum is further away than @p walkableClimb
/// from the current span's maximum.
/// This method removes the impact of the overestimation of conservative voxelization
/// so the resulting mesh will not have regions hanging in the air over ledges.
///
/// A span is a ledge if: <tt>rcAbs(currentSpan.smax - neighborSpan.smax) > walkableClimb</tt>
///
/// @see rcHeightfield, rcConfig
///
/// @ingroup recast
/// @param[in,out]	context				The build context to use during the operation.
/// @param[in]		walkableHeight	Minimum floor to 'ceiling' height that will still allow the floor area to
/// 								be considered walkable. [Limit: >= 3] [Units: vx]
/// @param[in]		walkableClimb	Maximum ledge height that is considered to still be traversable.
/// 								[Limit: >=0] [Units: vx]
/// @param[in,out]	heightfield			A fully built heightfield.  (All spans have been added.)
///
/// 将ledges的span标记为不可行走。
///
/// ledges是一个具有一个或多个邻居的span，其最大值距离当前跨度的最大值比 @p walkableClimb 更远。
/// 此方法消除了保守体素化高估的影响，因此生成的网格不会有悬挂在壁架上方空气中的区域。
///
/// 如果满足以下条件，则span是ledges：<tt>rcAbs(currentSpan.smax - neighborSpan.smax) > walkableClimb</tt>
void rcFilterLedgeSpans(rcContext* context, int walkableHeight, int walkableClimb, rcHeightfield& heightfield);

/// Marks walkable spans as not walkable if the clearance above the span is less than the specified walkableHeight.
///
/// For this filter, the clearance above the span is the distance from the span's
/// maximum to the minimum of the next higher span in the same column.
/// If there is no higher span in the column, the clearance is computed as the
/// distance from the top of the span to the maximum heightfield height.
///
/// @see rcHeightfield, rcConfig
/// @ingroup recast
///
/// @param[in,out]	context			The build context to use during the operation.
/// @param[in]		walkableHeight	Minimum floor to 'ceiling' height that will still allow the floor area to
/// 								be considered walkable. [Limit: >= 3] [Units: vx]
/// @param[in,out]	heightfield		A fully built heightfield.  (All spans have been added.)
///
/// 如果span上方的间隙小于指定的通行高度，则将可通行span标记为不可通行。
/// 对于此过滤器，跨度上方的间隙是从跨度的最大值到同一列中下一个较高跨度的最小值的距离。
/// 如果柱中没有更高的跨度，则间隙计算为从跨度顶部到最大高度场高度的距离。
void rcFilterWalkableLowHeightSpans(rcContext* context, int walkableHeight, rcHeightfield& heightfield);

/// Returns the number of spans contained in the specified heightfield.
///  @ingroup recast
///  @param[in,out]	context		The build context to use during the operation.
///  @param[in]		heightfield	An initialized heightfield.
///  @returns The number of spans in the heightfield.
///
/// 返回指定高度场(heightfield)中包含的span数。
int rcGetHeightFieldSpanCount(rcContext* context, const rcHeightfield& heightfield);

/// @}
/// @name Compact Heightfield Functions
/// @see rcCompactHeightfield
/// @{

/// Builds a compact heightfield representing open space, from a heightfield representing solid space.
///
/// This is just the beginning of the process of fully building a compact heightfield.
/// Various filters may be applied, then the distance field and regions built.
/// E.g: #rcBuildDistanceField and #rcBuildRegions
///
/// See the #rcConfig documentation for more information on the configuration parameters.
///
/// @see rcAllocCompactHeightfield, rcHeightfield, rcCompactHeightfield, rcConfig
/// @ingroup recast
///
/// @param[in,out]	context				The build context to use during the operation.
/// @param[in]		walkableHeight		Minimum floor to 'ceiling' height that will still allow the floor area
/// 									to be considered walkable. [Limit: >= 3] [Units: vx]
/// @param[in]		walkableClimb		Maximum ledge height that is considered to still be traversable.
/// 									[Limit: >=0] [Units: vx]
/// @param[in]		heightfield			The heightfield to be compacted.
/// @param[out]		compactHeightfield	The resulting compact heightfield. (Must be pre-allocated.)
/// @returns True if the operation completed successfully.
///
/// 从表示实体空间的高度场构建表示开放空间的紧凑高度场。
///
/// 这只是全面构建紧凑高度场过程的开始。 可以应用各种过滤器，然后构建距离场(distance field)和区域(regions)。
/// 例如：#rcBuildDistanceField 和 #rcBuildRegions
///
/// 有关配置参数的更多信息，请参阅#rcConfig 文档。
bool rcBuildCompactHeightfield(rcContext* context, int walkableHeight, int walkableClimb,
							   const rcHeightfield& heightfield, rcCompactHeightfield& compactHeightfield);

/// Erodes the walkable area within the heightfield by the specified radius.
///
/// Basically, any spans that are closer to a boundary or obstruction than the specified radius
/// are marked as un-walkable.
///
/// This method is usually called immediately after the heightfield has been built.
///
/// @see rcCompactHeightfield, rcBuildCompactHeightfield, rcConfig::walkableRadius
/// @ingroup recast
///
/// @param[in,out]	context				The build context to use during the operation.
/// @param[in]		erosionRadius		The radius of erosion. [Limits: 0 < value < 255] [Units: vx]
/// @param[in,out]	compactHeightfield	The populated compact heightfield to erode.
/// @returns True if the operation completed successfully.
///
/// 按指定半径侵蚀高度场内的可步行区域。
///
/// 基本上，任何距离边界或障碍物比指定半径更近的跨度都被标记为不可行走。
///
/// 通常在构建高度场后立即调用此方法。
bool rcErodeWalkableArea(rcContext* context, int erosionRadius, rcCompactHeightfield& compactHeightfield);

/// Applies a median filter to walkable area types (based on area id), removing noise.
///
/// This filter is usually applied after applying area id's using functions
/// such as #rcMarkBoxArea, #rcMarkConvexPolyArea, and #rcMarkCylinderArea.
///
/// @see rcCompactHeightfield
/// @ingroup recast
///
/// @param[in,out]	context		The build context to use during the operation.
/// @param[in,out]	compactHeightfield		A populated compact heightfield.
/// @returns True if the operation completed successfully.
///
/// 将中值过滤器应用于可步行区域类型（基于区域 ID），消除噪声。
/// 此过滤器通常在使用 #rcMarkBoxArea、#rcMarkConvexPolyArea 和 #rcMarkCylinderArea 等函数应用区域 id 之后应用。
bool rcMedianFilterWalkableArea(rcContext* context, rcCompactHeightfield& compactHeightfield);

/// Applies an area id to all spans within the specified bounding box. (AABB)
///
/// @see rcCompactHeightfield, rcMedianFilterWalkableArea
/// @ingroup recast
///
/// @param[in,out]	context				The build context to use during the operation.
/// @param[in]		boxMinBounds		The minimum extents of the bounding box. [(x, y, z)] [Units: wu]
/// @param[in]		boxMaxBounds		The maximum extents of the bounding box. [(x, y, z)] [Units: wu]
/// @param[in]		areaId				The area id to apply. [Limit: <= #RC_WALKABLE_AREA]
/// @param[in,out]	compactHeightfield	A populated compact heightfield.
///
/// 将区域 ID 应用于指定边界框内的所有span。 （AABB）
void rcMarkBoxArea(rcContext* context, const float* boxMinBounds, const float* boxMaxBounds, unsigned char areaId,
				   rcCompactHeightfield& compactHeightfield);

/// Applies the area id to the all spans within the specified convex polygon.
///
/// The value of spacial parameters are in world units.
///
/// The y-values of the polygon vertices are ignored. So the polygon is effectively
/// projected onto the xz-plane, translated to @p minY, and extruded to @p maxY.
///
/// @see rcCompactHeightfield, rcMedianFilterWalkableArea
/// @ingroup recast
///
/// @param[in,out]	context				The build context to use during the operation.
/// @param[in]		verts				The vertices of the polygon [For: (x, y, z) * @p numVerts]
/// @param[in]		numVerts			The number of vertices in the polygon.
/// @param[in]		minY				The height of the base of the polygon. [Units: wu]
/// @param[in]		maxY				The height of the top of the polygon. [Units: wu]
/// @param[in]		areaId				The area id to apply. [Limit: <= #RC_WALKABLE_AREA]
/// @param[in,out]	compactHeightfield	A populated compact heightfield.
///
/// 将area ID 应用于指定凸多边形内的所有span。
/// 空间参数的值采用世界单位。
///
/// 多边形顶点的 y 值将被忽略。 因此，多边形被有效地投影到 xz 平面上，转换为 @p minY，并拉伸为 @p maxY。
void rcMarkConvexPolyArea(rcContext* context, const float* verts, int numVerts,
						  float minY, float maxY, unsigned char areaId,
						  rcCompactHeightfield& compactHeightfield);

/// Expands a convex polygon along its vertex normals by the given offset amount.
/// Inserts extra vertices to bevel sharp corners.
///
/// Helper function to offset convex polygons for rcMarkConvexPolyArea.
///
/// @ingroup recast
///
/// @param[in]		verts		The vertices of the polygon [Form: (x, y, z) * @p numVerts]
/// @param[in]		numVerts	The number of vertices in the polygon.
/// @param[in]		offset		How much to offset the polygon by. [Units: wu]
/// @param[out]		outVerts	The offset vertices (should hold up to 2 * @p numVerts) [Form: (x, y, z) * return value]
/// @param[in]		maxOutVerts	The max number of vertices that can be stored to @p outVerts.
/// @returns Number of vertices in the offset polygon or 0 if too few vertices in @p outVerts.
///
/// 将凸多边形沿其顶点法向量扩展给定的偏移量。 插入额外的顶点以使尖角成斜角。
/// 用于偏移 rcMarkConvexPolyArea 凸多边形的辅助函数。
int rcOffsetPoly(const float* verts, int numVerts, float offset, float* outVerts, int maxOutVerts);

/// Applies the area id to all spans within the specified y-axis-aligned cylinder.
/// 将span id 应用到指定 y 轴对齐圆柱体内的所有span。
///
/// @see rcCompactHeightfield, rcMedianFilterWalkableArea
///
/// @ingroup recast
///
/// @param[in,out]	context				The build context to use during the operation.
/// @param[in]		position			The center of the base of the cylinder. [Form: (x, y, z)] [Units: wu]
/// @param[in]		radius				The radius of the cylinder. [Units: wu] [Limit: > 0]
/// @param[in]		height				The height of the cylinder. [Units: wu] [Limit: > 0]
/// @param[in]		areaId				The area id to apply. [Limit: <= #RC_WALKABLE_AREA]
/// @param[in,out]	compactHeightfield	A populated compact heightfield.
void rcMarkCylinderArea(rcContext* context, const float* position, float radius, float height,
						unsigned char areaId, rcCompactHeightfield& compactHeightfield);

/// Builds the distance field for the specified compact heightfield.
/// @ingroup recast
/// @param[in,out]	ctx		The build context to use during the operation.
/// @param[in,out]	chf		A populated compact heightfield.
/// @returns True if the operation completed successfully.
/// 构建指定紧凑高度场的距离场。
bool rcBuildDistanceField(rcContext* ctx, rcCompactHeightfield& chf);

/// Builds region data for the heightfield using watershed partitioning.
/// @ingroup recast
/// @param[in,out]	ctx				The build context to use during the operation.
/// @param[in,out]	chf				A populated compact heightfield.
/// @param[in]		borderSize		The size of the non-navigable border around the heightfield.
/// 								[Limit: >=0] [Units: vx]
/// @param[in]		minRegionArea	The minimum number of cells allowed to form isolated island areas.
/// 								[Limit: >=0] [Units: vx].
/// @param[in]		mergeRegionArea	Any regions with a span count smaller than this value will, if possible,
/// 								be merged with larger regions. [Limit: >=0] [Units: vx]
/// @returns True if the operation completed successfully.
///
/// 使用分水岭分区构建高度场的区域数据。
/// @param[in,out]	ctx				操作期间使用的构建上下文。
/// @param[in,out]	chf				紧凑高度场。
/// @param[in]		borderSize		高度场周围不可导航边框的大小。[Limit: >=0] [Units: vx]
/// @param[in]		minRegionArea	允许形成孤岛区域的最小单元数。[Limit: >=0] [Units: vx].
/// @param[in]		mergeRegionArea	如果可能，任何跨度计数小于此值的区域都将与更大的区域合并。 [限制：>=0] [单位：vx]
/// @returns True if the operation completed successfully.
bool rcBuildRegions(rcContext* ctx, rcCompactHeightfield& chf, int borderSize, int minRegionArea, int mergeRegionArea);

/// Builds region data for the heightfield by partitioning the heightfield in non-overlapping layers.
/// @ingroup recast
/// @param[in,out]	ctx				The build context to use during the operation.
/// @param[in,out]	chf				A populated compact heightfield.
/// @param[in]		borderSize		The size of the non-navigable border around the heightfield.
///  								[Limit: >=0] [Units: vx]
/// @param[in]		minRegionArea	The minimum number of cells allowed to form isolated island areas.
///  								[Limit: >=0] [Units: vx].
/// @returns True if the operation completed successfully.
/// 通过将高度场划分为非重叠层来构建高度场的区域数据。
bool rcBuildLayerRegions(rcContext* ctx, rcCompactHeightfield& chf, int borderSize, int minRegionArea);

/// Builds region data for the heightfield using simple monotone partitioning.
/// @ingroup recast
/// @param[in,out]	ctx				The build context to use during the operation.
/// @param[in,out]	chf				A populated compact heightfield.
/// @param[in]		borderSize		The size of the non-navigable border around the heightfield.
///  								[Limit: >=0] [Units: vx]
/// @param[in]		minRegionArea	The minimum number of cells allowed to form isolated island areas.
///  								[Limit: >=0] [Units: vx].
/// @param[in]		mergeRegionArea	Any regions with a span count smaller than this value will, if possible,
///  								be merged with larger regions. [Limit: >=0] [Units: vx]
/// @returns True if the operation completed successfully.
/// 使用简单的单调分区构建高度场的区域数据。
bool rcBuildRegionsMonotone(rcContext* ctx, rcCompactHeightfield& chf,
							int borderSize, int minRegionArea, int mergeRegionArea);

/// Sets the neighbor connection data for the specified direction.
/// @param[in]		span			The span to update.
/// @param[in]		direction		The direction to set. [Limits: 0 <= value < 4]
/// @param[in]		neighborIndex	The index of the neighbor span.
/// 设置指定方向的邻居连接数据。
/// @param[in,out]	span			The span to update.
/// @param[in]		direction		要设置的方向。 [限制：0 <= 值 < 4]
/// @param[in]		neighborIndex	邻居span的索引(1024)。
inline void rcSetCon(rcCompactSpan& span, int direction, int neighborIndex)
{
	// dir为 0左/1上/2右/3下，每个方向用6位存储，最高63层（6位都为1为不可行走,所以64-1=63）
    // 4个方向共24位，每6位存储一个方向，所有s.conn长度为24
	const unsigned int shift = (unsigned int)direction * 6;
	const unsigned int con = span.con;
	span.con = (con & ~(0x3f << shift)) | (((unsigned int)neighborIndex & 0x3f) << shift);
}

/// Gets neighbor connection data for the specified direction.
/// @param[in]		span		The span to check.
/// @param[in]		direction	The direction to check. [Limits: 0 <= value < 4]
/// @return The neighbor connection data for the specified direction, or #RC_NOT_CONNECTED if there is no connection.
/// 获取指定方向的邻居连接数据。
inline int rcGetCon(const rcCompactSpan& span, int direction)
{
	const unsigned int shift = (unsigned int)direction * 6;
	return (span.con >> shift) & 0x3f;
}

/// Gets the standard width (x-axis) offset for the specified direction.
/// @param[in]		direction		The direction. [Limits: 0 <= value < 4]
/// @return The width offset to apply to the current cell position to move in the direction.
/// 获取指定方向的标准宽度（x 轴）偏移量。
/// @return 应用于当前单元格位置以沿该方向移动的宽度偏移。
inline int rcGetDirOffsetX(int direction)
{
	static const int offset[4] = { -1, 0, 1, 0, };
	return offset[direction & 0x03];
}

// TODO (graham): Rename this to rcGetDirOffsetZ
/// Gets the standard height (z-axis) offset for the specified direction.
/// @param[in]		direction		The direction. [Limits: 0 <= value < 4]
/// @return The height offset to apply to the current cell position to move in the direction.
/// 获取指定方向的标准高度（z 轴）偏移量。
inline int rcGetDirOffsetY(int direction)
{
	static const int offset[4] = { 0, 1, 0, -1 };
	return offset[direction & 0x03];
}

/// Gets the direction for the specified offset. One of x and y should be 0.
/// @param[in]		offsetX		The x offset. [Limits: -1 <= value <= 1]
/// @param[in]		offsetZ		The z offset. [Limits: -1 <= value <= 1]
/// @return The direction that represents the offset.
/// 获取指定偏移的方向。 x 和 y 之一应为 0。
// x = 0,  y = -1 则 dirs 为 3 为 Z 轴负方向
// x = 0,  y = 1  则 dirs 为 1 为 Z 轴正方向
// x = -1, y = 0  则 dirs 为 0 为 X 轴负方向
// x = 1,  y = 0  则 dirs 为 2 为 X 轴正方向
inline int rcGetDirForOffset(int offsetX, int offsetZ)
{
	static const int dirs[5] = { 3, 0, -1, 2, 1 };
	return dirs[((offsetZ + 1) << 1) + offsetX];
}

/// @}
/// @name Layer, Contour, Polymesh, and Detail Mesh Functions
/// @see rcHeightfieldLayer, rcContourSet, rcPolyMesh, rcPolyMeshDetail
/// @{

/// Builds a layer set from the specified compact heightfield.
/// @ingroup recast
/// @param[in,out]	ctx				The build context to use during the operation.
/// @param[in]		chf				A fully built compact heightfield.
/// @param[in]		borderSize		The size of the non-navigable border around the heightfield. [Limit: >=0]
///  								[Units: vx]
/// @param[in]		walkableHeight	Minimum floor to 'ceiling' height that will still allow the floor area
///  								to be considered walkable. [Limit: >= 3] [Units: vx]
/// @param[out]		lset			The resulting layer set. (Must be pre-allocated.)
/// @returns True if the operation completed successfully.
/// 从指定的紧凑高度场构建图层集。
bool rcBuildHeightfieldLayers(rcContext* ctx, const rcCompactHeightfield& chf,
							  int borderSize, int walkableHeight,
							  rcHeightfieldLayerSet& lset);

/// Builds a contour set from the region outlines in the provided compact heightfield.
/// @ingroup recast
/// @param[in,out]	ctx			The build context to use during the operation.
/// @param[in]		chf			A fully built compact heightfield.
/// @param[in]		maxError	The maximum distance a simplified contour's border edges should deviate
/// 							the original raw contour. [Limit: >=0] [Units: wu]
/// @param[in]		maxEdgeLen	The maximum allowed length for contour edges along the border of the mesh.
/// 							[Limit: >=0] [Units: vx]
/// @param[out]		cset		The resulting contour set. (Must be pre-allocated.)
/// @param[in]		buildFlags	The build flags. (See: #rcBuildContoursFlags)
/// @returns True if the operation completed successfully.
/// 根据提供的紧凑高度场中的区域轮廓构建轮廓集。
bool rcBuildContours(rcContext* ctx, const rcCompactHeightfield& chf,
					 float maxError, int maxEdgeLen,
					 rcContourSet& cset, int buildFlags = RC_CONTOUR_TESS_WALL_EDGES);

/// Builds a polygon mesh from the provided contours.
/// @ingroup recast
/// @param[in,out]	ctx		The build context to use during the operation.
/// @param[in]		cset	A fully built contour set.
/// @param[in]		nvp		The maximum number of vertices allowed for polygons generated during the
/// 						contour to polygon conversion process. [Limit: >= 3]
/// @param[out]		mesh	The resulting polygon mesh. (Must be re-allocated.)
/// @returns True if the operation completed successfully.
/// 根据提供的轮廓构建多边形网格。
bool rcBuildPolyMesh(rcContext* ctx, const rcContourSet& cset, const int nvp, rcPolyMesh& mesh);

/// Merges multiple polygon meshes into a single mesh.
///  @ingroup recast
///  @param[in,out]	ctx		The build context to use during the operation.
///  @param[in]		meshes	An array of polygon meshes to merge. [Size: @p nmeshes]
///  @param[in]		nmeshes	The number of polygon meshes in the meshes array.
///  @param[in]		mesh	The resulting polygon mesh. (Must be pre-allocated.)
///  @returns True if the operation completed successfully.
/// 将多个多边形网格合并为单个网格。
bool rcMergePolyMeshes(rcContext* ctx, rcPolyMesh** meshes, const int nmeshes, rcPolyMesh& mesh);

/// Builds a detail mesh from the provided polygon mesh.
/// @ingroup recast
/// @param[in,out]	ctx				The build context to use during the operation.
/// @param[in]		mesh			A fully built polygon mesh.
/// @param[in]		chf				The compact heightfield used to build the polygon mesh.
/// @param[in]		sampleDist		Sets the distance to use when sampling the heightfield. [Limit: >=0] [Units: wu]
/// @param[in]		sampleMaxError	The maximum distance the detail mesh surface should deviate from
/// 								heightfield data. [Limit: >=0] [Units: wu]
/// @param[out]		dmesh			The resulting detail mesh.  (Must be pre-allocated.)
/// @returns True if the operation completed successfully.
/// 从提供的多边形网格构建细节网格。
bool rcBuildPolyMeshDetail(rcContext* ctx, const rcPolyMesh& mesh, const rcCompactHeightfield& chf,
						   float sampleDist, float sampleMaxError,
						   rcPolyMeshDetail& dmesh);

/// Copies the poly mesh data from src to dst.
/// @ingroup recast
/// @param[in,out]	ctx		The build context to use during the operation.
/// @param[in]		src		The source mesh to copy from.
/// @param[out]		dst		The resulting detail mesh. (Must be pre-allocated, must be empty mesh.)
/// @returns True if the operation completed successfully.
/// 将多边形网格数据从 src 复制到 dst。
bool rcCopyPolyMesh(rcContext* ctx, const rcPolyMesh& src, rcPolyMesh& dst);

/// Merges multiple detail meshes into a single detail mesh.
/// @ingroup recast
/// @param[in,out]	ctx		The build context to use during the operation.
/// @param[in]		meshes	An array of detail meshes to merge. [Size: @p nmeshes]
/// @param[in]		nmeshes	The number of detail meshes in the meshes array.
/// @param[out]		mesh	The resulting detail mesh. (Must be pre-allocated.)
/// @returns True if the operation completed successfully.
/// 将多个细节网格合并为一个细节网格。
bool rcMergePolyMeshDetails(rcContext* ctx, rcPolyMeshDetail** meshes, const int nmeshes, rcPolyMeshDetail& mesh);

/// @}

#endif // RECAST_H

///////////////////////////////////////////////////////////////////////////

// Due to the large amount of detail documentation for this file,
// the content normally located at the end of the header file has been separated
// out to a file in /Docs/Extern.
