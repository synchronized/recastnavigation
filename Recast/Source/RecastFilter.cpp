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

#include "Recast.h"
#include "RecastAssert.h"

#include <stdlib.h>

namespace
{
	const int MAX_HEIGHTFIELD_HEIGHT = 0xffff; // TODO (graham): Move this to a more visible constant and update usages.
}

/// 如果不可通行span的最大值在其下方span的 @p walkableClimb 范围内，则将不可通行span标记为可通行。
void rcFilterLowHangingWalkableObstacles(rcContext* context, const int walkableClimb, rcHeightfield& heightfield)
{
	rcAssert(context);

	rcScopedTimer timer(context, RC_TIMER_FILTER_LOW_OBSTACLES);

	const int xSize = heightfield.width;
	const int zSize = heightfield.height;

	for (int z = 0; z < zSize; ++z)
	{
		for (int x = 0; x < xSize; ++x)
		{
			rcSpan* previousSpan = NULL;
			bool previousWasWalkable = false;
			unsigned char previousAreaID = RC_NULL_AREA;

			// For each span in the column...
			// 对于列中的每个跨度...
			for (rcSpan* span = heightfield.spans[x + z * xSize]; span != NULL; previousSpan = span, span = span->next)
			{
				const bool walkable = span->area != RC_NULL_AREA;

				// If current span is not walkable, but there is walkable span just below it and the height difference
				// is small enough for the agent to walk over, mark the current span as walkable too.
				// 如果当前跨度不可步行，但其下方有可步行span，并且高度差足够小，可以让Agent走过，则也将当前span标记为可步行。
				if (!walkable && previousWasWalkable && (int)span->smax - (int)previousSpan->smax <= walkableClimb)
				{
					span->area = previousAreaID;
				}

				// Copy the original walkable value regardless of whether we changed it.
				// This prevents multiple consecutive non-walkable spans from being erroneously marked as walkable.
				// 复制原始的步行值，无论我们是否更改它。
				// 这可以防止多个连续的不可步行span被错误地标记为可步行。
				previousWasWalkable = walkable;
				previousAreaID = span->area;
			}
		}
	}
}

/// 将ledges的span标记为不可行走。
/// 过滤有效区间和陡峭区间
void rcFilterLedgeSpans(rcContext* context, const int walkableHeight, const int walkableClimb, rcHeightfield& heightfield)
{
	rcAssert(context);

	rcScopedTimer timer(context, RC_TIMER_FILTER_BORDER);

	const int xSize = heightfield.width;
	const int zSize = heightfield.height;

	// Mark spans that are adjacent to a ledge as unwalkable..
	// 将边界span标记为不可通行。
	for (int z = 0; z < zSize; ++z)
	{
		for (int x = 0; x < xSize; ++x)
		{
			for (rcSpan* span = heightfield.spans[x + z * xSize]; span; span = span->next)
			{
				// Skip non-walkable spans.
				if (span->area == RC_NULL_AREA)
				{
					continue;
				}

				const int floor = (int)(span->smax);
				const int ceiling = span->next ? (int)(span->next->smin) : MAX_HEIGHTFIELD_HEIGHT;

				// The difference between this walkable area and the lowest neighbor walkable area.
				// This is the difference between the current span and all neighbor spans that have
				// enough space for an agent to move between, but not accounting at all for surface slope.
				// 该步行区域与最低相邻步行区域之间的差异。
				// 这是当前span和所有具有足够空间供Agent移动的相邻span之间的差异，但根本不考虑表面坡度。
				int lowestNeighborFloorDifference = MAX_HEIGHTFIELD_HEIGHT;

				// Min and max height of accessible neighbours.
				int lowestTraversableNeighborFloor = span->smax;
				int highestTraversableNeighborFloor = span->smax;

				for (int direction = 0; direction < 4; ++direction)
				{
					const int neighborX = x + rcGetDirOffsetX(direction);
					const int neighborZ = z + rcGetDirOffsetY(direction);

					// Skip neighbours which are out of bounds.
					// 跳过超出范围的邻居。
					if (neighborX < 0 || neighborZ < 0 || neighborX >= xSize || neighborZ >= zSize)
					{
						lowestNeighborFloorDifference = -walkableClimb - 1;
						break;
					}

					const rcSpan* neighborSpan = heightfield.spans[neighborX + neighborZ * xSize];

					// The most we can step down to the neighbor is the walkableClimb distance.
					// Start with the area under the neighbor span
					// 我们可以下到邻居的最大距离是walkableClimb 距离。
					// 从邻居跨度下的区域开始
					int neighborCeiling = neighborSpan ? (int)neighborSpan->smin : MAX_HEIGHTFIELD_HEIGHT;

					// Skip neighbour if the gap between the spans is too small.
					// 如果span之间的间隙太小，则跳过邻居。
					if (rcMin(ceiling, neighborCeiling) - floor >= walkableHeight)
					{
						lowestNeighborFloorDifference = (-walkableClimb - 1);
						break;
					}

					// For each span in the neighboring column...
					//其余的span
					for (; neighborSpan != NULL; neighborSpan = neighborSpan->next)
					{
						const int neighborFloor = (int)neighborSpan->smax;
						neighborCeiling = neighborSpan->next ? (int)neighborSpan->next->smin : MAX_HEIGHTFIELD_HEIGHT;

						// Only consider neighboring areas that have enough overlap to be potentially traversable.
						// 仅考虑具有足够重叠以便可以穿越的相邻区域。
						// 如果Span之间的共通高度(指在高度上双方都是空心的高度)太小，则跳过 neightbour
						if (rcMin(ceiling, neighborCeiling) - rcMax(floor, neighborFloor) < walkableHeight)
						{
							// No space to traverse between them.
							// 它们之间没有空间可以穿越。
							continue;
						}

						const int neighborFloorDifference = neighborFloor - floor;
						lowestNeighborFloorDifference = rcMin(lowestNeighborFloorDifference, neighborFloorDifference);

						// Find min/max accessible neighbor height.
						// Only consider neighbors that are at most walkableClimb away.
						// 查找最小/最大可访问邻居高度。
						// 只考虑最多可步行爬走的邻居。
						// 查找最小/最大可访问邻居高度
						if (rcAbs(neighborFloorDifference) <= walkableClimb)
						{
							// There is space to move to the neighbor cell and the slope isn't too much.
							// 有空间可以移动到相邻的单元格，而且坡度也不是太大。
							lowestTraversableNeighborFloor = rcMin(lowestTraversableNeighborFloor, neighborFloor);
							highestTraversableNeighborFloor = rcMax(highestTraversableNeighborFloor, neighborFloor);
						}
						else if (neighborFloorDifference < -walkableClimb)
						{
							// We already know this will be considered a ledge span so we can early-out
							// 我们已经知道这将被视为一个边界 span，因此我们可以提前退出
							break;
						}
					}
				}

				// The current span is close to a ledge if the magnitude of the drop to any neighbour span is greater than the walkableClimb distance.
				// That is, there is a gap that is large enough to let an agent move between them, but the drop (surface slope) is too large to allow it.
				// (If this is the case, then biggestNeighborStepDown will be negative, so compare against the negative walkableClimb as a means of checking
				// the magnitude of the delta)
				// 如果任何相邻span的下降幅度大于 walkableClimb 距离，则当前span接近ledge。
				// 也就是说，间隙足够大，足以让Agent在它们之间移动，但落差（表面坡度）太大而不允许。
				// （如果是这种情况，那么 MaximumNeighborStepDown 将为负数，因此与负的 walkableClimb 进行比较，
				// 作为检查增量大小的方法）
				// 如果下降到任何邻居Span小于 walkableClimb，将Span标记为RC_NULL_AREA
				if (lowestNeighborFloorDifference < -walkableClimb)
				{
					span->area = RC_NULL_AREA;
				}
				// If the difference between all neighbor floors is too large, this is a steep slope, so mark the span as an unwalkable ledge.
				// 如果所有相邻楼层之间的差异太大，则这是一个陡坡，因此将span标记为无法行走的边界(ledge)。
				else if (highestTraversableNeighborFloor - lowestTraversableNeighborFloor > walkableClimb)
				{
					span->area = RC_NULL_AREA;
				}
			}
		}
	}
}

/// 如果span上方的间隙小于指定的通行高度，则将可通行span标记为不可通行。
void rcFilterWalkableLowHeightSpans(rcContext* context, const int walkableHeight, rcHeightfield& heightfield)
{
	rcAssert(context);
	rcScopedTimer timer(context, RC_TIMER_FILTER_WALKABLE);

	const int xSize = heightfield.width;
	const int zSize = heightfield.height;

	// Remove walkable flag from spans which do not have enough
	// space above them for the agent to stand there.
	// 从span上移除可通行标志，因为span上方没有足够的空间供Agent站在那里。
	for (int z = 0; z < zSize; ++z)
	{
		for (int x = 0; x < xSize; ++x)
		{
			for (rcSpan* span = heightfield.spans[x + z*xSize]; span; span = span->next)
			{
				const int floor = (int)(span->smax);
				const int ceiling = span->next ? (int)(span->next->smin) : MAX_HEIGHTFIELD_HEIGHT;
				if (ceiling - floor < walkableHeight)
				{
					span->area = RC_NULL_AREA;
				}
			}
		}
	}
}
