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

#ifndef DETOURNODE_H
#define DETOURNODE_H

#include "DetourNavMesh.h"

enum dtNodeFlags
{
	DT_NODE_OPEN = 0x01,
	DT_NODE_CLOSED = 0x02,
	// parent of the node is not adjacent. Found using raycast.
	// 节点的父节点不相邻。 使用光线投射发现。
	DT_NODE_PARENT_DETACHED = 0x04
};

typedef unsigned short dtNodeIndex;
static const dtNodeIndex DT_NULL_IDX = (dtNodeIndex)~0;

static const int DT_NODE_PARENT_BITS = 24;
static const int DT_NODE_STATE_BITS = 2;
struct dtNode
{
	///< Position of the node.
	///< 节点的位置
	float pos[3];
	///< Cost from previous node to current node.
	///< 从前一个节点到当前节点的成本。
	float cost;
	///< Cost up to the node.
	///< 成本取决于节点。
	float total;
	///< Index to parent node.
	///< 父节点的索引。
	unsigned int pidx : DT_NODE_PARENT_BITS;
	///< extra state information. A polyRef can have multiple nodes with different extra info. see DT_MAX_STATES_PER_NODE
	///< 额外的状态信息。 PolyRef 可以有多个带有不同额外信息的节点。 请参阅 DT_MAX_STATES_PER_NODE
	unsigned int state : DT_NODE_STATE_BITS;
	///< Node flags. A combination of dtNodeFlags.
	///< 节点标志。 dtNodeFlags 的组合。
	unsigned int flags : 3;
	///< Polygon ref the node corresponds to.
	///< 节点对应的多边形参考。
	dtPolyRef id;
};

// number of extra states per node. See dtNode::state
// 每个节点的额外状态数。 请参阅 dtNode::state
static const int DT_MAX_STATES_PER_NODE = 1 << DT_NODE_STATE_BITS;

class dtNodePool
{
  public:
	dtNodePool(int maxNodes, int hashSize);
	~dtNodePool();
	void clear();

	// Get a dtNode by ref and extra state information. If there is none then - allocate
	// There can be more than one node for the same polyRef but with different extra state information
	// 通过 ref 和额外的状态信息获取 dtNode。 如果没有则 - 分配
	// 同一 PolyRef 可以有多个节点，但具有不同的额外状态信息
	dtNode* getNode(dtPolyRef id, unsigned char state=0);
	//根据id和state查找单个节点，返回节点指针
	dtNode* findNode(dtPolyRef id, unsigned char state);
	//根据id查找节点列表,返回找到的节点数
	unsigned int findNodes(dtPolyRef id, dtNode** nodes, const int maxNodes);

	inline unsigned int getNodeIdx(const dtNode* node) const
	{
		if (!node) return 0;
		return (unsigned int)(node - m_nodes) + 1;
	}

	inline dtNode* getNodeAtIdx(unsigned int idx)
	{
		if (!idx) return 0;
		return &m_nodes[idx - 1];
	}

	inline const dtNode* getNodeAtIdx(unsigned int idx) const
	{
		if (!idx) return 0;
		return &m_nodes[idx - 1];
	}

	inline int getMemUsed() const
	{
		return sizeof(*this) +
			sizeof(dtNode)*m_maxNodes +
			sizeof(dtNodeIndex)*m_maxNodes +
			sizeof(dtNodeIndex)*m_hashSize;
	}

	inline int getMaxNodes() const { return m_maxNodes; }

	inline int getHashSize() const { return m_hashSize; }
	inline dtNodeIndex getFirst(int bucket) const { return m_first[bucket]; }
	inline dtNodeIndex getNext(int i) const { return m_next[i]; }
	inline int getNodeCount() const { return m_nodeCount; }

private:
	// Explicitly disabled copy constructor and copy assignment operator.
	// 显式禁用复制构造函数和复制赋值运算符。
	dtNodePool(const dtNodePool&);
	dtNodePool& operator=(const dtNodePool&);

	dtNode* m_nodes;
	dtNodeIndex* m_first;
	dtNodeIndex* m_next;
	const int m_maxNodes;
	const int m_hashSize;
	int m_nodeCount;
};

//队列，其实是小顶堆
class dtNodeQueue
{
public:
	dtNodeQueue(int n);
	~dtNodeQueue();

	inline void clear() { m_size = 0; }

	inline dtNode* top() { return m_heap[0]; }

	inline dtNode* pop()
	{
		dtNode* result = m_heap[0];
		m_size--;
		trickleDown(0, m_heap[m_size]);
		return result;
	}

	inline void push(dtNode* node)
	{
		m_size++;
		bubbleUp(m_size-1, node);
	}

	inline void modify(dtNode* node)
	{
		for (int i = 0; i < m_size; ++i)
		{
			if (m_heap[i] == node)
			{
				bubbleUp(i, node);
				return;
			}
		}
	}

	inline bool empty() const { return m_size == 0; }

	inline int getMemUsed() const
	{
		return sizeof(*this) +
		sizeof(dtNode*) * (m_capacity + 1);
	}

	inline int getCapacity() const { return m_capacity; }

private:
	// Explicitly disabled copy constructor and copy assignment operator.
	// 显式禁用复制构造函数和复制赋值运算符。
	dtNodeQueue(const dtNodeQueue&);
	dtNodeQueue& operator=(const dtNodeQueue&);

	void bubbleUp(int i, dtNode* node);
	void trickleDown(int i, dtNode* node);

	dtNode** m_heap;
	const int m_capacity;
	int m_size;
};


#endif // DETOURNODE_H
