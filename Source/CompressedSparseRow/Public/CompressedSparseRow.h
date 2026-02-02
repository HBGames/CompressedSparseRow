// Copyright Hitbox Games, LLC. All Rights Reserved.

#pragma once

#include "Algo/Sort.h"
#include "Containers/Array.h"
#include "Containers/ArrayView.h"
#include "Templates/IsIntegral.h"
#include "Templates/IsSigned.h"

/**
 * Compressed Sparse Row (CSR) graph adjacency container.
 *
 * CSR is a memory-efficient representation for sparse graphs that stores all neighbor
 * indices in a single contiguous array, with a separate offset array indicating where
 * each node's neighbors begin. This layout provides excellent cache locality for
 * neighbor iteration and minimal memory overhead compared to adjacency list approaches.
 *
 * Memory Layout:
 *   Offsets array:   [0, 2, 5, 5, 8]  (length = NumNodes + 1)
 *   Neighbors array: [1, 3, 0, 2, 4, 0, 1, 2]  (length = NumEdges)
 *
 *   Node 0's neighbors: Neighbors[Offsets[0]..Offsets[1]) = Neighbors[0..2) = {1, 3}
 *   Node 1's neighbors: Neighbors[Offsets[1]..Offsets[2]) = Neighbors[2..5) = {0, 2, 4}
 *   Node 2's neighbors: Neighbors[Offsets[2]..Offsets[3]) = Neighbors[5..5) = {} (no neighbors)
 *   Node 3's neighbors: Neighbors[Offsets[3]..Offsets[4]) = Neighbors[5..8) = {0, 1, 2}
 *
 * Performance Characteristics:
 *   - Neighbor iteration: O(degree) with optimal cache behavior
 *   - Memory usage: O(NumNodes + NumEdges) integers
 *   - Build time: O(NumEdges) without sorting, O(NumEdges * log(MaxDegree)) with sorting
 *   - Not designed for dynamic edge insertion/removal after Build()
 *
 * Template Parameters:
 * @tparam IndexType  Integer type for node indices. Must be signed integral (int32 recommended).
 *                    Determines maximum node count. Using int32 supports up to ~2 billion nodes.
 *
 * Usage:
 * @code
 *     TCSRGraph<int32> Graph;
 *     TArray<TCSRGraph<int32>::FEdge> Edges;
 *     Edges.Add({0, 1});  // Edge from node 0 to node 1
 *     Edges.Add({0, 2});
 *     Edges.Add({1, 2});
 *     Graph.Build(3, Edges, true, true);  // 3 nodes, sorted, deduplicated
 *
 *     for (int32 Neighbor : Graph.GetNeighbors(0))
 *     {
 *         // Process each neighbor of node 0
 *     }
 * @endcode
 *
 * @note Offsets use int32 regardless of IndexType, limiting total edges to MAX_int32.
 * @note Build() should be called infrequently (e.g., on topology changes). Iteration is the hot path.
 * @see TCSRGraphWithEdgeData for graphs requiring per-edge payload data.
 */
template <typename IndexType = int32>
class TCSRGraph
{
	static_assert(TIsIntegral<IndexType>::Value, "IndexType must be an integral type.");
	static_assert(TIsSigned<IndexType>::Value, "IndexType must be signed to support validity checks.");

public:
	/**
	 * Represents a directed edge in the graph.
	 * Used as input to Build() to define graph topology.
	 */
	struct FEdge
	{
		/** Source node index. Must be in range [0, NumNodes). */
		IndexType From = 0;

		/** Target node index. Must be in range [0, NumNodes). */
		IndexType To = 0;
	};

public:
	TCSRGraph() = default;

	/**
	 * Resets the graph to empty state, releasing memory.
	 * Equivalent to constructing a new graph.
	 */
	void Reset()
	{
		NumNodes = 0;
		Offsets.Reset();
		Neighbors.Reset();
	}

	/**
	 * Clears the graph but retains allocated memory for reuse.
	 * @param Slack  Optional slack to retain in internal arrays.
	 */
	void Empty(int32 Slack = 0)
	{
		NumNodes = 0;
		Offsets.Empty(Slack);
		Neighbors.Empty(Slack);
	}

	/** @return The number of nodes in the graph. */
	[[nodiscard]] FORCEINLINE IndexType GetNumNodes() const { return NumNodes; }

	/** @return The total number of directed edges in the graph. */
	[[nodiscard]] FORCEINLINE int32 GetNumEdges() const { return Neighbors.Num(); }

	/**
	 * Checks if a node index is valid for this graph.
	 * @param Node  The node index to validate.
	 * @return True if Node is in the range [0, NumNodes).
	 */
	[[nodiscard]] FORCEINLINE bool IsValidNode(IndexType Node) const
	{
		return Node >= 0 && Node < NumNodes;
	}

	/**
	 * Returns the out-degree (number of outgoing edges) for a node.
	 * @param Node  The node whose degree to retrieve. Must be a valid node index.
	 * @return The number of outgoing edges from this node.
	 */
	[[nodiscard]] FORCEINLINE int32 GetDegree(IndexType Node) const
	{
		check(IsValidNode(Node));
		return Offsets[static_cast<int32>(Node) + 1] - Offsets[static_cast<int32>(Node)];
	}

	/**
	 * Returns a view of all neighbors for a given node.
	 *
	 * The returned view is valid until the graph is modified via Build(), Reset(), or Empty().
	 * Iteration order depends on whether bSortNeighborsPerNode was set during Build().
	 *
	 * @param Node  The node whose neighbors to retrieve. Must be a valid node index.
	 * @return A const array view of neighbor indices. May be empty if the node has no outgoing edges.
	 */
	[[nodiscard]] FORCEINLINE TConstArrayView<IndexType> GetNeighbors(IndexType Node) const
	{
		check(IsValidNode(Node));
		const int32 Start = Offsets[static_cast<int32>(Node)];
		const int32 End = Offsets[static_cast<int32>(Node) + 1];
		return TConstArrayView<IndexType>(Neighbors.GetData() + Start, End - Start);
	}

	/**
	 * Returns the raw offsets array for direct access.
	 *
	 * The offsets array has length NumNodes + 1. For node i, its neighbors
	 * span indices [Offsets[i], Offsets[i+1]) in the neighbors array.
	 *
	 * @return Const view of the offsets array.
	 */
	[[nodiscard]] FORCEINLINE TConstArrayView<int32> GetOffsets() const
	{
		return TConstArrayView<int32>(Offsets);
	}

	/**
	 * Returns the raw neighbor indices array for direct access.
	 *
	 * Use with GetOffsets() for bulk processing or serialization.
	 * Neighbors for node i are at indices [Offsets[i], Offsets[i+1]).
	 *
	 * @return Const view of the neighbors array.
	 */
	[[nodiscard]] FORCEINLINE TConstArrayView<IndexType> GetNeighborArray() const
	{
		return TConstArrayView<IndexType>(Neighbors);
	}

	/**
	 * Builds the CSR graph from an edge list.
	 *
	 * This operation reconstructs the entire graph data structure. It should be called
	 * when the graph topology changes, not for frequent updates. The build process:
	 *   1. Counts the out-degree of each node
	 *   2. Computes prefix sums to determine offsets
	 *   3. Scatters edges into the neighbors array
	 *   4. Optionally sorts and/or deduplicates per-node neighbor lists
	 *
	 * @param InNumNodes              Total number of nodes. Nodes are indexed [0, InNumNodes).
	 * @param Edges                   Array of directed edges defining the graph topology.
	 * @param bSortNeighborsPerNode   If true, sorts each node's neighbors by index for deterministic iteration.
	 * @param bRemoveDuplicateEdges   If true, removes duplicate (From, To) pairs, keeping the first occurrence.
	 */
	void Build(
		IndexType InNumNodes,
		TConstArrayView<FEdge> Edges,
		bool bSortNeighborsPerNode = false,
		bool bRemoveDuplicateEdges = false)
	{
		check(InNumNodes >= 0);
		NumNodes = InNumNodes;

		// Offsets array uses int32, so total edge count must not exceed int32 max.
		checkf(Edges.Num() <= MAX_int32, TEXT("TCSRGraph::Build - Edge count exceeds int32 limit."));

		Offsets.SetNumZeroed(static_cast<int32>(NumNodes) + 1);

		// Pass 1: Count out-degree for each node by incrementing Offsets[From + 1].
		// This positions counts one slot ahead to prepare for prefix sum.
		for (const FEdge& Edge : Edges)
		{
			checkf(Edge.From >= 0 && Edge.From < NumNodes, TEXT("TCSRGraph::Build - Edge.From index out of range."));
			checkf(Edge.To >= 0 && Edge.To < NumNodes, TEXT("TCSRGraph::Build - Edge.To index out of range."));
			Offsets[static_cast<int32>(Edge.From) + 1]++;
		}

		// Pass 2: Convert degree counts to cumulative offsets via prefix sum.
		// After this, Offsets[i] = starting index for node i's neighbors.
		int64 RunningSum = 0;
		for (int32 i = 1; i < Offsets.Num(); ++i)
		{
			RunningSum += Offsets[i];
			checkf(RunningSum <= MAX_int32, TEXT("TCSRGraph::Build - Cumulative edge count overflows int32."));
			Offsets[i] = static_cast<int32>(RunningSum);
		}

		// Pass 3: Scatter edge targets into the neighbors array.
		// Use a cursor array to track the next write position for each node.
		Neighbors.SetNumUninitialized(Edges.Num());

		TArray<int32> WriteCursors;
		WriteCursors.SetNumUninitialized(Offsets.Num());
		FMemory::Memcpy(WriteCursors.GetData(), Offsets.GetData(), Offsets.Num() * sizeof(int32));

		for (const FEdge& Edge : Edges)
		{
			const int32 WriteIndex = WriteCursors[static_cast<int32>(Edge.From)]++;
			Neighbors[WriteIndex] = Edge.To;
		}

		// Pass 4: Optional post-processing for sorting and deduplication.
		if (bSortNeighborsPerNode || bRemoveDuplicateEdges)
		{
			SortAndOrDedup(bRemoveDuplicateEdges);
		}
	}

private:
	/**
	 * Sorts neighbor lists and optionally removes duplicate edges.
	 *
	 * Sorting is performed via an indirect permutation array to maintain stable ordering
	 * (earlier edges are kept when duplicates exist). The algorithm:
	 *   1. For each node, build a permutation array of local indices
	 *   2. Sort permutation by (neighbor index, original position) for stability
	 *   3. Emit neighbors in sorted order, skipping duplicates if requested
	 *
	 * @param bRemoveDuplicateEdges  If true, consecutive identical neighbors are collapsed to one.
	 */
	void SortAndOrDedup(bool bRemoveDuplicateEdges)
	{
		TArray<int32> NewOffsets;
		TArray<IndexType> NewNeighbors;
		NewOffsets.SetNumZeroed(static_cast<int32>(NumNodes) + 1);
		NewNeighbors.Reserve(Neighbors.Num());

		// Temporary buffers reused across nodes to minimize allocations.
		TArray<IndexType> NodeNeighbors;
		TArray<int32> SortPermutation;

		for (IndexType Node = 0; Node < NumNodes; ++Node)
		{
			const int32 RangeStart = Offsets[static_cast<int32>(Node)];
			const int32 RangeEnd = Offsets[static_cast<int32>(Node) + 1];
			const int32 NeighborCount = RangeEnd - RangeStart;

			if (NeighborCount == 0)
			{
				NewOffsets[static_cast<int32>(Node) + 1] = NewNeighbors.Num();
				continue;
			}

			// Copy this node's neighbors to temporary buffer.
			NodeNeighbors.SetNumUninitialized(NeighborCount);
			for (int32 i = 0; i < NeighborCount; ++i)
			{
				NodeNeighbors[i] = Neighbors[RangeStart + i];
			}

			// Build permutation array [0, 1, 2, ..., NeighborCount-1].
			SortPermutation.SetNumUninitialized(NeighborCount);
			for (int32 i = 0; i < NeighborCount; ++i)
			{
				SortPermutation[i] = i;
			}

			// Sort permutation by (neighbor value, original index) for stable ordering.
			// Using original index as tiebreaker ensures "first edge wins" semantics for dedup.
			Algo::Sort(SortPermutation, [&NodeNeighbors](int32 A, int32 B)
			{
				const IndexType NeighborA = NodeNeighbors[A];
				const IndexType NeighborB = NodeNeighbors[B];
				if (NeighborA != NeighborB)
				{
					return NeighborA < NeighborB;
				}
				return A < B;
			});

			// Emit sorted neighbors, optionally skipping duplicates.
			IndexType PreviousNeighbor = IndexType(-1);
			bool bHasPrevious = false;

			for (int32 PermIdx = 0; PermIdx < NeighborCount; ++PermIdx)
			{
				const IndexType CurrentNeighbor = NodeNeighbors[SortPermutation[PermIdx]];

				if (bRemoveDuplicateEdges && bHasPrevious && CurrentNeighbor == PreviousNeighbor)
				{
					continue;
				}

				PreviousNeighbor = CurrentNeighbor;
				bHasPrevious = true;
				NewNeighbors.Add(CurrentNeighbor);
			}

			NewOffsets[static_cast<int32>(Node) + 1] = NewNeighbors.Num();
		}

		Offsets = MoveTemp(NewOffsets);
		Neighbors = MoveTemp(NewNeighbors);
	}

private:
	/** Number of nodes in the graph. Nodes are indexed [0, NumNodes). */
	IndexType NumNodes = 0;

	/**
	 * Offset array defining neighbor ranges.
	 * Length is NumNodes + 1. Node i's neighbors are at Neighbors[Offsets[i]..Offsets[i+1]).
	 */
	TArray<int32> Offsets;

	/** Contiguous storage of all neighbor indices, partitioned by source node. */
	TArray<IndexType> Neighbors;
};


/**
 * Compressed Sparse Row (CSR) graph with per-edge payload data.
 *
 * Extends TCSRGraph with an additional parallel array storing arbitrary data for each edge.
 * The EdgeData array is aligned with the Neighbors array: EdgeData[i] corresponds to the
 * edge targeting Neighbors[i].
 *
 * This is useful for weighted graphs, edge labels, edge capacities, or any per-edge metadata.
 *
 * Memory Layout:
 *   Offsets array:   [0, 2, 4]  (length = NumNodes + 1)
 *   Neighbors array: [1, 2, 0, 2]  (length = NumEdges)
 *   EdgeData array:  [10, 20, 15, 5]  (length = NumEdges, parallel to Neighbors)
 *
 *   Node 0 -> Node 1 with data 10
 *   Node 0 -> Node 2 with data 20
 *   Node 1 -> Node 0 with data 15
 *   Node 1 -> Node 2 with data 5
 *
 * Template Parameters:
 * @tparam IndexType     Integer type for node indices. Must be signed integral (int32 recommended).
 * @tparam EdgeDataType  Type of per-edge payload data. Can be any copyable type (arithmetic, struct, etc.).
 *
 * Usage:
 * @code
 *     // Weighted graph example
 *     TCSRGraphWithEdgeData<int32, float> WeightedGraph;
 *     TArray<TCSRGraphWithEdgeData<int32, float>::FEdge> Edges;
 *     Edges.Add({0, 1, 1.5f});  // Edge 0->1 with weight 1.5
 *     Edges.Add({0, 2, 2.0f});  // Edge 0->2 with weight 2.0
 *     WeightedGraph.Build(3, Edges);
 *
 *     TConstArrayView<int32> Neighbors = WeightedGraph.GetNeighbors(0);
 *     TConstArrayView<float> Weights = WeightedGraph.GetEdgeData(0);
 *     for (int32 i = 0; i < Neighbors.Num(); ++i)
 *     {
 *         ProcessEdge(Neighbors[i], Weights[i]);
 *     }
 * @endcode
 *
 * @see TCSRGraph for graphs without edge data.
 */
template <typename IndexType = int32, typename EdgeDataType = uint16>
class TCSRGraphWithEdgeData
{
	static_assert(TIsIntegral<IndexType>::Value, "IndexType must be an integral type.");
	static_assert(TIsSigned<IndexType>::Value, "IndexType must be signed to support validity checks.");

public:
	/**
	 * Represents a directed edge with associated payload data.
	 * Used as input to Build() to define graph topology and edge data.
	 */
	struct FEdge
	{
		/** Source node index. Must be in range [0, NumNodes). */
		IndexType From = 0;

		/** Target node index. Must be in range [0, NumNodes). */
		IndexType To = 0;

		/** Arbitrary payload data associated with this edge. */
		EdgeDataType Data{};
	};

public:
	TCSRGraphWithEdgeData() = default;

	/**
	 * Resets the graph to empty state, releasing memory.
	 * Equivalent to constructing a new graph.
	 */
	void Reset()
	{
		NumNodes = 0;
		Offsets.Reset();
		Neighbors.Reset();
		EdgeData.Reset();
	}

	/**
	 * Clears the graph but retains allocated memory for reuse.
	 * @param Slack  Optional slack to retain in internal arrays.
	 */
	void Empty(int32 Slack = 0)
	{
		NumNodes = 0;
		Offsets.Empty(Slack);
		Neighbors.Empty(Slack);
		EdgeData.Empty(Slack);
	}

	/** @return The number of nodes in the graph. */
	[[nodiscard]] FORCEINLINE IndexType GetNumNodes() const { return NumNodes; }

	/** @return The total number of directed edges in the graph. */
	[[nodiscard]] FORCEINLINE int32 GetNumEdges() const { return Neighbors.Num(); }

	/**
	 * Checks if a node index is valid for this graph.
	 * @param Node  The node index to validate.
	 * @return True if Node is in the range [0, NumNodes).
	 */
	[[nodiscard]] FORCEINLINE bool IsValidNode(IndexType Node) const
	{
		return Node >= 0 && Node < NumNodes;
	}

	/**
	 * Returns the out-degree (number of outgoing edges) for a node.
	 * @param Node  The node whose degree to retrieve. Must be a valid node index.
	 * @return The number of outgoing edges from this node.
	 */
	[[nodiscard]] FORCEINLINE int32 GetDegree(IndexType Node) const
	{
		check(IsValidNode(Node));
		return Offsets[static_cast<int32>(Node) + 1] - Offsets[static_cast<int32>(Node)];
	}

	/**
	 * Returns a view of all neighbors for a given node.
	 *
	 * The returned view is parallel to GetEdgeData() for the same node:
	 * GetNeighbors(N)[i] and GetEdgeData(N)[i] describe the same edge.
	 *
	 * @param Node  The node whose neighbors to retrieve. Must be a valid node index.
	 * @return A const array view of neighbor indices. May be empty if the node has no outgoing edges.
	 */
	[[nodiscard]] FORCEINLINE TConstArrayView<IndexType> GetNeighbors(IndexType Node) const
	{
		check(IsValidNode(Node));
		const int32 Start = Offsets[static_cast<int32>(Node)];
		const int32 End = Offsets[static_cast<int32>(Node) + 1];
		return TConstArrayView<IndexType>(Neighbors.GetData() + Start, End - Start);
	}

	/**
	 * Returns a view of edge data for all outgoing edges from a node.
	 *
	 * The returned view is parallel to GetNeighbors() for the same node:
	 * GetEdgeData(N)[i] is the payload for the edge to GetNeighbors(N)[i].
	 *
	 * @param Node  The node whose edge data to retrieve. Must be a valid node index.
	 * @return A const array view of edge payloads. May be empty if the node has no outgoing edges.
	 */
	[[nodiscard]] FORCEINLINE TConstArrayView<EdgeDataType> GetEdgeData(IndexType Node) const
	{
		check(IsValidNode(Node));
		const int32 Start = Offsets[static_cast<int32>(Node)];
		const int32 End = Offsets[static_cast<int32>(Node) + 1];
		return TConstArrayView<EdgeDataType>(EdgeData.GetData() + Start, End - Start);
	}

	/** @return Const view of the offsets array (length = NumNodes + 1). */
	[[nodiscard]] FORCEINLINE TConstArrayView<int32> GetOffsets() const { return TConstArrayView<int32>(Offsets); }

	/** @return Const view of the raw neighbor indices array. */
	[[nodiscard]] FORCEINLINE TConstArrayView<IndexType> GetNeighborArray() const { return TConstArrayView<IndexType>(Neighbors); }

	/** @return Const view of the raw edge data array (parallel to neighbor array). */
	[[nodiscard]] FORCEINLINE TConstArrayView<EdgeDataType> GetEdgeDataArray() const { return TConstArrayView<EdgeDataType>(EdgeData); }

	/**
	 * Builds the CSR graph from an edge list with payloads.
	 *
	 * This operation reconstructs the entire graph data structure. Edge data is preserved
	 * in alignment with the neighbor array through all transformations including sorting
	 * and deduplication.
	 *
	 * @param InNumNodes              Total number of nodes. Nodes are indexed [0, InNumNodes).
	 * @param Edges                   Array of directed edges with payloads defining the graph.
	 * @param bSortNeighborsPerNode   If true, sorts each node's neighbors by index for deterministic iteration.
	 * @param bRemoveDuplicateEdges   If true, removes duplicate (From, To) pairs, keeping the first occurrence and its payload.
	 */
	void Build(
		IndexType InNumNodes,
		TConstArrayView<FEdge> Edges,
		bool bSortNeighborsPerNode = false,
		bool bRemoveDuplicateEdges = false)
	{
		check(InNumNodes >= 0);
		NumNodes = InNumNodes;

		checkf(Edges.Num() <= MAX_int32, TEXT("TCSRGraphWithEdgeData::Build - Edge count exceeds int32 limit."));

		Offsets.SetNumZeroed(static_cast<int32>(NumNodes) + 1);

		// Pass 1: Count out-degree for each node.
		for (const FEdge& Edge : Edges)
		{
			checkf(Edge.From >= 0 && Edge.From < NumNodes, TEXT("TCSRGraphWithEdgeData::Build - Edge.From index out of range."));
			checkf(Edge.To >= 0 && Edge.To < NumNodes, TEXT("TCSRGraphWithEdgeData::Build - Edge.To index out of range."));
			Offsets[static_cast<int32>(Edge.From) + 1]++;
		}

		// Pass 2: Convert degree counts to cumulative offsets via prefix sum.
		int64 RunningSum = 0;
		for (int32 i = 1; i < Offsets.Num(); ++i)
		{
			RunningSum += Offsets[i];
			checkf(RunningSum <= MAX_int32, TEXT("TCSRGraphWithEdgeData::Build - Cumulative edge count overflows int32."));
			Offsets[i] = static_cast<int32>(RunningSum);
		}

		// Pass 3: Scatter edge targets and payloads into respective arrays.
		Neighbors.SetNumUninitialized(Edges.Num());
		EdgeData.SetNumUninitialized(Edges.Num());

		TArray<int32> WriteCursors;
		WriteCursors.SetNumUninitialized(Offsets.Num());
		FMemory::Memcpy(WriteCursors.GetData(), Offsets.GetData(), Offsets.Num() * sizeof(int32));

		for (const FEdge& Edge : Edges)
		{
			const int32 WriteIndex = WriteCursors[static_cast<int32>(Edge.From)]++;
			Neighbors[WriteIndex] = Edge.To;
			EdgeData[WriteIndex] = Edge.Data;
		}

		// Pass 4: Optional post-processing for sorting and deduplication.
		if (bSortNeighborsPerNode || bRemoveDuplicateEdges)
		{
			SortAndOrDedup(bRemoveDuplicateEdges);
		}
	}

private:
	/**
	 * Sorts neighbor lists and optionally removes duplicate edges, preserving edge data alignment.
	 *
	 * @param bRemoveDuplicateEdges  If true, consecutive identical neighbors are collapsed, keeping first payload.
	 */
	void SortAndOrDedup(bool bRemoveDuplicateEdges)
	{
		TArray<int32> NewOffsets;
		TArray<IndexType> NewNeighbors;
		TArray<EdgeDataType> NewEdgeData;

		NewOffsets.SetNumZeroed(static_cast<int32>(NumNodes) + 1);
		NewNeighbors.Reserve(Neighbors.Num());
		NewEdgeData.Reserve(EdgeData.Num());

		// Temporary buffers reused across nodes to minimize allocations.
		TArray<IndexType> NodeNeighbors;
		TArray<EdgeDataType> NodeEdgeData;
		TArray<int32> SortPermutation;

		for (IndexType Node = 0; Node < NumNodes; ++Node)
		{
			const int32 RangeStart = Offsets[static_cast<int32>(Node)];
			const int32 RangeEnd = Offsets[static_cast<int32>(Node) + 1];
			const int32 NeighborCount = RangeEnd - RangeStart;

			if (NeighborCount == 0)
			{
				NewOffsets[static_cast<int32>(Node) + 1] = NewNeighbors.Num();
				continue;
			}

			// Copy this node's data to temporary buffers.
			NodeNeighbors.SetNumUninitialized(NeighborCount);
			NodeEdgeData.SetNumUninitialized(NeighborCount);
			for (int32 i = 0; i < NeighborCount; ++i)
			{
				NodeNeighbors[i] = Neighbors[RangeStart + i];
				NodeEdgeData[i] = EdgeData[RangeStart + i];
			}

			// Build and sort permutation array.
			SortPermutation.SetNumUninitialized(NeighborCount);
			for (int32 i = 0; i < NeighborCount; ++i)
			{
				SortPermutation[i] = i;
			}

			Algo::Sort(SortPermutation, [&NodeNeighbors](int32 A, int32 B)
			{
				const IndexType NeighborA = NodeNeighbors[A];
				const IndexType NeighborB = NodeNeighbors[B];
				if (NeighborA != NeighborB)
				{
					return NeighborA < NeighborB;
				}
				return A < B;
			});

			// Emit sorted data, optionally skipping duplicates.
			IndexType PreviousNeighbor = IndexType(-1);
			bool bHasPrevious = false;

			for (int32 PermIdx = 0; PermIdx < NeighborCount; ++PermIdx)
			{
				const int32 OriginalIdx = SortPermutation[PermIdx];
				const IndexType CurrentNeighbor = NodeNeighbors[OriginalIdx];

				if (bRemoveDuplicateEdges && bHasPrevious && CurrentNeighbor == PreviousNeighbor)
				{
					continue;
				}

				PreviousNeighbor = CurrentNeighbor;
				bHasPrevious = true;
				NewNeighbors.Add(CurrentNeighbor);
				NewEdgeData.Add(NodeEdgeData[OriginalIdx]);
			}

			NewOffsets[static_cast<int32>(Node) + 1] = NewNeighbors.Num();
		}

		Offsets = MoveTemp(NewOffsets);
		Neighbors = MoveTemp(NewNeighbors);
		EdgeData = MoveTemp(NewEdgeData);
	}

private:
	/** Number of nodes in the graph. Nodes are indexed [0, NumNodes). */
	IndexType NumNodes = 0;

	/** Offset array defining neighbor ranges. Length is NumNodes + 1. */
	TArray<int32> Offsets;

	/** Contiguous storage of all neighbor indices, partitioned by source node. */
	TArray<IndexType> Neighbors;

	/** Per-edge payload data, parallel to Neighbors array. EdgeData[i] corresponds to Neighbors[i]. */
	TArray<EdgeDataType> EdgeData;
};
