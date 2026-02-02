# Compact Sparse Row Graph for Unreal Engine

A high-performance, memory-efficient Compressed Sparse Row (CSR) graph data structure plugin for Unreal Engine 5.

## Overview

CSR is a compact representation for sparse graphs that stores all neighbor indices in a single contiguous array, with a separate offset array indicating where each node's neighbors begin. This layout provides:

- **Optimal cache locality** for neighbor iteration
- **Minimal memory overhead** compared to adjacency lists or matrices
- **O(1) access** to any node's neighbor list
- **Cache-friendly traversal** for graph algorithms

## Features

- `TCSRGraph<IndexType>` - Basic CSR graph for adjacency relationships
- `TCSRGraphWithEdgeData<IndexType, EdgeDataType>` - CSR graph with per-edge payload data (weights, labels, etc.)
- Optional neighbor sorting for deterministic iteration
- Optional duplicate edge removal with "first edge wins" semantics
- Header-only implementation - no linking required
- Zero external dependencies beyond Unreal Engine core

## Installation

### As a Project Plugin

1. Copy the `CompactSparseRow` folder to your project's `Plugins/` directory
2. Regenerate project files
3. Add `"CompactSparseRow"` to your module's `Build.cs`:

```csharp
PublicDependencyModuleNames.AddRange(new string[] { "CompactSparseRow" });
```

### As an Engine Plugin

1. Copy the `CompactSparseRow` folder to `Engine/Plugins/Runtime/`
2. Regenerate project files

## Usage

### Basic Graph

```cpp
#include "CompactSparseRow.h"

// Create a graph with 4 nodes
TCSRGraph<int32> Graph;
TArray<TCSRGraph<int32>::FEdge> Edges;

// Add directed edges
Edges.Add({0, 1});  // 0 -> 1
Edges.Add({0, 2});  // 0 -> 2
Edges.Add({1, 2});  // 1 -> 2
Edges.Add({2, 3});  // 2 -> 3
Edges.Add({3, 0});  // 3 -> 0 (cycle)

// Build the CSR structure
// Parameters: NumNodes, Edges, bSortNeighbors, bRemoveDuplicates
Graph.Build(4, Edges, true, true);

// Iterate neighbors of node 0
for (int32 Neighbor : Graph.GetNeighbors(0))
{
    UE_LOG(LogTemp, Log, TEXT("Node 0 -> Node %d"), Neighbor);
}
```

### Weighted Graph

```cpp
#include "CompactSparseRow.h"

// Graph with float edge weights
TCSRGraphWithEdgeData<int32, float> WeightedGraph;
TArray<TCSRGraphWithEdgeData<int32, float>::FEdge> Edges;

// Add weighted edges
Edges.Add({0, 1, 1.5f});   // 0 -> 1, weight 1.5
Edges.Add({0, 2, 2.0f});   // 0 -> 2, weight 2.0
Edges.Add({1, 2, 0.5f});   // 1 -> 2, weight 0.5

WeightedGraph.Build(3, Edges);

// Access neighbors and weights together
TConstArrayView<int32> Neighbors = WeightedGraph.GetNeighbors(0);
TConstArrayView<float> Weights = WeightedGraph.GetEdgeData(0);

for (int32 i = 0; i < Neighbors.Num(); ++i)
{
    UE_LOG(LogTemp, Log, TEXT("Edge 0 -> %d, Weight: %f"),
           Neighbors[i], Weights[i]);
}
```

### Custom Edge Data

```cpp
// Edge with multiple properties
struct FMyEdgeData
{
    float Weight;
    int32 Category;
    bool bEnabled;
};

TCSRGraphWithEdgeData<int32, FMyEdgeData> Graph;
TArray<TCSRGraphWithEdgeData<int32, FMyEdgeData>::FEdge> Edges;

Edges.Add({0, 1, {1.0f, 2, true}});
Edges.Add({0, 2, {0.5f, 1, false}});

Graph.Build(3, Edges);
```

## API Reference

### TCSRGraph

| Method | Description |
|--------|-------------|
| `Build(NumNodes, Edges, bSort, bDedup)` | Constructs the CSR from an edge list |
| `GetNeighbors(Node)` | Returns `TConstArrayView` of neighbor indices |
| `GetNumNodes()` | Returns the number of nodes |
| `GetNumEdges()` | Returns the total number of edges |
| `IsValidNode(Node)` | Checks if a node index is valid |
| `GetOffsets()` | Returns raw offset array for direct access |
| `GetNeighborArray()` | Returns raw neighbor array for direct access |
| `Reset()` | Clears the graph and releases memory |
| `Empty(Slack)` | Clears the graph but retains allocated memory |

### TCSRGraphWithEdgeData

All methods from `TCSRGraph`, plus:

| Method | Description |
|--------|-------------|
| `GetEdgeData(Node)` | Returns `TConstArrayView` of edge payloads (parallel to neighbors) |
| `GetEdgeDataArray()` | Returns raw edge data array for direct access |

## Memory Layout

```
Example: 4-node graph with edges 0->1, 0->3, 1->0, 1->2, 1->4, 3->0, 3->1, 3->2

Offsets:   [0, 2, 5, 5, 8]     (length = NumNodes + 1)
Neighbors: [1, 3, 0, 2, 4, 0, 1, 2]  (length = NumEdges)

Node 0: Neighbors[0..2) = {1, 3}
Node 1: Neighbors[2..5) = {0, 2, 4}
Node 2: Neighbors[5..5) = {} (no outgoing edges)
Node 3: Neighbors[5..8) = {0, 1, 2}
```

## Performance

| Operation | Complexity |
|-----------|------------|
| Build (unsorted) | O(E) |
| Build (sorted) | O(E log(MaxDegree)) |
| Get neighbors | O(1) |
| Iterate neighbors | O(degree) |
| Memory usage | O(V + E) |

Where V = number of nodes, E = number of edges.

## Use Cases

- **AI Navigation**: Store connectivity between navmesh polygons
- **Pathfinding**: Dijkstra, A*, BFS, DFS implementations
- **Dependency Graphs**: Module loading, task scheduling
- **Social Networks**: Follower/following relationships
- **Spatial Queries**: K-nearest neighbors, range queries
- **Flow Networks**: Max flow, min cut algorithms
- **Scene Graphs**: Parent-child relationships with edge metadata

## Limitations

- **Static structure**: Not designed for frequent edge insertion/removal after `Build()`
- **Directed edges**: For undirected graphs, add edges in both directions
- **int32 edge limit**: Total edges must fit in int32 (~2 billion)
- **Signed IndexType**: IndexType must be a signed integer type

## Requirements

- Unreal Engine 5.0+
- C++17 or later

## License

Copyright Hitbox Games, LLC. All Rights Reserved.

Licensed under the MIT License. See LICENSE file for details.
