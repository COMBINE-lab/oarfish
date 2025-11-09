//use mcmf::{Capacity, Cost, GraphBuilder, Vertex};
use crate::assignment::connected_components::ClusterSize;
//use ordered_float::OrderedFloat;
use network_simplex::graph::*;
use network_simplex::ns::network_simplex::solve_min_cost_flow;
use std::collections::{HashMap, HashSet, VecDeque, hash_map::OccupiedEntry};
use std::io::Write;
use tracing::{info, warn};
//use mcmf::*;

#[derive(Debug, Clone)]
pub struct Edge {
    pub from: usize,
    pub to: usize,
    pub capacity: i64,
    pub cost: f64,
    pub flow: i64,
}

#[derive(Debug)]
pub struct MinCostMaxFlow {
    n: usize,
    edges: Vec<Edge>,
    graph: Vec<Vec<usize>>, // adjacency list storing edge indices
}

impl MinCostMaxFlow {
    pub fn new(n: usize) -> Self {
        Self {
            n,
            edges: Vec::new(),
            graph: vec![Vec::new(); n],
        }
    }

    pub fn add_edge(
        &mut self,
        from: usize,
        to: usize,
        capacity: i64,
        cost: f64,
    ) -> anyhow::Result<()> {
        if from >= self.n {
            anyhow::bail!("from index {from} greater than num nodes = {}", self.n);
        }
        if to >= self.n {
            anyhow::bail!("to index {to} greater than num nodes = {}", self.n);
        }
        let edge_idx = self.edges.len();

        // Forward edge
        self.edges.push(Edge {
            from,
            to,
            capacity,
            cost,
            flow: 0,
        });
        self.graph[from].push(edge_idx);

        // Backward edge (reverse)
        self.edges.push(Edge {
            from: to,
            to: from,
            capacity: 0,
            cost: -cost,
            flow: 0,
        });
        self.graph[to].push(edge_idx + 1);
        Ok(())
    }

    /// Find minimum cost maximum flow from source to sink
    pub fn min_cost_max_flow(&mut self, source: usize, sink: usize) -> (i64, f64) {
        let mut total_flow = 0i64;
        let mut total_cost = 0f64;
        let mut iteration = 0;

        loop {
            iteration += 1;
            // Use SPFA (Shortest Path Faster Algorithm) to find shortest path
            let (dist, parent) = self.spfa(source, sink);

            if dist[sink] == f64::INFINITY {
                if self.edges.len() == 0 {
                    eprintln!(
                        "Iteration {}: No path found (dist[sink] = infinity)",
                        iteration
                    );
                }
                break; // No more augmenting paths
            }
            if self.edges.len() == 0 {
                eprintln!(
                    "Iteration {}: Found path with cost {}",
                    iteration, dist[sink]
                );
            }

            // Find minimum residual capacity along the path
            // Also check for cycles to prevent infinite loops
            let mut path_flow = i64::MAX;
            let mut node = sink;
            let mut visited = vec![false; self.n];
            let mut valid_path = true;
            let mut path_edges = Vec::new();

            while node != source {
                if visited[node] {
                    // Cycle detected - invalid path
                    valid_path = false;
                    if self.edges.len() == 0 {
                        eprintln!("  Cycle detected at node {}", node);
                    }
                    break;
                }
                visited[node] = true;

                if let Some(edge_idx) = parent[node] {
                    let edge = &self.edges[edge_idx];
                    let residual = edge.capacity - edge.flow;
                    path_edges.push((edge.from, edge.to, residual, edge_idx));
                    path_flow = path_flow.min(residual);
                    node = edge.from;
                } else {
                    // No parent - path doesn't reach source
                    valid_path = false;
                    if self.edges.len() == 0 {
                        eprintln!("  No parent for node {}", node);
                    }
                    break;
                }
            }

            if !valid_path || path_flow == 0 {
                if path_flow == 0 {
                    if self.edges.len() == 0 {
                        eprintln!("  Path found but bottleneck capacity is 0");
                    }
                    for (from, to, residual, idx) in path_edges.iter().rev() {
                        if self.edges.len() == 0 {
                            eprintln!(
                                "    {} -> {} (edge {}, residual: {})",
                                from, to, idx, residual
                            );
                        }
                    }
                }
                if !valid_path {
                    if self.edges.len() == 0 {
                        eprintln!("  Attempting to trace parent path from sink {}:", sink);
                    }
                    let mut trace_node = sink;
                    let mut trace_count = 0;
                    while trace_count < self.n * 2 {
                        if let Some(edge_idx) = parent[trace_node] {
                            let edge = &self.edges[edge_idx];
                            if self.edges.len() == 0 {
                                eprintln!(
                                    "    {} <- {} (edge {})",
                                    trace_node, edge.from, edge_idx
                                );
                            }
                            trace_node = edge.from;
                            trace_count += 1;
                            if trace_node == source {
                                if self.edges.len() == 0 {
                                    eprintln!("    Reached source!");
                                }
                                break;
                            }
                        } else {
                            if self.edges.len() == 0 {
                                eprintln!("    {} has no parent", trace_node);
                            }
                            break;
                        }
                    }
                }
                break; // No valid augmenting path
            }

            if self.edges.len() == 0 {
                eprintln!("  Sending {} units of flow", path_flow);
            }

            // Update flow along the path
            node = sink;
            while node != source {
                let edge_idx = parent[node].unwrap();
                self.edges[edge_idx].flow += path_flow;
                // Update reverse edge
                self.edges[edge_idx ^ 1].flow -= path_flow;
                total_cost += path_flow as f64 * self.edges[edge_idx].cost;
                node = self.edges[edge_idx].from;
            }

            total_flow += path_flow;
        }

        (total_flow, total_cost)
    }

    /// SPFA algorithm for finding shortest path with negative edge costs
    /// Includes cycle detection to avoid infinite loops
    fn spfa(&self, source: usize, sink: usize) -> (Vec<f64>, Vec<Option<usize>>) {
        let mut dist = vec![f64::INFINITY; self.n];
        let mut parent = vec![None; self.n];
        let mut in_queue = vec![false; self.n];
        let mut cnt = vec![0; self.n]; // Count of times each node is added to queue
        let mut queue = VecDeque::new();

        dist[source] = 0.0;
        queue.push_back(source);
        in_queue[source] = true;
        cnt[source] = 1;

        while let Some(u) = queue.pop_front() {
            in_queue[u] = false;

            for &edge_idx in &self.graph[u] {
                let edge = &self.edges[edge_idx];

                // Check if edge has residual capacity
                if edge.flow < edge.capacity {
                    let v = edge.to;
                    let new_dist = dist[u] + edge.cost;

                    // Use a small epsilon for floating point comparison
                    if new_dist < dist[v] {
                        dist[v] = new_dist;

                        // Check if setting parent[v] = u would create a cycle
                        // by checking if u is reachable from v through current parents
                        let mut creates_cycle = false;
                        let mut check_node = u;
                        let mut visited_check = vec![false; self.n];
                        while check_node != source {
                            if visited_check[check_node] || check_node == v {
                                creates_cycle = true;
                                break;
                            }
                            visited_check[check_node] = true;

                            if let Some(p_edge_idx) = parent[check_node] {
                                check_node = self.edges[p_edge_idx as usize].from;
                            } else {
                                break;
                            }
                        }

                        if creates_cycle {
                            // Skip this edge as it would create a cycle in the parent tree
                            continue;
                        }

                        parent[v] = Some(edge_idx);

                        if !in_queue[v] {
                            queue.push_back(v);
                            in_queue[v] = true;
                            cnt[v] += 1;

                            // If a node is added to queue more than n times, there's a negative cycle
                            if cnt[v] > self.n {
                                return (vec![f64::INFINITY; self.n], vec![None; self.n]);
                            }
                        }
                    }
                }
            }
        }

        (dist, parent)
    }

    /// Get all edges with non-zero flow
    pub fn get_flow_edges(&self) -> Vec<Edge> {
        self.edges
            .iter()
            .step_by(2) // Only forward edges (skip reverse edges)
            .filter(|e| e.flow > 0)
            .cloned()
            .collect()
    }

    /// Get all edges (including those with zero flow)
    pub fn get_all_edges(&self) -> Vec<Edge> {
        self.edges
            .iter()
            .step_by(2) // Only forward edges
            .cloned()
            .collect()
    }

    /// Get the flow paths as a list of sequences
    pub fn get_flow_paths(&self, source: usize, sink: usize) -> Vec<Vec<(usize, usize, i64, f64)>> {
        let mut paths = Vec::new();
        let mut remaining_flow: HashMap<(usize, usize), i64> = HashMap::new();

        // Build a map of remaining flows
        for edge in self.edges.iter().step_by(2) {
            if edge.flow > 0 {
                remaining_flow.insert((edge.from, edge.to), edge.flow);
            }
        }

        // Extract paths using DFS
        while let Some(&start_edge) = remaining_flow.keys().find(|&&(from, _)| from == source) {
            let mut path = Vec::new();
            let mut current = start_edge.0;
            let mut min_flow = i64::MAX;

            // Build path from source to sink
            let mut visited = vec![false; self.n];
            while current != sink {
                if visited[current] {
                    break; // Cycle detected, break out
                }
                visited[current] = true;

                if let Some(((from, to), flow)) = remaining_flow
                    .iter()
                    .find(|((f, _), _)| *f == current)
                    .map(|(k, v)| (*k, *v))
                {
                    // Find the edge to get its cost
                    let edge = self
                        .edges
                        .iter()
                        .step_by(2)
                        .find(|e| e.from == from && e.to == to)
                        .unwrap();

                    path.push((from, to, flow, edge.cost));
                    min_flow = min_flow.min(flow);
                    current = to;
                } else {
                    break;
                }
            }

            if current == sink && !path.is_empty() {
                // Subtract the minimum flow from all edges in the path
                for &(from, to, _, _) in &path {
                    let flow = remaining_flow.get_mut(&(from, to)).unwrap();
                    *flow -= min_flow;
                    if *flow == 0 {
                        remaining_flow.remove(&(from, to));
                    }
                }

                // Store the path with the actual flow value
                let path_with_flow: Vec<_> = path
                    .iter()
                    .map(|&(from, to, _, cost)| (from, to, min_flow, cost))
                    .collect();
                paths.push(path_with_flow);
            } else {
                break; // Can't form a complete path
            }
        }

        paths
    }
}

/*
#[derive(Clone, Debug)]
struct Edge {
    to: usize,
    cap: i64,
    cost: f64,
    rev: usize, // Index of reverse edge in adjacency list
    flow: i64,  // Current flow through this edge
}

#[derive(Clone, Debug)]
pub struct FlowEdge {
    pub from: usize,
    pub to: usize,
    pub flow: i64,
    pub capacity: i64,
    pub cost: f64,
}

pub struct MinCostFlow {
    graph: Vec<Vec<Edge>>,
    n: usize,
    original_n: usize,    // Number of original nodes (before node splitting)
    node_in: Vec<usize>,  // Maps original node to its "in" node
    node_out: Vec<usize>, // Maps original node to its "out" node
}

impl MinCostFlow {
    /// Creates a new min-cost flow graph with n nodes
    pub fn new(n: usize) -> Self {
        MinCostFlow {
            graph: vec![Vec::new(); n],
            n,
            original_n: n,
            node_in: (0..n).collect(),
            node_out: (0..n).collect(),
        }
    }

    /// Sets a capacity constraint on a node by splitting it into in/out nodes
    /// All existing edges TO this node will connect to the "in" node
    /// All existing edges FROM this node will connect from the "out" node
    /// A new edge with the given capacity connects in -> out
    /// Cost of 0.0 is used for the internal edge
    pub fn set_node_capacity(&mut self, node: usize, capacity: i64) {
        if node >= self.original_n {
            panic!("Node index out of bounds");
        }

        // If this node was already split, don't split again
        if self.node_in[node] != self.node_out[node] {
            // Update the capacity of the existing internal edge
            let in_node = self.node_in[node];
            let out_node = self.node_out[node];

            // Find and update the edge from in_node to out_node
            for edge in &mut self.graph[in_node] {
                if edge.to == out_node && edge.cost == 0.0 {
                    edge.cap = capacity;
                    // Also update the original capacity for extract_flow
                    let original_flow = edge.flow;
                    edge.cap = capacity - original_flow;
                    break;
                }
            }
            return;
        }

        let in_node = self.n;
        let out_node = self.n + 1;

        self.node_in[node] = in_node;
        self.node_out[node] = out_node;

        // Add two new nodes
        self.graph.push(Vec::new());
        self.graph.push(Vec::new());
        self.n += 2;

        // Redirect all incoming edges to point to in_node
        for u in 0..in_node {
            for edge in &mut self.graph[u] {
                if edge.to == node {
                    edge.to = in_node;
                }
            }
        }

        // Move all outgoing edges from node to out_node
        let outgoing = std::mem::take(&mut self.graph[node]);

        // Update reverse edge pointers for moved edges before assigning
        let to_update: Vec<(usize, usize)> =
            outgoing.iter().map(|edge| (edge.to, edge.rev)).collect();

        self.graph[out_node] = outgoing;

        // Now update the reverse edges
        for (to, rev_idx) in to_update {
            self.graph[to][rev_idx].to = out_node;
        }

        // Add internal edge with capacity constraint (0 cost)
        self.add_edge_internal(in_node, out_node, capacity, 0.0);
    }

    /// Adds a directed edge from u to v with capacity and cost
    /// Automatically handles node splitting for capacity-constrained nodes
    pub fn add_edge(&mut self, from: usize, to: usize, cap: i64, cost: f64) -> anyhow::Result<()> {
        if from >= self.original_n || to >= self.original_n {
            anyhow::bail!("Node index out of bounds");
        }

        let from_actual = self.node_out[from];
        let to_actual = self.node_in[to];

        self.add_edge_internal(from_actual, to_actual, cap, cost);
        Ok(())
    }

    /// Internal method to add edges directly between actual nodes
    fn add_edge_internal(&mut self, from: usize, to: usize, cap: i64, cost: f64) {
        let from_len = self.graph[from].len();
        let to_len = self.graph[to].len();

        self.graph[from].push(Edge {
            to,
            cap,
            cost,
            rev: to_len,
            flow: 0,
        });

        // Add reverse edge with 0 capacity and negative cost
        self.graph[to].push(Edge {
            to: from,
            cap: 0,
            cost: -cost,
            rev: from_len,
            flow: 0,
        });
    }

    /// Computes minimum cost flow from source to sink with given flow amount
    /// Returns (flow_achieved, total_cost, flow_edges) or None if flow cannot be satisfied
    /// flow_edges contains (from, to, flow_amount, cost) for each edge with positive flow
    pub fn min_cost_flow(
        &mut self,
        source: usize,
        sink: usize,
        mut flow: i64,
    ) -> Option<(i64, f64, Vec<FlowEdge>)> {
        if source >= self.original_n || sink >= self.original_n {
            panic!("Node index out of bounds");
        }

        let source_actual = self.node_out[source];
        let sink_actual = self.node_in[sink];

        let mut total_cost = 0.0;
        let mut total_flow = 0i64;

        while flow > 0 {
            // Find shortest path using Bellman-Ford (handles negative costs)
            let (dist, parent) = self.bellman_ford(source_actual, sink_actual)?;

            if dist[sink_actual] == f64::INFINITY {
                // No more augmenting paths
                break;
            }

            // Find minimum capacity along the path
            let mut path_flow = flow;
            let mut v = sink_actual;

            while v != source_actual {
                let (u, edge_idx) = parent[v];
                path_flow = path_flow.min(self.graph[u][edge_idx].cap);
                v = u;
            }

            if path_flow == 0 {
                break;
            }

            // Update flow along the path
            flow -= path_flow;
            total_flow += path_flow;
            total_cost += path_flow as f64 * dist[sink_actual];

            v = sink_actual;
            while v != source_actual {
                let (u, edge_idx) = parent[v];

                // Forward edge
                self.graph[u][edge_idx].cap -= path_flow;
                self.graph[u][edge_idx].flow += path_flow;

                // Reverse edge
                let rev_idx = self.graph[u][edge_idx].rev;
                self.graph[v][rev_idx].cap += path_flow;
                self.graph[v][rev_idx].flow -= path_flow;

                v = u;
            }
        }

        if total_flow == 0 {
            None
        } else {
            let flow_edges = self.extract_flow();
            Some((total_flow, total_cost, flow_edges))
        }
    }

    /// Extracts all edges with positive flow
    /// Maps internal nodes back to original node IDs
    pub fn extract_flow(&self) -> Vec<FlowEdge> {
        let mut result = Vec::new();

        // Create reverse mapping from actual nodes to original nodes
        let mut actual_to_original = vec![None; self.n];
        for original in 0..self.original_n {
            actual_to_original[self.node_in[original]] = Some(original);
            actual_to_original[self.node_out[original]] = Some(original);
        }

        for from_actual in 0..self.n {
            for edge in &self.graph[from_actual] {
                // Only include forward edges (non-negative cost) with positive flow
                if edge.cost >= 0.0 && edge.flow > 0 {
                    // Map back to original node IDs
                    let from_orig = actual_to_original[from_actual];
                    let to_orig = actual_to_original[edge.to];

                    // Skip internal capacity-constraint edges (where from == to in original graph)
                    if from_orig == to_orig {
                        continue;
                    }

                    if let (Some(from), Some(to)) = (from_orig, to_orig) {
                        result.push(FlowEdge {
                            from,
                            to,
                            flow: edge.flow,
                            capacity: edge.cap + edge.flow, // Original capacity
                            cost: edge.cost,
                        });
                    }
                }
            }
        }

        result
    }

    /// Bellman-Ford algorithm for finding shortest path with negative edge weights
    /// Returns (distances, parent) where parent[v] = (u, edge_index)
    fn bellman_ford(&self, source: usize, sink: usize) -> Option<(Vec<f64>, Vec<(usize, usize)>)> {
        let mut dist = vec![f64::INFINITY; self.n];
        let mut parent = vec![(0, 0); self.n];
        let mut in_queue = vec![false; self.n];
        let mut queue = VecDeque::new();

        dist[source] = 0.0;
        queue.push_back(source);
        in_queue[source] = true;

        // SPFA (Shortest Path Faster Algorithm) variant of Bellman-Ford
        let mut iterations = 0;
        while let Some(u) = queue.pop_front() {
            in_queue[u] = false;
            iterations += 1;

            // Detect negative cycles
            if iterations > self.n * self.n {
                return None;
            }

            for (idx, edge) in self.graph[u].iter().enumerate() {
                if edge.cap > 0 {
                    let new_dist = dist[u] + edge.cost;

                    if new_dist < dist[edge.to] {
                        dist[edge.to] = new_dist;
                        parent[edge.to] = (u, idx);

                        if !in_queue[edge.to] {
                            queue.push_back(edge.to);
                            in_queue[edge.to] = true;
                        }
                    }
                }
            }
        }

        Some((dist, parent))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple_flow() {
        let mut mcf = MinCostFlow::new(4);

        // Simple graph: 0 -> 1 -> 3
        //               0 -> 2 -> 3
        mcf.add_edge(0, 1, 10, 2.0);
        mcf.add_edge(0, 2, 10, 4.0);
        mcf.add_edge(1, 3, 10, 1.0);
        mcf.add_edge(2, 3, 10, 1.0);

        let result = mcf.min_cost_flow(0, 3, 10);
        assert!(result.is_some());

        let (flow, cost, edges) = result.unwrap();
        assert_eq!(flow, 10);
        assert!((cost - 30.0).abs() < 1e-9); // Path 0->1->3 costs 3 per unit

        // Verify flow conservation
        let total_flow: i64 = edges.iter().filter(|e| e.from == 0).map(|e| e.flow).sum();
        assert_eq!(total_flow, 10);
    }

    #[test]
    fn test_multiple_paths() {
        let mut mcf = MinCostFlow::new(4);

        mcf.add_edge(0, 1, 5, 1.0);
        mcf.add_edge(0, 2, 5, 2.0);
        mcf.add_edge(1, 3, 5, 3.0);
        mcf.add_edge(2, 3, 5, 1.0);

        let result = mcf.min_cost_flow(0, 3, 10);
        assert!(result.is_some());

        let (flow, _cost, edges) = result.unwrap();
        assert_eq!(flow, 10);

        // Verify each edge has non-negative flow
        for edge in &edges {
            assert!(edge.flow >= 0);
            assert!(edge.flow <= edge.capacity);
        }
    }

    #[test]
    fn test_node_capacity() {
        let mut mcf = MinCostFlow::new(4);

        // Graph: 0 -> 1 -> 3
        //        0 -> 2 -> 3
        // But node 1 has capacity 5
        mcf.add_edge(0, 1, 10, 1.0);
        mcf.add_edge(0, 2, 10, 2.0);
        mcf.add_edge(1, 3, 10, 1.0);
        mcf.add_edge(2, 3, 10, 1.0);

        mcf.set_node_capacity(1, 5); // Limit flow through node 1

        let result = mcf.min_cost_flow(0, 3, 10);
        assert!(result.is_some());

        let (flow, _cost, edges) = result.unwrap();
        assert_eq!(flow, 10);

        // Check that flow through node 1 doesn't exceed 5
        let flow_through_1: i64 = edges.iter().filter(|e| e.to == 1).map(|e| e.flow).sum();
        assert!(flow_through_1 <= 5);
    }
}
*/
/*
fn main() {
    // Example usage with node capacity constraints
    let mut mcf = MinCostFlow::new(5);

    // Build a sample graph
    mcf.add_edge(0, 1, 10, 1.5);
    mcf.add_edge(0, 2, 10, 2.0);
    mcf.add_edge(1, 2, 5, 0.5);
    mcf.add_edge(1, 3, 10, 3.0);
    mcf.add_edge(2, 3, 10, 1.0);
    mcf.add_edge(2, 4, 5, 2.5);
    mcf.add_edge(3, 4, 10, 1.0);

    // Add capacity constraint on node 2 (max 8 units can flow through it)
    mcf.set_node_capacity(2, 8);

    match mcf.min_cost_flow(0, 4, 15) {
        Some((flow, cost, edges)) => {
            println!("Maximum flow from 0 to 4: {}", flow);
            println!("Minimum cost: {:.2}", cost);
            println!("\nFlow decomposition:");
            for edge in &edges {
                println!(
                    "  {} -> {}: {}/{} units at cost {:.2} (total: {:.2})",
                    edge.from, edge.to, edge.flow, edge.capacity,
                    edge.cost, edge.flow as f64 * edge.cost
                );
            }

            // Verify node capacity constraint
            let flow_through_2: i64 = edges.iter()
                .filter(|e| e.to == 2)
                .map(|e| e.flow)
                .sum();
            println!("\nFlow through node 2: {} (capacity: 8)", flow_through_2);
        }
        None => {
            println!("Cannot achieve the desired flow");
        }
    }
}
*/

use super::connected_components::TranscriptConnectedComponentLabeling;
use crate::{
    assignment,
    util::oarfish_types::{AlnInfo, EMInfo},
};
use anyhow::Context;
use itertools::*;

#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum FlowVertex {
    Source,
    Sink,
    Read(usize),
    Transcript(usize),
    Penalty(usize),
}

pub fn solve<'a, I: Iterator<Item = (&'a [AlnInfo], &'a [f32], &'a [f64])> + 'a, F: Fn() -> I>(
    labels: &mut TranscriptConnectedComponentLabeling,
    em_info: &EMInfo,
    counts: &[f64],
    make_iter: F,
) -> anyhow::Result<Vec<u32>> {
    let make_iter = || em_info.eq_map.iter();
    let fops = &em_info.eq_map.filter_opts;

    let model_coverage = fops.model_coverage;
    let _density_fn = |x, y| -> f64 {
        match em_info.kde_model {
            Some(ref kde_model) => kde_model[(x, y)],
            _ => 1.,
        }
    };

    let total_reads = labels
        .cc_sizes
        .iter()
        .fold(0_usize, |acc, c| acc + c.1.nread as usize);

    const READ_NOT_MAPPED: u32 = u32::MAX;
    // will hold final read assignments
    let mut read_assignments = vec![READ_NOT_MAPPED; total_reads];
    let mut remaining_counts = counts.to_vec();

    // go over all alignments and, for each unique alignment,
    // assign this read to the corresponding transcript
    /*
    for (global_read_id, (alns, _probs, _coverage_probs)) in make_iter().enumerate() {
        if alns.len() == 1 {
            let target_id = alns[0].ref_id;
            let cluster = labels.labels[target_id as usize];
            // we took care of 1 read from this cluster
            labels.cc_sizes.entry(cluster).and_modify(|e| e.nread -= 1);

            // the read is assigned to this transcript
            read_assignments[global_read_id] = target_id;
            // the future capacity of this transcript decreases by 1
            remaining_counts[target_id as usize] -= 1.0f64;
        }
    }
    */

    // first, relabel the clusters so that we always see things in a consistent order
    // map from the component label (arbitrary transcript ID chosen by union find) to the
    // component's iteration order (rank).
    let nonempty_component_ranks = labels
        .cc_sizes
        .iter()
        .filter(|(_k, v)| v.nread > 0)
        .enumerate()
        .map(|(i, (k, _v))| (*k, i as u32))
        .collect::<std::collections::HashMap<u32, u32>>();
    // change the labels to match, giving an INVALID label id to empty components
    labels
        .labels
        .iter_mut()
        .for_each(|x| *x = *nonempty_component_ranks.get(x).unwrap_or(&u32::MAX));
    let cc_sizes = labels
        .cc_sizes
        .iter()
        .filter(|(_k, v)| v.nread > 0)
        .map(|(_k, v)| v.clone())
        .collect::<Vec<ClusterSize>>();

    // iterating the components in the prescribed order, how many total reads have we
    // seen at each component boundary
    let mut component_prefix_sums = cc_sizes
        .iter()
        .scan(0, |sum, i| {
            *sum += i.nread as usize;
            Some(*sum)
        })
        .collect::<Vec<_>>();

    // explicitly keep track of the total prefix sum, as we'll need it later once
    // we've filled in the read lists for each component
    let last_count = *component_prefix_sums.last().expect("should not be empty");
    let total_unassigned_reads = last_count;

    // will hold (contiguously) the ids for the reads for each component
    let mut component_read_lists = vec![u32::MAX; total_unassigned_reads];
    // set up read ids by component
    for (global_read_id, (alns, _probs, _coverage_probs)) in make_iter().enumerate() {
        // skip over the uniquely aligned reads because they are already dealt with
        //if alns.len() == 1 {
        //    continue;
        //}

        // get the component id
        let component_label = labels.labels[alns[0].ref_id as usize] as usize;
        // assert that there are still reads we can assign to this component!
        assert_ne!(cc_sizes[component_label].nread, 0_u32);

        // we start at the total number of reads in the component (one past the last index),
        // so we'd always like to place at the index one less. Hence, we pre-decremnt here
        component_prefix_sums[component_label] -= 1;
        let next_index = component_prefix_sums[component_label];
        component_read_lists[next_index] = global_read_id as u32;
    }

    // now, we've sutracted the appropriate numbers from each prefix sum, so that the
    // i-th entry now is where the (i-1)-st entry was before. This means we just need to
    // re-add the last item.
    component_prefix_sums.push(last_count);

    // because this could take some time.
    let bar = indicatif::ProgressBar::new(cc_sizes.len() as u64);
    bar.set_style(
        indicatif::ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:20.green/blue} {msg} {human_pos:>12}",
        )
        .unwrap()
        .tick_chars("⠁⠁⠉⠙⠚⠒⠂⠂⠒⠲⠴⠤⠄⠄⠤⠠⠠⠤⠦⠖⠒⠐⠐⠒⠓⠋⠉⠈⠈"),
    );
    bar.set_draw_target(indicatif::ProgressDrawTarget::stderr_with_hz(4));

    const COST_SCALE: f64 = 5f64;
    let mut component_txps = HashSet::<u32>::with_capacity(100);

    // now we have the component_read_lists vector that gives us a list of the specific reads that
    // are relevant to a component, allowing us to iterate over the reads relevant for a component
    for (component_root, size) in cc_sizes.iter().enumerate() {
        // if this component has no reads remaining
        // (namely, if it was a single transcript component so all
        // reads were unique) then we've nothing to do here.
        assert!(size.nread > 0);

        bar.set_message(format!("procesing {component_root}"));

        // number of reads remaining in this component graph
        let nread = size.nread as usize;

        // number of transcripts in this component graph
        let ntxp = size.ntxp as usize;

        // trying https://github.com/andryr/network-simplex-rust/
        let mut mcf_builder = network_simplex::graph::GraphBuilder::<FlowVertex>::new();
        mcf_builder.add_node(FlowVertex::Source, nread as f64);
        mcf_builder.add_node(FlowVertex::Sink, -(nread as f64));

        // build the flow graph -- iterate over the reads for this component
        // let label_rank = component_ranks[component_root] as usize;
        let start = component_prefix_sums[component_root];
        let end = component_prefix_sums[component_root + 1];
        assert_eq!(end - start, size.nread as usize);

        component_txps.clear();
        let reads_for_component = &component_read_lists[start..end];
        let mut distinct_edges = HashMap::<u32, (f64, bool)>::with_capacity(100);
        for (local_rank, global_read_id) in reads_for_component.iter().enumerate() {
            bar.set_message(format!(
                "procesing read {local_rank} of {nread} (# transcripts {ntxp})"
            ));
            let global_read_id = *global_read_id as usize;
            let (alns, probs, coverage_probs) = em_info
                .eq_map
                .get_alignments_for_read(global_read_id)
                .with_context(|| format!("expecting a valid read id, got {global_read_id}"))?;
            // we already dealt with all uniquely-mapped reads
            //assert_ne!(alns.len(), 1);

            distinct_edges.clear();
            let mut denom = 0.0f64;
            for (a, p, cp) in izip!(alns, probs, coverage_probs) {
                // Compute the probability of assignment of the
                // current read based on this alignment and the
                // target's estimated abundance.
                let target_id = a.ref_id as usize;

                //let txp_len = tinfo[target_id].lenf as usize;
                //let aln_len = a.alignment_span() as usize;
                let prob = *p as f64;
                let cov_prob = if model_coverage { *cp } else { 1.0 };
                let cond_prob = prob * cov_prob;
                let post_prob = counts[target_id] * cond_prob;

                distinct_edges
                    .entry(target_id as u32)
                    .and_modify(|p| p.0 += post_prob)
                    .or_insert((post_prob, false));

                denom += post_prob;
            }
            // edge from source to read node
            // network simplex
            mcf_builder.add_edge(
                FlowVertex::Source,
                FlowVertex::Read(global_read_id),
                1f64,
                0f64,
            );

            for (a, _p, _cp) in izip!(alns, probs, coverage_probs) {
                let target_id = a.ref_id as usize;

                // edge for this alignment
                // let prob = *p as f64;
                //let cov_prob = if model_coverage { *cp } else { 1.0 };
                //let post_prob = (counts[target_id] * prob * cov_prob) / denom;

                if let std::collections::hash_map::Entry::Occupied(mut o) =
                    distinct_edges.entry(target_id as u32)
                {
                    let ent = o.get_mut();
                    if !ent.1 {
                        let post_prob = ent.0 / denom;
                        let cost = COST_SCALE - (COST_SCALE * post_prob);
                        mcf_builder.add_edge(
                            FlowVertex::Read(global_read_id),
                            FlowVertex::Transcript(target_id),
                            1.0f64,
                            cost,
                        );
                        ent.1 = true;
                        component_txps.insert(target_id as u32);
                    }
                }
            }
        }

        // set transcript capacities
        for t in component_txps.iter() {
            // get the count; take the ceiling to ensure that we can assign all
            // reads using integral capacities
            let preferred_cap = remaining_counts[*t as usize].round() as isize;
            let max_cap = 1isize.max((remaining_counts[*t as usize] * 1.5).ceil() as isize);

            if preferred_cap > 0 {
                mcf_builder.add_edge(
                    FlowVertex::Transcript(*t as usize),
                    FlowVertex::Sink,
                    preferred_cap as f64,
                    0f64,
                );
            }

            // Second edge: overflow capacity with penalty
            let overflow_cap = max_cap.saturating_sub(preferred_cap);
            if overflow_cap > 0 {
                // Penalty for using overflow capacity
                // Higher penalty = less likely to use
                let penalty_cost = 10.0f64; // Tune this
                mcf_builder.add_edge(
                    FlowVertex::Transcript(*t as usize),
                    FlowVertex::Penalty(*t as usize),
                    overflow_cap as f64,
                    0.0f64,
                );
                mcf_builder.add_edge(
                    FlowVertex::Penalty(*t as usize),
                    FlowVertex::Sink,
                    overflow_cap as f64,
                    penalty_cost,
                );
            }
        }
        bar.set_message(format!(
            "running min cost flow for component {component_root} of size ({}, {})",
            size.ntxp, size.nread
        ));
        //network simplex version
        let mut total_cost = 0.0f64;
        let index_to_label = mcf_builder
            .node_label_to_index
            .iter()
            .map(|(k, v)| (*v, k.clone()))
            .collect::<std::collections::HashMap<usize, FlowVertex>>();
        let graph = mcf_builder.build();
        let flow_values = solve_min_cost_flow(&graph, 10e-9); //mcf_builder.get_flow_edges();
        for (e, f) in graph.edges.iter().zip(flow_values) {
            if let (&FlowVertex::Read(global_read_id), &FlowVertex::Transcript(global_txp_id)) = (
                index_to_label.get(&e.start).expect("valid start"),
                index_to_label.get(&e.end).expect("valid end"),
            ) {
                if f >= 0.999f64 {
                    total_cost += e.cost;
                    read_assignments[global_read_id] = global_txp_id as u32;
                    remaining_counts[global_txp_id] -= 1.0f64;
                }
            }
        }
        //info!("moving to next component");
        bar.set_message(format!(
            "finished flow problem for component {component_root}. Minimum cost was {total_cost}"
        ));

        // ensure that we assigned all of the reads in this component!
        for (local_rank, global_read_id) in reads_for_component.iter().enumerate() {
            if read_assignments[*global_read_id as usize] == u32::MAX {
                warn!(
                    "Read {global_read_id} was part of this component {component_root} with size {size:#?}, but had no assignment!"
                );
            }
        }
        bar.inc(1);
    }
    bar.finish();

    let unassigned_reads = read_assignments
        .iter()
        .enumerate()
        .filter_map(|(idx, txp_id)| if *txp_id == u32::MAX { Some(idx) } else { None })
        .collect::<Vec<usize>>();
    info!("unassigned read indices = {:#?}", unassigned_reads);
    Ok(read_assignments)
}
