use crate::assignment::connected_components::ClusterSize;
use network_simplex::ns::network_simplex::solve_min_cost_flow;
use std::collections::{HashMap, HashSet};
use tracing::{info, warn};

use super::connected_components::TranscriptConnectedComponentLabeling;
use crate::util::oarfish_types::{AlnInfo, EMInfo};
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
    _make_iter: F,
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
                let post_prob = (counts[target_id] + 0.5) * cond_prob;

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

        let mut total_cap = 0.0f64;
        // set transcript capacities
        for t in component_txps.iter() {
            total_cap += remaining_counts[*t as usize];
            // get the count; take the ceiling to ensure that we can assign all
            // reads using integral capacities
            let preferred_cap = remaining_counts[*t as usize].round() as isize;
            let max_cap = 1.max((remaining_counts[*t as usize] * 1.5).ceil() as isize);

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

        if total_cap + 0.1f64 < nread as f64 {
            warn!("total capacity : {total_cap} but nreads is {nread}");
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
        for global_read_id in reads_for_component {
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
