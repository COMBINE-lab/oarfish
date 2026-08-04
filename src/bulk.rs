use crate::NamedDigestVec;
use crate::alignment_parser;
use crate::em;
use crate::kde_utils;
use crate::prog_opts::Args;
use crate::util::constants::EMPTY_READ_NAME;
use crate::util::oarfish_types::AlnInfo;
use crate::util::oarfish_types::DiscardTable;
use crate::util::oarfish_types::{
    AlignmentFilters, EMInfo, InMemoryAlignmentStore, InputSourceType, ReadChunkWithNames,
    ReadSource, TranscriptInfo,
};
use crate::util::projection::{mapping_to_genomic_alignment, projected_to_records};
use crate::util::read_function::read_short_quant_vec;
use crate::util::write_function::{
    write_coverage_signals, write_infrep_file, write_out_prob, write_output,
};
use crate::{logistic_prob, normalize_read_probs};
use arrow2::{array::Float64Array, chunk::Chunk, datatypes::Field};
use bramble_rs::ProjectionConfig;
use bramble_rs::ProjectionContext;
use bramble_rs::g2t::G2TTree;
use bramble_rs::project_group_with;
use crossbeam::channel::Receiver;
use crossbeam::channel::Sender;
use crossbeam::channel::bounded;

use needletail::parse_fastx_file;
use noodles_bam as bam;
use num_format::{Locale, ToFormattedString};
use serde_json::json;
use std::io::BufRead;
use swapvec::{SwapVec, SwapVecConfig};
use tracing::{info, warn};

struct RunTiming {
    alignment: std::time::Duration,
    coverage: std::time::Duration,
    em: std::time::Duration,
}

/// Produce a [serde_json::Value] that encodes the relevant arguments and
/// parameters of the run that we wish to record to file. Ultimately, this
/// will be written to the corresponding `meta_info.json` file for this run.
fn get_json_info(
    args: &Args,
    emi: &EMInfo,
    em_result: &em::EMResult,
    seqcol_digest: &NamedDigestVec,
    timing: &RunTiming,
    coverage_diagnostics: &serde_json::Value,
) -> serde_json::Value {
    let prob = match args.coverage_model {
        crate::prog_opts::CoverageModel::None => "no_coverage",
        crate::prog_opts::CoverageModel::Logistic => "logistic_coverage",
        crate::prog_opts::CoverageModel::Endpoint => "endpoint_coverage",
        crate::prog_opts::CoverageModel::Hybrid => "hybrid_coverage",
        crate::prog_opts::CoverageModel::Adaptive => "adaptive_coverage",
        crate::prog_opts::CoverageModel::Degradation => "degradation_aware_coverage",
        crate::prog_opts::CoverageModel::Auto => "automatic_technology_coverage",
    };

    let source = if args.alignments.is_some() {
        "from_bam"
    } else {
        "from_raw_reads"
    };

    json!({
        "prob_model" : prob,
        "alignment_source" : source,
        "alignment_time" : {
            "comment" : "Time to parse (in alignment mode) or generate (in raw read mode) alignments, as well as apply filters, and compute conditional probabilities.",
            "human_time" : humantime::format_duration(timing.alignment).to_string(),
            "seconds" : timing.alignment.as_secs_f64()
        },
        "coverage_model_time": {
            "human_time": humantime::format_duration(timing.coverage).to_string(),
            "seconds": timing.coverage.as_secs_f64()
        },
        "em_time": {
            "human_time": humantime::format_duration(timing.em).to_string(),
            "seconds": timing.em.as_secs_f64()
        },
        "bin_width" : args.bin_width,
        "logistic_weight": args.logistic_weight,
        "endpoint_weight": args.endpoint_weight,
        "endpoint_support_scale": args.endpoint_support_scale,
        "coverage_folds": args.coverage_folds,
        "degradation_kernel": args.degradation_kernel,
        "coverage_max_bayes_factor": args.coverage_max_bayes_factor,
        "coverage_ablation": args.coverage_ablation,
        "coverage_warmup_iterations": args.coverage_warmup_iterations,
        "coverage_abundance_midpoint_per_million": args.coverage_abundance_midpoint_per_million,
        "sequencing_technology": args.seq_tech,
        "coverage_diagnostics": coverage_diagnostics,
        "censoring_model": args.censoring_model,
        "candidate_pruning": args.candidate_pruning,
        "dominance_bayes_factor": args.dominance_bayes_factor,
        "rank_blend": args.rank_blend,
        "rank_blend_floor": args.rank_blend_floor,
        "alignment_calibration": args.alignment_calibration,
        "filter_options" : &emi.eq_map.filter_opts,
        "discard_table" : &emi.eq_map.discard_table,
        "alignments": &args.alignments,
        "output": args.output.as_ref().expect("present"),
        "verbose": &args.verbose,
        "single_cell": &args.single_cell,
        "quiet": &args.quiet,
        "em_max_iter": &args.max_em_iter,
        "em_convergence_thresh": &args.convergence_thresh,
        "em_accel": &args.em_accel,
        "em_evaluations": em_result.evaluations,
        "em_converged": em_result.converged,
        "threads": &args.threads,
        "filter_group": &args.filter_group,
        "write_assignment_probs": &emi.eq_map.filter_opts.write_assignment_probs_type,
        "write_coverage_signals": args.write_coverage_signals,
        "coverage_signal_sample_rate": args.coverage_signal_sample_rate,
        "short_quant": &args.short_quant,
        "num_bootstraps": &args.num_bootstraps,
        "digest": seqcol_digest.to_json()
    })
}

#[allow(clippy::too_many_arguments)]
/// Connected components of transcripts linked by shared ambiguous reads.
///
/// Returns each transcript's component representative. Transcripts never seen
/// on a multi-candidate read are their own singleton component.
fn ambiguity_components<'a>(
    groups: impl Iterator<Item = &'a [AlnInfo]>,
    num_txps: usize,
) -> Vec<usize> {
    let mut parent: Vec<usize> = (0..num_txps).collect();

    fn find(parent: &mut [usize], mut node: usize) -> usize {
        while parent[node] != node {
            // Path halving keeps this near-constant without recursion.
            parent[node] = parent[parent[node]];
            node = parent[node];
        }
        node
    }

    for alns in groups {
        if alns.len() < 2 {
            continue;
        }
        let first = find(&mut parent, alns[0].ref_id as usize);
        for aln in &alns[1..] {
            let other = find(&mut parent, aln.ref_id as usize);
            if other != first {
                parent[other] = first;
            }
        }
    }
    (0..num_txps).map(|i| find(&mut parent, i)).collect()
}

fn perform_inference_and_write_output(
    header: &noodles_sam::header::Header,
    store: &mut InMemoryAlignmentStore,
    name_vec: Option<SwapVec<String>>,
    txps: &mut [TranscriptInfo],
    txps_name: &[String],
    seqcol_digest: NamedDigestVec,
    aln_time: std::time::Duration,
    args: &Args,
) -> anyhow::Result<()> {
    // print discard table information in which the user might be interested.
    info!("\ndiscard_table: \n{}\n", store.discard_table.to_table());

    // if we are using the KDE, create that here.
    let kde_opt: Option<kders::kde::KDEModel> = if args.use_kde {
        Some(kde_utils::get_kde_model(txps, store)?)
    } else {
        None
    };

    // Report tail detection so the rate can be compared against the offline
    // measurement before any modelling depends on it.
    {
        let with_tail = store.polya_tail.iter().filter(|t| **t > 0).count();
        let reads = store.len();
        if reads > 0 {
            info!(
                reads,
                with_tail,
                percent = format!("{:.1}", 100.0 * with_tail as f64 / reads as f64),
                "poly(A) tail detection"
            );
        }
    }
    let coverage_start = std::time::Instant::now();
    let mut warmup_abundances = None;
    let mut physical_creation_guard = None;
    let mut physical_warmup_diagnostics = None;
    let mut coverage_diagnostics = json!({});
    // Set before the dispatch chain, which reassigns `coverage_diagnostics`.
    let mut responsibility_diagnostics: Option<serde_json::Value> = None;
    let automatic_alignment_calibration = args.genome.is_none()
        && args.genome_alignments.is_none()
        && args.coverage_model == crate::prog_opts::CoverageModel::Auto
        && args.score_prob_denom.is_none();
    let alignment_calibration_diagnostics = if args.alignment_calibration
        == crate::prog_opts::AlignmentCalibration::Agreement
        || (args.alignment_calibration == crate::prog_opts::AlignmentCalibration::Auto
            && automatic_alignment_calibration)
    {
        let diagnostics =
            crate::util::alignment_calibration::apply_agreement_calibration(store, txps);
        info!(?diagnostics, "calibrated alignment-score likelihoods");
        Some(diagnostics)
    } else {
        None
    };
    // Optionally rebuild the per-transcript logistic coverage profile from
    // coverage-free abundance responsibilities. The default construction counts
    // every retained alignment with weight one, so an ambiguous read deposits a
    // full read of coverage on every candidate it touches (measured 5.09x
    // inflation on H69 cDNA). This reproduces the candidate rejected on
    // 2026-07-21, whose rejection rested on matched-Illumina comparator data.
    if args.coverage_ablation == crate::prog_opts::CoverageAblation::ResponsibilityProfiles {
        let start = std::time::Instant::now();
        let previous_model_coverage = store.filter_opts.model_coverage;
        store.filter_opts.model_coverage = false;
        let warmup_info = EMInfo {
            eq_map: store,
            txp_info: txps,
            max_iter: args.coverage_warmup_iterations.max(2),
            convergence_thresh: args.convergence_thresh,
            accel: crate::prog_opts::EmAccel::Squarem,
            init_abundances: None,
            kde_model: None,
            novel_locus: Vec::new(),
            novel_loci: 0,
            novel_odds_per_miss: 1.0,
        };
        let warm = if args.threads > 4 {
            em::em_par(&warmup_info, args.threads)
        } else {
            em::em(&warmup_info, args.threads)
        };
        store.filter_opts.model_coverage = previous_model_coverage;

        for txp in txps.iter_mut() {
            txp.clear_coverage_dist();
        }
        let mut reweighted = 0usize;
        for (alns, as_probs, _) in store.iter() {
            let denom: f64 = alns
                .iter()
                .zip(as_probs)
                .map(|(a, p)| warm.counts[a.ref_id as usize] * (*p as f64))
                .sum();
            if !(denom > 0.0) {
                // No abundance support: fall back to uniform mass over the
                // read's candidates rather than dropping it from the profile.
                let uniform = 1.0 / alns.len().max(1) as f64;
                for a in alns {
                    txps[a.ref_id as usize].add_interval(a.start, a.end, uniform);
                }
                continue;
            }
            for (a, p) in alns.iter().zip(as_probs) {
                let weight = warm.counts[a.ref_id as usize] * (*p as f64) / denom;
                txps[a.ref_id as usize].add_interval(a.start, a.end, weight);
            }
            reweighted += 1;
        }
        responsibility_diagnostics = Some(json!({
            "warmup_evaluations": warm.evaluations,
            "warmup_converged": warm.converged,
            "reweighted_reads": reweighted,
            "seconds": start.elapsed().as_secs_f64(),
        }));
        warmup_abundances = Some(warm.counts);
    }
    if args.coverage_model == crate::prog_opts::CoverageModel::Logistic {
        //obtaining the Cumulative Distribution Function (CDF) for each transcript
        logistic_prob(txps, args.growth_rate, &args.bin_width, args.threads);
        //Normalize the probabilities for the records of each read
        normalize_read_probs(store, txps, &args.bin_width);
    } else if args.coverage_model == crate::prog_opts::CoverageModel::Endpoint {
        crate::util::endpoint_probability::apply_endpoint_probabilities(store, txps);
        // Endpoint training does not need the legacy per-transcript coverage
        // bins, so this switch is deliberately enabled only after parsing.  It
        // tells the EM to consume the probabilities populated above.
        store.filter_opts.model_coverage = true;
    } else if args.coverage_model == crate::prog_opts::CoverageModel::Hybrid {
        logistic_prob(txps, args.growth_rate, &args.bin_width, args.threads);
        normalize_read_probs(store, txps, &args.bin_width);
        let endpoint = crate::util::endpoint_probability::endpoint_probabilities(
            store,
            txps,
            args.endpoint_support_scale,
        );
        info!(
            training_reads = endpoint.training_reads,
            logistic_weight = args.logistic_weight,
            endpoint_weight = args.endpoint_weight,
            endpoint_support_scale = args.endpoint_support_scale,
            "combining logistic and endpoint coverage models"
        );
        crate::util::hybrid_probability::apply_hybrid_probabilities(
            store,
            endpoint,
            args.logistic_weight,
            args.endpoint_weight,
        );
    } else if args.coverage_model == crate::prog_opts::CoverageModel::Adaptive
        || (args.coverage_model == crate::prog_opts::CoverageModel::Auto
            && !matches!(
                args.seq_tech,
                Some(crate::prog_opts::SequencingTech::OntDRNA)
            ))
    {
        let effective_ablation = if args.coverage_model == crate::prog_opts::CoverageModel::Auto
            && matches!(
                args.seq_tech,
                Some(crate::prog_opts::SequencingTech::PacBio)
                    | Some(crate::prog_opts::SequencingTech::PacBioHifi)
            )
            // The endpoint-geometry variants correct the fractional grid, which
            // PacBio `auto` does not use. Keeping them in this arm preserves the
            // promoted physical kernel so PacBio panels act as invariance
            // controls: a nonzero PacBio delta means this guard is wrong, not
            // that the geometry changed anything.
            && matches!(
                args.coverage_ablation,
                crate::prog_opts::CoverageAblation::Full
                    | crate::prog_opts::CoverageAblation::EndpointFeasibleSmoothing
                    | crate::prog_opts::CoverageAblation::EndpointFeasiblePrior
                    | crate::prog_opts::CoverageAblation::EndpointFeasibleGeometry
                    | crate::prog_opts::CoverageAblation::EndpointWidePriorGrid
                    | crate::prog_opts::CoverageAblation::EndpointNtMeasure
                    | crate::prog_opts::CoverageAblation::EndpointNtMeasureGeometry
                    | crate::prog_opts::CoverageAblation::ContinuousNestedGuard
                    | crate::prog_opts::CoverageAblation::ResponsibilityProfiles
                    | crate::prog_opts::CoverageAblation::ComponentMassConservation
                    | crate::prog_opts::CoverageAblation::PolyaThreePrime
            )
        {
            crate::prog_opts::CoverageAblation::PacbioPhysicalEndpoint
        } else {
            args.coverage_ablation
        };
        logistic_prob(txps, args.growth_rate, &args.bin_width, args.threads);
        normalize_read_probs(store, txps, &args.bin_width);
        let (endpoint, physical_diagnostics) = if effective_ablation
            == crate::prog_opts::CoverageAblation::PacbioPhysicalEndpoint
            && matches!(
                args.seq_tech,
                Some(crate::prog_opts::SequencingTech::PacBio)
                    | Some(crate::prog_opts::SequencingTech::PacBioHifi)
            ) {
            let warmup_start = std::time::Instant::now();
            let logistic_probabilities = store.coverage_probabilities.clone();
            let baseline_endpoint =
                crate::util::endpoint_probability::adaptive_endpoint_probabilities(
                    store,
                    txps,
                    args.coverage_folds,
                    args.endpoint_support_scale,
                    crate::prog_opts::CoverageAblation::NoEndpoint,
                );
            crate::util::hybrid_probability::apply_adaptive_probabilities(
                store,
                baseline_endpoint,
                args.logistic_weight,
                args.endpoint_weight,
                args.coverage_max_bayes_factor,
                crate::prog_opts::CoverageAblation::NoEndpoint,
            );
            let warmup_info = EMInfo {
                eq_map: store,
                txp_info: txps,
                max_iter: args.max_em_iter,
                convergence_thresh: args.convergence_thresh,
                // The guarded PacBio anchor commonly reaches the evaluation
                // limit. Ordinary EM has the cheapest evaluation and matches
                // the accelerated solutions in this regime.
                accel: crate::prog_opts::EmAccel::None,
                init_abundances: None,
                kde_model: None,
                novel_locus: Vec::new(),
                novel_loci: 0,
                novel_odds_per_miss: 1.0,
            };
            let warmup_result = if args.threads > 1 && store.len() >= 100_000 {
                em::em_par(&warmup_info, args.threads)
            } else {
                em::em(&warmup_info, args.threads)
            };
            physical_warmup_diagnostics = Some(json!({
                "kernel": "pacbio_auto_no_endpoint_baseline",
                "evaluations": warmup_result.evaluations,
                "converged": warmup_result.converged,
                "seconds": warmup_start.elapsed().as_secs_f64(),
            }));
            warmup_abundances = Some(warmup_result.counts);
            store.coverage_probabilities = logistic_probabilities;
            let (endpoint, diagnostics, protected) =
                crate::util::pacbio_endpoint_probability::physical_endpoint_probabilities(
                    store, txps,
                );
            physical_creation_guard = Some(protected);
            (endpoint, Some(diagnostics))
        } else {
            (
                crate::util::endpoint_probability::adaptive_endpoint_probabilities(
                    store,
                    txps,
                    args.coverage_folds,
                    args.endpoint_support_scale,
                    effective_ablation,
                ),
                None,
            )
        };
        let diagnostics = crate::util::hybrid_probability::apply_adaptive_probabilities(
            store,
            endpoint,
            args.logistic_weight,
            args.endpoint_weight,
            args.coverage_max_bayes_factor,
            effective_ablation,
        );
        info!(?diagnostics, "applied adaptive coverage model");
        coverage_diagnostics = if args.coverage_model == crate::prog_opts::CoverageModel::Auto {
            let technology_kernel = if matches!(
                args.seq_tech,
                Some(crate::prog_opts::SequencingTech::PacBio)
                    | Some(crate::prog_opts::SequencingTech::PacBioHifi)
            ) {
                if effective_ablation == crate::prog_opts::CoverageAblation::PacbioPhysicalEndpoint
                {
                    "pacbio_guarded_physical_endpoint"
                } else {
                    "pacbio_adaptive_logistic"
                }
            } else {
                "cross_fitted_adaptive_endpoint"
            };
            json!({
                "technology_kernel": technology_kernel,
                "effective_ablation": effective_ablation,
                "adaptive": diagnostics,
                "physical_endpoint": physical_diagnostics,
            })
        } else {
            serde_json::to_value(diagnostics)?
        };
    } else if args.coverage_model == crate::prog_opts::CoverageModel::Degradation
        || (args.coverage_model == crate::prog_opts::CoverageModel::Auto
            && matches!(
                args.seq_tech,
                Some(crate::prog_opts::SequencingTech::OntDRNA)
            ))
    {
        logistic_prob(txps, args.growth_rate, &args.bin_width, args.threads);
        normalize_read_probs(store, txps, &args.bin_width);
        let (endpoint, degradation_diagnostics) =
            crate::util::degradation_probability::degradation_adjusted_endpoint_probabilities(
                store,
                txps,
                args.coverage_folds,
                args.endpoint_support_scale,
                args.degradation_kernel,
                args.coverage_ablation,
            );
        let adaptive_diagnostics = crate::util::hybrid_probability::apply_adaptive_probabilities(
            store,
            endpoint,
            args.logistic_weight,
            args.endpoint_weight,
            args.coverage_max_bayes_factor,
            args.coverage_ablation,
        );
        info!(
            ?degradation_diagnostics,
            ?adaptive_diagnostics,
            "applied degradation-aware coverage model"
        );
        coverage_diagnostics = json!({
            "technology_kernel": "ont_drna_competing_risks",
            "degradation": degradation_diagnostics,
            "adaptive": adaptive_diagnostics,
        });
    }
    if args.coverage_ablation == crate::prog_opts::CoverageAblation::PolyaThreePrime {
        let diagnostics = crate::util::polya_probability::apply_polya_probabilities(store, txps);
        info!(?diagnostics, "applied poly(A) 3'-completeness likelihood");
        coverage_diagnostics["polya_three_prime"] = serde_json::to_value(diagnostics)?;
    }
    if let Some(diagnostics) = responsibility_diagnostics {
        coverage_diagnostics["responsibility_profiles"] = diagnostics;
    }
    if let Some(diagnostics) = alignment_calibration_diagnostics {
        coverage_diagnostics["alignment_calibration"] = serde_json::to_value(diagnostics)?;
    }
    let transcriptome_input = args.genome.is_none() && args.genome_alignments.is_none();
    let automatic_inference =
        transcriptome_input && args.coverage_model == crate::prog_opts::CoverageModel::Auto;
    const AUTO_RANK_BLEND_MIN_CENSOR_SCALE_NT: f64 = 37.0;
    const AUTO_RANK_BLEND_MIN_TRANSCRIPTS: usize = 1_000;
    let rank_blend_scale =
        if args.rank_blend == crate::prog_opts::RankBlend::Auto && automatic_inference {
            Some(crate::util::censoring_probability::estimate_censor_scale(
                store, txps,
            ))
        } else {
            None
        };
    let rank_blend_active = args.rank_blend == crate::prog_opts::RankBlend::Fixed
        || rank_blend_scale.is_some_and(|scale| {
            scale > AUTO_RANK_BLEND_MIN_CENSOR_SCALE_NT
                && txps.len() >= AUTO_RANK_BLEND_MIN_TRANSCRIPTS
        });
    coverage_diagnostics["rank_blend_selection"] = json!({
        "requested": args.rank_blend,
        "active": rank_blend_active,
        "learned_censor_scale_nt": rank_blend_scale,
        "minimum_censor_scale_nt": AUTO_RANK_BLEND_MIN_CENSOR_SCALE_NT,
        "minimum_transcripts": AUTO_RANK_BLEND_MIN_TRANSCRIPTS,
    });
    if (matches!(
        args.coverage_ablation,
        crate::prog_opts::CoverageAblation::AbundanceBlend
    ) || rank_blend_active)
        && warmup_abundances.is_none()
    {
        let warmup_start = std::time::Instant::now();
        let previous_model_coverage = store.filter_opts.model_coverage;
        store.filter_opts.model_coverage = false;
        let warmup_info = EMInfo {
            eq_map: store,
            txp_info: txps,
            max_iter: args.coverage_warmup_iterations.max(2),
            convergence_thresh: args.convergence_thresh,
            // SQUAREM is both faster and substantially less path-dependent
            // than DAAREM for the coverage-free warm start on highly
            // non-identifiable equivalence classes.
            accel: crate::prog_opts::EmAccel::Squarem,
            init_abundances: None,
            kde_model: None,
            novel_locus: Vec::new(),
            novel_loci: 0,
            novel_odds_per_miss: 1.0,
        };
        let warmup_result = if args.threads > 4 {
            em::em_par(&warmup_info, args.threads)
        } else {
            em::em(&warmup_info, args.threads)
        };
        store.filter_opts.model_coverage = previous_model_coverage;
        let warmup_seconds = warmup_start.elapsed().as_secs_f64();
        coverage_diagnostics["abundance_warmup"] = json!({
            "evaluations": warmup_result.evaluations,
            "converged": warmup_result.converged,
            "seconds": warmup_seconds,
        });
        warmup_abundances = Some(warmup_result.counts);
    }
    if let Some(diagnostics) = physical_warmup_diagnostics {
        coverage_diagnostics["abundance_warmup"] = diagnostics;
    }
    if args.censoring_model == crate::prog_opts::CensoringModel::Adaptive
        || (args.censoring_model == crate::prog_opts::CensoringModel::Auto && automatic_inference)
    {
        let diagnostics =
            crate::util::censoring_probability::apply_censoring_probabilities(store, txps);
        info!(?diagnostics, "applied read-level censoring likelihood");
        coverage_diagnostics["censoring"] = serde_json::to_value(diagnostics)?;
    }
    if args.candidate_pruning == crate::prog_opts::CandidatePruning::Dominance
        || (args.candidate_pruning == crate::prog_opts::CandidatePruning::Auto
            && automatic_inference)
    {
        let diagnostics = crate::util::dominance_pruning::apply_dominance_pruning(
            store,
            args.dominance_bayes_factor,
        );
        info!(?diagnostics, "applied candidate-dominance pruning");
        coverage_diagnostics["candidate_pruning"] = serde_json::to_value(diagnostics)?;
    }
    let coverage_time = coverage_start.elapsed();

    if store.junc_informative_reads > 0 {
        info!(
            "junction evidence: reads spanning >=1 internal boundary: {}; of those, disagreeing with EVERY candidate: {} ({:.2}%) -- candidate unannotated-isoform mass",
            store.junc_informative_reads.to_formatted_string(&Locale::en),
            store.junc_all_mismatch_reads.to_formatted_string(&Locale::en),
            100.0 * store.junc_all_mismatch_reads as f64 / store.junc_informative_reads as f64,
        );
    }
    info!(
        "Total number of alignment records : {}",
        store.total_len().to_formatted_string(&Locale::en)
    );
    info!(
        "number of aligned reads : {}",
        store.num_aligned_reads().to_formatted_string(&Locale::en)
    );
    info!(
        "number of unique alignments : {}",
        store.unique_alignments().to_formatted_string(&Locale::en)
    );

    // if we are seeding the quantification estimates with short read
    // abundances, then read those in here.
    let init_abundances = args
        .short_quant
        .as_ref()
        .map(|sr_path| read_short_quant_vec(sr_path, txps_name).unwrap_or_else(|e| panic!("{}", e)))
        .or_else(|| {
            if physical_creation_guard.is_none() {
                warmup_abundances.clone()
            } else {
                None
            }
        });
    let blend_reference = if matches!(
        args.coverage_ablation,
        crate::prog_opts::CoverageAblation::AbundanceBlend
    ) || physical_creation_guard.is_some()
        || rank_blend_active
    {
        warmup_abundances
    } else {
        None
    };

    // Reads whose splice structure disagrees with every annotated candidate get
    // a novel latent state, shared by all such reads at the same locus. The
    // locus is the read's ambiguity component, so no gene annotation is needed.
    let mut novel_locus: Vec<i32> = Vec::new();
    let mut novel_loci = 0usize;
    // locus id -> (member transcript ids, flagged read count), for the report.
    let mut novel_members: Vec<(Vec<u32>, usize)> = Vec::new();
    if args.models_unannotated_isoforms() && !store.min_junc_misses.is_empty() {
        let components =
            ambiguity_components(store.iter().map(|(alns, _, _)| alns), txps.len());
        let mut dense: std::collections::HashMap<usize, i32> = std::collections::HashMap::new();
        novel_locus = vec![-1; store.len()];
        for read_index in 0..store.len() {
            if store.min_junc_misses.get(read_index).copied().unwrap_or(0)
                < args.novel_min_misses.max(1)
            {
                continue;
            }
            if args.novel_require_hits
                && store.max_junc_hits.get(read_index).copied().unwrap_or(0) <= 0
            {
                continue;
            }
            let start = store.boundaries[read_index];
            let end = store.boundaries[read_index + 1];
            if end == start {
                continue;
            }
            let root = components[store.alignments[start].ref_id as usize];
            let next = dense.len() as i32;
            let id = *dense.entry(root).or_insert(next);
            novel_locus[read_index] = id;
        }
        // Drop loci that accumulated too little evidence, and re-densify. A
        // handful of disagreeing reads at a locus is alignment noise; a real
        // unannotated isoform accumulates many.
        if args.novel_min_locus_reads > 1 {
            let mut per_locus = vec![0usize; dense.len()];
            for &id in &novel_locus {
                if id >= 0 {
                    per_locus[id as usize] += 1;
                }
            }
            let mut remap = vec![-1i32; dense.len()];
            let mut next = 0i32;
            for (old_id, count) in per_locus.iter().enumerate() {
                if *count >= args.novel_min_locus_reads {
                    remap[old_id] = next;
                    next += 1;
                }
            }
            for id in novel_locus.iter_mut() {
                if *id >= 0 {
                    *id = remap[*id as usize];
                }
            }
            // Retain must be followed by remapping the surviving values: they
            // still hold pre-filter ids, which would index past the compacted
            // locus table.
            dense.retain(|_, v| remap[*v as usize] >= 0);
            for v in dense.values_mut() {
                *v = remap[*v as usize];
            }
            novel_loci = next as usize;
        } else {
            novel_loci = dense.len();
        }
        // Members of each surviving locus, so the report can name it. There is
        // no gene label to key on: bramble carries transcript names only, and a
        // genuinely novel locus has no gene by definition.
        novel_members = vec![(Vec::new(), 0usize); novel_loci];
        let mut root_of_locus: std::collections::HashMap<i32, usize> =
            std::collections::HashMap::new();
        for (root, id) in dense.iter() {
            root_of_locus.insert(*id, *root);
        }
        let mut by_root: std::collections::HashMap<usize, usize> =
            std::collections::HashMap::new();
        for (id, root) in root_of_locus.iter() {
            by_root.insert(*root, *id as usize);
        }
        for (txp, root) in components.iter().enumerate() {
            if let Some(&id) = by_root.get(root) {
                novel_members[id].0.push(txp as u32);
            }
        }
        for &id in novel_locus.iter() {
            if id >= 0 {
                novel_members[id as usize].1 += 1;
            }
        }
        let flagged = novel_locus.iter().filter(|v| **v >= 0).count();
        info!(
            novel_loci,
            flagged_reads = flagged,
            odds_per_miss = args.novel_odds_per_miss,
            min_misses = args.novel_min_misses,
            min_locus_reads = args.novel_min_locus_reads,
            require_hits = args.novel_require_hits,
            "annotation-omission model: novel latent states created"
        );
        coverage_diagnostics["annotation_omission"] = json!({
            "novel_loci": novel_loci,
            "flagged_reads": flagged,
            "odds_per_miss": args.novel_odds_per_miss,
        });
    }

    // wrap up all of the relevant information we need for estimation
    // in an EMInfo struct and then call the EM algorithm.
    let emi = EMInfo {
        eq_map: store,
        txp_info: txps,
        max_iter: args.max_em_iter,
        convergence_thresh: args.convergence_thresh,
        accel: args.em_accel,
        init_abundances,
        kde_model: kde_opt,
        novel_locus,
        novel_loci,
        novel_odds_per_miss: args.novel_odds_per_miss,
    };

    if args.use_kde {
        /*
        // run EM for model train iterations
        let orig_iter = emi.max_iter;
        emi.max_iter = 10;
        let counts = em::em(&emi, args.threads);
        // relearn the kde
        let new_model =
        kde_utils::refresh_kde_model(&txps, &store, &emi.kde_model.unwrap(), &counts);
        info!("refreshed KDE model");
        emi.kde_model = Some(new_model?);
        emi.max_iter = orig_iter;
        */
    }

    let em_start = std::time::Instant::now();
    // The parallel M-step does not implement the novel term.
    let use_parallel_em = novel_loci == 0
        && args.threads > 4
        || (physical_creation_guard.is_some() && args.threads > 1 && store.len() >= 100_000);
    let mut em_result = if use_parallel_em {
        em::em_par(&emi, args.threads)
    } else {
        em::em(&emi, args.threads)
    };
    let em_time = em_start.elapsed();
    if let Some(reference) = blend_reference {
        if let Some(protected) = physical_creation_guard {
            let additive = (store.num_aligned_reads() as f64 * 10.0 / 1_000_000.0).max(1.0);
            let minimum_extreme_count = store.num_aligned_reads() as f64 * 0.0025;
            let continuous = args.coverage_ablation
                == crate::prog_opts::CoverageAblation::ContinuousNestedGuard;
            // Continuous variant: instead of clamping the rare extreme case,
            // shrink every nested-candidate increase toward the coverage-free
            // anchor by `b / (b + m)`, with `m` the abundance-blend midpoint.
            let midpoint = (store.num_aligned_reads() as f64
                * args.coverage_abundance_midpoint_per_million
                / 1_000_000.0)
                .max(1.0);
            let mut clamped = 0usize;
            for ((count, &baseline), &is_protected) in
                em_result.counts.iter_mut().zip(&reference).zip(&protected)
            {
                if !is_protected {
                    continue;
                }
                if continuous {
                    if *count > baseline {
                        let anchor = baseline.max(0.0);
                        let gate = anchor / (anchor + midpoint);
                        *count = anchor + gate * (*count - anchor);
                        clamped += 1;
                    }
                } else {
                    let maximum = baseline.max(0.0) * 100.0 + additive;
                    if *count > maximum && *count > minimum_extreme_count {
                        *count = maximum;
                        clamped += 1;
                    }
                }
            }
            coverage_diagnostics["physical_creation_guard"] = json!({
                "protected_transcripts": protected.iter().filter(|&&value| value).count(),
                "clamped_transcripts": clamped,
                "mode": if args.coverage_ablation == crate::prog_opts::CoverageAblation::ContinuousNestedGuard { "continuous" } else { "hard_clamp" },
                "maximum_fold_creation": 100.0,
                "additive_cpm": 10.0,
                "minimum_extreme_library_fraction": 0.0025,
            });
        }
        if rank_blend_active
            || matches!(
                args.coverage_ablation,
                crate::prog_opts::CoverageAblation::AbundanceBlend
            )
        {
            let midpoint = (store.num_aligned_reads() as f64
                * args.coverage_abundance_midpoint_per_million
                / 1_000_000.0)
                .max(1.0);
            let conserve_components = args.coverage_ablation
                == crate::prog_opts::CoverageAblation::ComponentMassConservation;
            // Snapshot the corrected abundances so each ambiguity component's
            // total can be restored after the blend shrinks its members toward
            // the coverage-free anchor.
            let pre_blend = conserve_components.then(|| em_result.counts.clone());
            let mut gate_sum = 0.0;
            for (count, &baseline) in em_result.counts.iter_mut().zip(&reference) {
                let ratio = baseline.max(0.0) / midpoint;
                let ratio4 = ratio * ratio * ratio * ratio;
                let gate =
                    args.rank_blend_floor + (1.0 - args.rank_blend_floor) * ratio4 / (1.0 + ratio4);
                *count = baseline + gate * (*count - baseline);
                gate_sum += gate;
            }
            coverage_diagnostics["abundance_blend"] = json!({
                "midpoint_count": midpoint,
                "mean_gate": gate_sum / em_result.counts.len().max(1) as f64,
                "minimum_gate": args.rank_blend_floor,
            });
            if let Some(pre_blend) = pre_blend {
                let components =
                    ambiguity_components(store.iter().map(|(alns, _, _)| alns), em_result.counts.len());
                let minimum_mass = store.num_aligned_reads() as f64 * 0.0025;
                let mut before = vec![0.0f64; em_result.counts.len()];
                let mut after = vec![0.0f64; em_result.counts.len()];
                let mut members = vec![0usize; em_result.counts.len()];
                for (index, &root) in components.iter().enumerate() {
                    before[root] += pre_blend[index];
                    after[root] += em_result.counts[index];
                    members[root] += 1;
                }
                // Singleton components carry no ambiguity, and rescaling them
                // simply undoes the blend; the broad variant of this candidate
                // lost rank accuracy that way.
                let mut conserved = 0usize;
                for (index, &root) in components.iter().enumerate() {
                    if members[root] > 1 && before[root] >= minimum_mass && after[root] > 0.0 {
                        em_result.counts[index] *= before[root] / after[root];
                        if index == root {
                            conserved += 1;
                        }
                    }
                }
                coverage_diagnostics["component_mass_conservation"] = json!({
                    "conserved_components": conserved,
                    "minimum_component_mass": minimum_mass,
                });
            }
        }
        let total: f64 = em_result.counts.iter().sum();
        if total > 0.0 {
            let scale = store.num_aligned_reads() as f64 / total;
            em_result
                .counts
                .iter_mut()
                .for_each(|count| *count *= scale);
        }
    }
    // Per-locus unexplained mass. Keyed by ambiguity component rather than by
    // gene: bramble exposes transcript names only, and a genuinely unannotated
    // locus has no gene label by construction. A gene is reported when one can
    // be recovered from the transcript name (e.g. GENCODE's pipe-delimited
    // header), otherwise "." -- the member transcripts are the real identifier.
    if !novel_members.is_empty() {
        let mut path = args.output.as_ref().expect("output prefix").clone();
        let mut name = path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default();
        name.push_str(".unexplained.tsv");
        path.set_file_name(name);
        let mut w = std::io::BufWriter::new(std::fs::File::create(&path)?);
        use std::io::Write;
        writeln!(
            w,
            "locus\tgene\tn_transcripts\tflagged_reads\tunprojectable_reads\tunexplained_mass\tannotated_mass\tunexplained_fraction\ttranscripts"
        )?;
        let n_txps = txps.len();
        let mut reported = 0usize;
        for (id, (members, flagged)) in novel_members.iter().enumerate() {
            let unexplained = em_result.counts.get(n_txps + id).copied().unwrap_or(0.0);
            if unexplained <= 0.0 {
                continue;
            }
            let annotated: f64 = members
                .iter()
                .map(|t| em_result.counts.get(*t as usize).copied().unwrap_or(0.0))
                .sum();
            // A gene label if the reference name carries one; otherwise the
            // member transcripts identify the locus.
            let gene = members
                .first()
                .and_then(|t| txps_name.get(*t as usize))
                .and_then(|n| n.split('|').nth(1))
                .filter(|g| g.starts_with("ENSG") || g.starts_with("gene"))
                .unwrap_or(".");
            let mut names: Vec<&str> = members
                .iter()
                .filter_map(|t| txps_name.get(*t as usize))
                .map(|n| n.split('|').next().unwrap_or(n.as_str()))
                .take(8)
                .collect();
            if members.len() > names.len() {
                names.push("...");
            }
            // Unprojectable reads overlapping this locus. A read overlapping
            // several member transcripts increments each of them, so the max is
            // a lower bound on distinct reads while the sum would double-count;
            // report the lower bound.
            let unprojectable = members
                .iter()
                .filter_map(|t| store.unprojectable_per_txp.get(*t as usize))
                .copied()
                .max()
                .unwrap_or(0);
            writeln!(
                w,
                "{}\t{}\t{}\t{}\t{}\t{:.3}\t{:.3}\t{:.4}\t{}",
                id,
                gene,
                members.len(),
                flagged,
                unprojectable,
                unexplained,
                annotated,
                unexplained / (unexplained + annotated).max(f64::MIN_POSITIVE),
                names.join(",")
            )?;
            reported += 1;
        }
        info!(
            loci = reported,
            path = %path.display(),
            "wrote per-locus unexplained-mass report"
        );
    }

    // Per-transcript unprojectable overlap counts. Unlike the locus report above
    // this covers every annotated transcript with overlapping unprojectable
    // reads, including loci that produced no junction-mismatch evidence and so
    // have no novel state -- the mass there is invisible to the locus report.
    if store.unprojectable_per_txp.iter().any(|c| *c > 0) {
        let mut path = args.output.as_ref().expect("output prefix").clone();
        let mut name = path
            .file_name()
            .map(|n| n.to_string_lossy().into_owned())
            .unwrap_or_default();
        name.push_str(".unprojectable.tsv");
        path.set_file_name(name);
        let mut w = std::io::BufWriter::new(std::fs::File::create(&path)?);
        use std::io::Write;
        writeln!(w, "tname\tunprojectable_overlapping_reads")?;
        let mut rows = 0usize;
        for (tid, count) in store.unprojectable_per_txp.iter().enumerate() {
            if *count == 0 {
                continue;
            }
            let name = txps_name
                .get(tid)
                .map(|n| n.split('|').next().unwrap_or(n.as_str()))
                .unwrap_or(".");
            writeln!(w, "{}\t{}", name, count)?;
            rows += 1;
        }
        info!(
            transcripts = rows,
            path = %path.display(),
            "wrote per-transcript unprojectable-overlap report"
        );
    }

    let counts = &em_result.counts;

    let aux_txp_counts = crate::util::aux_counts::get_aux_counts(store, txps)?;

    // prepare the JSON object we'll write
    // to meta_info.json
    let json_info = get_json_info(
        args,
        &emi,
        &em_result,
        &seqcol_digest,
        &RunTiming {
            alignment: aln_time,
            coverage: coverage_time,
            em: em_time,
        },
        &coverage_diagnostics,
    );

    // write the output
    write_output(
        args.output.as_ref().expect("present"),
        json_info,
        header,
        counts,
        &aux_txp_counts,
    )?;

    // if the user requested bootstrap replicates,
    // compute and write those out now.
    if args.num_bootstraps > 0 {
        let breps = em::bootstrap(&emi, args.num_bootstraps, args.threads);

        let mut new_arrays = vec![];
        let mut bs_fields = vec![];
        for (i, b) in breps.into_iter().enumerate() {
            let bs_array = Float64Array::from_vec(b);
            bs_fields.push(Field::new(
                format!("bootstrap.{}", i),
                bs_array.data_type().clone(),
                false,
            ));
            new_arrays.push(bs_array.boxed());
        }
        let chunk = Chunk::new(new_arrays);
        write_infrep_file(args.output.as_ref().expect("present"), bs_fields, chunk)?;
    }

    if args.write_coverage_signals {
        let name_vec =
            name_vec.expect("cannot write coverage signals without a valid vector of read names");
        write_coverage_signals(
            args.output.as_ref().expect("present"),
            &emi,
            name_vec,
            txps_name,
            args.coverage_signal_sample_rate,
        )?;
    } else if args.write_assignment_probs.is_some() {
        let name_vec = name_vec
            .expect("cannot write assignment probabilities without valid vector of read names");
        write_out_prob(
            args.output.as_ref().expect("present"),
            &emi,
            counts,
            name_vec,
            txps_name,
            args.display_thresh,
        )?;
    }

    Ok(())
}

pub fn quantify_bulk_alignments_from_bam<R: BufRead>(
    header: &noodles_sam::Header,
    filter_opts: AlignmentFilters,
    reader: &mut bam::io::Reader<R>,
    txps: &mut [TranscriptInfo],
    txps_name: &[String],
    args: &Args,
    seqcol_digest: NamedDigestVec,
) -> anyhow::Result<()> {
    let mut name_vec = if filter_opts.write_assignment_probs {
        Some(SwapVec::<String>::with_config(SwapVecConfig {
            swap_after: Default::default(),
            batch_size: Default::default(),
            compression: Some(swapvec::Compression::Lz4),
        }))
    } else {
        None
    };
    // now parse the actual alignments for the reads and store the results
    // in our in-memory stor
    let mut store = InMemoryAlignmentStore::new(filter_opts, header);
    let read_aln_start = std::time::SystemTime::now();
    alignment_parser::parse_alignments(
        &mut store,
        &mut name_vec,
        header,
        reader,
        txps,
        args.sort_check_num,
        args.quiet,
    )?;
    let read_aln_time = read_aln_start.elapsed()?;
    info!(
        "Parsing of alignments from input took: {}",
        humantime::format_duration(read_aln_time).to_string()
    );

    perform_inference_and_write_output(
        header,
        &mut store,
        name_vec,
        txps,
        txps_name,
        seqcol_digest,
        read_aln_time,
        args,
    )
}

/// Genome-BAM mode: parse a name-collated, spliced genome BAM, project each
/// read's alignments onto the transcriptome with bramble, and quantify.
///
/// `txp_header` / `txps` / `txps_name` describe the *transcriptome* (built from
/// the annotation via the projection bridge); `genome_header` is the input BAM's
/// (chromosome) header used to read records and resolve reference ids.
#[allow(clippy::too_many_arguments)]
pub fn quantify_genome_alignments_from_bam<R: BufRead>(
    genome_header: &noodles_sam::Header,
    txp_header: &noodles_sam::Header,
    g2t: &G2TTree,
    proj_config: &ProjectionConfig,
    filter_opts: AlignmentFilters,
    reader: &mut bam::io::Reader<R>,
    txps: &mut [TranscriptInfo],
    txps_name: &[String],
    args: &Args,
    seqcol_digest: NamedDigestVec,
) -> anyhow::Result<()> {
    let mut name_vec = if filter_opts.write_assignment_probs {
        Some(SwapVec::<String>::with_config(SwapVecConfig {
            swap_after: Default::default(),
            batch_size: Default::default(),
            compression: Some(swapvec::Compression::Lz4),
        }))
    } else {
        None
    };

    // the store keys alignments by transcript, so it lives over the
    // transcriptome header.
    let mut store = InMemoryAlignmentStore::new(filter_opts, txp_header);
    let read_aln_start = std::time::SystemTime::now();
    alignment_parser::parse_genome_alignments(
        &mut store,
        &mut name_vec,
        genome_header,
        g2t,
        proj_config,
        args.projected_prob_beta,
        args.projected_prob_source,
        reader,
        txps,
        args.sort_check_num,
        args.quiet,
    )?;
    let read_aln_time = read_aln_start.elapsed()?;
    info!(
        "Parsing and projection of genome alignments took: {}",
        humantime::format_duration(read_aln_time).to_string()
    );

    perform_inference_and_write_output(
        txp_header,
        &mut store,
        name_vec,
        txps,
        txps_name,
        seqcol_digest,
        read_aln_time,
        args,
    )
}

/// Genome-read mode: spliced-align raw reads to the genome with the rammap
/// aligner, project each read's mappings onto the transcriptome with bramble,
/// and quantify.
///
/// This mirrors [`quantify_bulk_alignments_raw_reads`] (same producer → mapper →
/// consumer pipeline), but the mapper projects each read's genomic mappings onto
/// the transcripts (via bramble) rather than mapping directly to a transcriptome.
/// `txp_header` / `txps` / `txps_name` describe the transcriptome built from the
/// annotation; `aligner` is a spliced genome aligner; `g2t` the genome→
/// transcriptome index over the same reference order as the aligner targets.
#[allow(clippy::too_many_arguments)]
#[allow(unused_mut)]
pub fn quantify_genome_raw_reads(
    txp_header: &noodles_sam::Header,
    mut aligner: crate::util::mapper::Mapper,
    g2t: &G2TTree,
    proj_config: &ProjectionConfig,
    filter_opts: AlignmentFilters,
    read_paths: &[std::path::PathBuf],
    txps: &mut [TranscriptInfo],
    txps_name: &[String],
    args: &Args,
    seqcol_digest: NamedDigestVec,
) -> anyhow::Result<()> {
    // shared, read-only view of the transcript info (for the projected filters).
    let mut txp_info_view: Vec<TranscriptInfo> = Vec::with_capacity(txps.len());
    for ti in txps.iter() {
        txp_info_view.push(ti.clone());
    }

    let map_threads = args.threads.saturating_sub(2).max(1);

    let beta = args.projected_prob_beta;
    let use_fasta = proj_config.use_fasta;
    let prob_source = args.projected_prob_source;

    type ReadGroup = ReadChunkWithNames;
    type AlignmentGroupInfo = (
        Vec<AlnInfo>,
        Vec<f32>,
        Vec<usize>,
        Option<Vec<String>>,
        Vec<(i32, i32)>,
    );

    let (read_sender, read_receiver): (Sender<ReadGroup>, Receiver<ReadGroup>) =
        bounded(args.threads * 10);

    const READ_CHUNK_SIZE: usize = 200;
    let mut rpaths = vec![];
    read_paths.clone_into(&mut rpaths);

    // Producer thread: read sequences and send them to the channel.
    let producer = std::thread::spawn(move || {
        let mut ctr = 0_usize;
        let mut chunk_size = 0_usize;
        let mut read_chunk = ReadChunkWithNames::new();

        let mark_chunk = |chunk_size: &mut usize,
                          ctr: &mut usize,
                          read_chunk: &mut ReadGroup,
                          read_sender: &Sender<ReadGroup>| {
            *chunk_size += 1;
            *ctr += 1;
            if *chunk_size >= READ_CHUNK_SIZE {
                read_sender
                    .send(read_chunk.clone())
                    .expect("Error sending sequence");
                read_chunk.clear();
                *chunk_size = 0;
            }
        };

        for read_path in rpaths {
            match get_source_type(&read_path) {
                InputSourceType::Ubam => {
                    let mut reader = std::fs::File::open(read_path)
                        .map(bam::io::Reader::new)
                        .expect("could not create BAM reader");
                    let header = reader.read_header().expect("could not read BAM header");
                    for result in reader.record_bufs(&header) {
                        let record = result.expect("Error reading ubam record");
                        record.add_to_read_group(&mut read_chunk);
                        mark_chunk(&mut chunk_size, &mut ctr, &mut read_chunk, &read_sender);
                    }
                }
                s @ (InputSourceType::Fastx | InputSourceType::Unknown) => {
                    if matches!(s, InputSourceType::Unknown) {
                        warn!(
                            "could not determine input file type for {} from suffix; assuming (possibly gzipped) fastx",
                            read_path.display()
                        );
                    }
                    let mut reader =
                        parse_fastx_file(read_path).expect("valid path/file to read sequences");
                    while let Some(result) = reader.next() {
                        let record = result.expect("Error reading record");
                        record.add_to_read_group(&mut read_chunk);
                        mark_chunk(&mut chunk_size, &mut ctr, &mut read_chunk, &read_sender);
                    }
                }
            }
        }
        if chunk_size > 0 {
            read_sender
                .send(read_chunk)
                .expect("Error sending sequence");
        }
        ctr
    });

    let (mut store, name_vec, aln_time) = std::thread::scope(
        |s| -> anyhow::Result<(
            InMemoryAlignmentStore,
            Option<SwapVec<String>>,
            std::time::Duration,
        )> {
            const ALN_GROUP_CHUNK_LIMIT: usize = 100;

            let (aln_group_sender, aln_group_receiver): (
                Sender<AlignmentGroupInfo>,
                Receiver<AlignmentGroupInfo>,
            ) = bounded(args.threads * 100);

            // Read names are also required by the coverage-signal export; keep this
            // in sync with the `AlignmentFilters` builder in `main.rs`, which enables
            // name retention for either output.
            let write_assignment_probs: bool =
                args.write_assignment_probs.is_some() || args.write_coverage_signals;
            // Diagnostics: reads the aligner maps to the genome vs reads that
            // produce >=1 projected transcriptome alignment (projection loss).
            let n_genome_mapped = std::sync::Arc::new(std::sync::atomic::AtomicU64::new(0));
            let n_projected = std::sync::Arc::new(std::sync::atomic::AtomicU64::new(0));
            // Junction evidence: reads whose splice structure disagrees with
            // every candidate transcript are evidence that their true isoform is
            // missing from the annotation.
            let n_junc_informative = std::sync::Arc::new(std::sync::atomic::AtomicU64::new(0));
            let n_junc_all_mismatch = std::sync::Arc::new(std::sync::atomic::AtomicU64::new(0));
            // Unprojectable reads: mapped to the genome but compatible with no
            // annotated transcript's splice structure, so they are dropped before
            // quantification. Under `--model-unannotated-isoforms` we attribute
            // them to the loci whose exons they overlap so their mass is
            // reported rather than silently lost. Purely diagnostic: these reads
            // still never enter the EM, so no annotated estimate is perturbed.
            let track_unprojectable = args.models_unannotated_isoforms();
            let n_unprojectable = std::sync::Arc::new(std::sync::atomic::AtomicU64::new(0));
            let n_unprojectable_intergenic =
                std::sync::Arc::new(std::sync::atomic::AtomicU64::new(0));
            let unprojectable_per_txp: std::sync::Arc<Vec<std::sync::atomic::AtomicU32>> =
                std::sync::Arc::new(if track_unprojectable {
                    (0..txps.len())
                        .map(|_| std::sync::atomic::AtomicU32::new(0))
                        .collect()
                } else {
                    Vec::new()
                });
            // Mapper threads: align each read to the genome, then project its
            // mappings onto the transcriptome and filter.
            let consumers: Vec<_> = (0..map_threads)
                .map(|_| {
                    let receiver = read_receiver.clone();
                    let filter = filter_opts.clone();
                    let loc_aligner = aligner.clone();
                    let my_txp_info_view = &txp_info_view;
                    let aln_group_sender = aln_group_sender.clone();
                    let loc_junc_informative = n_junc_informative.clone();
                    let loc_junc_all_mismatch = n_junc_all_mismatch.clone();
                    let loc_unprojectable = n_unprojectable.clone();
                    let loc_unprojectable_intergenic = n_unprojectable_intergenic.clone();
                    let loc_unprojectable_per_txp = unprojectable_per_txp.clone();
                    let n_genome_mapped = std::sync::Arc::clone(&n_genome_mapped);
                    let n_projected = std::sync::Arc::clone(&n_projected);

                    s.spawn(move || {
                        let mut discard_table = DiscardTable::new();

                        let mut chunk_size = 0_usize;
                        let mut aln_group_alns: Vec<AlnInfo> = Vec::new();
                        let mut aln_group_probs: Vec<f32> = Vec::new();
                        let mut aln_group_boundaries: Vec<usize> = Vec::new();
                        let mut aln_group_read_names = write_assignment_probs.then(Vec::new);
                        let mut aln_group_junc: Vec<(i32, i32)> = Vec::new();
                        aln_group_boundaries.push(0);
                        // reused across all reads this worker projects (avoids
                        // per-read allocation of bramble's projection scratch).
                        let mut pctx = ProjectionContext::new();
                        // Reused across every read this worker projects, to avoid
                        // re-allocating the per-read alignment/score scratch vectors.
                        let mut galns: Vec<bramble_rs::GenomicAlignment> = Vec::new();
                        let mut src_scores: Vec<i32> = Vec::new();

                        for read_chunk in receiver {
                            for (name, seq) in read_chunk.iter() {
                                let map_res_opt =
                                    crate::util::mapper::map_read(&loc_aligner, name, seq);
                                let Ok(mappings) = map_res_opt else {
                                    warn!(
                                        "Error encountered mapping read: {}",
                                        map_res_opt.unwrap_err()
                                    );
                                    continue;
                                };

                                let query_name = String::from_utf8_lossy(name).into_owned();
                                // build GenomicAlignments and a parallel vector of
                                // their alignment scores (used by the
                                // score/combined probability sources).
                                galns.clear();
                                src_scores.clear();
                                for m in mappings.iter() {
                                    if let Some(ga) =
                                        mapping_to_genomic_alignment(m, &query_name, seq.len())
                                    {
                                        galns.push(ga);
                                        src_scores.push(crate::util::mapper::alignment_score(m));
                                    }
                                }
                                if galns.is_empty() {
                                    continue;
                                }
                                n_genome_mapped.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                                // attach the read sequence to the first alignment when
                                // soft-clip rescue is enabled (bramble shares one seq
                                // across the read group). bramble expects the seq in
                                // forward-reference orientation (as a BAM stores SEQ),
                                // so reverse-complement it when the canonical mapping is
                                // on the reverse strand — the read is held here in its
                                // original FASTQ orientation.
                                if use_fasta {
                                    galns[0].sequence = Some(if galns[0].is_reverse {
                                        crate::util::projection::revcomp(seq)
                                    } else {
                                        seq.to_vec()
                                    });
                                }

                                let projected =
                                    project_group_with(&galns, g2t, proj_config, &mut pctx);
                                if projected.is_empty() {
                                    if track_unprojectable {
                                        loc_unprojectable
                                            .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                                        let overlaps =
                                            crate::util::projection::overlapping_transcripts(
                                                &galns, g2t,
                                            );
                                        if overlaps.is_empty() {
                                            loc_unprojectable_intergenic
                                                .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                                        }
                                        for tid in overlaps {
                                            if let Some(slot) =
                                                loc_unprojectable_per_txp.get(tid as usize)
                                            {
                                                slot.fetch_add(
                                                    1,
                                                    std::sync::atomic::Ordering::Relaxed,
                                                );
                                            }
                                        }
                                    }
                                    continue;
                                }
                                n_projected.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                                let recs = projected_to_records(&projected, &src_scores);
                                // A read spanning no internal boundary says
                                // nothing about splice structure; only count
                                // reads that carry junction evidence.
                                if recs.iter().any(|r| r.junc_hits + r.junc_misses > 0) {
                                    loc_junc_informative
                                        .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                                    if recs.iter().all(|r| r.junc_misses > 0) {
                                        loc_junc_all_mismatch
                                            .fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                                    }
                                }
                                let (ag, aprobs) = filter.filter_projected(
                                    &mut discard_table,
                                    my_txp_info_view,
                                    &recs,
                                    seq.len(),
                                    beta,
                                    prob_source,
                                );

                                if !ag.is_empty() {
                                    // Junction evidence must ride along with the
                                    // group: the store is filled on the consumer
                                    // side, which never sees the projections.
                                    aln_group_junc.push((
                                        recs.iter().map(|r| r.junc_misses).min().unwrap_or(0),
                                        recs.iter().map(|r| r.junc_hits).max().unwrap_or(0),
                                    ));
                                    aln_group_alns.extend_from_slice(&ag);
                                    aln_group_probs.extend_from_slice(&aprobs);
                                    aln_group_boundaries.push(aln_group_alns.len());
                                    if let Some(ref mut names_vec) = aln_group_read_names {
                                        names_vec.push(query_name);
                                    }
                                    chunk_size += 1;
                                }

                                if chunk_size >= ALN_GROUP_CHUNK_LIMIT {
                                    aln_group_sender
                                        .send((
                                            aln_group_alns.clone(),
                                            aln_group_probs.clone(),
                                            aln_group_boundaries.clone(),
                                            aln_group_read_names,
                                            aln_group_junc.clone(),
                                        ))
                                        .expect("Error sending alignment group");
                                    aln_group_alns.clear();
                                    aln_group_probs.clear();
                                    aln_group_junc.clear();
                                    aln_group_boundaries.clear();
                                    aln_group_boundaries.push(0);
                                    aln_group_read_names = write_assignment_probs.then(Vec::new);
                                    chunk_size = 0;
                                }
                            }
                        }
                        if chunk_size > 0 {
                            aln_group_sender
                                .send((
                                    aln_group_alns,
                                    aln_group_probs,
                                    aln_group_boundaries,
                                    aln_group_read_names,
                                    aln_group_junc,
                                ))
                                .expect("Error sending alignment group");
                        }
                        discard_table
                    })
                })
                .collect();

            #[allow(clippy::useless_asref)]
            let txps_mut = txps.as_mut();
            let filter_opts_store = filter_opts.clone();
            let aln_group_consumer = s.spawn(move || {
                let mut name_vec = if filter_opts_store.write_assignment_probs {
                    Some(SwapVec::<String>::with_config(SwapVecConfig {
                        swap_after: Default::default(),
                        batch_size: Default::default(),
                        compression: Some(swapvec::Compression::Lz4),
                    }))
                } else {
                    None
                };

                let mut store = InMemoryAlignmentStore::new(filter_opts_store, txp_header);

                let pb = if args.quiet {
                    indicatif::ProgressBar::hidden()
                } else {
                    indicatif::ProgressBar::new_spinner().with_message("Number of reads mapped")
                };
                pb.set_style(
                    indicatif::ProgressStyle::with_template(
                        "[{elapsed_precise}] {spinner:4.green/blue} {msg} {human_pos:>12}",
                    )
                    .unwrap()
                    .tick_chars("⠁⠁⠉⠙⠚⠒⠂⠂⠒⠲⠴⠤⠄⠄⠤⠠⠠⠤⠦⠖⠒⠐⠐⠒⠓⠋⠉⠈⠈"),
                );
                pb.set_draw_target(indicatif::ProgressDrawTarget::stderr_with_hz(4));

                for (ags, aprobs, aln_boundaries, read_names, junc) in aln_group_receiver {
                    let mut junc_iter = junc.into_iter();
                    let mut reversed_read_names = if let Some(mut names_vec) = read_names {
                        names_vec.reverse();
                        Some(names_vec)
                    } else {
                        None
                    };

                    for window in aln_boundaries.windows(2) {
                        pb.inc(1);
                        let ag = &ags[window[0]..window[1]];
                        let as_probs = &aprobs[window[0]..window[1]];
                        let read_name_opt = if let Some(ref mut names_vec) = reversed_read_names {
                            names_vec.pop()
                        } else {
                            None
                        };

                        let (min_misses, max_hits) = junc_iter.next().unwrap_or((0, 0));
                        if store.add_filtered_group(ag, as_probs, txps_mut) {
                            store.min_junc_misses.push(min_misses);
                            store.max_junc_hits.push(max_hits);
                            if let Some(ref mut nvec) = name_vec {
                                let read_name =
                                    read_name_opt.unwrap_or(EMPTY_READ_NAME.to_string());
                                nvec.push(read_name)
                                    .expect("cannot push name to read name vector");
                            }
                            if ag.len() == 1 {
                                store.inc_unique_alignments();
                            }
                        }
                    }
                }
                pb.finish_with_message("Finished aligning and projecting reads.");
                (store, name_vec)
            });

            let aln_start = std::time::SystemTime::now();
            let total_reads = producer.join().expect("Producer thread panicked");

            let mut discard_tables: Vec<DiscardTable> = Vec::with_capacity(map_threads);
            for consumer in consumers {
                discard_tables.push(consumer.join().expect("Consumer thread panicked"));
            }

            drop(aln_group_sender);

            let (mut store, name_vec) = aln_group_consumer
                .join()
                .expect("Alignment group consumer panicked");

            info!(
                "Parsed {} total reads",
                total_reads.to_formatted_string(&Locale::en)
            );
            let n_gm = n_genome_mapped.load(std::sync::atomic::Ordering::Relaxed);
            let n_pr = n_projected.load(std::sync::atomic::Ordering::Relaxed);
            let n_ji = n_junc_informative.load(std::sync::atomic::Ordering::Relaxed);
            let n_jm = n_junc_all_mismatch.load(std::sync::atomic::Ordering::Relaxed);
            info!(
                "junction evidence: reads spanning >=1 internal boundary: {} ({:.2}% of projected);                  of those, disagreeing with EVERY candidate: {} ({:.2}%) -- candidate unannotated-isoform mass",
                n_ji.to_formatted_string(&Locale::en),
                if n_pr > 0 { 100.0 * n_ji as f64 / n_pr as f64 } else { 0.0 },
                n_jm.to_formatted_string(&Locale::en),
                if n_ji > 0 { 100.0 * n_jm as f64 / n_ji as f64 } else { 0.0 },
            );
            if track_unprojectable {
                store.unprojectable_reads =
                    n_unprojectable.load(std::sync::atomic::Ordering::Relaxed);
                store.unprojectable_intergenic =
                    n_unprojectable_intergenic.load(std::sync::atomic::Ordering::Relaxed);
                store.unprojectable_per_txp = unprojectable_per_txp
                    .iter()
                    .map(|c| c.load(std::sync::atomic::Ordering::Relaxed))
                    .collect();
                let attributed = store.unprojectable_reads - store.unprojectable_intergenic;
                info!(
                    "unprojectable reads: {} ({} attributable to an annotated locus, {} intergenic) -- reported, not quantified",
                    store.unprojectable_reads.to_formatted_string(&Locale::en),
                    attributed.to_formatted_string(&Locale::en),
                    store
                        .unprojectable_intergenic
                        .to_formatted_string(&Locale::en),
                );
            }
            info!(
                "reads aligned to genome: {}; of those, projected to >=1 transcript: {} ({:.2}%); projection-dropped: {}",
                n_gm.to_formatted_string(&Locale::en),
                n_pr.to_formatted_string(&Locale::en),
                if n_gm > 0 {
                    100.0 * (n_pr as f64) / (n_gm as f64)
                } else {
                    0.0
                },
                (n_gm - n_pr).to_formatted_string(&Locale::en),
            );
            let aln_time = aln_start.elapsed()?;
            info!(
                "Spliced genome alignment + projection of raw reads took: {}",
                humantime::format_duration(aln_time).to_string()
            );

            for dt in &discard_tables {
                store.aggregate_discard_table(dt);
            }
            Ok((store, name_vec, aln_time))
        },
    )?;

    perform_inference_and_write_output(
        txp_header,
        &mut store,
        name_vec,
        txps,
        txps_name,
        seqcol_digest,
        aln_time,
        args,
    )
}

fn get_source_type(pb: &std::path::Path) -> InputSourceType {
    let faq_endings = vec![
        ".fasta",
        ".fastq",
        ".FASTA",
        ".FASTQ",
        ".fa",
        ".fq",
        ".FA",
        ".FQ",
        ".fasta.gz",
        ".fastq.gz",
        ".FASTA.GZ",
        ".FASTQ.GZ",
        ".fa.gz",
        ".fq.gz",
        ".FA.GZ",
        ".FQ.GZ",
    ];
    let ubam_endings = vec![".bam", ".BAM", ".ubam", ".UBAM"];
    if let Some(ps) = pb.to_str() {
        for fe in faq_endings {
            if ps.ends_with(fe) {
                return InputSourceType::Fastx;
            }
        }
        for be in ubam_endings {
            if ps.ends_with(be) {
                return InputSourceType::Ubam;
            }
        }
    }

    InputSourceType::Unknown
}

#[allow(clippy::too_many_arguments)]
#[allow(unused_mut)]
pub fn quantify_bulk_alignments_raw_reads(
    header: &noodles_sam::Header,
    mut aligner: crate::util::mapper::Mapper,
    filter_opts: AlignmentFilters,
    read_paths: &[std::path::PathBuf],
    txps: &mut [TranscriptInfo],
    txps_name: &[String],
    args: &Args,
    seqcol_digest: NamedDigestVec,
) -> anyhow::Result<()> {
    // now parse the actual alignments for the reads and store the results
    // in our in-memory stor

    // we will take only shared refs to this
    let mut txp_info_view: Vec<TranscriptInfo> = Vec::with_capacity(txps.len());
    for ti in txps.iter() {
        txp_info_view.push(ti.clone());
    }

    // at least one mapping thread, otherwise everything but the fastx parser
    // and the in memory alignment store populator
    let map_threads = args.threads.saturating_sub(2).max(1);

    type ReadGroup = ReadChunkWithNames;
    type AlignmentGroupInfo = (
        Vec<AlnInfo>,
        Vec<f32>,
        Vec<usize>,
        Option<Vec<String>>,
        // No projection here, so no junction evidence.
        Vec<(i32, i32)>,
    );

    let (read_sender, read_receiver): (Sender<ReadGroup>, Receiver<ReadGroup>) =
        bounded(args.threads * 10);

    const READ_CHUNK_SIZE: usize = 200;
    let mut rpaths = vec![];
    read_paths.clone_into(&mut rpaths);

    // Producer thread: reads sequences and sends them to the channel
    let producer = std::thread::spawn(move || {
        let mut ctr = 0_usize;
        let mut chunk_size = 0_usize;
        let mut read_chunk = ReadChunkWithNames::new();

        // work shared between the two different
        // source types
        let mark_chunk = |chunk_size: &mut usize,
                          ctr: &mut usize,
                          read_chunk: &mut ReadGroup,
                          read_sender: &Sender<ReadGroup>| {
            *chunk_size += 1;
            *ctr += 1;
            if *chunk_size >= READ_CHUNK_SIZE {
                read_sender
                    .send(read_chunk.clone())
                    .expect("Error sending sequence");
                // prepare for the next chunk
                read_chunk.clear();
                *chunk_size = 0;
            }
        };

        // read from either a UBAM or (possibly compressed) FASTX file
        for read_path in rpaths {
            match get_source_type(&read_path) {
                InputSourceType::Ubam => {
                    let mut reader = std::fs::File::open(read_path)
                        .map(bam::io::Reader::new)
                        .expect("could not create BAM reader");
                    let header = reader.read_header().expect("could not read BAM header");
                    for result in reader.record_bufs(&header) {
                        let record = result.expect("Error reading ubam record");
                        record.add_to_read_group(&mut read_chunk);
                        mark_chunk(&mut chunk_size, &mut ctr, &mut read_chunk, &read_sender);
                    }
                }
                s @ (InputSourceType::Fastx | InputSourceType::Unknown) => {
                    if matches!(s, InputSourceType::Unknown) {
                        warn!(
                            "could not determine input file type for {} from suffix; assuming (possibly gzipped) fastx",
                            read_path.display()
                        );
                    }
                    let mut reader =
                        parse_fastx_file(read_path).expect("valid path/file to read sequences");
                    while let Some(result) = reader.next() {
                        let record = result.expect("Error reading record");
                        record.add_to_read_group(&mut read_chunk);
                        mark_chunk(&mut chunk_size, &mut ctr, &mut read_chunk, &read_sender);
                    }
                }
            }
        }
        // if any reads remain, send them off
        if chunk_size > 0 {
            read_sender
                .send(read_chunk)
                .expect("Error sending sequence");
        }
        ctr
    });

    // we need the scope here so we can borrow the relevant non-'static data
    let (mut store, name_vec, aln_time) = std::thread::scope(
        |s| -> anyhow::Result<(
            InMemoryAlignmentStore,
            Option<SwapVec<String>>,
            std::time::Duration,
        )> {
            const ALN_GROUP_CHUNK_LIMIT: usize = 100;

            let (aln_group_sender, aln_group_receiver): (
                Sender<AlignmentGroupInfo>,
                Receiver<AlignmentGroupInfo>,
            ) = bounded(args.threads * 100);

            // Consumer threads: receive sequences and perform alignment
            // Read names are also required by the coverage-signal export; keep this
            // in sync with the `AlignmentFilters` builder in `main.rs`, which enables
            // name retention for either output.
            let write_assignment_probs: bool =
                args.write_assignment_probs.is_some() || args.write_coverage_signals;
            let consumers: Vec<_> = (0..map_threads)
                .map(|_| {
                    let receiver = read_receiver.clone();
                    let mut filter = filter_opts.clone();
                    let loc_aligner = aligner.clone();

                    let my_txp_info_view = &txp_info_view;
                    let aln_group_sender = aln_group_sender.clone();
                    s.spawn(move || {
                        let mut discard_table = DiscardTable::new();

                        let mut chunk_size = 0_usize;
                        let mut aln_group_alns: Vec<AlnInfo> = Vec::new();
                        let mut aln_group_probs: Vec<f32> = Vec::new();
                        let mut aln_group_boundaries: Vec<usize> = Vec::new();
                        let mut aln_group_read_names = write_assignment_probs.then(Vec::new);
                        let mut aln_group_junc: Vec<(i32, i32)> = Vec::new();
                        let mut aln_group_junc: Vec<(i32, i32)> = Vec::new();
                        aln_group_boundaries.push(0);

                        // get the next chunk of reads
                        for read_chunk in receiver {
                            // iterate over every read
                            for (name, seq) in read_chunk.iter() {
                                // map the next read, with cigar string
                                let map_res_opt =
                                    crate::util::mapper::map_read(&loc_aligner, name, seq);
                                if let Ok(mut mappings) = map_res_opt {
                                    let (ag, aprobs, _polya) = filter.filter(
                                        &mut discard_table,
                                        header,
                                        my_txp_info_view,
                                        &mut mappings,
                                    );

                                    if !ag.is_empty() {
                                        aln_group_junc.push((0, 0));
                                        aln_group_alns.extend_from_slice(&ag);
                                        aln_group_probs.extend_from_slice(&aprobs);
                                        aln_group_boundaries.push(aln_group_alns.len());
                                        // if we are storing read names
                                        if let Some(ref mut names_vec) = aln_group_read_names {
                                            let name_str =
                                                String::from_utf8_lossy(name).into_owned();
                                            names_vec.push(name_str);
                                        }
                                        chunk_size += 1;
                                    }
                                    if chunk_size >= ALN_GROUP_CHUNK_LIMIT {
                                        aln_group_sender
                                            .send((
                                                aln_group_alns.clone(),
                                                aln_group_probs.clone(),
                                                aln_group_boundaries.clone(),
                                                aln_group_read_names,
                                                aln_group_junc.clone(),
                                            ))
                                            .expect("Error sending alignment group");
                                        aln_group_alns.clear();
                                        aln_group_probs.clear();
                                        aln_group_junc.clear();
                                        aln_group_boundaries.clear();
                                        aln_group_boundaries.push(0);
                                        aln_group_read_names =
                                            write_assignment_probs.then(Vec::new);
                                        chunk_size = 0;
                                    }
                                } else {
                                    warn!(
                                        "Error encountered mappread_ing read : {}",
                                        map_res_opt.unwrap_err()
                                    );
                                }
                            }
                        }
                        if chunk_size > 0 {
                            aln_group_sender
                                .send((
                                    aln_group_alns,
                                    aln_group_probs,
                                    aln_group_boundaries,
                                    aln_group_read_names,
                                    aln_group_junc,
                                ))
                                .expect("Error sending alignment group");
                        }
                        discard_table
                    })
                })
                .collect();

            #[allow(clippy::useless_asref)]
            let txps_mut = txps.as_mut();
            let filter_opts_store = filter_opts.clone();
            let aln_group_consumer = s.spawn(move || {
                let mut name_vec = if filter_opts_store.write_assignment_probs {
                    Some(SwapVec::<String>::with_config(SwapVecConfig {
                        swap_after: Default::default(),
                        batch_size: Default::default(),
                        compression: Some(swapvec::Compression::Lz4),
                    }))
                } else {
                    None
                };

                let mut store = InMemoryAlignmentStore::new(filter_opts_store, header);

                let pb = if args.quiet {
                    indicatif::ProgressBar::hidden()
                } else {
                    indicatif::ProgressBar::new_spinner().with_message("Number of reads mapped")
                };

                pb.set_style(
                    indicatif::ProgressStyle::with_template(
                        "[{elapsed_precise}] {spinner:4.green/blue} {msg} {human_pos:>12}",
                    )
                    .unwrap()
                    .tick_chars("⠁⠁⠉⠙⠚⠒⠂⠂⠒⠲⠴⠤⠄⠄⠤⠠⠠⠤⠦⠖⠒⠐⠐⠒⠓⠋⠉⠈⠈"),
                );
                pb.set_draw_target(indicatif::ProgressDrawTarget::stderr_with_hz(4));

                for (ags, aprobs, aln_boundaries, read_names, junc) in aln_group_receiver {
                    let mut junc_iter = junc.into_iter();
                    // if we are getting read names out then we are going to "reverse" them
                    // here so that we can simply pop the strings off the back to get them
                    // in order. We do this since we cannot otherwise "move" a string out of a
                    // Vec.
                    let mut reversed_read_names = if let Some(mut names_vec) = read_names {
                        names_vec.reverse();
                        Some(names_vec)
                    } else {
                        None
                    };

                    for window in aln_boundaries.windows(2) {
                        pb.inc(1);
                        let group_start = window[0];
                        let group_end = window[1];
                        let ag = &ags[group_start..group_end];
                        let as_probs = &aprobs[group_start..group_end];
                        let read_name_opt = if let Some(ref mut names_vec) = reversed_read_names {
                            names_vec.pop()
                        } else {
                            None
                        };

                        let (min_misses, max_hits) = junc_iter.next().unwrap_or((0, 0));
                        if store.add_filtered_group(ag, as_probs, txps_mut) {
                            store.min_junc_misses.push(min_misses);
                            store.max_junc_hits.push(max_hits);
                            if let Some(ref mut nvec) = name_vec {
                                let read_name =
                                    read_name_opt.unwrap_or(EMPTY_READ_NAME.to_string());
                                nvec.push(read_name)
                                    .expect("cannot push name to read name vector");
                            }
                            if ag.len() == 1 {
                                store.inc_unique_alignments();
                            }
                        }
                    }
                }
                pb.finish_with_message("Finished aligning reads.");
                (store, name_vec)
            });

            let aln_start = std::time::SystemTime::now();
            // Wait for the producer to finish reading
            let total_reads = producer.join().expect("Producer thread panicked");

            let mut discard_tables: Vec<DiscardTable> = Vec::with_capacity(map_threads);
            for consumer in consumers {
                let dt = consumer.join().expect("Consumer thread panicked");
                discard_tables.push(dt);
            }

            drop(aln_group_sender);

            let (mut store, name_vec) = aln_group_consumer
                .join()
                .expect("Alignment group consumer panicked");

            info!(
                "Parsed {} total reads",
                total_reads.to_formatted_string(&Locale::en)
            );

            let aln_time = aln_start.elapsed()?;
            info!(
                "Alignment of raw reads using rammap took: {}",
                humantime::format_duration(aln_time).to_string()
            );

            for dt in &discard_tables {
                store.aggregate_discard_table(dt);
            }
            Ok((store, name_vec, aln_time))
        },
    )?;

    perform_inference_and_write_output(
        header,
        &mut store,
        name_vec,
        txps,
        txps_name,
        seqcol_digest,
        aln_time,
        args,
    )
}

#[cfg(test)]
mod ambiguity_component_tests {
    use super::ambiguity_components;
    use crate::util::oarfish_types::AlnInfo;
    use bio_types::strand::Strand;

    fn group(ref_ids: &[u32]) -> Vec<AlnInfo> {
        ref_ids
            .iter()
            .map(|&ref_id| AlnInfo {
                ref_id,
                start: 1,
                end: 100,
                strand: Strand::Forward,
                left_clip: 0,
                right_clip: 0,
            })
            .collect()
    }

    fn components_of(groups: &[Vec<AlnInfo>], n: usize) -> Vec<usize> {
        ambiguity_components(groups.iter().map(|g| g.as_slice()), n)
    }

    #[test]
    fn unlinked_transcripts_are_singletons() {
        let c = components_of(&[], 4);
        assert_eq!(c, vec![0, 1, 2, 3]);
    }

    #[test]
    fn unique_reads_do_not_link_anything() {
        let groups = vec![group(&[0]), group(&[1])];
        let c = components_of(&groups, 3);
        assert_eq!(c[0], 0);
        assert_ne!(c[0], c[1]);
    }

    #[test]
    fn ambiguous_reads_merge_their_candidates() {
        let groups = vec![group(&[0, 2])];
        let c = components_of(&groups, 4);
        assert_eq!(c[0], c[2]);
        assert_ne!(c[0], c[1]);
        assert_ne!(c[0], c[3]);
    }

    #[test]
    fn components_merge_transitively_through_shared_reads() {
        // 0-1 and 1-2 are linked by different reads, so all three are one
        // component even though no read lists 0 and 2 together.
        let groups = vec![group(&[0, 1]), group(&[1, 2]), group(&[4, 5])];
        let c = components_of(&groups, 6);
        assert_eq!(c[0], c[1]);
        assert_eq!(c[1], c[2]);
        assert_eq!(c[4], c[5]);
        assert_ne!(c[0], c[4]);
        assert_ne!(c[0], c[3]);
        assert_eq!(c.iter().collect::<std::collections::HashSet<_>>().len(), 3);
    }
}
