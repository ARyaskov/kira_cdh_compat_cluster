use kira_cdh_compat_lsh::lsh::LshIndex;

/// Options controlling greedy clustering behavior.
#[derive(Clone, Debug)]
pub struct ClusterOptions {
    /// Minimum MinHash-based Jaccard estimate to accept a candidate.
    pub similarity_threshold: f32, // e.g., 0.97

    /// Maximum number of LSH candidates per seed.
    pub max_candidates: usize, // e.g., 256

    /// If true, sequences are processed in descending length order.
    /// Requires `seq_lengths` to be provided to the clusterer.
    pub sort_by_length: bool,

    /// If enabled (feature `parallel`), evaluate candidate similarities in parallel.
    pub parallel_sim_checks: bool,

    /// Optional coverage filter (CD-HIT-like).
    /// If set, both fields are interpreted in [0.0, 1.0].
    /// - Short-side coverage threshold (analog of -aS).
    pub min_coverage_short: Option<f32>,
    /// - Long-side coverage threshold (analog of -aL).
    pub min_coverage_long: Option<f32>,

    /// k-mer size used to approximate window counts (L - k + 1).
    /// Required if coverage thresholds are set.
    pub kmer_k: Option<u32>,
}

/// Result of greedy clustering: list of clusters of sequence indices.
#[derive(Clone, Debug)]
pub struct ClusterResult {
    pub clusters: Vec<Vec<usize>>,
}

/// Greedy clustering engine.
///
/// Assumptions:
/// - All `signatures` have identical length.
/// - `index` can retrieve candidate IDs for a given signature.
/// - If `sort_by_length` or coverage thresholds are enabled,
///   `seq_lengths` must be provided (base pairs / nt length).
pub struct GreedyClusterer<'a> {
    pub index: &'a LshIndex,
    pub signatures: &'a [Vec<u64>],
    pub seq_lengths: Option<&'a [u32]>,
    pub options: ClusterOptions,
}

impl<'a> GreedyClusterer<'a> {
    /// Run greedy clustering and return clusters of sequence indices.
    ///
    /// Algorithm:
    /// 1) Determine processing order (optionally by length).
    /// 2) For each unassigned seed `i`:
    ///    - Start a new cluster with `i`.
    ///    - Query LSH candidates.
    ///    - For each candidate `j` not yet assigned:
    ///      * Check MinHash-Jaccard threshold.
    ///      * If coverage thresholds are set, check coverage gating (short/long).
    ///      * If both pass, assign `j` to the cluster.
    pub fn run(&self) -> ClusterResult {
        self.validate_inputs();

        let n = self.signatures.len();
        if n == 0 {
            return ClusterResult {
                clusters: Vec::new(),
            };
        }

        let order = self.order_indices();
        let mut clusters: Vec<Vec<usize>> = Vec::new();
        let mut assigned = vec![false; n];

        for &i in &order {
            if assigned[i] {
                continue;
            }

            let mut cluster = Vec::with_capacity(64);
            cluster.push(i);
            assigned[i] = true;

            let cand_ids = self
                .index
                .query_candidates(&self.signatures[i], self.options.max_candidates);

            // Keep only not-yet-assigned and not self
            let cands: Vec<usize> = cand_ids
                .into_iter()
                .map(|(_, seq_id)| seq_id as usize)
                .filter(|&j| j != i && !assigned[j])
                .collect();

            let sig_i = &self.signatures[i];
            let thr = self.options.similarity_threshold;

            if self.options.parallel_sim_checks {
                #[cfg(feature = "parallel")]
                {
                    use rayon::prelude::*;
                    let passed: Vec<usize> = cands
                        .par_iter()
                        .copied()
                        .filter(|&j| {
                            let sim = jaccard_minhash(sig_i, &self.signatures[j]);
                            if sim < thr {
                                return false;
                            }
                            if self.coverage_filters_enabled() {
                                self.passes_coverage(i, j, sim)
                            } else {
                                true
                            }
                        })
                        .collect();

                    for j in passed {
                        if !assigned[j] {
                            assigned[j] = true;
                            cluster.push(j);
                        }
                    }
                }

                #[cfg(not(feature = "parallel"))]
                {
                    // Fallback to sequential when feature is disabled.
                    let mut to_add = Vec::new();
                    for j in cands {
                        if assigned[j] {
                            continue;
                        }
                        let sim = jaccard_minhash(sig_i, &self.signatures[j]);
                        if sim >= thr
                            && (!self.coverage_filters_enabled() || self.passes_coverage(i, j, sim))
                        {
                            to_add.push(j);
                        }
                    }
                    for j in to_add {
                        if !assigned[j] {
                            assigned[j] = true;
                            cluster.push(j);
                        }
                    }
                }
            } else {
                // Sequential path
                for j in cands {
                    if assigned[j] {
                        continue;
                    }
                    let sim = jaccard_minhash(sig_i, &self.signatures[j]);
                    if sim < thr {
                        continue;
                    }
                    if self.coverage_filters_enabled() && !self.passes_coverage(i, j, sim) {
                        continue;
                    }
                    assigned[j] = true;
                    cluster.push(j);
                }
            }

            clusters.push(cluster);
        }

        ClusterResult { clusters }
    }

    #[inline]
    fn coverage_filters_enabled(&self) -> bool {
        self.options.min_coverage_short.is_some() || self.options.min_coverage_long.is_some()
    }

    /// Check coverage constraints (CD-HIT-like -aS/-aL) using a MinHash-based approximation.
    ///
    /// Given MinHash Jaccard estimate `sim` and sequence lengths `L_i`, `L_j`,
    /// estimate the number of shared k-mer windows:
    ///   A = max(L_short - k + 1, 0)
    ///   B = max(L_long  - k + 1, 0)
    ///   I ≈ sim * (A + B) / (1 + sim)
    /// Then coverage for short = I / A, and for long = I / B.
    fn passes_coverage(&self, i: usize, j: usize, sim: f32) -> bool {
        let k = self
            .options
            .kmer_k
            .expect("Coverage thresholds require `kmer_k`");
        let lens = self
            .seq_lengths
            .expect("Coverage thresholds require `seq_lengths`");

        let (li, lj) = (lens[i] as i64, lens[j] as i64);
        let (short_idx, long_idx, ls, ll) = if li <= lj {
            (i, j, li, lj)
        } else {
            (j, i, lj, li)
        };

        // Window counts (non-negative).
        let k_i = k as i64;
        let a = (ls - k_i + 1).max(0) as f64; // short windows
        let b = (ll - k_i + 1).max(0) as f64; // long  windows

        if a == 0.0 || b == 0.0 {
            // degenerate case: treat as failing coverage
            return false;
        }

        let s = sim as f64;
        // I ≈ s * (A + B) / (1 + s)
        let inter = s * (a + b) / (1.0 + s);

        let cov_short = inter / a;
        let cov_long = inter / b;

        if let Some(min_s) = self.options.min_coverage_short {
            if cov_short < min_s as f64 {
                return false;
            }
        }
        if let Some(min_l) = self.options.min_coverage_long {
            if cov_long < min_l as f64 {
                return false;
            }
        }

        // Silence unused warnings for indices (kept for potential future logging)
        let _ = (short_idx, long_idx);
        true
    }

    fn validate_inputs(&self) {
        if self.signatures.is_empty() {
            return;
        }
        // Ensure equal signature length.
        let m = self.signatures[0].len();
        debug_assert!(
            self.signatures.iter().all(|s| s.len() == m),
            "All signatures must be of equal length"
        );

        // Sorting requires lengths.
        if self.options.sort_by_length {
            let lens = self
                .seq_lengths
                .expect("sort_by_length=true requires `seq_lengths`");
            debug_assert_eq!(
                lens.len(),
                self.signatures.len(),
                "`seq_lengths` must match signatures length"
            );
        }

        // Coverage thresholds require both lengths and k.
        if self.coverage_filters_enabled() {
            debug_assert!(
                self.options.kmer_k.is_some(),
                "Coverage thresholds require `kmer_k`"
            );
            let lens = self
                .seq_lengths
                .expect("coverage thresholds require `seq_lengths`");
            debug_assert_eq!(
                lens.len(),
                self.signatures.len(),
                "`seq_lengths` must match signatures length"
            );
        }
    }

    fn order_indices(&self) -> Vec<usize> {
        let n = self.signatures.len();
        let mut idx: Vec<usize> = (0..n).collect();

        if self.options.sort_by_length {
            let lens = self.seq_lengths.expect("seq_lengths is required");
            idx.sort_unstable_by_key(|&i| std::cmp::Reverse(lens[i]));
        }

        idx
    }
}

/// Estimate Jaccard similarity using MinHash signatures.
/// Returns a value in [0.0, 1.0].
///
/// Assumes equal-length signatures; compares position-wise equality.
pub fn jaccard_minhash(a: &[u64], b: &[u64]) -> f32 {
    debug_assert_eq!(a.len(), b.len(), "signatures must be equal length");
    let mut same = 0usize;
    for i in 0..a.len() {
        if a[i] == b[i] {
            same += 1;
        }
    }
    same as f32 / (a.len() as f32)
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::{Rng, SeedableRng, rngs::StdRng};

    fn make_sig(base: u64, len: usize) -> Vec<u64> {
        (0..len).map(|i| base.wrapping_add(i as u64)).collect()
    }

    #[test]
    fn test_jaccard_minhash() {
        let a = vec![1, 2, 3, 4];
        let b = vec![1, 2, 9, 9];
        assert!((jaccard_minhash(&a, &b) - 0.5).abs() < 1e-6);
    }

    #[test]
    fn test_order_with_lengths() {
        let idx = LshIndex::with_params(kira_cdh_compat_lsh::lsh::LshParams::new(8, 4).unwrap());
        let sigs = vec![make_sig(0, 4), make_sig(10, 4), make_sig(20, 4)];
        let lens = vec![100, 300, 200];
        let gc = GreedyClusterer {
            index: &idx,
            signatures: &sigs,
            seq_lengths: Some(&lens),
            options: ClusterOptions {
                similarity_threshold: 0.5,
                max_candidates: 10,
                sort_by_length: true,
                parallel_sim_checks: false,
                min_coverage_short: None,
                min_coverage_long: None,
                kmer_k: None,
            },
        };
        let order = gc.order_indices();
        assert_eq!(order, vec![1, 2, 0]); // 300, 200, 100
    }

    #[test]
    fn test_coverage_gate_basic() {
        // Two signatures, 75% equal positions => J ~ 0.75.
        let a = vec![1, 2, 3, 4];
        let b = vec![1, 2, 3, 9];
        let sim = jaccard_minhash(&a, &b);
        assert!((sim - 0.75).abs() < 1e-6);

        let idx = LshIndex::with_params(kira_cdh_compat_lsh::lsh::LshParams::new(8, 4).unwrap());
        let sigs = vec![a, b];
        let lens = vec![100, 120];

        let gc = GreedyClusterer {
            index: &idx,
            signatures: &sigs,
            seq_lengths: Some(&lens),
            options: ClusterOptions {
                similarity_threshold: 0.0, // ignore sim gate for this test
                max_candidates: 0,
                sort_by_length: false,
                parallel_sim_checks: false,
                min_coverage_short: Some(0.5),
                min_coverage_long: Some(0.4),
                kmer_k: Some(10),
            },
        };

        // Use the internal function indirectly via passes_coverage through run()
        // Quick sanity: just ensure validation does not panic and coverage gate is computable.
        gc.validate_inputs();
        // Manual check:
        let ok = gc.passes_coverage(0, 1, 0.75);
        assert!(ok);
    }
}
