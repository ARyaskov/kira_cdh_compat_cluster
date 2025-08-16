# kira_cdh_compat_cluster

Greedy clustering for CD-HIT–like pipelines in **Rust (edition 2024)**.

- **Input**: an LSH index and fixed-length MinHash/KMV signatures (one per sequence), plus optional sequence lengths.
- **Similarity**: MinHash-based Jaccard estimate (fraction of equal positions in two signatures).
- **Coverage gating**: optional CD-HIT–style `-aS/-aL` checks derived from Jaccard and sequence lengths (no alignment needed).
- **Deterministic**: stable order; optional length-first processing to mimic CD-HIT.
- **Parallel**: optional multi-threaded candidate checks via the `parallel` feature.

This crate is designed to “glue” a candidate retrieval stage (e.g., LSH) with a lightweight, deterministic 
greedy clustering engine that emits cluster partitions compatible with CD-HIT expectations.

---

## Installation

Add to your Cargo workspace (recommended):

```toml
[dependencies]
kira_cdh_compat_cluster = "*"

# If you need multi-threaded candidate checks:
[features]
default = []
# enable at build time: --features parallel
````

The engine expects LSH and MinHash/KMV signatures from relative crates, e.g.:

* `kira_cdh_compat_fastq_reader` — zero-copy FASTQ/FASTA reading (sync/async), with gzip support.
* `kira_cdh_compat_kmer_indexer` — k-mer indexer.
* `kira_cdh_compat_lsh` — KMV sketches and LSH index (insert/build/query).

---

## Quickstart

```rust
use kira_cdh_compat_cluster::{GreedyClusterer, ClusterOptions};
use kira_cdh_compat_lsh::{kmv::KmvSketch, lsh::{LshIndex, LshParams}};
use kira_cdh_compat_fastq_reader::FastqReader;

// 1) Build MinHash/KMV signatures and lengths
let mut signatures: Vec<Vec<u64>> = Vec::new(); // all same length (e.g., 128)
let mut lengths:    Vec<u32>      = Vec::new();

let k: usize = 11; // k-mer size used for hashing (tune per dataset)
let mut r = FastqReader::from_path("reads.fastq.gz", Default::default())?;
for item in &mut r {
    let rec = item?;
    let seq = rec.seq().as_bytes();
    lengths.push(seq.len() as u32);

    let mut kmv = KmvSketch::new(128); // signature length
    for win in seq.windows(k) {
        // Implement a fast, stable u64 hash for k-mers; pick a production hash (e.g., AHash, XXH3) or a 2-bit encoder + mix.
        let h = my_hash_u64(win);
        kmv.update(h);
    }
    signatures.push(kmv.finish());
}

// 2) Build an LSH index for candidate retrieval
let params = LshParams::new(32, 4)?; // bands, rows — tune for your target similarity
let mut index = LshIndex::with_params(params);
for (i, sig) in signatures.iter().enumerate() {
    index.insert(i as u64, sig)?;
}
index.build();

// 3) Run greedy clustering
let opts = ClusterOptions {
    similarity_threshold: 0.97,   // MinHash Jaccard gate
    max_candidates: 256,          // per-seed candidate cap from LSH
    sort_by_length: true,         // mimic CD-HIT (longer seeds first)
    parallel_sim_checks: true,    // requires building with --features parallel
    min_coverage_short: None,     // CD-HIT-like -aS (short-side coverage)
    min_coverage_long:  None,     // CD-HIT-like -aL (long-side coverage)
    kmer_k: None,                 // set when you enable coverage gates
};
let clusterer = GreedyClusterer {
    index: &index,
    signatures: &signatures,
    seq_lengths: Some(&lengths),
    options: opts,
};
let result = clusterer.run();

println!("clusters: {}", result.clusters.len());
// result.clusters is Vec<Vec<usize>> (sequence indices per cluster)
```

> To emit CD-HIT-compatible `.clstr`, pair this crate with a writer (e.g., a small `kira_cdh_compat_clstr` utility crate) and print:
>
> ```
> >Cluster 0
> 0\t150nt, >seqA... *
> 1\t140nt, >seqB...
> ...
> ```

---

## Coverage gating (CD-HIT-style `-aS/-aL`)

When sequences differ in length, a high MinHash Jaccard can still arise from partial overlap.
CD-HIT mitigates this with short/long coverage gates.
We provide a **fast approximation** that requires only Jaccard, lengths, and `k`.

* Let `L_s` and `L_l` be short/long sequence lengths; `k` is the k-mer size.
* Window counts:

  ```
  A = max(L_s - k + 1, 0)
  B = max(L_l - k + 1, 0)
  ```
* With `Ĵ` the MinHash Jaccard estimate (equal positions / signature length), approximate the number of shared windows:

  ```
  I ≈ Ĵ * (A + B) / (1 + Ĵ)
  ```
* Coverage for short/long:

  ```
  cov_short = I / A
  cov_long  = I / B
  ```
* Accept a candidate if:

  ```
  Ĵ ≥ similarity_threshold
  && (cov_short ≥ min_coverage_short  if set)
  && (cov_long  ≥ min_coverage_long   if set)
  ```

**Enabling coverage gates:**

```rust
let opts = ClusterOptions {
    similarity_threshold: 0.95,
    max_candidates: 256,
    sort_by_length: true,
    parallel_sim_checks: true,
    min_coverage_short: Some(0.90), // -aS
    min_coverage_long:  Some(0.40), // -aL (optional)
    kmer_k: Some(11),               // required when gating by coverage
};
```

### Practical profiles (nucleotide)

| Task                      | k     | `similarity_threshold` | `min_coverage_short` | `min_coverage_long` | Notes                                |
| ------------------------- | ----- | ---------------------- | -------------------- | ------------------- | ------------------------------------ |
| High-identity clustering  | 9–11  | 0.95–0.99              | 0.90–0.98            | 0.0–0.5             | Stable contigs/transcripts           |
| 97% OTU/ASV-like grouping | 9–12  | 0.90–0.97              | 0.85–0.95            | 0.0–0.5             | 16S/SSU; set `-aS` slightly below Ĵ |
| Aggressive deduplication  | 11–13 | 0.98–0.995             | 0.95–0.99            | 0.0–0.5             | Long reads or assembled contigs      |

**Hints**:

* Larger `k` tightens windows—great for high identity, too strict for noisy data.
* Use length-first ordering (`sort_by_length=true`) to reduce false merges.
* Tune LSH `(bands, rows)` to bend the collision curve near your target Jaccard.

---

## API overview

```rust
pub struct ClusterOptions {
    pub similarity_threshold: f32,
    pub max_candidates: usize,
    pub sort_by_length: bool,
    pub parallel_sim_checks: bool,
    pub min_coverage_short: Option<f32>, // like -aS
    pub min_coverage_long:  Option<f32>, // like -aL
    pub kmer_k: Option<u32>,             // required when coverage gates are used
}

pub struct GreedyClusterer<'a> {
    pub index: &'a LshIndex,
    pub signatures: &'a [Vec<u64>],  // fixed-length MinHash/KMV signatures
    pub seq_lengths: Option<&'a [u32]>,
    pub options: ClusterOptions,
}

impl<'a> GreedyClusterer<'a> {
    pub fn run(&self) -> ClusterResult;
}

pub struct ClusterResult {
    pub clusters: Vec<Vec<usize>>,  // sequence indices per cluster
}

/// MinHash Jaccard estimate (equal positions / signature length).
pub fn jaccard_minhash(a: &[u64], b: &[u64]) -> f32;
```

* The engine assumes all signatures share the same length.
* `LshIndex` should return candidate sequence IDs (not internal bucket IDs).
* Deterministic behavior: the processing order is fixed; with `sort_by_length=true`, longer sequences seed first.

---

## Performance tips

* Enable `--features parallel` to parallelize similarity checks within candidate sets.
* Batch insert/build in LSH; use “insert\_batch” APIs when available to reduce allocation churn.
* Keep signature lengths moderate (e.g., 128–256 u64s); increase only if collision rates are too high.
* Use a fast, stable k-mer hash (SIMD-friendly, no allocations).

---

## Examples

* `examples/greedy_basic.rs` — minimal run with synthetic signatures.
* `examples/greedy_coverage.rs` — same with coverage gating enabled.

Run:

```bash
cargo run --release --example greedy_basic
cargo run --release --example greedy_coverage --features parallel
```

---

## License

GPLv2

