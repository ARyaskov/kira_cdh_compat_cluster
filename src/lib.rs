//! Greedy clustering engine for CD-HIT-like workflows.
//!
//! - Input: LSH index, fixed-length MinHash/KMV signatures, optional sequence lengths.
//! - Similarity: MinHash-based Jaccard estimate.
//! - Optional coverage gating (CD-HIT-like -aS/-aL) derived from Jaccard + lengths + k.
//!
//! This crate aims to keep the API small, deterministic, and well-documented.

pub mod greedy;

pub use greedy::{ClusterOptions, ClusterResult, GreedyClusterer, jaccard_minhash};
