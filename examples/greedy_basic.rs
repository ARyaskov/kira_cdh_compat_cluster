use kira_cdh_compat_cluster::{ClusterOptions, GreedyClusterer};
use kira_cdh_compat_lsh::{
    kmv::KmvSketch,
    lsh::{LshIndex, LshParams},
};

fn main() -> anyhow::Result<()> {
    // Synthetic signatures for demo
    let make_sig = |seed: u64| {
        let mut kmv = KmvSketch::new(64);
        for i in 0..256u64 {
            kmv.update(seed.wrapping_mul(1_000_003).wrapping_add(i));
        }
        kmv.finish()
    };

    let sigs = vec![make_sig(1), make_sig(2), make_sig(1), make_sig(42)];
    let lengths = vec![1500u32, 900, 1400, 800];

    let params = LshParams::new(16, 4)?; // bands=16, rows=4
    let mut index = LshIndex::with_params(params);
    for (i, s) in sigs.iter().enumerate() {
        index.insert(i as u32, s)?;
    }
    index.build();

    let opts = ClusterOptions {
        similarity_threshold: 0.8,
        max_candidates: 64,
        sort_by_length: true,
        parallel_sim_checks: true, // enable parallel check if built with `--features parallel`
        min_coverage_short: Some(0.9), // analogous to -aS
        min_coverage_long: Some(0.0), // disable long-side gate
        kmer_k: Some(11),          // k used to approximate window counts
    };

    let clusterer = GreedyClusterer {
        index: &index,
        signatures: &sigs,
        seq_lengths: Some(&lengths),
        options: opts,
    };

    let res = clusterer.run();
    for (cid, members) in res.clusters.iter().enumerate() {
        println!(">Cluster {cid}");
        for &m in members {
            println!("  - seq#{m}");
        }
    }
    Ok(())
}
