use criterion::{Criterion, criterion_group, criterion_main};
use std::hint::black_box;
use std::path::PathBuf;
use std::process::Command;

fn bench_bed_complement(c: &mut Criterion) {
    let bin = env!("CARGO_BIN_EXE_rsomics-bed-complement");
    let manifest = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let bed = manifest.join("tests/golden/input.bed");
    let genome = manifest.join("tests/golden/genome.txt");
    c.bench_function("rsomics-bed-complement golden", |b| {
        b.iter(|| {
            let out = Command::new(black_box(bin))
                .args([bed.to_str().unwrap(), "-g", genome.to_str().unwrap()])
                .output()
                .unwrap();
            assert!(out.status.success());
        });
    });
}

criterion_group!(benches, bench_bed_complement);
criterion_main!(benches);
