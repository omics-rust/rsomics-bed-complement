use criterion::{Criterion, criterion_group, criterion_main};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::process::Command;

const N_RECORDS: usize = 50_000;
const CHROM_SIZE: u64 = 100_000_000;
const SEED: u64 = 0x00C0_FFEE;

fn xorshift(x: &mut u64) -> u64 {
    *x ^= *x << 13;
    *x ^= *x >> 7;
    *x ^= *x << 17;
    *x
}

fn synth_fixtures(bed: &PathBuf, genome: &PathBuf) {
    let chroms = [
        ("chr1", CHROM_SIZE),
        ("chr2", CHROM_SIZE / 2),
        ("chr3", CHROM_SIZE / 3),
    ];

    // Genome file.
    {
        let f = File::create(genome).expect("create genome");
        let mut w = BufWriter::new(f);
        for (c, s) in &chroms {
            writeln!(w, "{c}\t{s}").unwrap();
        }
    }

    // Sorted BED file.
    {
        let mut rows: Vec<(String, u64, u64)> = Vec::with_capacity(N_RECORDS);
        let mut rng = SEED;
        for _ in 0..N_RECORDS {
            let (chrom, size) = chroms[(xorshift(&mut rng) % chroms.len() as u64) as usize];
            let start = xorshift(&mut rng) % (size - 1000);
            let end = start + 100 + (xorshift(&mut rng) % 900);
            rows.push((chrom.to_string(), start, end));
        }
        rows.sort_unstable_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1)));
        let f = File::create(bed).expect("create bed");
        let mut w = BufWriter::new(f);
        for (c, s, e) in rows {
            writeln!(w, "{c}\t{s}\t{e}").unwrap();
        }
    }
}

fn ensure_fixtures() -> (PathBuf, PathBuf) {
    let mut bed = std::env::temp_dir();
    bed.push(format!("rsomics-bed-complement-bench-{N_RECORDS}.bed"));
    let mut genome = std::env::temp_dir();
    genome.push("rsomics-bed-complement-bench-genome.txt");
    if !bed.exists() || !genome.exists() {
        synth_fixtures(&bed, &genome);
    }
    (bed, genome)
}

fn bench(c: &mut Criterion) {
    let (bed, genome) = ensure_fixtures();
    let ours = env!("CARGO_BIN_EXE_rsomics-bed-complement");
    let mut group = c.benchmark_group(format!("bed_complement/{N_RECORDS}"));
    group.sample_size(10);

    group.bench_function("rsomics-bed-complement", |bm| {
        bm.iter(|| {
            let out = Command::new(ours)
                .arg(&bed)
                .arg("-g")
                .arg(&genome)
                .output()
                .expect("ours run");
            assert!(
                out.status.success(),
                "{}",
                String::from_utf8_lossy(&out.stderr)
            );
        });
    });

    if Command::new("bedtools")
        .arg("--version")
        .output()
        .is_ok_and(|o| o.status.success())
    {
        group.bench_function("bedtools-complement", |bm| {
            bm.iter(|| {
                let out = Command::new("bedtools")
                    .args(["complement", "-i"])
                    .arg(&bed)
                    .arg("-g")
                    .arg(&genome)
                    .output()
                    .expect("bedtools run");
                assert!(
                    out.status.success(),
                    "{}",
                    String::from_utf8_lossy(&out.stderr)
                );
            });
        });
    } else {
        eprintln!("bedtools not on PATH — skipping upstream comparison");
    }

    group.finish();
}

criterion_group!(benches, bench);
criterion_main!(benches);
