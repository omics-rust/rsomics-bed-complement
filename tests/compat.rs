use std::path::Path;
use std::process::Command;

use rsomics_bed_complement::{complement, read_genome};

fn golden(name: &str) -> std::path::PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn basic_complement_correctness() {
    let input = golden("input.bed");
    let genome_path = golden("genome.txt");
    let (genome_order, genome) = read_genome(&genome_path).unwrap();
    let mut out = Vec::new();
    complement(&input, &genome_order, &genome, &mut out).unwrap();
    let result = String::from_utf8(out).unwrap();
    // chr1 covered: [100,300) and [500,700)
    // complement: [0,100), [300,500), [700,1000)
    assert!(
        result.contains("chr1\t0\t100"),
        "chr1 leading gap: {result}"
    );
    assert!(
        result.contains("chr1\t300\t500"),
        "chr1 middle gap: {result}"
    );
    assert!(
        result.contains("chr1\t700\t1000"),
        "chr1 trailing gap: {result}"
    );
    // chr2 covered: [200,400)
    // complement: [0,200), [400,600)
    assert!(
        result.contains("chr2\t0\t200"),
        "chr2 leading gap: {result}"
    );
    assert!(
        result.contains("chr2\t400\t600"),
        "chr2 trailing gap: {result}"
    );
    // chr3 has no features → full interval [0,500)
    assert!(
        result.contains("chr3\t0\t500"),
        "chr3 full interval: {result}"
    );
    let lines: Vec<&str> = result.lines().filter(|l| !l.is_empty()).collect();
    assert_eq!(lines.len(), 6, "expected 6 complement intervals: {result}");
}

#[test]
fn empty_input_gives_full_chromosomes() {
    use std::io::Write;
    use tempfile::NamedTempFile;
    let mut f = NamedTempFile::new().unwrap();
    writeln!(f, "# empty").unwrap();
    let genome_path = golden("genome.txt");
    let (genome_order, genome) = read_genome(&genome_path).unwrap();
    let mut out = Vec::new();
    complement(f.path(), &genome_order, &genome, &mut out).unwrap();
    let result = String::from_utf8(out).unwrap();
    assert!(result.contains("chr1\t0\t1000"), "chr1 full: {result}");
    assert!(result.contains("chr2\t0\t600"), "chr2 full: {result}");
    assert!(result.contains("chr3\t0\t500"), "chr3 full: {result}");
}

#[test]
fn bedtools_compat() {
    let bedtools = Command::new("bedtools").arg("--version").output();
    if bedtools.is_err() || !bedtools.unwrap().status.success() {
        eprintln!("bedtools not available — skipping compat test");
        return;
    }

    let input = golden("input.bed");
    let genome_path = golden("genome.txt");
    let (genome_order, genome) = read_genome(&genome_path).unwrap();

    let mut ours = Vec::new();
    complement(&input, &genome_order, &genome, &mut ours).unwrap();
    let ours_str = String::from_utf8(ours).unwrap();

    let bt = Command::new("bedtools")
        .args(["complement", "-i"])
        .arg(&input)
        .arg("-g")
        .arg(&genome_path)
        .output()
        .expect("bedtools complement failed");
    let bt_str = String::from_utf8(bt.stdout).unwrap();

    let mut ours_lines: Vec<&str> = ours_str.lines().filter(|l| !l.is_empty()).collect();
    let mut bt_lines: Vec<&str> = bt_str.lines().filter(|l| !l.is_empty()).collect();
    ours_lines.sort_unstable();
    bt_lines.sort_unstable();

    assert_eq!(
        ours_lines, bt_lines,
        "output differs from bedtools complement"
    );
}
