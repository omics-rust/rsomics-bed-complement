//! BED complement — intervals not covered by any input feature.
//!
//! Equivalent to `bedtools complement -i input -g genome`. Requires a sorted
//! input (sorted by chrom, then start — same requirement as `bedtools
//! complement`). For each chromosome in the genome file, emits the uncovered
//! gaps: [0, first_start), [prev_end, next_start), …, [last_end, chrom_size).
//!
//! Chromosomes with no features in the input are emitted as a single interval
//! spanning the whole chromosome.
//!
//! Chromosomes in the genome file that never appear in the input are emitted
//! as a full-length interval.
//!
//! Algorithm: two-pass O(N):
//!   Pass 1 — collect all intervals per chrom into Vec<(start,end)>.
//!   Pass 2 — for each chrom in genome-file order, emit gaps.
//! No streaming (must collect to handle multi-record chroms correctly), but
//! typical BED files are in-RAM small compared to the analysis pipeline.

use std::collections::HashMap;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

/// Parse a two-column genome/chromsizes file into an ordered list and a map.
///
/// Returns `(ordered_chroms, chrom→size)` preserving genome-file order.
pub fn read_genome(path: &Path) -> Result<(Vec<String>, HashMap<String, u64>)> {
    let data = std::fs::read_to_string(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut order = Vec::new();
    let mut map = HashMap::new();
    for (lineno, line) in data.lines().enumerate() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut fields = line.splitn(2, '\t');
        let chrom = fields.next().unwrap_or("").to_owned();
        let size_str = fields.next().unwrap_or("");
        let size: u64 = size_str.parse().map_err(|_| {
            RsomicsError::InvalidInput(format!(
                "genome file {}: line {}: bad size {:?}",
                path.display(),
                lineno + 1,
                size_str
            ))
        })?;
        if !map.contains_key(&chrom) {
            order.push(chrom.clone());
        }
        map.insert(chrom, size);
    }
    Ok((order, map))
}

/// Compute the complement of `input` BED against `genome`, writing to `output`.
///
/// Input must be sorted by (chrom, start) — the same requirement as
/// `bedtools complement`. Unsorted input gives undefined output (gap
/// coordinates may be wrong), matching bedtools' own behaviour.
pub fn complement(
    input: &Path,
    genome_order: &[String],
    genome: &HashMap<String, u64>,
    output: &mut dyn Write,
) -> Result<()> {
    let file = std::fs::File::open(input)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", input.display())))?;
    complement_reader(BufReader::new(file), genome_order, genome, output)
}

/// Same as [`complement`] but reads from stdin.
pub fn complement_stdin(
    genome_order: &[String],
    genome: &HashMap<String, u64>,
    output: &mut dyn Write,
) -> Result<()> {
    complement_reader(
        BufReader::new(std::io::stdin()),
        genome_order,
        genome,
        output,
    )
}

fn complement_reader<R: Read>(
    reader: BufReader<R>,
    genome_order: &[String],
    genome: &HashMap<String, u64>,
    output: &mut dyn Write,
) -> Result<()> {
    // Collect intervals per chrom.
    let mut by_chrom: HashMap<String, Vec<(u64, u64)>> = HashMap::new();

    for (lineno_0, line) in reader.lines().enumerate() {
        let line = line.map_err(RsomicsError::Io)?;
        let bytes = line.as_bytes();
        if bytes.is_empty()
            || bytes[0] == b'#'
            || bytes.starts_with(b"track")
            || bytes.starts_with(b"browser")
        {
            continue;
        }
        let lineno = lineno_0 + 1;
        let mut fields = line.splitn(3, '\t');
        let chrom = fields.next().unwrap_or("").to_owned();
        let start_str = fields.next().unwrap_or("");
        let end_str = {
            let s = fields.next().unwrap_or("");
            // end field may have trailing columns — take only up to next tab.
            s.split('\t').next().unwrap_or(s)
        };
        let start: u64 = start_str.parse().map_err(|_| {
            RsomicsError::InvalidInput(format!("line {lineno}: bad start {start_str:?}"))
        })?;
        let end: u64 = end_str.parse().map_err(|_| {
            RsomicsError::InvalidInput(format!("line {lineno}: bad end {end_str:?}"))
        })?;
        by_chrom.entry(chrom).or_default().push((start, end));
    }

    let mut out = BufWriter::new(output);

    // Emit gaps in genome-file order.
    for chrom in genome_order {
        let chrom_size = match genome.get(chrom) {
            Some(&s) => s,
            None => continue,
        };

        let ivs = by_chrom.get(chrom).map(|v| v.as_slice()).unwrap_or(&[]);

        let mut cursor = 0u64;
        for &(start, end) in ivs {
            if start > cursor {
                writeln!(out, "{chrom}\t{cursor}\t{start}").map_err(RsomicsError::Io)?;
            }
            cursor = cursor.max(end);
        }
        if cursor < chrom_size {
            writeln!(out, "{chrom}\t{cursor}\t{chrom_size}").map_err(RsomicsError::Io)?;
        }
    }
    out.flush().map_err(RsomicsError::Io)?;
    Ok(())
}
