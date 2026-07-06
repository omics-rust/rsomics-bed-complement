//! BED complement — intervals not covered by any input feature.
//!
//! Equivalent to `bedtools complement -i input -g genome`. For each chromosome
//! in the genome file, emits the uncovered gaps: [0, first_start),
//! [prev_end, next_start), …, [last_end, chrom_size).
//!
//! Chromosomes with no features in the input are emitted as a single interval
//! spanning the whole chromosome; likewise chromosomes present in the genome
//! file but absent from the input.
//!
//! Input must be coordinate-sorted in genome-file chromosome order, then by
//! ascending start — the same requirement `bedtools complement` enforces. Like
//! bedtools, we reject rather than silently mis-handle malformed input: input
//! whose records are out of that order, a record on a chromosome absent from
//! the genome file, `start > end`, a negative or non-numeric coordinate, or
//! fewer than three columns all fail loud with a non-zero exit.
//!
//! A zero-length feature (`start == end`) is virtually widened to
//! `[start-1, end+1)` before gap computation, matching bedtools.

use std::collections::{HashMap, HashSet};
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
/// Input must be coordinate-sorted in genome-file chromosome order then by
/// ascending start; out-of-order, out-of-genome, or malformed records fail loud
/// (see the module docs), matching `bedtools complement`.
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
    let gidx: HashMap<&str, usize> = genome_order
        .iter()
        .enumerate()
        .map(|(i, c)| (c.as_str(), i))
        .collect();

    let mut by_chrom: HashMap<String, Vec<(u64, u64)>> = HashMap::new();

    let mut cur: Option<(String, usize)> = None;
    let mut last_start = 0u64;
    let mut closed: HashSet<String> = HashSet::new();

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
        let chrom = fields.next().unwrap_or("");
        let start_str = fields.next().unwrap_or("");
        let end_str = fields.next().unwrap_or("").split('\t').next().unwrap_or("");
        let start: u64 = start_str.parse().map_err(|_| {
            RsomicsError::InvalidInput(format!("line {lineno}: bad start {start_str:?}"))
        })?;
        let end: u64 = end_str.parse().map_err(|_| {
            RsomicsError::InvalidInput(format!("line {lineno}: bad end {end_str:?}"))
        })?;
        if start > end {
            return Err(RsomicsError::InvalidInput(format!(
                "line {lineno}: start {start} greater than end {end}"
            )));
        }

        let &idx = gidx.get(chrom).ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "line {lineno}: chromosome {chrom:?} does not exist in the genome file"
            ))
        })?;

        match &cur {
            Some((c, _)) if c == chrom => {
                if start < last_start {
                    return Err(RsomicsError::InvalidInput(format!(
                        "line {lineno}: input not coordinate-sorted (start {start} precedes {last_start} on {chrom})"
                    )));
                }
                last_start = start;
            }
            _ => {
                if closed.contains(chrom) {
                    return Err(RsomicsError::InvalidInput(format!(
                        "line {lineno}: input not coordinate-sorted (chromosome {chrom} reappears)"
                    )));
                }
                if let Some((c, i)) = cur.replace((chrom.to_owned(), idx)) {
                    if idx < i {
                        return Err(RsomicsError::InvalidInput(format!(
                            "line {lineno}: chromosome {chrom} out of genome-file order"
                        )));
                    }
                    closed.insert(c);
                }
                last_start = start;
            }
        }

        by_chrom
            .entry(chrom.to_owned())
            .or_default()
            .push((start, end));
    }

    let mut out = BufWriter::new(output);

    for chrom in genome_order {
        let &chrom_size = genome.get(chrom).unwrap();
        let ivs = by_chrom.get(chrom).map(|v| v.as_slice()).unwrap_or(&[]);

        let mut cursor = 0u64;
        for &(start, end) in ivs {
            // bedtools virtually widens a zero-length feature to [start-1, end+1).
            let (start, end) = if start == end {
                (start.saturating_sub(1), end + 1)
            } else {
                (start, end)
            };
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
