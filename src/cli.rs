use std::io;
use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Origin, Section};

use rsomics_bed_complement::{complement, complement_stdin, read_genome};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-bed-complement", disable_help_flag = true)]
pub struct Cli {
    /// Input sorted BED (default: stdin)
    input: Option<PathBuf>,
    /// Chromosome sizes file (required; two-column chrom\tsize TSV)
    #[arg(short = 'g', long, value_name = "FILE")]
    genome: PathBuf,
    /// Output BED (default: stdout)
    #[arg(short = 'o', long)]
    output: Option<PathBuf>,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        let (genome_order, genome) = read_genome(&self.genome)?;

        let mut stdout_lock;
        let mut file_out;
        let out: &mut dyn io::Write = if let Some(ref p) = self.output {
            file_out = std::fs::File::create(p).map_err(rsomics_common::RsomicsError::Io)?;
            &mut file_out
        } else {
            stdout_lock = io::stdout().lock();
            &mut stdout_lock
        };

        match self.input {
            Some(ref p) => complement(p.as_path(), &genome_order, &genome, out),
            None => complement_stdin(&genome_order, &genome, out),
        }
    }
}

pub const HELP: HelpSpec = HelpSpec {
    name: META.name,
    version: META.version,
    tagline: "Compute the complement of a BED file — uncovered genome intervals (bedtools complement).",
    origin: Some(Origin {
        upstream: "bedtools",
        upstream_license: "MIT",
        our_license: "MIT OR Apache-2.0",
        paper_doi: Some("10.1093/bioinformatics/btq033"),
    }),
    usage_lines: &["[OPTIONS] -g <genome> [INPUT]"],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: Some('g'),
                long: "genome",
                aliases: &[],
                value: Some("<FILE>"),
                type_hint: Some("Path"),
                required: true,
                default: None,
                description: "Chromosome sizes file (chrom\\tsize TSV)",
                why_default: None,
            },
            FlagSpec {
                short: Some('o'),
                long: "output",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("Path"),
                required: false,
                default: Some("stdout"),
                description: "Output BED path",
                why_default: None,
            },
            FlagSpec {
                short: Some('h'),
                long: "help",
                aliases: &[],
                value: None,
                type_hint: Some("bool"),
                required: false,
                default: None,
                description: "Show this help",
                why_default: None,
            },
        ],
    }],
    examples: &[Example {
        description: "Find uncovered genomic regions",
        command: "rsomics-bed-complement -g chrom.sizes peaks.bed",
    }],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use clap::CommandFactory;
    #[test]
    fn cli_definition_is_valid() {
        super::Cli::command().debug_assert();
    }
}
