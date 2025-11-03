use getopts::Options;
use search_primer::cbf_storage::load_cbf_and_extract_candidates;
use search_primer::config::Config;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::{env, process};

fn print_usage(program: &str, opts: &Options) {
    let brief = format!("Usage: {} CBF_FILE FASTA_FILE [options]", program);
    print!("{}", opts.usage(&brief));
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();

    let mut opts = Options::new();
    opts.optopt("o", "output", "output file path", "OUTPUT_FILE");
    opts.optopt("c", "config", "config file path", "CONFIG");
    opts.optopt("t", "threads", "number of threads (default: 8)", "THREADS");
    opts.optopt("a", "threshold", "threshold (default: 1000)", "THRESHOLD");
    opts.optflag("b", "binary", "output binary format");
    opts.optflag("h", "help", "print this help menu");

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => {
            eprintln!("Error: {}", f);
            process::exit(1);
        }
    };

    if matches.opt_present("h") {
        print_usage(&program, &opts);
        return;
    }

    if matches.free.len() < 2 {
        print_usage(&program, &opts);
        process::exit(1);
    }

    let cbf_file = &matches.free[0];
    let fasta_file = &matches.free[1];

    let config = if matches.opt_present("c") {
        let config_path = matches.opt_str("c").unwrap();
        Config::from_file(&config_path).unwrap_or_else(|e| {
            eprintln!("Failed to load config file: {}", e);
            process::exit(1);
        })
    } else {
        Config::default()
    };

    let threads = if matches.opt_present("t") {
        matches.opt_str("t").unwrap().parse::<usize>().unwrap_or(8)
    } else {
        8
    };

    let threshold = if matches.opt_present("a") {
        matches.opt_str("a").unwrap().parse::<u16>().unwrap_or(1000)
    } else {
        1000
    };

    let output_file = if matches.opt_present("o") {
        matches.opt_str("o").unwrap()
    } else {
        format!("{}_candidates_t{}.out", fasta_file, threshold)
    };

    let candidates = match load_cbf_and_extract_candidates(cbf_file, fasta_file, threshold, &config, threads) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("Error extracting candidates: {}", e);
            process::exit(1);
        }
    };

    let file = File::create(&output_file).unwrap();
    let mut writer = BufWriter::new(file);

    if matches.opt_present("b") {
        // Binary output
        for &candidate in &candidates {
            let bytes = candidate.to_le_bytes();
            writer.write_all(&bytes).unwrap();
        }
    } else {
        // Text output - just count for now
        writeln!(writer, "Found {} candidates with threshold {}", candidates.len(), threshold).unwrap();
        for &candidate in &candidates {
            writeln!(writer, "{:032x}", candidate).unwrap();
        }
    }

    writer.flush().unwrap();
    println!("Candidates extracted: {} -> {}", candidates.len(), output_file);
}
