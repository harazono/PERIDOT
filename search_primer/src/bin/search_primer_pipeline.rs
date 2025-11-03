use getopts::Options;
use search_primer::cbf_storage::{build_and_save_cbf, load_cbf_and_extract_candidates, remove_false_positives};
use search_primer::config::Config;
use search_primer::sequence_encoder_util::decode_u128_2_dna_seq;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::Path;
use std::{env, process};

fn print_usage(program: &str, opts: &Options) {
    let brief = format!("Usage: {} FASTA_FILE [options]", program);
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
    opts.optopt("", "cbf", "CBF file path (if exists, skip building)", "CBF_FILE");
    opts.optflag("", "save-cbf", "save CBF for reuse");
    opts.optflag("b", "binary", "output binary format");
    opts.optflag("r", "only-num", "output only count");
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

    let input_file = if !matches.free.is_empty() {
        matches.free[0].clone()
    } else {
        print_usage(&program, &opts);
        process::exit(1);
    };

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
        format!("{}_threshold{}_threads{}.out", input_file, threshold, threads)
    };

    // Determine CBF file path
    let cbf_file = if matches.opt_present("cbf") {
        matches.opt_str("cbf").unwrap()
    } else {
        format!("{}.cbf", input_file)
    };

    // Step 1: Build or load CBF
    if !Path::new(&cbf_file).exists() || matches.opt_present("save-cbf") {
        eprintln!("Building CBF...");
        if let Err(e) = build_and_save_cbf(&input_file, &cbf_file, &config, threads) {
            eprintln!("Error building CBF: {}", e);
            process::exit(1);
        }
    } else {
        eprintln!("Using existing CBF: {}", cbf_file);
    }

    // Step 2: Extract candidates
    eprintln!("Extracting candidates...");
    let candidates = match load_cbf_and_extract_candidates(&cbf_file, &input_file, threshold, &config, threads) {
        Ok(c) => c,
        Err(e) => {
            eprintln!("Error extracting candidates: {}", e);
            process::exit(1);
        }
    };

    // Step 3: Remove false positives
    eprintln!("Removing false positives...");
    let final_results = match remove_false_positives(&candidates, &input_file, threshold, &config, threads) {
        Ok(r) => r,
        Err(e) => {
            eprintln!("Error removing false positives: {}", e);
            process::exit(1);
        }
    };

    // Step 4: Output results
    let file = File::create(&output_file).unwrap();
    let mut writer = BufWriter::new(file);

    if matches.opt_present("r") {
        writeln!(writer, "lr_tuple count: {}\tthreshold: {}\tinput file: {}", 
                final_results.len(), threshold, input_file).unwrap();
    } else if matches.opt_present("b") {
        // Binary output
        for &lr_tuple in final_results.keys() {
            let bytes = lr_tuple.to_le_bytes();
            writer.write_all(&bytes).unwrap();
        }
    } else {
        // Text output with DNA sequences
        for &lr_tuple in final_results.keys() {
            let dna_seq = decode_u128_2_dna_seq(&lr_tuple, config.total_segment_len());
            writeln!(writer, "{}", String::from_utf8(dna_seq).unwrap()).unwrap();
        }
    }

    writer.flush().unwrap();
    
    eprintln!("Pipeline completed successfully!");
    eprintln!("Final results: {} LR-tuples", final_results.len());
    eprintln!("L: {}\tR: {}\tthreshold: {}\tthreads: {}", 
             config.l_len, config.r_len, threshold, threads);
    println!("Results saved to: {}", output_file);
}
