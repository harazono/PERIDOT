use getopts::Options;
use search_primer::cbf_storage::build_and_save_cbf;
use search_primer::config::Config;
use std::{env, process};

fn print_usage(program: &str, opts: &Options) {
    let brief = format!("Usage: {} FASTA_FILE [options]", program);
    print!("{}", opts.usage(&brief));
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();

    let mut opts = Options::new();
    opts.optopt("o", "output", "CBF output file path", "CBF_FILE");
    opts.optopt("c", "config", "config file path", "CONFIG");
    opts.optopt("t", "threads", "number of threads (default: 8)", "THREADS");
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

    let output_file = if matches.opt_present("o") {
        matches.opt_str("o").unwrap()
    } else {
        format!("{}.cbf", input_file)
    };

    if let Err(e) = build_and_save_cbf(&input_file, &output_file, &config, threads) {
        eprintln!("Error building CBF: {}", e);
        process::exit(1);
    }

    println!("CBF built successfully: {}", output_file);
}
