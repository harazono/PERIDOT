extern crate bio;
extern crate getopts;
extern crate rdxsort;
use crate::bio::io::fasta::FastaRead;
use bio::io::fasta::Reader as faReader;
use bio::io::fasta::Record as faRecord;
use getopts::Options;
use search_primer::counting_bloomfilter_util::BLOOMFILTER_TABLE_SIZE;
use search_primer::counting_bloomfilter_util::HASHSET_SIZE;
use search_primer::counting_bloomfilter_util::{
    build_counting_bloom_filter, count_lr_tuple_with_hashtable, number_of_high_occurence_lr_tuple,
};
use search_primer::sequence_encoder_util::DnaSequence;
use std::cmp::{max, min};
use std::collections::HashMap;
use std::collections::HashSet;
use std::fs;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::iter::zip;
use std::sync::{Arc, Mutex};
use std::thread;
use std::{env, process};

fn print_usage(program: &str, opts: &Options) {
    let brief = format!("Usage: {} FILE [options]", program);
    print!("{}", opts.usage(&brief));
    process::exit(0);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();

    let mut opts = Options::new();
    opts.optopt("x", "output", "set output file name(by CBF)", "NAME");
    opts.optopt(
        "t",
        "thread",
        "number of threads to use for radix sort. default value is 8.",
        "THREAD",
    );
    opts.optopt(
        "m",
        "bucket-size",
        "bucket size of CBF. default value is 1 << 30 as usize.",
        "BUCKET_SIZE",
    );

    opts.optflag("h", "help", "print this help menu");

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => {
            panic!("{}", f.to_string())
        }
    };
    if matches.opt_present("h") {
        print_usage(&program, &opts);
        return;
    }

    let input_file: String = if !matches.free.is_empty() {
        matches.free[0].clone()
    } else {
        print_usage(&program, &opts);
        return;
    };

    let threads: usize = if matches.opt_present("t") {
        matches.opt_str("t").unwrap().parse::<usize>().unwrap()
    } else {
        8
    };
    let bucket_size: usize = if matches.opt_present("m") {
        matches.opt_str("m").unwrap().parse::<usize>().unwrap()
    } else {
        BLOOMFILTER_TABLE_SIZE
    };

    let output_file: String = if matches.opt_present("output") {
        matches.opt_str("output").unwrap()
    } else {
        format!("{:?}_threads{}_CBF.out", input_file, threads)
    };

    eprintln!(
        "input  file: {:?}\t HASHSET_SIZE: {}\tBLOOMFILTER_TABLE_SIZE: {}",
        input_file, HASHSET_SIZE, bucket_size
    );
    let file: File = File::open(&input_file).expect("Error during opening the file");
    let mut reader: faReader<std::io::BufReader<File>> = faReader::new(file);
    let mut record: faRecord = faRecord::new();
    let mut sequences: Vec<DnaSequence> = Vec::new();
    eprintln!("loading {:?} done", input_file);
    'each_read: loop {
        reader.read(&mut record).unwrap();
        if record.is_empty() {
            break 'each_read;
        }
        let sequence_as_vec: Vec<u8> = record.seq().to_vec();
        let current_sequence: DnaSequence = DnaSequence::new(&sequence_as_vec);
        sequences.push(current_sequence);
    }

    /*
    ここにマルチスレッド処理を書く
    */

    let chunk_size: usize = max(sequences.len() / (threads - 1), 1);
    let sequences_ref = &sequences;
    let number_of_threads: usize = min(threads, sequences.len() + 1);
    let mut cbf_oyadama: Vec<u16> = vec![0; bucket_size];
    thread::scope(|scope| {
        let mut children_1 = Vec::new();
        for i in 1..number_of_threads {
            children_1.push(scope.spawn(move || {
                let start_idx: usize = (i - 1) * chunk_size;
                let end_idx: usize;
                if i != number_of_threads - 1 {
                    end_idx = i * chunk_size;
                } else {
                    end_idx = sequences_ref.len();
                }
                eprintln!(
                    "start calling build_counting_bloom_filter[{}], {}-{}",
                    i, start_idx, end_idx
                );
                let cbf: Vec<u16> = build_counting_bloom_filter(
                    sequences_ref,
                    start_idx,
                    end_idx,
                    bucket_size,
                    i,
                    0,
                );
                eprintln!(
                    "finish calling build_counting_bloom_filter[{}], {}-{}",
                    i, start_idx, end_idx
                );
                cbf
            }))
        }
        for child in children_1 {
            let cbf = child.join().unwrap();
            zip(cbf_oyadama.iter_mut(), cbf)
                .for_each(|(x, y)| *x = x.checked_add(y).unwrap_or(u16::MAX));
        }
    });
    eprintln!("finish building counting bloom filter");

    let mut cbf_hist: [u16; 65535] = [0u16; 65535];
    for each_bucket in cbf_oyadama {
        cbf_hist[each_bucket as usize] += 1;
    }

    let mut w: BufWriter<File> = BufWriter::new(fs::File::create(&output_file).unwrap());
    for each_cbf_hist in &cbf_hist {
        writeln!(&mut w, "{}", each_cbf_hist).unwrap();
    }

    eprintln!("finish writing to output file: {:?}", &output_file);
    eprintln!("threads: {}\tinput file {:?}", threads, &input_file);
}
