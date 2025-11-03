extern crate bio;
extern crate getopts;
extern crate rdxsort;
use crate::bio::io::fasta::FastaRead;
use bio::io::fasta::Reader as faReader;
use bio::io::fasta::Record as faRecord;
use getopts::Options;
use search_primer::counting_bloomfilter_util::BLOOMFILTER_TABLE_SIZE;
use search_primer::counting_bloomfilter_util::{
    build_counting_bloom_filter, count_lr_tuple_with_hashtable, number_of_high_occurence_lr_tuple,
};
use search_primer::counting_bloomfilter_util::{HASHSET_SIZE, L_LEN, R_LEN};
use search_primer::sequence_encoder_util::decode_u128_2_dna_seq;
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
    opts.optopt("x", "output1", "set output file name(by CBF)", "NAME");
    opts.optopt("y", "output2", "set output file name(by HS)", "NAME");
    opts.optopt(
        "t",
        "thread",
        "number of threads to use for radix sort. default value is 8.",
        "THREAD",
    );
    opts.optopt(
        "a",
        "threshold",
        "threshold for hyper log counter. default value is 8.",
        "THRESHOLD",
    );
    opts.optflag("b", "binary", "outputs binary file");
    opts.optflag("r", "only-num", "outputs only total number of lr-tuple");
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

    let input_file = if !matches.free.is_empty() {
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

    let threshold: u16 = if matches.opt_present("a") {
        matches.opt_str("a").unwrap().parse::<u16>().unwrap()
    } else {
        8
    };

    let output_file_1 = if matches.opt_present("output1") {
        matches.opt_str("output1").unwrap()
    } else {
        format!(
            "{:?}_threshold{}_threads{}_CBF.out",
            input_file, threshold, threads
        )
    };

    let output_file_2 = if matches.opt_present("output2") {
        matches.opt_str("output2").unwrap()
    } else {
        format!(
            "{:?}_threshold{}_threads{}_HS.out",
            input_file, threshold, threads
        )
    };

    eprintln!("input  file: {:?}", input_file);
    let file = File::open(&input_file).expect("Error during opening the file");
    let mut reader = faReader::new(file);
    let mut record = faRecord::new();
    let mut sequences: Vec<DnaSequence> = Vec::new();
    eprintln!("loading {:?} done", input_file);
    'each_read: loop {
        reader.read(&mut record).unwrap();
        if record.is_empty() {
            break 'each_read;
        }
        let sequence_as_vec: Vec<u8> = record.seq().to_vec();
        let current_sequence = DnaSequence::new(&sequence_as_vec);
        sequences.push(current_sequence);
    }

    /*
    ここにマルチスレッド処理を書く
    */

    let chunk_size: usize = max(sequences.len() / (threads - 1), 1);
    let sequences_ref = &sequences;
    let number_of_threads: usize = min(threads, sequences.len() + 1);
    let mut cbf_oyadama: Vec<u16> = vec![0; BLOOMFILTER_TABLE_SIZE];
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
                    BLOOMFILTER_TABLE_SIZE,
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
    let h_cbf_h_oyadama: Arc<Mutex<HashSet<u128>>> =
        Arc::new(Mutex::new(HashSet::with_capacity(HASHSET_SIZE)));
    let cbf_oyadama_ref = &cbf_oyadama;
    let h_cbf_h_oyadama_ref = &h_cbf_h_oyadama;

    eprintln!("start calling number_of_high_occurence_lr_tuple");
    thread::scope(|scope| {
        let mut children_2 = Vec::new();
        for i in 1..number_of_threads {
            children_2.push(scope.spawn(move || {
                let start_idx: usize = (i - 1) * chunk_size;
                let end_idx: usize;
                if i != number_of_threads - 1 {
                    end_idx = i * chunk_size;
                } else {
                    end_idx = sequences_ref.len();
                }
                eprintln!(
                    "thread [{}]: start calling number_of_high_occurence_lr_tuple",
                    i
                );
                let h_cbf_h: HashSet<u128> = number_of_high_occurence_lr_tuple(
                    cbf_oyadama_ref,
                    sequences_ref,
                    start_idx,
                    end_idx,
                    HASHSET_SIZE,
                    threshold,
                    BLOOMFILTER_TABLE_SIZE,
                    i,
                    0,
                );
                h_cbf_h_oyadama_ref.lock().unwrap().extend(&h_cbf_h);
                eprintln!(
                    "thread [{}]: finish calling number_of_high_occurence_lr_tuple",
                    i
                );
            }))
        }
        for child in children_2 {
            let _ = child.join();
        }
    });

    let mut high_occurence_lr_tuple: Vec<u128> =
        Vec::from_iter(h_cbf_h_oyadama.lock().unwrap().clone());
    high_occurence_lr_tuple.sort();
    let high_occurence_lr_tuple_ref: &HashSet<u128> =
        &HashSet::from_iter(high_occurence_lr_tuple.clone());
    let hashtable_count_result: HashMap<u128, u16> = count_lr_tuple_with_hashtable(
        &sequences,
        0,
        sequences.len(),
        high_occurence_lr_tuple_ref,
        1,
        0,
    );

    let mut sorted_hs_list: Vec<u128> = hashtable_count_result
        .iter()
        .filter_map(|(&key, &value)| if value >= threshold { Some(key) } else { None })
        .collect();
    sorted_hs_list.sort();

    eprintln!("length of sorted_hs_list: {}", sorted_hs_list.len());

    let mut w1 = BufWriter::new(fs::File::create(&output_file_1).unwrap());
    let mut w2 = BufWriter::new(fs::File::create(&output_file_2).unwrap());

    let mut previous_lr_tuple: u128 = 0;
    let mut cnt = 0;
    let mut buf_array: [u8; 16] = [0; 16];
    let mut buf_num: u128;

    if matches.opt_present("r") {
        eprintln!(
            "matches.opt_present('r'): {}\tmatches.opt_present('b'): {}",
            matches.opt_present("r"),
            matches.opt_present("b")
        );
        for each_lr_tuple in &high_occurence_lr_tuple {
            if previous_lr_tuple != *each_lr_tuple {
                cnt += 1;
            }
            previous_lr_tuple = *each_lr_tuple;
        }
        writeln!(
            &mut w1,
            "lr_tuple count: {}\tthreshold: {}\tinput file {:?}",
            cnt, threshold, &input_file
        )
        .unwrap();
    }
    if !matches.opt_present("r") && matches.opt_present("b") {
        eprintln!(
            "matches.opt_present('r'): {}\tmatches.opt_present('b'): {}",
            matches.opt_present("r"),
            matches.opt_present("b")
        );
        for each_lr_tuple in &high_occurence_lr_tuple {
            if previous_lr_tuple != *each_lr_tuple {
                cnt += 1;
                buf_num = *each_lr_tuple;
                for i in 0..16 {
                    buf_array[15 - i] = u8::try_from(buf_num & 0xFF).unwrap();
                    buf_num >>= 8;
                }
                w1.write(&buf_array).unwrap();
            }
            previous_lr_tuple = *each_lr_tuple;
        }
        for each_sorted_hs_list in &sorted_hs_list {
            buf_num = *each_sorted_hs_list;
            for i in 0..16 {
                buf_array[15 - i] = u8::try_from(buf_num & 0xFF).unwrap();
                buf_num >>= 8;
            }
            w2.write(&buf_array).unwrap();
        }
    }

    if !matches.opt_present("r") && !matches.opt_present("b") {
        eprintln!(
            "matches.opt_present('r'): {}\tmatches.opt_present('b'): {}",
            matches.opt_present("r"),
            matches.opt_present("b")
        );
        for each_lr_tuple in &high_occurence_lr_tuple {
            if previous_lr_tuple != *each_lr_tuple {
                cnt += 1;
                writeln!(
                    &mut w1,
                    "{:?}",
                    String::from_utf8(decode_u128_2_dna_seq(&each_lr_tuple, 54)).unwrap()
                )
                .unwrap();
            }
            previous_lr_tuple = *each_lr_tuple;
        }
    }

    eprintln!("finish writing to output file: {:?}", &output_file_1);
    eprint!(
        "L: {}\tR: {}\tthreshold:{}\tcardinarity: {}\t",
        L_LEN, R_LEN, threshold, cnt
    );
    eprintln!("threads: {}\tinput file {:?}", threads, &input_file);
}
