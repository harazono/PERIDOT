extern crate bio;
extern crate getopts;
use crate::bio::io::fasta::FastaRead;
use bio::io::fasta::Reader as faReader;
use bio::io::fasta::Record as faRecord;
use getopts::Options;
use search_primer::config::Config;
use search_primer::counting_bloomfilter_util::{
    build_counting_bloom_filter, count_lr_tuple_with_hashtable, number_of_high_occurence_lr_tuple,
};
use search_primer::sequence_encoder_util::decode_u128_2_dna_seq;
use search_primer::sequence_encoder_util::DnaSequence;
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
    let brief: String = format!("Usage: {} FILE [options]", program);
    print!("{}", opts.usage(&brief));
    process::exit(0);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();

    let mut opts = Options::new();
    opts.optopt("o", "output", "set output file name", "NAME");
    opts.optopt("c", "config", "config file path", "CONFIG");
    opts.optopt(
        "t",
        "thread",
        "number of threads to use for radix sort. default value is 8.",
        "THREAD",
    );
    opts.optopt(
        "a",
        "threshold",
        "threshold of occurence. default value is 1000.",
        "THRESHOLD",
    );
    opts.optflag("b", "binary", "outputs binary file");
    opts.optflag("r", "only-num", "outputs only total number of lr-tuple.");
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

    let config = if matches.opt_present("c") {
        let config_path = matches.opt_str("c").unwrap();
        Config::from_file(&config_path).unwrap_or_else(|e| {
            eprintln!("Failed to load config file: {}", e);
            process::exit(1);
        })
    } else {
        Config::default()
    };

    let threads: usize = if matches.opt_present("t") {
        std::cmp::max(matches.opt_str("t").unwrap().parse::<usize>().unwrap(), 2)
    } else {
        8
    };

    let threshold: u16 = if matches.opt_present("a") {
        matches.opt_str("a").unwrap().parse::<u16>().unwrap()
    } else {
        1000
    };

    let output_file: String = if matches.opt_present("o") {
        matches.opt_str("o").unwrap()
    } else {
        format!(
            "{:?}_threshold{}_threads{}.out",
            input_file, threshold, threads
        )
    };
    let mut w: BufWriter<File> = BufWriter::new(fs::File::create(&output_file).unwrap());

    eprintln!("input  file: {:?}", input_file);
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
    // CBFをマルチスレッドで作成する
    let chunk_size: usize = sequences.len() / (threads - 1);
    let sequences_ref: &Vec<DnaSequence> = &sequences;
    let mut cbf_oyadama: Vec<u16> = vec![0; config.bloomfilter_table_size];

    thread::scope(|scope: &thread::Scope<'_, '_>| {
        let mut children_1: Vec<thread::ScopedJoinHandle<'_, Vec<u16>>> = Vec::new();
        for i in 1..threads {
            children_1.push(scope.spawn(move || {
                let start_idx: usize = (i - 1) * chunk_size;
                let end_idx: usize;
                if i != threads - 1 {
                    end_idx = i * chunk_size;
                } else {
                    end_idx = sequences_ref.len() - 1;
                }
                eprintln!(
                    "start calling build_counting_bloom_filter[{}], {}-{}",
                    i, start_idx, end_idx
                );
                // CBFの構築
                let cbf: Vec<u16> = build_counting_bloom_filter(
                    sequences_ref,
                    start_idx,
                    end_idx,
                    &config,
                    i,
                );
                eprintln!(
                    "finish calling build_counting_bloom_filter[{}], {}-{}",
                    i, start_idx, end_idx
                );
                cbf
            }))
        }
        for child in children_1 {
            let cbf: Vec<u16> = child.join().unwrap();
            assert!(cbf.len() == config.bloomfilter_table_size);
            assert!(cbf_oyadama.len() == config.bloomfilter_table_size);
            zip(cbf_oyadama.iter_mut(), cbf)
                .for_each(|(x, y)| *x = x.checked_add(y).unwrap_or(u16::MAX));
        }
    });
    // CBFをマージする
    let h_cbf_h_oyadama: Arc<Mutex<HashSet<u128>>> =
        Arc::new(Mutex::new(HashSet::with_capacity(config.hashset_size)));

    //CBFを用いて高頻度のLR-tupleをマルチスレッドで列挙する
    let cbf_oyadama_ref: &Vec<u16> = &cbf_oyadama;
    let h_cbf_h_oyadama_ref: &Arc<Mutex<HashSet<u128>>> = &h_cbf_h_oyadama;
    thread::scope(|scope: &thread::Scope<'_, '_>| {
        let mut children_2: Vec<thread::ScopedJoinHandle<'_, ()>> = Vec::new();
        for i in 1..threads {
            children_2.push(scope.spawn(move || {
                let start_idx: usize = (i - 1) * chunk_size;
                let end_idx: usize;
                if i != threads - 1 {
                    end_idx = i * chunk_size;
                } else {
                    end_idx = sequences_ref.len() - 1;
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
                    threshold,
                    &config,
                    i,
                );
                h_cbf_h_oyadama_ref.lock().unwrap().extend(&h_cbf_h);
                eprintln!(
                    "thread [{}]: finish calling number_of_high_occurence_lr_tuple. h_cbf_h.len: {:?}",
                    i, &h_cbf_h.len()
                );
            }))
        }
        for child in children_2 {
            let _ = child.join();
        }
    });
    //高頻度のLR-tupleをマージする
    /*
        let mut high_occurence_lr_tuple: Vec<u128> =
            Vec::from_iter(h_cbf_h_oyadama.lock().unwrap().clone());
        high_occurence_lr_tuple.sort();

    */
    //高頻度のLR-tupleをハッシュテーブルを用いて数え直し、偽陽性を除去する
    let hashtable_count_result_oyadama: Arc<Mutex<HashMap<u128, u16>>> =
        Arc::new(Mutex::new(HashMap::with_capacity(config.hashset_size)));
    let hashtable_count_result_ref: &Arc<Mutex<HashMap<u128, u16>>> =
        &hashtable_count_result_oyadama;

    let high_occurence_lr_tuple: &HashSet<u128> = &*h_cbf_h_oyadama.lock().unwrap();
    eprintln!(
        "length of high_occurence_lr_tuple: {:?}",
        high_occurence_lr_tuple.len()
    );
    thread::scope(|scope| {
        let mut children_3 = Vec::new();
        for i in 1..threads {
            children_3.push(scope.spawn(move || {
                let start_idx: usize = (i - 1) * chunk_size;
                let end_idx: usize;
                if i != threads - 1 {
                    end_idx = i * chunk_size;
                } else {
                    end_idx = sequences_ref.len() - 1;
                }
                eprintln!(
                    "thread [{}]: start calling count_lr_tuple_with_hashtable",
                    i
                );
                let high_freq_ht: HashMap<u128, u16> = count_lr_tuple_with_hashtable(
                    &sequences_ref,
                    start_idx,
                    end_idx,
                    &high_occurence_lr_tuple,
                    &config,
                    i,
                );
                let mut hashtable_count_result: std::sync::MutexGuard<'_, HashMap<u128, u16>> =
                    hashtable_count_result_ref.lock().unwrap();
                for (key, &value) in &high_freq_ht {
                    *hashtable_count_result.entry(*key).or_insert(0) += value;
                }
                eprintln!(
                    "thread [{}]: finish calling count_lr_tuple_with_hashtable. high_freq_ht.len: {:?}",
                    &high_freq_ht.len(), i
                );
            }))
        }
        for child in children_3 {
            let _ = child.join();
        }
    });
    eprintln!(
        "count_lr_tuple_with_hashtable done: hashtable_count_result_ref.len() = {:?}",
        hashtable_count_result_ref.lock().unwrap().len()
    );

    let hashtable: std::sync::MutexGuard<'_, HashMap<u128, u16>> =
        hashtable_count_result_oyadama.lock().unwrap();

    if let Some((&key, &max_value)) = &hashtable.iter().max_by_key(|(_, &value)| value) {
        println!("Key with the maximum value: {}, Value: {}", key, max_value);
    } else {
        println!("The hashtable is empty.");
    }

    let mut sorted_hs_list: Vec<u128> = hashtable
        .iter()
        .filter_map(|(&key, &value)| if value >= threshold { Some(key) } else { None })
        .collect();
    sorted_hs_list.sort();

    let mut buf_array: [u8; 16] = [0; 16];
    let mut buf_num: u128;

    if matches.opt_present("r") {
        eprintln!(
            "matches.opt_present('r'): {}\tmatches.opt_present('b'): {}",
            matches.opt_present("r"),
            matches.opt_present("b")
        );
        writeln!(
            &mut w,
            "lr_tuple count: {}\tthreshold: {}\tinput file {:?}",
            sorted_hs_list.len(),
            threshold,
            &input_file
        )
        .unwrap();
    }
    if !matches.opt_present("r") && matches.opt_present("b") {
        eprintln!(
            "matches.opt_present('r'): {}\tmatches.opt_present('b'): {}",
            matches.opt_present("r"),
            matches.opt_present("b")
        );
        for each_lr_tuple in &sorted_hs_list {
            buf_num = *each_lr_tuple;
            for i in 0..16 {
                buf_array[15 - i] = u8::try_from(buf_num & 0xFF).unwrap();
                buf_num >>= 8;
            }
            w.write(&buf_array).unwrap();
        }
    }
    if !matches.opt_present("r") && !matches.opt_present("b") {
        eprintln!(
            "matches.opt_present('r'): {}\tmatches.opt_present('b'): {}",
            matches.opt_present("r"),
            matches.opt_present("b")
        );
        for each_lr_tuple in &sorted_hs_list {
            writeln!(
                &mut w,
                "{:?}",
                String::from_utf8(decode_u128_2_dna_seq(&each_lr_tuple, config.total_segment_len())).unwrap()
            )
            .unwrap();
        }
    }

    eprintln!("finish writing to output file: {:?}", &output_file);
    eprint!(
        "L: {}\tR: {}\tthreshold:{}\tcardinarity: {}\t",
        config.l_len,
        config.r_len,
        threshold,
        sorted_hs_list.len(),
    );
    eprintln!("threads: {}\tinput file {:?}", threads, &input_file);
}
