extern crate bio;
extern crate getopts;
extern crate rdxsort;
use crate::bio::io::fasta::FastaRead;
use bio::io::fasta::Reader as faReader;
use bio::io::fasta::Record as faRecord;
use getopts::Options;
use search_primer::counting_bloomfilter_util::aggregate_length_between_lr_tuple;
use search_primer::sequence_encoder_util::DnaSequence;
use std::fs::File;
use std::io::BufReader;
use std::io::{BufRead, Write};
use std::sync::Arc;
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
    opts.optopt("i", "input", "set input file name (TSV)", "NAME");
    opts.optopt("o", "output", "set output file name", "NAME");
    opts.optopt(
        "t",
        "thread",
        "number of threads to use for radix sort. default value is 8.",
        "THREAD",
    );
    opts.optopt("l", "length", "maximam size of target region", "LENGTH");

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

    let primer_filename: String = if matches.opt_present("i") {
        matches.opt_str("i").unwrap()
    } else {
        print_usage(&program, &opts);
        return;
    };

    let ngsread_input_file = if !matches.free.is_empty() {
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

    let max_offtarget_region_length: usize = if matches.opt_present("l") {
        matches.opt_str("l").unwrap().parse::<usize>().unwrap()
    } else {
        200
    };

    let output_file = if matches.opt_present("o") {
        matches.opt_str("o").unwrap()
    } else {
        "region.out".to_string()
    };

    eprintln!("Input  file name: {:?}", ngsread_input_file);
    eprintln!("Output file name: {:?}", output_file);
    eprintln!("Number of threads: {:?}", threads);

    let mut primer_tuple: Vec<(Vec<u8>, DnaSequence, DnaSequence)> = Vec::new();

    let f: File = File::open(&primer_filename).unwrap();
    let reader = BufReader::new(f);
    let mut lines = reader.lines();
    lines.next(); // ヘッダー行をスキップ

    lines.for_each(|line| {
        //eprintln!("line: {:?}", line);
        if let Ok(line) = line {
            // eprintln!("line: {:?}", line);
            let columns: Vec<&str> = line.split('\t').collect();
            assert!(columns.len() > 2, "column length is {}", columns.len());
            let primer_id: &[u8] = columns.get(0).unwrap_or(&"").as_bytes();
            let left_primer_seq: &str = columns.get(1).unwrap_or(&"");
            let right_primer_seq: &str = columns.get(2).unwrap_or(&"");
            let left_primer: DnaSequence = DnaSequence::new(&left_primer_seq.into());
            let right_primer: DnaSequence = DnaSequence::new(&right_primer_seq.into());
            let left_primer_revcomp: DnaSequence = left_primer.reverse_complement();
            let right_primer_revcomp: DnaSequence = right_primer.reverse_complement();
            let id_2: Vec<u8> =
                [&primer_id, &b"_LeftPrimerForward_RightPrimerRevcomp"[..]].concat();
            primer_tuple.push((id_2, left_primer.clone(), right_primer_revcomp.clone()));
        }
    });

    eprintln!("Number of primer_tuple: {:?}", &primer_tuple.len());
    for each_primer_tuple in &primer_tuple {
        eprint!(
            "primer_tuple: {}-",
            String::from_utf8(each_primer_tuple.1.decode(0, each_primer_tuple.1.len())).unwrap()
        );
        eprintln!(
            "{}",
            String::from_utf8(each_primer_tuple.2.decode(0, each_primer_tuple.2.len())).unwrap()
        );
    }
    let ngsread_file = File::open(&ngsread_input_file).expect("Error during opening the file");
    eprintln!("loading {:?}", &ngsread_input_file);
    let mut reader = faReader::new(ngsread_file);
    let mut record = faRecord::new();
    let mut sequences: Vec<DnaSequence> = Vec::new();
    eprintln!("loading {:?} done", &ngsread_input_file);
    'each_read: loop {
        reader.read(&mut record).unwrap();
        if record.is_empty() {
            break 'each_read;
        }
        let sequence_as_vec: Vec<u8> = record.seq().to_vec();
        let current_sequence = DnaSequence::new(&sequence_as_vec);
        sequences.push(current_sequence);
    }
    let chunk_size: usize = sequences.len() / (threads - 1);
    let sequences_ref = Arc::new(sequences);
    let primer_tuple_ref = Arc::new(primer_tuple);
    //let mut cbf_oyadama: Vec<u32> = vec![0;BLOOMFILTER_TABLE_SIZE];

    //let mut product_size_hashmap = HashMap::<u32, usize>::new();
    thread::scope(|scope| {
        let mut children_1 = Vec::new();
        for i in 1..threads {
            let sequences_ref = Arc::clone(&sequences_ref);
            let primer_tuple_ref = Arc::clone(&primer_tuple_ref);
            children_1.push(scope.spawn(move || {
                let start_idx: usize = (i - 1) * chunk_size;
                let end_idx: usize;
                if i != threads - 1 {
                    end_idx = i * chunk_size;
                } else {
                    end_idx = sequences_ref.len() - 1;
                }
                let primer_tuple_ref_mine = (*primer_tuple_ref).clone();
                let slice_sequences = Vec::from(sequences_ref[start_idx..end_idx].to_vec());

                eprintln!(
                    "start calling aggregate_length_between_primer_tuple[{}], # of sequence: {}",
                    i,
                    &slice_sequences.len()
                );
                let cbf: Vec<u8> = aggregate_length_between_lr_tuple(
                    &slice_sequences,
                    i,
                    &primer_tuple_ref_mine,
                    max_offtarget_region_length,
                );
                eprintln!(
                    "finish calling aggregate_length_between_primer_tuple[{}]",
                    i
                );
                return cbf;
            }))
        }
        eprintln!("start writing to output file: {:?}", &output_file);
        let mut file = File::create(&output_file).unwrap();
        for child in children_1 {
            let result = child.join(); // Resultを返さず、直接戻り値を取得
            file.write_all(&result.unwrap()).unwrap();
        }
        eprintln!("finish writing to output file: {:?}", &output_file);
    });
}
