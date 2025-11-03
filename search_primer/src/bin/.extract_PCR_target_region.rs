extern crate bio;
extern crate getopts;
extern crate rdxsort;
use crate::bio::io::fasta::FastaRead;
use bio::io::fasta::Reader as faReader;
use bio::io::fasta::Record as faRecord;
use core::any::Any;
use getopts::Options;
use search_primer::counting_bloomfilter_util::aggregate_length_between_lr_tuple;
use search_primer::sequence_encoder_util::DnaSequence;
// use std::any::Any;
use std::env;
use std::fs::File;
use std::io::BufReader;
use std::io::{BufRead, Write};
use std::sync::Arc;
use std::thread;

fn print_usage(program: &str, opts: &Options) {
    let brief = format!("Usage: {}", program);
    print!("{}", opts.usage(&brief));
}
fn main() -> Result<(), String> {
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();

    let mut opts = Options::new();
    opts.optopt("o", "output", "set output file name", "NAME");
    opts.optopt("r", "read", "set NGS read file name", "READ");
    opts.optopt("t", "tsv", "set TSV file name", "TSV");

    opts.optopt(
        "T",
        "thread",
        "number of threads to use for processing. Default is 8.",
        "THREAD",
    );
    opts.optopt("l", "length", "maximum size of target region", "LENGTH");

    opts.optflag("h", "help", "print this help menu");

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => {
            panic!("{}", f.to_string())
        }
    };

    if matches.opt_present("h") {
        print_usage(&program, &opts);
        return Ok(());
    }

    let output_file = matches
        .opt_str("o")
        .unwrap_or_else(|| "output.out".to_string());
    let ngsread_file = matches.opt_str("r").ok_or_else(|| {
        print_usage(&program, &opts);
        "NGS read file name is required".to_string()
    })?;
    let tsv_file = matches.opt_str("t").ok_or_else(|| {
        print_usage(&program, &opts);
        "TSV file name is required".to_string()
    })?;
    let threads: usize = matches
        .opt_str("T")
        .unwrap_or("8".to_string())
        .parse()
        .unwrap();
    let max_offtarget_region_length: usize = matches
        .opt_str("l")
        .unwrap_or("200".to_string())
        .parse()
        .unwrap();

    eprintln!("Output file name: {:?}", output_file);
    eprintln!("NGS read file name: {:?}", ngsread_file);
    eprintln!("TSV file name: {:?}", tsv_file);
    eprintln!("Number of threads: {:?}", threads);
    eprintln!(
        "Maximum target region length: {:?}",
        max_offtarget_region_length
    );

    let mut primer_tuple: Vec<(Vec<u8>, DnaSequence, DnaSequence)> = Vec::new();

    let f: File = File::open(&tsv_file).unwrap();
    let reader: BufReader<File> = BufReader::new(f);
    let mut lines: std::io::Lines<BufReader<File>> = reader.lines();
    lines.next(); // ヘッダー行をスキップ

    lines.for_each(|line| {
        //eprintln!("line: {:?}", line);
        if let Ok(line) = line {
            // eprintln!("line: {:?}", line);
            let columns: Vec<&str> = line.split('\t').collect();
            let primer_id: &[u8] = columns.get(0).unwrap_or(&"").as_bytes();
            let left_primer_seq: &str = columns.get(1).unwrap_or(&"");
            let right_primer_seq: &str = columns.get(2).unwrap_or(&"");
            /*
            入力TSVのLプライマーはリードと同じstrand、Rはrevcompとする。
            R側のプライマーの読み込み時にrevcompに変換する。
            リードの座標を基準に5to3か3to5か決める。
            5'---> 3'
             */
            let left_primer_5to3: DnaSequence = DnaSequence::new(&left_primer_seq.into());
            let right_primer_3to5: DnaSequence =
                DnaSequence::new(&right_primer_seq.into()).reverse_complement();
            let left_primer_3to5: DnaSequence = left_primer_5to3.reverse_complement();
            let right_primer_5to3: DnaSequence = right_primer_3to5.reverse_complement();
            let id_1: Vec<u8> = [&primer_id, &b"L-L"[..]].concat();
            let id_2: Vec<u8> = [&primer_id, &b"L-R"[..]].concat();
            let id_3: Vec<u8> = [&primer_id, &b"R-L"[..]].concat();
            let id_4: Vec<u8> = [&primer_id, &b"R-R"[..]].concat();

            primer_tuple.push((id_1, left_primer_5to3.clone(), left_primer_3to5.clone()));
            primer_tuple.push((id_2, left_primer_5to3.clone(), right_primer_3to5.clone()));
            primer_tuple.push((id_3, right_primer_5to3.clone(), left_primer_3to5.clone()));
            primer_tuple.push((id_4, right_primer_5to3.clone(), right_primer_3to5.clone()));
        }
    });

    eprintln!("Number of primer_tuple: {:?}", &primer_tuple.len());
    /*
    for each_primer_tuple in &primer_tuple {
        eprintln!(
            "primer_tuple: {}",
            String::from_utf8(each_primer_tuple.1.decode(0, each_primer_tuple.1.len())).unwrap()
        );
        eprintln!(
            "primer_tuple: {}",
            String::from_utf8(each_primer_tuple.2.decode(0, each_primer_tuple.2.len())).unwrap()
        );
    }
    */
    let ngsread_file_obj: File = File::open(&ngsread_file).expect("Error during opening the file");
    eprintln!("loading {:?}", &ngsread_file);
    let mut reader: faReader<BufReader<File>> = faReader::new(ngsread_file_obj);
    let mut record: faRecord = faRecord::new();
    let mut sequences: Vec<DnaSequence> = Vec::new();
    eprintln!("loading {:?} done", &ngsread_file);
    'each_read: loop {
        reader.read(&mut record).unwrap();
        if record.is_empty() {
            break 'each_read;
        }
        let sequence_as_vec: Vec<u8> = record.seq().to_vec();
        let current_sequence: DnaSequence = DnaSequence::new(&sequence_as_vec);
        sequences.push(current_sequence);
    }
    let chunk_size: usize = sequences.len() / (threads - 1);
    let sequences_ref: Arc<Vec<DnaSequence>> = Arc::new(sequences);
    let primer_tuple_ref: Arc<Vec<(Vec<u8>, DnaSequence, DnaSequence)>> = Arc::new(primer_tuple);
    //let mut cbf_oyadama: Vec<u32> = vec![0;BLOOMFILTER_TABLE_SIZE];

    //let mut product_size_hashmap = HashMap::<u32, usize>::new();
    thread::scope(|scope| {
        let mut children_1: Vec<thread::ScopedJoinHandle<'_, Vec<u8>>> = Vec::new();
        for i in 1..threads {
            let sequences_ref: Arc<Vec<DnaSequence>> = Arc::clone(&sequences_ref);
            let primer_tuple_ref: Arc<Vec<(Vec<u8>, DnaSequence, DnaSequence)>> =
                Arc::clone(&primer_tuple_ref);
            children_1.push(scope.spawn(move || {
                let start_idx: usize = (i - 1) * chunk_size;
                let end_idx: usize;
                if i != threads - 1 {
                    end_idx = i * chunk_size;
                } else {
                    end_idx = sequences_ref.len() - 1;
                }
                let primer_tuple_ref_mine: Vec<(Vec<u8>, DnaSequence, DnaSequence)> =
                    (*primer_tuple_ref).clone();
                let slice_sequences: Vec<DnaSequence> =
                    Vec::from(sequences_ref[start_idx..end_idx].to_vec());

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
        eprintln!("start  writing to output file: {:?}", &output_file);
        let mut file: File = File::create(&output_file).unwrap();
        for child in children_1 {
            let result: Result<Vec<u8>, Box<dyn Any + Send>> = child.join(); // Resultを返さず、直接戻り値を取得
            file.write_all(&result.unwrap()).unwrap();
        }
        eprintln!("finish writing to output file: {:?}", &output_file);
    });
    Ok(())
}
