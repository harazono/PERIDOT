extern crate getopts;
extern crate search_primer;
use std::fs::File;
use std::io::{BufReader, Read};
use std::{env, process};
//use std::io::prelude::*;
use getopts::Options;
use search_primer::sequence_encoder_util::{decode_u128_l, decode_u128_r};
use std::fs::OpenOptions;
use std::io::Write;
use std::mem;
use std::process::{Command, Output, Stdio};
use std::sync::{Arc, Mutex};
use std::thread;
use std::time::Instant;

extern crate rand;

use rand::distributions::Alphanumeric;
use rand::Rng;

/*

"SEQUENCE_ID={:0x}
SEQUENCE_TEMPLATE={}
PRIMER_TASK=pick_pcr_primers
PRIMER_OPT_SIZE=30
PRIMER_MIN_SIZE=16
PRIMER_MAX_SIZE=32
PRIMER_PRODUCT_SIZE_RANGE=101-200 201-301
P3_FILE_FLAG=0
PRIMER_EXPLAIN_FLAG=1
PRIMER_OPT_TM=66.0
PRIMER_MAX_TM=72.0
PRIMER_MAX_LIBRARY_MISPRIMING=11",

 */

fn primer3_core_input_sequences(sequences: &Vec<&u128>, primer3_config: &String) -> String {
    let mut ret_str: String = String::new();
    let many_n = "N".to_string().repeat(50);

    /*     eprintln!(
            "primer3_core_input_sequence: sequense length...{}",
            sequences.len()
        );
    */
    for each_seq in sequences {
        let l_u8_array: Vec<u8> = decode_u128_l(each_seq);
        let r_u8_array: Vec<u8> = decode_u128_r(each_seq);
        let l_str: &str = std::str::from_utf8(&l_u8_array).unwrap();
        let r_str: &str = std::str::from_utf8(&r_u8_array).unwrap();
        let sequence_with_internal_n: String = format!("{}{}{}", l_str, many_n, r_str);
        let primer3_fmt_str: String = format!(
            "SEQUENCE_ID={:0x}\nSEQUENCE_TEMPLATE={}\n{}=\n",
            each_seq, sequence_with_internal_n, primer3_config
        );
        // Check if library_file_name is Some or None

        /*         match library_file_name.as_ref() {
                    Some(file_name) => {
                       // If Some, append the file name to the string
                        primer3_fmt_str
                            .push_str(&format!("\nPRIMER_MISPRIMING_LIBRARY={}\n=\n", file_name));
                    }
                    None => {
                        primer3_fmt_str.push_str("\n=\n");
                    }
                }
        */
        ret_str.push_str(&primer3_fmt_str);
    }
    // eprintln!("mem::size_of_val(&ret_str): {}", mem::size_of_val(&ret_str));//これだとVecのアタマの24byteしか表示されない
    // eprintln!("(&ret_str).len(): {}", (&ret_str).len()); //function execute_primer3: formatted_string.len(): 483000
    return ret_str;
}

fn generate_random_string(length: usize) -> String {
    rand::thread_rng()
        .sample_iter(&Alphanumeric)
        .take(length)
        .map(char::from)
        .collect()
}

fn execute_primer3(
    formatted_string: String,
    temp_file_name_prefix: &String,
    thread_id: usize,
) -> String {
    /*
    /tmpなどに中間ファイルを書き込む
    let mut file = OpenOptions::new().write(true).create(true).open("/tmp/primer3_core_input.txt").unwrap();
    ファイルをクローズして、再び開ける
    Command::newにファイルオブジェクトを渡す（stdinは使わない）
    */
    // /tmpに中間ファイルを書き込む

    let timer: Instant = Instant::now();
    let temporary_file_name = format!(
        "{}_{}.txt",
        temp_file_name_prefix,
        generate_random_string(10)
    );
    eprintln!("temporary_file_name: {}", temporary_file_name);

    {
        let mut file = OpenOptions::new()
            .write(true)
            .create(true)
            .open(&temporary_file_name)
            .unwrap();
        file.write_all(formatted_string.as_bytes()).unwrap();
        file.flush().unwrap();
    } // ファイルをクローズする

    // ファイルを再度開き、Command::new("primer3_core")の標準入力として利用する
    eprintln!("start opening {}", temporary_file_name);
    let output: Output = Command::new("primer3_core")
        .stdin(Stdio::from(
            OpenOptions::new()
                .read(true)
                .open(&temporary_file_name)
                .unwrap(),
        ))
        .output()
        .expect("Failed to execute command");

    // ファイルを削除する
    eprintln!("start removing {}", temporary_file_name);
    let _ = std::fs::remove_file(temporary_file_name);

    // 結果を返す

    let end_time: std::time::Duration = timer.elapsed();
    eprintln!(
        "primer3_caller[{:02}] takes {}.{} sec",
        thread_id,
        end_time.as_secs(),
        end_time.subsec_millis()
    );

    String::from_utf8(output.stdout).unwrap()
}

fn print_usage(program: &str, opts: &Options) {
    let brief = format!("Usage: {} FILE ", program);
    print!("{}", opts.usage(&brief));
    process::exit(0);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();

    let mut opts = Options::new();
    opts.optflag("h", "help", "print this help menu");
    opts.optopt(
        "t",
        "thread",
        "number of thread to use for radix sort. default value is 8.",
        "THREAD",
    );
    opts.optopt("c", "config", "config file for primer3_core.", "CONFIG");
    opts.optopt(
        "o",
        "output",
        "output file name for primer3 results",
        "OUTPUT",
    ); // New option for output file
    opts.optopt("m", "tmpfile", "set temporary file name prefix", "TEMP");

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

    let thread_number: usize = if matches.opt_present("t") {
        matches.opt_str("t").unwrap().parse::<usize>().unwrap()
    } else {
        4
    };
    let config_file_name: String = if matches.opt_present("c") {
        matches.opt_str("c").unwrap()
    } else {
        print_usage(&program, &opts);
        return;
    };

    let temporary_file_name_prefix = if matches.opt_present("m") {
        matches.opt_str("m").unwrap()
    } else {
        "/work/primer3_tmp".to_string()
    };

    // let library_file_name: Option<String> = matches.opt_str("l");

    let input_file = if !matches.free.is_empty() {
        matches.free[0].clone()
    } else {
        print_usage(&program, &opts);
        return;
    };

    eprintln!("start  loading {:?}", &config_file_name);
    let mut f_c: File = File::open(&config_file_name).unwrap();
    let mut primer3_config = String::new();
    f_c.read_to_string(&mut primer3_config).unwrap();
    eprintln!("finish loading {:?}", &config_file_name);
    eprintln!("File primer3_config:\n{}", primer3_config);

    eprintln!("start  loading {:?}", &input_file);
    let f: File = File::open(&input_file).unwrap();
    let mut reader = BufReader::new(f);
    let mut buffer: [u8; 16] = [0; 16];
    let mut tmp_seq_as_u128: u128;
    let mut candidates: Vec<u128> = Vec::new();
    loop {
        let result = reader.by_ref().take(16).read_exact(&mut buffer);
        match result {
            Ok(_val) => {}
            Err(_err) => break,
        }
        tmp_seq_as_u128 = u128::from_be_bytes(buffer);
        //println!("{:?}", String::from_utf8(decode_u128_2_dna_seq(&tmp_seq_as_u128, 64)).unwrap());
        candidates.push(tmp_seq_as_u128);
    }
    eprintln!("finish loading {:?}", &input_file);

    //eprintln!("start formatting string");
    //let primer3_fmt_string: Vec<String> = primer3_core_input_sequence(&candidates, &library_file_name);
    //let bunch_of_50000_fmt_string: Vec<Vec<String>> = primer3_fmt_string.chunks(1).map(|chunk| chunk.to_vec()).collect();
    //eprintln!("finish formatting string");
    let mut chunks_of_input: Vec<Vec<u128>> = Vec::new();

    for _i in 0..thread_number {
        chunks_of_input.push(Vec::new()); //スレッド分のVecを用意する
    }
    for (index, bunch) in candidates.iter().enumerate() {
        chunks_of_input[index % thread_number].push(*bunch);
    }

    let arc_chunks_of_input: Arc<Vec<Vec<u128>>> = Arc::new(chunks_of_input);
    let final_result: Arc<Mutex<Vec<String>>> = Arc::new(Mutex::new(Vec::new()));
    let mut children: Vec<thread::JoinHandle<()>> = Vec::new();
    let output_file_name: String = matches
        .opt_str("o")
        .unwrap_or("default_output.txt".to_string()); // Get output file name
    let file_mutex: Arc<Mutex<File>> = Arc::new(Mutex::new(
        OpenOptions::new()
            .append(true)
            .create(true)
            .open(&output_file_name)
            .expect("Unable to open file"),
    ));

    for i in 0..thread_number {
        let chunks_of_input: Arc<Vec<Vec<u128>>> = Arc::clone(&arc_chunks_of_input);
        let arc_final_result: Arc<Mutex<Vec<String>>> = Arc::clone(&final_result);
        let thread_file_mutex: Arc<Mutex<File>> = Arc::clone(&file_mutex);
        // let library_file_name_clone: Option<String> = library_file_name.clone(); // Clone the library_file_name
        let temporary_file_name_prefix_clone: String = temporary_file_name_prefix.clone(); // Clone the temporary_file_name_prefix
        let primer3_config_clone = primer3_config.clone();
        children.push(thread::spawn(move || {
            let mut primer3_results: String = String::new();
            let total_elements: usize = chunks_of_input[i].len();
            let mut processed_elements: usize = 0;
            // chunks_of_input[i]を500個の要素ごとのチャンクに分割
            for bunch in chunks_of_input[i].chunks(100) {
                let start: Instant = Instant::now();
                let sequences: Vec<_> = bunch.iter().collect();
                processed_elements += sequences.len();
                primer3_results += &execute_primer3(
                    primer3_core_input_sequences(
                        &sequences,
                        // &library_file_name_clone,
                        &primer3_config_clone,
                    ),
                    &temporary_file_name_prefix_clone,
                    i,
                );
                //eprintln!("{}", &primer3_results);
                if mem::size_of_val(&primer3_results) > 2 * 1024 * 1024 * 1024 {
                    let mut file: std::sync::MutexGuard<'_, File> =
                        thread_file_mutex.lock().unwrap(); // Use the cloned mutex
                    file.write_all(primer3_results.as_bytes())
                        .expect("Unable to write to file");
                    primer3_results.clear();
                }
                let elapsed: std::time::Duration = start.elapsed();

                let percentage_done = (processed_elements as f64 / total_elements as f64) * 100.0;
                eprintln!(
                    "Thread[{:02}] - Processing chunk took {}.{} sec, {:.2}% done",
                    i,
                    elapsed.as_secs(),
                    elapsed.subsec_millis(),
                    percentage_done
                );
            }
            // After the loop, write remaining primer3_results to the file
            if !primer3_results.is_empty() {
                let mut file: std::sync::MutexGuard<'_, File> = thread_file_mutex.lock().unwrap();
                file.write_all(primer3_results.as_bytes())
                    .expect("Unable to write to file");
                primer3_results.clear();
            }

            arc_final_result.lock().unwrap().push(primer3_results);
        }));
    }

    for child in children {
        let _ = child.join();
    }
    eprintln!("finish waiting all threads");
    for i in final_result.lock().unwrap().iter() {
        println!("{}", i);
    }
}
