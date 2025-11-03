use getopts::Options;
use std::collections::HashSet;
use std::env;
use std::error::Error;
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};

const U128_SIZE: usize = 16; // u128 is 16 bytes

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();

    let mut opts = Options::new();
    opts.optflag("h", "help", "print this help menu");
    opts.optmulti("i", "inputfiles", "set input file names", "NAME");
    opts.optopt("o", "outputfile", "set output file name", "NAME");

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => {
            panic!("{}", f.to_string())
        }
    };

    if matches.opt_present("h") {
        println!("{}", opts.usage(&format!("Usage: {} [options]", args[0])));
        return Ok(());
    }

    let output_file = if matches.opt_present("o") {
        matches.opt_str("o").unwrap()
    } else {
        "u128_binary_merge_out.bin".to_string()
    };

    // Create a binary heap for merging

    let mut u128_counter: HashSet<u128> = HashSet::new();
    let files = matches.opt_strs("i");

    // Initialize the heap with the first element from each reader
    let mut u128_cnt: Vec<usize> = Vec::new();
    for file in &files {
        let f: File = File::open(&file).unwrap();
        let mut reader: BufReader<File> = BufReader::new(f);
        let mut buffer: [u8; 16] = [0u8; 16];
        let mut cnt_in_this_file: usize = 0;
        loop {
            let result = reader.by_ref().take(16).read_exact(&mut buffer);
            match result {
                Ok(_val) => {}
                Err(_err) => break,
            }
            let tmp_val: u128 = u128::from_be_bytes(buffer);
            u128_counter.insert(tmp_val);
            cnt_in_this_file += 1;
        }
        u128_cnt.push(cnt_in_this_file);
    }
    let output_count: usize = u128_counter.len();
    let total_input_cnt: usize = u128_cnt.iter().sum::<usize>();
    // Flush the writer to ensure all data is written to the output file

    // Calculate the total number of u128 integers read from input files

    // Output results to stderr in markdown format
    eprintln!("| Input File | Number of u128 integers |");
    eprintln!("|------------|------------------------:|");
    for file in &files {
        eprintln!("| {} | {} |", file, count_u128_integers(file)?);
    }

    eprintln!("\n| Output File |");
    eprintln!("|-------------|");
    eprintln!("| {} |", output_file);

    eprintln!("\n| Sum of u128 integers |");
    eprintln!("|-----------------------:|");
    eprintln!("| {} |", total_input_cnt);

    eprintln!("\n| Cardinality of u128 integers |");
    eprintln!("|-------------------------------:|");
    eprintln!("| {} |", output_count);

    // Calculate and print the number of duplicates removed
    let duplicates_removed = total_input_cnt - output_count;
    eprintln!("\n| Duplicates removed |");
    eprintln!("|-------------------:|");
    eprintln!("| {} |", duplicates_removed);

    eprintln!(
        "{}\t{}\t{}",
        total_input_cnt, output_count, duplicates_removed
    );

    let mut u128_counter_vec: Vec<u128> = u128_counter.into_iter().collect();
    u128_counter_vec.sort();
    let mut w: BufWriter<File> = BufWriter::new(File::create(&output_file).unwrap());
    let mut buf_array: [u8; 16] = [0; 16];
    let mut buf_num: u128;
    for each_u128_counter_vec in &u128_counter_vec {
        buf_num = *each_u128_counter_vec;
        for i in 0..16 {
            buf_array[15 - i] = u8::try_from(buf_num & 0xFF).unwrap();
            buf_num >>= 8;
        }
        w.write_all(&buf_array).unwrap();
    }

    Ok(())
}

/// Count the number of u128 integers in a file
fn count_u128_integers(file_name: &str) -> io::Result<usize> {
    let file: File = File::open(file_name)?;
    let metadata: std::fs::Metadata = file.metadata()?;
    Ok(metadata.len() as usize / U128_SIZE)
}
