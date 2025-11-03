extern crate search_primer;
use getopts::Options;
use search_primer::sequence_encoder_util::{decode_u128_l, decode_u128_r};
use std::collections::HashSet;
use std::env;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};

fn main() {
    let args: Vec<String> = env::args().collect();

    let mut opts: Options = Options::new();
    opts.optflag("h", "help", "print this help menu");
    opts.optmulti("i", "inputfiles", "set input file names", "NAME");
    opts.optopt("o", "outputfile", "set output file name", "NAME");

    let matches: getopts::Matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => {
            panic!("{}", f.to_string())
        }
    };

    if matches.opt_present("h") {
        println!("{}", opts.usage(&format!("Usage: {} [options]", args[0])));
        return;
    }

    let mut u128_primer_candidates_set: HashSet<u128> = HashSet::with_capacity(10000);
    let files: Vec<String> = matches.opt_strs("i");
    for file in files {
        let f: File = File::open(&file).unwrap();
        let mut reader: BufReader<File> = BufReader::new(f);
        let mut buffer: [u8; 16] = [0u8; 16];

        loop {
            let result: Result<(), std::io::Error> =
                reader.by_ref().take(16).read_exact(&mut buffer);
            match result {
                Ok(_val) => {}
                Err(_err) => break,
            }
            let tmp_val: u128 = u128::from_be_bytes(buffer);
            u128_primer_candidates_set.insert(tmp_val);
        }
    }
    let mut u128_primer_candidate_vec: Vec<u128> = u128_primer_candidates_set.into_iter().collect();
    u128_primer_candidate_vec.sort();

    let output_file = if matches.opt_present("o") {
        matches.opt_str("o").unwrap()
    } else {
        "u128_binary_merge_out.fa".to_string()
    };

    let mut w: BufWriter<File> = BufWriter::new(File::create(&output_file).unwrap());

    for each_primer_candidate in &u128_primer_candidate_vec {
        // let full_seq = decode_u128_2_dna_seq(&each_primer_candidate, 64);
        let l_seq = decode_u128_l(&each_primer_candidate);
        let r_seq = decode_u128_r(&each_primer_candidate);
        let hex_candidate = format!("{:X}", each_primer_candidate);

        w.write(b">").expect("Failed to write to file");
        w.write(hex_candidate.as_bytes())
            .expect("Failed to write to file");
        w.write(b"_l\n").expect("Failed to write to file");
        w.write(&l_seq).expect("Failed to write to file");
        w.write(b"\n\n").expect("Failed to write to file");

        // Write the right sequence in FASTA format
        w.write(b">").expect("Failed to write to file");
        w.write(hex_candidate.as_bytes())
            .expect("Failed to write to file");
        w.write(b"_r\n").expect("Failed to write to file");
        w.write(&r_seq).expect("Failed to write to file");
        w.write(b"\n\n").expect("Failed to write to file");
        // let dna_seq = decode_u128_2_dna_seq(&each_primer_candidate, 64);
        // w.write(&dna_seq).unwrap(); // Vec<u8>から&[u8]への参照を渡す
        // w.write(b"\n").unwrap();
    }
}
