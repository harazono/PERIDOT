extern crate search_primer;
extern crate getopts;
use std::{env, process};
use std::io::{Write, BufWriter};
use std::io::{Read,  BufReader};
use std::fs;
use std::fs::File;
use getopts::Options;
use search_primer::sequence_encoder_util::{decode_u128_l, decode_u128_r};


fn print_usage(program: &str, opts: &Options) {
    let brief = format!("Usage: {} FILE [options]", program);
    print!("{}", opts.usage(&brief));
    process::exit(0);
}

fn blast_formatter(sequence: &u128) -> String{
    let l_u8_array = decode_u128_l(sequence);
    let r_u8_array = decode_u128_r(sequence);
    let l_str: &str = std::str::from_utf8(&l_u8_array).unwrap();
    let r_str: &str = std::str::from_utf8(&r_u8_array).unwrap();
    let fasta_fmt = format!(">{:0x}-L\n{}\n>{:0x}-R\n{}\n", sequence, l_str,  sequence, r_str);
    return fasta_fmt;
}

fn main(){
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();

    let mut opts = Options::new();
    opts.optopt("a", "output1", "sequence.fasta", "FILENAME");
    opts.optopt("b", "output2", "namelist.txt", "FILENAME");

    opts.optflag("h", "help", "print this help menu");

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => { m }
        Err(f) => { panic!("{}", f.to_string()) }
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

    let output_file_1 = if matches.opt_present("a") {
        matches.opt_str("a").unwrap()
    }else{
        "out_1.fasta".to_string()
    };
    let output_file_2 = if matches.opt_present("b") {
        matches.opt_str("b").unwrap()
    }else{
        "out_2.txt".to_string()
    };

    let mut w1 = BufWriter::new(fs::File::create(&output_file_1).unwrap());
    let mut w2 = BufWriter::new(fs::File::create(&output_file_2).unwrap());


    let f: File = File::open(&input_file).unwrap();
    let mut reader = BufReader::new(f);
    let mut buf: [u8; 16] = [0; 16];
    let mut tmp_seq_as_u128: u128 = 0;

    loop {
        match reader.read(&mut buf).unwrap() {
            0 => break,
            n => {
                let buf = &buf[..n];
                for i in 0..16 {
                    tmp_seq_as_u128 <<= 8;
                    tmp_seq_as_u128 += u128::from(buf[i]);
                }
                writeln!(&mut w1, "{}", blast_formatter(&tmp_seq_as_u128)).unwrap();
                writeln!(&mut w2, "{:0x}", tmp_seq_as_u128).unwrap();
            }
        }
    }
}