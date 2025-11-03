use getopts::Options;
use std::env;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;

fn main() {
    let args: Vec<String> = env::args().collect();
    let program = args[0].clone();

    let mut opts = Options::new();
    opts.optopt("x", "output1", "set output file name (odd numbers)", "NAME");
    opts.optopt(
        "y",
        "output2",
        "set output file name (even numbers)",
        "NAME",
    );
    opts.optflag("h", "help", "print this help menu");

    let matches = match opts.parse(&args[1..]) {
        Ok(m) => m,
        Err(f) => panic!("{}", f.to_string()),
    };
    if matches.opt_present("h") {
        print_usage(&program, &opts);
        return;
    }

    let output_file_odd = matches.opt_str("output1").unwrap_or("odd.bin".to_string());
    let output_file_even = matches.opt_str("output2").unwrap_or("even.bin".to_string());

    let mut writer_odd = BufWriter::new(File::create(output_file_odd).unwrap());
    let mut writer_even = BufWriter::new(File::create(output_file_even).unwrap());

    for i in 0..100 {
        let number = i as u128;
        let buf = number.to_be_bytes(); // Convert to big-endian byte array

        if i % 2 == 0 {
            // Write even numbers (2x) to even file
            writer_even.write_all(&buf).unwrap();
        } else {
            // Write odd numbers (2x+1) to odd file
            writer_odd.write_all(&buf).unwrap();
        }
    }

    writer_odd.flush().unwrap();
    writer_even.flush().unwrap();
}

fn print_usage(program: &str, opts: &Options) {
    let brief = format!("Usage: {} FILE [options]", program);
    print!("{}", opts.usage(&brief));
}
