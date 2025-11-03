use std::env;
use std::error::Error;
use std::fs::File;
use std::io::{Read, Write};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();
    let mut input_file: String = String::new();
    let mut output_file_base: String = String::new();
    let mut split_count: usize = 20; // Default split count

    // Parse command line arguments
    let mut args_iter = args.iter().skip(1);
    while let Some(arg) = args_iter.next() {
        match arg.as_str() {
            "-o" => output_file_base = args_iter.next().unwrap().clone(),
            "-i" => input_file = args_iter.next().unwrap().clone(),
            "-n" => split_count = args_iter.next().unwrap().parse()?,
            _ => {}
        }
    }

    if input_file.is_empty() || output_file_base.is_empty() {
        eprintln!("Usage: merge_u128 -i <input_file> -o <output_file_base> [-n <split_count>]");
        return Ok(());
    }
    eprintln!("{}\n{}", input_file, output_file_base);
    // Process each input file
    let mut file = File::open(input_file)?;
    let mut buffer = [0u8; 16]; // u128 is 16 bytes
    let mut outputs = vec![Vec::new(); split_count];
    let mut index = 0;

    // Read and split the data
    while let Ok(_) = file.read_exact(&mut buffer) {
        // let number = u128::from_be_bytes(buffer);
        outputs[index].extend_from_slice(&buffer);
        index = (index + 1) % split_count;
    }

    // Write to output files
    for (i, data) in outputs.into_iter().enumerate() {
        let mut output_file = File::create(format!(
            "{}_{}.bin",
            output_file_base,
            format!("{:03}", i + 1)
        ))?;
        output_file.write_all(&data)?;
    }

    Ok(())
}
