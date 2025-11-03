use binary_heap_plus::BinaryHeap;
use std::env;
use std::error::Error;
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};

const U128_SIZE: usize = 16; // u128 is 16 bytes

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();
    let mut input_files: Vec<String> = Vec::new();
    let mut output_file: String = String::new();

    // Parse command line arguments
    let mut args_iter: std::iter::Skip<std::slice::Iter<'_, String>> = args.iter().skip(1);
    while let Some(arg) = args_iter.next() {
        match arg.as_str() {
            "-o" => output_file = args_iter.next().unwrap().clone(),
            "-i" => input_files.push(args_iter.next().unwrap().clone()),
            _ => {}
        }
    }

    if input_files.is_empty() || output_file.is_empty() {
        eprintln!("Usage: merge_u128 -i <input_files> -o <output_file>");
        return Ok(());
    }

    // Open all input files and create BufReaders
    let mut readers: Vec<_> = input_files
        .iter()
        .map(|file| BufReader::new(File::open(file).expect("Unable to open input file")))
        .collect();
    eprintln!("file count: {}", readers.len());
    // Create a BufWriter for the output file
    let mut writer: BufWriter<File> = BufWriter::new(File::create(&output_file)?);

    // Create a binary heap for merging
    let mut heap = BinaryHeap::new_by(|a: &(u128, usize), b: &(u128, usize)| {
        a.0.cmp(&b.0).reverse() // We want a min-heap
    });

    // Initialize the heap with the first element from each reader
    for (index, reader) in readers.iter_mut().enumerate() {
        let mut buffer: [u8; 16] = [0u8; U128_SIZE];
        if reader.read_exact(&mut buffer).is_ok() {
            let num: u128 = u128::from_be_bytes(buffer);
            // eprintln!("heap.push1: {:X}, {}", num, index);
            heap.push((num, index));
        }
    }

    // Track the last number written to avoid duplicates
    let mut last_written = None;
    // Track the number of unique u128 integers written to the output file
    // let mut output_count = 0;

    // Iterate over the heap and write unique values to the output file
    while let Some((number, index)) = heap.pop() {
        // eprintln!("loop  head: {:X}, {}", number, index);
        if last_written != Some(number) {
            writer.write_all(&number.to_be_bytes())?;
            last_written = Some(number);
            // output_count += 1; // Increment the count for each unique number written
        }
        // Read the next number from the reader that provided the last number
        let mut buffer = [0u8; U128_SIZE];
        if readers[index].read_exact(&mut buffer).is_ok() {
            // eprintln!("heap.push2: {:X}, {}", u128::from_be_bytes(buffer), index);
            heap.push((u128::from_be_bytes(buffer), index));
        }
    }

    // Flush the writer to ensure all data is written to the output file
    writer.flush()?;

    /*

    // Calculate the total number of u128 integers read from input files
    let total_input_count: usize = input_files
        .iter()
        .map(|file| count_u128_integers(file).unwrap_or(0))
        .sum();

    // Output results to stderr in markdown format
    eprintln!("| Input File | Number of u128 integers |");
    eprintln!("|------------|------------------------:|");
    for file in &input_files {
        eprintln!("| {} | {} |", file, count_u128_integers(file)?);
    }

    eprintln!("\n| Output File |");
    eprintln!("|-------------|");
    eprintln!("| {} |", output_file);

    eprintln!("\n| Sum of u128 integers |");
    eprintln!("|-----------------------:|");
    eprintln!("| {} |", total_input_count);

    eprintln!("\n| Cardinality of u128 integers |");
    eprintln!("|-------------------------------:|");
    eprintln!("| {} |", output_count);

    // Calculate and print the number of duplicates removed
    let duplicates_removed = total_input_count - output_count;
    eprintln!("\n| Duplicates removed |");
    eprintln!("|-------------------:|");
    eprintln!("| {} |", duplicates_removed);

    eprintln!(
        "{}\t{}\t{}",
        total_input_count, output_count, duplicates_removed
    );


    */

    Ok(())
}

/*
// Count the number of u128 integers in a file
fn count_u128_integers(file_name: &str) -> io::Result<usize> {
    let file: File = File::open(file_name)?;
    let metadata: std::fs::Metadata = file.metadata()?;
    Ok(metadata.len() as usize / U128_SIZE)
}
*/
