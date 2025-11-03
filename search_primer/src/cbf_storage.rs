use crate::config::Config;
use crate::counting_bloomfilter_util::{
    build_counting_bloom_filter, count_lr_tuple_with_hashtable, number_of_high_occurence_lr_tuple,
};
use crate::sequence_encoder_util::DnaSequence;
use bio::io::fasta::{FastaRead, Reader as faReader, Record as faRecord};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufReader, BufWriter, Read, Write};
use std::sync::{Arc, Mutex};
use std::thread;

const CBF_MAGIC: &[u8; 8] = b"CBFV0001";

#[derive(Debug)]
pub struct CbfHeader {
    pub magic: [u8; 8],
    pub config_hash: u64,
    pub table_size: u64,
    pub reserved: u64,
}

impl CbfHeader {
    fn new(config: &Config) -> Self {
        Self {
            magic: *CBF_MAGIC,
            config_hash: Self::hash_config(config),
            table_size: config.bloomfilter_table_size as u64,
            reserved: 0,
        }
    }

    fn hash_config(config: &Config) -> u64 {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};
        
        let mut hasher = DefaultHasher::new();
        config.l_len.hash(&mut hasher);
        config.r_len.hash(&mut hasher);
        config.chunk_max.hash(&mut hasher);
        config.margin_size.hash(&mut hasher);
        hasher.finish()
    }

    fn write_to<W: Write>(&self, writer: &mut W) -> std::io::Result<()> {
        writer.write_all(&self.magic)?;
        writer.write_all(&self.config_hash.to_le_bytes())?;
        writer.write_all(&self.table_size.to_le_bytes())?;
        writer.write_all(&self.reserved.to_le_bytes())?;
        Ok(())
    }

    fn read_from<R: Read>(reader: &mut R) -> std::io::Result<Self> {
        let mut magic = [0u8; 8];
        let mut config_hash_bytes = [0u8; 8];
        let mut table_size_bytes = [0u8; 8];
        let mut reserved_bytes = [0u8; 8];

        reader.read_exact(&mut magic)?;
        reader.read_exact(&mut config_hash_bytes)?;
        reader.read_exact(&mut table_size_bytes)?;
        reader.read_exact(&mut reserved_bytes)?;

        Ok(Self {
            magic,
            config_hash: u64::from_le_bytes(config_hash_bytes),
            table_size: u64::from_le_bytes(table_size_bytes),
            reserved: u64::from_le_bytes(reserved_bytes),
        })
    }

    fn is_valid(&self) -> bool {
        self.magic == *CBF_MAGIC
    }

    fn is_compatible(&self, config: &Config) -> bool {
        self.config_hash == Self::hash_config(config) && 
        self.table_size == config.bloomfilter_table_size as u64
    }
}

fn load_sequences(fasta_path: &str) -> Result<Vec<DnaSequence>, Box<dyn std::error::Error>> {
    let file = File::open(fasta_path)?;
    let mut reader = faReader::new(file);
    let mut record = faRecord::new();
    let mut sequences = Vec::new();

    loop {
        reader.read(&mut record)?;
        if record.is_empty() {
            break;
        }
        let sequence_as_vec: Vec<u8> = record.seq().to_vec();
        let current_sequence = DnaSequence::new(&sequence_as_vec);
        sequences.push(current_sequence);
    }
    Ok(sequences)
}

pub fn build_and_save_cbf(
    fasta_path: &str,
    cbf_output_path: &str,
    config: &Config,
    threads: usize,
) -> Result<(), Box<dyn std::error::Error>> {
    eprintln!("Loading sequences from {}", fasta_path);
    let sequences = Arc::new(load_sequences(fasta_path)?);
    
    eprintln!("Building CBF with {} threads", threads);
    let chunk_size = sequences.len() / threads.max(1);
    let mut cbf_result = vec![0u16; config.bloomfilter_table_size];

    thread::scope(|scope| {
        let mut children = Vec::new();
        
        for i in 0..threads {
            let sequences = Arc::clone(&sequences);
            children.push(scope.spawn(move || {
                let start_idx = i * chunk_size;
                let end_idx = if i == threads - 1 {
                    sequences.len()
                } else {
                    (i + 1) * chunk_size
                };
                
                if start_idx >= sequences.len() {
                    return vec![0u16; config.bloomfilter_table_size];
                }
                
                build_counting_bloom_filter(&sequences, start_idx, end_idx, config, i + 1)
            }));
        }

        for child in children {
            let cbf = child.join().unwrap();
            for (i, &value) in cbf.iter().enumerate() {
                cbf_result[i] = cbf_result[i].saturating_add(value);
            }
        }
    });

    eprintln!("Saving CBF to {}", cbf_output_path);
    let file = File::create(cbf_output_path)?;
    let mut writer = BufWriter::new(file);
    
    let header = CbfHeader::new(config);
    header.write_to(&mut writer)?;
    
    for &value in &cbf_result {
        writer.write_all(&value.to_le_bytes())?;
    }
    writer.flush()?;
    
    eprintln!("CBF saved successfully");
    Ok(())
}

pub fn load_cbf_and_extract_candidates(
    cbf_path: &str,
    fasta_path: &str,
    threshold: u16,
    config: &Config,
    threads: usize,
) -> Result<HashSet<u128>, Box<dyn std::error::Error>> {
    eprintln!("Loading CBF from {}", cbf_path);
    let file = File::open(cbf_path)?;
    let mut reader = BufReader::new(file);
    
    let header = CbfHeader::read_from(&mut reader)?;
    if !header.is_valid() {
        return Err("Invalid CBF file format".into());
    }
    if !header.is_compatible(config) {
        return Err("CBF file incompatible with current config".into());
    }

    let mut cbf_table = vec![0u16; header.table_size as usize];
    for value in &mut cbf_table {
        let mut bytes = [0u8; 2];
        reader.read_exact(&mut bytes)?;
        *value = u16::from_le_bytes(bytes);
    }

    eprintln!("Loading sequences from {}", fasta_path);
    let sequences = Arc::new(load_sequences(fasta_path)?);
    let cbf_table = Arc::new(cbf_table);
    
    eprintln!("Extracting candidates with threshold {}", threshold);
    let chunk_size = sequences.len() / threads.max(1);
    let candidates = Arc::new(Mutex::new(HashSet::with_capacity(config.hashset_size)));

    thread::scope(|scope| {
        let mut children = Vec::new();
        
        for i in 0..threads {
            let sequences = Arc::clone(&sequences);
            let cbf_table = Arc::clone(&cbf_table);
            let candidates = Arc::clone(&candidates);
            
            children.push(scope.spawn(move || {
                let start_idx = i * chunk_size;
                let end_idx = if i == threads - 1 {
                    sequences.len()
                } else {
                    (i + 1) * chunk_size
                };
                
                if start_idx >= sequences.len() {
                    return;
                }
                
                let local_candidates = number_of_high_occurence_lr_tuple(
                    &cbf_table,
                    &sequences,
                    start_idx,
                    end_idx,
                    threshold,
                    config,
                    i + 1,
                );
                
                candidates.lock().unwrap().extend(local_candidates);
            }));
        }

        for child in children {
            child.join().unwrap();
        }
    });

    let result = Arc::try_unwrap(candidates).unwrap().into_inner().unwrap();
    eprintln!("Found {} candidates", result.len());
    Ok(result)
}

pub fn remove_false_positives(
    candidates: &HashSet<u128>,
    fasta_path: &str,
    threshold: u16,
    config: &Config,
    threads: usize,
) -> Result<HashMap<u128, u16>, Box<dyn std::error::Error>> {
    eprintln!("Loading sequences from {}", fasta_path);
    let sequences = Arc::new(load_sequences(fasta_path)?);
    
    eprintln!("Removing false positives from {} candidates", candidates.len());
    let chunk_size = sequences.len() / threads.max(1);
    let result = Arc::new(Mutex::new(HashMap::with_capacity(candidates.len())));

    thread::scope(|scope| {
        let mut children = Vec::new();
        
        for i in 0..threads {
            let sequences = Arc::clone(&sequences);
            let result = Arc::clone(&result);
            
            children.push(scope.spawn(move || {
                let start_idx = i * chunk_size;
                let end_idx = if i == threads - 1 {
                    sequences.len()
                } else {
                    (i + 1) * chunk_size
                };
                
                if start_idx >= sequences.len() {
                    return;
                }
                
                let local_counts = count_lr_tuple_with_hashtable(
                    &sequences,
                    start_idx,
                    end_idx,
                    candidates,
                    config,
                    i + 1,
                );
                
                let mut global_result = result.lock().unwrap();
                for (key, value) in local_counts {
                    *global_result.entry(key).or_insert(0) += value;
                }
            }));
        }

        for child in children {
            child.join().unwrap();
        }
    });

    let mut final_result = Arc::try_unwrap(result).unwrap().into_inner().unwrap();
    final_result.retain(|_, &mut count| count >= threshold);
    
    eprintln!("Final result: {} LR-tuples above threshold", final_result.len());
    Ok(final_result)
}
