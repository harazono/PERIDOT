use crate::config::Config;
use crate::counting_bloomfilter_util::{
    build_counting_bloom_filter, count_lr_tuple_with_hashtable, number_of_high_occurence_lr_tuple,
    hash_from_u128, count_occurence_from_counting_bloomfilter_table,
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

pub fn save_cbf_to_file(cbf_data: &[u16], config: &Config, file_path: &str) -> Result<(), Box<dyn std::error::Error>> {
    let file = File::create(file_path)?;
    let mut writer = BufWriter::new(file);
    
    let header = CbfHeader::new(config);
    header.write_to(&mut writer)?;
    
    for &value in cbf_data {
        writer.write_all(&value.to_le_bytes())?;
    }
    writer.flush()?;
    
    Ok(())
}

pub fn load_cbf_from_file(file_path: &str, config: &Config) -> Result<Vec<u16>, Box<dyn std::error::Error>> {
    let file = File::open(file_path)?;
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

    Ok(cbf_table)
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

    save_cbf_to_file(&cbf_result, config, cbf_output_path)?;
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
    let cbf_table = load_cbf_from_file(cbf_path, config)?;

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cbf_query_functionality() -> Result<(), Box<dyn std::error::Error>> {
        let config = Config {
            l_len: 8,
            r_len: 8,
            chunk_max: 40,
            margin_size: 4,
            bloomfilter_table_size: 1024,
            hashset_size: 512,
        };

        // テスト用CBFデータを作成（特定のパターン）
        let mut cbf_data = vec![0u16; config.bloomfilter_table_size];
        // いくつかの位置に既知の値を設定
        cbf_data[100] = 500;
        cbf_data[200] = 1000;
        cbf_data[300] = 1500;

        let temp_cbf = NamedTempFile::new()?;
        let cbf_path = temp_cbf.path().to_str().unwrap();
        
        save_cbf_to_file(&cbf_data, &config, cbf_path)?;
        
        // CBFQueryを作成
        let cbf_query = CbfQuery::from_file(cbf_path, &config)?;
        
        // テスト配列を作成
        let test_seq = "ATCGATCGTTTTGCTAGCTA"; // L="ATCGATCG", margin="TTTT", R="GCTAGCTA"
        let dna_seq = DnaSequence::new(&test_seq.as_bytes().to_vec());
        
        // 問い合わせ実行
        let results = cbf_query.query_sequence(&dna_seq);
        
        // 結果検証
        assert!(!results.is_empty(), "Should find at least one LR-tuple");
        
        // 各結果の形式確認
        for (l_start, l_end, r_start, r_end, count) in &results {
            assert_eq!(l_end - l_start, config.l_len, "L-segment length mismatch");
            assert_eq!(r_end - r_start, config.r_len, "R-segment length mismatch");
            assert!(*r_start >= l_end + config.margin_size, "Margin size violation");
            println!("Found LR-tuple: L[{}:{}] R[{}:{}] count={}", l_start, l_end, r_start, r_end, count);
        }

        Ok(())
    }

    #[test]
    fn test_cbf_query_lr_tuple_direct() -> Result<(), Box<dyn std::error::Error>> {
        let config = Config {
            l_len: 4,
            r_len: 4,
            chunk_max: 20,
            margin_size: 2,
            bloomfilter_table_size: 256,
            hashset_size: 128,
        };

        // 既知のLR-tupleでテスト
        let test_seq = "ATCGTTGCTA"; // L="ATCG", margin="TT", R="GCTA"
        let dna_seq = DnaSequence::new(&test_seq.as_bytes().to_vec());
        let lr_tuple = dna_seq.subsequence_as_u128(vec![[0, 4], [6, 10]]);

        // テスト用CBFを作成
        let mut cbf_data = vec![0u16; config.bloomfilter_table_size];
        let hash_indices = hash_from_u128(lr_tuple, config.bloomfilter_table_size);
        
        // 特定の値を設定
        for &idx in &hash_indices {
            cbf_data[idx as usize] = 42;
        }

        let temp_cbf = NamedTempFile::new()?;
        let cbf_path = temp_cbf.path().to_str().unwrap();
        save_cbf_to_file(&cbf_data, &config, cbf_path)?;

        let cbf_query = CbfQuery::from_file(cbf_path, &config)?;
        
        // 直接問い合わせ
        let count = cbf_query.query_lr_tuple(lr_tuple);
        assert_eq!(count, 42, "Direct LR-tuple query failed");

        // 配列問い合わせ
        let results = cbf_query.query_sequence(&dna_seq);
        assert!(!results.is_empty(), "Sequence query should find the LR-tuple");
        
        let found = results.iter().any(|(_, _, _, _, c)| *c == 42);
        assert!(found, "Should find the LR-tuple with count 42");

        Ok(())
    }

    #[test]
    fn test_cbf_query_threshold_filtering() -> Result<(), Box<dyn std::error::Error>> {
        let config = Config {
            l_len: 6,
            r_len: 6,
            chunk_max: 30,
            margin_size: 3,
            bloomfilter_table_size: 512,
            hashset_size: 256,
        };

        // 複数のLR-tupleを含む配列
        let test_seq = "ATCGATTTGCTAGCTTTAGCTAGC"; // 複数のLR-tupleが生成される
        let dna_seq = DnaSequence::new(&test_seq.as_bytes().to_vec());

        // ランダムなCBFデータ（0-2000の範囲）
        let cbf_data: Vec<u16> = (0..config.bloomfilter_table_size)
            .map(|i| (i % 2001) as u16)
            .collect();

        let temp_cbf = NamedTempFile::new()?;
        let cbf_path = temp_cbf.path().to_str().unwrap();
        save_cbf_to_file(&cbf_data, &config, cbf_path)?;

        let cbf_query = CbfQuery::from_file(cbf_path, &config)?;
        
        // 全結果取得
        let all_results = cbf_query.query_sequence(&dna_seq);
        
        // 閾値フィルタリング
        let threshold = 1000;
        let filtered_results = cbf_query.query_above_threshold(&dna_seq, threshold);
        
        // 検証
        assert!(filtered_results.len() <= all_results.len(), "Filtered results should be subset");
        
        for (_, _, _, _, count) in &filtered_results {
            assert!(*count >= threshold, "All filtered results should be above threshold");
        }

        println!("All results: {}, Above threshold ({}): {}", 
                all_results.len(), threshold, filtered_results.len());

        Ok(())
    }

    #[test]
    fn test_cbf_query_various_k_values() -> Result<(), Box<dyn std::error::Error>> {
        let test_cases = vec![
            (4, 4, 2, "small k"),
            (12, 12, 6, "medium k"),
            (20, 20, 10, "large k"),
            (8, 16, 4, "asymmetric k"),
        ];

        for (l_len, r_len, margin_size, description) in test_cases {
            println!("Testing CBF query with: {}", description);
            
            let config = Config {
                l_len,
                r_len,
                chunk_max: 100,
                margin_size,
                bloomfilter_table_size: 256,
                hashset_size: 128,
            };

            // 十分長いテスト配列
            let test_seq = "ATCGATCGATCGATCGATCGTTTTTTTTTTTTTTTTTTGCTAGCTAGCTAGCTAGCTAGC";
            let dna_seq = DnaSequence::new(&test_seq.as_bytes().to_vec());

            // CBFデータ作成
            let cbf_data = vec![100u16; config.bloomfilter_table_size];
            
            let temp_cbf = NamedTempFile::new()?;
            let cbf_path = temp_cbf.path().to_str().unwrap();
            save_cbf_to_file(&cbf_data, &config, cbf_path)?;

            let cbf_query = CbfQuery::from_file(cbf_path, &config)?;
            let results = cbf_query.query_sequence(&dna_seq);

            // 基本検証
            for (l_start, l_end, r_start, r_end, count) in &results {
                assert_eq!(l_end - l_start, l_len, "L-segment length mismatch for {}", description);
                assert_eq!(r_end - r_start, r_len, "R-segment length mismatch for {}", description);
                assert!(*r_start >= l_end + margin_size, "Margin violation for {}", description);
                assert!(*count <= 100, "Count should not exceed CBF values for {}", description);
            }

            println!("  Found {} LR-tuples", results.len());
        }

        Ok(())
    }
}

pub struct CbfQuery {
    cbf_table: Vec<u16>,
    config: Config,
}

impl CbfQuery {
    pub fn from_file(cbf_path: &str, config: &Config) -> Result<Self, Box<dyn std::error::Error>> {
        let cbf_table = load_cbf_from_file(cbf_path, config)?;
        Ok(CbfQuery {
            cbf_table,
            config: *config,
        })
    }

    pub fn query_sequence(&self, sequence: &DnaSequence) -> Vec<(usize, usize, usize, usize, u16)> {
        let mut results = Vec::new();
        
        if sequence.len() < self.config.l_len + self.config.r_len + self.config.margin_size {
            return results;
        }

        let mut l_window_start = 0;
        while l_window_start + self.config.l_len <= sequence.len() {
            let l_window_end = l_window_start + self.config.l_len;
            
            // L-segmentのリピートチェック
            let (l_has_repeat, l_repeat_offset) = sequence.has_repeat(l_window_start, l_window_end);
            if l_has_repeat {
                l_window_start += l_repeat_offset + 1;
                continue;
            }

            let mut r_window_start = l_window_end + self.config.margin_size;
            while r_window_start + self.config.r_len <= sequence.len() {
                let r_window_end = r_window_start + self.config.r_len;
                
                // 最大距離チェック
                if r_window_end - l_window_start > self.config.chunk_max {
                    break;
                }

                // R-segmentのリピートチェック
                let (r_has_repeat, r_repeat_offset) = sequence.has_repeat(r_window_start, r_window_end);
                if r_has_repeat {
                    r_window_start += r_repeat_offset + 1;
                    continue;
                }

                // LR-tupleを生成してCBFに問い合わせ
                let lr_tuple = sequence.subsequence_as_u128(vec![
                    [l_window_start, l_window_end],
                    [r_window_start, r_window_end],
                ]);
                
                let hash_indices = hash_from_u128(lr_tuple, self.config.bloomfilter_table_size);
                let count = count_occurence_from_counting_bloomfilter_table(&self.cbf_table, hash_indices);
                
                results.push((l_window_start, l_window_end, r_window_start, r_window_end, count));
                
                r_window_start += 1;
            }
            l_window_start += 1;
        }
        
        results
    }

    pub fn query_lr_tuple(&self, lr_tuple: u128) -> u16 {
        let hash_indices = hash_from_u128(lr_tuple, self.config.bloomfilter_table_size);
        count_occurence_from_counting_bloomfilter_table(&self.cbf_table, hash_indices)
    }

    pub fn query_sequences_batch(&self, sequences: &[DnaSequence]) -> Vec<Vec<(usize, usize, usize, usize, u16)>> {
        sequences.iter().map(|seq| self.query_sequence(seq)).collect()
    }

    pub fn query_above_threshold(&self, sequence: &DnaSequence, threshold: u16) -> Vec<(usize, usize, usize, usize, u16)> {
        self.query_sequence(sequence)
            .into_iter()
            .filter(|(_, _, _, _, count)| *count >= threshold)
            .collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence_encoder_util::decode_u128_2_dna_seq;
    use std::fs;
    use tempfile::NamedTempFile;

    fn create_test_config() -> Config {
        Config {
            l_len: 10,
            r_len: 10,
            chunk_max: 50,
            margin_size: 5,
            bloomfilter_table_size: 1024, // 小さなサイズでテスト
            hashset_size: 512,
        }
    }

    fn create_test_cbf_data(size: usize) -> Vec<u16> {
        (0..size).map(|i| (i % 65536) as u16).collect()
    }

    fn create_test_fasta(sequences: &[&str]) -> NamedTempFile {
        let temp_file = NamedTempFile::new().unwrap();
        let mut content = String::new();
        for (i, seq) in sequences.iter().enumerate() {
            content.push_str(&format!(">seq{}\n{}\n", i, seq));
        }
        std::fs::write(temp_file.path(), content).unwrap();
        temp_file
    }

    #[test]
    fn test_cbf_header_creation() {
        let config = create_test_config();
        let header = CbfHeader::new(&config);
        
        assert_eq!(header.magic, *CBF_MAGIC);
        assert_eq!(header.table_size, config.bloomfilter_table_size as u64);
        assert!(header.is_valid());
        assert!(header.is_compatible(&config));
    }

    #[test]
    fn test_cbf_header_compatibility() {
        let config1 = create_test_config();
        let mut config2 = create_test_config();
        config2.l_len = 20; // 異なる設定
        
        let header = CbfHeader::new(&config1);
        
        assert!(header.is_compatible(&config1));
        assert!(!header.is_compatible(&config2));
    }

    #[test]
    fn test_cbf_save_and_load() -> Result<(), Box<dyn std::error::Error>> {
        let config = create_test_config();
        let test_data = create_test_cbf_data(config.bloomfilter_table_size);
        
        let temp_file = NamedTempFile::new()?;
        let file_path = temp_file.path().to_str().unwrap();
        
        // 保存
        save_cbf_to_file(&test_data, &config, file_path)?;
        
        // 読み込み
        let loaded_data = load_cbf_from_file(file_path, &config)?;
        
        // 検証
        assert_eq!(test_data.len(), loaded_data.len());
        assert_eq!(test_data, loaded_data);
        
        Ok(())
    }

    #[test]
    fn test_cbf_incompatible_config() -> Result<(), Box<dyn std::error::Error>> {
        let config1 = create_test_config();
        let mut config2 = create_test_config();
        config2.l_len = 20; // 異なる設定
        
        let test_data = create_test_cbf_data(config1.bloomfilter_table_size);
        
        let temp_file = NamedTempFile::new()?;
        let file_path = temp_file.path().to_str().unwrap();
        
        // config1で保存
        save_cbf_to_file(&test_data, &config1, file_path)?;
        
        // config2で読み込み（エラーになるはず）
        let result = load_cbf_from_file(file_path, &config2);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("incompatible"));
        
        Ok(())
    }

    #[test]
    fn test_cbf_invalid_magic() -> Result<(), Box<dyn std::error::Error>> {
        let temp_file = NamedTempFile::new()?;
        let file_path = temp_file.path().to_str().unwrap();
        
        // 無効なマジックナンバーでファイル作成
        let mut file = File::create(file_path)?;
        file.write_all(b"INVALID!")?; // 8バイトの無効なマジック
        file.write_all(&[0u8; 24])?;  // 残りのヘッダー
        file.flush()?;
        
        let config = create_test_config();
        let result = load_cbf_from_file(file_path, &config);
        
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string().contains("Invalid CBF file format"));
        
        Ok(())
    }

    #[test]
    fn test_cbf_file_size_validation() -> Result<(), Box<dyn std::error::Error>> {
        let config = create_test_config();
        let test_data = create_test_cbf_data(config.bloomfilter_table_size);
        
        let temp_file = NamedTempFile::new()?;
        let file_path = temp_file.path().to_str().unwrap();
        
        save_cbf_to_file(&test_data, &config, file_path)?;
        
        // ファイルサイズ確認
        let metadata = fs::metadata(file_path)?;
        let expected_size = 32 + (config.bloomfilter_table_size * 2); // ヘッダー + データ
        assert_eq!(metadata.len(), expected_size as u64);
        
        Ok(())
    }

    #[test]
    fn test_cbf_empty_data() -> Result<(), Box<dyn std::error::Error>> {
        let config = create_test_config();
        let test_data = vec![0u16; config.bloomfilter_table_size];
        
        let temp_file = NamedTempFile::new()?;
        let file_path = temp_file.path().to_str().unwrap();
        
        save_cbf_to_file(&test_data, &config, file_path)?;
        let loaded_data = load_cbf_from_file(file_path, &config)?;
        
        assert_eq!(test_data, loaded_data);
        assert!(loaded_data.iter().all(|&x| x == 0));
        
        Ok(())
    }

    #[test]
    fn test_cbf_max_values() -> Result<(), Box<dyn std::error::Error>> {
        let config = create_test_config();
        let test_data = vec![u16::MAX; config.bloomfilter_table_size];
        
        let temp_file = NamedTempFile::new()?;
        let file_path = temp_file.path().to_str().unwrap();
        
        save_cbf_to_file(&test_data, &config, file_path)?;
        let loaded_data = load_cbf_from_file(file_path, &config)?;
        
        assert_eq!(test_data, loaded_data);
        assert!(loaded_data.iter().all(|&x| x == u16::MAX));
        
        Ok(())
    }

    // 複数のk,mの組み合わせテスト
    #[test]
    fn test_cbf_various_k_m_combinations() -> Result<(), Box<dyn std::error::Error>> {
        let test_cases = vec![
            // (l_len, r_len, margin_size, description)
            (10, 10, 0, "small k, no margin"),
            (15, 15, 5, "medium k, small margin"),
            (20, 20, 10, "large k, medium margin"),
            (8, 12, 3, "asymmetric k, small margin"),
            (25, 15, 20, "large asymmetric k, large margin"),
        ];

        for (l_len, r_len, margin_size, description) in test_cases {
            println!("Testing: {}", description);
            
            let config = Config {
                l_len,
                r_len,
                chunk_max: 100,
                margin_size,
                bloomfilter_table_size: 512, // 小さなサイズでテスト
                hashset_size: 256,
            };

            // テストデータ作成
            let test_data = (0..config.bloomfilter_table_size)
                .map(|i| (i % 1000) as u16)
                .collect::<Vec<u16>>();

            let temp_file = NamedTempFile::new()?;
            let file_path = temp_file.path().to_str().unwrap();

            // 保存・読み込みテスト
            save_cbf_to_file(&test_data, &config, file_path)?;
            let loaded_data = load_cbf_from_file(file_path, &config)?;

            assert_eq!(test_data, loaded_data, "Failed for: {}", description);

            // ファイルサイズ検証
            let metadata = fs::metadata(file_path)?;
            let expected_size = 32 + (config.bloomfilter_table_size * 2);
            assert_eq!(metadata.len(), expected_size as u64, "File size mismatch for: {}", description);
        }

        Ok(())
    }

    #[test]
    fn test_cbf_binary_format_verification() -> Result<(), Box<dyn std::error::Error>> {
        let config = Config {
            l_len: 16,
            r_len: 16,
            chunk_max: 80,
            margin_size: 8,
            bloomfilter_table_size: 256,
            hashset_size: 128,
        };

        // 特定のパターンでテストデータ作成
        let test_data: Vec<u16> = (0..256).map(|i| {
            match i % 4 {
                0 => 0x0000,
                1 => 0xFFFF,
                2 => 0xAAAA,
                3 => 0x5555,
                _ => unreachable!(),
            }
        }).collect();

        let temp_file = NamedTempFile::new()?;
        let file_path = temp_file.path().to_str().unwrap();

        save_cbf_to_file(&test_data, &config, file_path)?;

        // バイナリファイルの内容を直接検証
        let file_content = fs::read(file_path)?;
        
        // ヘッダー検証
        assert_eq!(&file_content[0..8], CBF_MAGIC, "Magic number mismatch");
        
        // テーブルサイズ検証
        let table_size_bytes = &file_content[16..24];
        let table_size = u64::from_le_bytes(table_size_bytes.try_into().unwrap());
        assert_eq!(table_size, 256, "Table size mismatch");

        // データ部分検証
        let data_start = 32;
        for (i, &expected_value) in test_data.iter().enumerate() {
            let offset = data_start + i * 2;
            let actual_bytes = &file_content[offset..offset + 2];
            let actual_value = u16::from_le_bytes(actual_bytes.try_into().unwrap());
            assert_eq!(actual_value, expected_value, "Data mismatch at index {}", i);
        }

        Ok(())
    }

    #[test]
    fn test_cbf_with_real_sequences() -> Result<(), Box<dyn std::error::Error>> {
        let config = Config {
            l_len: 12,
            r_len: 12,
            chunk_max: 60,
            margin_size: 6,
            bloomfilter_table_size: 1024,
            hashset_size: 512,
        };

        // 実際のDNA配列でテスト
        let sequences = vec![
            "ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG",
            "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT",
            "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT",
            "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        ];

        let fasta_file = create_test_fasta(&sequences);
        let fasta_path = fasta_file.path().to_str().unwrap();

        let temp_cbf = NamedTempFile::new()?;
        let cbf_path = temp_cbf.path().to_str().unwrap();

        // CBF構築・保存
        build_and_save_cbf(fasta_path, cbf_path, &config, 1)?;

        // CBF読み込み検証
        let loaded_cbf = load_cbf_from_file(cbf_path, &config)?;
        assert_eq!(loaded_cbf.len(), config.bloomfilter_table_size);

        // 非ゼロ要素の存在確認（実際にデータが処理されたことを確認）
        let non_zero_count = loaded_cbf.iter().filter(|&&x| x > 0).count();
        assert!(non_zero_count > 0, "CBF should contain non-zero elements");

        println!("CBF non-zero elements: {} / {}", non_zero_count, loaded_cbf.len());

        Ok(())
    }

    #[test]
    fn test_lr_tuple_encoding_format() {
        // LR-tupleのエンコーディング形式を検証
        let config = Config {
            l_len: 8,
            r_len: 8,
            chunk_max: 40,
            margin_size: 4,
            bloomfilter_table_size: 256,
            hashset_size: 128,
        };

        // テスト配列: L="AAAAAAAA", margin="TTTT", R="CCCCCCCC"
        let test_seq = "AAAAAAAATTTTCCCCCCCC";
        let dna_seq = DnaSequence::new(&test_seq.as_bytes().to_vec());

        // LR-tupleを抽出
        let lr_tuple = dna_seq.subsequence_as_u128(vec![[0, 8], [12, 20]]);

        // デコードして検証
        let decoded = decode_u128_2_dna_seq(&lr_tuple, 16);
        let decoded_str = String::from_utf8(decoded).unwrap();
        
        assert_eq!(decoded_str, "AAAAAAAACCCCCCCC", "LR-tuple encoding/decoding mismatch");

        // ビット表現の検証
        // A=00, C=01, G=10, T=11
        // "AAAAAAAA" = 0x0000, "CCCCCCCC" = 0x5555
        // 結合: 0x00005555
        assert_eq!(lr_tuple, 0x5555, "LR-tuple bit representation mismatch");
    }

    #[test]
    fn test_cbf_edge_cases() -> Result<(), Box<dyn std::error::Error>> {
        // エッジケースのテスト
        let edge_cases = vec![
            // (l_len, r_len, margin_size, description)
            (1, 1, 0, "minimum k"),
            (32, 32, 0, "maximum k for u128"),
            (16, 16, 100, "large margin"),
            (30, 30, 2, "near maximum k"),
        ];

        for (l_len, r_len, margin_size, description) in edge_cases {
            println!("Testing edge case: {}", description);
            
            let config = Config {
                l_len,
                r_len,
                chunk_max: 200,
                margin_size,
                bloomfilter_table_size: 128,
                hashset_size: 64,
            };

            let test_data = vec![42u16; config.bloomfilter_table_size];
            let temp_file = NamedTempFile::new()?;
            let file_path = temp_file.path().to_str().unwrap();

            save_cbf_to_file(&test_data, &config, file_path)?;
            let loaded_data = load_cbf_from_file(file_path, &config)?;

            assert_eq!(test_data, loaded_data, "Failed for edge case: {}", description);
        }

        Ok(())
    }
}
