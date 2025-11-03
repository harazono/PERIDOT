use serde::{Deserialize, Serialize};
use std::fs;
use std::path::Path;

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct Config {
    pub l_len: usize,
    pub r_len: usize,
    pub chunk_max: usize,
    pub margin_size: usize,
    pub bloomfilter_table_size: usize,
    pub hashset_size: usize,
}

impl Default for Config {
    fn default() -> Self {
        Self {
            l_len: 30,
            r_len: 30,
            chunk_max: 200,
            margin_size: 50,
            bloomfilter_table_size: 1 << 30,
            hashset_size: 1 << 29,
        }
    }
}

impl Config {
    pub fn from_file<P: AsRef<Path>>(path: P) -> Result<Self, Box<dyn std::error::Error>> {
        let content = fs::read_to_string(path)?;
        let config: Config = toml::from_str(&content)?;
        Ok(config)
    }

    pub fn save_to_file<P: AsRef<Path>>(&self, path: P) -> Result<(), Box<dyn std::error::Error>> {
        let content = toml::to_string_pretty(self)?;
        fs::write(path, content)?;
        Ok(())
    }

    pub fn total_segment_len(&self) -> usize {
        self.l_len + self.r_len
    }
}
