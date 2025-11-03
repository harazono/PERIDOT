use std::collections::HashSet;
use std::fs;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::Path;
use std::process::Command;
use std::sync::mpsc;
use std::thread;
use flate2::write::GzEncoder;
use flate2::Compression;

fn get_db_taxids(db: &str) -> HashSet<String> {
    let output = Command::new("blastdbcmd")
        .args(&["-db", db, "-entry", "all", "-outfmt", "%T"])
        .output()
        .expect("failed to execute process");
    let stdout = String::from_utf8_lossy(&output.stdout);
    stdout.lines().map(|s| s.to_string()).collect()
}

fn call(
    tgt_db_list_filename: &str,
    taxon_id_to_be_ignored: &str,
    fasta_file: &str,
    output_dir: &str,
    slots: usize,
) {
    // negative_taxon_idsの読み込み
    let negative_taxon_ids: HashSet<String> = BufReader::new(File::open(taxon_id_to_be_ignored).unwrap())
        .lines()
        .map(|l| l.unwrap().trim().to_string())
        .collect();

    // 出力ディレクトリとログディレクトリの作成
    fs::create_dir_all(output_dir).expect("failed to create output directory");
    let log_dir = Path::new(output_dir).join("logs");
    fs::create_dir_all(&log_dir).expect("failed to create log directory");

    // 出力ファイルの設定
    let output_file_path = Path::new(output_dir).join("blastn_results.gz");

    // データベースリストの読み込み
    let db_list: Vec<String> = BufReader::new(File::open(tgt_db_list_filename).unwrap())
        .lines()
        .filter_map(|l| {
            let line = l.unwrap();
            if line.trim().is_empty() {
                None
            } else {
                Some(line.trim().to_string())
            }
        })
        .collect();

    // スレッドプール
    let (tx, rx) = mpsc::channel();
    let mut handles = vec![];

    // 各データベースに対してblastnを実行
    for (i, db) in db_list.iter().enumerate() {
        let tx = tx.clone();
        let db = db.clone();
        let fasta_file = fasta_file.to_string();
        let output_dir = output_dir.to_string();
        let log_dir = log_dir.clone();
        let negative_taxon_ids = negative_taxon_ids.clone();

        let handle = thread::spawn(move || {
            // ログファイルの設定
            let log_file = log_dir.join(format!("blastn_{}.log", i));
            let exec_cmd_logfile = log_dir.join(format!("blastn_{}_exec_cmd.txt", i));

            let temp_output = Path::new(&output_dir).join(format!("{}.temp.gz", db));

            let taxon_ids_in_db = get_db_taxids(&db);
            let negative_taxon_ids_for_this_db: Vec<String> = taxon_ids_in_db
                .intersection(&negative_taxon_ids)
                .cloned()
                .collect();

            let negative_taxids_arg = if !negative_taxon_ids_for_this_db.is_empty() {
                format!("-negative_taxids {}", negative_taxon_ids_for_this_db.join(","))
            } else {
                "".to_string()
            };

            let execute_cmd = format!(
                r#"echo {negative_taxids_arg} > {log_file:?}; \
                /home/harazono/miniconda3/bin/blastn \
                -task blastn-short \
                -query {fasta_file:?} \
                -db {db:?} \
                -outfmt "6 qseqid sseqid sacc qlen qstart qend slen sstart send qseq sseq evalue length staxid staxids ssciname scomname" \
                -num_threads {slots} \
                -dust no \
                -soft_masking false \
                -word_size 15 \
                -best_hit_overhang 0.01 \
                -best_hit_score_edge 0.49 \
                -perc_identity 100 \
                {negative_taxids_arg} \
                2>> {log_file:?} | gzip > {temp_output:?}
            "#
            );

            // 実行コマンドのログ
            let mut f = File::create(exec_cmd_logfile).unwrap();
            writeln!(f, "{}", execute_cmd).unwrap();

            // blastnの実行
            let output = Command::new("sh")
                .arg("-c")
                .arg(execute_cmd)
                .output()
                .expect("failed to execute process");

            // blastnの出力
            tx.send((i, output)).unwrap();
        });
        handles.push(handle);
    }

    // 各スレッドの終了を待つ
    for handle in handles {
        handle.join().unwrap();
    }

    // blastnの結果を結合
    let f_out = File::create(output_file_path).unwrap();
    let mut f_out = BufWriter::new(GzEncoder::new(f_out, Compression::default()));
    for _ in 0..db_list.len() {
        let (i, output) = rx.recv().unwrap();
        println!("{} in {}", i, db_list.len());

        if !output.status.success() {
            let stderr = String::from_utf8_lossy(&output.stderr);
            panic!("blastn failed: {}", stderr);
        }

        // 一時ファイルの削除
        let temp_output = Path::new(output_dir).join(format!("{}.temp.gz", db_list[i]));
        let mut f_in = BufReader::new(File::open(temp_output.clone()).unwrap());
        std::io::copy(&mut f_in, &mut f_out).unwrap();
        fs::remove_file(temp_output).unwrap();
    }
}

fn main() {
    let args: Vec<String> = std::env::args().collect();

    if args.len() != 6 {
        eprintln!(
            "Usage: {} <tgt_db_list_filename> <taxon_id_to_be_ignored> <fasta_file> <slots> <output_dir>",
            args[0]
        );
        std::process::exit(1);
    }

    let tgt_db_list_filename = &args[1];
    let taxon_id_to_be_ignored = &args[2];
    let fasta_file = &args[3];
    let slots: usize = args[4].parse().unwrap();
    let output_dir = &args[5];

    call(
        tgt_db_list_filename,
        taxon_id_to_be_ignored,
        fasta_file,
        output_dir,
        slots,
    );
}