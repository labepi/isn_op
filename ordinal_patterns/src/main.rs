use std::{fs::ReadDir, io::ErrorKind};
use std::fs;
use ordinal_patterns::{FloatWrapper, OrdinalPatternData};

fn main() {

    let paths = fs::read_dir("../data/new_csv/windows/").unwrap();

    // setup series
    let d = 7;
    let tau = 1;

    // generate_random_results(d, tau);
    generate_isn_results(paths, d, tau);
}
fn compute_data(
    filename: &str,
    d: i32,
    tau: i32,
    column: usize 
) -> OrdinalPatternData {

    let full_vec = ordinal_patterns::read_csv(filename).unwrap();
    let series = ordinal_patterns::get_column(full_vec, column);
    let series: Vec<FloatWrapper> = series
        .iter()
        .map(|x| FloatWrapper::new(*x).unwrap()).collect();

    let result = ordinal_patterns::ordinal_patterns(&series, d, tau);
    let freq = ordinal_patterns::get_frequency(&result);

    // probability
    let prob = ordinal_patterns::get_probability(
        &freq, 
        series.len() as i32,
        d,
        tau
    );

    // entropy & normalized entropy
    // let entropy = ordinal_patterns::get_entropy(&prob);
    // let vec_entropy: Vec<f64> = entropy.values().cloned().collect();

    // let _norm_entropy = ordinal_patterns::normalized_entropy(
    //     &vec_entropy, 
    //     d
    // );


    // complexity
    let vec_prob: Vec<f64> = prob.values().cloned().collect();
    let shannon_entropy = ordinal_patterns::shannon_entropy(&vec_prob, true);
    let complexity = ordinal_patterns::complexity(&vec_prob, None);

    let data_struct = OrdinalPatternData {
        entropy: shannon_entropy,
        complexity
    };

    data_struct

}
fn generate_random_results(d: i32, tau: i32) {
    let filename = "../dados_random.csv";

    fs::create_dir_all(format!("../results/random_{d}"))
        .expect("couldn't create the directory");

    for i in 0..4 {
        let data_struct = compute_data(filename, d, tau, i);

        let filename = filename.split('/').collect::<Vec<&str>>()[1];

        let result_path = format!(
            "../results/random_{d}/{i}_{filename}"
        );

        match ordinal_patterns::write_csv(&result_path, data_struct) {
            Err(e) if e.kind() == ErrorKind::AlreadyExists => {
                println!("file already exists")
            }
            _ => println!("Wrote {result_path}"),
        }
    }
}
fn generate_isn_results(paths: ReadDir, d: i32, tau: i32) {
    for path in paths {
        let file = path.unwrap().path().to_str().unwrap().to_owned();
        let data_struct = compute_data(&file, d, tau, 0);
        
        // set the result file with the same name as the original file
        let file = file.split('/').collect::<Vec<&str>>()
            .last()
            .unwrap()
            .to_string()
            .to_owned();
        fs::create_dir_all(format!("../results/d_{d}"))
            .expect("couldn't create the directory");

        let result_path = format!("../results/d_{d}/{file}");
        match ordinal_patterns::write_csv(&result_path, data_struct) {
            Err(e) if e.kind() == ErrorKind::AlreadyExists => {
                continue;
            }
            _ => println!("Wrote {result_path}"),
        }

    }
}
