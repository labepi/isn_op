use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::io;
use std::{error::Error, fs::File, fs::OpenOptions};
use serde::Serialize;

#[derive(Debug, Serialize)]
pub struct OrdinalPatternData {
    pub entropy: f64,
    pub complexity: f64,
}

#[derive(PartialEq, PartialOrd)]
pub struct FloatWrapper(f64);

impl FloatWrapper {
    pub fn new(val: f64) -> Option<FloatWrapper> {
        if val.is_nan() {
            None
        } else {
            Some(FloatWrapper(val))
        }
    }
}

impl Eq for FloatWrapper {}
impl Ord for FloatWrapper {
    fn cmp(&self, other: &FloatWrapper) -> Ordering {
        self.partial_cmp(other).unwrap()
    }
}

pub fn factorial(num: u64) -> u64 {
    match num {
        0 => 1,
        1 => 1,
        _ => factorial(num - 1) * num,
    }
}

fn arg_sort<T>(data: &[T]) -> Vec<i32> 
where
    T: Ord
{

    let mut indices: Vec<_> = (0..data.len()).collect();
    indices.sort_by_key(|&i| &data[i]);
    indices.into_iter().map(|x| (x as i32) + 1).collect()
}

pub fn ordinal_patterns<T>(series: &Vec<T>, d: i32, tau: i32) -> Vec<Vec<i32>>
where
    T: Ord
{

    let n: i32 = series.len() as i32;
    let mut patterns: Vec<Vec<i32>> = Vec::new();

    let m = n - (d - 1) * tau;

    for i in 0..m {
        let w = i + (d - 1) * tau;
        let subseq = &series[i as usize..=w as usize];
        let pattern = arg_sort(subseq);

        patterns.push(pattern);
    }

    patterns

}

pub fn read_csv(filename: &str) -> Result<Vec<Vec<f64>>, Box<dyn Error>> {
    let mut full_vec = Vec::new();

    let file = File::open(filename)?;
    let mut reader = csv::Reader::from_reader(file);

    for result in reader.records() {
        let record: Vec<f64> = result.unwrap()
            .iter()
            .map(|x| x.parse::<f64>().unwrap()).collect();
        full_vec.push(record);
    }

    Ok(full_vec)
}

pub fn write_csv(
    filename: &str, 
    data: OrdinalPatternData) 
-> io::Result<()> {
    let file = OpenOptions::new().write(true)
        .append(false)
        .create(true)
        .open(filename);
    let mut wtr = csv::Writer::from_writer(file.unwrap());

    wtr.serialize(data)?;
    wtr.flush()?;
    Ok(())
}

pub fn get_column(vec: Vec<Vec<f64>>, col: usize) -> Vec<f64> {
    let mut extracted_columns: Vec<f64> = Vec::new();

    for row in &vec {
        extracted_columns.push(row[col])
    }

    extracted_columns
}

pub fn get_frequency(patterns: &Vec<Vec<i32>>) -> BTreeMap<Vec<i32>, i32> {
    let mut frequency_map: BTreeMap<Vec<i32>, i32> = BTreeMap::new();

    for s in patterns {
        *frequency_map.entry(s.to_vec()).or_default() += 1;
    }

    frequency_map
}

pub fn get_probability(
    pattern_freq: &BTreeMap<Vec<i32>, i32>,
    len: i32,
    d: i32,
    tau: i32
) -> BTreeMap<Vec<i32>, f64> {

    let mut pattern_prob: BTreeMap<Vec<i32>, f64> = BTreeMap::new();
    let m: f64 = (len - (d - 1)*tau).into();
    for (k, v) in pattern_freq {
        *pattern_prob.entry(k.to_vec()).or_default() = (*v as f64/m).into();
    }

    pattern_prob

}

pub fn get_entropy(
    pattern_prob: &BTreeMap<Vec<i32>, f64>,
) -> BTreeMap<Vec<i32>, f64> {

    let mut pattern_entropy: BTreeMap<Vec<i32>, f64> = BTreeMap::new();

    for (k, v) in pattern_prob {
        *pattern_entropy.entry(k.to_vec()).or_default() = v * v.ln();
    }

    pattern_entropy
}

pub fn normalized_entropy(
    pattern_entropy: &Vec<f64>,
    d: i32,
) -> f64 {

    let mut sum = 0.0;

    for v in pattern_entropy {
        sum += v;
    }

    sum = sum * -1.0;

    let d_fact: f64 = factorial(d as u64) as f64;
    let norm_entropy = sum/d_fact.ln();

    norm_entropy

}

pub fn shannon_entropy(
    pattern_prob: &Vec<f64>,
    norm: bool,
) -> f64 {

    let mut sum = 0.0;

    let n = pattern_prob.len() as f64;
    for v in pattern_prob {
        sum += v * v.ln();
    }

    let mut entropy = sum * -1.0;

    if norm {
        entropy = entropy/n.ln();
    }

    entropy

}

pub fn complexity(
    vec_prob: &Vec<f64>,
    mut entropy: Option<f64>
) -> f64 {
    if let None = entropy {
        entropy = Some(shannon_entropy(&vec_prob, true));
    }

    let n: f64 = vec_prob.len() as f64;

    let p_u: Vec<f64> = vec![1.0/n; n as usize];

    let sum_vecs: Vec<f64>= vec_prob
        .iter()
        .zip(&p_u)
        .map(|(a, b)| (a + b) / 2.0)
        .collect();

    let entropy_sum: f64 = shannon_entropy(&sum_vecs, false);
    let entropy_prob = shannon_entropy(&vec_prob, false);
    let entropy_uni = shannon_entropy(&p_u, false);

    let js = (entropy_sum) - (entropy_prob/2.0) - (entropy_uni/2.0);

    let aux = 
        ((n+1.0)/n) * (n + 1.0).ln() - 2.0*((2.0*n).ln()) + n.ln();

    let q_0: f64;
    if aux != 0.0 {
        q_0 = -2.0 * (1.0/aux);
    } else {
        q_0 = -2.0
    }

    let q = q_0 * js;
    let c = q * entropy.unwrap();

    c
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_arg_sort() {
        let data = vec![4, 7, 9];

        assert_eq!(vec![1, 2, 3], arg_sort(&data))
    }

    #[test]
    fn test_ordinal_patterns() {
        let d = 3;
        let tau = 1;
        let series = vec![4, 7, 9 , 10, 6, 11, 3];

        let test_vector: Vec<Vec<i32>> = vec![
            vec![1, 2, 3],
            vec![1, 2, 3],
            vec![3, 1, 2],
            vec![2, 1, 3],
            vec![3, 1, 2],
        ];

        assert_eq!(test_vector, ordinal_patterns(&series, d, tau))
    }

    #[test]
    fn test_float_ordinal_patterns() {
        let d = 3;
        let tau = 1;
        let series: Vec<_> = vec![4.2, 7.1, 9.3, 10.3, 6.2, 11.3, 3.2]
            .iter()
            .map(|x| FloatWrapper::new(*x).unwrap()).collect();

        let test_vector: Vec<Vec<i32>> = vec![
            vec![1, 2, 3],
            vec![1, 2, 3],
            vec![3, 1, 2],
            vec![2, 1, 3],
            vec![3, 1, 2],
        ];
        assert_eq!(test_vector, ordinal_patterns(&series, d, tau))
    }

}
