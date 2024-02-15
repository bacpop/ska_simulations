use std::fs::File;
use std::io::prelude::*;
use std::io::{self, BufReader, BufWriter};

use rand::distributions::WeightedIndex;
use rand::prelude::*;
use rand_distr::{Distribution, Gamma};

use clap::Parser;

// These can be hard coded
const GAMMA_SCALE: f64 = 5.0;
const GAMMA_SHAPE: f64 = 1.0;
const INVAR_RATE: f64 = 0.01;

const ALPHA: f64 = 0.9850577;
const BETA: f64 = 0.8185017;
const GAMMA: f64 = 0.3737484;
const DELTA: f64 = 0.5641096;
const EPS: f64 = 0.2259800;
const ETA: f64 = 1.0;

const PI_A: f64 = 0.301807845032044;
const PI_C: f64 = 0.197996188588821;
const PI_G: f64 = 0.199163647868774;
const PI_T: f64 = 0.301032318510361;

const BASE_VEC: [u8; 4] = [b'A', b'C', b'G', b'T'];
const A_MUT: [char; 3] = ['C', 'G', 'T'];
const C_MUT: [char; 3] = ['A', 'G', 'T'];
const G_MUT: [char; 3] = ['A', 'C', 'T'];
const T_MUT: [char; 3] = ['A', 'C', 'G'];

const INDEL_DIST: [f64; 20] = [
    0.320238687,
    0.140726007,
    0.134261561,
    0.085529587,
    0.055693685,
    0.050223769,
    0.032819493,
    0.026852312,
    0.032819493,
    0.019890602,
    0.018896072,
    0.020885132,
    0.012928891,
    0.007956241,
    0.008950771,
    0.008950771,
    0.006464446,
    0.006961711,
    0.004972650,
    0.003978120,
];

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
struct Args {
    /// File to mutate (fasta with single contig)
    #[arg(long, default_value = "pneumo.fa")]
    ref_file: String,

    /// Output file for split-kmers
    #[arg(long, default_value = "split_kmers.txt")]
    kmer_file: String,

    /// K-mer size
    #[arg(short, default_value_t = 31)]
    k: usize,

    /// Branch length (mutations per site)
    #[arg(short, long, default_value_t = 0.0001)]
    end_dist: f64,

    /// Rate of INDELs (compared to SNPs at rate 1.0)
    #[arg(short, long, default_value_t = 0.175)]
    indel_rate: f64,
}

struct GtrDraws {
    a_rates: WeightedIndex<f64>,
    c_rates: WeightedIndex<f64>,
    g_rates: WeightedIndex<f64>,
    t_rates: WeightedIndex<f64>,
}

// Linear approx to exp(Ql) as l << 1
impl GtrDraws {
    pub fn new(
        alpha: f64,
        beta: f64,
        gamma: f64,
        delta: f64,
        eps: f64,
        eta: f64,
        pi_a: f64,
        pi_c: f64,
        pi_g: f64,
        pi_t: f64,
    ) -> Self {
        let a_rates = WeightedIndex::new(vec![alpha * pi_g, beta * pi_c, gamma * pi_t]).unwrap();
        let c_rates = WeightedIndex::new(vec![alpha * pi_a, delta * pi_c, eps * pi_t]).unwrap();
        let g_rates = WeightedIndex::new(vec![beta * pi_a, delta * pi_g, eta * pi_t]).unwrap();
        let t_rates = WeightedIndex::new(vec![gamma * pi_a, eps * pi_g, eta * pi_c]).unwrap();

        Self {
            a_rates,
            c_rates,
            g_rates,
            t_rates,
        }
    }

    pub fn draw_mutation(&self, base: char, rng: &mut ThreadRng) -> char {
        match base {
            'a' | 'A' => A_MUT[self.a_rates.sample(rng)],
            'c' | 'C' => C_MUT[self.c_rates.sample(rng)],
            'g' | 'G' => G_MUT[self.g_rates.sample(rng)],
            't' | 'T' | 'u' | 'U' => T_MUT[self.t_rates.sample(rng)],
            _ => panic!("Mutating non ACGT(U) base"),
        }
    }
}

// Generate gamma heterogeneity
fn site_rates(n_sites: usize, gamma_dist: &Gamma<f64>, rng: &mut ThreadRng) -> Vec<f64> {
    let mut weights = Vec::with_capacity(n_sites);
    for _idx in 0..n_sites {
        let u1: f64 = rng.gen();
        let weight = if u1 > INVAR_RATE {
            gamma_dist.sample(rng)
        } else {
            0.0
        };
        weights.push(weight);
    }
    weights
}

fn open_ref(file: &str) -> io::Result<String> {
    let f = File::open(file)?;
    let f = BufReader::new(f);

    let mut seq = String::new();
    for line in f.lines().skip(1) {
        seq.push_str(&line.unwrap());
    }

    Ok(seq)
}

fn main() {
    let args = Args::parse();

    // Read input and open split-k output file
    let mut start_seq = open_ref(&args.ref_file).unwrap().as_bytes().to_vec();
    let mut k_file = BufWriter::new(File::create(args.kmer_file).unwrap());

    // RNG and simulation progress
    let mut rng = rand::thread_rng();
    let split_k = (args.k - 1) / 2;
    let mut sites = start_seq.len();
    let mut dist = 0.0;

    // Set up GTR+GAMMA+I model
    let gtr_matrix = GtrDraws::new(ALPHA, BETA, GAMMA, DELTA, EPS, ETA, PI_A, PI_C, PI_G, PI_T);
    let gamma = Gamma::new(GAMMA_SHAPE, GAMMA_SCALE).unwrap();
    let mut weights = site_rates(sites, &gamma, &mut rng);
    let mut pos_dist = WeightedIndex::new(&weights).unwrap();
    let indel_size_dist = WeightedIndex::new(INDEL_DIST).unwrap();
    let base_dist = WeightedIndex::new(vec![PI_A, PI_C, PI_G, PI_T]).unwrap();

    // Gillespie algorithm, sum of rates is 1
    let r_indel = args.indel_rate / (1.0 + args.indel_rate);
    // not needed as only two possibilities
    // let r_sub = 1.0 / (1.0 + INDEL_RATE);
    while dist < (args.end_dist * start_seq.len() as f64) {
        let mutated_pos = pos_dist.sample(&mut rng);
        //eprintln!("{mutated_pos}");
        let u1: f64 = rng.gen();
        dist += -u1.ln();

        let u2: f64 = rng.gen();
        if u2 < r_indel {
            // Make an INDEL (in/del equally likely)
            let size = indel_size_dist.sample(&mut rng) + 1;
            let u3: f64 = rng.gen();
            if u3 > 0.5 {
                // Insertion
                sites += size;
                let mut new_seq = vec![start_seq[mutated_pos]];
                for _idx in 0..size {
                    new_seq.push(BASE_VEC[base_dist.sample(&mut rng)]);
                }
                /*
                eprintln!(
                    "IN {mutated_pos}:{}",
                    String::from_utf8(new_seq.clone()).unwrap()
                );
                */
                start_seq.splice(mutated_pos..(mutated_pos + 1), new_seq);

                let mut new_weights = vec![weights[mutated_pos]];
                new_weights.append(&mut site_rates(size, &gamma, &mut rng));
                weights.splice(mutated_pos..(mutated_pos + 1), new_weights);
            } else {
                // Deletion
                sites -= size;
                /*
                eprintln!(
                    "DEL {mutated_pos}:{}",
                    String::from_utf8(start_seq[mutated_pos..(mutated_pos + size)].to_vec())
                        .unwrap()
                );
                */
                start_seq.drain(mutated_pos..(mutated_pos + size));
                weights.drain(mutated_pos..(mutated_pos + size));
            }
            // Need new gamma heterogeneity
            pos_dist = WeightedIndex::new(&weights).unwrap();
        } else {
            // Make a substitution
            let old_base = start_seq[mutated_pos] as char;
            start_seq[mutated_pos] =
                gtr_matrix.draw_mutation(start_seq[mutated_pos] as char, &mut rng) as u8;

            // Write out the split k-mer in ska nk --full-info format
            let k_start = mutated_pos.saturating_sub(split_k);
            let k_end = (mutated_pos + split_k).min(sites - 1);
            writeln!(
                k_file,
                "{}\t{}\t{},{}",
                String::from_utf8(start_seq[k_start..mutated_pos].to_vec()).unwrap(),
                String::from_utf8(start_seq[(mutated_pos + 1)..=k_end].to_vec()).unwrap(),
                old_base,
                start_seq[mutated_pos] as char
            )
            .unwrap();
        }
    }

    println!(">mutated");
    println!("{}", String::from_utf8(start_seq).unwrap());
}
