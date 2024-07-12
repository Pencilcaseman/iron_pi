use clap::Parser;
use colored::Colorize;
use iron_pi::{
    bsplit::binary_split,
    fact::{PrimeFactorSieve, PrimeFactors},
};
use rug::{Float, Integer};
use std::io::Write;
use std::ops::{AddAssign, MulAssign};

const BITS_PER_DIGIT: f64 = std::f64::consts::LOG2_10;
const DIGITS_PER_ITER: f64 = 14.181647462725478;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Number of digits to calculate
    #[clap(short, long, default_value_t = 1000)]
    digits: usize,

    /// File to write the result to (default is `pi.txt`)
    #[clap(short, long, default_value = "pi.txt")]
    out_file: String,

    /// Number of digits per block (default is 10)
    #[clap(short, long, default_value_t = 10)]
    block_size: usize,

    /// Number of blocks per line (default is 5)
    #[clap(short, long, default_value_t = 5)]
    num_blocks: usize,
}

fn format_with_commas(num: usize) -> String {
    let mut s = num.to_string();
    let mut i = s.len() as isize - 3;
    while i > 0 {
        s.insert(i as usize, ',');
        i -= 3;
    }
    s
}

fn main() {
    println!("{}", "PI Calculator".red().bold());

    let Args {
        digits,
        out_file,
        block_size,
        num_blocks,
    } = Args::parse();

    let prec = (digits as f64 * BITS_PER_DIGIT) as u32 + 1;
    let iters = ((digits as f64) / DIGITS_PER_ITER) as usize + 1;

    println!(
        "{} {}",
        "Digits    : ".green(),
        format_with_commas(digits).cyan().bold(),
    );
    println!(
        "{} {}",
        "Precision : ".green(),
        format!("{} bits", format_with_commas(prec as usize))
            .cyan()
            .bold()
    );
    println!(
        "{} {}",
        "Iterations: ".green(),
        format_with_commas(iters).cyan().bold()
    );
    println!();

    print!("{}", "Generating factor sieve... ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();
    let sieve = PrimeFactorSieve::new((3 * 5 * 23 * 29 + 1).max(6 * iters));
    let end = std::time::Instant::now();
    println!(
        "{} {}",
        "Done in".green(),
        format!("{:?}", end - start).cyan()
    );

    print!("{}", "Binary splitting...        ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();
    let (_, q, r) = binary_split(1, iters, 0, &sieve);
    let q = q.num;
    let r = r.num;
    let end = std::time::Instant::now();
    println!(
        "{} {}",
        "Done in".green(),
        format!("{:?}", end - start).cyan()
    );

    print!("{}", "Calculating numerator...   ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();
    let mut num = Integer::from(426880);
    num.mul_assign(&q);
    let end = std::time::Instant::now();
    println!(
        "{} {}",
        "Done in".green(),
        format!("{:?}", end - start).cyan()
    );

    print!("{}", "Calculating denominator... ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();
    let mut den = Integer::from(13591409);
    den.mul_assign(q);
    den.add_assign(r);
    let end = std::time::Instant::now();
    println!(
        "{} {}",
        "Done in".green(),
        format!("{:?}", end - start).cyan()
    );

    print!("{}", "Dividing...                ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();
    let div = Float::with_val(prec, num) / Float::with_val(prec, den);
    let end = std::time::Instant::now();
    println!(
        "{} {}",
        "Done in".green(),
        format!("{:?}", end - start).cyan()
    );

    print!("{}", "Square root...             ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();
    let sqrt = Float::with_val(prec, 10005).sqrt();
    let end = std::time::Instant::now();
    println!(
        "{} {}",
        "Done in".green(),
        format!("{:?}", end - start).cyan()
    );

    print!("{}", "Final multiplication...    ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();
    let pi = div * sqrt;
    let end = std::time::Instant::now();
    println!(
        "{} {}",
        "Done in".green(),
        format!("{:?}", end - start).cyan()
    );
    println!();

    print!("{}", "Converting to string...    ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();
    let pi_bytes: Vec<u8> = pi
        .to_string_radix(10, Some(digits + 1))
        .into_bytes()
        .into_iter()
        .skip(2)
        .collect();
    let end = std::time::Instant::now();
    println!(
        "{} {}",
        "Done in".green(),
        format!("{:?}", end - start).cyan()
    );

    print!("{}", "Formatting string...       ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();

    let formatted_pi: Vec<u8> = pi_bytes
        .into_iter()
        .enumerate()
        .flat_map(|(i, c)| {
            // let pos = i + 1;
            let pos = i;
            if pos % (block_size * num_blocks) == 0 {
                // Start of a new line
                vec![b'\n', b' ', b' ', c]
            } else if pos % block_size == 0 {
                // Start of a new block (but not a new line)
                vec![b' ', c]
            } else {
                // Regular digit
                vec![c]
            }
        })
        .skip(1) // Skip the initial newline
        .collect();

    let end = std::time::Instant::now();
    println!(
        "{} {}",
        "Done in".green(),
        format!("{:?}", end - start).cyan()
    );

    print!("{}", "Writing to file...         ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();
    let mut file = std::fs::File::create(out_file).unwrap();
    file.write_all(b"3.\n").unwrap();
    // file.write_all(&pi_bytes).unwrap();
    file.write_all(&formatted_pi).unwrap();
    let end = std::time::Instant::now();
    println!(
        "{} {}",
        "Done in".green(),
        format!("{:?}", end - start).cyan()
    );
}
