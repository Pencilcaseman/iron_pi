use clap::Parser;
use colored::Colorize;
use iron_pi::{bsplit::par_split, fact::PrimeFactorSieve};
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

    /// Depth of the parallel binary splitting (default is 4)
    #[clap(short, long, default_value_t = 4)]
    lookup_depth: usize,

    /// File to write the result to (default is `pi.txt`)
    #[clap(short, long, default_value = "pi.txt")]
    out_file: String,

    /// Number of digits per block (default is 10)
    #[clap(short, long, default_value_t = 10)]
    block_size: usize,

    /// Number of blocks per line (default is 5)
    #[clap(short, long, default_value_t = 5)]
    num_blocks: usize,

    /// Whether to divide by the GCD (default is true)
    #[clap(short, long)]
    gcd: bool,
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
        lookup_depth,
        out_file,
        block_size,
        num_blocks,
        gcd: div_gcd,
    } = Args::parse();

    let prec = (digits as f64 * BITS_PER_DIGIT) as u32 + 16; // Double the precision
    let iters = ((digits as f64) * 1.25 / DIGITS_PER_ITER) as usize + 16;
    let max_depth = iters.ilog2();

    unsafe {
        // We need the most precision possible with MPFR. Without this,
        // MPFR cannot convert the numerator and denominator into floats
        // after more than ~100,000,000 digits
        let max = gmp_mpfr_sys::mpfr::get_emax_max();
        gmp_mpfr_sys::mpfr::set_emax(max);
    }

    println!(
        "{} {}",
        "Digits        : ".green(),
        format_with_commas(digits).cyan().bold(),
    );

    println!(
        "{} {}",
        "Precision     : ".green(),
        format!("{} bits", format_with_commas(prec as usize))
            .cyan()
            .bold()
    );

    println!(
        "{} {}",
        "Iterations    : ".green(),
        format_with_commas(iters).cyan().bold()
    );

    println!(
        "{} {}",
        "Max depth     : ".green(),
        format_with_commas(max_depth as usize).cyan().bold()
    );

    println!(
        "{} {}",
        "Parallel depth: ".green(),
        format_with_commas(lookup_depth).cyan().bold()
    );

    println!();

    let global_start = std::time::Instant::now();

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
    // let (_, q_full, r_full) = binary_split(1, iters, 0, &sieve, usize::MAX, &mut HashMap::new());
    // let (_, q_full, r_full) = binary_split_par(1, iters, 0, &sieve, lookup_depth, None);

    let mut lookup = par_split(1, iters, 0, &sieve, lookup_depth);
    let (_, q_full, r_full) = lookup.remove(&[1, iters]).unwrap();
    let mut q = q_full.num;
    let mut r = r_full.num;
    let end = std::time::Instant::now();
    println!(
        "{} {}",
        "Done in".green(),
        format!("{:?}", end - start).cyan()
    );

    if div_gcd {
        print!("{}", "Finding GCD...             ".green());
        std::io::stdout().flush().unwrap();
        let start = std::time::Instant::now();
        let mut gcd = q.clone();
        gcd.gcd_mut(&r);
        let end = std::time::Instant::now();
        println!(
            "{} {}",
            "Done in".green(),
            format!("{:?}", end - start).cyan()
        );

        print!("{}", "Dividing by GCD...         ".green());
        std::io::stdout().flush().unwrap();
        let start = std::time::Instant::now();
        q.div_exact_mut(&gcd);
        r.div_exact_mut(&gcd);
        let end = std::time::Instant::now();
        println!(
            "{} {}",
            "Done in".green(),
            format!("{:?}", end - start).cyan()
        );
    } else {
        println!("{}", "Skipping GCD...".purple().bold());
    }

    print!("{}", "Calculating numerator...   ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();
    let mut num = Integer::from(426880);
    num.mul_assign(&q);
    let end = std::time::Instant::now();
    print!(
        "{} {}",
        "Done in".green(),
        format!("{:?}", end - start).cyan()
    );
    println!(
        "\t {} {}",
        format_with_commas(unsafe { gmp_mpfr_sys::gmp::mpz_sizeinbase(num.as_raw(), 10) })
            .truecolor(255, 47, 106)
            .bold(),
        "digits".truecolor(255, 47, 106),
    );

    print!("{}", "Calculating denominator... ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();
    let mut den = Integer::from(13591409);
    den.mul_assign(q);
    den.add_assign(r);
    let end = std::time::Instant::now();
    print!(
        "{} {}",
        "Done in".green(),
        format!("{:?}", end - start).cyan()
    );
    println!(
        "\t {} {}",
        format_with_commas(unsafe { gmp_mpfr_sys::gmp::mpz_sizeinbase(den.as_raw(), 10) })
            .truecolor(255, 47, 106)
            .bold(),
        "digits".truecolor(255, 47, 106),
    );

    print!("{}", "Dividing...                ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();
    // let div = Float::with_val(prec, num) / Float::with_val(prec, den);

    let num = Float::with_val(prec, num);
    let den = Float::with_val(prec, den);

    if num.is_infinite() || num.is_nan() {
        println!("{}", "Numerator is 'infinite' or 'NaN'.".red().bold());
        return;
    }

    if den.is_infinite() || den.is_nan() {
        println!("{}", "Denominator is 'infinite' or 'NaN'.".red().bold());
        return;
    }

    let div = num / den;

    // let div = if div_gcd {
    //     unsafe { Rational::from_canonical(num, den) }
    // } else {
    //     Rational::from((num, den))
    // };

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

    println!();

    let global_end = std::time::Instant::now();
    println!(
        "{} {}",
        "Total time: ".green(),
        format!("{:?}", global_end - global_start).cyan()
    );

    println!();
}
