use std::io::Write;

use clap::Parser;
use colored::Colorize;
use iron_pi::bsplit::{binary_split, binary_split_arb};

const BITS_PER_DIGIT: f64 = std::f64::consts::LOG2_10;
const DIGITS_PER_ITER: f64 = 14.181647462725478;

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Number of digits to calculate
    #[clap(short, long, default_value_t = 1000)]
    digits: u64,

    #[clap(short, long, default_value_t = 0)]
    max_parallel_depth: usize,

    /// File to write the result to
    #[clap(short, long, default_value = "pi.txt")]
    out_file: String,

    /// Number of digits per block
    #[clap(short, long, default_value_t = 10)]
    block_size: usize,

    #[clap(short = 'B', long, default_value_t = 10)]
    base: i32,

    /// Number of blocks per line
    #[clap(short, long, default_value_t = 5)]
    num_blocks: usize,

    /// Number of threads to use (defaults to all available threads)
    #[clap(short, long, default_value_t = 0)]
    threads: usize,
}

fn format_with_commas(num: u64) -> String {
    let mut s = num.to_string();
    let mut i = s.len() as i64 - 3;
    while i > 0 {
        s.insert(i as usize, ',');
        i -= 3;
    }
    s
}

fn main() {
    println!("{}", "============== IronPI ==============".red().bold());

    let Args {
        digits,
        mut max_parallel_depth,
        out_file,
        mut base,
        block_size,
        num_blocks,
        threads,
    } = Args::parse();

    // Set the number of threads
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();

    // let prec = (digits as f64 * BITS_PER_DIGIT) as u32 + 16;
    let prec = (digits as f64 * BITS_PER_DIGIT) as u64 + 16;
    let iters = ((digits as f64) * 1.25 / DIGITS_PER_ITER) as u64 + 16;
    let max_depth = iters.ilog2();

    if max_parallel_depth == 0 {
        max_parallel_depth = (threads.ilog2() + 1) as usize;
    }

    unsafe {
        // We need the most precision possible with MPFR. Without this,
        // MPFR cannot convert the numerator and denominator into floats
        // after more than ~100,000,000 digits
        let max = gmp_mpfr_sys::mpfr::get_emax_max();
        gmp_mpfr_sys::mpfr::set_emax(max);
    }

    if !(2..=36).contains(&base) {
        println!(
            "{}",
            "Warning: Base must be in the range [2, 36]".red().bold()
        );

        base = base.clamp(2, 36);
        println!("{}", format!("Warning: Setting base to {base}").red().bold());
    }

    println!(
        "{} {}",
        "Digits          : ".green(),
        format_with_commas(digits).cyan().bold(),
    );

    println!(
        "{} {}",
        "Precision       : ".green(),
        format!("{} bits", format_with_commas(prec)).cyan().bold()
    );

    println!(
        "{} {}",
        "Iterations      : ".green(),
        format_with_commas(iters).cyan().bold()
    );

    println!(
        "{} {}",
        "Max depth       : ".green(),
        format_with_commas(max_depth as u64).cyan().bold()
    );

    println!(
        "{} {}",
        "Parallel depth  : ".green(),
        format_with_commas(max_parallel_depth as u64).cyan().bold()
    );

    println!(
        "{} {}",
        "Threads         : ".green(),
        format_with_commas(threads as u64).cyan().bold()
    );

    println!();

    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .expect("Failed to build pool");

    unsafe {
        flint3_sys::flint_set_num_threads(threads as i32);
    }

    // unsafe {
    //     let start = std::time::Instant::now();
    //     let mut tmp_pi = iron_pi::util::new_arb();
    //     flint3_sys::arb_const_pi(&mut tmp_pi[0], prec as i64);
    //     println!("Elapsed: {:?}\n", start.elapsed());
    //     flint3_sys::arb_clear(&mut tmp_pi[0]);
    // }

    let global_start = std::time::Instant::now();

    print!("{}", "Binary splitting...        ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();

    // let (_, q_full, r_full) =
    //     binary_split_work_stealing(1, iters, &sieve, &pool);

    // let (_, q, mut r) = binary_split(1, iters, &pool, max_parallel_depth);

    let (_, q, mut r) =
        binary_split_arb(1, iters, &pool, max_parallel_depth, prec as i64);

    // let q = q_full.num;
    // let mut r = r_full.num;

    let end = std::time::Instant::now();
    println!("{} {}", "Done in".green(), format!("{:?}", end - start).cyan());

    print!("{}", "Calculating numerator...   ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();
    // let mut num = Integer::from(426880);

    let mut num = iron_pi::util::new_arb_with_u64(426_880);

    // num.mul_assign(&q);
    unsafe {
        // flint3_sys::fmpz_mul(&mut num[0], &num[0], &q[0]);

        flint3_sys::arb_mul(&mut num[0], &num[0], &q[0], prec as i64);
    }

    let end = std::time::Instant::now();
    println!("{} {}", "Done in".green(), format!("{:?}", end - start).cyan());

    // println!(
    //     "\t {} {}",
    //     format_with_commas(unsafe {
    //         // gmp_mpfr_sys::gmp::mpz_sizeinbase(num.as_raw(), 10)
    //         flint3_sys::fmpz_sizeinbase(&num[0], 10) as u64
    //     })
    //     .truecolor(255, 47, 106)
    //     .bold(),
    //     "digits".truecolor(255, 47, 106),
    // );

    print!("{}", "Calculating denominator... ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();

    // den = 13591409 * q + r = r + q * 13591409
    // let mut den = Integer::from(13591409);
    // den.mul_assign(q);
    // den.add_assign(r);

    // let mut den = iron_pi::util::new_fmpz_with(13_591_409);
    unsafe {
        // flint3_sys::fmpz_addmul_ui(&mut r[0], &q[0], 13_591_409);

        flint3_sys::arb_addmul_ui(&mut r[0], &q[0], 13_591_409, prec as i64);
    }

    let end = std::time::Instant::now();
    println!("{} {}", "Done in".green(), format!("{:?}", end - start).cyan());

    // println!(
    //     "\t {} {}",
    //     format_with_commas(unsafe {
    //         // gmp_mpfr_sys::gmp::mpz_sizeinbase(den.as_raw(), 10)
    //         flint3_sys::fmpz_sizeinbase(&r[0], 10) as u64
    //     })
    //     .truecolor(255, 47, 106)
    //     .bold(),
    //     "digits".truecolor(255, 47, 106),
    // );

    print!("{}", "Dividing...                ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();

    // let num = Float::with_val_64(prec, num);
    // let den = Float::with_val_64(prec, den);
    // let num = iron_pi::util::new_arb_with_fmpz(&num);
    // let den = iron_pi::util::new_arb_with_fmpz(&r);

    // if num.is_infinite() || num.is_nan() {
    //     println!("{}", "Numerator is 'infinite' or 'NaN'.".red().bold());
    //     return;
    // }

    // if den.is_infinite() || den.is_nan() {
    //     println!("{}", "Denominator is 'infinite' or 'NaN'.".red().bold());
    //     return;
    // }

    // let div = num / den;

    let mut div = iron_pi::util::new_arb();
    unsafe {
        // flint3_sys::arb_fmpz_div_fmpz(&mut div[0], &num[0], &r[0], prec as
        // i64);

        flint3_sys::arb_div(&mut div[0], &num[0], &r[0], prec as i64);
    }

    let end = std::time::Instant::now();
    println!("{} {}", "Done in".green(), format!("{:?}", end - start).cyan());

    print!("{}", "Square root...             ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();

    // let sqrt = Float::with_val_64(prec, 10005).sqrt();

    let mut sqrt = iron_pi::util::new_arb();

    unsafe {
        flint3_sys::arb_sqrt_ui(&mut sqrt[0], 10_005, prec as i64);
    }

    let end = std::time::Instant::now();
    println!("{} {}", "Done in".green(), format!("{:?}", end - start).cyan());

    print!("{}", "Final multiplication...    ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();

    // let pi = div * sqrt;
    unsafe {
        flint3_sys::arb_mul(&mut div[0], &div[0], &sqrt[0], prec as i64);
    }

    let end = std::time::Instant::now();
    println!("{} {}", "Done in".green(), format!("{:?}", end - start).cyan());
    println!();

    println!(
        "{} {}\n",
        "Calculated pi in".green(),
        format!("{:?}", global_start.elapsed()).cyan()
    );

    print!("{}", "Converting to string...    ".green());
    std::io::stdout().flush().unwrap();
    let str_start = std::time::Instant::now();

    let base_digits = (digits as f64 / (base as f64).log10()).round() as usize;

    // let pi_bytes: Vec<u8> = pi
    //     .to_string_radix(base, Some(base_digits + 1))
    //     .into_bytes()
    //     .into_iter()
    //     .skip(2)
    //     .collect();

    let pi_bytes = unsafe {
        let bytes = flint3_sys::arb_get_str(&div[0], base_digits as i64 + 1, 2);
        &&std::ffi::CStr::from_ptr(bytes)
            .to_str()
            .expect("Failed to convert to string")
            .as_bytes()[2..]
    };

    let end = std::time::Instant::now();
    println!(
        "{} {}",
        "Done in".green(),
        format!("{:?}", end - str_start).cyan()
    );

    // print!("{}", "Formatting string...       ".green());
    // std::io::stdout().flush().unwrap();
    // let start = std::time::Instant::now();
    //
    // let formatted_pi: Vec<u8> = pi_bytes
    //     .par_iter()
    //     .enumerate()
    //     .flat_map(|(pos, &c)| {
    //         if pos % (block_size * num_blocks) == 0 {
    //             // Start of a new line
    //             vec![b'\n', b' ', b' ', c]
    //         } else if pos % block_size == 0 {
    //             // Start of a new block (but not a new line)
    //             vec![b' ', c]
    //         } else {
    //             // Regular digit
    //             vec![c]
    //         }
    //     })
    //     .collect();
    //
    // let end = std::time::Instant::now();
    // println!("{} {}", "Done in".green(), format!("{:?}", end -
    // start).cyan());
    //
    // print!("{}", "Writing to file...         ".green());
    // std::io::stdout().flush().unwrap();
    // let start = std::time::Instant::now();
    // let mut file = std::fs::File::create(out_file).unwrap();
    // file.write_all(b"3.\n").unwrap();
    // file.write_all(&formatted_pi).unwrap();
    // let end = std::time::Instant::now();
    // println!("{} {}", "Done in".green(), format!("{:?}", end -
    // start).cyan());

    // ========================================================================
    // ========================================================================
    // ========================================================================

    print!("{}", "Writing string...          ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();

    let file = std::fs::File::create(out_file).expect("Failed to create file");
    let mut writer = std::io::BufWriter::new(file);

    writer.write_all(b"3.").expect("Failed to write prefix");

    pi_bytes.iter().enumerate().for_each(|(pos, &c)| {
        let result = if pos % (block_size * num_blocks) == 0 {
            writer.write_all(b"\n  ")
        } else if pos % block_size == 0 {
            writer.write_all(b" ")
        } else {
            Ok(())
        };

        result
            .and_then(|_| writer.write_all(&[c]))
            .expect("Failed to write to file");
    });

    writer.flush().expect("Failed to flush writer");

    println!(
        "{} {}",
        "Done in".green(),
        format!("{:?}", start.elapsed()).cyan()
    );

    // ========================================================================
    // ========================================================================
    // ========================================================================

    println!(
        "\n{} {}\n",
        "IO completed in".green(),
        format!("{:?}", str_start.elapsed()).cyan()
    );

    let global_end = std::time::Instant::now();
    println!(
        "{} {}",
        "Total time: ".green(),
        format!("{:?}", global_end - global_start).cyan()
    );
}
