use std::{io::Write, num::NonZero, thread::available_parallelism};

use clap::Parser;
use colored::Colorize;
use iron_pi::bsplit::binary_split;

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

// unsafe fn arb_const_pi_chudnovsky_eval(s: &mut flint3_sys::arb_t, prec: i64)
// {     let mut series = iron_pi::util::new_hypgeom();
//
//     let mut t = iron_pi::util::new_arb();
//     let mut u = iron_pi::util::new_arb();
//
//     flint3_sys::arb_init(&mut t[0]);
//     flint3_sys::arb_init(&mut u[0]);
//     flint3_sys::hypgeom_init(&mut series[0]);
//
//     flint3_sys::fmpz_poly_set_str(
//         &mut series[0].A[0],
//         std::ffi::CString::new("2  13591409 545140134").unwrap().as_ptr(),
//     );
//     flint3_sys::fmpz_poly_set_str(
//         &mut series[0].B[0],
//         std::ffi::CString::new("1  1").unwrap().as_ptr(),
//     );
//     flint3_sys::fmpz_poly_set_str(
//         &mut series[0].P[0],
//         std::ffi::CString::new("4  5 -46 108 -72").unwrap().as_ptr(),
//     );
//     flint3_sys::fmpz_poly_set_str(
//         &mut series[0].Q[0],
//         std::ffi::CString::new("4  0 0 0
// 10939058860032000").unwrap().as_ptr(),     );
//
//     // prec += FLINT_CLOG2(prec) + 5;
//     flint3_sys::arb_hypgeom_infsum(
//         &mut s[0],
//         &mut t[0],
//         &mut series[0],
//         prec,
//         prec,
//     );
//
//     let mut s_thing = iron_pi::util::new_arb();
//     flint3_sys::arb_log_base_ui(&mut s_thing[0], &s[0], 10, 128);
//     let s_string = std::ffi::CString::from_raw(flint3_sys::arb_get_str(
//         &s_thing[0],
//         prec,
//         0,
//     ));
//
//     let mut t_thing = iron_pi::util::new_arb();
//     flint3_sys::arb_log_base_ui(&mut t_thing[0], &t[0], 10, 128);
//     let t_string = std::ffi::CString::from_raw(flint3_sys::arb_get_str(
//         &t_thing[0],
//         prec,
//         0,
//     ));
//
//     println!("S THING: {s_string:?}");
//     println!("T THING: {t_string:?}");
//
//     flint3_sys::arb_rsqrt_ui(&mut u[0], 640320, prec);
//     flint3_sys::arb_mul(&mut s[0], &s[0], &u[0], prec);
//
//     flint3_sys::arb_mul_ui(&mut t[0], &t[0], 640320 / 12, prec);
//     flint3_sys::arb_div(&mut s[0], &t[0], &s[0], prec);
//
//     flint3_sys::hypgeom_clear(&mut series[0]);
//     flint3_sys::arb_clear(&mut t[0]);
//     flint3_sys::arb_clear(&mut u[0]);
// }

fn main() {
    println!("{}", "============== IronPI ==============".red().bold());

    let Args {
        digits,
        mut max_parallel_depth,
        out_file,
        mut base,
        block_size,
        num_blocks,
        mut threads,
    } = Args::parse();

    // Set the number of threads
    rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build_global()
        .unwrap();

    let prec = (digits as f64 * BITS_PER_DIGIT) as u64 + 16;
    // let iters = ((digits as f64) * 1.25 / DIGITS_PER_ITER) as u64 + 16;
    let iters = ((digits as f64) / DIGITS_PER_ITER) as u64 + 16;
    let max_depth = iters.ilog2();

    if threads == 0 {
        threads = available_parallelism()
            .unwrap_or(NonZero::new(1usize).unwrap())
            .get();
    }

    if max_parallel_depth == 0 {
        max_parallel_depth = (threads.ilog2() + 1) as usize;
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

    let global_start = std::time::Instant::now();

    print!("{}", "Binary splitting...        ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();

    let (mut tmp, mut q, mut r) =
        binary_split(1, iters, &pool, max_parallel_depth, prec as i64);

    let end = std::time::Instant::now();
    println!("{} {}", "Done in".green(), format!("{:?}", end - start).cyan());

    let mut u = iron_pi::util::new_arb();

    print!("{}", "Calculating numerator...   ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();

    unsafe {
        flint3_sys::arb_mul_ui(&mut tmp[0], &q[0], 426_880, prec as i64);
    }

    let end = std::time::Instant::now();
    println!("{} {}", "Done in".green(), format!("{:?}", end - start).cyan());

    print!("{}", "Calculating denominator... ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();

    unsafe {
        flint3_sys::arb_addmul_ui(&mut r[0], &q[0], 13_591_409, prec as i64);
    }

    let end = std::time::Instant::now();
    println!("{} {}", "Done in".green(), format!("{:?}", end - start).cyan());

    print!("{}", "Dividing...                ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();

    unsafe {
        flint3_sys::arb_div(&mut tmp[0], &tmp[0], &r[0], prec as i64);
    }

    let end = std::time::Instant::now();
    println!("{} {}", "Done in".green(), format!("{:?}", end - start).cyan());

    print!("{}", "Square root...             ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();

    unsafe {
        flint3_sys::arb_sqrt_ui(&mut q[0], 10_005, prec as i64);
    }

    let end = std::time::Instant::now();
    println!("{} {}", "Done in".green(), format!("{:?}", end - start).cyan());

    print!("{}", "Final multiplication...    ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();

    unsafe {
        flint3_sys::arb_mul(&mut tmp[0], &tmp[0], &q[0], prec as i64);
    }

    let end = std::time::Instant::now();
    println!("{} {}", "Done in".green(), format!("{:?}", end - start).cyan());
    println!();

    println!(
        "{} {}\n",
        "Calculated pi in".green(),
        format!("{:?}", global_start.elapsed()).cyan()
    );

    print!("{}", "Freeing Memory...          ".green());
    std::io::stdout().flush().unwrap();
    let start = std::time::Instant::now();

    unsafe {
        flint3_sys::arb_clear(&mut q[0]);
        flint3_sys::arb_clear(&mut r[0]);
    }

    let end = std::time::Instant::now();
    println!("{} {}", "Done in".green(), format!("{:?}", end - start).cyan());
    println!();

    print!("{}", "Converting to string...    ".green());
    std::io::stdout().flush().unwrap();
    let str_start = std::time::Instant::now();

    let base_digits = (digits as f64 / (base as f64).log10()).round() as usize;

    let pi_bytes = unsafe {
        let bytes = flint3_sys::arb_get_str(&tmp[0], base_digits as i64 + 1, 2);
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
