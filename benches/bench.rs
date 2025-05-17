use criterion::{black_box, criterion_group, criterion_main, Criterion};
use iron_pi::fact::*;

fn bench_prime_sieve(crit: &mut Criterion) {
    crit.bench_function("PrimeFactorSieve::new(1_000_000)", |bench| {
        bench.iter(|| PrimeFactorSieve::new(1_000_000))
    });
}

fn bench_prime_factors(crit: &mut Criterion) {
    let sieve = PrimeFactorSieve::new(10_000_000);
    let a = PrimeFactors::new(&sieve, 2);
    let b = PrimeFactors::new(&sieve, 3);
    let c = PrimeFactors::new(&sieve, 4);
    let d = PrimeFactors::new(&sieve, 5);
    let e = PrimeFactors::new(&sieve, 6);

    crit.bench_function("a * b * c * d * e", |bench| {
        bench.iter(|| {
            black_box(&a)
                * black_box(&b)
                * black_box(&c)
                * black_box(&d)
                * black_box(&e)
        })
    });

    crit.bench_function("PrimeFactors::new(3_628_800)", |b| {
        b.iter(|| PrimeFactors::new(&sieve, black_box(3_628_800)))
    });
}

fn bench_split_mul(crit: &mut Criterion) {
    let sieve = PrimeFactorSieve::new(100_000_000);

    let numbers = vec![
        PrimeFactors::new(&sieve, 37618973),
        PrimeFactors::new(&sieve, 20571526),
        PrimeFactors::new(&sieve, 24775278),
        PrimeFactors::new(&sieve, 73066629),
        PrimeFactors::new(&sieve, 81616626),
        PrimeFactors::new(&sieve, 40049046),
        PrimeFactors::new(&sieve, 53335747),
        PrimeFactors::new(&sieve, 23004181),
        PrimeFactors::new(&sieve, 54313362),
        PrimeFactors::new(&sieve, 12116873),
        PrimeFactors::new(&sieve, 40534212),
        PrimeFactors::new(&sieve, 21841054),
        PrimeFactors::new(&sieve, 24501285),
        PrimeFactors::new(&sieve, 27515545),
        PrimeFactors::new(&sieve, 44154554),
        PrimeFactors::new(&sieve, 44886639),
        PrimeFactors::new(&sieve, 73567365),
        PrimeFactors::new(&sieve, 32437521),
    ];

    let m = numbers.iter().fold(PrimeFactors::new(&sieve, 1), |a, b| a * b);

    crit.bench_function("bench_split_mul", |bench| {
        bench.iter(|| black_box(&m).split_mul(0, black_box(m.0.len())))
    });
}

fn bench_div_gcd(crit: &mut Criterion) {
    let sieve = PrimeFactorSieve::new(100_000_000);

    let sieve = PrimeFactorSieve::new(100_000_000);

    let numbers = vec![
        PrimeFactors::new(&sieve, 37618973),
        PrimeFactors::new(&sieve, 20571526),
        PrimeFactors::new(&sieve, 24775278),
        PrimeFactors::new(&sieve, 73066629),
        PrimeFactors::new(&sieve, 81616626),
        PrimeFactors::new(&sieve, 40049046),
        PrimeFactors::new(&sieve, 53335747),
        PrimeFactors::new(&sieve, 23004181),
        PrimeFactors::new(&sieve, 54313362),
        PrimeFactors::new(&sieve, 12116873),
        PrimeFactors::new(&sieve, 40534212),
        PrimeFactors::new(&sieve, 21841054),
        PrimeFactors::new(&sieve, 24501285),
        PrimeFactors::new(&sieve, 27515545),
        PrimeFactors::new(&sieve, 44154554),
        PrimeFactors::new(&sieve, 44886639),
        PrimeFactors::new(&sieve, 73567365),
        PrimeFactors::new(&sieve, 32437521),
    ];

    let mut m = numbers.iter().fold(PrimeFactors::new(&sieve, 1), |a, b| a * b);

    crit.bench_function("bench_div_gcd", |bench| {
        bench.iter(|| {
            numbers.iter().for_each(|a| {
                let _ = black_box(&m).div_gcd(a);
            })
        });
    });
}

criterion_group!(
    benches,
    bench_prime_sieve,
    bench_prime_factors,
    bench_split_mul,
    bench_div_gcd
);
criterion_main!(benches);
