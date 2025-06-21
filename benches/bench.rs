use criterion::{criterion_group, criterion_main, Criterion};

fn bench_pi(_crit: &mut Criterion) {
    todo!()
}

criterion_group!(benches, bench_pi);
criterion_main!(benches);
