# Iron Pi

Iron Pi is able to calculate billions\* of digits of pi extremely quickly.

It applies binary splitting to the Chudnovsky algorithm, combined with
some reasonably complicated multithreading.

## Usage

```none
PI Calculator
Usage: iron_pi [OPTIONS]

Options:
  -d, --digits <DIGITS>              Number of digits to calculate [default: 1000]
  -l, --lookup-depth <LOOKUP_DEPTH>  Depth of the parallel binary splitting [default: 4]
  -o, --out-file <OUT_FILE>          File to write the result to [default: pi.txt]
  -b, --block-size <BLOCK_SIZE>      Number of digits per block [default: 10]
  -n, --num-blocks <NUM_BLOCKS>      Number of blocks per line [default: 5]
  -g, --gcd                          Whether to divide by the GCD
  -t, --threads <THREADS>            Number of threads to use (defaults to all available threads) [default: 0]
  -h, --help                         Print help
  -V, --version                      Print version
```
