# Chudnovsky

Requires GMP and MPFR to be installed and included while compiling. To compile:

`
gcc filename.c -lm -lgmp -lmpfr;
`

Then, run with first parameter as number of digits. For example: 

`
./a.out 100 
`

will produce 100 digits of pi and output to ./data/out.txt (make sure this file in this folder exists). This is not a final version and their may be small bugs and inconveniences.
