#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <gmp.h>
#include <mpfr.h>
#include <pthread.h>

#define BITS_PER_DIGIT 3.32192809488736234787
#define C 10939058860032000
#define MAX_SPLIT_DEPTH 3

struct PQT
{
    mpz_t P;
    mpz_t Q;
    mpz_t T;
};

struct bs_args {
    mpz_t a;
    mpz_t b;
    unsigned int depth;
};

void calc_n_iterations( mpz_t n_iterations, mpz_t n_digits )
{
    mpz_div_ui(n_iterations, n_digits, log( C/64 )/log(10));
    mpz_add_ui(n_iterations, n_iterations, 1);
}

void *binary_split(void *void_bs_args)
{
    // Create P, Q, and T struct
    struct PQT *PQT_ab = malloc(sizeof(struct PQT));
    struct bs_args *args = (struct bs_args *)void_bs_args;
    mpz_init(PQT_ab->P);
    mpz_init(PQT_ab->Q);
    mpz_init(PQT_ab->T);

    mpz_t diff;
    mpz_init(diff);
    mpz_sub(diff, args->b, args->a);

    if( mpz_cmp_ui( diff, 1 ) == 0 )
    {
        if( mpz_cmp_ui( args->a, 0 ) == 0)
        {
            // Set P and Q to 1
            mpz_set_ui(PQT_ab->P, 1);
            mpz_set_ui(PQT_ab->Q, 1);
        } else
        {
            // Calculate P
            mpz_set(PQT_ab->Q, args->a);

            mpz_mul_ui(PQT_ab->Q, PQT_ab->Q, 2);
            mpz_sub_ui(PQT_ab->P, PQT_ab->Q, 1);

            mpz_mul_ui(PQT_ab->Q, PQT_ab->Q, 3);
            mpz_sub_ui(PQT_ab->Q, PQT_ab->Q, 1);
            mpz_mul(PQT_ab->P, PQT_ab->P, PQT_ab->Q);

            mpz_sub_ui(PQT_ab->Q, PQT_ab->Q, 4);
            mpz_mul(PQT_ab->P, PQT_ab->P, PQT_ab->Q);

            // Calculate Q
            mpz_set(PQT_ab->Q, args->a);
            mpz_pow_ui(PQT_ab->Q, PQT_ab->Q, 3);
            mpz_mul_ui(PQT_ab->Q, PQT_ab->Q, C);
        }

        // Compute T
        mpz_set_ui(PQT_ab->T, 545140134);
        mpz_mul(PQT_ab->T, PQT_ab->T, args->a);
        mpz_add_ui(PQT_ab->T, PQT_ab->T, 13591409);
        mpz_mul(PQT_ab->T, PQT_ab->T, PQT_ab->P);

        // Switch sign(T) if necessary

        if( mpz_even_p( args->a ) == 0 )
        {
            mpz_neg(PQT_ab->T, PQT_ab->T);
        }
    } else
    {
        // Calculate m
        mpz_t m;
        mpz_init(m);
        mpz_add(m, args->a, args->b);
        mpz_div_ui(m, m, 2);

        struct PQT *PQT_am;
        struct PQT *PQT_mb;

        struct bs_args args_A;
        struct bs_args args_B;
        mpz_init_set(args_A.a, args->a);
        mpz_init_set(args_A.b, m);
        args_A.depth = args->depth+1;

        mpz_init_set(args_B.a, m);
        mpz_init_set(args_B.b, args->b);
        args_B.depth = args->depth+1;

        if( args->depth < MAX_SPLIT_DEPTH )
        {
            pthread_t pid_A;
            pthread_t pid_B;
            void *res_A;
            void *res_B;

            pthread_create(&pid_A, NULL, binary_split, &args_A);
            pthread_create(&pid_B, NULL, binary_split, &args_B);

            pthread_join(pid_A, &res_A);
            pthread_join(pid_B, &res_B);

            PQT_am = (struct PQT *)res_A;
            PQT_mb = (struct PQT *)res_B;
        } else
        {
            PQT_am = (struct PQT *)binary_split(&args_A);
            PQT_mb = (struct PQT *)binary_split(&args_B);
        }

        mpz_gcd(m, PQT_am->P, PQT_mb->Q);
        mpz_div(PQT_am->P, PQT_am->P, m);
        mpz_div(PQT_mb->Q, PQT_mb->Q, m);

        // Calculate PQT
        mpz_mul(PQT_ab->P, PQT_am->P, PQT_mb->P);
        mpz_mul(PQT_ab->Q, PQT_am->Q, PQT_mb->Q);

        mpz_mul(PQT_ab->T, PQT_mb->Q, PQT_am->T);
        mpz_mul(m, PQT_am->P, PQT_mb->T);
        mpz_add(PQT_ab->T, PQT_ab->T, m);

        mpz_clear(PQT_am->P);
        mpz_clear(PQT_am->Q);
        mpz_clear(PQT_am->T);
        mpz_clear(PQT_mb->P);
        mpz_clear(PQT_mb->Q);
        mpz_clear(PQT_mb->T);
        mpz_clear(m);
        free(PQT_am);
        free(PQT_mb);

        // Call Function On Two Halves
        /*
        struct PQT *PQT_am = (struct PQT *)binary_split(a, m, depth+1);
        struct PQT *PQT_mb = (struct PQT *)binary_split(m, b, depth+1);
        */
    }

    return (void *)PQT_ab;
}



int main(int argc, char *argv[])
{
    mpz_t N_DIGITS, ITERATIONS;
    mpfr_prec_t PRECISION;

    mpz_init(N_DIGITS);
    mpz_init(ITERATIONS);

    if( argc == 2 )
    {

        mpz_set_str(N_DIGITS, argv[1], 10);
        PRECISION = (mpfr_prec_t)(strtoll(argv[1], NULL, 10)*BITS_PER_DIGIT + 16);
    }
    else
    {
        mpfr_printf("Invalid Arguments\n");
        return 1;
    }

    calc_n_iterations( ITERATIONS, N_DIGITS );
    mpfr_set_default_prec(PRECISION);

    mpfr_printf("Number Of Digits: %Zd\n", N_DIGITS);
    mpfr_printf("Number Of Bits Per Float: %Pd\n", mpfr_get_default_prec());
    mpfr_printf("Number Of Iterations Required: %Zd\n", ITERATIONS);

    mpfr_printf("Seting Up Required Variables\n");

    mpz_t S;

    mpz_init_set_ui(S, 0);

    mpfr_printf("Starting Computation\n");

    struct bs_args initial_args;
    mpz_init_set(initial_args.a, S);
    mpz_init_set(initial_args.b, ITERATIONS);
    initial_args.depth = 0;

    struct PQT *R =(struct PQT *)binary_split((void *)&initial_args);

    mpfr_printf("Freeing Some Memory\n");
    mpz_clear(S);

    mpfr_t H1;
    mpfr_init(H1);
    mpfr_sqrt_ui(H1, 10005, MPFR_RNDN);
    mpfr_mul_ui(H1, H1, 426880, MPFR_RNDN);
    mpfr_mul_z(H1, H1, R->Q, MPFR_RNDN);
    mpfr_div_z(H1, H1, R->T, MPFR_RNDN);

    FILE *fp;
    fp = fopen("./data/out.txt", "w");
    mpfr_out_str(fp, 10, 0, H1, MPFR_RNDN);

    mpfr_printf("Cleaning Up\n");
    mpz_clear(R->P);
    mpz_clear(R->Q);
    mpz_clear(R->T);
    mpz_clear(N_DIGITS);

    mpz_clear(ITERATIONS);
    mpfr_clear(H1);
    free(R);
}

// gcc chudnovsky.c -lgmp -lpthread -lm -lmpfr -o chudnovsky.out; /usr/bin/time -v ./chudnovsky.out 10000000000; python ./check/calc.py;
//
