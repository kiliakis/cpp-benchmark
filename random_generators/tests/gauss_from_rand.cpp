#include <iostream>
#include <chrono>
#include <vector>
#include <cmath>
#include <configuration.h>
#include <utilities.h>
#include <optionparser.h>
#include <omp.h>
#include <math_functions.h>
// #include <rng64.h>
#include <cstdlib>
#include <limits>

unsigned long long  seed;

const double epsilon = std::numeric_limits<double>::min();
const double two_pi = 2.0 * 3.14159265358979323846;

double generateGaussianNoise(double mu, double sigma)
{

    static double z0, z1;
    static bool generate;
    generate = !generate;

    if (!generate)
        return z1 * sigma + mu;

    double u1, u2;
    do {
        u1 = rand() * (1.0 / RAND_MAX);
        u2 = rand() * (1.0 / RAND_MAX);
    } while (u1 <= epsilon);

    z0 = sqrt(-2.0 * log(u1)) * cos(two_pi * u2);
    z1 = sqrt(-2.0 * log(u1)) * sin(two_pi * u2);
    return z0 * sigma + mu;
}


unsigned N = 100000;
unsigned ITERS = 10;
unsigned n_threads  = 1;
int dump_on = 0;
void parse_args(int argc, char **argv);

int main(int argc, char **argv)
{
    parse_args(argc, argv);

    std::cout << "Num of Iterations : " << ITERS << "\n";
    std::cout << "Num of Elems/Iteration : " << N << "\n";
    std::cout << "Num of Threads : " << n_threads << "\n";
    omp_set_num_threads(n_threads);

    std::cout << "\n\n";
    std::cout.precision(4);

    timespec start;
    auto elapsed = 0.0L;
    // auto sum = 0.0L;
    auto mean = 0.0, std = 0.0;
    f_vector_t out(N);

    // seed = 1;
    // initrand(&seed, 1, 32, 16, 16);

    for (unsigned iter = 0; iter < ITERS; ++iter) {

        util::get_time(start);

        for (uint i = 0; i < out.size(); i++)
            out[i] = generateGaussianNoise(0, 1);

        elapsed += util::time_elapsed(start);//end - start;

        mean += mymath::mean(out.data(), out.size());
        std += mymath::standard_deviation(out.data(), out.size());
        if (dump_on)
            util::dump_to_file(out, "gauss_from_rand.dat");
    }



    std::cout << "Generation of " << N << " elems\n";
    std::cout << "Elapsed Time : " << elapsed << " s\n";

    std::cerr << "Throughput : " << (N * ITERS * sizeof(double)) / (elapsed * 1000000) << " MB/s\n";

    std::cout << "Mean : " << mean / ITERS << std::endl;
    std::cout << "std : " << std / ITERS << std::endl;
    std::cout << "\n\n";

    return 0;
}


void parse_args(int argc, char **argv)
{
    using namespace std;
    using namespace option;

    enum optionIndex {UNKNOWN, HELP, N_ELEMS, N_ITERS, N_THREADS, DUMP, OPTIONS_NUM
                     };

    const option::Descriptor usage[] = {
        {
            UNKNOWN, 0, "", "", Arg::None,                  "USAGE: ./my_prog [options]\n\n"
            "Options:"
        },
        {  HELP, 0, "h", "help", Arg::None, "--help, -h  Print usage and exit." },
        {N_ELEMS, 0, "n", "elems", util::Arg::Numeric, "--elems=<num>, -n <num>  Number of elems (default: 10k)" },
        {N_ITERS, 0, "i", "iters", util::Arg::Numeric, "--iters=<num>, -i <num>  Number of iterations (default: 10k)" },
        {N_THREADS, 0, "t", "threads", util::Arg::Numeric, "--threads=<num>, -t <num>  Number of Threads (default: 1)" },
        {DUMP, 0, "d", "dump", util::Arg::Numeric, "--dump=<0|1>, -d <0|1>  Dump to file (default: 0)" },

        {
            UNKNOWN, 0, "", "", Arg::None, "\nExamples:\n"
            "\t./my_prog\n"
            "\t./my_prog -n 1000\n"
        },
        {0, 0, 0, 0, 0, 0}
    };

    argc -= (argc > 0);
    argv += (argc > 0); // skip program name argv[0] if present
    Stats stats(usage, argc, argv);
    vector<Option> options(stats.options_max);
    vector<Option> buffer(stats.buffer_max);
    Parser parse(usage, argc, argv, &options[0], &buffer[0]);

    if (options[HELP]) {
        printUsage(cout, usage);
        exit(0);
    }

    for (int i = 0; i < parse.optionsCount(); ++i) {
        Option &opt = buffer[i];
        //fprintf(stdout, "Argument #%d is ", i);
        switch (opt.index()) {
            case HELP:
            // not possible, because handled further above and exits the program
            case N_ELEMS:
                N = atoi(opt.arg);
                //fprintf(stdout, "--numeric with argument '%s'\n", opt.arg);
                break;
            case N_ITERS:
                ITERS = atoi(opt.arg);
                //fprintf(stdout, "--numeric with argument '%s'\n", opt.arg);
                break;
            case N_THREADS:
                n_threads = atoi(opt.arg);
                //fprintf(stdout, "--numeric with argument '%s'\n", opt.arg);
                break;
            case DUMP:
                dump_on = atoi(opt.arg);
                //fprintf(stdout, "--numeric with argument '%s'\n", opt.arg);
                break;
            case UNKNOWN:
                // not possible because Arg::Unknown returns ARG_ILLEGAL
                // which aborts the parse with an error
                break;
        }
    }


}

