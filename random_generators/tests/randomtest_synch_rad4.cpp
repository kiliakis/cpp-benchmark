
#include <iostream>
#include <fstream>
#include <blond/vector_math.h>
#include <blond/utilities.h>
#include <blond/math_functions.h>
#include <blond/optionparser.h>
#include <random>
#include <chrono>
#include <math.h>
#include <stdlib.h>
#include <thread>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <omp.h>

// #include <boost/random.hpp>
// #include <boost/random/normal_distribution.hpp>

using namespace std;
long int N_t = 1;
long int N_p = 100000;
int N_threads = 1;
long int n_kicks = 1;
string prog_name = "synch_rad4";

void parse_args(int argc, char **argv);

int main(int argc, char *argv[])
{
    parse_args(argc, argv);
    omp_set_num_threads(N_threads);
    std::hash<std::thread::id> hash;
    
    ofstream files[N_threads];
    for (int i = 0; i < N_threads; ++i) {
        files[i].open(prog_name+"_"+std::to_string(i) + "_numbers.txt");
    }

    for (int i = 0; i < N_t; i++) {
        #pragma omp parallel
        {
            int tid = omp_get_thread_num();
            // __thread ofstream file(prog_name+"_"+std::to_string(tid) + "_numbers.txt");
            // static __thread ofstream file(prog_name+"_"+std::to_string(tid) + "_numbers.txt");
            static __thread mt19937_64 *gen = nullptr;
            if(!gen) gen = new mt19937_64(clock() + hash(this_thread::get_id()));    
            static __thread std::normal_distribution<> dist(0.0, 1.0);
            #pragma omp for
            for (int j = 0; j < N_p; ++j) {
                files[tid] << dist(*gen) << "\n";
            }
        }
    }

    for (int i = 0; i < N_threads; ++i) {
        files[i].close();
    }


    return 0;
}


void parse_args(int argc, char **argv)
{
    using namespace std;
    using namespace option;

    enum optionIndex {
        UNKNOWN,
        HELP,
        N_THREADS,
        N_TURNS,
        N_PARTICLES,
        N_KICKS,
        OPTIONS_NUM
    };

    const option::Descriptor usage[] = {
        {
            UNKNOWN, 0, "", "", Arg::None, "USAGE: ./myprog [options]\n\n"
            "Options:"
        },
        {
            HELP, 0, "h", "help", Arg::None,
            "  --help,              -h        Print usage and exit."
        },
        {
            N_TURNS, 0, "t", "turns", util::Arg::Numeric,
            "  --turns=<num>,       -t <num>  Number of turns (default: 500)"
        },
        {
            N_KICKS, 0, "k", "n_kicks", util::Arg::Numeric,
            "  --kicks=<num>,       -k <num>  Number of kicks (default: 1)"
        },
        {
            N_PARTICLES, 0, "p", "particles", util::Arg::Numeric,
            "  --particles=<num>,   -p <num>  Number of particles (default: "
            "100k)"
        },
        {
            N_THREADS, 0, "m", "threads", util::Arg::Numeric,
            "  --threads=<num>,     -m <num>  Number of threads (default: 1)"
        },
        {
            UNKNOWN, 0, "", "", Arg::None,
            "\nExamples:\n"
            "\t./myprog\n"
            "\t./myprog -t 1000 -p 10000 -m 4\n"
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
        // fprintf(stdout, "Argument #%d is ", i);
        switch (opt.index()) {
            case HELP:
            // not possible, because handled further above and exits the program
            case N_TURNS:
                N_t = atoi(opt.arg);
                // fprintf(stdout, "--numeric with argument '%s'\n", opt.arg);
                break;
            case N_KICKS:
                n_kicks = atoi(opt.arg);
                // fprintf(stdout, "--numeric with argument '%s'\n", opt.arg);
                break;
            case N_THREADS:
                N_threads = atoi(opt.arg);
                // fprintf(stdout, "--numeric with argument '%s'\n", opt.arg);
                break;
            case N_PARTICLES:
                N_p = atoi(opt.arg);
                // fprintf(stdout, "--numeric with argument '%s'\n", opt.arg);
                break;
            case UNKNOWN:
                // not possible because Arg::Unknown returns ARG_ILLEGAL
                // which aborts the parse with an error
                break;
        }
    }
}
