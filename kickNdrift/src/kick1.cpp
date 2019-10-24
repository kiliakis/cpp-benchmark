
#include <iostream>
#include <blond/vector_math.h>
#include <blond/utilities.h>
#include <blond/math_functions.h>
#include <blond/optionparser.h>
#include <random>
#include <chrono>

using namespace std;

long int N_t = 1;
long int N_p = 100000;
int N_threads = 1;
long int N_rf = 1;


inline void kick(const double *__restrict beam_dt,
                 double *__restrict beam_dE,
                 const int n_rf,
                 const double *__restrict voltage,
                 const double *__restrict omega_rf,
                 const double *__restrict phi_rf,
                 const int n_macroparticles,
                 const double acc_kick)
{
    // KICK
    //#pragma omp parallel for collapse(2)
    for (int j = 0; j < n_rf; ++j) {
        #pragma omp parallel for
        for (int i = 0; i < n_macroparticles; ++i) {
            const double a = omega_rf[j] * beam_dt[i] + phi_rf[j];
            beam_dE[i] += voltage[j] * mymath::fast_sin(a);
        }
    }

// SYNCHRONOUS ENERGY CHANGE
    #pragma omp parallel for
    for (int i = 0; i < n_macroparticles; ++i)
        beam_dE[i] += acc_kick;
}


void parse_args(int argc, char **argv);

int main(int argc, char *argv[])
{
    parse_args(argc, argv);
    omp_set_num_threads(N_threads);
    cout << "Number of turns: " <<  N_t << "\n";
    cout << "Number of points: " << N_p << "\n";
    cout << "Number of rf sections: " << N_rf << "\n";
    cout << "Number of openmp threads: " <<  N_threads << "\n";

    vector<double> dE(N_p), dt(N_p);
    vector<double> omega_rf(N_rf), voltage(N_rf), phi_rf(N_rf);

    default_random_engine generator;
    uniform_real_distribution<double> d(0.0, 1.0);
    // apply_f_in_place(dt, [d, generator](double x) {return d(generator);});
    // apply_f_in_place(dE, [d, generator](double x) {return d(generator);});
    // apply_f_in_place(omega_rf, [d, generator](double x) {return d(generator);});
    // apply_f_in_place(voltage, [d, generator](double x) {return d(generator);});
    // apply_f_in_place(phi_rf, [d, generator](double x) {return d(generator);});
    
    for (int i = 0; i < N_p; i++) {
        dt[i] = 10e-6 * d(generator);
        dE[i] = 10e6 * d(generator);
    }

    for (int i = 0; i < N_rf; i++) {
        omega_rf[i] = d(generator);
        voltage[i] = 10000 * d(generator);
        phi_rf[i] = d(generator);
    }


    double sum = 0.0;
    chrono::time_point<chrono::high_resolution_clock>start;
    chrono::duration<double> elapsed_time(0.0);
    start = chrono::system_clock::now();

    for (int i = 0; i < N_t; i++) {
        kick(dt.data(), dE.data(), N_rf, voltage.data(),
             omega_rf.data(), phi_rf.data(), N_p, 0.0);
    }

    elapsed_time = chrono::system_clock::now() - start;
    sum = mymath::sum(dE);

    double throughput = N_t * N_p * N_rf / elapsed_time.count() / 1e6;
    cout << "kick with fast_sin\n";
    cout << "Elapsed time: " << elapsed_time.count() << " sec\n";
    cout << "Throughput: " << throughput << " MP/sec\n";
    cout << "Throughput/thread: " << throughput/N_threads << " MP/sec\n";
    cout << "Sum: " << sum << "\n";

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
        N_RF,
        OPTIONS_NUM
    };

    const option::Descriptor usage[] = {
        {
            UNKNOWN, 0, "", "", Arg::None, "USAGE: ./kick1 [options]\n\n"
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
            N_RF, 0, "r", "rf_sections", util::Arg::Numeric,
            "  --rf_sections=<num>,       -r <num>  Number of rf_sections (default: 1)"
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
            "\t./kick1\n"
            "\t./kick1 -t 1000 -p 10000 -m 4\n"
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
            case N_RF:
                N_rf = atoi(opt.arg);
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
