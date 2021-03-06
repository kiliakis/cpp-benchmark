
#include <iostream>
#include <blond/vector_math.h>
#include <blond/utilities.h>
#include <blond/math_functions.h>
#include <blond/optionparser.h>
#include <random>
#include <chrono>
#include <math.h>
#include <algorithm>    // std::fill_n
#include <stdlib.h>
#include <string.h>     // memset()

// #include <boost/random.hpp>
// #include <boost/random/normal_distribution.hpp>

using namespace std;
string prog_name = "histo2";
long int N_t = 1;
int N_p = 100000;
int N_threads = 1;
int n_slices = 100;
chrono::duration<double> set_zero_d(0.0);
chrono::duration<double> histo_calc_d(0.0);
chrono::duration<double> free_mem_d(0.0);


// This function calculates and applies only the synchrotron radiation damping term
extern "C" void histogram(const double *__restrict__ input,
                                 double *__restrict__ output, const double cut_left,
                                 const double cut_right, const int n_slices,
                                 const int n_macroparticles)
{

    // Number of Iterations of the inner loop
    // const int STEP = 16;
    const double inv_bin_width = n_slices / (cut_right - cut_left);

    chrono::time_point<chrono::high_resolution_clock> start_t;
    start_t = chrono::system_clock::now();

    double **histo = (double **) malloc(omp_get_max_threads() * sizeof(double *));
    histo[0] = (double *) malloc (omp_get_max_threads() * n_slices * sizeof(double));
    for (int i = 0; i < omp_get_max_threads(); i++)
        histo[i] = (*histo + n_slices * i);
    
    set_zero_d += chrono::system_clock::now() - start_t; 
    
    start_t = chrono::system_clock::now();

    #pragma omp parallel
    {
        const int id = omp_get_thread_num();
        const int threads = omp_get_num_threads();
        memset(histo[id], 0., n_slices * sizeof(double));
        // float fbin[STEP];
        #pragma omp for
        for (int i = 0; i < n_macroparticles; i ++) {
            int bin = floor((input[i] - cut_left) * inv_bin_width);
            if (bin < 0 || bin >= n_slices) continue;
            histo[id][bin] += 1.0;
            // const int loop_count = n_macroparticles - i > STEP ?
            //                        STEP : n_macroparticles - i;

            // First calculate the index to update
            // for (int j = 0; j < loop_count; j++) {
            //     fbin[j] = floor((input[i + j] - cut_left) * inv_bin_width);
            // }
            // // Then update the corresponding bins
            // for (int j = 0; j < loop_count; j++) {
            //     const int bin  = (int) fbin[j];
            //     if (bin < 0 || bin >= n_slices) continue;
            //     histo[id][bin] += 1.;
            // }
        }

        // Reduce to a single histogram
        #pragma omp for
        for (int i = 0; i < n_slices; i++) {
            output[i] = 0.;
            for (int t = 0; t < threads; t++)
                output[i] += histo[t][i];
        }
    }
    
    histo_calc_d += chrono::system_clock::now() - start_t; 

    start_t = chrono::system_clock::now();
    // free memory
    free(histo[0]);
    free(histo);
    free_mem_d += chrono::system_clock::now() - start_t; 



}


void parse_args(int argc, char **argv);

int main(int argc, char *argv[])
{
    parse_args(argc, argv);
    omp_set_num_threads(N_threads);
    cout << "Number of turns: " <<  N_t << "\n";
    cout << "Number of points: " << N_p << "\n";
    cout << "Number of slices: " << n_slices << "\n";
    cout << "Number of openmp threads: " <<  N_threads << "\n";

    vector<double> dt(N_p);
    vector<double> profile(n_slices, 0.0);

    default_random_engine generator;
    uniform_real_distribution<double> d(0.0, 1.0);

    for (int i = 0; i < N_p; i++) {
        dt[i] = 1e-6 * d(generator);
    }

    double cut_left = 1.01 * (*min_element(dt.begin(), dt.end()));
    double cut_right = 0.99 * (*max_element(dt.begin(), dt.end()));

    double sum = 0.0;
    chrono::time_point<chrono::high_resolution_clock>start;
    chrono::duration<double> elapsed_time(0.0);
    start = chrono::system_clock::now();

    for (int i = 0; i < N_t; i++) {
        histogram(dt.data(), profile.data(),
                         cut_left, cut_right, n_slices, N_p);
    }

    elapsed_time = chrono::system_clock::now() - start;
    double other_d = elapsed_time.count() - set_zero_d.count() - histo_calc_d.count() - free_mem_d.count();
    sum = mymath::sum(profile);

    double throughput = N_t * N_p / elapsed_time.count() / 1e6;
    cout << "Sum: " << sum << "\n";
    cout << prog_name << "\n";
    cout << "Elapsed time:\t\t" << elapsed_time.count() << " sec\t" << 100 * elapsed_time.count() / elapsed_time.count() << "%\n";
    cout << "init:\t" << set_zero_d.count() << " sec\t" << 100 * set_zero_d.count() / elapsed_time.count() << "%\n";
    cout << "calc:\t" << histo_calc_d.count() << " sec\t" << 100 * histo_calc_d.count() / elapsed_time.count() << "%\n";
    cout << "free:\t" << free_mem_d.count() << " sec\t" << 100 * free_mem_d.count() / elapsed_time.count() << "%\n";
    cout << "other:\t" << other_d << " sec\t" << 100 * other_d / elapsed_time.count() << "%\n";
    cout << "Throughput: " << throughput << " MP/sec\n";
    cout << "Throughput/thread: " << throughput/N_threads << " MP/sec\n";

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
        N_SLICES,
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
            N_SLICES, 0, "s", "slices", util::Arg::Numeric,
            "  --slices=<num>,       -s <num>  Number of slices (default: 100)"
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
        case N_SLICES:
            n_slices = atoi(opt.arg);
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
