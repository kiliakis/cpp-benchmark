
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
string prog_name = "fastres3";
long int N_t = 1;
int N_p = 100000;
int N_threads = 1;
int n_res = 100;
int n_freqs = 100;
chrono::duration<double> set_zero_d(0.0);
chrono::duration<double> histo_calc_d(0.0);
chrono::duration<double> free_mem_d(0.0);


extern "C" void fast_resonator_real_imag(float *__restrict__ impedanceReal,
        float *__restrict__ impedanceImag,
        const float *__restrict__ frequencies,
        const float *__restrict__ shunt_impedances,
        const float *__restrict__ Q_values,
        const float *__restrict__ resonant_frequencies,
        const int n_resonators,
        const int n_frequencies)
        
{   /*
    This function takes as an input a list of resonators parameters and 
    computes the impedance in an optimised way.
    
    Parameters
    ---------- 
    frequencies : float array
        array of frequency in Hz
    shunt_impedances : float array
        array of shunt impedances in Ohm
    Q_values : float array
        array of quality factors
    resonant_frequencies : float array
        array of resonant frequency in Hz
    n_resonators : int
        number of resonantors
    n_frequencies : int
        length of the array 'frequencies'
    
    Returns
    -------
    impedanceReal : float array
        real part of the impedance
    impedanceImag : float array
        imaginary part of the impedance
      */


    for (int res = 0; res < n_resonators; res++) {
        const float Qsquare = Q_values[res] * Q_values[res];
        #pragma omp parallel for
        for (int freq = 1; freq < n_frequencies; freq++) {
            const float commonTerm = (frequencies[freq]
                                       / resonant_frequencies[res]
                                       - resonant_frequencies[res]
                                       / frequencies[freq]);

            impedanceReal[freq] += shunt_impedances[res]
                                   / (1.0 + Qsquare * commonTerm * commonTerm);

            impedanceImag[freq] -= shunt_impedances[res]
                                   * (Q_values[res] * commonTerm)
                                   / (1.0 + Qsquare * commonTerm * commonTerm);
        }
    }

}


void parse_args(int argc, char **argv);

int main(int argc, char *argv[])
{
    parse_args(argc, argv);
    omp_set_num_threads(N_threads);
    cout << "Number of turns: " <<  N_t << "\n";
    // cout << "Number of points: " << N_p << "\n";
    cout << "Number of resonantors: " << n_res << "\n";
    cout << "Number of frequencies: " << n_freqs << "\n";
    cout << "Number of openmp threads: " <<  N_threads << "\n";

    vector<float> impedReal(n_freqs);
    vector<float> impedImag(n_freqs);
    vector<float> freqs(n_freqs);

    vector<float> Q_values(n_res);
    vector<float> shunt_impedances(n_res);
    vector<float> resonant_frequencies(n_res);

    default_random_engine generator;
    uniform_real_distribution<float> d(0.0, 1.0);

    for (int i = 0; i < n_freqs; i++) {
        impedReal[i] = d(generator);
        impedImag[i] = d(generator);
        freqs[i] = d(generator);
    }

    for (int i = 0; i < n_res; i++) {
        Q_values[i] = d(generator);
        shunt_impedances[i] = d(generator);
        resonant_frequencies[i] = d(generator);
    }

    // float sumR = 0.0;
    // float sumI = 0.0;
    chrono::time_point<chrono::high_resolution_clock>start;
    chrono::duration<double> elapsed_time(0.0);
    start = chrono::system_clock::now();

    for (int i = 0; i < N_t; i++) {
        fast_resonator_real_imag(impedReal.data(), impedImag.data(),
                                 freqs.data(), shunt_impedances.data(), Q_values.data(),
                                 resonant_frequencies.data(), n_res, n_freqs);
    }

    elapsed_time = chrono::system_clock::now() - start;
    // double other_d = elapsed_time.count() - set_zero_d.count() - histo_calc_d.count() - free_mem_d.count();
    // sumR = mymath::sum(impedReal);
    // sumI = mymath::sum(impedImag);

    double throughput = N_t * n_freqs * n_res  / elapsed_time.count() / 1e6;
    cout << "Sum: " << mymath::sum(impedReal) + mymath::sum(impedImag)
         << "\tAvg: " << mymath::mean(impedReal) + mymath::mean(impedImag) 
         << "\tSTD: " << mymath::standard_deviation(impedReal) + mymath::standard_deviation(impedImag) 
         << "\n";
    cout << prog_name << "\n";
    cout << "Elapsed time:\t\t" << elapsed_time.count() << " sec\t" << 100 * elapsed_time.count() / elapsed_time.count() << "%\n";
    // cout << "init:\t" << set_zero_d.count() << " sec\t" << 100 * set_zero_d.count() / elapsed_time.count() << "%\n";
    // cout << "calc:\t" << histo_calc_d.count() << " sec\t" << 100 * histo_calc_d.count() / elapsed_time.count() << "%\n";
    // cout << "free:\t" << free_mem_d.count() << " sec\t" << 100 * free_mem_d.count() / elapsed_time.count() << "%\n";
    // cout << "other:\t" << other_d << " sec\t" << 100 * other_d / elapsed_time.count() << "%\n";
    cout << "Throughput: " << throughput << " MP/sec\n";
    cout << "Throughput/thread: " << throughput / N_threads << " MP/sec\n";

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
        N_RES,
        N_FREQS,
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
            N_RES, 0, "r", "resonators", util::Arg::Numeric,
            "  --resonators=<num>,       -r <num>  Number of resonators (default: 100)"
        },
        {
            N_FREQS, 0, "f", "frequencies", util::Arg::Numeric,
            "  --frequencies=<num>,       -f <num>  Number of frequencies (default: 1000)"
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
        case N_RES:
            n_res = atoi(opt.arg);
            // fprintf(stdout, "--numeric with argument '%s'\n", opt.arg);
            break;
        case N_FREQS:
            n_freqs = atoi(opt.arg);
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
