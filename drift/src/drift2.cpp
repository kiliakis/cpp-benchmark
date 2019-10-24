/*
Copyright 2016 CERN. This software is distributed under the
terms of the GNU General Public Licence version 3 (GPL Version 3),
copied verbatim in the file LICENCE.md.
In applying this licence, CERN does not waive the privileges and immunities
granted to it by virtue of its status as an Intergovernmental Organization or
submit itself to any jurisdiction.
Project website: http://blond.web.cern.ch/
*/
// Optimised C++ routine that calculates the drift.
// Author: Danilo Quartullo, Helga Timko, Alexandre Lasheen

#include <iostream>
#include <blond/vector_math.h>
#include <blond/utilities.h>
#include <blond/math_functions.h>
#include <blond/optionparser.h>
#include <random>
#include <chrono>
#include <string.h>

using namespace std;

long int N_t = 1;
long int N_p = 100000;
int N_threads = 1;


void drift(double * __restrict__ beam_dt,
           const double * __restrict__ beam_dE,
           const char * __restrict__ solver,
           const double T0, const double length_ratio,
           const double alpha_order, const double eta_zero,
           const double eta_one, const double eta_two,
           const double alpha_zero, const double alpha_one,
           const double alpha_two,
           const double beta, const double energy,
           const int n_macroparticles) {

    int i;
    double T = T0 * length_ratio;

    if ( strcmp (solver, "simple") == 0 )
    {
        double coeff = eta_zero / (beta * beta * energy);
        #pragma omp parallel for
        for (int i = 0; i < n_macroparticles; i++)
            beam_dt[i] += T * coeff * beam_dE[i];
    }

    else if ( strcmp (solver, "legacy") == 0 )
    {
        const double coeff = 1. / (beta * beta * energy);
        const double eta0 = eta_zero * coeff;
        const double eta1 = eta_one * coeff * coeff;
        const double eta2 = eta_two * coeff * coeff * coeff;

        if ( alpha_order == 1 )
            for ( i = 0; i < n_macroparticles; i++ )
                beam_dt[i] += T * (1. / (1. - eta0 * beam_dE[i]) - 1.);
        else if (alpha_order == 2)
            for ( i = 0; i < n_macroparticles; i++ )
                beam_dt[i] += T * (1. / (1. - eta0 * beam_dE[i]
                                         - eta1 * beam_dE[i] * beam_dE[i]) - 1.);
        else
            for ( i = 0; i < n_macroparticles; i++ )
                beam_dt[i] += T * (1. / (1. - eta0 * beam_dE[i]
                                         - eta1 * beam_dE[i] * beam_dE[i]
                                         - eta2 * beam_dE[i] * beam_dE[i] * beam_dE[i]) - 1.);
    }

    else
    {

        const double invbetasq = 1. / (beta * beta);
        const double invenesq = 1. / (energy * energy);
        double beam_delta;

        #pragma omp parallel for
        for ( i = 0; i < n_macroparticles; i++ )
        {

            beam_delta = sqrt(1. + invbetasq *
                              (beam_dE[i] * beam_dE[i] * invenesq + 2.*beam_dE[i] / energy)) - 1.;

            beam_dt[i] += T * ((1. + alpha_zero * beam_delta +
                                alpha_one * (beam_delta * beam_delta) +
                                alpha_two * (beam_delta * beam_delta * beam_delta)) *
                               (1. + beam_dE[i] / energy) / (1. + beam_delta) - 1.);

        }

    }

}




void parse_args(int argc, char **argv);

int main(int argc, char *argv[])
{
    parse_args(argc, argv);
    omp_set_num_threads(N_threads);
    cout << "Number of turns: " <<  N_t << "\n";
    cout << "Number of points: " << N_p << "\n";
    cout << "Number of openmp threads: " <<  N_threads << "\n";

    vector<double> dE(N_p), dt(N_p);
    // vector<double> omega_rf(N_rf), voltage(N_rf), phi_rf(N_rf);

    default_random_engine gen;
    uniform_real_distribution<double> d(0.0, 1.0);
    // apply_f_in_place(dt, [d, gen](double x) {return d(gen);});
    // apply_f_in_place(dE, [d, gen](double x) {return d(gen);});
    // apply_f_in_place(omega_rf, [d, gen](double x) {return d(gen);});
    // apply_f_in_place(voltage, [d, gen](double x) {return d(gen);});
    // apply_f_in_place(phi_rf, [d, gen](double x) {return d(gen);});

    for (int i = 0; i < N_p; i++) {
        dt[i] = 10e-6 * d(gen);
        dE[i] = 10e6 * d(gen);
    }

    const char *solver = "full";
    const double T0 = d(gen);
    const double length_ratio = 1.0;
    const double alpha_order = 3;
    const double eta_zero = d(gen);
    const double eta_one = d(gen);
    const double eta_two = d(gen);
    const double alpha_zero = d(gen);
    const double alpha_one = d(gen);
    const double alpha_two = d(gen);
    const double beta = d(gen);
    const double energy = d(gen);


    double sum = 0.0;
    chrono::time_point<chrono::high_resolution_clock>start;
    chrono::duration<double> elapsed_time(0.0);
    start = chrono::system_clock::now();

    for (int i = 0; i < N_t; i++) {
        drift(dt.data(), dE.data(), solver, T0,
              length_ratio, alpha_order, eta_zero, eta_one,
              eta_two, alpha_zero, alpha_one, alpha_two,
              beta, energy, N_p);
    }



    elapsed_time = chrono::system_clock::now() - start;
    sum = mymath::sum(dE);

    double throughput = N_t * N_p / elapsed_time.count() / 1e6;
    cout << "Drift new\n";
    cout << "Elapsed time: " << elapsed_time.count() << " sec\n";
    cout << "Throughput: " << throughput << " MP/sec\n";
    cout << "Throughput/thread: " << throughput / N_threads << " MP/sec\n";
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
        // N_RF,
        OPTIONS_NUM
    };

    const option::Descriptor usage[] = {
        {
            UNKNOWN, 0, "", "", Arg::None, "USAGE: ./drift2 [options]\n\n"
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
        // {
        //     N_RF, 0, "r", "rf_sections", util::Arg::Numeric,
        //     "  --rf_sections=<num>,       -r <num>  Number of rf_sections (default: 1)"
        // },
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
            "\t./drift2\n"
            "\t./drift2 -t 1000 -p 10000 -m 4\n"
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
        // case N_RF:
        //     N_rf = atoi(opt.arg);
        //     // fprintf(stdout, "--numeric with argument '%s'\n", opt.arg);
        //     break;
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
