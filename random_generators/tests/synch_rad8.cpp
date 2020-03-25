
#include <iostream>
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

using namespace std;
string prog_name = "synch_rad8";
long int N_t = 1;
long int N_p = 100000;
int N_threads = 1;
long int n_kicks = 1;

// std::random_device rd;
// std::mt19937_64 gen(rd());
// std::normal_distribution<> dist(0.0, 1.0);

chrono::duration<double> synch_rad_d(0.0);
chrono::duration<double> rand_gen_d(0.0);
chrono::duration<double> quantum_d(0.0);

// This function calculates and applies only the synchrotron radiation damping term
extern "C" void synchrotron_radiation(double * __restrict__ beam_dE, const double U0, 
                            const int n_macroparticles, const double tau_z, 
                            const int n_kicks){

    // SR damping constant
    const double const_synch_rad = 2.0 / tau_z;

    for (int j=0; j<n_kicks; j++){
        // SR damping term due to energy spread
        #pragma omp parallel for
        for (int i = 0; i < n_macroparticles; i++)
            beam_dE[i] -= const_synch_rad * beam_dE[i];
    
        // Average energy change due to SR
        #pragma omp parallel for
        for (int i = 0; i < n_macroparticles; i++)
            beam_dE[i] -= U0;
    }
}


// This function calculates and applies synchrotron radiation damping and
// quantum excitation terms
extern "C" void synchrotron_radiation_full(double * __restrict__ beam_dE, const double U0,
                                        const int n_macroparticles, const double sigma_dE,
                                        const double tau_z,const double energy,
                                        double * __restrict__ random_array,
                                        const int n_kicks){
    
    chrono::time_point<chrono::high_resolution_clock> start_t;
    std::hash<std::thread::id> hash;
    // Quantum excitation constant
    const double const_quantum_exc = 2.0 * sigma_dE / sqrt(tau_z) * energy;
    const double const_synch_rad = 1 - 2.0 / tau_z;
    const int STEP = 16;
    // Random number generator for the quantum excitation term
    for (int j=0; j<n_kicks; j++){
        // Compute synchrotron radiation damping term
        // start_t = chrono::system_clock::now();
        // synchrotron_radiation(beam_dE, U0, n_macroparticles, tau_z, 1);
        // synch_rad_d += chrono::system_clock::now() - start_t; 

        start_t = chrono::system_clock::now();
        // Re-calculate the random (Gaussian) number array
        #pragma omp parallel
        {
            static __thread boost::mt19937_64 *gen = nullptr;
            if(!gen) gen = new boost::mt19937_64(clock() + hash(this_thread::get_id()));    
            static __thread boost::normal_distribution<> dist(0.0, 1.0);
            
            double temp[STEP];
            #pragma omp for
            for (int i = 0; i < n_macroparticles; i+= STEP) {
                const int loop_count = n_macroparticles - i > STEP ?
                       STEP : n_macroparticles - i;
                for (int j = 0; j < loop_count; j++) {
                    temp[j] = dist(*gen);
                }
                // random_array[i] = dist(*gen);
                for (int j = 0; j < loop_count; j++) {
                    beam_dE[i+j] = beam_dE[i+j] * const_synch_rad + const_quantum_exc * temp[j] - U0;
                }
            }
        }
        rand_gen_d += chrono::system_clock::now() - start_t; 

        // start_t = chrono::system_clock::now();
        // Applies the quantum excitation term
        // #pragma omp parallel for
        // for (int i = 0; i < n_macroparticles; i++){
        // }
        // quantum_d += chrono::system_clock::now() - start_t; 
    }
}



void parse_args(int argc, char **argv);

int main(int argc, char *argv[])
{
    parse_args(argc, argv);
    omp_set_num_threads(N_threads);
    cout << "Number of turns: " <<  N_t << "\n";
    cout << "Number of points: " << N_p << "\n";
    cout << "Number of kicks: " << n_kicks << "\n";
    cout << "Number of openmp threads: " <<  N_threads << "\n";

    vector<double> dE(N_p), random_array(N_p);

    default_random_engine generator;
    uniform_real_distribution<double> d(0.0, 1.0);
    
    for (int i = 0; i < N_p; i++) {
        dE[i] = 10e6 * d(generator);
    }

    double U0 = d(generator);
    double sigma_dE = d(generator);
    double tau_z = d(generator);
    double energy = d(generator);

    double sum = 0.0;
    chrono::time_point<chrono::high_resolution_clock>start;
    chrono::duration<double> elapsed_time(0.0);
    start = chrono::system_clock::now();

    for (int i = 0; i < N_t; i++) {
        synchrotron_radiation_full(dE.data(), U0, N_p, sigma_dE, 
             tau_z, energy, random_array.data(), n_kicks);
    }

    // elapsed_time = chrono::system_clock::now() - start;
    elapsed_time = synch_rad_d + rand_gen_d + quantum_d;
    sum = mymath::sum(dE);

    double throughput = N_t * N_p * n_kicks / elapsed_time.count() / 1e6;
    cout << "Sum: " << sum << "\n";
    cout << prog_name << "\n";
    cout << "Elapsed time:\t\t" << elapsed_time.count() << " sec\t" << 100 * elapsed_time.count() / elapsed_time.count() << "%\n";
    cout << "synchrotron_radiation:\t" << synch_rad_d.count() << " sec\t" << 100 * synch_rad_d.count() / elapsed_time.count() << "%\n";
    cout << "random_generation:\t" << rand_gen_d.count() << " sec\t" << 100 * rand_gen_d.count() / elapsed_time.count() << "%\n";
    cout << "quantum_excitation:\t" << quantum_d.count() << " sec\t" << 100 * quantum_d.count() / elapsed_time.count() << "%\n";
    // cout << "Throughput: " << throughput << " MP/sec\n";
    // cout << "Throughput/thread: " << throughput/N_threads << " MP/sec\n";

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
