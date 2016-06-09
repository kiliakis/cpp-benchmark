#include <iostream>
#include <chrono>
#include <vector>
#include <cmath>
#include "math_functions.h"
#include "configuration.h"
#include "utilities.h"
#include "optionparser.h"

unsigned N = 100000;
unsigned ITERS = 10;
void parse_args(int argc, char **argv);

int main(int argc, char **argv)
{
   parse_args(argc, argv);

   std::cout << "Num of Iterations : " << ITERS << "\n";
   std::cout << "Num of Elems/Iteration : " << N << "\n";
   std::cout << "\n\n";
   std::cout.precision(4);

   timespec start;
   auto elapsed = 0.0L;
   auto sum = 0.0L;

   // if (not fftw_init_threads()) {
   //    std::cerr << "Failed to initialize fftw multithreaded version\n";
   // }

   complex_vector_t v(N);
   complex_vector_t out(N);
   //auto p = mymath::init_fft(N, v.data(), out.data(), FFTW_FORWARD);

   for (unsigned iter = 0; iter < ITERS; ++iter) {

      for (unsigned i = 0; i < N; i++) {
         //float r2 = static_cast <float>(rand()) / (static_cast <float>(RAND_MAX / 100.0));
         v[i] = complex_t(i, i);
      }

      util::get_time(start);
      //auto p = mymath::init_fft(N, v.data(), out.data(), FFTW_FORWARD);
      //mymath::run_fft(p);
      //mymath::destroy_fft(p);
      mymath::fft(v, out, v.size());
      elapsed += util::time_elapsed(start);//end - start;

      for (const auto &z : out)
         sum += std::abs(z);

   }

//   mymath::destroy_fft(p);
   //fftw_cleanup_threads();
   std::cout << "FFT of " << N << " elems\n";
   std::cout << "Elapsed Time : " << elapsed << " s\n";
   std::cout << "Throughput : " << (N * ITERS) / (elapsed * 1000000) << " M/s\n";
   std::cout << "Sum : " << sum / ITERS << std::endl;
   std::cout << "\n\n";

   return 0;
}


void parse_args(int argc, char **argv)
{
   using namespace std;
   using namespace option;

   enum optionIndex {UNKNOWN, HELP, N_ELEMS, N_ITERS, OPTIONS_NUM
                    };

   const option::Descriptor usage[] = {
      {
         UNKNOWN, 0, "", "", Arg::None,                  "USAGE: ./fft-bench [options]\n\n"
         "Options:"
      },
      {HELP, 0, "h", "help", Arg::None, "--help, -h  Print usage and exit." },
      {N_ELEMS, 0, "n", "elems", util::Arg::Numeric, "--elems=<num>, -n <num>  Number of elems (default: 10k)" },
      {N_ITERS, 0, "i", "iters", util::Arg::Numeric, "--iters=<num>, -i <num>  Number of iterations (default: 10k)" },
      {
         UNKNOWN, 0, "", "", Arg::None, "\nExamples:\n"
         "\t./fft-bench\n"
         "\t./fft-bench -n 1000\n"
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
         case UNKNOWN:
            // not possible because Arg::Unknown returns ARG_ILLEGAL
            // which aborts the parse with an error
            break;
      }
   }


}

