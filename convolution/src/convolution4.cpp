#include <iostream>
#include <chrono>
#include <vector>
#include "convolution.h"
#include "configuration.h"
#include "utilities.h"
#include "optionparser.h"
#include <omp.h>


unsigned N = 100000;
unsigned ITERS = 10;
unsigned n_threads  = 1;
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
   auto sum = 0.0L;

   f_vector_t in(N);
   f_vector_t kernel(N / 2);
   double *out;
   int out_size = N + N / 2 - 1;

   for (unsigned iter = 0; iter < ITERS; ++iter) {

      for (unsigned i = 0; i < in.size(); i++) {
         //auto r2 = static_cast <double>(rand()) / (static_cast <double>(RAND_MAX / 100.0));
         in[i] = (i * 1.0) / N;
      }

      for (unsigned i = 0; i < kernel.size(); i++) {
         //auto r2 = static_cast <double>(rand()) / (static_cast <double>(RAND_MAX / 100.0));
         kernel[i] = (i * 1.0) / (10 * N);
      }

      util::get_time(start);
      out = convolution4(in.data(), in.size(),
                         kernel.data(), kernel.size());

      elapsed += util::time_elapsed(start);//end - start;
      for (int i = 0; i < out_size; ++i)
         sum += std::abs(out[i]);
   }
   std::cout << "Convolution4 of " << N << " elems\n";
   std::cout << "Elapsed Time : " << elapsed << " s\n";

   std::cerr << "Throughput : " << (N * ITERS * 8) / (elapsed * 1000000) << " MB/s\n";

   std::cout << "Sum : " << sum / ITERS << std::endl;
   std::cout << "\n\n";

   if (out) free(out);
   return 0;
}


void parse_args(int argc, char **argv)
{
   using namespace std;
   using namespace option;

   enum optionIndex {UNKNOWN, HELP, N_ELEMS, N_ITERS, N_THREADS, OPTIONS_NUM
                    };

   const option::Descriptor usage[] = {
      {
         UNKNOWN, 0, "", "", Arg::None,                  "USAGE: ./convolution4 [options]\n\n"
         "Options:"
      },
      {  HELP, 0, "h", "help", Arg::None, "--help, -h  Print usage and exit." },
      {N_ELEMS, 0, "n", "elems", util::Arg::Numeric, "--elems=<num>, -n <num>  Number of elems (default: 10k)" },
      {N_ITERS, 0, "i", "iters", util::Arg::Numeric, "--iters=<num>, -i <num>  Number of iterations (default: 10k)" },
      {N_THREADS, 0, "t", "threads", util::Arg::Numeric, "--threads=<num>, -t <num>  Number of Threads (default: 1)" },

      {
         UNKNOWN, 0, "", "", Arg::None, "\nExamples:\n"
         "\t./convolution4\n"
         "\t./convolution4 -n 1000\n"
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
         case UNKNOWN:
            // not possible because Arg::Unknown returns ARG_ILLEGAL
            // which aborts the parse with an error
            break;
      }
   }


}

