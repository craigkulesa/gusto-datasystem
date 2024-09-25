#include "corrspec.h"
#include "callback.h"

#define handle_error(msg) \
	do { perror(msg); exit(EXIT_FAILURE); } while (0)

extern struct coeffs c;

// load up the QC coefficients from file and store as struct coeffs c
int readQCfile()
{
  FILE *fp;
  int idx=0;
  fp = fopen("coeffs.txt", "r");
  while(fscanf(fp,"%lf %d %d", &c.a[idx], &c.i[idx], &c.j[idx]) != EOF)
      idx++;
  fclose(fp);
  return idx;
}

struct Spectrum spec[4];

static void handler(int sig, siginfo_t *si, void *unused)
{
   printf("Got SIGSEGV at address: 0x%lx\n", (long) si->si_addr);
   printf("Gracefully exiting\n");
   sleep(3);
}

int main(int argc, char **argv) {

   // Set up SIGSEGV handler
   FILE *fp;
   char fileName[128];
   struct sigaction sa;
   sa.sa_flags = SA_SIGINFO;
   sigemptyset(&sa.sa_mask);
   sa.sa_sigaction = handler;
   if (sigaction(SIGSEGV, &sa, NULL) == -1)
	   handle_error("sigaction");

   c.len = readQCfile();

   // Setup all possible FFTW array lengths
   printf("readying fft\n");
   for(int i=0; i<4; i++){
     int N=(i+1)*128;
     spec[i].in  = (fftw_complex *) fftw_malloc((4*N) *  sizeof(fftw_complex));
     spec[i].out = (fftw_complex *) fftw_malloc((4*N) *  sizeof(fftw_complex));
     spec[i].p = fftw_plan_dft_1d((4*N-1), spec[i].in, spec[i].out, FFTW_FORWARD, FFTW_PATIENT|FFTW_PRESERVE_INPUT);
   }
   fftw_import_system_wisdom();
   printf("ready to start\n");

   fp = fopen(argv[1], "r");
   while(fscanf(fp,"%s", fileName) != EOF)
      callback(fileName);
   fclose(fp);
   return 0;
}

