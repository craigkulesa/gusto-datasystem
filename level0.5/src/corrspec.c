#include "corrspec.h"
#include "callback.h"

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

int main(int argc, char **argv) {
   // Set up SIGSEGV handler
   FILE *fp;
   char fileName[128];

   c.len = readQCfile();
   
   // Setup all possible FFTW array lengths
   for(int i=0; i<4; i++){
     int N=(i+1)*128;
     spec[i].in  = (fftw_complex *) fftw_malloc((4*N) *  sizeof(fftw_complex));
     spec[i].out = (fftw_complex *) fftw_malloc((4*N) *  sizeof(fftw_complex));
     spec[i].p = fftw_plan_dft_1d((4*N-1), spec[i].in, spec[i].out, FFTW_FORWARD, FFTW_PATIENT|FFTW_PRESERVE_INPUT);
   }
   fftw_import_system_wisdom();

   fp = fopen(argv[1], "r");
   while(fscanf(fp,"%s", fileName) != EOF)
      callback(fileName);
   fclose(fp);
   return 0;
}

