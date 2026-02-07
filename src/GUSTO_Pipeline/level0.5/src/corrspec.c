#include "corrspec.h"
#include "callback.h"
#include <sys/stat.h>

extern struct coeffs c;
char commit_info[256];

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

void get_git_commit_info(const char *filename, char *commit_info) {
    char command[BUFSIZ];
    FILE *fp;
    char hash[BUFSIZ];

    // Construct the command to get the commit hash
    snprintf(command, sizeof(command), "git log -1 --format=%%cd --date=format-local:'%%Y-%%m-%%d %%H:%%M:%%S %%Z' --pretty=format:\"Level 0.5 commit %%h by %%an %%ad\" -- %s", filename);
    fp = popen(command, "r");
    if (fp == NULL) {
        perror("popen");
        exit(EXIT_FAILURE);
    }
    fgets(hash, sizeof(hash), fp);
    pclose(fp);
    hash[strcspn(hash, "\n")] = 0; // Remove trailing newline

    // Combine hash and date into the commit_info string
    snprintf(commit_info, BUFSIZ, "%s", hash);
}


int main(int argc, char **argv) {
   // Set up SIGSEGV handler
   FILE *fp;
   char fileName[128], outputPath[128]="./build";

   if(argc < 2)
     {
       printf("Syntax:\n\t./corrspec <file list> <optional output directory>\n");
       return 1;
     }
   else if(argc == 3)
     sprintf(outputPath, "%s", argv[2]);
   
   printf("Using output directory: %s\n", outputPath);
   
   c.len = readQCfile();
   get_git_commit_info("./src/callback.c", commit_info);
   
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
     callback(fileName, outputPath);
   fclose(fp);
   return 0;
}
