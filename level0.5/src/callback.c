#include "callback.h"
#include "corrspec.h"
#include "influx.h"

#define PI 3.14159
#define DEBUG 0
#define REF 0.025
#define CIIFREQ 1900536.9
#define NIIFREQ 1461132.0

struct coeffs c;
int CALID = -1, lastBand = 0, lastScanID = 0;
float dacV[4][4]; // DEV and then 0=VIhi, 1=VQhi, 2=VIlo, 3=VQlo

// perform the polynomial fit to the QC correction
// using the coefficients loaded from coeffs.txt at runtime
float polyfit(double x, double y)
{
  double z = 0.0;
  int idx;
  // These limits ensure that the returned vals are ceiling'd to 1.0 or floored to -1.0
  // even if the polyfit diverges outside of the fit region.  For safety; derived from omnisys dll.
  float yMax[50]={0.00,0.02,0.03,0.04,0.06,0.07,0.08,0.09,0.11,0.12,0.14,0.15,0.17,0.19,0.20,0.22,0.24,0.26,0.28,0.30,0.32,0.34,0.37,0.39,0.41,0.43,0.45,0.47,0.49,0.51,0.53,0.55,0.57,0.59,0.60,0.62,0.64,0.65,0.66,0.68,0.69,0.70,0.71,0.72,0.73,0.73,0.74,0.74,0.75,0.75};
  for(idx = 0; idx < c.len+1; idx++)
      z += c.a[idx] * pow(x, (double)c.i[idx]) * pow(y, (double)c.j[idx]);
  if(z>1.0 || y>yMax[(int)(x*100)]) // x*100 turns the power level into the corresponding array index)
    z=1.0;
  if(z<-1.0 || y<-1.0*yMax[(int)(x*100)])
    z=-1.0;
  return (float)z;
}

// Single SELECT for CORRELATOR DACS from nearest previous Correlator Cal
int getDACVfromInflux(int band, char *scanIDregex) {
   int DEV, UNIT;
   CURL *curl;
   char *query = malloc(BUFSIZ);

   for(DEV=0; DEV<4; DEV++) {
      if(band==1)
         UNIT=6;
      else
         UNIT=4;
      curl = init_influx();
      sprintf(query, "&q=SELECT * FROM /^ACS%d_DEV%d_*/ WHERE \"scanID\"=~/^%s/ ORDER BY time DESC LIMIT 5", UNIT-1, DEV+1, scanIDregex);
      influxReturn = influxWorker(curl, query);
      //      Order is reverse alphabetical from InfluxDB SELECT
      dacV[DEV][0] = influxReturn->value[3]; //VIhi  
      dacV[DEV][1] = influxReturn->value[1]; //VQhi
      dacV[DEV][2] = influxReturn->value[2]; //VIlo
      dacV[DEV][3] = influxReturn->value[0]; //VQlo
      CALID = influxReturn->scanID;
      // Free Influx struct from ACS_DEV_VDAC
      freeinfluxStruct(influxReturn);
   }
   free(query);
   return(CALID);
}

void get_proctime(char *proctime) {
   char command[BUFSIZ];
   FILE *fp;
   char date[BUFSIZ];
   
   snprintf(command, sizeof(command), "date +%%Y%%m%%d_%%H%%M%%S");
   fp = popen(command, "r");
   if (fp == NULL) {
      perror("popen");
      exit(EXIT_FAILURE);
   }
   fgets(date, sizeof(date), fp);
   pclose(fp);
   date[strcspn(date, "\n")] = 0; // Remove trailing newline

   snprintf(proctime, BUFSIZ, "%s", date);
}


void append_to_fits_table(const char *filename, struct s_header *fits_header, double *array) {
    fitsfile *fptr;  // FITS file pointer
    int status = 0;  // CFITSIO status value MUST be initialized to zero!
    int array_length, band, npix, seqflag = 0;
    long nrows;
    char extname[] = "DATA_TABLE", proctime[256], line[4];
    float linefreq;
    
    // Try to open the FITS file in read/write mode. If it doesn't exist, create a new one.
    if (fits_open_file(&fptr, filename, READWRITE, &status)) {
        if (status == FILE_NOT_OPENED) {
            status = 0;  // Reset the status
            if (fits_create_file(&fptr, filename, &status)) {
                fits_report_error(stderr, status);  // Print any error message
                return;
            }

            // Create a primary array image (needed before any extensions can be created)
            if (fits_create_img(fptr, BYTE_IMG, 0, NULL, &status)) {
                fits_report_error(stderr, status);  // Print any error message
                return;
            }

            // Construct the primary FITS HEADER
	    // Various indices and keywords in the primary header that depend on Band #
            if (fits_header->unit == 6){ //ACS5 B1
	       strcpy(line, "NII");
	       band = 1;
	       npix = 512;
	       linefreq = NIIFREQ;
            }
            if (fits_header->unit == 4){ //ACS3 B2
	       strcpy(line, "CII");
	       band = 2;
	       npix = 1024;
	       linefreq = CIIFREQ;
            }

	    // Create some Primary header keyword value pairs and fill them from the current fits_header struct
            fits_write_key(fptr, TINT,      "CALID",   &fits_header->CALID,  "ID of correlator calibration", &status);
            fits_write_key(fptr, TSTRING, "TELESCOP",  "GUSTO",   "Observatory Name", &status);
            fits_write_key(fptr, TSTRING, "LINE",      &line,     "Line Name", &status);
            fits_write_key(fptr, TFLOAT,  "LINEFREQ",  &linefreq, "Line freq in GHz", &status);
            fits_write_key(fptr, TINT,    "BAND",      &band,     "GUSTO band #",     &status);
            fits_write_key(fptr, TINT,    "NPIX",      &npix,     "N spec pixels",    &status);
            fits_write_key(fptr, TSTRING, "DLEVEL",    "0.5",      "data level",      &status);
            get_proctime(proctime);
            fits_write_key(fptr, TSTRING, "Proctime",  proctime,  "processing time",  &status);
            fits_write_key(fptr, TINT,    "SEQ_FLAG",   &seqflag,       "SEQUENCE FLAG",    &status);


            // Define the column parameters
            char *ttype[] ={"MIXER", "NINT", "UNIXTIME", "NBYTES", "CORRTIME", "INTTIME", "ROW_FLAG", "Ihigh", \
			    "Qhigh", "Ilow", "Qlow", "Ierr", "Qerr", "VIhi", "VQhi", "VIlo", "VQlo", "scanID", \
		            "subScan", "scan_type", "RA", "DEC", "filename", "PSat", "Vmon", "Imon", "Gmon", \
			    "DATA", "CHANNEL_FLAG"};

            char *tunit[] ={" ", " ", "sec", " ", " ", "sec", " ", "counts", "counts", "counts",   \
			    "counts", " ", " ", "Volts", "Volts", "Volts", "Volts", " ", " ", " ", \
		            "degrees", "degrees", "text", "Amps", "mV", "uA", "Amps", " ", " "};

	    // All header values are signed 32-bit except UNIXTIME which is uint64_t
            char *tform[29];
	    tform[0]  = "1J"; //int     mixer
	    tform[1]  = "1J"; //int	nint
	    tform[2]  = "1D"; //double	64 bit unixtime+fractional
	    tform[3]  = "1J"; //int	nbytes
	    tform[4]  = "1J"; //int     corrtime
	    tform[5]  = "1E"; //float   integration time
            tform[6]  = "1J"; //32 bit row flag
            tform[7]  = "1J"; //int	ihi
            tform[8]  = "1J"; //int	qhi
            tform[9]  = "1J"; //int	ilo
            tform[10] = "1J"; //int	qlo
            tform[11] = "1J"; //int	ierr
            tform[12] = "1J"; //int	qerr
            tform[13] = "1E"; //float	Vdac
            tform[14] = "1E"; //float	Vdac
            tform[15] = "1E"; //float	Vdac
            tform[16] = "1E"; //float	Vdac
            tform[17] = "1J"; //int	scanID
            tform[18] = "1J"; //int	subScan
            tform[19] = "6A"; //char    scan type
            tform[20] = "1E"; //float	RA
            tform[21] = "1E"; //float	DEC
	    tform[22] = "48A";//char    filename
            tform[23] = "1E"; //float	PSat
            tform[24] = "1E"; //float	Vmon
	    tform[25] = "1E"; //float	Imon
	    tform[26] = "1E"; //float   Gmon
	    // Various indices and keywords in the per row header that depend on Band #
            if (fits_header->unit==6){ //ACS5 B1
	       tform[27] = "512E";  //32 bit float
	       tform[28] = "512I";  //16 bit int
	    }
            if (fits_header->unit==4){ //ACS3 B2
	       tform[27] = "1024E";  //32 bit float
	       tform[28] = "1024I";  //16 bit int
	    }

	    int tfields = 29;
            // Create a binary table
            if (fits_create_tbl(fptr, BINARY_TBL, 0, tfields ,ttype, tform, tunit, extname, &status)) {
                fits_report_error(stderr, status);  // Print any error message
                return;
            }

        } else {
            fits_report_error(stderr, status);  // Print any error message
            return;
        }
    } ////////// END PRIMARY HEADER SECTION //////////

    // need to do this *again* because the one several lines back is only done for the fits file creation.
    // that bug took a day to find.
    if (fits_header->unit==6){ //ACS5 B1
        array_length=512;
    }
    if (fits_header->unit==4){ //ACS3 B2
        array_length=1024;
    }

    // Move to the named HDU (where the table is stored)
    if (fits_movnam_hdu(fptr, BINARY_TBL, extname, 0, &status)) {
        fits_report_error(stderr, status);  // Print any error message
        return;
    }

    // Get the current number of rows in the table
    if (fits_get_num_rows(fptr, &nrows, &status)) {
        fits_report_error(stderr, status);  // Print any error message
        return;
    }

    // insert a single empty row at the end of the output table
    if (fits_insert_rows(fptr, nrows, 1, &status)) {
        fits_report_error(stderr, status);  // Print any error message
        return;
    }

    // Write the header data
    fits_write_col(fptr, TINT32BIT,  1, nrows+1, 1, 1, &fits_header->mixer,    &status);
    fits_write_col(fptr, TINT32BIT,  2, nrows+1, 1, 1, &fits_header->nint,     &status);
    fits_write_col(fptr, TDOUBLE,    3, nrows+1, 1, 1, &fits_header->fulltime, &status);
    fits_write_col(fptr, TINT32BIT,  4, nrows+1, 1, 1, &fits_header->nbytes,   &status);
    fits_write_col(fptr, TINT32BIT,  5, nrows+1, 1, 1, &fits_header->corrtime, &status);
    fits_write_col(fptr, TFLOAT,     6, nrows+1, 1, 1, &fits_header->inttime,  &status);
    fits_write_col(fptr, TINT32BIT,  7, nrows+1, 1, 1, &fits_header->row_flag, &status);
    fits_write_col(fptr, TINT32BIT,  8, nrows+1, 1, 1, &fits_header->Ihi,      &status);
    fits_write_col(fptr, TINT32BIT,  9, nrows+1, 1, 1, &fits_header->Qhi,      &status);
    fits_write_col(fptr, TINT32BIT, 10, nrows+1, 1, 1, &fits_header->Ilo,      &status);
    fits_write_col(fptr, TINT32BIT, 11, nrows+1, 1, 1, &fits_header->Qlo,      &status);
    fits_write_col(fptr, TINT32BIT, 12, nrows+1, 1, 1, &fits_header->Ierr,     &status);
    fits_write_col(fptr, TINT32BIT, 13, nrows+1, 1, 1, &fits_header->Qerr,     &status);
    fits_write_col(fptr, TFLOAT,    14, nrows+1, 1, 1, &fits_header->VIhi,     &status);
    fits_write_col(fptr, TFLOAT,    15, nrows+1, 1, 1, &fits_header->VQhi,     &status);
    fits_write_col(fptr, TFLOAT,    16, nrows+1, 1, 1, &fits_header->VIlo,     &status);
    fits_write_col(fptr, TFLOAT,    17, nrows+1, 1, 1, &fits_header->VQlo,     &status);
    fits_write_col(fptr, TINT32BIT, 18, nrows+1, 1, 1, &fits_header->scanID,   &status);
    fits_write_col(fptr, TINT32BIT, 19, nrows+1, 1, 1, &fits_header->subScan,  &status);
    fits_write_col(fptr, TSTRING,   20, nrows+1, 1, 1, &fits_header->type,     &status);
    fits_write_col(fptr, TSTRING,   23, nrows+1, 1, 1, &fits_header->filename, &status);

    // Write the spectra as a single 2*N column
    if (fits_write_col(fptr, TDOUBLE, 28, nrows+1, 1, 1 * array_length, array, &status)) {
        fits_report_error(stderr, status);  // Print any error message
        return;
    }
    // Write the channel flag as a single 2*N column
    if (fits_write_col(fptr, TDOUBLE, 29, nrows+1, 1, 1 * array_length, array, &status)) {
        fits_report_error(stderr, status);  // Print any error message
        return;
    }

    // Close the FITS file
    if (fits_close_file(fptr, &status)) {
        fits_report_error(stderr, status);  // Print any error message
        return;
    }

    if (DEBUG)
        printf("Array appended as a new row in the FITS table successfully.\n");
}

void printDateTimeFromEpoch(time_t ts)
{
   struct tm *tm = gmtime(&ts);

   char buffer[26];
   strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S", tm);
   printf("UTC Date and Time: %s\n", buffer);
}

// Callback function to process the file
void callback(char *filein){
   char *fullpath= malloc(128*sizeof(char));
   strcpy(fullpath, filein); // make a copy leaving filein intact for later tokenization

   char *datafile;	// datafile is filename with no path - used in fits header
   datafile = strrchr(fullpath, '/');
   if (datafile != NULL) {
	   datafile++;
   } else {
	   datafile = fullpath;
   }

   //timing
   struct timeval begin, end;
   gettimeofday(&begin, 0);

   //correlator file objects
   int N = 0;
   FILE *fp;
   double P_I=0;
   double P_Q=0;
   struct corrType corr;

   // correlator variables from datafile
   uint64_t UNIXTIME=0;
   int NINT=0, UNIT, DEV, NBYTES, MIXER;
   float FRAC=0, FS_FREQ; // Full scale frequency, B1==5000MHz || B2==5000MHz
   float VIhi = 0.0, VQhi = 0.0, VIlo = 0.0, VQlo = 0.0;

   // For normalized float correlator lags from Quant Corr
   float *Rn, *Rn2;

   // file open notification
   printf("opened file: %s\n", filein);
   fp = fopen(filein, "r");

   // tokenize scanID from filename
   char *token, *prefix=malloc(8*sizeof(char));
   int position = 0, band = -1, scanID = -1, subScan = -1;
   bool error  = FALSE;	// status of error which may end processing early

   // Find file type from filename
   int i=0;
   char *ptr = NULL;
   const char *prefix_names[]={"SRC", "REF", "OTF", "HOT", "COLD", "FOC", "UNK"};

   // Use strtok to tokenize the filename using underscores as delimiters
   token = strtok(filein, "_");

   // Iterate through the tokens until reaching the 2nd position
   while (token != NULL ) {

      if (position == 0 ) {      //get band
	 if (strstr(token, "ACS3"))
            band = 2;
	 if (strstr(token, "ACS5"))
            band = 1;
      }
      
      if (position == 1 ) {      //get scan type
         while ( ptr==NULL ){
            ptr = strstr(token, prefix_names[i]);
            i++;
         }
         int len = strlen(prefix_names[i-1]);
         strncpy(prefix, ptr, len);
         prefix[len] = '\0';         
      }

      if (position == 2 ) 
         if (atoi(token)>0) scanID = atoi(token);

      if (position == 3 ) 
         subScan = atoi(token);

      token = strtok(NULL, "_");
      position++;
   }

   if(DEBUG) {
	printf("Band is: %d\n", band);
	printf("The type is %s\n", prefix);
	printf("The scanID is: %d\n", scanID);
	printf("The data file index # is: %05d\n", subScan);
   }
   // Build a regex with the range of the previous 16 scanID #s for Correlator DACs
   char scanIDregex[512];
   int pos = 0;
   pos += sprintf(&scanIDregex[pos], "^(");
   for (int k=0; k<15; k++){
      pos += sprintf(&scanIDregex[pos], "%d|", scanID-k);
   }
   pos += sprintf(&scanIDregex[pos], "%d)$", scanID-15);
   if (DEBUG)
      printf("%s\n", scanIDregex);

   // integer variables for fread from datafile
   int32_t value;
   uint32_t value1;	// for first 32 bits of UNIXTIMRE
   uint32_t value2;	//

   // figure out how many spectra in the file
   fseek(fp, 24, SEEK_SET);     // go to header position
   fread(&value, 4, 1, fp);
   int32_t bps = value;         // get bytes per spectrum
   fseek(fp, 0L, SEEK_END);
   int sz = ftell(fp);          // get total bytes in file
   fseek(fp, 0L, SEEK_SET);     // go back to beginning of file
   if(DEBUG)
     printf("File has %.1f spectra\n", (float)sz/bps);

   int32_t header[22];

   corr.corrtime=0;
//////////////////////////////  LOOP OVER ALL SPECTRA IN FILE  ///////////////////////////////////

   // Start at beginning of data file
   for (int j=0; j<(int)sz/bps; j++)
   {
   if (DEBUG)
      printf("The type is %s\n", prefix);
      // Loop over header location
   for (int i=0; i<22; i++){
     if (i==3){
       //UNIXTIME is 64 bits
       fread(&value1, 4, 1, fp); //Least significant 32 bits
       fread(&value2, 4, 1, fp); //Most significant 32 bis
       UNIXTIME = (((uint64_t)value2 << 32) | value1 ) / 1000.; //unixtime is to msec, store as 1sec int
       FRAC     = (((uint64_t)value2 << 32) | value1 ) % 1000 ; //fractional part is 1msec
     }
     else
       fread(&value, 4, 1, fp);
     
     header[i] = (value);
   }

      // fill variables from header array
      UNIT          = header[0];
      DEV           = header[1];
      NINT          = header[2];
      NBYTES        = header[5];
      corr.corrtime = header[6];
      corr.Ihi      = header[8];
      corr.Qhi      = header[9];
      corr.Ilo      = header[10];
      corr.Qlo      = header[11];
      corr.Ierr     = header[12];
      corr.Qerr     = header[13];

      // Since we know UNIT and NLAGS, set correlator frequency and fits spectrum for later use
      // And also secret decoder ring (unit,dev) -> (Mixer #)
      int mixerTable[2][4] = {
	      {2, 3, 4, 6},	// band 1 mixers
	      {2, 3, 5, 8}	// band 2 mixers
      };
      MIXER = mixerTable[band-1][DEV-1];
      FS_FREQ = 5000.;

      // Indication that this is a broken header file from correlator STOP signal
      if (corr.Ierr!=0 || corr.Qerr!=0 || \
           corr.Ihi==0 ||  corr.Qhi==0 || corr.Ilo==0 || corr.Qlo==0 || \
           (corr.corrtime*256.)/(FS_FREQ*1000000.)<0.1 || (corr.corrtime*256.)/(FS_FREQ*1000000.)>10.0 )
      {
         printf("######################## ERROR ###########################\n");
         printf("#                Error, data is no good!                 #\n");
         printf("#                        Exiting!                        #\n");
         printf("######################## ERROR ###########################\n");
         printf("CORRTIME was %.6f\n", (corr.corrtime*256.)/(FS_FREQ*1000000.));
	 error = TRUE;
         break;
      }
      // go to influx now that we are more sure the data are OK
      if(j == 0 && (band != lastBand || scanID != lastScanID)) {  // but only do it once per band or scanID
	if(DEBUG)
	  printf("ScanID or band changed, doing Influxdb lookup\n");
	CALID = getDACVfromInflux(band, scanIDregex);
	lastScanID = scanID;  // we now have data for this scanID 
	lastBand = band;      // and band
      }
	// just copy from vector into floats
      VIhi = dacV[DEV-1][0];
      VQhi = dacV[DEV-1][1];
      VIlo = dacV[DEV-1][2];
      VQlo = dacV[DEV-1][3];

	// this section unfuck-ifys special cases when ICE was off by one
      if (VQlo==0.){
	VIhi=VIhi-(VIlo-VQhi);  //make up this lost data, it'l be close enough
	VQhi = dacV[DEV-1][0];
	VIlo = dacV[DEV-1][1];
	VQlo = dacV[DEV-1][2];
      }
      
      if (VIhi==0.){ //Still no values?  bail and don't make spectra
	printf("######################## ERROR ###########################\n");
	printf("#                  Error, no DAC values!                 #\n");
	printf("#                        Exiting!                        #\n");
	printf("######################## ERROR ###########################\n");
	break;
      }
      if (DEBUG)
	printf("VIhi %.3f\tVQhi %.3f\tVIlo %.3f\tVQlo %.3f\n", VIhi, VQhi, VIlo, VQlo);
      
      if (NBYTES==8256)
         N = 512;
      else if (NBYTES==6208)
         N = 384;
      else if (NBYTES==4160)
         N = 256;
      else if (NBYTES==2112)
         N = 128;
      int specA = (int) N/128 - 1;

      //We don't know the lag # until we open the file, so malloc now
      corr.II   = malloc(N*sizeof(int32_t));   //Uncorrected ints
      corr.QI   = malloc(N*sizeof(int32_t));
      corr.IQ   = malloc(N*sizeof(int32_t));
      corr.QQ   = malloc(N*sizeof(int32_t));
      corr.IIqc = malloc(N*sizeof(float));     //Normalized Quantization Corrected floats
      corr.QIqc = malloc(N*sizeof(float));
      corr.IQqc = malloc(N*sizeof(float));
      corr.QQqc = malloc(N*sizeof(float));
      Rn  = malloc(2*N*sizeof(float));         //Rn,Rn2 floats for normalized, quantization corrected
      Rn2 = malloc(4*N*sizeof(float));

      // Read lags in from file in order after header
      for (int i=0; i<N; i++){
         fread(&value, 4, 1, fp);
         corr.II[i] = value;
         fread(&value, 4, 1, fp);
         corr.QI[i] = value;
         fread(&value, 4, 1, fp);
         corr.IQ[i] = value;
         fread(&value, 4, 1, fp);
         corr.QQ[i] = value;
      }

      double IIlevel = 0.5*(((double)corr.Ilo/(double)corr.corrtime)+(double)corr.Ihi/(double)corr.corrtime);
      double QQlevel = 0.5*(((double)corr.Qlo/(double)corr.corrtime)+(double)corr.Qhi/(double)corr.corrtime);
      double IQlevel = 0.25*(((double)corr.Ilo/(double)corr.corrtime)+(double)corr.Ihi/(double)corr.corrtime+(double)corr.Qhi/(double)corr.corrtime+(double)corr.Qlo/(double)corr.corrtime);
      double Idiff = fabs((double)corr.Ilo/(double)corr.corrtime - (double)corr.Ihi/(double)corr.corrtime) ;
      double Qdiff = fabs((double)corr.Qlo/(double)corr.corrtime - (double)corr.Qhi/(double)corr.corrtime) ;
      float asymCorr = (Idiff/REF * Qdiff/REF)*0.0015;

      for (int i=0; i<N; i++){
	 corr.IIqc[i] = polyfit(IIlevel, (double)corr.II[i]/(double)corr.corrtime) - asymCorr;
         corr.IQqc[i] = polyfit(IQlevel, (double)corr.IQ[i]/(double)corr.corrtime) - asymCorr;
         corr.QIqc[i] = polyfit(IQlevel, (double)corr.QI[i]/(double)corr.corrtime) - asymCorr;
         corr.QQqc[i] = polyfit(QQlevel, (double)corr.QQ[i]/(double)corr.corrtime) - asymCorr;
      }

      // GET THE QUANTIZATION CORRECTED LAGS BACK AND  CONTINUE ON WITH LAGS->SPECTRA
      // Combine IQ lags into R[]
         for (int i=0; i<(2*N)-1; i++){
            if (i%2 == 0) Rn[i] = 0.5* (corr.IIqc[i/2] + corr.QQqc[i/2]);           // Even indicies
            if (i%2 == 1) Rn[i] = 0.5* (corr.IQqc[(i-1)/2] + corr.QIqc[(i-1)/2+1]); // Odd  indicies
         }
         int le = 2*N-1; //last element
         Rn[le] = 0.5* (corr.IQqc[(le-1)/2] + corr.QIqc[(le-1)/2]);         // Last element

      // Mirror R[] symmetrically
         for (int i=0; i<4*N-1; i++){
           if ( i<2*N ) Rn2[2*N-1-i] = Rn[i];
           else        Rn2[i]       = Rn[i-(2*N-1)];
         }

      // Apply Hann window
         for (int i=0; i<4*N; i++){
           Rn2[i] = Rn2[i]*0.5*(1-cos(2*PI*i/(4*N)));     //real with Hann window
         }

      // Fill fft struct
         for (int i=0; i<4*N; i++){
           spec[specA].in[i][0] = Rn2[i];  //real
           spec[specA].in[i][1] = 0.;      //imag
         }

      // Do FFT and print
         fftw_execute(spec[specA].p);

      // Compute Power Coefficients
         P_I = pow((VIhi-VIlo),2) * 1.8197 / pow((erfinv(1-2*(double)corr.Ihi/(double)corr.corrtime)) + \
                                                 (erfinv(1-2*(double)corr.Ilo/(double)corr.corrtime)),2);
         P_Q = pow((VQhi-VQlo),2) * 1.8197 / pow((erfinv(1-2*(double)corr.Qhi/(double)corr.corrtime)) + \
                                                 (erfinv(1-2*(double)corr.Qlo/(double)corr.corrtime)),2);

      // Construct the per row FITS HEADER
      s_header *fits_header = (s_header *)malloc(sizeof(s_header));
      fits_header->type     = (char *)malloc(8);
      fits_header->filename = (char *)malloc(48);

      fits_header->unit     = UNIT;
      fits_header->dev      = DEV;
      fits_header->mixer    = MIXER;
      fits_header->nint     = NINT;
      fits_header->fulltime = (double) UNIXTIME + 0.001*FRAC;
      fits_header->nbytes   = NBYTES;
      fits_header->corrtime = (int) corr.corrtime;
      fits_header->inttime  =(corr.corrtime*256.)/(FS_FREQ*1000000.);
      fits_header->row_flag = 0;
      fits_header->channel_flag = 0;

      fits_header->Ihi      = corr.Ihi;
      fits_header->Qhi      = corr.Qhi;
      fits_header->Ilo      = corr.Ilo;
      fits_header->Qlo      = corr.Qlo;
      fits_header->Ierr     = corr.Ierr;
      fits_header->Qerr     = corr.Qerr;

      fits_header->VIhi     = VIhi;
      fits_header->VQhi     = VQhi;
      fits_header->VIlo     = VIlo;
      fits_header->VQlo     = VQlo;

      fits_header->scanID   = scanID;
      fits_header->subScan  = subScan;
      fits_header->CALID    = CALID;

      strcpy(fits_header->type, prefix);
      strcpy(fits_header->filename, datafile);
      
      // Construct the FITS DATA
      double *array = malloc(2*N*sizeof(double));
      for (int i=0; i<2*N; i++){
         array[i] = sqrt(P_I*P_Q) * sqrt(pow(spec[specA].out[i][0],2)+pow(spec[specA].out[i][1],2));
      }

      char fitsfile[32];
      if(UNIT == 4) // band 2
	 sprintf(fitsfile, "./build/B2/ACS%d_%s_%05d.fits", UNIT-1, prefix, scanID);
      else if (UNIT == 6) // band 1
	 sprintf(fitsfile, "./build/B1/ACS%d_%s_%05d.fits", UNIT-1, prefix, scanID);
      if (DEBUG)
         printf("%s\n", fitsfile);
      append_to_fits_table(fitsfile, fits_header, array); 

      // Free items before next spectra within this file
      // All of these objects are malloced at the start of every spectrum
      free(corr.II);      //free all mallocs
      free(corr.QI);
      free(corr.IQ);
      free(corr.QQ);
      free(corr.IIqc);    //free all mallocs
      free(corr.QIqc);
      free(corr.IQqc);
      free(corr.QQqc);
      free(Rn);
      free(Rn2);
      free(array);

      free(fits_header->type);
      free(fits_header->filename);
      free(fits_header);
   } // END SPECTRA LOOP
/////////////////////////////  END LOOP OVER ALL SPECTRA IN FILE  ///////////////////////////////////
   
   fclose(fp);   //close input data file		    
   if (!error && DEBUG){
      //ouput stats for last spectra in file
      printf("UNIXTIME is %" PRIu64 "\n", UNIXTIME);
      printf("CORRTIME is %.6f\n", (corr.corrtime*256.)/(FS_FREQ*1000000.));
      printDateTimeFromEpoch((long long) UNIXTIME);
      printf("UNIT is %d\n", UNIT);
      printf("DEV  is %d\n",  DEV); 
      printf("NINT is %d\n", NINT);
      printf("nlags=%d\n", N);
      printf("etaQ\t%.3f\n", 1/sqrt(P_I*P_Q));
      printf("The cal    is from scanID: %d\n", CALID);
      printf("The data   is from scanID: %d\n", scanID);
   }

   gettimeofday(&end, 0);
   int seconds = end.tv_sec - begin.tv_sec;
   int microseconds = end.tv_usec - begin.tv_usec;
   double elapsed = seconds + microseconds*1e-6;
   if (!error)
      printf("AVG FFTW %.1f ms in %d spectra\n", 1000.*elapsed/(sz/bps), sz/bps);
   fflush(stdout);

   free(prefix);
   free(fullpath);
}

