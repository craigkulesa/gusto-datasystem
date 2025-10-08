#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <curl/curl.h>
#include <sys/param.h>
#include "corrspec.h"
#include "influx.h"

#define MAXLENGTH 5
#define DEBUG 0

influxStruct *influxReturn=NULL;
int numNames=0;

void extract_all_keyword_values(const char *source, const char *keyword) {
    const char *found = source;
    char result[64];
    int result_size = sizeof(result);
    int i=0;
    char ID[64]="";

    // Loop until no more keywords are found
    while ((found = strstr(found, keyword)) != NULL && i<MAXLENGTH) {
        // Move the pointer to the position after the keyword
        found += strlen(keyword);

        // Check if the next character is a quote
        if (*found == '"') {
            found++; // Move past the starting quote

            // Find the ending quote
            const char *end = strchr(found, '"');
            if (end) {
                // Calculate the length of the string inside the quotes
                size_t length = end - found;

                // Ensure we don't copy more than result_size - 1 characters
                if (length >= result_size) {
                    length = result_size - 1;
                }

                // Copy the string inside the quotes to the result buffer
                strncpy(result, found, length);
                result[length] = '\0'; // Null-terminate the result

                // Print the extracted value
	        memset(ID, '\0', sizeof(ID));
		strncpy(ID, result, sizeof(ID));
		strncpy(influxReturn->name[i], ID, 64);
		if(DEBUG) {
		  printf("Extracted value: %d %s\n", i, result);
		  printf("stored value: %d %s\n", i, influxReturn->name[i]);
		}
            }
        }
        i++;
    }
}


int countString(const char *haystack, const char *needle){
	int count = 0;
	const char *tmp = haystack;
	while( (tmp = strstr(tmp, needle)) )
	{
		count++;
		tmp++;
	}
	return count;
}

// extract substring starting from nth occurence of start_str to end_char
void extract_substring(char *source, char *start_str, char end_char, int occurrence, char *result) {
    char *start_pos = source;
    for (int i = 0; i < occurrence; i++) {
        start_pos = strstr(start_pos, start_str);
        if (start_pos == NULL) {
            // nth occurrence of starting string not found
            result[0] = '\0';
            return;
        }
        start_pos += strlen(start_str); // Move to the end of the current occurrence
    }

    char *end_pos = strchr(start_pos, end_char);
    if (end_pos == NULL) {
        // Ending character not found
        result[0] = '\0';
        return;
    }

    size_t length = end_pos - start_pos;
    strncpy(result, start_pos, length);
    result[length] = '\0'; // Null-terminate the result

}



// Callback function to handle the response from the InfluxDB server
size_t write_callback(char *contents, size_t size, size_t nmemb, void *userp) {

    size_t realsize = size * nmemb;

    // Parser column indicies
    int time_indx = 0;  // column containing time
    int scan_indx = 0;  // column containing scanID
    int text_indx = 0;  // column containing any text values (TARGET NAME)

    int pos = 0;	//token count
    char *token;
    char ID[64]="";
    char data[64]="";

    // vars to extract number of columns substring
    char *start_str = "\"columns\":\[";
    char *start_str2 = "\"values\":\[\[";
    char end_char = ']';
    char result[256];

    int nmeas = 0;	// num measurements
			// ACS3_DEV4_VQlo, ACS3_DEV4_VQhi, ACS3_DEV4_VIlo, ACS3_DEV4_VIhi  would be meas=4
			// udpPointing would be meas=1
			// HK_TEMP11 would be meas=1
			
    int ncols = 0;	// num influx columns
			// "time", "DEC", "RA", "scanID"
			// "time", "temp"
			// "time", "scanID", "volt"
			
    //get the number of measurements
    nmeas = countString(contents, "name");

    //get the number of columns
    extract_substring(contents, start_str, end_char, 1, result);
    ncols = countString(result, ",") + 1;

    // find the dynamically reported time and scanID columns
    token = strtok(result, ",");
    while (token != NULL ){

	memset(ID, '\0', sizeof(ID));
        strncpy(ID, token+1, MIN(strlen(token)-2, sizeof(ID)-1));

        if(!strcmp(ID, "time"))
            time_indx = pos;

        if(!strcmp(ID, "scanID"))
            scan_indx = pos;

        if(!strcmp(ID, "volt"))
            text_indx = pos;

        token = strtok(NULL, ",");
        pos++;
    }


    // dynamically allocate structure
    influxReturn = malloc(sizeof(*influxReturn));
    influxReturn->length = nmeas * (ncols-2);
    influxReturn->value = (float *)malloc(nmeas*(ncols-2) * sizeof(float));

    influxReturn->name   = (char **)malloc(nmeas*(ncols-2) * sizeof(char *));

    for(int i=0; i<nmeas*(ncols-2); i++)
	    influxReturn->name[i] = (char *)malloc(64*sizeof(char));

    numNames = nmeas*(ncols-2); // for use with free
    // load names into allocated space for names
    extract_all_keyword_values(contents, "name\":");

    for (int i=0; i<influxReturn->length; i++){
	    influxReturn->value[i] = 0.;
    }

    // Parse the InfluxDB return string 
    int data_indx=0;
    //    printf("%s\n", contents);
    for(int i=0; i<nmeas; i++){
        extract_substring(contents, start_str2, end_char, i+1, result);

        pos=0;
        token = strtok(result, ",");
        while (token != NULL ){

	        if(pos == time_indx){
	                memset(ID, '\0', sizeof(ID));
		        strncpy(ID, token+1, MIN(strlen(token)-2, sizeof(ID)-1));
			strcpy(influxReturn->time, ID);
                        if(DEBUG)
			  printf("time string = %s\n", ID);
	        }
    
	        else if(pos == scan_indx){
	                memset(ID, '\0', sizeof(ID));
		        strncpy(ID, token+1, MIN(strlen(token)-2, sizeof(ID)-1));
		        influxReturn->scanID = atoi(ID);
                        if(DEBUG)
			  printf("scanID string = %s\n", ID);
	        }

	        else if(pos == text_indx){
	                memset(data, '\0', sizeof(data));
		        strncpy(data, token, MIN(strlen(token), sizeof(data)-1));
		        influxReturn->value[data_indx] = atof(data);
                        if(DEBUG)
			  printf("data string %d = %s\n", data_indx, data);
		        data_indx++;
	        }
	        pos++;
                token = strtok(NULL, ",");
        }
    }
    return realsize;

}


void freeinfluxStruct(influxStruct *influxReturn) {
   free(influxReturn->value);
   for(int i=0; i<numNames; i++) {
     free(influxReturn->name[i]);
   }
   free(influxReturn->name);
   free(influxReturn);
}


// worker for the Influx DB query
// Takes: curl handle
// Operates: makes the function callback( )
influxStruct* influxWorker(CURL *curl, char *query)
{
    if (curl) {
        // Set the InfluxDB URL
        curl_easy_setopt(curl, CURLOPT_URL, INFLUXDB_URL);

        // Set the POST data (query)
        curl_easy_setopt(curl, CURLOPT_POSTFIELDS, query);

        // Set the callback function to handle the response
        curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_callback);

        // Make the HTTP POST request
        curl_easy_perform(curl);
	curl_easy_cleanup(curl);
    }
   // Cleanup Influx DB
   curl_global_cleanup();

   return influxReturn;      //return the struct
}

// Initialization function to connect to the Influx DB.
// Returns: The curl handle to pass to the handler function
CURL *init_influx()
{
    CURL *curl;

    // Initialize libcurl
    curl_global_init(CURL_GLOBAL_DEFAULT);

    // Create a CURL handle
    curl = curl_easy_init();
    return curl;
}

