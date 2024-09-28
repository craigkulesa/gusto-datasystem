#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <curl/curl.h>
#include <stdint.h>


// global variable to be returned by worker from callback fn

// struct for the influxDB worker to return
typedef struct 
{
  char     time[64];	// time string
  char     text[64];	// text string for TARGETNAMES
			
  int16_t  scanID;	// scanID
  float    *value;	// values of measurements
  char     **name;	// names of measurements
			
  char     **cols;	// how many columns of values are returned from influx
  size_t   length;	// how many time are a measurement returned from influx
} influxStruct;
extern influxStruct *influxReturn;


// InfluxDB parameters
#define INFLUXDB_URL "http://localhost:8086/query?&db=gustoDBlp"

// Callback function to handle the response from the InfluxDB server
size_t write_callback(char *contents, size_t size, size_t nmemb, void *userp);

// Initialization function to connect to the Influx DB.
// Returns: The curl handle to pass to the handler function
CURL *init_influx();

influxStruct* influxWorker(CURL *curl, char *query);

void freeinfluxStruct(influxStruct *influxReturn);
