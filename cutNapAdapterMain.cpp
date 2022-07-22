/********************************************************
 * cutNapAdapter: remove the adapters in NAP-seq data
 * jianhua yang yangjh7@mail.sysu.edu.cn
 * $ 2019/12/09
 *******************************************************/
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<ctype.h>
#include<math.h>
#include <getopt.h>
#include <map>
#include <algorithm>
#include <ios>
#include <iostream>
#include <string>
#include <ostream>
#include <fstream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <vector>
#include <zlib.h>

using namespace std;

#include "cutNapAdapter.h"

void usage(void);
char version[] = "cutNapAdapter version 0.1 [2015/11/06]";

int main(int argc, char **argv)
{
  int i = 0;
  int c = 0;
  char *read1File = NULL;
  char *read2File = NULL;
  gzFile read1fp;
  gzFile read2fp;
  char *outdir = NULL;
  char *prefix = NULL;
  const char *defpre = "cutNapAdapter";
  const char *defout  = "cutNapAdapterResults";
  char createDir[255];
  char outputDir[255];
  char *adapter1 = NULL;
  char *adapter2 = NULL;
  int showVersion = 0;
  int showHelp = 0;
  parameters paraInfo;

  if (argc == 1)
  {
    usage();
  }

  const char *shortOptions = "vhVda:o:m:l:e:1:2:b:f:c:C:";

  const struct option longOptions[] =
  {
    { "verbose" , no_argument , NULL, 'v' },
    { "help" , no_argument , NULL, 'h' },
    { "version" , no_argument , NULL, 'V' },
    { "discard" , no_argument , NULL, 'd' },
    { "output" , required_argument , NULL, 'o' },
    { "prefix" , required_argument , NULL, 'f' },
    { "adapter1" , required_argument , NULL, 'a' },
    { "adapter2" , required_argument , NULL, 'b' },
    { "min-match" , required_argument , NULL, 'm' },
    { "min-len" , required_argument , NULL, 'l' },
    { "p5-len" , required_argument , NULL, 'c' },
    { "p3-len" , required_argument , NULL, 'C' },
    { "error-rate" , required_argument , NULL, 'e' },
    { "read1" , required_argument , NULL, '1' },
    { "read2" , required_argument , NULL, '2' },
    {NULL, 0, NULL, 0} ,  /* Required at end of array. */
  };

  paraInfo.verbose = 0;
  paraInfo.minMatchLen = 5;
  paraInfo.minSeqLen = 5;
  paraInfo.maxMismatch = 3;
  paraInfo.keepAll = 1;
  paraInfo.maxErrorRate = 0.1;
  paraInfo.p5BarcodeLen = 6;
  paraInfo.p3BarcodeLen = 6;

  while ((c = getopt_long(argc, argv, shortOptions, longOptions, NULL)) >= 0)
  {
    switch (c)
    {
    case 'v':
      paraInfo.verbose = 1;
      break;
    case 'h':
      showHelp = 1;
      break;
    case 'V':
      showVersion = 1;
      break;
    case 'd':
      paraInfo.keepAll = 0;
      break;
    case 'o':
      outdir  = optarg;
      break;
    case 'f':
      prefix  = optarg;
      break;
    case 'a':
      adapter1 = optarg;
      break;
    case 'b':
      adapter2 = optarg;
      break;
    case '1':
      read1File = optarg;
      break;
    case '2':
      read2File = optarg;
      break;
    case 'm':
      paraInfo.minMatchLen = atoi(optarg);
      break;
    case 'l':
      paraInfo.minSeqLen = atoi(optarg);
      break;
    case 'c':
      paraInfo.p5BarcodeLen = atoi(optarg);
      break;
    case 'C':
      paraInfo.p3BarcodeLen = atoi(optarg);
      break;
    case 'e':
      paraInfo.maxErrorRate = atof(optarg);
      break;
    case '?':
      showHelp = 1;
      break;
    default:
      usage();
    }
  }

  if (adapter1 == NULL)
  {
    fprintf(stderr, "you must provide adapter sequence using -a option\n");
    usage();
  }
  if (adapter2 == NULL)
  {
    fprintf(stderr, "you must provide adapter sequence for read2 using -b option\n");
    usage();
  }

  if (paraInfo.minMatchLen <= 0)
  {
    fprintf(stderr, "-m option: minimum matches bewteen adapter and sequence must be large than 1 nt\n");
    usage();
  }
  if (paraInfo.minSeqLen <= 0)
  {
    fprintf(stderr, "-s option:  minimum sequence length after trim adapter must be large than 1 nt\n");
    usage();
  }
  if (read1File == NULL)
  {
    fprintf(stderr, "please read1 seqence file with option: -1 \n");
    usage();
  }

  if (strstr(read1File, ".gz") == NULL)
  {
    fprintf(stderr, "sequence file must be gzip format.\n");
    usage();
  }

  read1fp = gzopen(read1File, "rb");
  if (read1fp == NULL)
  {
    fprintf(stderr, "ERROR: Can't open gzip file: %s\n", read1File);
    usage();
  }

  if (read2File == NULL)
  {
    fprintf(stderr, "please read2 seqence file with option: -2 \n");
    usage();
  }
  if (strstr(read2File, ".gz") == NULL)
  {
    fprintf(stderr, "sequence file must be gzip format.\n");
    usage();
  }
  read2fp = gzopen(read2File, "r");
  if (read2fp == NULL)
  {
    fprintf(stderr, "ERROR: Can't open gzip file: %s\n", read2File);
    usage();
  }

  if (showVersion)
  {
    fprintf(stderr, "%s", version);
    exit(1);
  }

  if (showHelp)
  {
    usage();
    exit(1);
  }

  if (paraInfo.verbose) fprintf(stderr, "#clip sequencing reads...\n");

  strcpy(createDir, "mkdir -p ");
  if (outdir != NULL) {
    strcpy(outputDir, outdir);
    strcat(createDir, outdir);
  }
  else {
    strcpy(outputDir, defout);
    strcat(createDir, defout);
  }
  strcat(outputDir, "/");
  if (prefix != NULL) {
    strcat(outputDir, prefix);
  }
  else {
    strcat(outputDir, defpre);
  }
  if (paraInfo.verbose)  fprintf(stderr, "#create dir \"%s\"\n", createDir);
  system(createDir);

  cutNapAdapter(&paraInfo, read1fp, read2fp, adapter1, adapter2, outputDir);

  if (paraInfo.verbose)  fprintf(stderr, "close read files\n");
  gzclose(read1fp);
  if (read2fp != NULL) gzclose(read2fp);
  return 0;
}

void usage(void)
{
  fprintf(stderr, "%s",
          "Usage:  cutNapAdapter [options] -a <adapter5> -b <adapter3> -1 <read1 fiile, gzip format> -2 <read2 File, gzip format>\n\
[options]\n\
-o/--output <string>       : output dir[default=cutNapAdapterResults]\n\
-f/--prefix <string>       : prefix for output files[default=cutNapAdapter]\n\
-a/--adapter1 <string>     : adapter1 sequence [required]\n\
-b/--adapter2 <string>     : adapter2 sequence [required for paired-end sequencing]\n\
-c                         : barcode length of 5'-end[default=6]\n\
-C                         : barcode length of 3'-end[default=6]\n\
-m/--min-match             : minimum matches bewteen adapter and sequences[default=5]\n\
-l/--min-len               : minimum sequence length after trim adapter[default=5]\n\
-e/--error-rate            : maximum error rate [default<=0.1]\n\
-v/--verbose               : verbose information\n\
-V/--version               : cutNapAdapter version\n\
-h/--help                  : help informations\n\
");
  exit(1);
}
