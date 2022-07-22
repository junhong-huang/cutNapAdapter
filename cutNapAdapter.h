/********************************************************
 * cutNapAdapter: remove the adapters
 * jianhua yang yangjh7@mail.sysu.edu.cn
 * $ 2015/11/6
 *******************************************************/

#ifndef cutNapAdapter_HEAD_H
#define cutNapAdapter_HEAD_H

#ifndef MIN
#define MIN(x, y) (x)<(y)?(x):(y)
#endif

#ifndef MAX
#define MAX(x, y) (x)>(y)?(x):(y)
#endif

#define GAP_OPEN -2.0
#define GAP_CONT -2.0
#define MATCH 1.0
#define MISMATCH -2.0

#define H 3

struct parameterInfo
{
	int verbose;
	int minMatchLen;
	int minSeqLen;
	int maxMismatch;
	int keepAll;
	//int barcodeLen;
	int p5BarcodeLen;
	int p3BarcodeLen;
	double maxErrorRate;
};

typedef struct parameterInfo parameters;

struct readInfo
{
	char *readName;
	char *readSeq;
	char *qualityName;
	char *quality;
};

typedef struct readInfo readInfo;

void cutNapAdapter(parameters *paraInfo, gzFile fp1, gzFile fp2, char *adapter1, char *adapter2, char *outdir);


void outputPEresults(readInfo *read1, readInfo *read2, gzFile outfp1, gzFile outfp2, int p5Flag, int p3Flag, char *bcSeq);

void cutPEAdapter(parameters *paraInfo, gzFile fp1, gzFile fp2, char *adapter5, char *adapter3,
                  char *read1OutFile, char *read1UnpairOutFile, char *read2OutFile, char *read2UnpairOutFile);

char *getGzipLine(gzFile gzf);

readInfo *getGzipOneRead(gzFile gzf, char *line, char format);

void freeFloatMatrix(double **matrix);

void freeIntMatrix(int **matrix);

int argmax(double m[], int len);

double max(double m[], int len);

double score(char a, char b);

int scanAdapter(char *adapter, char *seq, parameters *paraInfo);

readInfo *getOneRead(FILE *fp, char *line, char format);

void freeRead(readInfo *read);

void copyReadInfo(readInfo *tar, readInfo *org);

void reverseCompPEread(readInfo *read1, readInfo *read2);

int cut5pAdapter(parameters *paraInfo, char *adapter, char *seq, int barcodeLen);

int cut3pAdapter(parameters *paraInfo, char *adapter, char *seq, int barcodeLen);

#endif /* End cutNapAdapter_HEAD_H */
