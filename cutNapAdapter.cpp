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
#include "bioUtils.h"
#include "cutNapAdapter.h"

void cutNapAdapter(parameters *paraInfo, gzFile fp1, gzFile fp2, char *adapter1, char *adapter2, char *outdir)
{
	char read1OutFile[512];
	char read2OutFile[512];
	char read1UnpairOutFile[512];
	char read2UnpairOutFile[512];
	strcpy(read1OutFile, outdir);
	strcpy(read1UnpairOutFile, outdir);
	strcat(read1OutFile, "_R1.clip5p3pAdapter.pair.fq.gz");
	strcat(read1UnpairOutFile, "_R1.clip5p3pAdapter.unpair.fq.gz");
	strcpy(read2OutFile, outdir);
	strcat(read2OutFile, "_R2.clip5p3pAdapter.pair.fq.gz");
	strcpy(read2UnpairOutFile, outdir);
	strcat(read2UnpairOutFile, "_R2.clip5p3pAdapter.unpair.fq.gz");

	if (paraInfo->verbose) fprintf(stderr, "cut adapter in read2\n");
	fprintf(stderr, "###statistical information###\n");
	cutPEAdapter(paraInfo, fp1, fp2, adapter1, adapter2, read1OutFile, read1UnpairOutFile, read2OutFile, read2UnpairOutFile);
}

void cutPEAdapter(parameters *paraInfo, gzFile fp1, gzFile fp2, char *adapter5, char *adapter3,
                  char *read1OutFile, char *read1UnpairOutFile, char *read2OutFile, char *read2UnpairOutFile)
{
	char *line1 = NULL;
	char *line2 = NULL;
	int lessLenNum = 0;
	int lessOneNum = 0;
	int largeLenNum = 0;
	int noAdapterNum = 0;

	gzFile read1ofp = (gzFile) gzopen(read1OutFile, "w6");
	if (read1ofp == NULL)
	{
		fprintf(stderr, "ERROR: Can't open %s\n", read1OutFile);
		exit(1);
	}
	gzFile read2ofp = (gzFile) gzopen(read2OutFile, "w6");
	if (read2ofp == NULL)
	{
		fprintf(stderr, "ERROR: Can't open %s\n", read2OutFile);
		exit(1);
	}
	gzFile read1ufp = (gzFile) gzopen(read1UnpairOutFile, "w6");
	if (read1ufp == NULL)
	{
		fprintf(stderr, "ERROR: Can't open %s\n", read1UnpairOutFile);
		exit(1);
	}
	gzFile read2ufp = (gzFile) gzopen(read2UnpairOutFile, "w6");
	if (read2ufp == NULL)
	{
		fprintf(stderr, "ERROR: Can't open %s\n", read2UnpairOutFile);
		exit(1);
	}

	if (paraInfo->verbose) fprintf(stderr, "clip reads\n");
	while (line1 = getGzipLine(fp1))
	{
		line2 = getGzipLine(fp2);
		if (gzeof(fp1) || gzeof(fp2) || line1 == NULL || line2 == NULL)
		{
			safeFree(line1);
			safeFree(line2);
			break;
		}
		if ((line1[0] != '@' || line2[0] != '@') && (line1[0] != '>' || line2[0] != '>'))
		{
			fprintf(stderr, "error read format: %c\n", line1[0]);
			fprintf(stderr, "error read format: %c\n", line2[0]);
			fflush(stderr);
			safeFree(line1);
			safeFree(line2);
			break;
		}
		if ((line1[0] == '@' && line2[0] == '@') || (line1[0] == '>' && line2[0] == '>'))
		{
			readInfo *readPtr1 = getGzipOneRead(fp1, line1, line1[0]);
			readInfo *readPtr2 = getGzipOneRead(fp2, line2, line2[0]);
			// set one ptr to one struct
			readInfo read1 = *readPtr1;
			readInfo read2 = *readPtr2;

			int minusFlag = 0;
			int p5aFlag = 0;
			int p3aFlag = 0;

			int bcLen = paraInfo->p3BarcodeLen;
			if (bcLen < 1) bcLen = 1;
			char *bcSeq = (char *) safeMalloc(bcLen + 1);
			bcSeq[0] = 'N';
			bcSeq[1] = '\0';

			int read1p5aPos = cut5pAdapter(paraInfo, adapter5, read1.readSeq, paraInfo->p5BarcodeLen);
			if (read1p5aPos >= paraInfo->minSeqLen)
			{
				int pos = read1p5aPos + paraInfo->p5BarcodeLen + 1;
				read1.readSeq = read1.readSeq + pos;
				if (read1.quality != NULL)
					read1.quality = read1.quality + pos;
				p5aFlag = 5;
			}
			int read2p5aPos = cut5pAdapter(paraInfo, adapter5, read2.readSeq, paraInfo->p5BarcodeLen);
			if (read2p5aPos >= paraInfo->minSeqLen)
			{
				int pos = read2p5aPos + paraInfo->p5BarcodeLen + 1;
				read2.readSeq = read2.readSeq + pos;
				if (read2.quality != NULL)
					read2.quality = read2.quality + pos;
				minusFlag = 1;
				p5aFlag = 5;
			}
			int read1p3aPos = cut3pAdapter(paraInfo, adapter3, read1.readSeq, paraInfo->p3BarcodeLen);
			if (read1p3aPos >= paraInfo->minSeqLen)
			{
				int pos = read1p3aPos - paraInfo->p3BarcodeLen;

				strncpy(bcSeq, read1.readSeq + pos, paraInfo->p3BarcodeLen);
				bcSeq[paraInfo->p3BarcodeLen] = '\0';

				read1.readSeq[pos] = '\0';
				if (read1.quality != NULL)
					read1.quality[pos] = '\0';
				p3aFlag = 3;
			}
			int read2p3aPos = cut3pAdapter(paraInfo, adapter3, read2.readSeq, paraInfo->p3BarcodeLen);
			if (read2p3aPos >= paraInfo->minSeqLen)
			{
				int pos = read2p3aPos - paraInfo->p3BarcodeLen;

				strncpy(bcSeq, read2.readSeq + pos, paraInfo->p3BarcodeLen);
				bcSeq[paraInfo->p3BarcodeLen] = '\0';

				read2.readSeq[pos] = '\0';
				if (read2.quality != NULL)
					read2.quality[pos] = '\0';
				minusFlag = 1;
				p3aFlag = 3;
			}
			char *rcAdapter5 = strClone(adapter5);
			reverseComp(rcAdapter5);
			char *rcAdapter3 = strClone(adapter3);
			reverseComp(rcAdapter3);
			int read1p5aRcPos = cut3pAdapter(paraInfo, rcAdapter5, read1.readSeq, paraInfo->p5BarcodeLen);
			if (read1p5aRcPos >= paraInfo->minSeqLen)
			{
				int pos = read1p5aRcPos - paraInfo->p5BarcodeLen;
				read1.readSeq[pos] = '\0';
				if (read1.quality != NULL)
					read1.quality[pos] = '\0';
				minusFlag = 1;
				p5aFlag = 5;
			}
			int read2p5aRcPos = cut3pAdapter(paraInfo, rcAdapter5, read2.readSeq, paraInfo->p5BarcodeLen);
			if (read2p5aRcPos >= paraInfo->minSeqLen)
			{
				int pos = read2p5aRcPos - paraInfo->p5BarcodeLen;
				read2.readSeq[pos] = '\0';
				if (read2.quality != NULL)
					read2.quality[pos] = '\0';
				p5aFlag = 5;
			}
			int read1p3aRcPos = cut5pAdapter(paraInfo, rcAdapter3, read1.readSeq, paraInfo->p3BarcodeLen);
			if (read1p3aRcPos >= paraInfo->minSeqLen)
			{
				int pos = read1p3aRcPos + paraInfo->p3BarcodeLen + 1;

				strncpy(bcSeq, read1.readSeq + pos - paraInfo->p3BarcodeLen, paraInfo->p3BarcodeLen);
				bcSeq[paraInfo->p3BarcodeLen] = '\0';
				reverseComp(bcSeq);

				read1.readSeq = read1.readSeq + pos;
				if (read1.quality != NULL)
					read1.quality = read1.quality + pos;
				minusFlag = 1;
				p3aFlag = 3;
			}
			int read2p3aRcPos = cut5pAdapter(paraInfo, rcAdapter3, read2.readSeq, paraInfo->p3BarcodeLen);
			if (read2p3aRcPos >= paraInfo->minSeqLen)
			{
				int pos = read2p3aRcPos + paraInfo->p3BarcodeLen + 1;

				strncpy(bcSeq, read2.readSeq + pos - paraInfo->p3BarcodeLen, paraInfo->p3BarcodeLen);
				bcSeq[paraInfo->p3BarcodeLen] = '\0';
				reverseComp(bcSeq);

				read2.readSeq = read2.readSeq + pos;
				if (read2.quality != NULL)
					read2.quality = read2.quality + pos;
				p3aFlag = 3;
			}

			safeFree(rcAdapter5);
			safeFree(rcAdapter3);

			int read1Len = strlen(read1.readSeq);
			int read2Len = strlen(read2.readSeq);
			if ( read1Len >= paraInfo->minSeqLen && read2Len >= paraInfo->minSeqLen)
			{
				largeLenNum++;
				if (minusFlag == 1)
				{
					reverseCompPEread(&read1, &read2);
				}
				outputPEresults(&read1, &read2, read1ofp, read2ofp, p5aFlag, p3aFlag, bcSeq);
			}
			else
			{
				outputPEresults(&read1, &read2, read1ufp, read2ufp, p5aFlag, p3aFlag, bcSeq);
			}

			safeFree(bcSeq);
			freeRead(readPtr1);
			freeRead(readPtr2);
		}
	} // gofile while loops
	gzclose(read1ofp);
	gzclose(read2ofp);
	gzclose(read1ufp);
	gzclose(read2ufp);
	fprintf(stderr, "#The number of reads with length >= %d is %d(%.2f)\n", paraInfo->minSeqLen, largeLenNum, largeLenNum / (double)(largeLenNum + noAdapterNum + lessLenNum) * 100);
	fprintf(stderr, "#The number of reads with length <  %d is %d(%.2f)\n", paraInfo->minSeqLen, lessLenNum, lessLenNum / (double)(largeLenNum + noAdapterNum + lessLenNum) * 100);
	fprintf(stderr, "#The number of reads with length <=  1 is %d(%.2f)\n", lessOneNum, lessOneNum / (double)(largeLenNum + noAdapterNum + lessLenNum) * 100);
	fprintf(stderr, "#The number of reads without adapter   is %d(%.2f)\n", noAdapterNum, noAdapterNum / (double)(largeLenNum + noAdapterNum + lessLenNum) * 100);
}

void outputPEresults(readInfo *read1, readInfo *read2, gzFile outfp1, gzFile outfp2, int p5Flag, int p3Flag, char *bcSeq)
{
	int fieldNum = 0;
	char **fields = NULL;
	fields = splitWhitespace(read1->readName, &fieldNum);
	gzprintf(outfp1, "%s:%s:%d\n", fields[0], bcSeq, (p5Flag + p3Flag));
	gzprintf(outfp1, "%s\n", read1->readSeq);
	if (read1->qualityName != NULL && read1->quality != NULL)
	{
		gzprintf(outfp1, "%s\n", read1->qualityName);
		gzprintf(outfp1, "%s\n", read1->quality);
	}
	freeWords(fields, fieldNum);
	fields = splitWhitespace(read2->readName, &fieldNum);
	gzprintf(outfp2, "%s:%s:%d\n", fields[0], bcSeq, (p5Flag + p3Flag));
	gzprintf(outfp2, "%s\n", read2->readSeq);
	if (read2->qualityName != NULL && read2->quality != NULL)
	{
		gzprintf(outfp2, "%s\n", read2->qualityName);
		gzprintf(outfp2, "%s\n", read2->quality);
	}
	freeWords(fields, fieldNum);
}

char *getGzipLine(gzFile gzf)
/* get whole line from gzip
*/
{
	const int lineLen = 512;
	char s[lineLen];
	char *line = NULL;
	char *cp = NULL;
	int done = FALSE;
	line = NULL;

	do {
		if (gzgets(gzf, s, lineLen) == Z_NULL)
			break; /* EOF */
		/* for unix OS */
		cp = strchr(s, '\n');
		if (cp != NULL) {
			*cp   = '\0';
			done  = TRUE;
		}
		/* for window OS */
		cp = strchr(s, '\r');
		if (cp != NULL) {
			*cp  = '\0';
			done = TRUE;
		}
		if (line == NULL)
			line = (char *)safeMalloc(strlen(s) + 1); /* don't use malloc, for we will using strcat function, so we must initilized the line with '0' */
		else
			line = (char *)safeRealloc(line, strlen(s) + strlen(line) + 1);
		strcat(line, s);
	} while (!done);

	return line;
}

readInfo *getGzipOneRead(gzFile gzf, char *line, char format)
{
	char *seq = NULL;
	char *qualityName = NULL;
	char *quality = NULL;
	seq = getGzipLine(gzf);
	if (format == '@')
	{
		qualityName = getGzipLine(gzf);
		quality = getGzipLine(gzf);
	}
	readInfo *read = (readInfo *)safeMalloc(sizeof(readInfo));
	read->readName    = strClone(line);
	read->readSeq     = seq;
	read->qualityName = qualityName;
	read->quality     = quality;
	return read;
}

void reverseCompPEread(readInfo *read1, readInfo *read2)
{
	reverseComp(read1->readSeq);
	reverseComp(read2->readSeq);
	reverseBytes(read1->quality, strlen(read1->quality));
	reverseBytes(read2->quality, strlen(read2->quality));
}

int cut5pAdapter(parameters *paraInfo, char *adapter, char *seq, int barcodeLen)
{
	int i = 0;
	int j = 0;
	int headMatchLen = 6;
	int minLen = paraInfo->minMatchLen;
	int maxMismatchNum = paraInfo->maxMismatch;
	int start = 0;
	int adapterLen = strlen(adapter);
	int seqLen = strlen(seq);
	int idxPos = adapterLen + barcodeLen;
	int pos = 0;
	if (idxPos > seqLen - barcodeLen)  idxPos = seqLen - barcodeLen;
	double errorRate = 100;
	for (start = idxPos; start >= minLen; start--)
	{
		i = adapterLen - 1;
		j = start;
		int match = 0;
		int mismatch = 0;
		int idx = 0; // skip the reads with mismatch at first six bases
		for (; i >= 0 && j >= 0; i--, j--)
		{
			if (adapter[i] == seq[j])
			{
				match += 1;
			}
			else
			{
				mismatch += 1;
			}
			/*if (idx < headMatchLen && mismatch > 0)
			{
				// skip the reads with mismatch at first six bases
				mismatch = 100;
				break;
			}*/
			if (mismatch > paraInfo->maxMismatch) break;
			//idx++;
		}
		int totalLen = match + mismatch;
		errorRate = (double)mismatch / (double)totalLen;
		if (errorRate <= paraInfo->maxErrorRate && j <= 0)
		{
			pos = start;
			break;
		}
	}
	return pos;
}

int cut3pAdapter(parameters *paraInfo, char *adapter, char *seq, int barcodeLen)
{
	int i = 0;
	int j = 0;
	int headMatchLen = 6;
	int minLen = paraInfo->minMatchLen;
	int maxMismatchNum = paraInfo->maxMismatch;
	int start = 0;
	int adapterLen = strlen(adapter);
	int seqLen = strlen(seq);
	int idxPos = seqLen - adapterLen - barcodeLen;
	int pos = 0;
	if (idxPos < minLen + barcodeLen)  idxPos = minLen + barcodeLen;
	double errorRate = 100;
	for (start = idxPos; start < seqLen - minLen; start++)
	{
		i = 0;
		j = start;
		int match = 0;
		int mismatch = 0;
		int idx = 0; // skip the reads with mismatch at first six bases
		for (; i < adapterLen && j < seqLen; i++, j++)
		{
			if (adapter[i] == seq[j])
			{
				match += 1;
			}
			else
			{
				mismatch += 1;
			}
			/*if (idx < headMatchLen && mismatch > 0)
			{
				// skip the reads with mismatch at first six bases
				mismatch = 100;
				break;
			}*/
			if (mismatch > paraInfo->maxMismatch) break;
			//idx++;
		}
		int totalLen = match + mismatch;
		errorRate = (double)mismatch / (double)totalLen;
		if (errorRate <= paraInfo->maxErrorRate && j >= seqLen - 1)
		{
			pos = start;
			break;
		}
	}
	return pos;
}

void freeRead(readInfo *read)
{
	safeFree(read->readName);
	safeFree(read->readSeq);
	if (read->qualityName != NULL && read->quality != NULL)
	{
		safeFree(read->qualityName);
		safeFree(read->quality);
	}
	safeFree(read);
}

void copyReadInfo(readInfo *tar, readInfo *org)
{
	tar->readName    = org->readName;
	tar->readSeq     = org->readSeq;
	tar->qualityName = org->qualityName;
	tar->quality     = org->quality;

}
readInfo *getOneRead(FILE *fp, char *line, char format)
{
	char *seq = NULL;
	char *qualityName = NULL;
	char *quality = NULL;
	seq = getLine(fp);
	if (format == '@')
	{
		qualityName = getLine(fp);
		quality = getLine(fp);
	}
	readInfo *read = (readInfo *)safeMalloc(sizeof(readInfo));
	read->readName    = strClone(line);
	read->readSeq     = seq;
	read->qualityName = qualityName;
	read->quality     = quality;
	return read;
}

void freeFloatMatrix(double **matrix)
/* free float matrix */
{
	safeFree(matrix[0]);
	safeFree(matrix);
}

void freeIntMatrix(int **matrix)
/*free int matrix */
{
	safeFree(matrix[0]);
	safeFree(matrix);
}

double max(double m[], int len)
/* return maxinum value */
{
	int i = 0;
	double maxScore = m[0];

	for (i = 1; i < len; i++)
	{
		if (m[i] > maxScore)
		{
			maxScore = m[i];
		}
	}
	return maxScore;
}

int argmax(double m[], int len)
/* return the index for maximum value */
{
	int i = 0;
	int best_i = 0;
	double maxScore = m[0];
	for (i = 1; i < len; i++)
	{
		if (m[i] > maxScore)
		{
			maxScore = m[i];
			best_i = i;
		}
	}
	return best_i;
}

double score(char a, char b)
/* score match and mismatch */
{
	if (a == b) return MATCH;
	else if (a != b) return MISMATCH;
	return 0;
}
