#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define strequals(src,dst) (!strcmp(src,dst))
#define null NULL
#define DEBUG false

const int READ_TO_EOF = 0;
const int READ_TO_TAB = 1;
const int READ_TO_NL  = 10;

const int MODE_COUNT_ONLY = 0;
const int MODE_PGM        = 1;
const int MODE_BIN        = 2;

const unsigned char MAJOR_A =(unsigned char)128;
const unsigned char MAJOR_C =64;
const unsigned char MAJOR_G =32;
const unsigned char MAJOR_T =16;

const unsigned char MAJOR_MASK= (unsigned char)240;

const unsigned char MINOR_A =8;
const unsigned char MINOR_C =4;
const unsigned char MINOR_G =2;
const unsigned char MINOR_T =1;

const unsigned char MINOR_MASK= 15;

unsigned char * DB_ALLELLS=null;
unsigned int  * DB_ADDRESS=null;

const int DATA_OUT_MODE_ENCODED = 0;
const int DATA_OUT_MODE_DIFF    = 1;
const int DATA_OUT_MODE_ASCII   = 4;

const int PATIENT_ALL = -1;
int targetpatient=0;

int dataoutmode= DATA_OUT_MODE_ENCODED;

char readbuffer[1024];

int outmode=MODE_COUNT_ONLY;

FILE * fout=stdout;

unsigned char allell(char cmajor,char cminor) {
    unsigned char c=0;
         if (cmajor=='A') c=MAJOR_A;
    else if (cmajor=='C') c=MAJOR_C;
    else if (cmajor=='G') c=MAJOR_G;
    else if (cmajor=='T') c=MAJOR_T;
         if (cminor=='A') c+=MINOR_A;
    else if (cminor=='C') c+=MINOR_C;
    else if (cminor=='G') c+=MINOR_G;
    else if (cminor=='T') c+=MINOR_T;
    return c;
}

unsigned char headerallell(char *s) {
    int i=0;
    while (s[i]!=':') i++;
    i++;
    while (s[i]!=':') i++;
    i++;
    int j=i;
    while (s[j]!=':') j++;
    j++;
    // destructive debug hack
    s[0]=s[i];
    s[1]=s[j];
    return allell(s[i],s[j]);
}


int tabber(char *s) {
  char c=getchar();
  int i=0;
  s[i++]=c;

  while (c!='\n' && c!='\t' && c!=EOF)
    s[i++]=c=getchar();
  s[i]=0;

  if (c=='\t') return READ_TO_TAB;
  if (c=='\n') return READ_TO_NL;

  return READ_TO_EOF;
}

int  count=0;
long missingdata=0;

void readDBAllells() {
    int index=count-5;
    DB_ALLELLS[index]=headerallell(readbuffer);
    if (DEBUG)
        fprintf(stderr,"%i: %c %c : %03u\n",
                       index,
                       readbuffer[0],
                       readbuffer[1],
                       (unsigned int)DB_ALLELLS[index]);
}

void processPatientAllells() {
    if (dataoutmode==DATA_OUT_MODE_ASCII) {
        fputc(readbuffer[0],fout);
        fputc(readbuffer[1],fout);
        return;
    }
    if (readbuffer[0]!=0 && readbuffer[1]!=0) {
      char c=allell(readbuffer[0],readbuffer[0]);
      if (dataoutmode==DATA_OUT_MODE_ENCODED)
        fputc(c,fout);
      else {
          bool majorMatch=(readbuffer[0]&MAJOR_MASK)==(DB_ALLELLS[count]&MAJOR_MASK);
          bool minorMatch=(readbuffer[0]&MINOR_MASK)==(DB_ALLELLS[count]&MINOR_MASK);
          if (majorMatch && minorMatch)
            fputc(255,fout);
          else if (majorMatch && !minorMatch)
            fputc(127,fout);
          else if (               minorMatch)
            fputc( 63,fout);
          else
            fputc( 31,fout);
      }
    } else {
      fputc(0,fout);
      missingdata++;
    }
}


int main (int nargs, char **args) {
  int  marker=0;
  int  nlines=0;
  int  w=-1;
  int  nchar=0;

  for (int i=1; i<nargs; i++) {
    if (strequals(args[i],"-h")) {
       printf("usage : %s -counts | -pgm | -bin\n",args[0]);
       return 0;
    } else if (strequals(args[i],"-counts")) {
       outmode=MODE_COUNT_ONLY;
    } else if (strequals(args[i],"-pgm")) {
       outmode=MODE_PGM;
    } else if (strequals(args[i],"-bin")) {
       outmode=MODE_BIN;
    } else if (strequals(args[i],"-nchar")) {
       sscanf(args[++i],"%i",&nchar);
    } else if (strequals(args[i],"-allpatients")) {
        targetpatient=PATIENT_ALL;
    } else if (strequals(args[i],"-patient")) {
       sscanf(args[++i],"%i",&targetpatient);
    } else if (strequals(args[i],"-ascii")) {
       dataoutmode=DATA_OUT_MODE_ASCII;
    } else if (strequals(args[i],"-encoded")) {
       dataoutmode=DATA_OUT_MODE_ENCODED;
    } else if (strequals(args[i],"-diff")) {
       dataoutmode=DATA_OUT_MODE_DIFF;
    }
  }

  if (outmode!=MODE_COUNT_ONLY)
    DB_ALLELLS= (unsigned char *)malloc(524288);

  while ((marker=tabber(readbuffer))!=READ_TO_EOF) {
    if (marker==READ_TO_TAB) {
      if (outmode!=MODE_COUNT_ONLY && w>1) {
        if (targetpatient<0 || targetpatient==nlines-1)
          if (count>4)
            processPatientAllells();
      } else if (w<1){
        if (count>4) readDBAllells();
      }
      count++;
    } else { // NEWLINE
      if (outmode==MODE_COUNT_ONLY) printf("%i\n",count);
      if (w<0) {
         w=count;
         readDBAllells();
         // should create target file
         // sprintf(outfile,"patient.%04i.pgm",nlines);
         // outfile=fopen(outfilename,"wb");
         if (outmode!=MODE_COUNT_ONLY) {
           if (fout!=null && fout!=stdout)
             fclose(fout);
           char filename[32];
           snprintf(filename,32,"p%04i.c%03i.%s",nlines,nchar,outmode==MODE_PGM?"pgm":"bin");
           //fout=fopen(finename,"wb");
         }
         if (outmode==MODE_PGM) {
             int sqrt_w=(int)sqrt(w);
             int sqrt_w_reminder = (int)(w/sqrt_w);
             // NOTE: EXTRA REMINDER: w%sqrt_w_reminder -> may need to add row
             fprintf(fout,"P5\n#CRIC %i\n%i %i\n255\n",w-5,sqrt_w,sqrt_w_reminder);
         }
      } else if (outmode!=MODE_COUNT_ONLY) {
        if (targetpatient<0 || targetpatient==nlines-1)
           processPatientAllells();
        // should start new output file?
        if (targetpatient==nlines-1) return 0;
      } else {
      }
      count=0;
      nlines++;
    }
  }
  if (missingdata>0) fprintf(stderr,"Warning: missing data at %li places.\n",missingdata);
  if (outmode==MODE_COUNT_ONLY) printf("\n%i\n",nlines);
  return 0;
}
