#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int nocolumns=2;
int minbins=1;
int maxbins=1000;
int validateinput(FILE *input);
void checkwindowsize(int w, int inputsize);
double findS_t(int inputsize, double mtx[][nocolumns]);
double findS_ybar(int w, int nosegments, double boxcardata[], int remainder);
double findS_tybar(int inputsize, double mtx[][nocolumns], double boxcardata[], int w);
double findS_tt(int inputsize, double mtx[][nocolumns]);
double findmax(double detrendmtx[][nocolumns],int inputsize);
double findmin(double detrendmtx[][nocolumns],int inputsize);
double findmean(double detrendmtx[][nocolumns],int inputsize);
double findsampleSD(double detrendmtx[][nocolumns],int inputsize, double mean);

int main(int argc, char* argv[])
{
/* open file for reading */
FILE *input;
char fname[]="proj3input.dat";
input = fopen(fname, "r");
/*check file has opened */
if(input == (FILE*) NULL)
    {
    printf("Can not open file %s\n", fname);
    exit(EXIT_FAILURE);
    }
/*finds size of/checks input*/
int inputsize = validateinput(input);
printf("Input size is %d\n", inputsize);
/*checks window size is appropriate for size of input
(say between 1/20th and 1/200th size of input).
Checks numbers of arguments, puts window size into variable w */
int w;
int nobins;
int validarguments;
validarguments = (argc == 3);
validarguments = validarguments && sscanf(argv[1], "%d", &w);
validarguments = validarguments && sscanf(argv[2], "%d", &nobins);
validarguments = ((1<=w)&&(w<=inputsize));
validarguments = ((minbins<=nobins)&&(nobins<=maxbins));
if (validarguments ==0)
{
    printf("*******************************\n");
    printf("* Input error.                *\n");
    printf("*Execution: 1205556_proj3 a b *\n");
    printf("*a is int between %d and %d   *\n",(inputsize/200),(inputsize/20));
    printf("*b is int between %d and %d   *\n",minbins, maxbins);
    printf("*******************************\n");
    exit(EXIT_FAILURE);
}
checkwindowsize(w,inputsize);
/*reads file into 2D array*/
double mtx[inputsize][nocolumns];
int i, j;
for (i=0;i<inputsize;i++)
    {
    for (j=0;j<nocolumns ;j++)
        {
    fscanf(input, "%lf ", &mtx[i][j]);
        }
    }
/******BOXCAR AVERAGE******/
int nosegments = ceil(((double)inputsize)/((double)w));
double boxcardata[nosegments];
double summation;
for (i =0;i<(nosegments-1);i++)
{
            summation = 0;
    for (j=0;j<w;j++)
    {

        summation = summation+mtx[j+(w*i)][1];
    }
    boxcardata[i]=(summation/((double)w));
}
/*average for remainder of data*/
int remainder = inputsize%w;
if (remainder ==0)
{
    remainder=w;
}
i=nosegments-1;
summation = 0;
for (j=0;j<(remainder);j++)
{
    summation = summation+mtx[j+(w*i)][1];
}
boxcardata[i]=(double)summation/((double)(remainder));
/****find requisite values for equations****/
double S_t= findS_t(inputsize, mtx);
double S_ybar=findS_ybar(w, nosegments, boxcardata, remainder);
double S_tybar=findS_tybar(inputsize, mtx, boxcardata, w);
double S_tt=findS_tt(inputsize, mtx);
/***find m ***/
double N=inputsize;
double m;
m= (((N*S_tybar) - (S_t*S_ybar))/(N*S_tt - (S_t*S_t)));
/***find c ***/
double c;
c= (((S_ybar*S_tt)-(S_t*S_tybar))/((N*S_tt)-(S_t*S_t)));
/********Detrend the data******/
double detrendmtx[inputsize][nocolumns];
for (i=0;i<inputsize;i++)
    {
    detrendmtx[i][0]=mtx[i][0];
    detrendmtx[i][1]=mtx[i][1]-(mtx[i][0]*m)-c;
    }
/*******find minimum, maximum, mean and standard deviation***/
double maximum=findmax(detrendmtx,inputsize);
double minimum=findmin(detrendmtx,inputsize);
double mean=findmean(detrendmtx,inputsize);
double sampleSD=findsampleSD(detrendmtx,inputsize,mean);
/***Choose histogram bins and count number of samples in each bin***/
double histogram[nobins][3];
/*initialise histogram with zero values*/
for (i=0;i<nobins;i++)
{
    histogram[i][0]=0;
}
/*______________*/
double binsize=(maximum-minimum)/nobins;
/*write bin boundaries*/
for (i=0;i<nobins;i++)
{
    histogram[i][1]=minimum + i*binsize;
    histogram[i][2]=minimum + (i+1)*binsize;

}
double temp;
for(i=0;i<inputsize;i++)
{
    temp=((detrendmtx[i][1]-minimum)/binsize);
    for (j=0;j<nobins;j++)
    {
        if ((temp>=j)&&(temp<j+1)) {histogram[j][0]++;}
    }
}
/*open first file for writing */
FILE *studnum_proj3_1;
char fname2[]="1205556_proj3_1.out";
studnum_proj3_1 = fopen(fname2, "w");
/*check file has opened*/
if(studnum_proj3_1 == (FILE*) NULL)
    {
    printf("Can not open file %s\n", fname2);
    exit(EXIT_FAILURE);
    }
/*prints statistics into studnum_proj3_1*/
fprintf(studnum_proj3_1, "Maximum: %lf ",maximum);
fprintf(studnum_proj3_1, "Minimum: %lf ",minimum);
fprintf(studnum_proj3_1, "Mean: %lf ",mean);
fprintf(studnum_proj3_1, "Sample Standard Deviation: %lf \n",sampleSD);
fprintf(studnum_proj3_1, "Window size: %d ",w);
fprintf(studnum_proj3_1, "Slope: %lf ",m);
fprintf(studnum_proj3_1, "c: %lf \n",c);
        /*prints data into studnum_proj3_1*/
for (i=0;i<inputsize;i++)
    {
        if (i==0)
            {
            fprintf(studnum_proj3_1, "\n | ");
            }
        else
            {
            fprintf(studnum_proj3_1, "|\n | ");
            }
        for (j=0;j<nocolumns;j++)
            {
            fprintf(studnum_proj3_1, "%lf ", detrendmtx[i][j]);
            }
    }
fclose(studnum_proj3_1);
fprintf(studnum_proj3_1, "|\n");
/*open second file for writing */
FILE *studnum_proj3_2;
char fname3[]="1205556_proj3_2.out";
studnum_proj3_2 = fopen(fname3, "w");
/*check file has opened*/
if(studnum_proj3_2 == (FILE*) NULL)
    {
    printf("Can not open file %s\n", fname3);
    exit(EXIT_FAILURE);
    }
for (i=0;i<nobins;i++)
    {
        if (i==0)
            {
            fprintf(studnum_proj3_2, "\n | ");
            }
        else
            {
            fprintf(studnum_proj3_2, "|\n | ");
            }
        for (j=0;j<3;j++)
            {
            if (j==0)
                {
                    fprintf(studnum_proj3_2, "%d samples in range:",
                             (int)histogram[i][j]); /*(code overflows line)*/
                }
            else if (j==1)
                {
                    fprintf(studnum_proj3_2, "%lf to ", histogram[i][j]);
                }
            else {fprintf(studnum_proj3_2, "%lf", histogram[i][j]); }
            }
    }
fprintf(studnum_proj3_2, "|\n");
fclose(studnum_proj3_2);
/****EXTENSION: Create matlab files ***/
/*histogram*/
FILE *studnum_proj3_3;
char fname4[]="histogram.m";
studnum_proj3_3 = fopen(fname4, "w");
/*check file has opened*/
if(studnum_proj3_3 == (FILE*) NULL)
    {
    printf("Can not open file %s\n", fname4);
    exit(EXIT_FAILURE);
    }
fprintf(studnum_proj3_3,"x=%lf:%lf:%lf;\n",minimum+(binsize/2),binsize,
        maximum); /*code overflows line*/
fprintf(studnum_proj3_3, "y = [");
for (i=0;i<nobins-1;i++)
{
    fprintf(studnum_proj3_3, "%d, ", (int)histogram[i][0]);
}
fprintf(studnum_proj3_3, "%d];\n", (int)histogram[nobins-1][0]);
fprintf(studnum_proj3_3, "figure;\n");
fprintf(studnum_proj3_3, "bar(x,y) ;");
fclose(studnum_proj3_3);
/**Detrend visualisation plot**/
FILE *studnum_proj3_4;
char fname5[]="detrend.m";
studnum_proj3_4 = fopen(fname5, "w");
/*check file has opened*/
if(studnum_proj3_4 == (FILE*) NULL)
    {
    printf("Can not open file %s\n", fname5);
    exit(EXIT_FAILURE);
    }
fprintf(studnum_proj3_4, "t = [");
for (i=0;i<inputsize-1;i++)
{
    fprintf(studnum_proj3_4, "%lf, ", mtx[i][0]);
}
fprintf(studnum_proj3_4, "%lf];\n", mtx[inputsize-1][0]);
fprintf(studnum_proj3_4, "x = [");
for (i=0;i<inputsize-1;i++)
{
    fprintf(studnum_proj3_4, "%lf, ", mtx[i][1]);
}
fprintf(studnum_proj3_4, "%lf];\n", mtx[inputsize-1][1]);
fprintf(studnum_proj3_4, "z = [");
for (i=0;i<inputsize-1;i++)
{
    fprintf(studnum_proj3_4, "%lf, ", detrendmtx[i][1]);
}
fprintf(studnum_proj3_4, "%lf];\n", detrendmtx[inputsize-1][1]);
fprintf(studnum_proj3_4, "y = [");
for (i=0;i<inputsize-1;i++)
{
    fprintf(studnum_proj3_4, "%lf, ", c+(m*mtx[i][0]));
}
fprintf(studnum_proj3_4, "%lf];\n", c+(m*mtx[inputsize-1][0]));
fprintf(studnum_proj3_4, "b = [");
for (i=0;i<inputsize-1;i++)
{
    fprintf(studnum_proj3_4, "%lf, ", boxcardata[i/w]);
}
fprintf(studnum_proj3_4, "%lf];\n", boxcardata[(inputsize-1)/w]);
fprintf(studnum_proj3_4, "m = [");
for (i=0;i<inputsize-1;i++)
{
    fprintf(studnum_proj3_4, "%lf, ", mean);
}
fprintf(studnum_proj3_4, "%lf];\n", mean);
fprintf(studnum_proj3_4, "figure;\n");
fprintf(studnum_proj3_4, "plot(t,x,t,z,t,y,t,m,t,b) ;");
fclose(studnum_proj3_4);
fclose(input);
return 0;
}

/******FUNCTIONS*****/

int validateinput(FILE *input)
{
/* checks input of correct form, outputs size of input*/
double temp1=0;
char temp2=0;
int fileend = 0;
int rowcounter=0;
int columncounter=0;
int inputsize=0;
while (fileend == 0)
    {
    temp1=-1;
    fscanf(input, "%lf", &temp1);
    if (temp1==-1) /*fscanf failed:end of file*/
    {
        fileend =1;
               if (rowcounter == 0)
                {
                    rowcounter = 0;
                    inputsize=columncounter;
                    fseek(input, 0, SEEK_SET);/* reset fscanf to beginning */
                    return inputsize;
                }
            else
            {
                printf ("Error in input file: row not correct length");
                exit(EXIT_FAILURE);
            }

    }
    rowcounter++;
    temp2=-1;
    fscanf(input, "%c", &temp2);
    if (temp2==-1) /*fscanf failed:end of file*/
        {
               fileend =1;
               if (rowcounter == nocolumns)
                {
                    rowcounter = 0;
                    columncounter++;
                    inputsize=columncounter;
                    fseek(input, 0, SEEK_SET);/* reset fscanf to beginning */
                    return inputsize;
                }
            else
            {
                printf ("Error in input file: row not correct length");
                exit(EXIT_FAILURE);
            }

        }
    else if (temp2 == '\n')
        {
            if (rowcounter == nocolumns)
                {
                    rowcounter = 0;
                    columncounter++;
                }
            else
            {
                printf ("Error in input file: row not correct length");
                exit(EXIT_FAILURE);
            }
        }

    }
}

void checkwindowsize(int w, int inputsize)
{
    /* gives warning for window sizes that cause many or
    few windows to be produced. Does not end program */
if (w < (inputsize/200))
    {
    printf("Warning: window size may be too small. \n");
    printf("Choose int between %d and %d. \n",(inputsize/200),(inputsize/20));
    }
else if (w > (inputsize/20))
    {
    printf("Warning: window size may be too large. \n");
    printf("Choose int between %d and %d. \n",(inputsize/200),(inputsize/20));
    }
}

double findS_t(int inputsize, double mtx[][nocolumns])
{
/* find S_t */
double S_t=0;
int i;
for (i=0; i<inputsize;i++)
{
    S_t+=(mtx[i][0]);
}
return S_t;
}

double findS_ybar(int w, int nosegments, double boxcardata[],int remainder)
{
/* find S_ybar */
double S_ybar=0;
int i;
for (i=0;i<nosegments-1;i++)
{
    S_ybar=S_ybar+ (w*boxcardata[i]);
}
i=nosegments-1;
S_ybar=S_ybar+ (remainder*boxcardata[i]);
return S_ybar;
}

double findS_tybar(int inputsize, double mtx[][nocolumns],
                    double boxcardata[], int w) /*code overflows line*/
{
/*find S_tybar */
double S_tybar=0;
int i;
for (i=0;i<inputsize;i++)
{
    S_tybar=S_tybar+ (mtx[i][0]*boxcardata[i/w]);
}
return S_tybar;
}

double findS_tt(int inputsize, double mtx[][nocolumns])
{
double S_tt=0;
int i;
for (i=0; i<inputsize;i++)
{
    S_tt=S_tt+mtx[i][0]*mtx[i][0];
}
return S_tt;
}

double findmax(double detrendmtx[][nocolumns],int inputsize)
{
    /*finds maximum of detrended data*/
double currentmax=detrendmtx[0][1];
int i;
for(i=1;i<inputsize;i++)
{
    if( currentmax< detrendmtx[i][1]) {currentmax = detrendmtx[i][1];}
}
return currentmax;
}

double findmin(double detrendmtx[][nocolumns],int inputsize)
{
/*finds minimum of detrended data*/
double currentmin=detrendmtx[0][1];
int i;
for(i=1;i<inputsize;i++)
{
    if( currentmin> detrendmtx[i][1]) {currentmin = detrendmtx[i][1];}
}
return currentmin;
}

double findmean(double detrendmtx[][nocolumns],int inputsize)
{
/*finds mean of detrended data */
double detrendsummation=0;
int i;
for(i=0;i<inputsize;i++)
{
    detrendsummation=detrendsummation+detrendmtx[i][1];
}
double mean=detrendsummation/((double)inputsize);
return mean;
}

double findsampleSD(double detrendmtx[][nocolumns],int inputsize, double mean)
{
/*finds sample standard deviation of detrended data*/
double sumsquareerror=0;
int i;
for (i=0;i<inputsize;i++)
{
sumsquareerror=sumsquareerror+(detrendmtx[i][1]-mean)*(detrendmtx[i][1]-mean);
}
double sampleSD = sqrt((sumsquareerror/(inputsize)));
return sampleSD;
}
