/****************************************
 ** 390-P1
 ** Tom Beadman
 ** J C
 ** M C
 ** 1 Nov 2011
 ** 
 ** Peakfinder:
 ** This program find the peaks in the spectrum data provided
 ** by file 'spectrum.dat'. The user is prompted for parameters 
 ** 'm' (which is currently not used), and 'n' which governs
 ** the splitting of the distribution.
 ****************************************/

#include <stdio.h>
#include <stdlib.h>

#define FILE_CHUNK_BYTES 4
#define true 1
#define false 0


typedef struct
{
    double x;
    double y;
} Point;

typedef struct
{
    Point *Data;
    unsigned int Len; /* Number of Points */
} Series;


/************************************************************************
 * PROTOTYPES
 ************************************************************************/
void intro(); /*J*/
void SPrint(Series This); /*M*/
Series *NPMA(Series Src, unsigned int n); /*M*/
Series *SRead(char Path[]); /*J*/
Series *SCreate(unsigned int Len); /*M*/
double SwapEnd(float Src); /*J*/
unsigned int *SPeaks(Series Src); /*M and Tom*/
char Sgn(double x); /*M*/
int getm(void); /*Tom*/
unsigned int getn(unsigned int noofx); /*Tom*/
Series *Mean(unsigned int array[], int m, Series *MA, Series *Src); /*Tom*/
unsigned int *Filter(unsigned int darray[], Series ma, Series src); /*Tom*/

/************************************************************************
 * CORE FUNCTIONS
 ************************************************************************/
int main(int argc, char** argv)
{


    void intro();
    int m = getm();
    Series *Src = SRead("src/390/spectrum.dat");
    Series *MA = NPMA(*Src, getn(Src->Len));
    unsigned int *dum = SPeaks(*MA);
    unsigned int *dum2 = Filter(dum, *MA, *Src);
    // Series *ME = Mean(dum2,m,MA,Src);
    printf("--------MA-------\n\n\n\n\n\n");
    SPrint(*Src);
    printf("--------Src-------\n\n\n\n\n\n");
    SPrint(*MA);
    printf("-------Speaks-----\n\n\n\n\n\n");
    int i;
    for (i = 0; i < 10; i++)
        printf("dum(%i)=%i, dum2(%i)=%i\n", i, dum[i], i, dum2[i]);

    return (EXIT_SUCCESS);
}
/**** J FUNCTIONS *****/

/**** M FUNCTIONS*****/

/**** TOM FUNCTIONS ****/

/* Reads a little endian file on a big-endian machine into a Series struct and
 * returns it */
Series *SRead(char Path[])
{
    FILE *File;
    /* Read file and check for NULL pointer */
    if ((File = fopen(Path, "rb")) != NULL) {

        /* Check file size > 0 and divisible by the size of each (x, y) pair */
        fseek(File, 0, SEEK_END);
        unsigned int Bytes = ftell(File);
        if (Bytes && (Bytes % (2 * FILE_CHUNK_BYTES) == 0)) {
            rewind(File); /* Reset file pointer */

            /* Allocate Series struct to hold data (of size equal to #{(x, y) points}) */
            Series *Src = SCreate(Bytes / (2 * FILE_CHUNK_BYTES));
            if (Src != NULL) {
                /* Loop through file reading in each (x, y) pair */
                unsigned int i;
                float FlBuf;
                for (i = 0; i < Src->Len; i++) {
                    /* Read in x value to float buffer, then swap endian and store */
                    fread(&FlBuf, FILE_CHUNK_BYTES, 1, File);
                    Src->Data[i].x = SwapEnd(FlBuf);

                    /* Do likewise for y value */
                    fread(&FlBuf, FILE_CHUNK_BYTES, 1, File);
                    Src->Data[i].y = SwapEnd(FlBuf);
                }
            }
            fclose(File);
            return Src;
        }
        else
            fclose(File);
    }
}

/* Returns an n-point running average Series of Src - notice it may be of a
 * different length to Src
 *
 * The numerator sum of the running average is constructed by summing in the 
 * new y-term and subtracting the previous y-term. That is, for a previous
 * sum of y[i] + ... + y[i + n - 1], the new sum y[i + 1] + ... + y[i + n]
 * is constructed as follows:
 * 
 * y[i] + ... + y[i + n - 1] +   y[i + n]    -   y[i]
 * |---- previous sum -----| + |-new term-| - |old term|
 * 
 * The sum is then divided by n for the n-average.
 */
Series *NPMA(Series Src, unsigned int n)
{
    if (Src.Len && (n & 1) && (Src.Len >= n)) { /* n must be odd to preserve peaks */

        /* Allocate return Series of appropriate size */
        Series *MA = SCreate(Src.Len - n + 1);
        if (MA != NULL) {

            /* i = index to first y-val in current running average,
             * j = index to first y-val in next running average
             * k = x index*/
            unsigned int i, j = n - 1, k = j / 2; /* n - 1 since 0-based index */
            double s = 0; /* Holds the partial sums of running averages */

            /* Build a partial sum to start */
            for (i = 0; i < n - 1; i++)
                s += Src.Data[i].y;

            for (i = 0; i < MA->Len; i++, j++, k++) {
                MA->Data[i].x = k; /* Re-index of x */
                s += Src.Data[j].y; /* Sum the next y-value */
                MA->Data[i].y = s / n; /* Calculate and store the average*/
                s -= Src.Data[i].y; /* Drop the old obsolete y-value */
            }
        }
        return MA;
    }
}

/* Creates a Series pointer struct via malloc */
Series *SCreate(unsigned int Len)
{
    Series *This = (Series*) malloc(sizeof (Series));
    if (This != NULL) {
        This->Data = (Point*) malloc((This->Len = Len) * sizeof (Point));
        if (This->Data == NULL)
            /* Failed to allocate memory for data so free the pointer s. This 
             allows the caller to test for (This == NULL) rather than also having
             to check (This->Data == NULL) */
            free(This);
    }
    return This;
}

unsigned int *SPeaks(Series Src)
{
    /* s = sign flag (sign of y[i + 1] - y[i])
     * g = "gradient" flag (1 if d >= 0) */
    char s = -1;
    unsigned char g = 0, *Flags = (unsigned char*) malloc(Src.Len * sizeof (unsigned char));

    unsigned int i, c = 0; /* c = count */
    if (Flags != NULL) {
        for (i = 0; i < Src.Len - 2; i++) {
            if (Flags[i] = (((s = Sgn(Src.Data[i + 1].y - Src.Data[i].y)) <= 0) && g))
                c++;
            g = (s >= 0);
        }

        Flags[Src.Len - 1] = (Src.Data[Src.Len - 1].y == Src.Data[Src.Len - 2].y);

        unsigned int *Idx = (unsigned int*) malloc(c * sizeof (unsigned int));
        if (Idx != NULL) {
            for (i = 0, c = 0; i < Src.Len - 1; i++) {
                if (Flags[i])
                    Idx[c++] = (unsigned int) Src.Data[i].x;
            }
            return Idx;
        }
    }
}

char Sgn(double x)
{
    return x ? (x > 0 ? 1 : -1) : 0;
}

void SPrint(Series S)
{
    int i;
    printf("Printing Array of length %d\n", S.Len);
    for (i = 0; i < S.Len; i++)
        printf("%lf %lf\n", S.Data[i].x, S.Data[i].y);
    printf("\n");
}

double SwapEnd(float Src)
{
    unsigned int *SrcInt = (unsigned int *) (&Src);

    union
    {
        unsigned int Int;
        float Float;
    } Dest;
    Dest.Int = (*SrcInt >> 24) + ((*SrcInt & 0xFF0000) >> 8) + ((*SrcInt & 0xFF00) << 8) + (*SrcInt << 24);
    return (double) Dest.Float;
}

void intro()
{
    puts("-------------------------------------------");
    puts("-------------------------------------------");
    puts("-          PX390 - P1   PEAKFINDER            -");
    puts("- will ask for m and n values at a later stage -");
    puts("--------------------------------------------");
    puts("--------------------------------------------");
}

/*Function to get m value from user input- TOM*/
int getm(void)
{

    int maxmval = 99999, mt;
    char buf1[BUFSIZ], *p1;
    int input_completem = false;


    while (input_completem == false) {
        printf("Please enter a sensible, positive, non-zero even value for the variable m:");

        if (fgets(buf1, sizeof (buf1), stdin) != NULL) {
            mt = strtol(buf1, &p1, 10); /* Use of strtol assisted by http://www.mkssoftware.com/docs/man3/strtol.3.asp */

            /*Checks for initial newline character, non numerical characters and appropriate total value */
            if (buf1[0] != '\n' && (*p1 == '\n' || *p1 == '\0') && (mt > 0 && mt < maxmval && (mt & 1) == 0)) {
                printf("Valid number of %ld entered\n\n", mt);
                input_completem = true;
            }
            else {
                printf("\n***Invalid number entered***\n");
                printf("***Restarting...***\n\n");
            }
        }
    }
    return (mt);
}

/*Function to get n value from input, upper limit passed in as integar- TOM*/
unsigned int getn(unsigned int noofx)
{ /*noofx should be set to Src->Len in main()*/

    unsigned int maxnval = 0.1 * noofx, nt;
    char buf2[BUFSIZ];
    char *p2;
    int input_completen = false;
    while (input_completen == false) {
        printf("Please enter a sensible, positive, non-zero odd value for the variable n:");
        if (fgets(buf2, sizeof (buf2), stdin) != NULL) {
            nt = strtol(buf2, &p2, 10); /*string to integer function*/
            /*Checks for initial newline character, non numerical characters and appropriate total value*/
            if (buf2[0] != '\n' && (*p2 == '\n' || *p2 == '\0') && (nt > 0 && nt < maxnval && (nt & 1) != 0)) {
                printf("Valid number of %ld entered\n", nt);
                input_completen = true;
            }
            else {
                printf("\n***Invalid number entered for variable 'm'***\n");
                printf("***Restarting...***\n\n");
            }
        }
    }
    return (nt);
}

Series *Mean(unsigned int array[], int m, Series *MA, Series *Src)
{

    unsigned int arraylength = sizeof (array) / sizeof (unsigned int);

    /*Create Series of length equal to that of the inputted filter array*/
    Series *ME = SCreate(arraylength); /*Divide by 4 because array contains int type*/
    int i, j, temp; /*Temporary variables for loops*/
    int count; /*Loop counter variable*/
    if (ME != NULL) {
        for (i = 0; i < ME->Len; i++) {
            count = 0; /* Initializing count value*/
            temp = array[i]; /*Takes the value from input array*/
            ME->Data[i].x = Src->Data[temp].x; /*Takes the X value from original Series*/
            j = temp - (m / 2); /*Initializes j value the lower limit*/

            while ((j - m / 2) <= temp <= (j + m / 2)) { /* Cumulative addition of y values from 2nd Series*/
                ME->Data[i].y = ME->Data[i].y + MA->Data[j].y;
                j++;
                count++;
            }
            ME->Data[i].y = ME->Data[i].y / count; /*Division of total value by number of values to get Mean*/
        }
        return (ME); /*Return Series*/
    }
    else {
        printf("Unable to create mean value series, filtered array may be NULL\n");
    }
}

unsigned int *Filter(unsigned int darray[], Series ma, Series src)
{/*Filter - Tom*/

    int i = 0;
    unsigned int tp; /*Temp variables*/
    int n = 10; /* use for span either side of testing point, 20 chosen for this dataset*/
    int length = (sizeof (darray) / sizeof (unsigned int));
    unsigned int *farray = (unsigned int *) malloc(sizeof (darray)); /*Allocate memory for new Array*/

    while (i < length) {
        tp = (unsigned int) darray[i];
        if (tp < 10)
            farray[i] = tp;

        i++;

    }
    return farray;
}

