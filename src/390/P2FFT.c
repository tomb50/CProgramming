/****************************************
 ** 390-P3
 ** Tom Beadman
 ** J C
 ** M C
 ** 18 Nov 2011
 ** 
 ** Fast Fourier Transform:
 ** This program calculates the FFT on a supplied spectrum using the implementation
 ** detailed in Numerical Recipes in C 2nd edition.
 ****************************************/

/************************************************************************
 * PREPROCESSOR
 ************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h> /* used in myfft function for size_t type */
#include <errno.h>    
#include <math.h>     /* MUST USE -lm to compile because of math.h library*/
#include <string.h>   
#include <limits.h>
#include <float.h>

/***************************************************
 * TOM BEADMAN */
#define TRUE 1
#define FALSE 0

/***************************************************
 * M C */
#define SPEC_FILEPATH "src/390/spectrum.dat"

typedef struct {
    double *DataX;
    double *DataY;
    unsigned int Pts; /* Number of Points */
    unsigned int LenX; /* Size of array after padding and E2 rounding */
    unsigned int LenY;
} FTArr;

/***************************************************
 * J C */

#define out_fn1 "src/390/reft.txt"            /*name real output file here */
#define out_fn2 "src/390/imft.txt"    	/*name imaginary output file here */


/************************************************************************
 * PROTOTYPES
 ************************************************************************/

/***************************************************
 * J C */
void intro();
void PrintHeader(void);
void PrintFooter(void);
double swapend(float floatin);
double* myfft(double *, size_t, int, size_t *);
FTArr *sintestplot(unsigned int npoints); /*for myfft testing purposes*/
void PrintResult(FTArr S);
void textprint(FTArr);

/***************************************************
 * M C */
char GetEnd();
void ExitFailure(char FuncName[], char ErrorMsg[]);
FTArr *FTACreate(unsigned int Pts);
void FTAPrint(FTArr FTA);
FTArr *FTARead(char Path[]);
unsigned int RoundE2(unsigned int v);

/***************************************************
 * M C + TOM BEADMAN*/
double GetParam(char Name[], char Prompt[], double Min, double Max, char IntFlg);

/***************************************************
 * TOM BEADMAN */
double GaussianValue(double x, double stddev);
double GaussianInverseValue(double y, double std);
FTArr *GaussianPlot(unsigned int npoints, double stddev);

/************************************************************************
 * MAIN - J C + M C
 ************************************************************************/
int main() {
    intro();

    /* Get operation choice - 0 is FT, 1 is inverse FT */
    char OpFlg = (char) GetParam("operation choice",
            "Please enter the number corresponding to the desired operation",
            1, 2, 1) - 1;
    printf("%s selected\n\n", OpFlg ? "Inverse Fourier Transform" : "Fourier Transform");

    /* Get data source choice - 0 is spectrum.dat, 1 is Gaussian */
    char DatFlg = (char) GetParam("data choice",
            "Please enter the number corresponding to the desired data set",
            1, 2, 1) - 1;
    printf("%s selected as data source\n\n", DatFlg ? "Gaussian distribution" : "File spectrum.dat");

    /* Read/generate source data to FTArr struct */
    FTArr *Src;

    if (DatFlg) {
        unsigned int pts = (unsigned int) GetParam("number of points", "", 0, UINT_MAX, 1);
        double sd = GetParam("standard deviation", "", 0, 99999, 0);
        Src = GaussianPlot(pts, sd);
        FTAPrint(*Src);
    } else {
        Src = FTARead(SPEC_FILEPATH);
        FTAPrint(*Src);
    }

    size_t myarray;
    double *FT = myfft(Src->DataY, (size_t) (Src->LenX), 1, &myarray);

    switch (OpFlg + (DatFlg << 1)) { /* Read/generate source data to FTArr struct */
        case 0:
            puts("FT+Spec");
            break;

        case 1:
            puts("Inverse FT+Spec");
            break;

        case 2:
            puts("FT+Gaussian");
            break;

        case 3:
            puts("Inverse FT+Gaussian");
            break;
    }

    return (EXIT_SUCCESS);
}

/************************************************************************
 * J C */

double swapend(float floatin)
/*  4 byte swap function, form inspired by M.Dowsett Px390 lecture course
     university of warwick, lecture 2 slide 6 - J */ {
    unsigned int *floatin_int = (unsigned int *) (&floatin);

    union {
        unsigned int bp1; /*using int's in unions similar to struct's to allow bitwise operations */
        float floatout;
    } floaty;
    /* individual byte manipulation here */
    floaty.bp1 = (*floatin_int >> 24) + ((*floatin_int & 0xFF0000) >> 8) + ((*floatin_int & 0xFF00) << 8) + (*floatin_int << 24);
    return (double) floaty.floatout; /* returns input float as a double via casting */
}

void textprint(FTArr finalarray)
/* function to print peak data to reft.txt and imft.txt files */
/* need to get in numpoints number of entries to be written*/
/* this will need to be the size of the final array??*/
/* in summary the variable finalpeak.Len needs changing or ok?? */
/* names of output files are located in preprocessor */ {
    unsigned int num_points = (finalarray.LenX); /* calculating num of points to be written to file*/
    unsigned int i; /* for looping through input series */
    FILE *outfile; /* output text file pointer, only one needed: reallocated pointer for second file*/


    if ((outfile = fopen(out_fn1, "w")) != NULL) /* open real text file output.txt */ {

        /*write real data to file*/
        for (i = 0; i < (num_points) - 1; i += 2) /*start loop at zero for real elements */ {
            fprintf(outfile, "%g\n\r", finalarray.DataY[i]);
            /* are the rest of the group ok with using %g instead of %f  ? */
        }
        fclose(outfile); /* close reft.txt after writing*/
    } else {
        printf("Cannot write real text file: \t imft.txt. \n");
        exit(0);
    }

    /* second file to be printed now */

    if ((outfile = fopen(out_fn2, "w")) != NULL) /* open imaginary text file output.txt */ {

        /*write imaginary data to file*/
        for (i = 1; i < num_points; i += 2) /*start loop at one for imaginary elements*/ {
            fprintf(outfile, "%g\n\r", finalarray.DataY[i]);
            /* are the rest of the group ok with using %g instead of %f  ? */

        }
        fclose(outfile); /* close imft.txt after writing*/
    } else {
        printf("Cannot write imaginary text file: \t reft.txt. \n");
        exit(0);
    }
}

double* myfft(double paddeddatain[], size_t nn, int isign, size_t myarray[0])
/*

   Decimation in time function adapted from four1 code in
   NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING, (ISBN 0-521-43108-5) 1988-1992, Cambridge University Press,
   William H. Press, Saul A. Teukolsky, William T. Vetterling, Polaroid Corporation, Brian P. Flannery.

   and Elster's algorithm incorporated to handle bit reversed subscripts from 
   "A Parallel Bit Reversal Algorithm and it's CILK Implementation", D�Á�niza C. Morales Berrios, Advisor: Dr. Jaime Seguel, Mathematics Department, University of Puerto Rico
   Mathematics Department, University of Puerto Rico, Mayag�Á�ez Campus, Mayag�Á�ez, Puerto Rico 00681-5000, danizam@cs.uprm.edu 
   ( is one of the top google hits so can be found easily...)

  - - Note: size_t type defined in sys/stat.h library. Group members please do not delete it!

   - - Takes in from left to right: (input data array which must be pre padded with complex elements as zero and padded to next power 2, 
       nn as described above, isign distinguishes between  FFT or FFT^-1 operations to be performed.) 

  - - Inverse is done by reversing sign of exponent. i.e if isign is -1.

  - - Second swap array used when handling bit reversed subscripts is created within the myfft function and eventually this array becomes basis of output.

  - - Haven't used size_t myarray[0]. Could I use this for the bit reversed subscript array? - J note to self
  
Source's description:

  %"Calculates the discrete Fourier transform of datain[0...2*nn], if isign is input as 1; or calculates the
  inverse discrete Fourier transform of datain[1..2*nn], if isign is input as -1. datain must be a complex array 
  of length nn or, equivalently, a real array of length 2*nn. nn MUST be an integer power of 2 
  (this is not checked for!).%"  --note that it refers to "datain" which is "paddeddatain" in this modified version.
  This is to emphasise that the input data must go through the padding process.
  
  Modifications to original NRiC four1 code:
  - -	No need to decrement array for  C's zero offset as indexes replaced with index-1.
  - -	No SWAP as there are both input and output arrays being used in the myfft function. It is inefficient to loop through and swap
          in a single array as shown in fig. 12.2.1 on page 506 of NRiC.
   - -    Rather than using the trigonometric relationships in four1, used "sin" and  "cos" look ups for 
          accuracy as in my opinion they will not make the FFT slow. Code in book is old, lookup is faster now! 
  - -      In addition: instead of standard "sin", "cos"  i'm using GNU extension because calculating both together, according to the 
        GNU C manual, offers a speed advantage. See pg 400 of manual for reference.
  - -     Used M_PI from math library instead of writing out pi to several decimal places as was done in four1.
 */ {
    size_t n, mmax, m, istep, i, l, j;
    double *dataout;
    double wr, wi, ntheta, theta, tempr, tempi;

    static int init = 1;
    static unsigned long *B, *arr_ptr;

    dataout = (double*) malloc(sizeof (double)*2 * nn);


    if (init == 1) /*On first call store bit reversed array subscripts in array*/ { /*doin this instead passing via size_t array to function */
        l = nn / 2;

        arr_ptr = (unsigned long*) malloc(nn * sizeof (unsigned long));

        B = arr_ptr;

        B[0] = 0;

        for (i = 1; i < nn; i = 2 * i, l = l / 2) /*Elster's Algorithm- instead of using SWAP, this handles bitwise swapping*/ {
            for (j = 0; j < i; j++) {
                B[i + j] = B[j] + l;

            }
        }
    }

    for (i = 0; i < nn; i++) /*Place data into bit reversed array subscripts*/ {
        dataout[2 * B[i]] = paddeddatain[2 * i];
        dataout[2 * B[i] + 1] = paddeddatain[2 * i + 1];
    }

    if (init == 3) {
        free(arr_ptr); /*Free on the last call*/
    }
    /*Danielson-Lanczos section beginning here*/
    n = nn << 1;
    mmax = 2;

    while (n > mmax) /*Outer loop executed log2 nn times.*/ {
        istep = mmax << 1;

        theta = isign * ((2 * M_PI) / mmax); /*Initialize the trigonometric recurrence.*/
        /* used pi from math library instead of writing out pi as NRiC had it*/
        ntheta = 0;
        /* wtemp=sin(0.5*theta);
             wpr = -2.0*wtemp*wtemp; */

        wr = 1.0; /*initial "sin" and "cos" values defined here*/
        wi = 0.0;

        for (m = 1; m < mmax; m += 2) /*Here are the two nested inner loops.*/ {
            for (i = m; i <= n; i += istep) {
                j = i + mmax;
                tempr = wr * dataout[j - 1] - wi * dataout[j]; /*All loop indexes decreased by 1*/
                tempi = wr * dataout[j] + wi * dataout[j - 1];
                /*                                                           */
                dataout[j - 1] = dataout[i - 1] - tempr;
                dataout[j] = dataout[i] - tempi;
                dataout[i - 1] += tempr;
                dataout[i] += tempi;
            }

            ntheta += theta;
            fsincos(ntheta, &wi, &wr); /* "sin" and "cos" lookup discussed in function preamble is here */
            /* Trigonometric recurrence no longer here */
        }

        mmax = istep;
    }

    if (isign == -1) /*If isign==-1, divide all points by nn*/ { /* which is for the inverse transform */
        for (i = 0; i < 2 * nn; i += 2) {
            dataout[i] = dataout[i] / nn;
            dataout[i + 1] = dataout[i + 1] / nn;
        }
    }

    init++; /*Count the amount of times function has been called*/

    return dataout;
}

FTArr *sintestplot(unsigned int npoints)
/* Function for generating test array for myfft function - Based on Gaussianplot*/ {

    FTArr *Sinetest = FTACreate(npoints);
    double xlim = (2 * M_PI)*10;
    double x;
    double step = (xlim / npoints)*2; /*Increment in x value*/
    int i = 0; /*Array index & loop count*/
    int j = 0;
    for (x = (-1) * xlim; x <= xlim; x += step) { /*x initialised at the negative xbound*/
        printf("%f,     ", x);
        printf("%i,     ", i);
        Sinetest->DataX[j] = x; /*assign array value*/
        Sinetest->DataY[i] = sin(x); /*assign array value*/
        printf("%.50f,     \n", Sinetest->DataY[i]);
        i += 2;
        j++;
    }
    return (Sinetest); /*returns a pointer to the array*/
}

void intro()
/* basic intro function to console display */ {
    puts("\n\n************************************************************\n"
            "PX390 Project 2 - FFT \n"
            "************************************************************\n\n"

            "Authors: J C, M C, Tom Beadman\n\n"

            "Usage:\n"
            "This program will ask for whether you wish to perform\n"
            "\t1) A Fourier Transform\n"
            "\t2) An inverse Fourier Transform\n"
            "followed by your choice of data to which you wish to apply the operation:\n"
            "\t1) The specific data file spectrum.dat\n"
            "\t2) A Gaussian distribution defined by your choice of 'p' and 'sigma'\n"
            "\t\twhere 'p' is number of points and 'sigma' is standard\n"
            "\t\tdeviation for which you will be asked to specify at runtime.\n\n"

            "Once entered, the program will run and output data files for imaginary\n"
            "and real parts separately for operations performed.\n\n\n"

            "************************************************************\n"
            "EXECUTION \n"
            "************************************************************\n");
}

void PrintResult(FTArr S) {
    printf("\n\n\n************************************************************\n"
            "RESULT \n"
            "************************************************************\n\n"
            "All analysis successfully completed.\n\n");
    FTAPrint(S);
    /*textprint(S); */
    PrintFooter();
}

/* Defines and prints program footer*/
void PrintFooter(void) {
    printf("\n\n\nThe program will now exit.\n"
            "************************************************************\n\n\n\n");
}









/************************************************************************
 * M C */

/* Defines failure routine. Takes a parameter FuncName which corresponds to the
 * name of the function in which the failure occurred and ErrorMsg, the error
 * message to explain the failure. Both are printed to stdout to inform the user. */
void ExitFailure(char FuncName[], char ErrorMsg[]) {
    printf("An exception occured in function %s():\n%s\n\n"
            "The program cannot continue.\n", FuncName, ErrorMsg);
    PrintFooter();
    exit(EXIT_FAILURE);
}

/*------------------------------------------------------------------------------*/

/* Creates a FTArr pointer struct via calloc */
FTArr *FTACreate(unsigned int Pts) {
    FTArr *FTA = (FTArr*) malloc(sizeof (FTArr));
    if (FTA != NULL) {
        FTA->Pts = Pts,
                FTA->DataX = (double*) calloc(FTA->LenX = RoundE2(Pts), sizeof (double)),
                FTA->DataY = (double*) calloc((FTA->LenY = 2 * FTA->LenX), sizeof (double));

        /* If failed to allocate memory for data, free the pointer. This 
           allows the caller to test for (EA == NULL) rather than also having
           to check (EA->Data == NULL) */
        if ((FTA->DataX == NULL) || (FTA->DataY == NULL)) free(FTA), FTA = NULL;
    }
    return FTA;
}

/*------------------------------------------------------------------------------*/

/* Loop through Series and print element by element */
void FTAPrint(FTArr FTA) {
    unsigned int i, j;
    for (i = 0, j = 0; i < FTA.Pts; i++, j += 2)
        printf("%g, %g + %gi\n", FTA.DataX[i], FTA.DataY[j], FTA.DataY[j + 1]);
    printf("%d points in array of length %d.\n\n", FTA.Pts, FTA.LenX);
}

/*------------------------------------------------------------------------------*/

/* Reads a little endian file on a big-endian machine into a Series struct and
 * returns it. Largely lifted from what I coded in previous proj.
 * 
 * Changes: support for FTA struct and added interrogation to endianess */
#define FILE_CHUNK_BYTES 4      /* # of bytes taken up by floats in file */

FTArr *FTARead(char Path[]) {
    FILE *File = fopen(Path, "rb");

    /* Read file and check for NULL pointer */
    if (File != NULL) {

        /* Check file size > 0 and divisible by the size of each (x, y) pair */
        fseek(File, 0, SEEK_END);
        unsigned int Bytes = ftell(File);
        if (Bytes && (Bytes % (2 * FILE_CHUNK_BYTES) == 0)) {
            rewind(File); /* Reset file pointer */

            /* Allocate FTA struct to hold data (of size equal to #{y-points}) */
            FTArr *FTA = FTACreate(Bytes / (2 * FILE_CHUNK_BYTES));
            if (FTA != NULL) {

                /* Loop through file reading in each (x, y) pair */
                char EndFlg = GetEnd();
                unsigned int i, j;
                float FlBuf;
                for (i = 0, j = 0; i < FTA->Pts; i++, j += 2) {
                    /* Read in x value to float buffer, then swap endian and store */
                    fread(&FlBuf, FILE_CHUNK_BYTES, 1, File);
                    FTA->DataX[i] = EndFlg ? ((double) FlBuf) : swapend(FlBuf);

                    /* Do likewise for y value */
                    fread(&FlBuf, FILE_CHUNK_BYTES, 1, File);
                    FTA->DataY[j] = EndFlg ? ((double) FlBuf) : swapend(FlBuf);
                }
            }
            fclose(File);
            return FTA;
        } else
            fclose(File);
    }
    return NULL;
}

FTArr FTAFFT() {

}

/*------------------------------------------------------------------------------*/

/* compute the next highest power of 2 of 32-bit int */
unsigned int RoundE2(unsigned int x) {
    return x--, x |= x >> 1, x |= x >> 2, x |= x >> 4, x |= x >> 8, x |= x >> 16, x++, x += (x == 0);
}

/*------------------------------------------------------------------------------*/

/* Interrogates endianness - from
 * http://www.ibm.com/developerworks/aix/library/au-endianc/index.html?ca=drs-#list5
 * 
 * If return is 1, the running platform is assumed to be little-endian.
 * If it is 0, it is assumed to be big-endian. */
char GetEnd() {
    unsigned int i = 1;
    return ((*(char*) &i) == 0);
}









/************************************************************************
 * M C + TOM BEADMAN */

/*Generic function to get and validate user input, validates against datatypes, max and min restriction*/
/*,invalid characters eg newline etc, informs user to specific issue with validity*/

/*Function takes in the variable name, prompt specificed by the calling function, min-max values and a flag determining if the data is an integer*/
double GetParam(char Name[], char Prompt[], double Min, double Max, char IntFlg) {
    char Buf[BUFSIZ], *EndChar, /* String buffer variables */
            ValidFlg = FALSE, NLFlg, MinFlg, MaxFlg, EndFlg; /* Validation flags */
    double Input;

    /* Check function parameters - revert to defaults if necessary */
    if (!strlen(Name)) Name = "input"; /*If no name specified, default to input*/
    if (Min > Max) { /*Switching Max by Min if Min > Max*/
        double TmpSwp;
        TmpSwp = Min, Min = Max, Max = TmpSwp;
    }

    do {
        /* If prompt specified, print. Else print default prompt. */
        if (strlen(Prompt)) printf(Prompt);
        else printf("Please specify a number between %g and %g for the %s", Min, Max, Name);
        printf(": ");

        if (fgets(Buf, sizeof (Buf), stdin) != NULL) { /*Gets string from user through input stream */
            Input = IntFlg ? ((double) strtol(Buf, &EndChar, 0)) : strtod(Buf, &EndChar); /*Runs String to int function or string to double function depending on Integer flag*/
            ValidFlg = !(NLFlg = (Buf[0] == '\n')) /* Check for initial new line char*/
                    && (MinFlg = (Min <= Input)) && (MaxFlg = (Input <= Max)) /* Check for min and max */
                    && (EndFlg = (*EndChar == '\n' || *EndChar == '\0')); /* Check illegal termination character */

            if (!(ValidFlg || NLFlg)) { /*Output dependant on validity flags*/
                printf("Invalid value specified");
                if (!EndFlg) printf(": the %s must be a %s.", Name, IntFlg ? "integer" : "number"); /*Specific promts to user*/
                else if (!MinFlg) printf(": the %s must be at least %g.", Name, Min);
                else if (!MaxFlg) printf(": the %s must be at most %g.", Name, Max);
                printf(" Please try again.\n\n");
            }
        } else puts("Failed to read from input. Please respecify your input.\n");
    } while (ValidFlg == FALSE);

    return Input; /*returns required value*/
}

/************************************************************************
 * TOM BEADMAN

/*Create an array containing Gaussian distribution, calls GaussianInverse and Gaussian Value functions */
FTArr *GaussianPlot(unsigned int npoints, double stddev) {


    FTArr *Gau = FTACreate(npoints); /*Creates a FTArr struct, of length 2*npoints + necessery value to pad up to power of two*/

    if (Gau != NULL) {
        double xlim = GaussianInverseValue(DBL_MIN, stddev); /*define the interval limit(s) by seeing at which abscisa value the orinate is as close to zero as representable*/
        double x, step = (xlim / (npoints - 1))*2; /*Determine increment in abscissa value*/
        int i = 0; /*Array index & loop count*/
        int j = 0;
        for (x = -xlim; x <= xlim; x += step) { /*abscissa initialised at the negative xbound, increasing by step*/
            Gau->DataX[j] = x; /*Assign abscissa value, not specifically needed in the project but kept for good practice*/
            Gau->DataY[i] = GaussianValue(x, stddev); /*Assign array value from the Gaussian*/
            i += 2; /*Move to next real field*/
            j++;
        }
    }

    return (Gau); /*returns a pointer to the array*/
}

/*This function calculates the value of the Gaussian for a given x and stddev */

/*Called from GaussianPlot to populate the array structure*/
double GaussianValue(double x, double stddev) {

    double a = 1 / (stddev * sqrt(2 * M_PI)); /*Calculate coefficient term*/
    double temp = x / stddev; /*Temporary expression to reduce cpu expense in calculation*/
    double exponent = -0.5 * temp*temp; /*Calculate exponent*/
    double Gval = a * exp(exponent); /*Calculate Value*/
    return (Gval); /*return value*/
}

/*This function calculates the abscissa value from the ordinate argument y for the Gaussian function*/

/*Called from GaussianPlot, used in conjunction with lowest representable floating point value FLT_MIN to determine abscissa limit*/
double GaussianInverseValue(double y, double stddev) {

    double a = 1 / (stddev * sqrt(2 * M_PI)); /*Calculate coefficient term*/
    double yovera = y / a; /*Temporary expression to reduce cpu expense in calculation*/
    double tmp = ((-2) * stddev * stddev); /*Temporary expression*/
    double xsquared = tmp * log(yovera); /*Value for x^2*/
    double x = sqrt(xsquared); /*Calculate x value (C will only take the positive value*/
    return (x); /*Return value*/
}
