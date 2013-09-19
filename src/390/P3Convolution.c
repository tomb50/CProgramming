/****************************************
 ** 390-P3
 ** Tom Beadman
 ** J C
 ** M C
 ** 12 Dec 2011
 ** 
 ** Convolution / Deconvolution:
 ** This program can convolute or de-convolute a supplied spectrum "spectrum.dat"
 ** with a response function - a gaussian distribution.
 ****************************************/



/************************************************************************
 * PREPROCESSOR
 ************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <sys/stat.h> /* used in myfft function is necessary for size_t type - J needs this library!*/
#include <errno.h>    
#include <math.h>     /* MUST USE -lm to compile because of math.h library*/
#include <string.h>   
#include <limits.h>
#include <float.h>

/***************************************************
 * J C */
#define TRUE 1
#define FALSE 0

/***************************************************
 * M C */
#define DECNV_FILT_MAXPROP 0.5  /* Maximum proportional width of filter to response */
#define SCreateDatX(Pts, DatX) SCreateDat(Pts, DatX, (double*)NULL)
#define SCreateDatY(Pts, DatY) SCreateDat(Pts, (double*)NULL, DatY)
#define SCreate(Pts) SCreateDat(Pts, (double*)NULL, (double*)NULL)
#define SDiv(S1, S2) SBinOp(S1, S2, CDiv)
#define SProd(S1, S2) SBinOp(S1, S2, CProd)

/* An fourier transform data structure of {(x), (y-real,y-imag)}
 * including associated metadata required for padding */
typedef struct {
    double *DatX;
    double *DatY;
    unsigned int Pts;           /* Number of Points */
    unsigned int Pad;           /* Padding */
    unsigned int OffDat;
    double OffX;             /* if t != 0, this is the offset */
} Series;

/* A struct to store metadata of response Series*/
typedef struct {
    unsigned char Class;
    unsigned char OOM;
    double RParam;
} RMeta;

/***************************************************
 * TOM BEADMAN */
#define PPSD_MIN 2.0            /*Min points per stddev suggested by M.Dowset in Example Gaussian code*/
#define PPHWHM_MIN 3.0          /*Min points per HWHM for use in Lorentzian, following from HWHM for gaussian
                                 *is aproximattely 1.17741*stddev, and that min points per HWHM logic can be applied
                                 *similarly to the Lorentzian, which would the give 2.35482 minpts per HWHM, which is then rounded up*/
#define POINTS_COEF_MAX 0.5     /*ratio coefficent for total number points*/ 
#define RPARAM_MAX FLT_MAX      /*Safer than dbl max*/
#define OOM_MIN 1               /*Minimum order of magnitude, used when creating response function*/
#define OOM_MAX 8               /*Maximum order of magnitude, used when creating response function*/









/************************************************************************
 * PROTOTYPES
 ************************************************************************/

/***************************************************
 * J C */
double       EndSwp(float floatin);                                                     /* Swaps 4-byte floats, used for endianess conversion */
double*      myfft(double paddeddatain[], size_t nn, int isign, size_t myarray[0]);     /* FFT function*/
void         PrintHead(void);                                                           /* Basic printing of information useful to user on console */
Series*      SConv(Series *S, Series *R, Series **fftdat, Series **fftfunct);           /* Instructions to perform convolution using myfft */
Series*      SConvFT(Series *S, Series *R, Series **fftdat, Series **fftfunct);         /* Partial protion of convolution calculation (blur) */
int          SWrite(Series *S, char Path[], int CFlg);                                  /* Generalised printing function for abscissae and ordinates encountered in this assignment*/

/***************************************************
 * J C + M C */
void         CDiv(double *x, double *y, double *r);
void         CProd(double *x, double *y, double *r);
void         PrintFoot(void);
void         PrintRes(void);
Series*      SBinOp(Series *S1, Series *S2, void (*Op)(double*, double*, double*));
void         SMeta(Series *S, double OffX, unsigned int Pad, unsigned OffDat);


/***************************************************
 * M C */
char         EndGet(void);
void         ExitFail(char FName[], char ErrMsg[]);
double       Max(double x, double y);
double       Min(double x, double y);
unsigned int ModP(int a, unsigned int n);
double       RFWHM(RMeta *RM);
Series*      RGet(char Name[], Series *Src, RMeta *RM);
unsigned int RoundE2(unsigned int v);
Series*      SCreateDat(unsigned int Pts, double *DatX, double *DatY);
Series*      SDeconv(Series *S, Series *R, Series *Filt, Series **FTS, Series **FTR);
double       SDelta(Series *S);
void         SFree(Series **SP);
Series*      SFT(Series *Src, unsigned int PadN, int OffDat);
Series*      SFTI(Series *Src);
Series*      SPad(Series *Src, unsigned int N, int OffDat);
Series*      SPadCpy(Series *Src, Series *Dest);
void         SPrint(Series *S);
Series*      SRead(char Path[]);
Series*      SUnpad(double *DatY, Series *Meta);

/***************************************************
 * M C + TOM BEADMAN */
double       GetParam(char Name[], char Prmpt[], double Min, double Max, unsigned char IntFlg);    /*General function to get numerical input - written in Project2*/

/***************************************************
 * TOM BEADMAN */
double        GauLim(unsigned char OOM);                       /* Calculate abscissa limit for Gaussian based on required range of OOM*/
Series*       GauPlot(Series *spec, RMeta *meta);              /* Creates series of Gaussian distribution*/
double        GauVal(double x, double stddev);                 /* Calculates Gaussian value for given abscissa value*/
double        LorLim(unsigned char OOM);                       /* Calculates abscissa value for Lorentzian based on OOm*/    
Series*       LorPlot(Series *spec, RMeta *meta);              /* Creates series of Lorentzian distribution*/
double        LorVal(double x, double HWHM);                   /* Calculates Lorentzian value for a given abscissa value*/
void          Normalise(Series *resp);                         /* Normalises the elements of a series containing a response function so the area is = 1*/
unsigned char OOM(Series *Arr);                                /* Displays the order of magnitude for a given series, for use on spectrum data to suggest 
                                                               * appropriate OOM range for response function for convolution/deconvolution*/





/************************************************************************
 * MAIN - M C & J C
 ************************************************************************/
int main() {
    PrintHead();
    
    Series *Src, *R, *FTS, *FTR, *Dest;
    /* read in little endian binary data file */
    if ((Src = SRead("src/390/spectrum.dat")) == NULL)
            ExitFail("SRead", "Failed to read from local file. Please ensure "
                              "that the file exists in the working directory.");
    
    
    /* Get operation choice - 0 is Conv 1 is Deconv */
    char *OpDesc[2] = {"Convolution", "Deconvolution"};
    unsigned char OpFlg = (unsigned char)GetParam("operation choice",
                                "Please enter the desired choice of 1)convolution or 2)deconvolution",
                                1, 2, 1) - 1;
    printf("%s selected\n\n", OpDesc[OpFlg]);
    
    /* Get response function */
    RMeta RM;
    char RDesc[256];
    sprintf(RDesc, "%s %s", OpDesc[OpFlg], "response function choice");
    if ((R = RGet(RDesc, Src, &RM)) == NULL)
            ExitFail("SRead", "Failed to read from local file. Please ensure "
                              "that the file exists in the working directory.");
    
    if (OpFlg) { /*perform deconvolution process*/
        RMeta FiltM;
        Series *Filt;
        unsigned char FiltFlg = 0;
        
        /* Read in filtering function and heck that it is not too wide (i.e.
         * wide in comparison to the response so as to make deconvolution ineffective)*/
        do {
            puts("Please choose a smoothing (noise filter):");
            Filt = RGet("deconvolution noise filter", Src, &FiltM);
            if (!(FiltFlg = RFWHM(&FiltM) < DECNV_FILT_MAXPROP*RFWHM(&RM))) {
                puts("The filter specified is too wide in relation to the "
                     "response function for an adequate deconvolution. Please "
                     "reduce the width of the filtering function. See constant "
                     "DECNV_FILT_MAXPROP for more.\n");
            }
        } while (!FiltFlg);
        puts("************************************************************\n\n"
             "Performing deconvolution. Please wait...");
        Dest = SDeconv(Src, R, Filt, &FTS, &FTR);
    }
    else { /* Perform convolution process*/
        puts("************************************************************\n\n"
             "Performing convolution. Please wait...");
        Dest = SConv(Src, R, &FTS, &FTR);   /* not sure if i use fftdat and fftfunct here in these???M!??>*/
    }
    
    /* Check operation is successful */
    if ((FTS == NULL) || (FTR == NULL))
        ExitFail("SRead", "Failed to perform specified operation on data set.");
    
    /* Write to file */
    puts("Writing output to files...");
    if (!(SWrite(Src, "src/390/s.txt", 3) && SWrite(FTS,"src/390/fts.txt", 7) && SWrite(FTR, "src/390/ftr.txt", 7)
            && SWrite(Dest, OpFlg ? "src/390/Deconv.txt" : "src/390/Conv.txt", 3)))
                ExitFail("SWrite", "Failed to write to one or more output files. "
                                   "Please ensure that the disk is not full and "
                                   "that you have sufficient write-access.");
    
    PrintRes();
    return (EXIT_SUCCESS);
}









/************************************************************************
 * TOM BEADMAN
*************************************************************************/

/*Calculates the abscissa limit for a Gaussian distribution given a require range in order of magnitude*/
double GauLim(unsigned char OOM){
    
    double mag=10;
    unsigned int i;
    for(i=0;i<OOM;i++){
        mag=mag*10;
    }
    
    double temp = 2*log(mag);
    double HWXM = sqrt(temp); /*HWXM denotes the width at for a given order of magnitude*/
    return (HWXM);
}

/*------------------------------------------------------------------------------*/

/*This function calculates the value of the Gaussian for a given x and stddev */
/*Called from GaussianPlot to populate the array structure*/
double GauVal(double x, double stddev){
    
    double a=1/(stddev*sqrt(2*M_PI)); /*Calculate coefficient term*/
    double temp = x/stddev; /*Temporary expression to reduce cpu expense in calculation*/
    double exponent = -0.5*temp*temp; /*Calculate exponent*/
    double Gval = a*exp(exponent); /*Calculate Value*/
    return (Gval);    /*return value*/
}

/*------------------------------------------------------------------------------*/

/*Create an array containing Gaussian distribution, calls GaussianInverse and Gaussian Value functions */
/*Suggested improvements based on M. Dowsett's Gaussian example presented in feedback section denoted */
Series *GauPlot(Series *spec, RMeta *meta){

    double stddev = meta->RParam;
    double x_delta = SDelta(spec);
    double xlim = GauLim(meta->OOM)*stddev; /*define the interval limit(s) by seeing at which abscisa value the orinate is as close to zero as representable*/
    unsigned int npoints = floor(2*xlim/x_delta);
    unsigned int invalidOOMflag=FALSE;
    
    if(!(npoints & 1)){ /*Taken from M.Dowsett example Gaussian source code in project2 FAQ feedback page 2011*/ 
        npoints += 1;
    }
        
    /* This check in the generated number of points is greater than the specified limit, if so reduction in the order of magnitude is performed automatically to reduce point count
     * If the order of magnitude is reduced to zero the loop is exited through break and user prompted to reduce HWHM or increase count limit.*/
    while(npoints>spec->Pts*POINTS_COEF_MAX){ /*Loop checks whether or not the number of points generated is greater than the minimum number specified by POINTS_COEF_MAX in preprosessor*/
        meta->OOM--; 
        if(meta->OOM<OOM_MIN){ /*This would imply OOM is being reduced passed a threshold, implicitly OOM=0, therefore user must select smaller HWHM to get sufficent points*/
            invalidOOMflag=TRUE;
            break;
        }
        else{    
            xlim = GauLim(meta->OOM)*stddev;
            npoints = floor(2*xlim/x_delta);
        }
    }
    
    Series *Gau = SCreate(npoints); /*Creates a Series struct*/
    if(invalidOOMflag==FALSE){
    
        double x;
        int i=0; /*Array index & loop count*/
        int j=0;
        x=-xlim;
        for(j=0;j<Gau->Pts/2;j++){ /*Abscissa initialised at the negative xbound, increasing by x_delta from spectrum data*/    
            Gau->DatX[j]=x; /*Assign abscissa value*/ 
            Gau->DatY[i] = GauVal(x,stddev); /*Assign array value from the Gaussian*/
            Gau->DatX[Gau->Pts-1-j]=-x;
            Gau->DatY[2*(Gau->Pts-1)-i]=Gau->DatY[i];/*Symmetry is used to populate the distribution from both ends to split function calling by half*/
            i+=2; /*Move to next real field*/
            x+=x_delta;
        }
        Gau->DatY[i]=GauVal(0,stddev); /*Centroid value assigned*/
    }
    
    else{
        SFree(&Gau); /*Destroys Gau to return NULL for error handling (OOM<OOM_MIN*/
    }
    return(Gau); /*returns a pointer to the array*/
}    

/*------------------------------------------------------------------------------*/

/*Calculates the abscissa limit for a Lorentzian given a required range of order of magnitude*/
double LorLim(unsigned char OOM){
    
    double mag=10;
    unsigned int i;
    for(i=0;i<OOM;i++){ /*This loop calculates the raw value of an order of magnitude ie OOM=3 -> mag=1000*/
        mag=mag*10;        
    }
    double HWXM = sqrt(mag-1); /*Width calculated and returned*/
    return HWXM;
}

/*------------------------------------------------------------------------------*/

/*Create an array containing Lorentzian distribution*/
/*Based on the GaussisnPlot template with appropriate variable change*/
Series *LorPlot(Series *spec, RMeta *meta){
    
    double HWHM = meta->RParam;
    double x_delta = SDelta(spec);
    double xlim = LorLim(meta->OOM)*HWHM; /*define the interval limit*/
    unsigned int npoints = floor(2*xlim/x_delta);
    unsigned int invalidOOMflag=FALSE;
    
    
    if (!(npoints & 1)){ /*Taken from M.Dowsett example Gaussian source code example on project2 faq/feedback page 2011*/ 
        npoints += 1;
    }
    
    /* This check in the generated number of points is greater than the specified limit, if so reduction in the order of magnitude is performed automatically to reduce point count
     * If the order of magnitude is reduced to zero the loop is exited through break and user prompted to reduce HWHM or increase count limit.*/
    
    while(npoints>spec->Pts*POINTS_COEF_MAX){/*Loop checks whether or not the number of points generated is greater than the minimum number specified by POINTS_COEF_MAX in preprosessor*/
        meta->OOM--; 
        if(meta->OOM<OOM_MIN){/*This would imply OOM is being reduced passed a threshold, implicitly OOM=0, therefore user must select smaller HWHM to get sufficent points*/
            invalidOOMflag=TRUE;
            break;
        }
        else{    
            xlim = LorLim(meta->OOM)*HWHM; 
            npoints = floor(2*xlim/x_delta);
        }
    }
    
    Series *Lor = SCreate(npoints); /*Creates a Series struct*/
    if(invalidOOMflag==FALSE){
    
        double x; /*Determine increment in abscissa value*/
        int i=0; /*Array index & loop count*/
        int j=0;

        x=-xlim;
        for(j=0;j<Lor->Pts/2;j++){ /*Abscissa initialised at the negative xbound, increased by x_delta*/    
            Lor->DatX[j]=x; /*Assign abscissa value*/ 
            Lor->DatY[i] = LorVal(x,HWHM); /*Assign array value from the Lorentzian*/
            Lor->DatX[Lor->Pts-1-j]=-x;
            Lor->DatY[2*(Lor->Pts-1)-i]=Lor->DatY[i]; /*Symmetry is used to populate the distribution from both ends to split function calling by half*/
            i+=2; /*Move to next real field*/
            x+=x_delta;
        }
        Lor->DatY[i]=LorVal(0,HWHM); /*Centroid value assigned*/
    }
    else{
        SFree(&Lor);
    }
    return(Lor); /*returns a pointer to the array*/
}

/*------------------------------------------------------------------------------*/

/*This function calculates the value of the Lorentzian for a given x and HWHM */
/*Called from LorentzianPlot to populate the array structure*/
double LorVal(double x, double HWHM){
    
    double tmp = x*x+HWHM*HWHM;
    double denom = M_PI*tmp;
    double Lz;
    if(denom==0 || HWHM==0){ /*If the denominator is zero or HWHM zero*/
        Lz = 0;
    }
    else{    
        Lz=HWHM/denom;
    }
    return (Lz);    /*return value*/
}
 
/*------------------------------------------------------------------------------*/

/*Series to normalise a response function*/
void Normalise(Series *resp){
    double area=0;
    double delta = SDelta(resp); /*Gets x_delta from spectrum data*/
    unsigned int i,j;
    
    for(i=0;i<resp->Pts;i++){
        area+=resp->DatY[2*i]*delta; /*calculates the area of the input series*/   
    }
    double corfact=1/area;  /*Determines a corrective factor to be applied*/
    
    if(area!=1){
        printf("Normalisation being applied...\n");
        for(j=0;j<resp->Pts;j++){
            resp->DatY[2*j]*=corfact; /*scaled each ordinate value by the correction factor*/
        }
    }
}
   
/*------------------------------------------------------------------------------*/

/* This function calculates the order of magnitude for the inputted dataset, measuring the values of Series->DatY
   and taking the maximum and minimum values of the dataset to generate an OOM, Which is then passed onto the Gaussian/Lorentzian 
   functions*/ 
unsigned char OOM(Series *Arr){
    
    double xmax,xmin,maxval,minval,diff; 
    unsigned char oomspec;
    unsigned int position=0,index1,index2;

    maxval=Arr->DatY[0],xmax=Arr->DatX[0], index1=position; /*initializing the max values to the first values of the array*/
    minval=Arr->DatY[0],xmin=Arr->DatX[0], index2=position;/*initializing the min values to the first values of the array*/

    for(position=0;position<Arr->Pts;position++){  /*Main loop, moving through entire length of dataset by step of 2 (due to Real an Im parts for each DatY)*/

        if(Arr->DatY[2*position]>maxval){     /*Basic max filter using condition expression to assign value*/ 
            maxval=Arr->DatY[2*position],index1=position,xmax=Arr->DatX[position];  /*index and corresponding abscissa values are reported (for testing puropses)*/
        }

        if(Arr->DatY[2*position]<minval){ /*Basic max filter using condition expression to assign value*/
             minval=Arr->DatY[2*position],index2=position,xmin=Arr->DatX[position];/*index and corresponding abscissa values are reported (for testing purposes)*/
        }              
    }
    
    if(minval==0){
        diff=fabs(maxval);
    }
    else{
        diff=fabs(maxval/minval);
    }
    oomspec=(unsigned char)floor(fabs(log10(diff)));
    return oomspec;
}









/************************************************************************
 * TOM BEADMAN & M C */
 
 /* Generic function to get and validate numeric user input, validates against datatypes,
  * max and min restriction, invalid characters (eg newline etc), informs user to
  * specific issue with validity. Function takes in the variable name, prompt specified 
  * by the calling function, min-max values and a flag determining if the data is
  * an integer
  * 
  * Params:
  * Name is the string name of param.
  * 
  * Prompt is the string prompt shown to user. If left empty, a default prompt
  * will be used.
  * 
  * Min, Max is bound restrictions.
  *
  * If IntFlg = 0, then will validate for double. If IntFlg = 1, then validates
  * for integers. If IntFlg = 2, then will validate for EVEN integers. If 
  * IntFlg = 3 then will validate for ODD integers.
  */
double GetParam(char Name[], char Prmpt[], double Min, double Max, unsigned char IntFlg) {
    char Buf[BUFSIZ], *EndChar,                                 /* String buffer variables */
         ValidFlg = FALSE, NLFlg, MinFlg, MaxFlg, Parity, ParFlg, EndFlg,       /* Validation flags */
         *Desc[4] = {"a number", "an integer", "an even integer", "an odd integer"};
    double Input;
    
    /* Check function parameters - revert to defaults if necessary */
    if (!strlen(Name)) Name = "input"; /*If no name specified, default to input*/
    if (Min > Max) { /*Switching Max by Min if Min > Max*/
        double TmpSwp;
        TmpSwp = Min, Min = Max, Max = TmpSwp;
    }
    
    do {
        /* If prompt specified, print. Else print default prompt. */
        if (strlen(Prmpt)) printf(Prmpt);
        else printf("Please specify %s between %g and %g for the %s", Desc[IntFlg], Min, Max, Name);
        printf(": ");
        
        /*Gets string from user through input stream */  
        if (fgets(Buf, sizeof(Buf), stdin) != NULL){
            
            /*Runs string to int function or string to double function depending on IntFlg */
            Input = IntFlg ? ((double)strtol(Buf, &EndChar, 0)) : strtod(Buf, &EndChar),
            Parity = (char)((int)Input) & 1,
            
            /* Builds validity flag */
            ValidFlg = !(NLFlg = (Buf[0] == '\n'))                                 /* Check for initial new line char*/
                        && (ParFlg = !(((IntFlg == 2) && Parity) || ((IntFlg == 3) && !Parity))) /* Check parity */
                        && (MinFlg = (Min <= Input)) && (MaxFlg = (Input <= Max))  /* Check for min and max */
                        && (EndFlg = (*EndChar == '\n' || *EndChar == '\0'));      /* Check illegal termination character */
            
            if (!(ValidFlg || NLFlg)) {
                printf("Invalid value specified");
                
                /* Specific warning dependant on validity flags */
                if (!(EndFlg && ParFlg)) printf(": the %s must be %s.", Name, Desc[IntFlg]); /*Specific prompts to user*/
                else if (!MinFlg) printf(": the %s must be at least %g.", Name, Min);
                else if (!MaxFlg) printf(": the %s must be at most %g.", Name, Max);
                printf(" Please try again.\n\n");
            }
        }
        else puts("Failed to read from input. Please respecify your input.\n");
    } while (!ValidFlg);
    
    return Input; /*returns required value*/
}









/************************************************************************
 * M C */ 

/* Interrogates endianness - from
 * http://www.ibm.com/developerworks/aix/library/au-endianc/index.html?ca=drs-#list5
 * 
 * If return is 1, the running platform is assumed to be little-endian.
 * If it is 0, it is assumed to be big-endian. */
char EndGet() {
    unsigned int i = 1;
    return ((*(char*)&i) == 0);
}

/*------------------------------------------------------------------------------*/

/* Defines failure routine. Takes a parameter FName which corresponds to the
 * name of the function in which the failure occurred and ErrorMsg, the error
 * message to explain the failure. Both are printed to stdout to inform the user. */
void ExitFail(char FName[], char ErrMsg[]) {
    printf("An exception occurred in function %s():\n%s\n\n"
           "The program cannot continue.\n", FName, ErrMsg);
    PrintFoot();
    exit(EXIT_FAILURE);
}

/*------------------------------------------------------------------------------*/

/* Finds max of x and y */

double Max(double x, double y) {
    return x > y ? x : y;
}

/*------------------------------------------------------------------------------*/

/* Finds min of x and y */

double Min(double x, double y) {
    return x < y ? x : y;
}

/*------------------------------------------------------------------------------*/

/* Finds smallest positive solution to a modulo n */

unsigned int ModP(int a, unsigned int n) {
    return ((a % n) + n) % n;
}

/*------------------------------------------------------------------------------*/

/* Find FWHM of a response function from its metadata.
 * FWHM = 2*sqrt(2*ln(2))*SD */
double RFWHM(RMeta *RM) {
    return RM->Class ? RM->RParam : 2*sqrt(2*log(2))*RM->RParam;
}


/*------------------------------------------------------------------------------*/
/* Asks for user input to populate a response Series struct of either Gaussian
 * or Lorentzian. Name is the full name of R, Src is the corresponding Series data,
 * RMeta is the associated metadata struct for the response function
 */

Series *RGet(char Name[], Series *Src, RMeta *RM) {
    Series *R = NULL;
    
    /* Get choice of R - either 0 for Gau or 1 for Lor */
    static char *RDesc[2] = {"Gaussian", "Lorentzian"};
    RM->Class = (unsigned char)GetParam(Name,
                                        "Please enter the desired choice of 1)Gaussian or 2)Lorentzian",
                                        1, 2, 1) - 1;
    printf("%s selected\n", RDesc[RM->Class]);
    
    /* Get corresponding parameter (either SD or FWHM) and OOM. If the combination
     * of function parameter and OOM produces an illegal R Series, then ask the 
     * use again for different values */
    char *ParDesc[2] = {"Gaussian SD", "Lorentzian FWHM"},
          RParamPrmpt[BUFSIZ], OOMPrmpt[BUFSIZ];
    sprintf(RParamPrmpt, "%s %s", "Please enter the desired", ParDesc[RM->Class]),
    sprintf(OOMPrmpt, "%s %u%s", "Please enter the desired OOM (the approx. "
                                 "OOM for specified data is", OOM(Src), ")");
    
    do {
        RM->RParam = GetParam("response function parameter choice", RParamPrmpt,
                              DBL_EPSILON, RPARAM_MAX, 0);
        
        /* Get Order of magnitude user input for response function */
        RM->OOM = (unsigned char)GetParam("response function order of magnitude (OOM)",
                                          OOMPrmpt, OOM_MIN, OOM_MAX, 1);
        printf("\n");
        
        /* Populate plot function pointer depending on RFlg */
        Series* (*RPlot)(Series*, RMeta*) = RM->Class ? LorPlot : GauPlot;
        if ((R = RPlot(Src, RM)) == NULL) {   /* Notify if plot failed and loop back */
            puts("Failed to create response function. Please ensure that "
                "the OOM and/or HWHM are not too large as to cause the response "
                "function to be extremely long or short in length (please see the "
                "POINTS_COEF_MAX constant) for more.\n");
        }
    } while (R == NULL);
    
    return R;
}

/*------------------------------------------------------------------------------*/

/* Compute the next highest power of 2 of 32-bit int from 
 * http://graphics.stanford.edu/~seander/bithacks.html#RoundUpPowerOf2
 * 
 * Essentially, sets all all bits to right of the first bit equal to 1, to 1.
 */

unsigned int RoundE2(unsigned int x) {
     return x |= --x >> 1, x |= x >> 2, x |= x >> 4, x |= x >> 8, x |= x >> 16, x += (++x == 0);
}


/*------------------------------------------------------------------------------*/

/* Creates a Series struct pointer using DatX and DatY as the data arrays. If 
 * DatX or DatY are missing (i.e. passed as NULL), calloc is used to create the
 * required array. */

Series *SCreateDat(unsigned int Pts, double *DatX, double *DatY) {
    Series *S = (Series*)malloc(sizeof(Series));
    if (S != NULL) {
        S->Pts = Pts,
        S->DatX = (DatX == NULL) ? (double*)calloc(S->Pts, sizeof(double)) : DatX,
        S->DatY = (DatY == NULL) ? (double*)calloc(2*S->Pts, sizeof(double)) : DatY;
        
        /* If failed to allocate memory for data, free the pointer. */
        if ((S->DatX == NULL) || S->DatY == NULL) SFree(&S);
        else SMeta(S, 0, 0, 0);
    }
    return S;
}

/*------------------------------------------------------------------------------*/

/* Deconvolution of S using response R. Smooths by filter Filt then divides by
 * the fourier transform of the response to obtain a frouier transform of the
 * deconvolved Series. Finally performs an inverse fourier transform. */

Series *SDeconv(Series *S, Series *R, Series *Filt, Series **FTS, Series **FTR) {
    double Q = (R->Pts - 1)/2;  /* The # of points at end of file */
    
    Series *FTFilt, *SmthFT = SConvFT(S, Filt, FTS, &FTFilt);
    SFree(&FTFilt);
    if (SmthFT != NULL) {
        *FTR = SFT(R, S->Pts - R->Pts + Q + 1, -Q);
        if (*FTR != NULL) {
            Series *DcnvFT = SDiv(SmthFT, *FTR);
            if (DcnvFT != NULL) {
                SFree(&SmthFT);
                return SFTI(DcnvFT);
            }
            SFree(FTR);
        }
        SFree(&SmthFT);
    }
    return NULL;
}

/*------------------------------------------------------------------------------*/

/* Calculates the value of delta (x-value spacing) in a Series containing equally
 * spaced x-values */

double SDelta(Series *S) {
    return (S->DatX[S->Pts - 1] - S->DatX[0])/(S->Pts - 1);
}

/*------------------------------------------------------------------------------*/

/* Frees a Series Struct by freeing DatX and DatY and finally of the struct
 * itself. Notice a pointer to the Series pointer is required in order to set
 * the Series pointer to NULL. */

void SFree(Series **SP) {
    Series *S = *SP;
    if (S->DatX != NULL) free(S->DatX);
    if (S->DatY != NULL) free(S->DatY);
    free(S), *SP = NULL;
}

/*------------------------------------------------------------------------------*/

/* Wrapper for myfft to integrate it with Series struct. This will also deal with
 * x-values in the struct so scatter plots can be created even though only y-values
 * are required for output.
 * 
 * Notice LenX, LenY is assumed to be a power of 2 and hence even. Also, it is 
 * assumed that the x-values are equally spaced 
 * 
 * Some general comments about myfft wrapper functions:
 * Note that padding/unpadding was designed to be performed locally within myfft
 * wrapper functions so as to make the whole padding process invisible to and
 * independent from other functions - all other functions only need to work with
 * unpadded arrays, and hence don't have to worry whether to take into account
 * padding at centre on at the end etc.
 * 
 * SFT and SFTI are separate functions as it eliminates having to pass a switch
 * parameter and acts as syntax sugar.
 */

Series *SFT(Series *Src, unsigned int PadN, int OffDat) {
    Series *FT = SPad(Src, PadN, OffDat);
    if (FT != NULL) {
        /* Create fft y array and replace existing DatY */
        size_t *Dummy;
        double *DatY = FT->DatY;        /* Cache pointer to free later */
        FT->DatY = myfft(FT->DatY, (size_t)(FT->Pts), 1, Dummy),
        free(DatY);

        /* Calculate x-values.
         * d = delta, sn = partial sum cache - faster than to keep multiplying.
         * Notice Pts is assumed to be a power of 2 and hence even */
        double d = SDelta(Src), nd = 1/(d*FT->Pts), sn = 0;
        unsigned int i, j=FT->Pts-1, u = FT->Pts/2;

        /* Calculate x-values from both ends of the array for efficiency */
        for(i = 1; i < u; i++, j--) {
            FT->DatX[i] = (sn += nd),
            FT->DatX[j] = -FT->DatX[i];
        }
        FT->DatX[0] = 0, FT->DatX[u] = 1/(2*d);
    }
    return FT;
}

/*------------------------------------------------------------------------------*/

/* Wrapper for myfft for inverse FT to integrate myfft with Series struct. The
 * function calculates and unpads the inverse FT, storing in a Series struct.
 * To regenerate the abscissa values it calculates delta (interval between samples)
 * and progressively sums to rebuild all abscissa values.
 * 
 * Also note that the 1/(2*delta) frequency value is not used to calculate delta
 * in case it's set to -1/2d manually by user 
 * 
 * For more, please also see the general comments about myfft wrapper functions
 * in the comments to function SFT().
 */

Series *SFTI(Series *Src) {
    /* First, perform inverse FT and unpad */
    size_t *Dummy;
    Series *IFT = SUnpad(myfft(Src->DatY, (size_t)(Src->Pts), -1, Dummy), Src);
    if (IFT != NULL) {
        /* Regenerate abscissae values
         * d = delta, sn = partial sum cache - faster than to keep multiplying. */
        double d = (1/Src->DatX[1])/(Src->Pts), sn = IFT->DatX[0] = Src->OffX;
        unsigned int i;
        for(i = 1; i < IFT->Pts; i++)
            IFT->DatX[i] = (sn += d);
    }
    return IFT;
}

/*------------------------------------------------------------------------------*/

/* A generalised padding function: pads the Src Series S by an extra N, then to 
 * the next power of 2. OffDat gives the offset of the data from the beginning 
 * of the array in the padded array so formats such as wrap-around can be achieved.
 * 
 * Negative OffDat offsets it counting backwards from the end of the file instead.
 */

Series *SPad(Series *Src, unsigned int N, int OffDat) {
    Series *Dest = SCreate(RoundE2(Src->Pts + N));      /* Create padded array */
    if (Dest != NULL) {
        /* Calculate metadata */
        Dest->Pad = Dest->Pts - Src->Pts;
        Dest->OffDat = ModP(OffDat, Dest->Pts);  /* Find lowest positive offset equivalent */
        Dest->OffX = Src->OffX;
        SPadCpy(Src, Dest);       /* Transfer data from Src to Dest */
    }
    return Dest;
}

/*------------------------------------------------------------------------------*/

/* Copies y-values from a source Series to a destination Series, taking into
 * account padding metadata in both arrays. Thus essentially, it transfers data
 * from Src Series (either padded or unpadded) to Dest Series (of opposite padding
 * type). memcpy() is used in lieu of a loop for speed (as it is written in optimised
 * assembly code) and arrays in C always occupy contiguous space in memory. See
 * http://stackoverflow.com/questions/1696074/how-can-i-concatenate-two-arrays-in-c
 * for more.
 * 
 * The transfer in done in 2 stages, corresponding to the 2 possible places the data
 * be located in a padded array (i.e. if in wrap-around format). */

Series *SPadCpy(Series *Src, Series *Dest) {
    
    /* Determine which is the padded and unpadded Series */
    Series *Pad, *Upad;
    if (Src->Pad) Pad = Src, Upad = Dest;
    else          Pad = Dest, Upad = Src;
    
    /* Calculate limits for padded and unpadded series. Define (dummy) offsets for
     * Src and Dest offsets (Src may be padded or unpadded and vice/versa so the
     * dummy variables will be linked to the true offsets later).
     * 
     * Variable naming: U/P for Unpadded/Padded, O for Offset, 1/2 for copy stage,
     * N for number of points to copy, S for Src and D for Dest */
    unsigned int PO1 = Pad->OffDat, N1 = (unsigned int)Min((double)Upad->Pts, (double)(Pad->Pts - Pad->OffDat)), /* Stage 1 calculations */
                 UO2 = N1, N2 = Upad->Pts - N1,          /* Stage 2 calculations */
                 SO1 = 0, DO1 = 0, SO2 = 0, DO2 = 0;     /* Src/Dest dummy data offsets */
    
    /* Link up Src/Dest offsets with Padded/unpadded offsets depending on which
     * one of Src/Dest is padded */
    if (Src->Pad) SO1 = PO1, DO2 = UO2;
    else          DO1 = PO1, SO2 = UO2;
    
    /* Copy first chunk of data (from the padded Series offset to as far as
     * either the padded or unpadded Series will allow) */
    memcpy(&Dest->DatX[DO1], &Src->DatX[SO1], N1 *= sizeof(double));
    memcpy(&Dest->DatY[2*DO1], &Src->DatY[2*SO1], 2*N1);
    
    /* If unpadded Series still has points not copied, wrap around the padded
     * Series and keep copying. */
    if (N2) {
        memcpy(&Dest->DatX[DO2], &Src->DatX[SO2], N2 *= sizeof(double));
        memcpy(&Dest->DatY[2*DO2], &Src->DatY[2*SO2], 2*N2);
    }
    
    return Dest;
}

/*------------------------------------------------------------------------------*/

/* Loop through Series and print element by element */
void SPrint(Series *S) {
    unsigned int i, j;
    for(i = 0, j = 0; i < S->Pts; i++, j += 2)
        printf("%g, %g + %gi\n", S->DatX[i], S->DatY[j], S->DatY[j + 1]);
    printf("%d points in array.\n\n", S->Pts);
}

/*------------------------------------------------------------------------------*/

/* Reads a little endian file from Path[] into a Series struct and
 * returns it. */
#define FILE_CHUNK_BYTES 4      /* # of bytes taken up by floats in file */
Series *SRead(char Path[]) {
    
    /* Read file and check for NULL pointer */
    FILE *File = fopen(Path, "rb");
    if (File != NULL) {
        
        /* Check file size > 0 and divisible by the size of each (x, y) pair */
        fseek(File, 0, SEEK_END);
        unsigned int Bytes = ftell(File);
        if (Bytes && (Bytes % (2*FILE_CHUNK_BYTES) == 0)) {
            rewind(File);       /* Reset file pointer */
            
            /* Allocate S struct to hold data (of size equal to #{(x, y) pairs) */
            Series *S = SCreate(Bytes/(2*FILE_CHUNK_BYTES));
            if (S != NULL) {
                
                /* Loop through file reading in each (x, y) pair */
                char EndFlg = EndGet();
                unsigned int i, j;
                float FlBuf;
                for (i = 0, j = 0; i < S->Pts; i++, j += 2) {
                    
                    /* Read in x value to float buffer, then swap endian and store */
                    fread(&FlBuf, FILE_CHUNK_BYTES, 1, File),
                    S->DatX[i] = EndFlg ? EndSwp(FlBuf) : ((double)FlBuf),
                    
                    /* Do likewise for y value */
                    fread(&FlBuf, FILE_CHUNK_BYTES, 1, File),
                    S->DatY[j] = EndFlg ? EndSwp(FlBuf) : ((double)FlBuf);
                }
                S->OffX = S->DatX[0];
            }
            fclose(File);
            return S;
        }
        else
            fclose(File);
    }
    return NULL;
}

/*------------------------------------------------------------------------------*/

/* Unpads a padded double array DatY and stores in a Series Struct. Uses
 * underlying metadata of Meta Series to base the new Series struct upon. Frees
 * the padded input array upon success */
Series *SUnpad(double *DatY, Series *Meta) {
    
    /* Create source and destination arrays for SPadTrans */
    Series *Src = SCreateDatY(Meta->Pts, DatY),
           *Dest = SCreate(Meta->Pts - Meta->Pad);
    
    /* If Series failed to create, ensure both are destroyed, else send it off
     * to SPadTrans to transfer the data */
    char SrcN = (Src == NULL), DestN = (Dest == NULL);
    if (SrcN || DestN) {
        if (!SrcN)  SFree(&Src);
        if (!DestN) SFree(&Dest);
    }
    else
        Dest->OffX = Src->OffX, SPadCpy(Src, Dest), SFree(&Src);    /* Also frees DatY */
    
    return Dest;
}









/************************************************************************
 * M C & J C */ 

/* Function to perform complex division of two complex inputs */
void CDiv(double *x, double *y, double *r) {
    double d = y[0]*y[0] + y[1]*y[1];
    r[0] = ((x[0]*y[0]) + (x[1]*y[1]))/d,
    r[1] = ((x[1]*y[0]) - (x[0]*y[1]))/d; 
}

/*------------------------------------------------------------------------------*/

/* Function to perform complex multiplication */
void CProd(double *x, double *y, double *r) {
    r[0] = ((x[0]*y[0]) - (x[1]*y[1])),
    r[1] = ((x[0]*y[1]) + (x[1]*y[0]));
}

/*------------------------------------------------------------------------------*/

/* Defines and prints program footer*/
void PrintFoot(void) {
   puts("The program will now exit.\n\n"
        "************************************************************\n\n\n");
}
/*------------------------------------------------------------------------------*/

/* Prints result UI */
void PrintRes(void) {
    puts("\n\n\n************************************************************\n"
         "RESULT\n"
         "************************************************************\n\n"
         "All analysis successfully completed.");
    PrintFoot();
}

/*------------------------------------------------------------------------------*/

/* Function to perform point-by-point binary operation on two input Series struts.
 * Return NULL if NaN or Inf detected*/
Series *SBinOp(Series *S1, Series *S2, void (*Op)(double*, double*, double*)) {
    
    Series *DomS = (S1->Pts > S2->Pts) ? S2 : S1,       /* Find dominant Series (i.e. takes precedence with metatdata)*/
           *Dest = SCreateDatX(DomS->Pts, DomS->DatX);
    
    SMeta(Dest, DomS->OffX, DomS->Pad, DomS->OffDat);
    
    unsigned int i;
    for (i = 0; i < Dest->Pts; i++) { /*point-by-point mathematics performed here */
        Op(&S1->DatY[i], &S2->DatY[i], &Dest->DatY[i]);
    
        /* NaN Inf detection */
        if (isnan(Dest->DatY[i]) || isinf(Dest->DatY[i])
                || isnan(Dest->DatY[++i]) || isinf(Dest->DatY[i])) {
                    SFree(&Dest);
                    break;
        }
    }
    return Dest;
}

/*------------------------------------------------------------------------------*/

/* Sets the metadata of Series struct S */
void SMeta(Series *S, double OffX, unsigned int Pad, unsigned OffDat) {
    S->OffX = OffX, S->Pad = Pad, S->OffDat = OffDat;
}









/************************************************************************
 * J C */

double EndSwp(float floatin) 
/*  (As written in project 2) 4 byte swap function, form inspired by M.Dowsett Px390 lecture course
     university of warwick, lecture 2 slide 6 - J */
{
    unsigned int *floatin_int = (unsigned int *)(&floatin);
    union{
        unsigned int bp1; /*using int's in unions similar to struct's to allow bitwise operations */
        float floatout;                                              
    } floaty;
    /* individual byte manipulation here */
    floaty.bp1 = (*floatin_int >> 24) + ((*floatin_int & 0xFF0000) >> 8) + ((*floatin_int & 0xFF00) << 8) + (*floatin_int << 24);
    return (double)floaty.floatout; /* returns input float as a double via casting */
}

/*------------------------------------------------------------------------------*/

double* myfft(double paddeddatain[], size_t nn, int isign, size_t myarray[0])
/*
   (AS I WRITTEN FOR PROJECT 2 BUT MODIFIED HEAVILY FOR CONVOLUTION AND DECONVOLUTION PURPOSES)
   
   Decimation in time function adapted from four1 code in
   NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING, (ISBN 0-521-43108-5) 1988-1992, Cambridge University Press,
   William H. Press, Saul A. Teukolsky, William T. Vetterling, Polaroid Corporation, Brian P. Flannery.

   and Elster's algorithm incorporated to handle bit reversed subscripts from 
   "A Parallel Bit Reversal Algorithm and it's CILK Implementation", Dániza C. Morales Berrios, Advisor: Dr. Jaime Seguel, Mathematics Department, University of Puerto Rico
   Mathematics Department, University of Puerto Rico, Mayagüez Campus, Mayagüez, Puerto Rico 00681-5000, danizam@cs.uprm.edu 
   ( is one of the top google hits so can be found easily...)

- - Takes in from left to right: (input data array which must be pre padded with complex elements as zero and padded to next power 2, 
       nn as described above, isign distinguishes between  FFT or FFT^-1 operations to be performed.) 

  - - Inverse is done by reversing sign of exponent. i.e if isign is -1.

  - - Second swap array used when handling bit reversed subscripts is created within the myfft function and eventually this array becomes basis of output.
   
Source's description:

  %"Calculates the discrete Fourier transform of datain[0...2*nn], if isign is input as 1; or calculates the
  inverse discrete Fourier transform of datain[1..2*nn], if isign is input as -1. datain must be a complex array 
  of length nn or, equivalently, a real array of length 2*nn. nn MUST be an integer power of 2 
  (this is not checked for!).%"  --note that it refers to "datain" which is "paddeddatain" in this modified version.
  This is to emphasise that the input data must go through the padding process.
  
  Modifications to original NRiC four1 code:
  - -    No need to decrement array for  C's zero offset as indexes replaced with index-1.
  - -    No SWAP as there are both input and output arrays being used in the myfft function. It is inefficient to loop through and swap
          in a single array as shown in fig. 12.2.1 on page 506 of NRiC.
   - -    Rather than using the trigonometric relationships in four1, used "sin" and  "cos" look ups for 
          accuracy as in my opinion they will not make the FFT slow. Code in book is old, lookup is faster now! 
  - -     In addition: instead of standard "sin", "cos"  i'm using GNU extension because calculating both together, according to the 
    GNU C manual, offers a speed advantage. See pg 400 of manual for reference.
  - -     Used M_PI from math library instead of writing out pi to several decimal places as was done in four1.
  - -     The bit reversed subscript section is executed only upon the first call OR when the function is detected to be
     using a new input array length. i.e bit reversed subscripts are calculated for the first call to the function for each array size.
     This is for convolution and deconvolution purposes. Callcounter at end of function keeps count of how many times
     function is called and stores in a static variable so that the information is retained throughout entire
     duration of program run. 

  - -     I have recently become aware that Dr. Dowsett dislikes global variables but I could not see another way of utilising counters in the
           manner I have.
  */
{
    size_t n,mmax,m,istep,i, l, j;
    double *dataout;
    double wr,wi,ntheta,theta,tempr,tempi; 
    
    unsigned char truthflag;                  /*Used to state whether input array size is a new length or not to previous function call*/
    static int callcounter=1; 
    static size_t initialarraysize;
    size_t newarraysize;    
    
    static unsigned long *B,*arr_ptr;   
    
    if(callcounter==1) /* Set initial value of initialarraysize with first call to myfft */
    {
        initialarraysize=nn;   
        truthflag=TRUE;    
    }
    else
    {
        newarraysize=nn;     /*Used to check to see if array being performed on with this call of myfft function */
                                /* is the same length as the previous time myfft function was called. If different length  */
                                /* is detected then bit reversed subscripts must be recalculated for new array length.  */
                                /* If array is same length as previous call then subscripts can remain unchanged. */    
        if(initialarraysize==newarraysize)
        {
            truthflag=FALSE; /* If new input array is of same length as previous call then bit reversed subscript section of code will not execute */
        }
        else
        {
            truthflag=TRUE;             /* If new input array length differs from previous call then bit reversed subscript code must be executed */
            initialarraysize=newarraysize;       /* Load new array length into memory as being the initial array length to be checked against for the*/
            free(arr_ptr);                        /*next time the function is called.*/
        }
    }
    dataout=(double*)malloc(sizeof(double)*2*nn);
        
    if(truthflag)                /*On first call store bit reversed array subscripts in array*/
    {                        /*doing this instead passing via size_t array to function */
        l=nn/2;                              /*Also if new input array size is detected then this section of code will execute as truthflag is set to true */
        arr_ptr=(unsigned long*)malloc(nn*(sizeof(unsigned long))*2);    
        B=arr_ptr;    
        B[0]=0;
        for(i=1; i < nn ;i=2*i, l=l/2)                 /*Elster's Algorithm- instead of using SWAP, this handles bitwise swapping*/
        {                                 
            for( j=0; j<i; j++)                
            {                            
                B[i+j] = B[j] + l;                           
            }
        }
    }                                    
    for(i=0; i<nn; i++)                        /*Place data into bit reversed array subscripts*/
    {
        dataout[2*B[i]]=paddeddatain[2*i];
        dataout[2*B[i]+1]=paddeddatain[2*i+1];
    }
    
    /*Danielson-Lanczos section beginning here*/
    n=nn << 1;                                    
    mmax=2;
    while (n > mmax)                             /*Outer loop executed log2 nn times.*/
    { 
        istep=mmax << 1;                                             
        theta=isign*((2*M_PI)/mmax);     /*Initialize the trigonometric recurrence.*/
                                /* used pi from math library instead of writing out pi as NRiC had it*/
        ntheta=0;
        /* wtemp=sin(0.5*theta);
        wpr = -2.0*wtemp*wtemp; */
        wr=1.0; /*initial "sin" and "cos" values defined here*/
        wi=0.0; 
        for (m=1;m<mmax;m+=2)                     /*Here are the two nested inner loops.*/
        { 
            for (i=m;i<=n;i+=istep) 
            {
                j=i+mmax;                         
                tempr=wr*dataout[j-1]-wi*dataout[j];    /*All loop indexes decreased by 1*/
                tempi=wr*dataout[j]+wi*dataout[j-1];
                /*                                                           */
                dataout[j-1]=dataout[i-1]-tempr;
                dataout[j]=dataout[i]-tempi;
                dataout[i-1] += tempr;
                dataout[i] += tempi;
            }
            ntheta+=theta;            
            sincos (ntheta, &wi, &wr);     /* "sin" and "cos" lookup discussed in function preamble is here */
                                                     /* Trigonometric recurrence no longer here */
        }
        mmax=istep;
    }
    if(isign==-1)                                /*If isign==-1, divide all points by nn*/
    {                                                                              /* which is for the inverse transform.*/
        for(i=0; i < 2*nn ; i+=2)
        {
            dataout[i]=dataout[i]/nn;
            dataout[i+1]=dataout[i+1]/nn;
        }
    }
    
    callcounter++;                        /*Count the amount of times function has been called*/
    return dataout;
}

/*------------------------------------------------------------------------------*/

void PrintHead(void)
/* basic intro function to console display */
{
    puts("\n\n************************************************************\n"
         "PX390 Project 3 - FFT Convolution Deconvolution \n"
         "************************************************************\n\n"

         "Authors: Tom Beadman, M C, J C\n\n"
       
         "Usage:\n"
         "This program will convolve/deconvolve data with a specified response\n"
         "function from the binary file spectrum.dat. For your choice of operation,\n"
         "you may select from the following:\n"
             "\t1) Convolve the data\n"
             "\t2) Deconvolve the data\n\n"
             
         "For the response function you may select from the following:\n"
             "\t1) A Gaussian\n"
             "\t2) A Lorentzian\n\n"
            
         "Once entered, you will also be asked for the SD/FWHM and order of\n"
         "magnitude (OOM) depending on your choice of response function.\n"
         "For deconvolution, you will also be prompted for a smoothing (noise\n"
         "filter) function.\n\n"   
         
         "Output data will be in the following files\n"
             "\ts.txt\t\tRaw data from spectrum.dat\n"
             "\tftr.txt\t\tFT of response function\n"
             "\tfts.txt\t\tFT of spectrum.dat\n"
             "\tconv/deconv.txt\tResulting convolution/deconvolution\n\n\n\n"

         "************************************************************\n"
         "EXECUTION \n"
         "************************************************************\n");
}

/*------------------------------------------------------------------------------*/

/*  Function to produce output text files.
 *
 *  - - Names of output file are taken in when calling the function. This function will serve as a generic print
 *   function to increase versatility. f_name must be entered including quotation marks.
 *  - - Toggleparam settings: Enter as number below to perform the listed task when printing:
 *
 *  0: Does nothing                  4: imaginary ordinates only
 *  1: abscissae only                5: abscissae and imaginary ordinates
 *  2: real ordinates only           6: both imaginary and real ordinates
 *  3: abscissae and real ordinates  7: abscissae, real and imaginary ordinates        
 */

int SWrite(Series *S, char Path[], int CFlg) {
    FILE *File = fopen(Path, "w");       /* open input data text file */
    unsigned char RetFlg = (File != NULL) && CFlg;
    if (RetFlg) {
        unsigned int i, j = 0, k = 1;
        unsigned char C1Flg = CFlg & 1, C2Flg = CFlg & 2, C3Flg = CFlg & 4;
        for (i = 0; i < S->Pts; i++, j += 2, k += 2) {
            if (C1Flg) fprintf(File, "%g%s", S->DatX[i], C2Flg || C3Flg ? " " : "");
            if (C2Flg) fprintf(File, "%g%s", S->DatY[j], C3Flg ? " " : "");
            if (C3Flg) fprintf(File, "%g", S->DatY[k]);
            fprintf(File, "\r\n");
        }                
        fclose(File);
    }    
    return RetFlg;
}

/*------------------------------------------------------------------------------*/

Series *SConv(Series *S, Series *R, Series **fftdat, Series **fftfunct)
  /* Function inspired by Lecture 8, page 13 of Px390 2011 lectured by Mark Dowsett, Univeristy of Warwick.
       Function finds FFT of data and also the response (This is done in within SConvFT function). It then multiplies them together point-by-point (using SProd function) 
       and then takes the FFT^-1 of the subsequent combined array. 
  
      Input two arrays which are complex (pad complex elements to zero prior to entering this
      function if needed) and of equal length and padded to the nearest power of 2.  Returns
      a pointer to another array of the convoluted pair.
  
  - - S is the bit reversed array (after various padding also) created from the data
      file read in.      
  - - R is the array generated by choice of user of the Gaussian or Lorentzian functions
      present in this assignment. i.e the response function.      
   */  
  
{
    Series *temp = SConvFT(S,R,fftdat,fftfunct);
    
    if(temp!=NULL){
        Series *convout=SFTI(temp);
        return convout;
    }
    return NULL;      /*Function returns a pointer to the convolution array*/
} 

/*------------------------------------------------------------------------------*/
      
Series *SConvFT(Series *S, Series *R, Series **fftdat, Series **fftfunct) 
  /* This function is a small subsection of code used in convolution process present in myconv function, however it was taken 
     out to be a separate function as it is also utilised for blurring purposes in deconvolution function.*/
{
    *fftdat=SFT(S,(unsigned int)((R->Pts+1)/2),0);                                                /*Calculate fft of input data file*/
    
    if(*fftdat!=NULL)
    {
        *fftfunct=SFT(R,(unsigned int)(S->Pts-R->Pts+(R->Pts+1)/2),(int)(-(R->Pts-1)/2));    /*Calculating fft of Gaussian or Lorentzian */
        
        if(*fftfunct!=NULL)
        {
            Series *temp = SProd(*fftdat, *fftfunct);
            
            if(temp!=NULL)
            {
                return temp;
            }
        }
        SFree(fftfunct);
    }
    SFree(fftdat);
    return NULL;
}
