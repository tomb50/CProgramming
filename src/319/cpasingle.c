/****************************************
 ** Tom Beadman
 ** 25 Feb 2012
 ** 
 ** Count per Area (central density calculator):
 ** This program works with astronomical object field datasets in csv format 
 ** provided by the SDSS catalog, a single example input file is provided 
 **  "QSO0441I.csv". The column of interest is the last column which specifies
 ** the radial distance from the center of the field. From this a linear 
 ** regression is calculated for the relative object density against distance,
 ** confirming whether or not there is an over/under-density at INNERDIS compared
 ** with the value at OUTERDIS. An additional alternate value is calulated by
 ** following the same process but omitting the linear regression, taking the 
 ** required densities from the raw values. 
 ****************************************//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "nrutil.h"
#define ITMAX 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30
#define INNERDIS 0.5
#define MIDDIS 1.0
#define OUTERDIS 2.5
#define TRUE 1
#define FALSE 0

float fitval(float a, float b, float x);
void rerror(char error_text[]);
float gammln(float xx);
void gcf(float *gammcf, float a, float x, float *gln);
void gser(float *gamser, float a, float x, float *gln);
void rerror(char error_text[]);
void fit(float x[], float y[], int ndata, float sig[], int mwt, float *a, float *b, float *siga, float *sigb, float *chi2, float *q);
float gammq(float a, float x);
float interpol(float a, float fa, float b, float fb, float c);
void binpop(double dist[], int count[], float rad[], int nobjects, float bwidth, float upperlim);
void areacalc(float rad[], int nbins, float binwidth, float scalefactor, float centre[], float area[]);
void countperparea(int nbins, float cpa[], int count[], float area[], float error[]);
int closestindex(float pdist[], float value, int nbins);
void subset(float pdis[], float cpa[], float err[], float subpdis[], float subcpa[], float suberr[], int cutoff, int originalnbins);

int main()
{


    /*Wrapper*/

    /*--------------------------------------------------*/
    float scalefactor = 0.3408; /*z=0.44*/
    float bwidth = 1.0;
    float upperlim = 15.0;
    int nbins = (int) upperlim / bwidth;
    int* count = calloc(nbins, sizeof (int));
    float* rad = calloc(nbins, sizeof (float));
    float* pdis = calloc(nbins, sizeof (float));
    float* area = calloc(nbins, sizeof (float));
    float* cpa = calloc(nbins, sizeof (float));
    float* err = calloc(nbins, sizeof (float));

    /*fit variables*/
    int mwt = TRUE;
    float* fitgrad = malloc(sizeof (float));
    float* fitintercept = malloc(sizeof (float));
    float* fituncertgrad = malloc(sizeof (float));
    float* fituncertintercept = malloc(sizeof (float));
    float* fitchisquared = malloc(sizeof (float));
    float* fitgoodness = malloc(sizeof (float));

    char inp_fn[50] = "src/319/QSO044I.csv";
    // printf("Enter field file name:\n");
    // scanf("%s", &inp_fn);
    // printf("Set scalefactor(arcsec/mpc):\n");
    // scanf("%f",&scalefactor);

    FILE *input = fopen(inp_fn, "r");

    int noobj = 0;
    char check[550];
    char temp[550];
    int i = 0;
    int counter = 0, header = 0;
    char headertemp1[220];

    if ((input != (FILE*) NULL)) {


        while (fgets(temp, sizeof (temp), input) != NULL) {
            fscanf(input, "%s", check);
            if (check != '\n') noobj++;
        }
        rewind(input); /*returns to beginning of file*/

        double* distance = calloc(noobj, sizeof (double)); /*array created to store distance*/
        double* tempcol = calloc(noobj, sizeof (double));

        while (!feof(input)) {
            if (header == 0) {
                fscanf(input, "%s%*c \n", &headertemp1[header]);
                header++;
            }
            if (header != 0) {
                fscanf(input, "%lf%*c %lf%*c", &tempcol[counter], &tempcol[counter]);
                fscanf(input, "%lf%*c %lf%*c", &tempcol[counter], &tempcol[counter]);
                fscanf(input, "%lf%*c %lf%*c", &tempcol[counter], &tempcol[counter]);
                fscanf(input, "%lf%*c %lf%*c", &tempcol[counter], &tempcol[counter]);
                fscanf(input, "%lf%*c %lf%*c", &tempcol[counter], &tempcol[counter]);
                fscanf(input, "%lf%*c %lf%*c", &tempcol[counter], &tempcol[counter]);
                fscanf(input, "%lf%*c %lf%*c", &tempcol[counter], &tempcol[counter]);
                fscanf(input, "%lf%*c %lf%*c", &tempcol[counter], &tempcol[counter]);
                fscanf(input, "%lf%*c %lf%*c", &tempcol[counter], &tempcol[counter]);
                fscanf(input, "%lf%*c %lf%*c", &tempcol[counter], &tempcol[counter]);
                fscanf(input, "%lf%*c %lf%*c", &tempcol[counter], &tempcol[counter]);
                fscanf(input, "%lf%*c %lf%*c", &tempcol[counter], &tempcol[counter]);
                fscanf(input, "%lf%*c %lf%*c", &tempcol[counter], &tempcol[counter]);
                fscanf(input, "%lf%*c %lf%*c", &tempcol[counter], &tempcol[counter]);
                fscanf(input, "%lf%*c %lf*%c\n", &tempcol[counter], &distance[counter]);
                counter++;
            }
        }
        noobj--;
        fclose(input);
        float maxdis = 0;

        for (i = 0; i < noobj; i++) {
            if (distance[i] > maxdis) maxdis = distance[i];
        }


        binpop(distance, count, rad, noobj, bwidth, upperlim); /*Populate Bins*/

        areacalc(rad, nbins, bwidth, scalefactor, pdis, area); /*Calculate Area*/

        countperparea(nbins, cpa, count, area, err); /*Calculate CPA and error*/

        //    printf("\ncpa started:\n");        
        //    for(i=0;i<nbins;i++){
        //        printf("i=%i, rad[]=%f, count[]=%i, centre[]=%f, area=%f, cpa=%f, err=%f\n",i,rad[i],count[i],pdis[i],area[i],cpa[i],err[i]);
        //    }    

        int iinner = closestindex(pdis, INNERDIS, nbins); /*Determines index for greatest pdis less than INNERIS*/
        int iouter = closestindex(pdis, OUTERDIS, nbins); /*Determines index for greatest pdis less than OUTERDIS*/

        //    printf("\nclosest index = %i, outerindex = %i\n", iinner, iouter);
        //    printf("innerldis:%f, innerlcpa:%f\n",pdis[iinner],cpa[iinner]);
        //    printf("innerudis:%f, innerucpa:%f\n",pdis[iinner+1],cpa[iinner+1]);
        //    printf("outerldis:%f, outerlcpa:%f\n",pdis[iouter],cpa[iouter]);
        //    printf("outerudis:%f, outerucpa:%f\n",pdis[iouter+1],cpa[iouter+1]);

        float innercpa = interpol(pdis[iinner], cpa[iinner], pdis[iinner + 1], cpa[iinner + 1], INNERDIS); /*Interpolation for estimate of CPA at INNERDIS*/
        float outercpa2 = interpol(pdis[iouter], cpa[iouter], pdis[iouter + 1], cpa[iouter + 1], OUTERDIS); /*Interploation for estimate of CPA at OUTERDIS*/

        //   printf("\n innercpa = %f, outercpa = %f\n", innercpa, outercpa2);

        int icut = closestindex(pdis, MIDDIS, nbins); /*Determines index for greatest pdis less than MIDDIS*/
        int subbins = nbins - icut - 2; /*CNumber of bins for outer subset, misses of outermost bin*/
        float* subpdis = calloc(subbins, sizeof (float));
        float* subcpa = calloc(subbins, sizeof (float));
        float* suberr = calloc(subbins, sizeof (float));

        subset(pdis, cpa, err, subpdis, subcpa, suberr, icut, nbins); /*Populates outer subset for Linear Correlation fit*/
        //   printf("\n Subset:\n");

        /*   int j;
        for (j=0;j<subbins;j++){
            printf("j = %i, subpdis = %f, subcpa = %f, suberr=%f\n",j,subpdis[j],subcpa[j],suberr[j]);         
        }
         */

        /*Fit function*/
        fit(subpdis, subcpa, subbins, suberr, mwt, fitintercept, fitgrad, fituncertintercept, fituncertgrad, fitchisquared, fitgoodness);

        //   printf("\nfitting process for subset\n");
        //   printf("gradient=%f, intercept=%f, chi2 =%f, goodness=%f\n",*fitgrad,*fitintercept,*fitchisquared,*fitgoodness);


        float outercpa1 = fitval(*fitgrad, *fitintercept, OUTERDIS);
        float density1 = innercpa / outercpa1;
        float density2 = innercpa / outercpa2;

        //   printf("outercpaval = %f\n", outercpa1);

        printf("%s,density(LC),%f,density,%f\n", inp_fn, density1, density2);
    }

    return (EXIT_SUCCESS);
}

void subset(float pdis[], float cpa[], float err[], float subpdis[], float subcpa[], float suberr[], int cutoff, int originalnbins)
{
    /*Populates a subset array from a given data array, subset is left bounded by index cutoff up to the last bin, which is excluded*/

    int subbins = originalnbins - cutoff - 2;
    int j;
    for (j = 0; j < subbins; j++) {
        subpdis[j] = pdis[cutoff + 1];
        subcpa[j] = cpa[cutoff + 1];
        suberr[j] = err[cutoff + 1];
        cutoff++;
    }
}

int closestindex(float pdist[], float value, int nbins)
{
    /*returns index for which the data is the greatest value less than "value" set in second argument*/

    int i, cindex = 0;
    float temp = 1;
    for (i = 0; i < nbins; i++) {
        if (((value - pdist[i]) < temp)&&((value - pdist[i]) >= 0)) cindex = i;
    }
    return cindex;
}

void binpop(double dist[], int count[], float rad[], int nobjects, float bwidth, float upperlim)
{
    /*Populates bins given count and radius and upper lim*/

    int i = 0, j = 0, k = 0;
    /*printf("start:binpop\n");
    printf("nobjects: %i, bwidth: %f, upperlim: %f\n\n",nobjects,bwidth,upperlim);
     */
    /*upper radius annulus calculation*/
    for (j = 0; j < upperlim / bwidth; j++) {
        rad[j] = (j + 1) * bwidth;
    }
    /*population of count*/
    for (i = 0; i < nobjects; i++) {
        for (k = 0; k < upperlim / bwidth; k++) {
            if (k == 0) {
                if (dist[i] <= rad[k]) count[k]++;
            }
            else {
                if (dist[i] <= rad[k] && dist[i] > rad[k - 1]) count[k]++;
            }
        }
    }
}

void areacalc(float rad[], int nbins, float binwidth, float scalefactor, float centre[], float area[])
{
    /*Calculates area in Mpc^2 given scalefactor and radius, assumes scalefactor is adjusted for arcsecs/Mpc
     *Also converts distance into Mpc*/

    int i;
    for (i = 0; i < nbins; i++) {
        centre[i] = rad[i] - binwidth / 2;
        centre[i] *= scalefactor;
        if (i == 0) area[i] = M_PI * rad[i] * rad[i];
        else area[i] = M_PI * (rad[i] * rad[i] - rad[i - 1] * rad[i - 1]);
        area[i] *= scalefactor*scalefactor;
    }
}

void countperparea(int nbins, float cpa[], int count[], float area[], float error[])
{
    /*Calculates count per area and associated error given count and area*/
    int i;
    //printf("\n\nCPA:\n");
    for (i = 0; i < nbins; i++) {
        cpa[i] = count[i] / area[i];
        error[i] = sqrt(count[i]) / area[i];
        //printf("i=%i, count:%i, are:%f, cpa: %f err: %f\n",i,count[i],area[i],cpa[i],error[i]);
    }
}

void fit(float x[], float y[], int ndata, float sig[], int mwt, float *a, float *b, float *siga, float *sigb, float *chi2, float *q)
{
    /*NRIC - Fitting function takes abscissa values, function values with associated uncertainty (mwt flag !=0 => no uncertainty) returns gradient, intercept
     * associated abscissa and function uncertainties, chi^2 value and goodness of fit value page - page 665*/

    int i;
    float wt, t, sxoss, sx = 0.0, sy = 0.0, st2 = 0.0, ss, sigdat;

    *b = 0.0;
    if (mwt) {
        ss = 0.0;
        for (i = 0; i < ndata; i++) {
            if (sig[i] == 0) wt = 1;
            else wt = 1.0 / SQR(sig[i]);
            ss += wt;
            sx += x[i] * wt;
            sy += y[i] * wt;
        }

    }
    else {
        for (i = 0; i < ndata; i++) {
            sx += x[i];
            sy += y[i];
        }
        ss = ndata;
    }
    sxoss = sx / ss;
    if (mwt) {
        for (i = 0; i < ndata; i++) {
            if (sig[i] == 0) t = x[i] - sxoss;
            else t = (x[i] - sxoss) / sig[i];
            st2 += t*t;
            if (sig[i] == 0) *b += t * y[i];
            else *b += t * y[i] / sig[i];
            /*printf("grad %f\n", *b);*/
        }

    }
    else {
        for (i = 0; i < ndata; i++) {
            t = x[i] - sxoss;
            st2 += t*t;
            *b += t * y[i];
            /*printf("grad %f\n", *b);*/
        }
    }
    *b /= st2; /*grad here*/
    *a = (sy - sx * (*b)) / ss; /*cept here*/
    *siga = sqrt((1.0 + sx * sx / (ss * st2)) / ss);
    *sigb = sqrt(1.0 / st2);

    *chi2 = 0.0;
    *q = 1.0;
    if (mwt == 0) {
        for (i = 0; i < ndata; i++)
            *chi2 += SQR(y[i]-(*a)-(*b) * x[i]);
        sigdat = sqrt((*chi2) / (ndata - 2));
        *siga *= sigdat;
        *sigb *= sigdat;
    }
    else {
        for (i = 0; i < ndata; i++)
            if (sig[i] == 0) *chi2 += SQR(y[i]-(*a)-(*b) * x[i]);
            else *chi2 += SQR((y[i]-(*a)-(*b) * x[i]) / sig[i]);
        if (ndata > 2) *q = gammq(0.5 * (ndata - 2), 0.5 * (*chi2));
    }
}

float gammq(float a, float x)
{
    /*NRIC - return the incomplete Gamma Function Q(a,x)= 1-P(a,x) - page 218*/

    float gamser, gammcf, gln;

    if (x < 0.0 || a <= 0.0) rerror("Invalid arguments in routine gammq");
    if (x < (a + 1.0)) {
        gser(&gamser, a, x, &gln);
        return 1.0 - gamser;
    }
    else {
        gcf(&gammcf, a, x, &gln);
        return gammcf;
    }
}

void gser(float *gamser, float a, float x, float *gln)
{
    /*NRIC - Returns incomplete Gamma function as a series representation - page 218*/

    int n;
    float sum, del, ap;

    *gln = gammln(a);
    if (x <= 0.0) {
        if (x < 0.0) rerror("x less than 0 in routine gser");
        *gamser = 0.0;
        return;
    }
    else {
        ap = a;
        del = sum = 1.0 / a;
        for (n = 1; n <= ITMAX; n++) {
            ++ap;
            del *= x / ap;
            sum += del;
            if (fabs(del) < fabs(sum) * EPS) {
                *gamser = sum * exp(-x + a * log(x)-(*gln));
                return;
            }
        }
        rerror("a too large, ITMAX too small in routine gser");
        return;
    }
}

void gcf(float *gammcf, float a, float x, float *gln)
{
    /*NRIC - Returns incomplete Gamma function Q(a,x) evaluated by continued fraction representation - page 219*/

    int i;
    float an, b, c, d, del, h;

    *gln = gammln(a);
    b = x + 1.0 - a;
    c = 1.0 / FPMIN;
    d = 1.0 / b;
    h = d;
    for (i = 1; i <= ITMAX; i++) {
        an = -i * (i - a);
        b += 2.0;
        d = an * d + b;
        if (fabs(d) < FPMIN) d = FPMIN;
        c = b + an / c;
        if (fabs(c) < FPMIN) c = FPMIN;
        d = 1.0 / d;
        del = d*c;
        h *= del;
        if (fabs(del - 1.0) < EPS) break;
    }
    if (i > ITMAX) rerror("a too large, ITMAX too small on gcf");
    *gammcf = exp(-x + a * log(x)-(*gln)) * h;
}

float gammln(float xx)
{
    /*NRIC  - Logarithm of Gamma function - page 214*/

    double x, y, tmp, ser;
    static double cof[6] = {76.18009172947146, -86.50532032941677,
        24.01409824083091, -1.231739572450155,
        0.1208650973866179e-2, -0.5395239384953e-5};
    int j;

    y = x = xx;
    tmp = x + 5.5;
    tmp -= (x + 0.5) * log(tmp);
    ser = 1.000000000190015;
    for (j = 0; j <= 5; j++) ser += cof[j] / ++y;
    return -tmp + log(2.5066282746310005 * ser / x);
}

float interpol(float a, float fa, float b, float fb, float c)
{
    /*Interpolates a function value between to abscissa points - error thrown if extrapolation is implied*/

    if (a > c || b < c) rerror("Unsuitable points passed for interpolation in module: interpval");
    float grad = (fb - fa) / (b - a);
    float val = fa + (c - a) * grad;
    return val;
}

float fitval(float a, float b, float x)
{
    /*Calculates straight line function from gradient, intercept and abscissa*/

    float val = a * x + b;
    return val;
}

void rerror(char error_text[])
{
    /* Error handler */

    fprintf(stderr, "Run-time error...\n");
    fprintf(stderr, "%s\n", error_text);
    fprintf(stderr, "...now exiting to system...\n");
    exit(EXIT_FAILURE);
}

