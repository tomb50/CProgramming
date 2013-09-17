/****************************************
 ** Tom Beadman
 ** 16 Feb 2012
 ** 
 ** Kolmogorov-Smirnov Test:
 ** This program is a wrapper for the kstwo test implementation detailed
 ** in Numerical Recipes in C second edition. The program compares datasets in
 ** qso.txt and grb.txt files with that contained in blank.txt. 
 ****************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define EPS1 0.001
#define EPS2 1.0e-8

void kstwo(float data1[], unsigned long n1, float data2[], unsigned long n2, float *d, float *prob);
void sort(unsigned long n, float arr[]);
float probks(float alam);
int getval(void);

int main()
{

    FILE *input1; /*blank fields in separate file*/
    FILE *input2; /*Subject fields e.g. qso's or grbs*/
    int i = 0, j = 0;
    unsigned long nblank = 0, nqso = 0;
    const char inp_fn1[] = "src/319/blank.txt";
    const char inp_fn2[] = "src/319/grb.txt";
    const char inp_fn3[] = "src/319/qso.txt";
    float blanks[100], qsos[100];
    float* x = malloc(sizeof (float));
    float* y = malloc(sizeof (float));
    printf("\n\n**************************************\n");
    printf("Program to determine K-S statistics for two given datasets\n");
    printf("Blank fields should be kept in blank.txt in the source folder\n");
    printf("Qso objects to be kept in qso.txt, grbs to be kept in grb.txt in the same location\n");
    printf("The program sorts each array into ascending order before proceeding\n");
    printf("**************************************\n\n");
    int choice = getval();
    input1 = fopen(inp_fn1, "r");

    if (choice == 1) {
        input2 = fopen(inp_fn2, "r");
    }

    if (choice == 2) {
        input2 = fopen(inp_fn3, "r");
    }

    if ((input1 != (FILE*) NULL) && (input2 != (FILE*) NULL)) {
        while (!feof(input1)) {
            fscanf(input1, "%f \n", &blanks[i]);
            nblank++;
            i++;
        }
        fclose(input1);
        while (!feof(input2)) {
            fscanf(input2, "%f \n", &qsos[j]);
            nqso++;
            j++;
        }
        fclose(input2);
        sort(nblank, blanks);
        sort(nqso, qsos);
        printf("BLANK Fields:\n");
        for (i = 0; i < nblank; i++) {
            printf("%f \n", blanks[i]);
        }
        printf("\nQSO/GRB Fields:\n");
        for (j = 0; j < nqso; j++) {
            printf("%f \n", qsos[j]);
        }

        kstwo(blanks, nblank, qsos, nqso, x, y);

        printf("\nBlank fields imported: %lu\n", nblank);
        printf("QSO/GRB fields imported: %lu\n", nqso);
        printf("K-S statistic D=%f, K-S Probability=%f\n", *x, *y);
        printf("\n**************************************\n");
        printf("Program Complete\n");
        printf("**************************************\n");

        return (EXIT_SUCCESS);
    }
    printf("\n\n**************************************\n");
    printf("Program failed and closed - check input files\n");
    printf("**************************************\n\n");
    return (EXIT_FAILURE); /*Failure due to input files being incorrect*/
}

void kstwo(float data1[], unsigned long n1, float data2[], unsigned long n2, float *d, float *prob)
{

    unsigned long j1 = 1, j2 = 1;
    float d1, d2, dt, en1, en2, en, fn1 = 0.0, fn2 = 0.0;
    en1 = n1;
    en2 = n2;
    *d = 0.0;

    while (j1 <= n1 && j2 <= n2) {
        if ((d1 = data1[j1]) <= (d2 = data2[j2])) fn1 = j1++ / en1;
        if (d2 <= d1) fn2 = j2++ / en2;
        if ((dt = fabs(fn2 - fn1)) > *d) *d = dt;
    }
    en = sqrt(en1 * en2 / (en1 + en2));
    *prob = probks((en + 0.12 + 0.11 / en)*(*d));
}

/*The Kolmogorov-Smirnov Probability Function - NRIC pg 626*/
float probks(float alam)
{

    int j;
    float a2, fac = 2.0, sum = 0.0, term, termbf = 0.0;
    a2 = -2.0 * alam*alam;
    for (j = 1; j <= 100; j++) {
        term = fac * exp(a2 * j * j);
        sum += term;
        if (fabs(term) <= EPS1 * termbf || fabs(term) <= EPS2 * sum) return sum;
        fac = -fac;
        termbf = fabs(term);
    }
    return 1.0;
}

/*Sorts a specified array into increasing order*/
void sort(unsigned long n, float arr[])
{

    int i, j;
    float a;
    for (j = 1; j < n; j++) {
        a = arr[j];
        i = j;
        while ((i > 0) && (arr[i - 1] > a)) {
            arr[i] = arr[i - 1];
            i--;
        }
        arr[i] = a;
    }
}

int getval(void)
{

    int maxmval, mt;
    char buf1[BUFSIZ];
    char *p1;
    int input_completem, true, false;
    maxmval = 2;
    true = 1;
    false = 0;
    input_completem = false;

    while (input_completem == false) {
        printf("Please enter 1 to use grb.txt or 2 to use qso.txt:");

        if (fgets(buf1, sizeof (buf1), stdin) != NULL) {
            mt = strtol(buf1, &p1, 10); /* Use of strtol assisted by http://www.mkssoftware.com/docs/man3/strtol.3.asp */

            /*Checks for initial newline character, non numerical characters and appropriate total value */
            if (buf1[0] != '\n' && (*p1 == '\n' || *p1 == '\0') && (mt > 0 && mt <= maxmval)) {
                printf("Valid number of %d entered\n\n", mt);
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