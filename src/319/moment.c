/****************************************
 ** 319-1A
 ** Tom Beadman
 ** 20 Feb 2012
 ** 
 ** Statistical Moments:
 ** This program calculates the statistical moments from a 100 element 
 ** distribution contained in blank.txt.
 ** 
 ** Currently the program is producing nonsense results half of the time and
 ** correct results the other half, this needs to be fixed.  
 ****************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>


void moment(float data[], int n, float *ave, float *adev, float *sdev, float *var, float *skew, float *curt);
void sort(unsigned long n, float arr[]);


int main()
{

    unsigned int i = 0;
    unsigned int j = 0;
    FILE *input;
    const char inp_fn[] = "src/319/blank.txt";
    float* ave = malloc(sizeof (float));
    float* adev = malloc(sizeof (float));
    float* sdev = malloc(sizeof (float));
    float* var = malloc(sizeof (float));
    float* skew = malloc(sizeof (float));
    float* curt = malloc(sizeof (float));
    float data[100];

    printf("\n\n**************************************\n");
    printf("Program to determine statistical moments for a given dataset\n");
    printf("**************************************\n\n");


    input = fopen(inp_fn, "r");
    if ((input != (FILE*) NULL)) {
        while (!feof(input)) {
            fscanf(input, "%f \n", &data[i]);
            i++;
        }

        printf("\nTransferred:\n");
        for (j = 0; j < i; j++) {
            printf("%f\n", data[j]);
        }

        moment(data, j, ave, adev, sdev, var, skew, curt);
        printf("\nMoment Statistics\n");
        printf("********************************\n");
        printf("Number of points: %u\n", j);
        printf("Mean value: %f\n", *ave);
        printf("Average Deviation: %f\n", *adev);
        printf("Standard Deviation: %f\n", *sdev);
        printf("Variance: %f\n", *var);
        printf("Skewness: %f\n", *skew);
        printf("Kurtosis: %f\n", *curt);
        printf("********************************\n");


        //sort(j, data);
        printf("\nSorted:\n");
        for (j = 0; j < i; j++) {
            printf("%f\n", data[j]);
        }

        printf("\nMoment Statistics\n");
        printf("********************************\n");
        printf("Number of points: %u\n", j);
        printf("Mean value: %f\n", *ave);
        printf("Average Deviation: %f\n", *adev);
        printf("Standard Deviation: %f\n", *sdev);
        printf("Variance: %f\n", *var);
        printf("Skewness: %f\n", *skew);
        printf("Kurtosis: %f\n", *curt);
        printf("********************************\n");

        return (EXIT_SUCCESS);
    }


    return (EXIT_FAILURE);

}

void moment(float data[], int n, float *ave, float *adev, float *sdev, float *var, float *skew, float *curt)
{

    int j;
    float ep = 0.0, s, p;

    if (n <= 1) {
        printf("n must be at least 2 in moment");
        exit(EXIT_FAILURE);
    }
    s = 0.0;
    for (j = 1; j <= n; j++) s += data[j];
    *ave = s / n;
    *adev = (*var) = (*skew) = (*curt) = 0.0;
    for (j = 1; j <= n; j++) {
        *adev += fabs(s = data[j]-(*ave));
        ep += s;
        *var += (p = s * s);
        *skew += (p *= s);
        *curt += (p *= s);
    }
    *adev /= n;
    *var = (*var - ep * ep / n) / (n - 1);
    *sdev = sqrt(*var);
    if (*var) {
        *skew /= (n * (*var)*(*sdev));
        *curt = (*curt) / (n * (*var)*(*var)) - 3.0;
    }
    else {
        printf("No skew/kurtosis when variance = 0 (in moment)");
        exit(EXIT_FAILURE);
    }
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