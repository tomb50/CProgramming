/****************************************
 ** 270-3A
 ** Tom Beadman
 ** 8 Dec 2010
 ** 
 ** Numeric Integration Calculator:
 ** This program determines the the integral of the cosine function provided 
 ** in "prog3a.dat" from 0 up a value provided as a command line argument.
 ** The value must be less than Pi/2.   
 ****************************************/

#include <stdio.h> 
#include <stdlib.h>
#include <math.h>

#define PI (3.141592653589793)

float INTEGRAL(float dx[], float fx[], int end, float B, float A);

int main(int argc, char *argv[])
{
    FILE *input;
    int i, end, validInput;
    float tempx[10000], tempfx[10000], B, A, *x, *fx, *dx, gral;
    const char inp_fn[] = "src/270/prog3a.dat";

    validInput = (argc == 2);
    /*validate upper limit type and value*/
    validInput = validInput && sscanf(argv[1], "%f", &B); 
    validInput = validInput && (B >= 0) && (B <= PI / 2);

    if (!validInput) {
        printf("Error: B value must be between 0 and Pi/2 \n");

    }
    else {
        /* Opens file */
        input = fopen(inp_fn, "r");

        /* Check that file exists */
        if (input != (FILE*) NULL) {
            for (i = 0; tempx[i - 1] <= B; i++) {
                /* read file and makes two arrays slightly larger than the upper limit*/
                fscanf(input, "%f %f\n", &tempx[i], &tempfx[i]); 
                end = i;
            }
            fclose(input);
            /*allocate appropriate memory for x, fx and dx arrays*/
            x = malloc((end - 1) * sizeof (float));
            fx = malloc((end - 2) * sizeof (float)); 
            dx = malloc((end - 2) * sizeof (float));
            i = 0;
            for (i = 0; i <= (end - 1); i++) {
                x[i] = tempx[i];
                fx[i] = tempfx[i]; /*creates x and fx arrays */
            }
            i = 0;
            for (i = 0; i <= (end - 2); i++) {
                dx[i] = (x[i + 1] - x[i]); /*creates dx*/
            }

            gral = INTEGRAL(dx, fx, end, B, A);
            if (B == 0)
                printf("integral of f(x) between 0 and 0 is 0\n");
            else /*prints integral not needed but nice too know its doing something*/
                printf("integral of f(x) between 0 and %f is %f\n", B, gral);
        }
        else {
            printf("*** Couldn't open the input file! ***\n");
        }
        free(x);
        free(fx);
        free(dx);
        return (0);
    }
}

float INTEGRAL(float dx[], float fx[], int end, float B, float A) /*external function*/
{
    int i;
    float gral = 0;

    for (i = 0; i <= (end - 2); i++) {
        gral = (gral + (fx[i] * dx[i])); /*calculates intergral*/
    }
    return (gral);
}
