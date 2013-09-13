/****************************************
 ** 270-2A
 ** Tom Beadman
 ** 24 Nov 2010
 ** 
 ** 3X3 Determinant calculator:
 ** This programs calculates the determinant of a 3x3 matrix defined
 ** in the input file "matrix.dat"
 ****************************************/

#include <stdio.h> 
#include <stdlib.h>

float DTR(float A[3][3]);

int main(int argc, char* argv[])
{
    FILE *input;
    int i;
    float A[3][3], det;
    const char *inp_fn = "src/270/matrix.dat";

    /* Open file */
    input = fopen(inp_fn, "r");

    /* Check that file exists */
    if (input != (FILE*) NULL) {
        for (i = 0; i < 3; i++) {
            /*creates the 3x3 matrix*/
            fscanf(input, "%f %f %f\n", &A[i][0], &A[i][1], &A[i][2]); 
        }
        /* Close file */
        fclose(input);
        det = DTR(A); /* Assigns a value for the determinant */
        printf(" |%.4f  %.4f  %.4f| \n |%.4f  %.4f  %.4f| = %f \n |%.4f  %.4f  %.4f| \n",
               A[0][0], A[0][1], A[0][2], A[1][0], A[1][1], A[1][2], det, A[2][0], A[2][1], A[2][2]);
    }
    else {
        printf("*** Couldn't open the input file! ***\n");
        return(EXIT_FAILURE);
    }

    return (EXIT_SUCCESS);
}


float DTR(float A[3][3])
{
    /* Finds the value of the determinant using Laplace's formula */
    float det =
            (A[0][0]*(A[1][1] * A[2][2] - A[1][2] * A[2][1])
            - A[0][1]*(A[1][0] * A[2][2] - A[1][2] * A[2][0])
            + A[0][2]*(A[1][0] * A[2][1] - A[1][1] * A[2][0]));

    return (det);
}

