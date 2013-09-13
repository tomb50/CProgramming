/****************************************
 ** 270-4A
 ** Tom Beadman
 ** 27 Jan 2011
 ** 
 ** Restricted Three Body Problem:
 ** This program models the 3 body problem of interacting
 ** masses, the restriction is that the mass of the 3rd body is 
 ** small enough to not effect the remaining two.
 ** Initial masses, positions and velocities as passed as command line arguments
 ** as detailed in the prompt to the user. The bodies paths are printed to
 ** output file 'prog4A.out'.
 ****************************************/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

float Ang(float, float, float, float, float *, float *); /*External function*/

int main(int argc, char *argv[])
{

    FILE* output; /*Variable declarations*/
    float mneg, mpos, tempx, tempy, ANGx, ANGy;
    float x0, vx0, y0, vy0, delt/*mn, mp*/;
    float r[12000][2];
    float t[12000];
    int i, validInput;

    validInput = (argc == 7); /*Input validation*/
    validInput = validInput && sscanf(argv[1], "%f", &mneg);
    validInput = validInput && sscanf(argv[2], "%f", &mpos);
    validInput = validInput && sscanf(argv[3], "%f", &x0);
    validInput = validInput && sscanf(argv[4], "%f", &y0);
    validInput = validInput && sscanf(argv[5], "%f", &vx0);
    validInput = validInput && sscanf(argv[6], "%f", &vy0);
    validInput = validInput && mneg >= 0;
    validInput = validInput && mpos >= 0;

    if (!validInput) { /*Input error message*/
        printf("***Error in variable input.***\n");
        printf("***Enter program name followed by variable values: M-, M+, X0, Y0, Vx0, Vy0***\n");
        printf("***Values for M- and M+ must be non negative***\n");
        printf("***Program terminated***\n");
        exit(1);
    }

    else { /*Main program body*/

        output = fopen("prog4A.out", "w");
        delt = 0.005;

        for (i = 0; i < 1; i++) { /*Initial values*/
            r[i][0] = x0;
            r[i][1] = y0;
            t[i] = (i)*(delt);
            fprintf(output, "%f %f %f\n", t[i], r[i][0], r[i][1]);
            if ((r[i][1] == 0) && abs(r[i][0]) == 1) {
                fprintf(output, "********IMPACT********\n");
                exit(0);
            }
        }
        for (i == 1; i < 2; i++) { /*Using initial velocity parameters*/
            r[i][0] = x0 + (vx0) * delt;
            r[i][1] = y0 + (vy0) * delt;
            t[i] = (i)*(delt);
            fprintf(output, "%f %f %f\n", t[i], r[i][0], r[i][1]);
            if ((r[i][1] == 0) && abs(r[i][0]) == 1) {
                fprintf(output, "********IMPACT********\n");
                exit(0);
            }
        }
        for (i == 2; i < 12000; i++) { /*Main calculation loop*/
            tempx = r[i - 1][0];
            tempy = r[i - 1][1];
            if (Ang(tempx, tempy, mneg, mpos, &ANGx, &ANGy) == 1) {
                return (1);
            }

            r[i][0] = (delt)*(delt)*(ANGx)+(2)*(r[i - 1][0])-(r[i - 2][0]);
            r[i][1] = (delt)*(delt)*(ANGy)+(2)*(r[i - 1][1])-(r[i - 2][1]);
            t[i] = (i)*(delt);
            fprintf(output, "%f %f %f\n", t[i], r[i][0], r[i][1]);
            if ((r[i][1] == 0) && abs(r[i][0]) == 1) {
                fprintf(output, "********IMPACT********\n");
                exit(0);
            }
        }
    }

    return (0);
}

float Ang(float atempx, float atempy, float mn, float mp, float *pANGx, float *pANGy)/*External function*/
{
    // /*Variable declaration*/
    float cube1, cube2, squarexpos, squarexneg, squarey;
    float ang1x, ang1y, ang2x, ang2y;

    squarexpos = (atempx + 1)*(atempx + 1); /*Value assignments and calculation*/
    squarexneg = (atempx - 1)*(atempx - 1);
    squarey = (atempy)*(atempy);
    cube1 = pow(squarexpos + squarey, (float) 3 / 2);
    cube2 = pow(squarexneg + squarey, (float) 3 / 2);
    ang1x = (-1)*(mn * (atempx + 1)) / cube1;
    ang1y = (-1)*(mn * atempy) / cube1;
    ang2x = (-1)*(mp * (atempx - 1)) / cube2;
    ang2y = (-1)*(mp * atempy) / cube2;
    *pANGx = ang1x + ang2x;
    *pANGy = ang1y + ang2y;

    return (0);
}
