/****************************************
** 270-1A
** Tom Beadman
** 10 Nov 2010
** 
** Quadratic equation solver:
** This program determines the roots to a quadratic equation, 
** the coefficients of which are passed in as command line arguments.  
****************************************/


#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{

    float a, b, c, d;
    float q, q1, q2, root1, root2;
    int validInput;

    validInput = (argc == 4);
    validInput = validInput && sscanf(argv[1], "%f", &a);
    validInput = validInput && sscanf(argv[2], "%f", &b);
    validInput = validInput && sscanf(argv[3], "%f", &c);

    if (!validInput) {
        printf("Please ensure that valid numeric coefficients are entered as command line arguments \n");
        return ( EXIT_FAILURE);
    }

    else {

        printf("Calculating roots for %.3f X^2 + %.3f X + %.3f: \n", a, b, c);
        if (a == 0) { /*case where a equal to 0*/
            if (b == 0) {
                printf("Error: no X or X^2 terms entered \n");
                return EXIT_FAILURE;
            }
            else {
                q = -b;
                root2 = c / q;
                printf("Single root due to no X^2 term: %.3f \n", root2);
                return EXIT_SUCCESS;
            }
        }

        else { /*case when a does not equal to 0*/

            d = b * b - 4 * a*c;
            q1 = sqrt(d);
            q2 = (b < 0) ? b - q1 : b + q1;
            q = -(q2) / 2;

            root1 = q / a;
            root2 = c / q;
        }

        if (d > 0) { /*the discriminant is greater than 0*/
            printf("There are two real roots: \n");
            printf("Root1: %.6f \n", root1);
            printf("Root2: %.6f \n", root2);
        }

        else if (d == 0) { /*the discriminant is equal to 0*/
            printf("There is one repeated root: \n");
            printf("Root: %.6f \n", root1);
        }

        else { /*the discriminant is less than 0*/
            printf("There are no real roots \n");
        }
    }
}

