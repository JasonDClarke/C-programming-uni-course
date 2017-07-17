#include <stdio.h> /* error 1. corrected typo from b to h */
#include <stdlib.h>
#include <math.h>

/* Define useful quantities */
#define  ERROR  1E-6
#define  SOLUTION_INTERSECT     1
#define  SOLUTION_PARALLEL      2
#define  SOLUTION_COINCIDENT    3

/* Define first line segment, by points */
#define  A_X  -2.0 /*error 2. removed superfluous semi-colon */
#define  A_Y   3.0
#define  B_X   2.0
#define  B_Y  -1.0

int main(int argc, char *argv[]) /* error 3. changed type double to int */
{
    int    validInput, solution_type;
    /* Declare intersection point */
    double I[2];
    /* Declare all vectors, for clarity */
    double a[2], b[2], c[2], d[2], b_minus_a[2], d_minus_c[2], c_minus_a[2]; 
    double Denom, Numer1, Numer2, Line1_Par, Line2_Par;  

    /***********************
     *  Input / Validation
     ***********************/

    /* Check numbers of arguments, and read input */
    validInput = (argc == 5); /* error 4. changed = to == */
    validInput = validInput && sscanf(argv[1], "%lf", &c[0]);
    validInput = validInput && sscanf(argv[2], "%lf", &c[1]);
    validInput = validInput && sscanf(argv[3], "%lf", &d[0]);
    validInput = validInput && sscanf(argv[4], "%lf", &d[1]); 
	/* error 5. above changed %f to %lf */
    
    /* For clarity assign fixed segment ends to an array */
    a[0] = A_X; /*error 6 . array a indexes changed to start at zero */
    a[1] = A_Y; /*(error 6 continued to this line)*/
    b[0] = B_X;
    b[1] = B_Y;
    
    /* Check for intersection if input validated */
    if(validInput)
    {
       /* Vector needed for a parametric description of line 1 */ 
       b_minus_a[0] = b[0] - a[0]; b_minus_a[1] = b[1] - a[1];

       /* Vector needed for a parametric description of line 2 */ 
       d_minus_c[0] = d[0] - c[0]; d_minus_c[1] = d[1] - c[1];
       
       /* Other useful vector differences */
       c_minus_a[0] = c[0] - a[0]; c_minus_a[1] = c[1] - a[1];
       
       /* Denominator, common to both equations of line parameters */
       Denom = -b_minus_a[1] * d_minus_c[0] + b_minus_a[0] * d_minus_c[1];

       /* Numerators 1 and 2 */
       Numer1 =  -d_minus_c[0] * c_minus_a[1] + c_minus_a[0] * d_minus_c[1];
       Numer2 =   c_minus_a[0] * b_minus_a[1] - b_minus_a[0] * c_minus_a[1];
       
       if( (fabs(Denom) <= ERROR)&&(fabs(Numer1) <= ERROR) ) /*error 7. fabs(Denom <= ERROR) written twice */
          solution_type = SOLUTION_COINCIDENT;
	/***** NEW CODE: solution for parallel case ***/
       else if ( (fabs(Denom) <= ERROR)&&(fabs(Numer1) >= ERROR) )
          solution_type = SOLUTION_PARALLEL;
	/*** NEW CODE: solution for parallel case ***/
       else
       {
          /* Calculate line paremters */
          Line1_Par = Numer1 / Denom;
          Line2_Par = Numer2 / Denom;
          
          /* Point of intersection */
          if( (Line1_Par >= 0)&&(Line1_Par <= 1)&&(Line2_Par >= 0)&&(Line2_Par <= 1) )
          {
               I[0] = a[0] + Line1_Par * b_minus_a[0];
               I[1] = a[1] + Line1_Par * b_minus_a[1];
			   solution_type = SOLUTION_INTERSECT;
		  }
       }
              /* Print appropriate message */
       if(solution_type==SOLUTION_COINCIDENT)
          printf("Line segments are coincident\n");
       else if(solution_type==SOLUTION_PARALLEL)
          printf("Line segments are parallel\n");
       else if(solution_type==SOLUTION_INTERSECT)
          printf("Line segments intersect at (%d, %d)\n", I[0], I[1]); /*error 8. changed %lf to %d */
       else
          printf("No intersection for these segments\n");

    }
    else
    {
       printf("*****************************************\n");
       printf("* !!!  User input error detected !!!    *\n");
       printf("*   Execution method: prog2A a b c d    *\n");
       printf("* where a, b, c and d are real numbers  *\n");
       printf("*****************************************\n");
    }
    
    return(EXIT_SUCCESS);
}