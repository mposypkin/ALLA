#include <iostream>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <common/vec.hpp>
#include "tsoffadgrad.hpp"


//-----------------------------------------------------------------------------------

void Ini_Data(int I, TersoffParams *tsof, int *lb, int *ub, double *x); // Zadanie vhodnyh dannyh
//-----------------------------------------------------------------------------------

int main() {
    int i, I, *lb, *ub;
    double *x, *g, *z, *gz;
    TersoffParams tsof; // struktura, soderzhashchaya parametry potenciala Tersoffa
    double F; // ehnergiya sistemy atomov

    I = 4; //    kolichestvo sloev;
    const int n = 3 * I;

    if ((lb = (int *) malloc(I * sizeof (int))) == NULL) {
        printf("\n  Net pamyati \n");
        exit(10);
    }
    if ((ub = (int *) malloc(I * sizeof (int))) == NULL) {
        printf("\n  Net pamyati \n");
        exit(10);
    }
    if ((x = (double *) malloc(n * sizeof (double))) == NULL) {
        printf("\n  Net pamyati \n");
        exit(10);
    }
    if ((z = (double *) malloc(n * sizeof (double))) == NULL) {
        printf("\n  Net pamyati \n");
        exit(10);
    }
    if ((g = (double *) malloc(n * sizeof (double))) == NULL) {
        printf("\n  Net pamyati \n");
        exit(10);
    }
    if ((gz = (double *) malloc(n * sizeof (double))) == NULL) {
        printf("\n  Net pamyati \n");
        exit(10);
    }
    Ini_Data(I, &tsof, lb, ub, x);

    //-----------------------------------------------------------------------
    //  Vychislenie energhii i gradienta energhii sistemy po parametram koordinat atomov
    F = grad(I, tsof, lb, ub, x, g);

    printf("\n  F = %18.12e", F);
    for (i = 0; i < n; i++) {
        printf("\n i=%d  g = %e", i, g[i]);
    }

    double alpha = 0.01;
    int cnt = 0;
    
    while (true) {
        cnt ++;
        g[1] = 0;
        snowgoose::VecUtils::vecSaxpy(n, x, g, -alpha, z);
        double Fn = grad(I, tsof, lb, ub, z, gz);
        
        if (Fn < F) {
            F = Fn;
            printf("\n %d, Fn = %lf, alpha = %lf, gnorm = %lf\n", cnt,  Fn, alpha, snowgoose::VecUtils::vecNormTwoSqr(n, gz));
            snowgoose::VecUtils::vecCopy(n, z, x);
            snowgoose::VecUtils::vecCopy(n, gz, g);
            alpha *= 2;
        } else if (alpha > 0.00001) {
            alpha *= 0.5;
        } else {
            break;
        }
//        printf("alpha = %lf\n", alpha);

        /*        char c;
                std::cout << "Continue?\n";
                std::cin >> c;
         */

    }


    //-----------------------------------------------------------------------
    if (lb != NULL) free(lb);
    if (ub != NULL) free(ub);
    if (x != NULL) free(x);
    if (z != NULL) free(z);
    if (g != NULL) free(g);
}



/**************************************************************************************************************/
/**************************************************************************************************************/
/**************************************************************************************************************/
// Zadanie vhodnyh dannyh
//  N - kolichestvo izmenyaemyh parametrov potenciala Tersoffa,
//  I - kolichestvo sloev atomov v sisteme,
//  U - parametry potenciala Tersoffa.
//  d,s,h - parametry koordinat

void Ini_Data(int I, TersoffParams *tsof, int *lb, int *ub, double *x) {
    //--------------------------------------------------------------------
    //   lb[i] --- nomer pervogo atoma na sloe i;   ub[i] --- nomer poslednego atoma na sloe i;
    lb[0] = 0;
    ub[0] = 7;
    lb[1] = 0;
    ub[1] = 7;
    lb[2] = 0;
    ub[2] = 6;
    lb[3] = 0;
    ub[3] = 6;
    //----------------------------------------------------------------------
    double z[] = {0.7125554220,0.0000000000,2.2857142816,1.4536654278,0.0000000567,2.2857142775,0.7125554121,1.1217308149,2.2927564072,1.4533956522,1.1217307869,2.2927564059};

    for (int i = 0; i < 12; i++) x[i] = z[i];

    //-----------------------------------------------------------------------------------------------
    // Zapolnenie vsekh parametrov potenciala Tersoffa
    tsof->mDe = 5.166159;
    tsof->mS = 1.576879;
    tsof->mB = 1.964037;
    tsof->mRe = 1.4471159;
    tsof->mR = 1.95;
    tsof->mRcut = 0.15;
    tsof->mC = 38049;
    tsof->mD = 4.3484;
    tsof->mH = -0.57058;
    tsof->mNu = 0.72751;
    tsof->mGamma = 1.5724E-7;
    tsof->mLambda = 0;
}

void Ini_Data_Old(int I, TersoffParams *tsof, int *lb, int *ub, double *x) {
    //--------------------------------------------------------------------
    //   lb[i] --- nomer pervogo atoma na sloe i;   ub[i] --- nomer poslednego atoma na sloe i;
    lb[0] = 0;
    ub[0] = 6;
    lb[1] = 0;
    ub[1] = 5;
    lb[2] = 0;
    ub[2] = 5;
    lb[3] = 0;
    ub[3] = 6;
    //----------------------------------------------------------------------
    double z[] = {1.46, 0, 2.5, 0.73, 1.25, 2.5, 1.46, 1.25, 2.5, 0.73, 0, 2.5};

    for (int i = 0; i < 12; i++) x[i] = z[i];

    //-----------------------------------------------------------------------------------------------
    // Zapolnenie vsekh parametrov potenciala Tersoffa
    tsof->mDe = 5.166159;
    tsof->mS = 1.576879;
    tsof->mB = 1.964037;
    tsof->mRe = 1.4471159;
    tsof->mR = 1.95;
    tsof->mRcut = 0.15;
    tsof->mC = 38049;
    tsof->mD = 4.3484;
    tsof->mH = -0.57058;
    tsof->mNu = 0.72751;
    tsof->mGamma = 1.5724E-7;
    tsof->mLambda = 0;
}
/***********************************************************************************************/
