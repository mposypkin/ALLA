/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   tsoffadgrad.hpp
 * Author: posypkin
 *
 * Created on November 13, 2016, 1:50 PM
 */

#ifndef TSOFFADGRAD_HPP
#define TSOFFADGRAD_HPP

#define VARIABLE_TERSOFF_PARAMS 10




//-----------------------------------------------------------------------------------
// Ob"yavlenie struktury Param_Tersoff

struct TersoffParams {
    double mDe, mRe, mB, mS, mNu, mGamma, mLambda, mC, mD, mH, mR, mRcut; // izmenyaemye parametry - pervye 10
};


/**************************************************************************************************************/
/**************************************************************************************************************/

/**************************************************************************************************************/
double expnew(double XX) {
    double xz, az, enew;
    xz = 36.0; // 40.0;
    az = exp(xz)*(0.5 * XX * XX + (1.0 - xz) * XX + (1.0 - xz + 0.5 * XX * XX));
    if ((XX <= xz)&&(XX >= (-xz))) enew = exp(XX);
    if (XX > xz) enew = az;
    if (XX < (-xz)) enew = 1.0 / az;

    return (enew);
}

int m_index(int I, int *Jz, int *Jzz) {
    int i, m;
    m = 0;
    for (i = 1; i <= I; i++) m = m + Jzz[i] - Jz[i] + 1;
    return (m);
}

class TersoffFAD {
private:
    double **r; // rasstoyaniya 
    double **rr; // rasstoyaniya 
    double **rrr; // rasstoyaniya 
    double **z; // fazovye peremennye
    double **v; // fazovye peremennye
    double **w; // fazovye peremennye
    double **R; // fazovye peremennye
    double **T; // fazovye peremennye
    double **Q; // fazovye peremennye
    double **V; // fazovye peremennye
    double **f; // Sglazhivayushchaya funkciya f
    double **df; // Proizvodnye sglazhivayushchej funkcii f
    double **pf; // Sopryazhennye peremennye
    double **pr; // Sopryazhennye peremennye
    double **pz; // Sopryazhennye peremennye
    double **pv; // Sopryazhennye peremennye
    double **pw; // Sopryazhennye peremennye
    double **pR; // Sopryazhennye peremennye
    double **pT; // Sopryazhennye peremennye
    double **pQ; // Sopryazhennye peremennye
    double **pV; // Sopryazhennye peremennye
    double **d; // Rabochij massiv


    double ***ar; // rasstoyaniya 
    double ***br; // rasstoyaniya 
    double ***q; // ������� ����
    double ***g; // fazovye peremennye
    double ***x; // fazovye peremennye
    double ***y; // fazovye peremennye
    double ***ff; // Sglazhivayushchaya funkciya f
    double ***dff; // Proizvodnye sglazhivayushchej funkcii f
    double ***pq; // Sopryazhennye peremennye
    double ***pg; // Sopryazhennye peremennye
    double ***px; // Sopryazhennye peremennye
    double ***py; // Sopryazhennye peremennye
    double ***pff; // Sopryazhennye peremennye
    double ***par; // Sopryazhennye peremennye
    double ***pbr; // Sopryazhennye peremennye
    double ***tau; // rij - rik

    int *Jz, *Jzz, *im, *jm, *g_s, *g_h;
    double *x1, *x2, *gr_x1, *gr_x2, F, *grad_h, *grad_d, *grad_s, *fad_d, *fad_s, *h, *dm, *dp;

    int mMM = 0;
    int mII = 0;

    void lazyAllocArrays(int MM) {

        if (MM != mMM) {
            if (mMM > 0) {
                deleteArrays(mMM);
            }
            allocArrays(MM);
            mMM = MM;
        }
    }

    void deleteArrays(int MM) {
        // Osvobozhdenie vspomogatel'noj pamyati dlya dvumernyh massivov

        for (int i = 0; i < MM; i++) {
            delete[] r[i];
            delete[] rr[i];
            delete[] rrr[i];
            delete[] z[i];
            delete[] v[i];
            delete[] w[i];
            delete[] R[i];
            delete[] T[i];
            delete[] Q[i];
            delete[] V[i];
            delete[] f[i];
            delete[] df[i];
            delete[] pf[i];
            delete[] pr[i];
            delete[] pz[i];
            delete[] pv[i];
            delete[] pw[i];
            delete[] pR[i];
            delete[] pT[i];
            delete[] pQ[i];
            delete[] pV[i];
            delete[] d[i];
        }
        delete[] r;
        delete[] rr;
        delete[] rrr;
        delete[] z;
        delete[] v;
        delete[] w;
        delete[] R;
        delete[] T;
        delete[] Q;
        delete[] V;
        delete[] f;
        delete[] df;
        delete[] pf;
        delete[] pr;
        delete[] pz;
        delete[] pv;
        delete[] pw;
        delete[] pR;
        delete[] pT;
        delete[] pQ;
        delete[] pV;
        delete[] d;

        //-----------------------------------------------------------------------
        // Osvobozhdenie vspomogatel'noj pamyati dlya trekhmernyh massivov

        for (int i = 0; i < MM; i++) {
            for (int j = 0; j < MM; j++) {
                delete[] ar[i][j];
                delete[] br[i][j];
                delete[] q[i][j];
                delete[] g[i][j];
                delete[] x[i][j];
                delete[] y[i][j];
                delete[] ff[i][j];
                delete[] dff[i][j];
                delete[] pq[i][j];
                delete[] pg[i][j];
                delete[] px[i][j];
                delete[] py[i][j];
                delete[] pff[i][j];
                delete[] par[i][j];
                delete[] pbr[i][j];
                delete[] tau[i][j];
            }

            delete[] ar[i];
            delete[] br[i];
            delete[] q[i];
            delete[] g[i];
            delete[] x[i];
            delete[] y[i];
            delete[] ff[i];
            delete[] dff[i];
            delete[] pq[i];
            delete[] pg[i];
            delete[] px[i];
            delete[] py[i];
            delete[] pff[i];
            delete[] par[i];
            delete[] pbr[i];
            delete[] tau[i];
        }

        delete[] ar;
        delete[] br;
        delete[] q;
        delete[] g;
        delete[] x;
        delete[] y;
        delete[] ff;
        delete[] dff;
        delete[] pq;
        delete[] pg;
        delete[] px;
        delete[] py;
        delete[] pff;
        delete[] par;
        delete[] pbr;
        delete[] tau;

        if (x1 != NULL) free(x1);
        if (x2 != NULL) free(x2);
        if (gr_x1 != NULL) free(gr_x1);
        if (gr_x2 != NULL) free(gr_x2);
        if (im != NULL) free(im);
        if (jm != NULL) free(jm);
        if (g_h != NULL) free(g_h);
        if (g_s != NULL) free(g_s);

        // Konec osvobozhdeniya vspomogatel'noj pamyati
        //-----------------------------------------------------------------------
    }

    void allocArrays(int MM) {
        //---------------------------------------------------------------------------
        // Vydelenie vspomogatel'noj pamyati dlya dvumernyh massivov

        r = new double*[MM]; // rasstoyaniya 
        rr = new double*[MM]; // rasstoyaniya 
        rrr = new double*[MM]; // rasstoyaniya 
        z = new double*[MM]; // fazovye peremennye
        v = new double*[MM]; // fazovye peremennye
        w = new double*[MM]; // fazovye peremennye
        R = new double*[MM]; // fazovye peremennye
        T = new double*[MM]; // fazovye peremennye
        Q = new double*[MM]; // fazovye peremennye
        V = new double*[MM]; // fazovye peremennye
        f = new double*[MM]; // Sglazhivayushchaya funkciya f
        df = new double*[MM]; // Proizvodnye sglazhivayushchej funkcii f
        pf = new double*[MM]; // Sopryazhennye peremennye
        pr = new double*[MM]; // Sopryazhennye peremennye
        pz = new double*[MM]; // Sopryazhennye peremennye
        pv = new double*[MM]; // Sopryazhennye peremennye
        pw = new double*[MM]; // Sopryazhennye peremennye
        pR = new double*[MM]; // Sopryazhennye peremennye
        pT = new double*[MM]; // Sopryazhennye peremennye
        pQ = new double*[MM]; // Sopryazhennye peremennye
        pV = new double*[MM]; // Sopryazhennye peremennye
        d = new double*[MM]; // Rabochij massiv

        for (int i = 0; i < MM; i++) {
            r[i] = new double[MM];
            rr[i] = new double[MM];
            rrr[i] = new double[MM];
            z[i] = new double[MM];
            v[i] = new double[MM];
            w[i] = new double[MM];
            R[i] = new double[MM];
            T[i] = new double[MM];
            Q[i] = new double[MM];
            V[i] = new double[MM];
            f[i] = new double[MM];
            df[i] = new double[MM];
            pf[i] = new double[MM];
            pr[i] = new double[MM];
            pz[i] = new double[MM];
            pv[i] = new double[MM];
            pw[i] = new double[MM];
            pR[i] = new double[MM];
            pT[i] = new double[MM];
            pQ[i] = new double[MM];
            pV[i] = new double[MM];
            d[i] = new double[MM];
        }

        //------------------------------------------------------------------------------
        // Vydelenie vspomogatel'noj pamyati dlya trekhmernyh massivov

        ar = new double**[MM]; // rasstoyaniya 
        br = new double**[MM]; // rasstoyaniya 
        q = new double**[MM]; // ������� ����
        g = new double**[MM]; // fazovye peremennye
        x = new double**[MM]; // fazovye peremennye
        y = new double**[MM]; // fazovye peremennye
        ff = new double**[MM]; // Sglazhivayushchaya funkciya f
        dff = new double**[MM]; // Proizvodnye sglazhivayushchej funkcii f
        pq = new double**[MM]; // Sopryazhennye peremennye
        pg = new double**[MM]; // Sopryazhennye peremennye
        px = new double**[MM]; // Sopryazhennye peremennye
        py = new double**[MM]; // Sopryazhennye peremennye
        pff = new double**[MM]; // Sopryazhennye peremennye
        par = new double**[MM]; // Sopryazhennye peremennye
        pbr = new double**[MM]; // Sopryazhennye peremennye
        tau = new double**[MM]; // rij - rik

        for (int i = 0; i < MM; i++) {
            ar[i] = new double*[MM];
            br[i] = new double*[MM];
            q[i] = new double*[MM];
            g[i] = new double*[MM];
            x[i] = new double*[MM];
            y[i] = new double*[MM];
            ff[i] = new double*[MM];
            dff[i] = new double*[MM];
            pq[i] = new double*[MM];
            pg[i] = new double*[MM];
            px[i] = new double*[MM];
            py[i] = new double*[MM];
            pff[i] = new double*[MM];
            par[i] = new double*[MM];
            pbr[i] = new double*[MM];
            tau[i] = new double*[MM];

            for (int j = 0; j < MM; j++) {
                ar[i][j] = new double[MM];
                br[i][j] = new double[MM];
                q[i][j] = new double[MM];
                g[i][j] = new double[MM];
                x[i][j] = new double[MM];
                y[i][j] = new double[MM];
                ff[i][j] = new double[MM];
                dff[i][j] = new double[MM];
                pq[i][j] = new double[MM];
                pg[i][j] = new double[MM];
                px[i][j] = new double[MM];
                py[i][j] = new double[MM];
                pff[i][j] = new double[MM];
                par[i][j] = new double[MM];
                pbr[i][j] = new double[MM];
                tau[i][j] = new double[MM];
            }
        }

        if ((x1 = (double *) malloc(MM * sizeof (double))) == NULL) {
            printf("\n  Net pamyati \n");
            exit(10);
        }
        if ((x2 = (double *) malloc(MM * sizeof (double))) == NULL) {
            printf("\n  Net pamyati \n");
            exit(10);
        }
        if ((gr_x1 = (double *) malloc(MM * sizeof (double))) == NULL) {
            printf("\n  Net pamyati \n");
            exit(10);
        }
        if ((gr_x2 = (double *) malloc(MM * sizeof (double))) == NULL) {
            printf("\n  Net pamyati \n");
            exit(10);
        }
        if ((im = (int *) malloc(MM * sizeof (int))) == NULL) {
            printf("\n  Net pamyati \n");
            exit(10);
        }
        if ((jm = (int *) malloc(MM * sizeof (int))) == NULL) {
            printf("\n  Net pamyati \n");
            exit(10);
        }
        if ((g_s = (int *) malloc(MM * sizeof (int))) == NULL) {
            printf("\n  Net pamyati \n");
            exit(10);
        }
        if ((g_h = (int *) malloc(MM * sizeof (int))) == NULL) {
            printf("\n  Net pamyati \n");
            exit(10);
        }
    }

    double GRAD(int MP, int M, int N, TersoffParams tsof, double *x1, double *x2, double *gr_x1, double *gr_x2) {
        int i, j, k, I, IM;
        double a, b, c, hlp, pE, E;
        int NN, MM;
        double RR, RC; // zadannye parametry 

        double r1, r2, r3, r4, r5, r6, r7, r8, r9, r10, r12, r14, r15;
        double PI = 3.141592653589793238462643;

        double u[VARIABLE_TERSOFF_PARAMS]; // vektor izmenyayushchihsya parametrov potenciala Tersoffa
        //double *u;

        MM = MP + 1;
        I = MP;
        IM = M;
        NN = N + 1;
        //u = new double[NN];


        //------------------------------------------------------------------------------
        //  Obnulenie vspomogatel'nyh massivov

        for (i = 0; i <= I; i++)
            for (j = 0; j <= I; j++) {
                r[i][j] = 0.0;
                rr[i][j] = 0.0;
                rrr[i][j] = 0.0;
                z[i][j] = 0.0;
                v[i][j] = 0.0;
                w[i][j] = 0.0;
                R[i][j] = 0.0;
                T[i][j] = 0.0;
                Q[i][j] = 0.0;
                V[i][j] = 0.0;
                pz[i][j] = 0.0;
                pv[i][j] = 0.0;
                pw[i][j] = 0.0;
                pR[i][j] = 0.0;
                pT[i][j] = 0.0;
                pQ[i][j] = 0.0;
                pV[i][j] = 0.0;
                pf[i][j] = 0.0;
                pr[i][j] = 0.0;
                df[i][j] = 0.0;
                f[i][j] = 0.0;
                d[i][j] = 0.0;
                pr[i][j] = 0.0;
                for (k = 0; k <= I; k++) {
                    q[i][j][k] = 0.0;
                    g[i][j][k] = 0.0;
                    x[i][j][k] = 0.0;
                    y[i][j][k] = 0.0;
                    tau[i][j][k] = 0.0;
                    ar[i][j][k] = 0.0;
                    br[i][j][k] = 0.0;
                    ff[i][j][k] = 0.0;
                    pg[i][j][k] = 0.0;
                    px[i][j][k] = 0.0;
                    py[i][j][k] = 0.0;
                    pq[i][j][k] = 0.0;
                    par[i][j][k] = 0.0;
                    pbr[i][j][k] = 0.0;
                    pff[i][j][k] = 0.0;
                    dff[i][j][k] = 0.0;
                }
            }
        //------------------------------------------------------------------------------

        //------------------------------------------------------------------------------
        //  Zapolnenie massiva parametrov
        u[1] = tsof.mDe;
        u[2] = tsof.mRe;
        u[3] = tsof.mB;
        u[4] = tsof.mS;
        u[5] = tsof.mNu;
        u[6] = tsof.mGamma;
        u[7] = tsof.mLambda;
        u[8] = tsof.mC;
        u[9] = tsof.mD;
        u[10] = tsof.mH;
        RR = tsof.mR;
        RC = tsof.mRcut;
        //------------------------------------------------------------------------------
        r1 = u[8] * u[8];
        r2 = u[9] * u[9];
        r3 = 1.0 + r1 / r2;
        r4 = u[7] * u[7] * u[7];
        r5 = -u[3] * sqrt(2.0 * u[4]);
        r6 = u[1] / (u[4] - 1);
        r7 = u[1] * u[4] / (u[4] - 1);
        r8 = -u[3] * sqrt(2.0 / u[4]);
        r9 = 1.0 / u[4] - 1.0;
        r15 = -1.0 / 2.0 / u[5];

        E = 0.0;
        for (i = 1; i <= IM; i++)
            for (j = 1; j <= I; j++) {
                if (i != j) {
                    r[i][j] = sqrt((x1[i] - x1[j])*(x1[i] - x1[j]) + (x2[i] - x2[j])*(x2[i] - x2[j]));

                    if ((r[i][j]) <= (RR + RC)) {

                        c = r[i][j];
                        if ((c) < (RR - RC)) f[i][j] = 1.0;
                        if (((c) >= (RR - RC)) && ((c) <= (RR + RC))) f[i][j] = 0.5 * (1.0 - sin(PI * (c - RR) / (2.0 * RC)));
                        z[i][j] = 0.0;
                        //-----------------------------------------------------------------------------------------------
                        for (k = 1; k <= I; k++)
                            if ((k != i) && (k != j)) {
                                br[i][j][k] = sqrt((x1[i] - x1[k])*(x1[i] - x1[k]) + (x2[i] - x2[k])*(x2[i] - x2[k]));
                                ar[i][j][k] = sqrt((x1[j] - x1[k])*(x1[j] - x1[k]) + (x2[j] - x2[k])*(x2[j] - x2[k]));
                                if ((br[i][j][k]) <= (RR + RC)) {

                                    a = ar[i][j][k];
                                    b = br[i][j][k];
                                    q[i][j][k] = (b * b + c * c - a * a) / (2.0 * b * c);
                                    tau[i][j][k] = (c - b) * (c - b) * (c - b);

                                    g[i][j][k] = r3 - r1 / (r2 + (u[10] - q[i][j][k])*(u[10] - q[i][j][k]));
                                    x[i][j][k] = expnew(tau[i][j][k] * r4);

                                    if ((b) < (RR - RC)) ff[i][j][k] = 1.0;
                                    if (((b) >= (RR - RC)) && ((b) <= (RR + RC))) ff[i][j][k] = 0.5 * (1.0 - sin(PI * (b - RR) / (2.0 * RC)));
                                    y[i][j][k] = ff[i][j][k] * x[i][j][k] * g[i][j][k];

                                    z[i][j] = z[i][j] + y[i][j][k];
                                }
                            } //         k
                        //-------------------------------------------------------------------------------------------------
                        d[i][j] = expnew(r5 * (c - u[2]));
                        v[i][j] = u[6] * z[i][j];
                        w[i][j] = pow(v[i][j], u[5]);
                        R[i][j] = pow((1.0 + w[i][j]), r15);
                        T[i][j] = r6 * d[i][j];
                        Q[i][j] = r7 * expnew(r8 * (c - u[2]));
                        V[i][j] = f[i][j] * r6 * d[i][j]*(1.0 - R[i][j] * u[4] * pow(d[i][j], (r9)));
                        E = E + V[i][j];
                    }
                }
            } //         i, j

        r10 = -1.0 - 0.5 / u[5];
        pE = 1.0;
        for (i = 1; i <= IM; i++)
            for (j = 1; j <= I; j++) {
                if (i != j) {

                    if ((r[i][j]) <= (RR + RC)) {
                        pV[i][j] = 1.0;
                        pQ[i][j] = -f[i][j] * R[i][j];
                        pT[i][j] = f[i][j];
                        pR[i][j] = -f[i][j] * Q[i][j];
                        pw[i][j] = 0.5 * f[i][j] * Q[i][j] * pow((1.0 + w[i][j]), (r10)) / u[5];
                        if (v[i][j] < 1.0e-25) hlp = pow((1.0e-25), (u[5] - 1));
                        else hlp = pow((v[i][j]), (u[5] - 1));
                        pv[i][j] = u[5] * hlp * pw[i][j];
                        pz[i][j] = u[6] * pv[i][j];
                        pf[i][j] = r6 * d[i][j]*(1 - R[i][j] * u[4] * pow(d[i][j], r9));
                        c = r[i][j];
                        if ((c) < (RR - RC)) df[i][j] = 0.0;
                        if (((c) >= (RR - RC)) && ((c) <= (RR + RC)))
                            df[i][j] = -PI * cos(PI * (c - RR) / (2.0 * RC)) / (4.0 * RC);

                        pr[i][j] = Q[i][j] * r8 * pQ[i][j] + T[i][j] * r5 * pT[i][j] + df[i][j] * pf[i][j];
                        //--------------------------------------------------------------------------------
                        for (k = 1; k <= I; k++)
                            if ((k != i) && (k != j)) {
                                if ((br[i][j][k]) <= (RR + RC)) {

                                    py[i][j][k] = pz[i][j];
                                    px[i][j][k] = ff[i][j][k] * g[i][j][k] * py[i][j][k];
                                    pg[i][j][k] = ff[i][j][k] * x[i][j][k] * py[i][j][k];
                                    r12 = u[10] - q[i][j][k];
                                    r14 = (r2 + (r12)*(r12))*(r2 + (r12)*(r12));
                                    pq[i][j][k] = -2.0 * r1 * (r12) * pg[i][j][k] / (r14);
                                    a = ar[i][j][k];
                                    b = br[i][j][k];
                                    if ((b) < (RR - RC)) dff[i][j][k] = 0.0;
                                    if (((b) >= (RR - RC)) && ((b) <= (RR + RC))) dff[i][j][k] = -PI * cos(PI * (b - RR) / (2.0 * RC)) / (4.0 * RC);

                                    pff[i][j][k] = x[i][j][k] * g[i][j][k] * py[i][j][k];

                                    pr[i][j] = pr[i][j] + (3.0 * x[i][j][k] * r4 * (c - b)*(c - b)) * px[i][j][k] +
                                            ((4.0 * b * c * c - 2.0 * b * (b * b + c * c - a * a)) / (4.0 * b * b * c * c)) * pq[i][j][k];

                                    par[i][j][k] = par[i][j][k] - ((4.0 * a * b * c) / (4.0 * b * b * c * c)) * pq[i][j][k];

                                    pbr[i][j][k] = pbr[i][j][k] - (3.0 * x[i][j][k] * r4 * (r[i][j] - b)*(c - b)) * px[i][j][k] +
                                            ((4.0 * b * b * c - 2.0 * c * (b * b + c * c - a * a)) /
                                            (4.0 * b * b * c * c)) * pq[i][j][k] + dff[i][j][k] * pff[i][j][k];
                                }
                            }
                        //-----------------------------------------------------------------------------------------------
                    }
                }
            }
        for (i = 0; i <= I; i++) {
            gr_x1[i] = 0.0;
            gr_x2[i] = 0.0;
        }

        for (i = 1; i <= IM; i++)
            for (j = 1; j <= I; j++) {
                if (i != j) {

                    if ((r[i][j]) <= (RR + RC)) {
                        //				if (i<=IM)
                        {
                            gr_x1[i] = gr_x1[i] + (x1[i] - x1[j]) * pr[i][j] / r[i][j];
                            gr_x2[i] = gr_x2[i] + (x2[i] - x2[j]) * pr[i][j] / r[i][j];
                        }
                        //				if (j<=IM)
                        {
                            gr_x1[j] = gr_x1[j] - (x1[i] - x1[j]) * pr[i][j] / r[i][j];
                            gr_x2[j] = gr_x2[j] - (x2[i] - x2[j]) * pr[i][j] / r[i][j];
                        }
                        for (k = 1; k <= I; k++)
                            if ((k != i) && (k != j)) {
                                if ((br[i][j][k]) <= (RR + RC)) {
                                    //				if (j<=IM)
                                    {
                                        gr_x1[j] = gr_x1[j] + (x1[j] - x1[k]) * par[i][j][k] / ar[i][j][k];
                                        gr_x2[j] = gr_x2[j] + (x2[j] - x2[k]) * par[i][j][k] / ar[i][j][k];
                                    }
                                    //				if (k<=IM)
                                    {
                                        gr_x1[k] = gr_x1[k] - (x1[j] - x1[k]) * par[i][j][k] / ar[i][j][k];
                                        gr_x2[k] = gr_x2[k] - (x2[j] - x2[k]) * par[i][j][k] / ar[i][j][k];
                                    }
                                    //				if (i<=IM)
                                    {
                                        gr_x1[i] = gr_x1[i] + (x1[i] - x1[k]) * pbr[i][j][k] / br[i][j][k];
                                        gr_x2[i] = gr_x2[i] + (x2[i] - x2[k]) * pbr[i][j][k] / br[i][j][k];
                                    }
                                    //				if (k<=IM)
                                    {
                                        gr_x1[k] = gr_x1[k] - (x1[i] - x1[k]) * pbr[i][j][k] / br[i][j][k];
                                        gr_x2[k] = gr_x2[k] - (x2[i] - x2[k]) * pbr[i][j][k] / br[i][j][k];
                                    }
                                }
                            }
                    }
                }
            }

        //deleteArrays(MM);

        return (E);
    }

    
    void lazyAllocGradArrays(int II) {
        if(II != mII) {
            if(mII != 0) {
                freeGradArrays();
            }
            allocGradArrays(II);
            mII = II;
        }
    }
    void allocGradArrays(int II) {
        if ((Jz = (int *) malloc(II * sizeof (int))) == NULL) {
            printf("\n  Net pamyati \n");
            exit(10);
        }
        if ((Jzz = (int *) malloc(II * sizeof (int))) == NULL) {
            printf("\n  Net pamyati \n");
            exit(10);
        }
        if ((fad_d = (double *) malloc(II * sizeof (double))) == NULL) {
            printf("\n  Net pamyati \n");
            exit(10);
        }
        if ((fad_s = (double *) malloc(II * sizeof (double))) == NULL) {
            printf("\n  Net pamyati \n");
            exit(10);
        }
        if ((h = (double *) malloc(II * sizeof (double))) == NULL) {
            printf("\n  Net pamyati \n");
            exit(10);
        }
        if ((dm = (double *) malloc(II * sizeof (double))) == NULL) {
            printf("\n  Net pamyati \n");
            exit(10);
        }
        if ((dp = (double *) malloc(II * sizeof (double))) == NULL) {
            printf("\n  Net pamyati \n");
            exit(10);
        }
        if ((grad_h = (double *) malloc(II * sizeof (double))) == NULL) {
            printf("\n  Net pamyati \n");
            exit(10);
        }
        if ((grad_d = (double *) malloc(II * sizeof (double))) == NULL) {
            printf("\n  Net pamyati \n");
            exit(10);
        }
        if ((grad_s = (double *) malloc(II * sizeof (double))) == NULL) {
            printf("\n  Net pamyati \n");
            exit(10);
        }

    }

    void freeGradArrays() {
        if (grad_h != NULL) free(grad_h);
        if (grad_s != NULL) free(grad_s);
        if (grad_d != NULL) free(grad_d);
        if (Jz != NULL) free(Jz);
        if (Jzz != NULL) free(Jzz);
        if (h != NULL) free(h);
        if (fad_s != NULL) free(fad_s);
        if (fad_d != NULL) free(fad_d);
        if (dm != NULL) free(dm);
        if (dp != NULL) free(dp);
    }


public:

    double grad(int nlay, const struct TersoffParams& tsof, const int* lb, const int* ub, const double* x, double* g) {
        int ki, N, M, LM, MP, MM, II, I, k, n, i, j, m, num, ind, K;
        double aaa, xmin, xmax, Sum, RO;
        double RR, RC; // zadannye parametry

        I = nlay;
        II = I + 1;
        N = VARIABLE_TERSOFF_PARAMS; //   kolichestvo izmenyaemyh parametrov potenciala;

        lazyAllocGradArrays(II);

        RR = tsof.mR;
        RC = tsof.mRcut;

        h[0] = 0.0;
        ind = 0;
        for (i = 1; i <= I; i++) {
            Jz[i] = lb[i - 1];
            Jzz[i] = ub[i - 1];
            h[i] = x[ind];
            ind = ind + 1;
            fad_d[i] = x[ind];
            ind = ind + 1;
            fad_s[i] = x[ind];
            ind = ind + 1;
        }

        M = m_index(I, Jz, Jzz); //   kolichestvo atomov sistemy
        LM = M;
        //-------------------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------------
        //                              kolichestvo atomov ocrujeniya:
        m = M;
        RO = RR + RC;
        xmin = fad_d[1] - RO;
        xmax = fad_d[1] + (Jzz[1] - Jz[1]) * fad_s[1] + RO;
        for (i = 2; i <= I; i++) {
            if ((fad_d[i] - RO) < xmin) xmin = fad_d[i] - RO;
            if ((fad_d[i] + (Jzz[i] - Jz[i]) * fad_s[i] + RO) > xmax) xmax = fad_d[i] + (Jzz[i] - Jz[i]) * fad_s[i] + RO;
        }
        //----------------------------------------------------------------------------
        //                              OCRUJENIE po gorizontali
        for (i = 1; i <= I; i++) {
            dm[i] = fad_d[i];
            do {
                dm[i] = dm[i] - fad_s[i];
                m = m + 2;
            } while (dm[i] > xmin);
        }
        //------------------------------------------------------------------------------------
        Sum = 0.0;
        for (i = 1; i <= I; i++) Sum = Sum + h[i];
        ki = 0;
        do {
            ki = ki + 1;
        } while ((ki * Sum) < RO);
        K = ki;
        //------------------------------------------------------------------------------------

        for (k = 1; k <= K; k++) {
            for (i = 1; i <= I; i++) {
                dm[i] = fad_d[i];
                do {
                    dm[i] = dm[i] - fad_s[i];
                    m = m + 4;
                } while (dm[i] > xmin);
                for (j = Jz[i]; j <= Jzz[i]; j++) {
                    m = m + 2;
                }
            }
        }
        //-------------------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------------
        MP = m;
        MM = m + 1;

        lazyAllocArrays(MM);



        //-------------------------------------------------------------------------------
        m = 0;
        for (i = 1; i <= I; i++)
            for (j = Jz[i]; j <= Jzz[i]; j++) {
                m = m + 1;
                x1[m] = fad_d[i] + (j - Jz[i]) * fad_s[i];
                aaa = 0.0;
                for (k = 2; k <= i; k++) aaa = aaa + h[k];
                x2[m] = aaa;
                im[m] = i;
                jm[m] = j;
            }

        //-------------------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------------
        //                              PARAMETRY ATOMOV OCRUJENIYA:
        //----------------------------------------------------------------------------
        //                              OCRUJENIE po gorizontali
        for (i = 1; i <= I; i++) {
            aaa = 0.0;
            for (k = 2; k <= i; k++) aaa = aaa + h[k];

            dm[i] = fad_d[i];
            dp[i] = fad_d[i] + (Jzz[i] - Jz[i]) * fad_s[i];
            num = 0;
            do {
                dm[i] = dm[i] - fad_s[i];
                dp[i] = dp[i] + fad_s[i];
                num = num + 1;
                m = m + 1;
                x1[m] = dm[i];
                x2[m] = aaa;
                im[m] = i;
                g_s[m] = -num;
                g_h[m] = i;
                m = m + 1;
                x1[m] = dp[i];
                x2[m] = aaa;
                im[m] = i;
                g_s[m] = (Jzz[i] - Jz[i]) + num;
                g_h[m] = i;
            } while (dm[i] > xmin);
        }

        //------------------------------------------------------------------------------------
        for (k = 1; k <= K; k++) {
            for (i = 1; i <= I; i++) {
                aaa = -Sum*k;
                for (n = 2; n <= i; n++) aaa = aaa + h[n];

                dm[i] = fad_d[i];
                dp[i] = fad_d[i] + (Jzz[i] - Jz[i]) * fad_s[i];
                num = 0;
                do {
                    dm[i] = dm[i] - fad_s[i];
                    dp[i] = dp[i] + fad_s[i];
                    num = num + 1;
                    m = m + 1;
                    x1[m] = dm[i];
                    x2[m] = aaa;
                    im[m] = i;
                    g_s[m] = -num;
                    g_h[m] = -(k - 1) * I - (I - i + 1);
                    m = m + 1;
                    x1[m] = dp[i];
                    x2[m] = aaa;
                    im[m] = i;
                    g_s[m] = (Jzz[i] - Jz[i]) + num;
                    g_h[m] = -(k - 1) * I - (I - i + 1);
                } while (dm[i] > xmin);
                for (j = Jz[i]; j <= Jzz[i]; j++) {
                    m = m + 1;
                    x1[m] = fad_d[i] + (j - Jz[i]) * fad_s[i];
                    x2[m] = aaa;
                    im[m] = i;
                    g_s[m] = (j - Jz[i]);
                    g_h[m] = -(k - 1) * I - (I - i + 1);
                }
            }

            for (i = 1; i <= I; i++) {
                aaa = Sum * k - h[1];
                for (n = 1; n <= i; n++) aaa = aaa + h[n];

                dm[i] = fad_d[i];
                dp[i] = fad_d[i] + (Jzz[i] - Jz[i]) * fad_s[i];
                num = 0;
                do {
                    dm[i] = dm[i] - fad_s[i];
                    dp[i] = dp[i] + fad_s[i];
                    num = num + 1;
                    m = m + 1;
                    x1[m] = dm[i];
                    x2[m] = aaa;
                    im[m] = i;
                    g_s[m] = -num;
                    g_h[m] = k * I + i;
                    m = m + 1;
                    x1[m] = dp[i];
                    x2[m] = aaa;
                    im[m] = i;
                    g_s[m] = (Jzz[i] - Jz[i]) + num;
                    g_h[m] = k * I + i;
                } while (dm[i] > xmin);
                for (j = Jz[i]; j <= Jzz[i]; j++) {
                    m = m + 1;
                    x1[m] = fad_d[i] + (j - Jz[i]) * fad_s[i];
                    x2[m] = aaa;
                    im[m] = i;
                    g_s[m] = (j - Jz[i]);
                    g_h[m] = k * I + i;
                }
            }
        }
        //                             END OCRUJENIE;
        //-------------------------------------------------------------------------------------------------------
        //-------------------------------------------------------------------------------------------------------
        F = GRAD(MP, M, N, tsof, x1, x2, gr_x1, gr_x2);

        for (i = 1; i <= I; i++) grad_h[i] = 0.0;
        for (m = 1; m <= M; m++) {
            i = im[m];
            for (n = 2; n <= i; n++)
                grad_h[n] = grad_h[n] + gr_x2[m];
        }

        for (i = 1; i <= I; i++) {
            grad_d[i] = 0.0;
            for (j = Jz[i]; j <= Jzz[i]; j++) {
                for (m = 1; m <= M; m++) {
                    if ((im[m] == i)&&(jm[m] == j)) aaa = gr_x1[m];
                }
                grad_d[i] = grad_d[i] + aaa;
            }
        }

        for (i = 1; i <= I; i++) {
            grad_s[i] = 0.0;
            for (j = Jz[i]; j <= Jzz[i]; j++) {
                for (m = 1; m <= M; m++) {
                    if ((im[m] == i)&&(jm[m] == j)) aaa = gr_x1[m];
                }
                grad_s[i] = grad_s[i] + (j - Jz[i]) * aaa;
            }
        }

        for (m = M + 1; m <= MP; m++) {
            i = im[m];
            grad_d[i] = grad_d[i] + gr_x1[m];
            grad_s[i] = grad_s[i] + (g_s[m]) * gr_x1[m];
        }

        for (m = M + 1; m <= MP; m++) {
            i = im[m];

            if ((g_h[m] >= 2)&&(g_h[m] <= I)) {
                for (n = 2; n <= i; n++)
                    grad_h[n] = grad_h[n] + gr_x2[m];
            }
            if (g_h[m] > I) {
                grad_h[1] = grad_h[1] - gr_x2[m];
                for (n = 1; n <= I; n++)
                    grad_h[n] = grad_h[n] + ((g_h[m] - i) / I) * gr_x2[m];
                for (n = 1; n <= i; n++)
                    grad_h[n] = grad_h[n] + gr_x2[m];
            }
            if (g_h[m] < 0) {
                grad_h[1] = grad_h[1] - gr_x2[m];
                for (n = 1; n <= I; n++)
                    grad_h[n] = grad_h[n] + ((g_h[m] + I - i + 1) / I) * gr_x2[m];
                for (n = i + 1; n <= I; n++)
                    grad_h[n] = grad_h[n] - gr_x2[m];
            }
        }

        ind = 0;
        for (i = 1; i <= I; i++) {
            g[ind] = grad_h[i];
            ind = ind + 1;
            g[ind] = grad_d[i];
            ind = ind + 1;
            g[ind] = grad_s[i];
            ind = ind + 1;
        }

        //freeGradArrays();

        return (F);

    }

};


#endif /* TSOFFADGRAD_HPP */

