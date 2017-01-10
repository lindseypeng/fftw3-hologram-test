#include "complexnumber.h"
#include <complex>
#include <fftw3.h>
#include <iostream>
#include <math.h>

using namespace std;

int complexnumber(int step1,int width1, int height1)
{
    const complex<double> J (0.0,1.0);
    const complex<double> d (130.0,0.0);
    const complex<double> PI(M_PI,0.0);//distance in mm


    double dx=0.005195;
    double dy=0.005203125;
    double lambda0=0.000488;

    double XX[width1*height1];
    double YY[width1*height1];

    complex <double> num[width1*height1];
    complex <double> den[width1*height1];
    complex <double> g[width1*height1];

//for YY row major loop
    for (int j=0,k=0;j<height1;j++){
        for (int i=0;i<width1;i++){
        YY[k]=j*dy;
        k++;
        }
           }
//for XX row major loop

    for (int j=0,k=0;j<height1;j++){
        for (int i=0;i<width1;i++){
        XX[k]=i*dx;
        k++;
        }
    }

//creating row major array of transfer function g
    for (int j=0;j<height1*width1;j++){
        num[j]=exp(J*2.0*M_PI/lambda0*sqrt(pow(d,2.0)+pow(XX[j],2.0)+pow(YY[j],2.0)));
        den[j]=sqrt(pow(d,2.0)+pow(XX[j],2.0)+pow(YY[j],2.0));
        g[j]=-J/lambda0*num[j]/den[j];
    }

//converting complex(double) to fftw_complex double to *in
        for (int j=0;j<height1*width1;j++)
    {
        fftw_complex* in = reinterpret_cast<fftw_complex*>(&num[j]);
    }

     return 0;
    }

