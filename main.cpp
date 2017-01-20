#include <stdlib.h>
#include <math.h>
#include <complex>
#include <fftw3.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include <iomanip>


using namespace cv;
using namespace std;
/////////////allocating space y defining gloal variales for transfer function//////////////

/////////////////////////////////


int main()
{
    cv::Mat img1;
    cv::Mat img3;

    uchar *img1_data;
    uchar *img3_data;

    fftw_complex    *data_in;
    fftw_complex    *fft;
    fftw_complex    *transfer;
    fftw_complex    *products;
    fftw_complex    *inverse;

    //try memcpy tranfer g to another variale
    fftw_complex    *gg;

    fftw_plan       plan_f;
    fftw_plan       plan_g;
    fftw_plan       plan_I;

    int             width, height, step;
    int             i, j, k;

    const complex<double> J (0.0,1.0);
    const complex<double> d (130.0,0.0);
    const complex<double> PI(M_PI,0.0);//distance in mm
    double dx=5.32/1024;
    double dy=6.66/1280;
    double lambda0=0.000488;

    complex <double> *num = new complex<double>[1024*1280];
    complex <double> *den = new complex<double>[1024*1280] ;
    complex <double> *g = new complex<double>[1024*1280];

    double *XX =new double[1024*1280] ;
    double *YY = new double[1024*1280] ;
    double *final2=new double[1024*1280] ;
    // load original image
    img1= imread("/home/alinsi/codeblocks/fft2/microbead1.tif", CV_LOAD_IMAGE_GRAYSCALE);

    // create new image for IFFT result
    img3 = img1.clone();

    // get image properties
    width  	  = img1.size().width;
    height 	  = img1.size().height;
    step	  = img1.step;
    img1_data = ( uchar* ) img1.data;
    img3_data = ( uchar* ) img3.data;

    // initialize arrays for fftw operations
    data_in = fftw_alloc_complex(width * height);
    fft     = fftw_alloc_complex(width * height);
    transfer = fftw_alloc_complex(width * height);
    products =fftw_alloc_complex(width * height);
    inverse = fftw_alloc_complex(width * height);
    gg = fftw_alloc_complex(width*height);

    // create plans
    plan_f = fftw_plan_dft_2d( height, width, data_in, fft,  FFTW_FORWARD,  FFTW_ESTIMATE );

    // load img1's data to fftw input
    for( i = 0, k = 0 ; i < height ; i++ ) {
        for( j = 0 ; j < width ; j++ ) {
            data_in[k][0] = ( double )img1_data[i * step + j];
            data_in[k][1] = 0.0;
            k++;
        }
    }

    ////calculating transfer function//
//for YY row major loop
    for (int j=-height/2,k=0;j<height/2;j++){
        for (int i=-width/2;i<width/2;i++){
        YY[k]=j*dy;
        k++;
        }
           }
//for XX row major loop
    for (int j=-height/2,k=0;j<height/2;j++){
        for (int i=-width/2;i<width/2;i++){
        XX[k]=(i)*dx;
        k++;
        }
    }

//creating row major array of transfer function g
    for (int j=0;j<height*width;j++){
        num[j]=-exp(J*2.00000000*M_PI/lambda0*sqrt(pow(d,2.00000)+pow(XX[j],2.00000)+pow(YY[j],2.00000)));
        den[j]=sqrt(pow(d,2.0)+pow(XX[j],2.0)+pow(YY[j],2.0));
        g[j]=-J/lambda0*num[j]/den[j];

    }
//copy the result from complex<doule> to fftw_complex
memcpy( &gg, &g, sizeof (fftw_complex));

 // perform FFT
    fftw_execute( plan_f );

//plan_g = fftw_plan_dft_2d( width, height, reinterpret_cast<fftw_complex*>(g), transfer,FFTW_FORWARD, FFTW_ESTIMATE);
plan_g = fftw_plan_dft_2d( height, width, gg, transfer,FFTW_FORWARD, FFTW_ESTIMATE);

    //perform another FFT
    fftw_execute( plan_g );

//now multiple FFT1  FFT2/////////double check multiplication of copmlex numbers(xu-yv)+(xv+yu)i
     for( i=0,k = 0 ; i < (height*width); i++ ) {

            products[k][0] = (fft[k][0]*transfer[k][0]-fft[k][1]*transfer[k][1]);
            products[k][1] = (fft[k][0]*transfer[k][1]+fft[k][1]*transfer[k][0]);
            k++;
        }

//create another plan
plan_I = fftw_plan_dft_2d( height, width, products,inverse,FFTW_BACKWARD, FFTW_ESTIMATE);

    fftw_execute( plan_I );
    //    // normalize final inverse result
    for( i = 0 ; i < ( width * height ) ; i++ ) {
        inverse[i][0] /= (( double )( width * height ));
        inverse[i][1] /= (( double )( width * height ));

    }

//obtain magnitude of inverse
    for( i = 0 ; i < ( width * height ) ; i++ ) {
    final2[i]=sqrt(inverse[i][0]*inverse[i][0]+inverse[i][1]*inverse[i][1]);

    }
//otain maximum value from inverse array
    double max2=final2[0];
    double min2=final2[0];

    for( i = 0; i < ( width * height ) ; i++ ) {
        if (final2[i]>max2)
            max2=final2[i];
        if (final2[i]<min2)
            min2=final2[i];
    }
//normalize final2 into range of 0 and 255
    for (i=0,k=0;i<width*height;i++){
        final2[k] /= max2;
        final2[k] *= 255;
        k++;
        }
	    //copy  transfer data into img3
    for( i = 0, k = 0 ; i < height ; i++ ) {
		for( j = 0 ; j < width ; j++ ) {
			img3_data[i * step + j] = ( uchar )(final2[k++]);
		}
	}

// rearrange the quadrants of Fourier image
// so that the origin is at the image center
int cx = width/2;
int cy = height/2;

Mat tmp;
Mat q0(img3, Rect(0, 0, cx, cy));
Mat q1(img3, Rect(cx, 0, cx, cy));
Mat q2(img3, Rect(0, cy, cx, cy));
Mat q3(img3, Rect(cx, cy, cx, cy));
q0.copyTo(tmp);
q3.copyTo(q0);
tmp.copyTo(q3);

Mat tmp1;
q1.copyTo(tmp1);
q2.copyTo(q1);
tmp1.copyTo(q2);


//
//
    // display images
    cv::namedWindow( "original_image", CV_WINDOW_NORMAL );
    cv::namedWindow( "Inverse", CV_WINDOW_NORMAL );
    //cv::namedWindow("transfer",CV_WINDOW_AUTOSIZE);
    cv::imshow( "original_image", img1 );
    cv::imshow( "Inverse", img3 );
  //  cv::imshow("transfer",img4);
    char key;
    while (true) {
        key = cv::waitKey( 0 );
        if (27 == key)
            break;
    }

    // free memory
    cv::destroyWindow( "original_image" );
    cv::destroyWindow( "IFFT" );
    cv::destroyWindow( "Inverse");
    //cv::destroyWindow( "transfer");

    delete [] num;
    delete [] den;
    delete [] g;
    delete [] XX;
    delete [] YY;
    delete [] final2;
   //delete [] temporary;

    img1.release();
    img3.release();

    fftw_destroy_plan( plan_f );
    fftw_destroy_plan( plan_I );
    fftw_destroy_plan( plan_g );

    fftw_free( data_in );
    fftw_free( fft );
    fftw_free( products );
    fftw_free( inverse);
    fftw_free(transfer);

    void fftw_cleanup_threads(void);

    return 0;
}
