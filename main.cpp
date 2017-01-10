#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>
#include"complexnumber.h"

using namespace cv;
using namespace std;


int main(void)
{
    cv::Mat img1;
    cv::Mat img2;
    uchar *img1_data;
    uchar *img2_data;


    fftw_complex    *data_in;
    fftw_complex    *fft1;
    fftw_complex    *fft2;//added
    fftw_complex    *products;//added
    fftw_complex    *ifft;
    fftw_plan       plan_f;
    fftw_plan       plan_g;//added
    fftw_plan       plan_b;

    int             width, height, step;
    int             i, j, k;



    // load original image
    img1= imread("/home/alinsi/codeblocks/fft2/3beads2.tif", CV_LOAD_IMAGE_GRAYSCALE);


    // create new image for IFFT result
    img2 = img1.clone();

    // get image properties
    width  	  = img1.size().width;
    height 	  = img1.size().height;
    step	  = img1.step;
    img1_data =  img1.data;
    img2_data =  img2.data;



    // initialize arrays for fftw operations
    data_in = fftw_alloc_complex(width * height);//raw image data
    fft1    = fftw_alloc_complex(width * height);//fourier transform of image
    fft2    = fftw_alloc_complex(width * height);//fourier transform of transfer function.added
    products= fftw_alloc_complex(width * height);//product of fft1 and fft2.added
    ifft    = fftw_alloc_complex(width * height);//inverse fft of products



    // create plans
    plan_f = fftw_plan_dft_1d( width * height, data_in, fft1,  FFTW_FORWARD,  FFTW_ESTIMATE );
    plan_g = fftw_plan_dft_1d( width * height, in, fft2,  FFTW_FORWARD,  FFTW_ESTIMATE );//added.in
    plan_b = fftw_plan_dft_1d( width * height, products,     ifft, FFTW_BACKWARD, FFTW_ESTIMATE );//changed fft to products

    // load img1's data to fftw input
    for( i = 0, k = 0 ; i < height ; i++ ) {
        for( j = 0 ; j < width ; j++ ) {
            data_in[k][0] = ( double )img1_data[i * step + j];
            data_in[k][1] = 0.0;
            k++;
        }
    }

    // perform FFT
    fftw_execute( plan_f );
    //perform FFT on products
    fftw_execute( plan_g ); //added


    //now multiple FFT1  FFT2/////////double check multiplication
    for( i = 0, k = 0 ; i < height ; i++ ) {
        for( j = 0 ; j < width ; j++ ) {
           products[k][0] = ( double )img1_data[i * step + j];

            k++;
        }
    }


    // perform IFFT
    fftw_execute( plan_b );

    // normalize IFFT result
    for( i = 0 ; i < ( width * height ) ; i++ ) {
        ifft[i][0] /= ( double )( width * height );
    }

    // copy IFFT result to img2's data
    for( i = 0, k = 0 ; i < height ; i++ ) {
		for( j = 0 ; j < width ; j++ ) {
			img2_data[i * step + j] = ( uchar )ifft[k++][0];
		}
	}

    // display images
    cv::namedWindow( "original_image", CV_WINDOW_AUTOSIZE );
    cv::namedWindow( "IFFT", CV_WINDOW_AUTOSIZE );
    cv::imshow( "original_image", img1 );
    cv::imshow( "IFFT", img2 );

    char key;
    while (true) {
        key = cv::waitKey( 0 );
        if (27 == key)
            break;
    }

    // free memory
    cv::destroyWindow( "original_image" );
    cv::destroyWindow( "IFFT" );
    img1.release();
    img2.release();
    fftw_destroy_plan( plan_f );
    fftw_destroy_plan( plan_b );
    fftw_free( data_in );
    fftw_free( fft );
    fftw_free( ifft );

    return 0;
}
