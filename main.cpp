


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
    cv::Mat img2;
    cv::Mat img3;
    cv::Mat img4;
    uchar *img1_data;
    uchar *img2_data;
    uchar *img3_data;
    uchar *img4_data;

    fftw_complex    *data_in;
    fftw_complex    *fft;
    fftw_complex    *ifft;
    fftw_complex    *transfer;
    fftw_complex    *products;
    fftw_complex    *inverse;

    //try memcpy tranfer g to another variale
    fftw_complex    *gg;

    fftw_plan       plan_f;
    fftw_plan       plan_b;
    fftw_plan       plan_g;
    fftw_plan       plan_I;


    int             width, height, step;
    int             i, j, k;


    const complex<double> J (0.0,1.0);
    const complex<double> d (30.0,0.0);
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
    img2 = img1.clone();
    img3 = img1.clone();
    img4 = img4.clone();

    // get image properties
    width  	  = img1.size().width;
    height 	  = img1.size().height;
    step	  = img1.step;
    img1_data = ( uchar* ) img1.data;
    img2_data = ( uchar* ) img2.data;
    img3_data = ( uchar* ) img3.data;
    img4_data = ( uchar* ) img4.data;
    // initialize arrays for fftw operations
    data_in = fftw_alloc_complex(width * height);
    fft     = fftw_alloc_complex(width * height);
    ifft    = fftw_alloc_complex(width * height);
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
//checking gg values--good to go
//for (int i=1306870; i<1306878; i++){
//
//    cout<<"tranfer function gg  "<<gg[i][0]<<" + "<<gg[i][1]<<endl;
//
//}
 // perform FFT
    fftw_execute( plan_f );

//normalize fft

//    for( i = 0 ; i < ( width * height ) ; i++ ) {
//        fft[i][0] /=(( double )( width * height ));
//        fft[i][1] /= (( double )( width * height ));
//
//    }

//checking values of original in fourier space

    for (int i=1306870; i<1306878; i++){
        cout<<" fft "<<fft[i][0]<<" + "<<fft[i][1]<<endl;
    }

//create another plan

//plan_g = fftw_plan_dft_2d( width, height, reinterpret_cast<fftw_complex*>(g), transfer,FFTW_FORWARD, FFTW_ESTIMATE);
plan_g = fftw_plan_dft_2d( height, width, gg, transfer,FFTW_FORWARD, FFTW_ESTIMATE);
plan_b = fftw_plan_dft_2d( height, width, transfer,     ifft, FFTW_BACKWARD, FFTW_ESTIMATE );
    //perform another FFT
    fftw_execute( plan_g );
//normalize transfer results

//    for( i = 0 ; i < ( width * height ) ; i++ ) {
//        transfer[i][0] /= (( double )( width * height ));
//        transfer[i][1] /= (( double )( width * height ));
//
//    }
//checking fourier transfer numers and size
for (int i=1306870; i<1306878; i++){

    cout<<"tranfer function gg  "<<transfer[i][0]<<" + "<<transfer[i][1]<<endl;

}



//now multiple FFT1  FFT2/////////double check multiplication of copmlex numbers(xu-yv)+(xv+yu)i
     for( k = 0 ; k < (height*width); k++ ) {

            products[k][0] = fft[k][0]*transfer[k][0]-fft[k][1]*transfer[k][1];
            products[k][1] = fft[k][0]*transfer[k][1]+fft[k][1]*transfer[k][0];
            k++;
        }

//create another plan
plan_I = fftw_plan_dft_2d( height, width, products,inverse,FFTW_BACKWARD, FFTW_ESTIMATE);


//test products//
for (int i=0; i<10; i++){
        cout<<" products "<<products[i][0]<<" + "<<products[i][1]<<endl;
    }
//perform IFFT inverse on products

    fftw_execute( plan_I );

    //    // normalize final inverse result
    for( i = 0 ; i < ( width * height ) ; i++ ) {
        inverse[i][0] /= (( double )( width * height ));
        inverse[i][1] /= (( double )( width * height ));

    }

//for (int i=0; i<10; i++){
//
//        cout<<" transfer fourier "<<ifft[i][0]<<ifft[i][0]<<endl;
//        cout<<" original transfer "<<ifft[i][1]<<" + "<<ifft[i][1]<<endl;
//
//
//    }

//	double *floatPtr = img4.ptr<double>();
//    //add g into img4
//    for( i = 0, k = 0 ; i < height ; i++ ) {
//		for( j = 0 ; j < width ; j++ ) {
//			*floatPtr++ = (uchar)(g[k].real());
//			k++;
//		}
//	}


	    // normalize inverse result
//    for( i = 0 ; i < ( width * height ) ; i++ ) {
//        inverse[i][0] /= ( double )( width * height );
//        inverse[i][1] /= ( double )( width * height );
//        //log scale//
//    }

////calculating absolute values
////    for( i = 0 ; i < ( width * height ) ; i++ ) {
////        final2[i]=pow(inverse[i][0],2.0)+pow(inverse[i][1],2.0);//log scale//
//// }
//// normalize inverse result
////    for( i = 0 ; i < ( width * height ) ; i++ ) {
////            final2[i] /= (double) (width*height);//log scale//
////    }
    // perform IFFT 
    fftw_execute( plan_b );



//     normalize IFFT result
    for( i = 0 ; i < ( width * height ) ; i++ ) {
        ifft[i][0] /= (( double )( width * height ));
        ifft[i][1] /= (( double )( width * height ));

    }

    // copy IFFT result to img2's data
    for( i = 0, k = 0 ; i < height ; i++ ) {
		for( j = 0 ; j < width ; j++ ) {
			img2_data[i * step + j] = ( uchar )ifft[k++][1];
		}
	}
	    //copy real components of transfer data into img3
    for( i = 0, k = 0 ; i < height ; i++ ) {
		for( j = 0 ; j < width ; j++ ) {
			img3_data[i * step + j] = ( uchar )inverse[k++][0];
		}
	}

//		uchar c = 1/log(1+255);
//		cv::Mat fg;
//        cv::Mat j2;
//fg=img1.clone;
//    for( i = 0, k = 0 ; i < height ; i++ ) {
//		for( j = 0 ; j < width ; j++ ) {
//			fg[i * step + j]=( uchar )inverse[k++][0]*255;
//			j2=c*log*(1+fg);
//		}
//	}

//    cv::Mat fg;
//    img3.convertTo(fg,CV_32F);
//    fg=fg+1;
//    log(fg,fg);
//    convertScaleAbs(fg,fg);
	 //print out the minmax of ifft before normalization
    double minv;
    double maxv;
    Point minL;
    Point maxL;

    minMaxLoc(img3,&minv,&maxv,&minL,&maxL);
    cout<<"min val of img3:"<<minv<<endl;
    cout<<"max val of img3:"<<maxv<<endl;



//
//
    // display images
    cv::namedWindow( "original_image", CV_WINDOW_AUTOSIZE );
    cv::namedWindow( "Transfer function", CV_WINDOW_AUTOSIZE );
    cv::namedWindow( "Inverse", CV_WINDOW_AUTOSIZE );
    //cv::namedWindow("transfer",CV_WINDOW_AUTOSIZE);
    cv::imshow( "original_image", img1 );
    cv::imshow( "Transfer function", img2 );
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

    img1.release();
    img2.release();
    img3.release();
    img4.release();
    fftw_destroy_plan( plan_f );
    fftw_destroy_plan( plan_b );
    fftw_destroy_plan( plan_I );
    fftw_destroy_plan( plan_g );

    fftw_free( data_in );
    fftw_free( fft );
    fftw_free( ifft );

    fftw_free( products );
    fftw_free( inverse);
    fftw_free(transfer);

    void fftw_cleanup_threads(void);

    return 0;
}
