#include "ffttools.hpp"

namespace FFTTools{
    
cv::Mat fftd1d(cv::Mat img, int dim, bool backwards)
{
    if (img.channels() == 1)
    {
        cv::Mat planes[] = {cv::Mat_<float> (img), cv::Mat_<float>::zeros(img.size())};
        //cv::Mat planes[] = {cv::Mat_<double> (img), cv::Mat_<double>::zeros(img.size())};
        cv::merge(planes, 2, img);
    }
    if (0==dim)
        cv::dft(img, img, backwards ? (cv::DFT_ROWS|cv::DFT_INVERSE | cv::DFT_SCALE) : cv::DFT_ROWS );
    else
    {
        img = img.t();
        cv::dft(img, img, backwards ? (cv::DFT_ROWS|cv::DFT_INVERSE | cv::DFT_SCALE) : cv::DFT_ROWS );
        img = img.t();
    }
    return img;
}
 //   CvMat * matCircshift( CvMat *mat,int rowMove,int colMove) //矩阵平移
//{
// int i,j,k;
 
// if(rowMove>=0) //rowMove>=0 矩阵下移
// {
//  CvMat * temp_rowMat=cvCreateMat(1,mat->cols,CV_32FC1);
//  CvMat * rowMat=cvCreateMat(1,mat->cols,CV_32FC1);
//  while(rowMove--)
//  {
//   cvGetRow(mat,temp_rowMat,(mat->rows)-1);
//   cvCopy(temp_rowMat,rowMat,NULL );
//   for(i=mat->rows-1;i>0;i--)
//   {
//    for(j=0;j<mat->cols;j++)
//    {
//     cvmSet(mat,i,j,cvmGet(mat,i-1,j));
//    }
//    }
//   }
//   for(k=0;k<mat->cols;k++)
//   {
//    cvmSet(mat,0,k,cvmGet(rowMat,0,k));
//   }
//  }
// }
// else  //rowMove<0 矩阵上移
// {
//  rowMove=-rowMove;
//  CvMat* temp_rowMat=cvCreateMat(1,mat->cols,CV_32FC1);
//  CvMat* rowMat=cvCreateMat(1,mat->cols,CV_32FC1);
//  while(rowMove--)/
//  {
//   cvGetRow(mat,temp_rowMat,0);
// cvCopy(temp_rowMat,rowMat,NULL);
//   for(i=0;i<mat->rows-1;i++)
//   {
//    for(j=0;j<mat->cols;j++)
//    {
//     cvmSet(mat,i,j,cvmGet(mat,i+1,j));/
//    }
//   }
//   for(k=0;k<mat->cols;k++)
//   {
//    cvmSet(mat,mat->rows-1,k,cvmGet(rowMat,0,k));
//   }
//  }
// }
   
// if(colMove>=0) //colMove>=0 右移
// {
//  CvMat* temp_colMat=cvCreateMat(mat->rows,1,CV_32FC1);
//  CvMat* colMat=cvCreateMat(mat->rows,1,CV_32FC1);
//  while(colMove--)
//  {
//   cvGetCol(mat,temp_colMat,(mat->cols)-1);
//   cvCopy(temp_colMat,colMat,NULL);
//   for(i=mat->cols-1;i>0;i--)
//   {
//    for(j=0;j<mat->rows;j++)
//    {
//     cvmSet(mat,j,i,cvmGet(mat,j,i-1));
//    }
//   }
//   for(k=0;k<mat->rows;k++)
//   {
//    cvmSet(mat,k,0,cvmGet(colMat,k,0));
 //  }
 // }
// }
// else //colMove<0 左移
// {
//  colMove=-colMove;
// CvMat* colMat=cvCreateMat(mat->rows,1,CV_32FC1);
//  CvMat* temp_colMat=cvCreateMat(mat->rows,1,CV_32FC1);
//  while(colMove--)
//  {
//   cvGetCol(mat,temp_colMat,0);
//   cvCopy(temp_colMat,colMat,NULL);/
//   for(i=0;i<mat->cols-1;i++)
//   {
//    for(j=0;j<mat->rows;j++)
//    {
//     cvmSet(mat,j,i,cvmGet(mat,j,i+1));
//    }
//   }
//   for(k=0;k<mat->rows;k++)
//   {
//    cvmSet(mat,k,mat->cols-1,cvmGet(colMat,k,0));
//   }
//  }
// }
// return mat;/
//}
	
	
cv::Mat fftd(cv::Mat img, bool backwards)
{
/*
#ifdef USE_FFTW

    fftw_complex * fm = (fftw_complex*) fftw_malloc(sizeof (fftw_complex) * img.cols * img.rows);

    fftw_plan p = fftw_plan_dft_2d(img.rows, img.cols, fm, fm, backwards ? 1 : -1, 0 * FFTW_ESTIMATE);


    if (img.channels() == 1)
    {
        for (int i = 0; i < img.rows; i++)
            for (int j = 0; j < img.cols; j++)
            {
                fm[i * img.cols + j][0] = img.at<float>(i, j);
                fm[i * img.cols + j][1] = 0;
            }
    }
    else
    {
        assert(img.channels() == 2);
        for (int i = 0; i < img.rows; i++)
            for (int j = 0; j < img.cols; j++)
            {
                fm[i * img.cols + j][0] = img.at<cv::Vec2d > (i, j)[0];
                fm[i * img.cols + j][1] = img.at<cv::Vec2d > (i, j)[1];
            }
    }
    fftw_execute(p);
    cv::Mat res(img.rows, img.cols, CV_64FC2);


    for (int i = 0; i < img.rows; i++)
        for (int j = 0; j < img.cols; j++)
        {
            res.at<cv::Vec2d > (i, j)[0] = fm[i * img.cols + j][0];
            res.at<cv::Vec2d > (i, j)[1] = fm[i * img.cols + j][1];

            //  _iout(fm[i * img.cols + j][0]);
        }

    if (backwards)res *= 1.d / (float) (res.cols * res.rows);

    fftw_free(p);
    fftw_free(fm);
    return res;

#else
*/
    if (img.channels() == 1)
    {
        cv::Mat planes[] = {cv::Mat_<float> (img), cv::Mat_<float>::zeros(img.size())};
        //cv::Mat planes[] = {cv::Mat_<double> (img), cv::Mat_<double>::zeros(img.size())};
        cv::merge(planes, 2, img);
    }
    cv::dft(img, img, backwards ? (cv::DFT_INVERSE | cv::DFT_SCALE) : 0 );

    return img;

/*#endif*/

}

cv::Mat real(cv::Mat img)
{
    std::vector<cv::Mat> planes;
    cv::split(img, planes);
    return planes[0];
}

cv::Mat imag(cv::Mat img)
{
    std::vector<cv::Mat> planes;
    cv::split(img, planes);
    return planes[1];
}

cv::Mat magnitude(cv::Mat img)
{
    cv::Mat res;
    std::vector<cv::Mat> planes;
    cv::split(img, planes); // planes[0] = Re(DFT(I), planes[1] = Im(DFT(I))
    if (planes.size() == 1) res = cv::abs(img);
    else if (planes.size() == 2) cv::magnitude(planes[0], planes[1], res); // planes[0] = magnitude
    else assert(0);
    return res;
}

cv::Mat complexMultiplication(cv::Mat a, cv::Mat b)
{
    std::vector<cv::Mat> pa;
    std::vector<cv::Mat> pb;
    cv::split(a, pa);
    cv::split(b, pb);

    std::vector<cv::Mat> pres;
    pres.push_back(pa[0].mul(pb[0]) - pa[1].mul(pb[1]));
    pres.push_back(pa[0].mul(pb[1]) + pa[1].mul(pb[0]));

    cv::Mat res;
    cv::merge(pres, res);

    return res;
}

cv::Mat complexConjMult(cv::Mat a, cv::Mat b)
{
    // b is conj
    std::vector<cv::Mat> pa;
    std::vector<cv::Mat> pb;
    cv::split(a, pa);
    cv::split(b, pb);
    std::vector<cv::Mat> pres;
    pres.push_back(pa[0].mul(pb[0]) + pa[1].mul(pb[1]));
    pres.push_back(-pa[0].mul(pb[1]) + pa[1].mul(pb[0]));
    cv::Mat res;
    cv::merge(pres, res);
    return res;
}

cv::Mat complexSelfConjMult(cv::Mat a)
{
    std::vector<cv::Mat> pa;
    cv::split(a, pa);
    cv::Mat res = pa[0].mul(pa[0]) + pa[1].mul(pa[1]);
    return res;
}

cv::Mat complexDivision(cv::Mat a, cv::Mat b)
{
    std::vector<cv::Mat> pa;
    std::vector<cv::Mat> pb;
    cv::split(a, pa);
    cv::split(b, pb);

    cv::Mat divisor = 1. / (pb[0].mul(pb[0]) + pb[1].mul(pb[1]));

    std::vector<cv::Mat> pres;

    pres.push_back((pa[0].mul(pb[0]) + pa[1].mul(pb[1])).mul(divisor));
    pres.push_back((pa[1].mul(pb[0]) + pa[0].mul(pb[1])).mul(divisor));

    cv::Mat res;
    cv::merge(pres, res);
    return res;
}

cv::Mat complexDivReal(cv::Mat cMat, cv::Mat rMat)
{
    cv::Mat res = cMat.clone();
    float *pRes = (float *)(res.data);
    float *pR = (float *)(rMat.data);
    for (int h=0; h<cMat.rows; h++)
    {
        for (int w=0; w<rMat.cols; w++)
        {
            pRes[0]= pRes[0]/pR[0];
            pRes[1]= pRes[1]/pR[0]; 
            pRes += 2;
            pR += 1;
        }
    }
    return res;
}

void rearrange(cv::Mat &img)
{
    // img = img(cv::Rect(0, 0, img.cols & -2, img.rows & -2));
    int cx = img.cols / 2;
    int cy = img.rows / 2;

    cv::Mat q0(img, cv::Rect(0, 0, cx, cy)); // Top-Left - Create a ROI per quadrant
    cv::Mat q1(img, cv::Rect(cx, 0, cx, cy)); // Top-Right
    cv::Mat q2(img, cv::Rect(0, cy, cx, cy)); // Bottom-Left
    cv::Mat q3(img, cv::Rect(cx, cy, cx, cy)); // Bottom-Right

    cv::Mat tmp; // swap quadrants (Top-Left with Bottom-Right)
    q0.copyTo(tmp);
    q3.copyTo(q0);
    tmp.copyTo(q3);
    q1.copyTo(tmp); // swap quadrant (Top-Right with Bottom-Left)
    q2.copyTo(q1);
    tmp.copyTo(q2);
}
/*
template < typename type>
cv::Mat fouriertransFull(const cv::Mat & in)
{
    return fftd(in);

    cv::Mat planes[] = {cv::Mat_<type > (in), cv::Mat_<type>::zeros(in.size())};
    cv::Mat t;
    assert(planes[0].depth() == planes[1].depth());
    assert(planes[0].size == planes[1].size);
    cv::merge(planes, 2, t);
    cv::dft(t, t);

    //cv::normalize(a, a, 0, 1, CV_MINMAX);
    //cv::normalize(t, t, 0, 1, CV_MINMAX);

    // cv::imshow("a",real(a));
    //  cv::imshow("b",real(t));
    // cv::waitKey(0);

    return t;
}*/

void normalizedLogTransform(cv::Mat &img)
{
    img = cv::abs(img);
    img += cv::Scalar::all(1);
    cv::log(img, img);
    // cv::normalize(img, img, 0, 1, CV_MINMAX);
}

}
