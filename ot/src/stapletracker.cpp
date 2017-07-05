#include <math.h>
#include <slapi-plugin.h>
#include "stapletracker.hpp"
#include "ffttools.hpp"
#include "recttools.hpp"
#include "fhog.hpp"
#include "dollar_hog.hpp"
#include "labdata.hpp"
#include "comdef.h"
#include "integral.hpp"
#include "labdata.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "shift.hpp"


static float subPixelPeak(float left, float center, float right)
{
    float divisor = 2*center - right - left;

    if (divisor == 0)
        return 0;
    return 0.5 * (right-left)/divisor;
}

StapleTracker::StapleTracker(int maxSide, int minSide,
                             int maxObjNum)
    :Tracker(maxSide, minSide, maxObjNum)
{
    m_trans_cell_size=4; //HoG cell size
    m_trans_inner_padding = 0.2f;
    m_trans_fixed_area=150*150; //
    m_trans_color_bins = 16;
    m_trans_lr_pwp = 0.04f;
    m_trans_lr_cf = 0.01f;
    m_trans_lambda = 1e-3f; //regularization factor
    m_trans_merge_factor = 0.3f;
    m_trans_y_sigma = 1.f/16;

    //Scale parameter
    m_scale_num = 17;//orginal 33 ,m_scale_num = 5;
    m_scale_cell_size = 4;
    m_scale_max_area = 32*16;
    m_scale_lr = 0.025;
    m_scale_lambda = 1e-3f;
    m_scale_step = 1.02; //orginal 1.02 1.01
    m_scale_y_sigma = 1.f/4;
    m_scale_hann = 0;
    m_scale_factors = 0;
    // by zz in 2017/05/15
    BACF_lr=0.0125;
     pfcell_size = 4;

    MaxItr=2;
    frameW=480;
    frameH=640;

    //
    //by zhangzheng in 2017/5/22
    //pvisualzation=1;
     visualization=0;
    debug=0;
    search_area_scale=4;
    filter_size=1.2000;
    output_sigma_factor=0.0625;
    ini_imgs=8;
    etha=0.0125;// 0.0125;
    upResize=100;
    lowResize=50;
    term=0.000001;
    lambda=0.001;
    slambda=0.0100;
    mu=1;
    maxMu=1000;
    search_area_shape="sequre";
    search_area_scale=4;
    min_image_sample_size=40000;
    max_image_sample_size=90000;
    fix_model_size=6400;
     pfcolorUpdataRate=0.0100;
     pft[0]="greyHoG";
     pft[1]="grey";
     pfhog_orientation=9;
     pfnbin=10;
     cell_size=4;

     pfinterPatchRate=0.3;
     pfgrey=1;
     pfgreyHoG=1;
     pfcolorProb=0;
     pfcolorProbHoG=0;
     pfcolorName=0;
     pfgreyProb=0;
     pflbp=0;
     nScales=17;
     nScalesInterp=33;
     scale_step=1.02;
     scale_sigma_factor=0.0625;
     scale_model_factor=1;
    scale_model_max_area=512;
     s_num_compressed_dim="MAX";
    // interp_factor=0.0250;
     search_size[0] =1;
     search_size[1]=0.9850;
     search_size[2]=0.9900;
     search_size[3]=0.9950;
     search_size[4]=1.0050;
     search_size[5]=1.0100;
     search_size[6]=1.0150;
     resize_image=0;
     resize_scale=1;
    pmax_scale_dim=1;
}

int StapleTracker::init()
{
    if (0!=Tracker::init())
        return -1;

    if (0 != initTransCF())
        return -1;

    if (0 != initScaleCF())
        return -1;

    m_cfs.resize(m_maxObjNum);
    return 0;
}

int StapleTracker::initTransCF()
{
    return 0;
}

int StapleTracker::initScaleCF()
{
    if (0==m_scale_hann)
    {
        m_scale_hann = new float[m_scale_num];
        if (0==m_scale_hann)
            return -1;
        for (int i = 0; i < m_scale_num; i++)
            m_scale_hann[i] = 0.5*(1-std::cos(2*PI_T*i/(m_scale_num-1)));
    }

    if (0==m_scale_factors) {
        m_scale_factors = new float[m_scale_num];
        for (int i = 0; i < m_scale_num; i++) {

            m_scale_factors[i] = pow(m_scale_step,
                                     (float)((i-floor((m_scale_num-1) / 2))*(nScalesInterp/(1.0f*nScales)) ));
    }
    }
// interp scale factors which reduce calcuated amount
    interp_factor= new float[nScalesInterp];
    for (int i = 0; i < nScalesInterp; i++) {
        if (i <= floor(nScalesInterp/2)) {
            interp_factor[i] = pow(m_scale_step,
                                     i);
        } else {
            interp_factor[i] = pow(m_scale_step,
                                     i - nScalesInterp);
        }
    }

    cv::Mat y = cv::Mat(cv::Size(m_scale_num, 1),
                        CV_32F, cv::Scalar(0));

    float sigma = std::sqrt(m_scale_num)*m_scale_y_sigma;
    float scale_sigma=nScalesInterp*scale_sigma_factor;
    float mult = -0.5f/(sigma * sigma);
    mult= -0.5f/(scale_sigma*scale_sigma);
    // circleshift model_y;
    for (int i = 0; i < m_scale_num; i++)
    {
        int dist = (i-int(m_scale_num/2));
        float dist1;
        if(i<=floor(m_scale_num/2)) {
            dist1 = ((float)i ) * (float)nScalesInterp / (float)nScales;
        }
        else
        {
            dist1 = ((float)i - (float)m_scale_num) * (float)nScalesInterp / (float)nScales;
        }
        y.at<float>(0, i)=std::exp(mult*dist1*dist1);
    }
    m_scale_y = FFTTools::fftd(y);
    return 0;
}

int StapleTracker::add(Rect_T &roi, int cate_id) {
    Rect_T roi_s, winR;
    float model_area=roi.w*roi.h*filter_size*filter_size;
    float model_currentScaleFactor= sqrt(model_area/fix_model_size);
    roi_s.x = (int) (roi.x );
    roi_s.y = (int) (roi.y );
    roi_s.w= (int) (roi.w);
    roi_s.h =(int) (roi.h);
    int idx = getIdleIdx();
    if (-1 == idx)
        return 0;

    if (1 == isAlreadyIn(roi_s))
        return 0;

    float pad = (roi_s.w + roi_s.h) / 2.f;
    float bg_w, bg_h, fg_w, fg_h;
    bg_w = round(roi_s.w + pad);
    bg_h = round(roi_s.h + pad);
    float scale = sqrt(m_trans_fixed_area / (bg_w * bg_h));
    bg_w = (int(int(round(bg_w * scale)) / m_trans_cell_size) / 2 * 2 + 1) * m_trans_cell_size / scale;
    bg_h = (int(int(round(bg_h * scale)) / m_trans_cell_size) / 2 * 2 + 1) * m_trans_cell_size / scale;
    fg_w = roi_s.w - pad * m_trans_inner_padding;
    fg_h = roi_s.h - pad * m_trans_inner_padding;

    m_cfs[idx].pos[0] = roi_s.x + roi_s.w / 2;
    m_cfs[idx].pos[1] = roi_s.y + roi_s.h / 2;
    m_cfs[idx].target_size[0] = roi_s.w;
    m_cfs[idx].target_size[1] = roi_s.h;
    m_cfs[idx].base_tz[0] = roi_s.w;
    m_cfs[idx].base_tz[1] = roi_s.h;
    m_cfs[idx].bg_size[0] = bg_w;
    m_cfs[idx].bg_size[1] = bg_h;
    m_cfs[idx].fg_size[0] = fg_w;
    m_cfs[idx].fg_size[1] = fg_h;
    m_cfs[idx].scale = scale;
    frameW;
    float radius = MIN_T((bg_w - roi_s.w) * m_cfs[idx].scale,
                         (bg_h - roi_s.h) * m_cfs[idx].scale);
    m_cfs[idx].norm_delta_size = 2 * floor(radius / 2) + 1;

    m_objs[idx].roi = roi_s;
    m_objs[idx].cate_id = cate_id;
    m_objs[idx].obj_id = m_accumObjId++;
    m_objs[idx].status = 1;
    //Set the scale parameters
    float factor = 1;
    int norm_tw = round(roi_s.w * m_cfs[idx].scale);
    int norm_th = round(roi_s.h * m_cfs[idx].scale);
    if (norm_tw * norm_th > m_scale_max_area)
        factor = sqrt(m_scale_max_area / (norm_tw * norm_th));
    m_cfs[idx].scale_norm_tz[0] = floor(norm_tw * factor);
    m_cfs[idx].scale_norm_tz[1] = floor(norm_th * factor);
    m_cfs[idx].scale_adapt = 1;

    float t = sqrt(m_cfs[idx].target_size[0] * m_cfs[idx].target_size[1]);
    bool  resize_image1;
    int resize_scale1;
    if (t >= lowResize && t < upResize)
    {
        resize_image1 = 0;
        resize_scale1=1;
    } else if(sqrt(m_cfs[idx].target_size[0] * m_cfs[idx].target_size[1])>= upResize){
        resize_image1 = 0;
        resize_scale1 = 1;
    }
    else
    {
        resize_image1 =0;
        resize_scale1 = 1 ;
    }
    resize_image=resize_image1;
    resize_scale=resize_scale1;

    if(resize_image)
    {
        //cv::Mat img=cv::resize(m_img,1/resize_scale);
    }

    float search_area= m_cfs[idx].target_size[0]*m_cfs[idx].target_size[1]*filter_size*filter_size;

    float currentScaleFactor = sqrt(search_area/fix_model_size);
    float target_sz[2];
    target_sz[0]=m_cfs[idx].target_size[0]/currentScaleFactor;
    target_sz[1]=m_cfs[idx].target_size[1]/currentScaleFactor;

    float sz1[2];
    sz1[0]=sqrt(target_sz[0]*target_sz[1])*search_area_scale;
    sz1[1]=sqrt(target_sz[0]*target_sz[1])*search_area_scale;
    int sz[2];
    sz[0]=round(sz1[0]);
    sz[1]=round(sz1[1]);
    int tsz= round(target_sz[0]);
    sz[0]=sz[0] - (sz[0]-tsz)%2;
    tsz = round(target_sz[1]);
    sz[1]=sz[1] - (sz[1]-tsz)%2;

    m_cfs[idx].scale_max_factor = pow(m_scale_step,
                                      floor(log(MIN_T(m_img.width * 1.f /  target_sz[0], m_img.height * 1.f /  target_sz[1])) /
                                            log(m_scale_step)));

    m_cfs[idx].scale_min_factor = pow(m_scale_step,
                                      ceil(log(MAX_T(5.f / sz[0], 5.f / sz[1])) / log(m_scale_step)));
    m_cfs[idx].s_filt_sz[0] = floor(target_sz[0]*filter_size);
    m_cfs[idx].s_filt_sz[1] = floor(target_sz[1]*filter_size);
    m_cfs[idx].s_filt_sz[0] = floor(m_cfs[idx].s_filt_sz[0]/pfcell_size);
    m_cfs[idx].s_filt_sz[1] = floor(m_cfs[idx].s_filt_sz[1]/pfcell_size);

    int b_filt_sz[2];
    b_filt_sz[0] = floor(sz[0]/pfcell_size);
    b_filt_sz[1] = floor(sz[1]/pfcell_size);
    m_cfs[idx].b_filt_sz[0]=b_filt_sz[0];
    m_cfs[idx].b_filt_sz[1]=b_filt_sz[1];

    float output_sigma= sqrt(m_cfs[idx].s_filt_sz[0]*m_cfs[idx].s_filt_sz[1])*output_sigma_factor;
    m_cfs[idx].scale=1;
            int patch_sz[2];
    patch_sz[0]=floor(sz[0]*currentScaleFactor);
    patch_sz[1]=floor(sz[1]*currentScaleFactor);

    if(patch_sz[0]<1) patch_sz[0]=2;
    if(patch_sz[1]<1) patch_sz[1]=2;
    m_cfs[idx].currentScaleFactor=currentScaleFactor;

    m_cfs[idx].window_sz[0]=sz[0];
    m_cfs[idx].window_sz[1]=sz[1];
    m_cfs[idx].target_size1[0]=(float)(target_sz[0]);
    m_cfs[idx].target_size1[1]=(float)(target_sz[1]);

    m_cfs[idx].last_target_sz[0]=(float)(target_sz[0]);
    m_cfs[idx].last_target_sz[1]=(float)(target_sz[1]);
    m_cfs[idx].last_pos[0]=m_cfs[idx].pos[0] ;
    m_cfs[idx].last_pos[1]=m_cfs[idx].pos[1] ;
    pbase_target_sz[0]=(float)(target_sz[0]);
    pbase_target_sz[1]=(float)(target_sz[1]);
    pfeatures_sz[0]=b_filt_sz[0];
    pfeatures_sz[1]=b_filt_sz[1];
    psz[0]=sz[0];
    psz[1]=sz[1];
    float scale_sigma = nScalesInterp * scale_sigma_factor;
    m_scale_num=nScales;
    // above code initialize  params/

    // corresponding to initBACF.m in matlab edition of BACF
    cv::Mat roiImg;
    roiImg = getSubWin(idx);
    trainTransCF(idx, roiImg, 1.f, true,0.0125);//BACF_lr=0.0125
    trainScaleCF(idx, 1.f, true);
    m_curObjNum ++;
    return 0;
}

int StapleTracker::update()
{
    m_curObjNum = 0;
#ifdef __BENCHMARK
    float train_th=0.0f, update_th=0.0f;
#else
    float train_th=0.2f, update_th=0.15f;
#endif
    float confT, confS;
    for (int i=0 ;i<m_maxObjNum; i++)
    {
        if (-1 == m_objs[i].status)
            continue;

        detectTrans(i, confT);
        if (confT < update_th)
        {
            printf("trans conf %f\n", confT);
            m_objs[i].status = -1;
            continue;
        }

        detectScaleCF(i, confS);
        if (confS < update_th)
        {
            printf("scale conf %f\n", confS);
            m_objs[i].status = -1;
            continue;
        }

        if (m_cfs[i].pos[0]<0 ||
            m_cfs[i].pos[0]>m_img.width-1 ||
            m_cfs[i].pos[1]<0 ||
            m_cfs[i].pos[1]>m_img.height-1)
        {
            printf("img size: %d, %d pos: %f, %f\n",
                   m_img.width, m_img.height,
                   m_cfs[i].pos[0],
                   m_cfs[i].pos[1]);
            m_objs[i].status = -1;
            continue;
        }
        m_objs[i].roi = getShowRect(i);
        if (confT > train_th)
        {
            cv::Mat roiImg;
            roiImg = getSubWin(i);
            trainTransCF(i, roiImg, m_trans_lr_cf, false,0.0125);
            //trainTransPWP(i, roiImg, m_trans_lr_pwp,false);
            if (confS > train_th)
                trainScaleCF(i, m_scale_lr, false);
        }
        m_curObjNum ++;
    }
    return m_curObjNum;
}

Rect_T StapleTracker::getShowRect(int idx)
{
    Rect_T roi;
    roi.x = m_cfs[idx].pos[0]-m_cfs[idx].target_size[0]/2;
    roi.y = m_cfs[idx].pos[1]-m_cfs[idx].target_size[1]/2;
    roi.w = m_cfs[idx].target_size[0];
    roi.h = m_cfs[idx].target_size[1];

    roi.x = MAX_T(MIN_T(m_img.width, roi.x), 0);
    roi.y = MAX_T(MIN_T(m_img.height, roi.y), 0);
    roi.w = MIN_T(m_img.width-roi.x-1, roi.w);
    roi.h = MIN_T(m_img.height-roi.y-1, roi.h);
    return roi;
}
int StapleTracker::detectTrans(int idx, float &conf)
{

    cv::Mat roiImg = getSubWin(idx);

    cv::Mat resCF1 = detectTransCF(idx, roiImg);

    double maxcf,maxpwp;
    cv::minMaxLoc(resCF1,NULL,&maxcf,NULL,NULL);
   // cv::minMaxLoc(resPWP,NULL,&maxpwp,NULL,NULL);
    //cv::Mat resCF = 0.7f*maxpwp/maxcf*resCF1;//+0.3f*resPWP;
    cv::Mat resCF = resCF1;
    cv::Point2i pi,pipwp;
    double pv,pvpwp;
    cv::minMaxLoc(resCF, NULL, &pv, NULL, &pi);
   // cv::minMaxLoc(resPWP, NULL, &pvpwp, NULL, &pipwp);
    cv::Point2f pf((float)pi.x, (float)pi.y);

    int center = 0;
    pf.x=pi.x;
    pf.y = pi.y;

    if(pf.x>= ceil(resCF.cols/2)){pf.x= pf.x-resCF.cols;}
    if(pf.y>= ceil(resCF.rows/2)){pf.y= pf.y - resCF.rows;}
    float pos[2];
    // position update
    pos[0]=m_cfs[idx].pos[0] + (pf.x-center)*m_cfs[idx].currentScaleFactor*m_scale_cell_size;
    pos[1]=m_cfs[idx].pos[1] + (pf.y-center)*m_cfs[idx].currentScaleFactor*m_scale_cell_size;
    if(pos[0]-floor(m_cfs[idx].target_size[0]/2)>0&&pos[0]+floor(m_cfs[idx].target_size[0]/2)<m_img.width-1&&pos[1]-floor(m_cfs[idx].target_size[1]/2)>0&&pos[1]+floor(m_cfs[idx].target_size[1]/2)<m_img.height-1) {
        m_cfs[idx].pos[0] += (pf.x - center) * m_cfs[idx].currentScaleFactor * m_scale_cell_size;
        m_cfs[idx].pos[1] += (pf.y - center) * m_cfs[idx].currentScaleFactor * m_scale_cell_size;
    }

    conf = (float)pv;
    return 0;
}

cv::Mat StapleTracker::detectTransCF(int idx,
                                     cv::Mat &roiImg)
{
    cv::Mat feaCF;
    getTransFeaCF(roiImg, feaCF);
    applyHann2D(feaCF, m_cfs[idx].hann);
    cv::Mat  denzz = cv::Mat::zeros(feaCF.rows, feaCF.cols, CV_64F);
    std::vector<cv::Mat> feaSplit;
    std::vector<cv::Mat> feaSplit1;
    cv::split(feaCF, feaSplit);
    cv::split(m_cfs[idx].df, feaSplit1);
    int ch = feaSplit.size();
    int stepCh = ch*2;

    std::vector<cv::Mat> splidf;
    cv::split(m_cfs[idx].df,splidf);
    cv::Mat resF = cv::Mat::zeros(feaCF.rows, feaCF.cols,
                                  CV_32FC2);
    for (int c=0; c<ch; c++) {
        cv::Mat feaFFT = FFTTools::fftd(feaSplit[c]);
        float *pFea = (float *) (feaFFT.data);
            float *pRes = (float *) (resF.data);
            for (int h = 0; h < feaCF.rows; h++) {
                for (int w = 0; w < feaCF.cols; w++) {
                    pRes[0] += splidf[2*c].at<float>(h,w) * pFea[0] - splidf[2*c+1].at<float>(h,w) * pFea[1];
                    pRes[1] += splidf[2*c].at<float>(h,w) * pFea[1] + splidf[2*c+1].at<float>(h,w) * pFea[0];
                    pRes += 2;
                    pFea += 2;
                }
            }
    }
    return FFTTools::real(FFTTools::fftd(resF, true));
}

cv::Mat StapleTracker::detectTransPWP(int idx,
                                      cv::Mat &roiImg)
{
    cv::Mat res1 = cv::Mat::zeros(roiImg.rows,roiImg.cols,
                                 CV_32F);
    cv::Mat roiImg1= getSubWin(idx);
    cv::Mat res = cv::Mat::zeros(roiImg1.rows,roiImg1.cols,
                                  CV_32F);
    unsigned char *pImg = (unsigned char*)(roiImg1.data);
    float *pRes = (float *)(res.data);
    float *pBG = (float *)(m_cfs[idx].bgHist.data);
    float *pFG = (float *)(m_cfs[idx].fgHist.data);

    int range = 256/m_trans_color_bins;
    for (int h=0; h<res.rows; h++)
    {
        for (int w=0; w<res.cols; w++)
        {
            int idx1 = pImg[w*3]/range;
            int idx2 = pImg[w*3+1]/range;
            int idx3 = pImg[w*3+2]/range;
            int idx = idx3*m_trans_color_bins*m_trans_color_bins+idx2*m_trans_color_bins+idx1;
            float fg = pFG[idx];
            float bg = pBG[idx];
            if ((fg+bg)<0.000001f)
                *(pRes++) = 0;
            else
                *(pRes++) = fg/(fg+bg);
        }
        pImg += roiImg.step[0];
    }
    return res;
}

cv::Mat StapleTracker::getSubWinPWP(int idx)
{
    //Get the search region of PWP
    int norm_pwp_bg_w = round(m_cfs[idx].target_size[0]*m_cfs[idx].scale)+m_cfs[idx].norm_delta_size-1;
    int norm_pwp_bg_h = round(m_cfs[idx].target_size[1]*m_cfs[idx].scale)+m_cfs[idx].norm_delta_size-1;
    int pwp_bg_w = round(norm_pwp_bg_w/m_cfs[idx].scale);
    int pwp_bg_h = round(norm_pwp_bg_h/m_cfs[idx].scale);

    cv::Rect roi;

    float cx = m_cfs[idx].pos[0];
    float cy = m_cfs[idx].pos[1];
    roi.width = pwp_bg_w;
    roi.height = pwp_bg_h;
    roi.x = round(cx - roi.width/2.f);
    roi.y = round(cy - roi.height/2.f);

    // Get sub Image
    if(roi.width<=0||roi.height<=0)
    {
        printf("roi.width or roi.height no more than 0 No.2\n");
    }
    cv::Mat image = cv::Mat(m_img.height,m_img.width,
                            CV_8UC3, m_img.data[0]);
    cv::Mat z = RectTools::subwindow(image, roi,
                                     cv::BORDER_REPLICATE);
    if (z.cols!=norm_pwp_bg_w || z.rows!=norm_pwp_bg_h)
        cv::resize(z, z, cv::Size(norm_pwp_bg_w,
                                  norm_pwp_bg_h));
    return z;
}

cv::Mat StapleTracker::getSubWin(int idx)
{
    cv::Rect roi;
    float pos0=m_cfs[idx].pos[0];///m_cfs[idx].rate2img[0];   // code by zhangzheng in 2017/6/22
    float pos1=m_cfs[idx].pos[1];///m_cfs[idx].rate2img[1];

    float cx = m_cfs[idx].pos[0];
    float cy = m_cfs[idx].pos[1];

    roi.width = floor(m_cfs[idx].window_sz[0]*m_cfs[idx].currentScaleFactor);
    roi.height = floor(m_cfs[idx].window_sz[1]*m_cfs[idx].currentScaleFactor);

    roi.x = floor(cx - floor(roi.width/2.f)+1);
    roi.y = floor(cy - floor(roi.height/2.f)+1);
    int bg_w = round(m_cfs[idx].bg_size[0]*m_cfs[idx].scale);
    int bg_h = round(m_cfs[idx].bg_size[1]*m_cfs[idx].scale);

    cv::Mat image = cv::Mat(m_img.height,m_img.width,
                            CV_8UC3, m_img.data[0]);

    if(roi.width<=0||roi.height<=0)
    {
        printf("roi.width or roi.height no more than 0 No.3\n");
    }
    cv::Mat ztemp = RectTools::subwindow(image, roi,
                                         cv::BORDER_REPLICATE);

    int izz=0;
    cv::Mat z=cv::Mat::zeros(m_cfs[idx].window_sz[0], m_cfs[idx].window_sz[1],CV_32FC3);
    if ((ztemp.cols != m_cfs[idx].window_sz[0]) || (ztemp.rows !=  m_cfs[idx].window_sz[1])) {
        if(m_cfs[idx].window_sz[0]>ztemp.rows) {
            cv::resize(ztemp, ztemp, cv::Size(m_cfs[idx].window_sz[0], m_cfs[idx].window_sz[1]),0.0,0.0,1);
          }
        else //if(m_cfs[idx].window_sz[0]<z.rows)
        {
            cv::resize(ztemp, ztemp, cv::Size(m_cfs[idx].window_sz[0], m_cfs[idx].window_sz[1]),0.0,0.0,3);
        }
        return ztemp;
    }
    else{
       return ztemp;
    }
}

int StapleTracker::getTransFeaCF(cv::Mat &roiImg1,
                                 cv::Mat &feaHog1)
{
    // Extract HOG Feature
    cv::Mat roiGray;
    cv::cvtColor(roiImg1, roiGray, CV_BGR2GRAY);
    cv::Mat feaHog = fhog(roiGray, m_trans_cell_size,9);

    std::vector<cv::Mat> splitzz,splitzz1;
    std::vector<cv::Mat> res;
    cv::split(feaHog,splitzz);
    for(int i=0;i<feaHog.channels();i++)
    {
        res.push_back(splitzz[i]);
    }

    cv::Mat roiImg;
    if(feaHog.cols>roiImg1.rows)
    {
        cv::resize(roiImg1,roiImg,cv::Size(feaHog.rows,feaHog.cols),0,0,1);
    }
    else
    {
        cv::resize(roiImg1,roiImg,cv::Size(feaHog.rows,feaHog.cols),0,0,3);
    }
    cv::split(roiImg,splitzz1);
    for(int i=0;i<roiImg.channels();i++)
    {
        res.push_back(splitzz1[i]);
    }
    int p2323=0;
    int p1212=0;
    int a=0;

    cv::Mat feaHogzz(roiImg.rows,roiImg.cols,CV_32FC(res.size()),cv::Scalar(0));
    float * pa=(float *)(feaHogzz.data);
    float * p = (float *)(feaHog.data );
    std::vector<cv::Mat> spl;
    cv::split(roiImg,spl);

    // hog features is 31 channels while the BACF merge RGB feature to form 34 channles features.
    for(int i=0;i<roiImg.rows;i++)
    {
        for(int j=0;j <roiImg.cols;j++)
        {
            for(int n=0;n<roiImg.channels();n++)
            {
                pa[0]=(((float)roiImg.data[p2323])/255)-0.5;
                a++;
                // p1++;
                p2323++;
                pa++;
            }
            for(int j=0;j<feaHog.channels();j++)
            {
                pa[0]=p[0];
                p++;
                p1212++;
                pa++;
            }
        }
    }
    feaHog1=feaHogzz.clone();

#ifdef ENABLE_LAB_TRANS
                // Extract LAB Feature
    int cell_sizeQ = m_trans_cell_size*m_trans_cell_size;
    cv::Mat imgLab;
    cv::cvtColor(roiImg, imgLab, CV_BGR2Lab);
    unsigned char *pLAB = (unsigned char*)(imgLab.data);
    cv::Mat labFea = cv::Mat::zeros(cv::Size(feaHog.cols,
                                             feaHog.rows),
                                    CV_32FC(nLABCentroid));

    float *pFea = (float *)(labFea.data);
    for (int cY=0; cY<imgLab.rows; cY+=m_trans_cell_size)
    {
        for (int cX=0;cX<imgLab.cols;cX+=m_trans_cell_size)
        {
            for(int y=cY; y<cY+m_trans_cell_size; ++y)
            {
                for(int x=cX; x<cX+m_trans_cell_size; ++x)
                {
                    int idx = (imgLab.cols * y + x) * 3;
                    float l = (float)pLAB[idx];
                    float a = (float)pLAB[idx + 1];
                    float b = (float)pLAB[idx + 2];

                    // Iterate trough each centroid
                    float minDist = FLT_MAX;
                    int minIdx = 0;
                    for(int k=0; k<nLABCentroid; ++k)
                    {
                        float ld=(l-pLABCentroids[3*k]);
                        float ad=(a-pLABCentroids[3*k+1]);
                        float bd=(b-pLABCentroids[3*k+2]);

                        float dist =ld*ld + ad*ad + bd*bd;
                        if(dist < minDist){
                            minDist = dist;
                            minIdx = k;
                        }
                    }
                    pFea[minIdx] += 1.f/cell_sizeQ;
                }
            }
            pFea += nLABCentroid;
        }
    }
    std::vector<cv::Mat> fv0;
    std::vector<cv::Mat> fv1;
    cv::split(feaHog, fv0);
    cv::split(labFea, fv1);
    fv0.insert(fv0.end(), fv1.begin(), fv1.end());
    cv::merge(fv0, feaHog);
#endif
    return 0;

}


int StapleTracker:: trainTransCF(int idx, cv::Mat &roiImg,
                                float lr, bool isInit,double BACF_lr)
{
    cv::Mat feaCF;
    getTransFeaCF(roiImg, feaCF);//hog ,whole real data
    if (isInit)
    {
        m_cfs[idx].hann = createHann2D(feaCF.rows,
                                       feaCF.cols);
        m_cfs[idx].y=createGaussianPeak(idx, feaCF.rows,
                                        feaCF.cols,
                                        m_trans_y_sigma);
    }

    cv::Mat num, den, numVZZ, ZX,tempnumV,tempZX,tempZZ;//by zz in 2017/5/18
    applyHann2D(feaCF, m_cfs[idx].hann);

    int azz=0;
    float * feazz=(float *)(feaCF.data);
    std::vector<cv::Mat>feaCF1;
    cv::split(feaCF,feaCF1);
    cv::Mat feaFFT;
    solveTransCF(num, den, numVZZ, ZX, feaCF,  m_cfs[idx].y,feaFFT,idx);

    if (isInit)
    {
        //by zz in 2017/5/18 init for some new values, ZZ, ZX , X , df,sf, Ldsf.
        m_cfs[idx].ADMM_iteration=2;
        std::vector<cv::Mat> df1;
        std::vector<cv::Mat> sf1;
        std::vector<cv::Mat> Ldsf1;
        for(int n=0;n<feaFFT.channels();n++)
        {
            df1.push_back(cv::Mat::zeros(feaFFT.rows,feaFFT.cols,CV_32FC1));
            sf1.push_back(cv::Mat::zeros(feaFFT.rows,feaFFT.cols,CV_32FC1));
            Ldsf1.push_back(cv::Mat::zeros(feaFFT.rows,feaFFT.cols,CV_32FC1));
         }
        cv::merge(df1,m_cfs[idx].df);
        cv::merge(sf1,m_cfs[idx].sf);
        cv::merge(Ldsf1,m_cfs[idx].Ldsf);

        m_cfs[idx].num = num.clone();
        m_cfs[idx].den = den.clone();
        m_cfs[idx].alpha = num.clone();
        m_cfs[idx].X=feaFFT.clone();
        m_cfs[idx].numVZZ=numVZZ.clone();
        m_cfs[idx].ZX=ZX.clone();

    }
    else
    {
        m_cfs[idx].num = (1-lr)*m_cfs[idx].num + lr*num;
        m_cfs[idx].den = (1-lr)*m_cfs[idx].den + lr*den;
        m_cfs[idx].den = (1-lr)*m_cfs[idx].den + lr*den;
        m_cfs[idx].numVZZ = (1-lr)*m_cfs[idx].numVZZ + lr*numVZZ;
        m_cfs[idx].ZX = (1-lr)*m_cfs[idx].ZX + lr*ZX;
        m_cfs[idx].X = (1-lr)*m_cfs[idx].X + lr*feaFFT;
        //by zz in 2017/5/18
        tempZX=numVZZ.clone();
        tempZZ=ZX.clone();
        //by zz in 2017/5/18
    }

    //Compute the alpha
  float *pA1 = (float *)(m_cfs[idx].alpha.data);
  float *pNum1 = (float *)(m_cfs[idx].num.data);
   float *pDen1 = (float *)(m_cfs[idx].den.data);
   int channels =  feaCF.channels();
   for(int iz1=0; iz1<feaCF.rows;iz1++)
  {
       for (int w=0; w<feaCF.cols; w++)
        {
            float factor = 1.0f/(*(pDen1++)+m_trans_lambda);
            for (int c=0; c<channels; c++)
            {
                pA1[0]=pNum1[0]*factor;
              pA1[1]=pNum1[1]*factor;
                pA1 += 2;
                pNum1 += 2;
            }
        }
    }


    int minItr = 1;
    int maxItr=2;
    int term = 1;
    int visfilt=0;
    cv::Mat X=feaFFT;
    int  muo=1;
    int Mx0;
    int Mx1;
    Mx0=X.rows;
    Mx1=X.cols;
    int Nf=feaCF.channels();
    int Mf[2];
    Mf[0]=m_cfs[idx].bg_size[0];
    Mf[1]=m_cfs[idx].bg_size[1];
    ECF(idx,muo,X,Nf,term,minItr,maxItr,visfilt);  // main filter df learning precession

    return 0;
}

bool StapleTracker::ECF(int idx, int muo, cv::Mat &X,int Nf,int term,int minItr,int maxItr,int visfilt) {
    //  int MMx=Mx[0]*Mx[1];
    int MMx = X.rows * X.cols;
    int MMx0 = X.rows;
    int MMx1 = X.cols;
    int Nx = Nf;
    float lambda = m_trans_lambda;
    int i = 1;
    m_cfs[idx].ADMM_iteration=2;
    m_cfs[idx].muo=muo;
    m_cfs[idx].beta=10;
    m_cfs[idx].mumax=1000;
    while (i <= m_cfs[idx].ADMM_iteration) { //ADMM
        std::vector<cv::Mat> split1;
        cv::split(m_cfs[idx].df,split1);
        std::vector<cv::Mat> split2;
        cv::split(m_cfs[idx].numVZZ,split2);
        std::vector<cv::Mat> split3;
        cv::split(m_cfs[idx].ZX,split3);
        std::vector<cv::Mat> split4;
        cv::split(m_cfs[idx].Ldsf,split4);
        std::vector<cv::Mat> split5;
        cv::split(m_cfs[idx].sf,split5);

        // calculate sf (g in the paper): the same as  argmin_s.m in matlab implementation of BACF
        float * pasf = (float *)(m_cfs[idx].sf.data);
        float * pas2 = (float *)(m_cfs[idx].numVZZ.data);
        int izz1 = 0;
        for (int h = 0; h < X.rows; h++) {
            for (int w = 0; w < X.cols; w++) {
                for (int n = 0; n < Nf; n++) {

                    split2[2 * n].at<float>(h, w) = split2[2 * n].at<float>(h, w) + m_cfs[idx].muo;
                    split3[2 * n].at<float>(h, w) = split3[2 * n].at<float>(h, w) + float(m_cfs[idx].muo) * split1[2 * n].at<float>(h, w) - split4[2 * n].at<float>(h, w);
                    split3[2 * n + 1].at<float>(h, w) = split3[2 * n + 1].at<float>(h, w) + float(m_cfs[idx].muo) * split1[2 * n + 1].at<float>(h, w) - split4[2 * n + 1].at<float>(h, w);
                    pasf[0] = split3[2 * n].at<float>(h, w) / split2[2 * n].at<float>(h, w);
                    pasf[1] = split3[2 * n + 1].at<float>(h, w) / split2[2 * n].at<float>(h, w);
                    pasf += 2;
                }
            }
        }

        int N = Nf;
        int prodMx = MMx;
        int M2 = std::floor(m_cfs[idx].s_filt_sz[0] );
        int M1 = std::floor(m_cfs[idx].s_filt_sz[1] );

        //calculate df: (mu*sf+Lf)/(mu+(lambda/sqrt(prod(Md)))) in matlab implementation of BACF
        float at = (m_cfs[idx].muo + lambda / (std::sqrt(prodMx)));
        izz1=0;
        float * pzzzdf= (float *)(m_cfs[idx].df.data);
       /* for (int h = 0; h < X.rows; h++) {
            for (int w = 0; w < X.cols; w++) {
                for (int n = 0; n < Nf; n++)
                {
                    split1[2*n].at<float>(h,w) = (split5[2*n].at<float>(h,w) * m_cfs[idx].muo + split4[2*n].at<float>(h,w) ) / at;
                    split1[2*n+1].at<float>(h,w)  = ( split5[1+2*n].at<float>(h,w)  * m_cfs[idx].muo + split4[2*n+1].at<float>(h,w)) / at;
                    pzzzdf[0]=split1[2 * n].at<float>(h, w);
                    pzzzdf[1]=split1[2 * n + 1].at<float>(h, w);
                    pzzzdf+=2;
                    izz1+=2;
                }
            }
        }*/
        m_cfs[idx].df=(m_cfs[idx].sf*m_cfs[idx].muo+m_cfs[idx].Ldsf)/at;
        // ifftvec.m: ifft2(x2)  //by zhangzheng in 2017 6.16
        cv::Mat xtem;
        std::vector<cv::Mat> xtempzz;
        std::vector<cv::Mat> dftemp;
        float *pdfz2=(float *)(m_cfs[idx].df.data);
        cv::split(m_cfs[idx].df,dftemp);
        for (int c=0; c<Nf; c++)
        {
            cv::Mat resF= cv::Mat::zeros(m_cfs[idx].df.rows, m_cfs[idx].df.cols,
                                         CV_32FC2);
            float *pRes3 = (float *)(resF.data);
            for (int h=0; h<m_cfs[idx].df.rows; h++)
            {
                for (int w=0; w<m_cfs[idx].df.cols; w++) {
                    pRes3[0]=dftemp[2*c].at<float>(h,w);
                    pRes3[1]=dftemp[2*c+1].at<float>(h,w);
                    pRes3 += 2;
                    pdfz2 += 2;
                }
            }
            xtempzz.push_back(FFTTools::fftd(resF, true));
        }
        cv::merge(xtempzz,xtem);
        std::vector<cv::Mat> splitx;
        cv::split(xtem,splitx);

        // ifftvec.m:  cropping:
        cv::Mat azztemp = m_cfs[idx].df.clone(); //A.t() transpose of A
        int delta[2];
        delta[0]=(M1-m_cfs[idx].df.rows>0)?std::floor((M1-m_cfs[idx].df.rows)/2):std::ceil((M1-m_cfs[idx].df.rows)/2);
        delta[1]=(M2-m_cfs[idx].df.cols>0)?std::floor((M2-m_cfs[idx].df.cols)/2):std::ceil((M2-m_cfs[idx].df.cols)/2);
        cv::Mat matzz=cv::Mat(azztemp.cols,azztemp.rows,CV_32FC(azztemp.channels()));
        matzz=xtem;
        std::vector<cv::Mat> splazz;
        cv::split(matzz,splazz);
        float * pmat = (float *)(matzz.data);

        cv::Mat matemp=cv::Mat(azztemp.rows,azztemp.cols, CV_32FC(azztemp.channels()));
        cv::Point2f P2f(delta[1],delta[0]);
        std::vector<cv::Mat> splpma1;
        cv::split(matzz,splpma1);
        std::vector<cv::Mat> splmatzz;
        cv::split(matemp,splmatzz);
        for(int i=0 ; i< splpma1.size();i ++)
        {
            shift(splpma1[i], splmatzz[i], P2f);
        }
        cv::merge(splmatzz,matemp);
        float * pmatemp= (float *)(matemp.data);
        std::vector<cv::Mat> splma,splitazz;
        cv::split(matemp,splma);
        cv::split(matemp,splitazz);
        cv::Mat r=cv::Mat(M1,M2,CV_32FC(matemp.channels()));
        r=matemp(cv::Rect(0,0,M2,M1));

        // fftvec.m: padding
        std::vector<cv::Mat> splitma;
        cv::split(m_cfs[idx].df,splitma);
        std::vector<cv::Mat> splir;
        cv::split(r,splir);
        cv::Mat zztemp;
        cv::merge(splitma,zztemp);
        int ch = Nf;
        std::vector<cv::Mat> xtemp1;
        cv::Mat xtemp;
        cv::Mat padded=cv::Mat::zeros(m_cfs[idx].df.rows,m_cfs[idx].df.cols,CV_32FC(m_cfs[idx].df.channels()));
        cv::split(padded,xtemp1);
        cv::Rect roi(-delta[1],-delta[0],r.cols,r.rows);
        for(int i=0; i< r.channels();i++) {
            splir[i].copyTo(xtemp1[i](roi));
        }
        cv::merge(xtemp1,padded);
        std::vector<cv::Mat> splipa;
        // fftvec: fft2(xpad)
        std::vector<cv::Mat> xtem3;
        float *pdf2=(float *)(m_cfs[idx].df.data);
        for (int c=0; c<ch; c++)
        {
            cv::Mat resF= cv::Mat::zeros(m_cfs[idx].df.rows, m_cfs[idx].df.cols,
                                    CV_32FC2);
                float *pRes3 = (float *)(resF.data);
                for (int h=0; h<m_cfs[idx].df.rows; h++)
                {
                    for (int w=0; w<m_cfs[idx].df.cols; w++) {
                        pRes3[0]=xtemp1[2*c].at<float>(h,w);
                        pRes3[1]=xtemp1[2*c+1].at<float>(h,w);
                        pdf2[0]=xtemp1[2*c].at<float>(h,w);
                        pdf2[1]=xtemp1[2*c+1].at<float>(h,w);
                        pRes3 += 2;
                        pdf2 += 2;
                    }
                }
            xtem3.push_back(FFTTools::fftd(resF));
        }
        cv::merge(xtem3,xtemp);
        xtemp1.clear();
        xtem3.clear();
        m_cfs[idx].df=xtemp.clone(); // get df

        //get Ldsf
        m_cfs[idx].Ldsf=m_cfs[idx].Ldsf+(m_cfs[idx].sf-m_cfs[idx].df)* m_cfs[idx].muo;
        // get mu
            if((m_cfs[idx].muo * m_cfs[idx].beta) > m_cfs[idx].mumax) {
                m_cfs[idx].muo = m_cfs[idx].mumax;
            } else {
                m_cfs[idx].muo = m_cfs[idx].muo * m_cfs[idx].beta;
            }
            //float *DFZ = (float *)(dfo.data);
    i = i + 1;
    }
    return true;
}


int StapleTracker::trainTransPWP(int idx,
                                 cv::Mat &roiImg,
                                 float lr, bool isInit)
{
    cv::Mat histBg = cv::Mat::zeros(1, m_trans_color_bins*m_trans_color_bins*m_trans_color_bins, CV_32F);
    cv::Mat histFg = cv::Mat::zeros(1, m_trans_color_bins*m_trans_color_bins*m_trans_color_bins, CV_32F);

    int bg_h = roiImg.rows, bg_w = roiImg.cols;
    //Get PWP Histgram
    cv::Mat bgMask=cv::Mat::ones(bg_h, bg_w, CV_8U);
    int offsetX = (m_cfs[idx].bg_size[0]-m_cfs[idx].target_size[0])*m_cfs[idx].scale/2;
    int offsetY = (m_cfs[idx].bg_size[1]-m_cfs[idx].target_size[1])*m_cfs[idx].scale/2;
    if(offsetX<=0){ offsetX=0;}
    if(offsetY<=0){offsetY=0;}
    if(offsetY>=0.5*bg_w){offsetY=0.5*bg_w;}
    if(offsetX>=0.5*bg_h){offsetX=0.5*bg_h;}
    bgMask(cv::Rect(offsetX, offsetY,
                    bg_w-2*offsetX, bg_h-2*offsetY))=cv::Scalar::all(0.0);

    cv::Mat fgMask=cv::Mat::zeros(bg_h, bg_w, CV_8U);
    offsetX = (m_cfs[idx].bg_size[0]-m_cfs[idx].fg_size[0])*m_cfs[idx].scale/2;
    offsetY = (m_cfs[idx].bg_size[1]-m_cfs[idx].fg_size[1])*m_cfs[idx].scale/2;
    if(offsetX<=0){ offsetX=0;}
    if(offsetY<=0){offsetY=0;}
    if(offsetY>=0.5*bg_w){offsetY=0.5*bg_w;}
    if(offsetX>=0.5*bg_h){offsetX=0.5*bg_h;}
    fgMask(cv::Rect(offsetX, offsetY,
                    bg_w-2*offsetX, bg_h-2*offsetY))=1;

    unsigned char *pBGMask=(unsigned char *)(bgMask.data);
    unsigned char *pFGMask=(unsigned char *)(fgMask.data);
    unsigned char *pImg=(unsigned char *)(roiImg.data);
    float *pBG=(float *)(histBg.data);
    float *pFG=(float *)(histFg.data);
    int range = 256/m_trans_color_bins;
    for (int h=0; h<bg_h; h++)
    {
        for (int w=0; w<bg_w; w++)
        {
            int idx1 = pImg[w*3]/range;
            int idx2 = pImg[w*3+1]/range;
            int idx3 = pImg[w*3+2]/range;
            int idx = idx3*m_trans_color_bins*m_trans_color_bins+idx2*m_trans_color_bins+idx1;
            pBG[idx] += *(pBGMask++);
            pFG[idx] += *(pFGMask++);
        }
        pImg += roiImg.step[0];

    }
    histBg = histBg/(cv::sum(histBg)[0]);
    histFg = histFg/(cv::sum(histFg)[0]);
    if (isInit)
    {
        m_cfs[idx].bgHist = histBg.clone();
        m_cfs[idx].fgHist = histFg.clone();
    }
    else
    {
        m_cfs[idx].bgHist=(1-lr)*m_cfs[idx].bgHist+lr*histBg;
        m_cfs[idx].fgHist=(1-lr)*m_cfs[idx].fgHist+lr*histFg;
    }
    return 0;
}

int StapleTracker::solveTransCF(cv::Mat &num, cv::Mat &den, cv::Mat &numVZZ, cv::Mat &ZX,
                                cv::Mat &fea, cv::Mat y, cv::Mat &feaFFT,int idx)
{
    float norm = 1.0f/(fea.cols*fea.rows);
    std::vector<cv::Mat> feaSplit;
    cv::split(fea, feaSplit);
    cv::Mat temp;

    std::vector<cv::Mat> numV;
    std::vector<cv::Mat> numVZZ21;
    std::vector<cv::Mat> numVZX;
    std::vector<cv::Mat> feaFFT1;

    den = cv::Mat::zeros(fea.rows, fea.cols, CV_32F);
    // debug
    for (int i=0; i<feaSplit.size(); i++)
    {
        cv::Mat feaFFT2 = FFTTools::fftd(feaSplit[i]);

        feaFFT1.push_back(FFTTools::fftd(feaSplit[i]));
        numV.push_back(FFTTools::complexConjMult(y,feaFFT2));


        numVZZ21.push_back(FFTTools::complexConjMult(feaFFT2,
                                                     feaFFT2));
        den=den + FFTTools::complexSelfConjMult(feaFFT2)*norm;
    }
    cv::merge(numV, num);
    cv::merge(numVZZ21, numVZZ);
    cv::merge(feaFFT1,feaFFT);
    feaSplit.clear();
    numV.clear();
    numVZZ21.clear();
    numVZX.clear();
    feaFFT1.clear();
    ZX=num;

    return 0;
}
int StapleTracker::applyHann2D(cv::Mat& fea, cv::Mat& hann)
{
    int channels = fea.channels();
    int width = fea.cols;
    int height = fea.rows;
    float *pFea = (float *)(fea.data);
    float *pHann = (float *)(hann.data);
    for (int h=0; h<height; h++)
    {
        for (int w=0; w<width; w++)
        {
            float factor = *(pHann++);
            for (int c=0; c<channels; c++) {
                pFea[c] = pFea[c] * factor;
              }
            pFea += channels;
        }
    }
    return 0;
}

cv::Mat StapleTracker::createHann2D(int height,
                                    int width)
{
    cv::Mat hann1t = cv::Mat(cv::Size(width,1), CV_32F, cv::Scalar(0));
    cv::Mat hann2t = cv::Mat(cv::Size(1,height), CV_32F, cv::Scalar(0));

    for (int i = 0; i < hann1t.cols; i++)
        hann1t.at<float > (0, i) = 0.5 * (1 - std::cos(2 * PI_T * i / (hann1t.cols - 1)));
    for (int i = 0; i < hann2t.rows; i++)
        hann2t.at<float > (i, 0) = 0.5 * (1 - std::cos(2 * PI_T * i / (hann2t.rows - 1)));

    return hann2t*hann1t;

}

cv::Mat StapleTracker::createGaussianPeak(int idx, int sizey, int sizex, float sigma)
{
    cv::Mat_<float> res(sizey, sizex);


    int syh = (sizey) / 2;
    int sxh = (sizex) / 2;
    //printf("m_target_size 66666: %f , %f\n",m_cfs[idx].target_size[0],m_cfs[idx].target_size[1]);
    int tw = round(m_cfs[idx].target_size[0]);//*m_cfs[idx].scale);
    int th = round(m_cfs[idx].target_size[1]);//*m_cfs[idx].scale);
    //float output_sigma = std::sqrt((float)th*tw)*sigma/m_trans_cell_size;
    float output_sigma = std::sqrt((float)m_cfs[idx].s_filt_sz[0]*m_cfs[idx].s_filt_sz[1])*output_sigma_factor;
    float mult = -0.5 / (output_sigma * output_sigma);

    for (int i = 0; i < sizey; i++)
        for (int j = 0; j < sizex; j++)
        {
            int ih = i - syh;
            int jh = j - sxh;
            res.at<float>(i, j) = std::exp(mult* (ih*ih+jh*jh));
        }

    // circshift the data
    cv::Mat resCS = res.clone();
    for (int h=0; h<sizey; h++)
    {
        int idy = (h+syh)%sizey;
        for (int w=0; w<sizex; w++)
        {
            int idx = (w+sxh)%(sizex);
            resCS.at<float>(idy,idx) = res.at<float>(h,w);
        }
    }
    cv::Mat resF = FFTTools::fftd(resCS);
    return resF;
}

cv::Mat StapleTracker::cropTransResponseCF(int idx,
                                           cv::Mat &res)
{
    int size = floor((res.cols + res.rows)/2 );
    int sizeH = int(size/2);
    cv::Mat newRes = cv::Mat(cv::Size(size, size),
                             CV_32F);
    float *pDst = (float *)(newRes.data);
    for (int h=0; h<newRes.rows; h++)
    {
        int idy = (res.rows+h-sizeH-1)%res.rows;
        for (int w=0; w<newRes.cols; w++)
        {
            int idx = (res.cols+w-sizeH-1)%res.cols;
            *(pDst++) = res.at<float>(idy, idx);
        }
    }

    return newRes;
}

cv::Mat StapleTracker::cropTransResponsePWP(int idx,
                                            cv::Mat &res)
{

    cv::Mat II = cv::Mat(cv::Size(res.cols+1,
                                  res.rows+1),
                         CV_32F);

    integral<float, float>((float *)(res.data),
                           res.cols, res.cols, res.rows,
                           (float *)(II.data));
    int tw = round(m_cfs[idx].target_size1[0]);
    int th = round(m_cfs[idx].target_size1[1]);
    float factor = 1.f/(tw*th);

    cv::Mat newRes=cv::Mat(cv::Size(0-(tw-res.rows)+1,0-(th-res.cols)+1),CV_32F);
    float *pDst = (float *)(newRes.data);
    for (int h=0; h<newRes.rows; h++)
    {
        for (int w=0; w<newRes.cols; w++)
        {
            float val = II.at<float>(h,w)+            \
                II.at<float>(h+th,w+tw)-              \
                II.at<float>(h,w+tw)   -              \
                II.at<float>(h+th,w);
            *(pDst++) = val*factor;
        }
    }
    return newRes;
}

int StapleTracker::detectScaleCF(int idx, float &conf)
{
    cv::Mat fea = getScaleFeaCF(idx);
    fea = FFTTools::fftd1d(fea, 1);

    cv::Mat resH = cv::Mat(cv::Size(m_scale_num, 1),
                           CV_32FC2);
    float *ptr1, *ptr2, *ptr3, *ptr4;
    ptr1 = (float *)(resH.data);
    ptr2 = (float *)(fea.data);
    ptr3 = (float *)(m_cfs[idx].num_scale.data);
    ptr4 = (float *)(m_cfs[idx].den_scale.data);

    for (int i=0; i<m_scale_num; i++)
    {
        double realV=0, complV=0;
        for (int j=0; j<fea.cols; j++)
        {
            realV += ptr2[0]*ptr3[0] - ptr2[1]*ptr3[1];
            complV+= ptr2[0]*ptr3[1] + ptr2[1]*ptr3[0];
            ptr2 += 2;
            ptr3 += 2;
        }
        ptr1[0] = realV / (ptr4[i] + m_scale_lambda);
        ptr1[1] = complV/ (ptr4[i] + m_scale_lambda);
        ptr1 += 2;
    }
    int len=nScales;
    int minsz=std::min(len,nScalesInterp);

    float scaling=(float)nScalesInterp/len;
    int newSize=nScalesInterp;
    cv::Mat resH1 = cv::Mat(cv::Size(nScalesInterp, 1),
                           CV_32FC2);
    float * presH1=(float *)(resH1.data);
    float * presH= (float *)(resH.data);
    std::vector<cv::Mat> splitres;
    cv::split(resH,splitres);
    int mids=ceil(minsz/2)+1;
    int mide=floor((minsz-1)/2)-1;

    presH1=(float *)(resH1.data);

    presH= (float *)(resH.data);
    int a=0;
    // 17 scale response interp to 33 scales.
    for(int i=0;i<nScalesInterp;i++)
    {
        if(i<mids)
        {
            presH1[0]=scaling*presH[0];
            presH++;
            presH1[1]=scaling*presH[0];
            presH++;
        }
        else if((i<nScalesInterp-mide-1)&&(i>=mids))
        {
            presH1[0]=0;
            presH1[1]=0;
        }
        else //if(i>=nScalesInterp-mide-1)
        {
            presH1[0]=scaling*presH[0];
            presH++;
            presH1[1]=scaling*presH[0];
            presH++;
        }
        presH1+=2;
    }

    //cv::Mat res1 = FFTTools::real(FFTTools::fftd(resH,
    //                                            true));
    cv::Mat res = FFTTools::real(FFTTools::fftd(resH1,
                                                true));

    cv::Point2i pi,pit;
    double pv,pvt;
    //cv::minMaxLoc(res1, NULL, &pv, NULL, &pi);
    cv::minMaxLoc(res, NULL, &pvt, NULL, &pit);
    conf = (float)pvt;

   // float bestScale = m_scale_factors[pi.x];
    float bestScale=interp_factor[pit.x];
    m_cfs[idx].scale_adapt *= bestScale;

    if (m_cfs[idx].scale_adapt < m_cfs[idx].scale_min_factor)
        m_cfs[idx].scale_adapt=m_cfs[idx].scale_min_factor;
    else if(m_cfs[idx].scale_adapt > m_cfs[idx].scale_max_factor)
        m_cfs[idx].scale_adapt= m_cfs[idx].scale_max_factor;

    float tz_w, tz_h, bg_w, bg_h;//, fg_w, fg_h;
    tz_w = m_cfs[idx].base_tz[0]*m_cfs[idx].scale_adapt;
    tz_h = m_cfs[idx].base_tz[1]*m_cfs[idx].scale_adapt;
    float pad = (tz_w+tz_h)/2.f;

    bg_w = tz_w+pad;
    bg_h = tz_h+pad;
    float scale = sqrt(m_trans_fixed_area/(bg_w*bg_h));

    //printf("m_target_size 77777: %d , %d\n",m_cfs[idx].target_size[0],m_cfs[idx].target_size[1]);
    int  bg_area=round(sqrt(m_cfs[idx].target_size[0]*m_cfs[idx].target_size[1])*search_area_scale);
    int pre_norm_bg_w=round(m_cfs[idx].bg_size[0]*m_cfs[idx].scale);
    int pre_norm_bg_h=round(m_cfs[idx].bg_size[1]*m_cfs[idx].scale);
    bg_w = pre_norm_bg_w/scale;
    bg_h = pre_norm_bg_h/scale;
    int fg_w = round(tz_w-pad*m_trans_inner_padding);
    int fg_h = round(tz_h-pad*m_trans_inner_padding);

    m_cfs[idx].target_size[0] = floor(m_cfs[idx].target_size1[0]*m_cfs[idx].currentScaleFactor);
    m_cfs[idx].target_size[1] = floor(m_cfs[idx].target_size1[1]*m_cfs[idx].currentScaleFactor);
    m_cfs[idx].bg_size[0] = bg_area-round((bg_area-m_cfs[idx].base_tz[0])%2);
    m_cfs[idx].bg_size[1] = bg_area-(bg_area-m_cfs[idx].base_tz[1])%2;
    m_cfs[idx].fg_size[0] = (fg_w+((m_cfs[idx].bg_size[0])-fg_w)%2);
    m_cfs[idx].fg_size[1] = fg_h+((m_cfs[idx].bg_size[1])-fg_w)%2;
    m_cfs[idx].currentScaleFactor = bestScale* m_cfs[idx].currentScaleFactor;

    return 0;
}

int StapleTracker::getOneScaleFeaCF(cv::Mat &roiImg,
                                    cv::Mat &feaHog)
{
    //Get HOG Feature
    cv::Mat roiGray;
    cv::cvtColor(roiImg, roiGray, CV_BGR2GRAY);
    feaHog = fhog(roiGray, m_scale_cell_size);


#ifdef ENABLE_LAB_SCALE
    // Extract LAB Feature
    int cell_sizeQ = m_scale_cell_size*m_scale_cell_size;
    cv::Mat imgLab;
    cv::cvtColor(roiImg, imgLab, CV_BGR2Lab);
    unsigned char *pLAB = (unsigned char*)(imgLab.data);
    cv::Mat feaLab = cv::Mat::zeros(cv::Size(feaHog.cols,
                                             feaHog.rows),
                                    CV_32FC(nLABCentroid));

    float *pFea = (float *)(feaLab.data);
    for (int cY=0; cY<imgLab.rows-m_scale_cell_size; cY+=m_scale_cell_size)
    {
        for (int cX=0;cX<imgLab.cols-m_scale_cell_size;cX+=m_scale_cell_size)
        {
            for(int y=cY; y<cY+m_scale_cell_size; ++y)
            {
                for(int x=cX; x<cX+m_scale_cell_size; ++x)
                {
                    int idx = (imgLab.cols * y + x) * 3;
                    float l = (float)pLAB[idx];
                    float a = (float)pLAB[idx + 1];
                    float b = (float)pLAB[idx + 2];

                    // Iterate trough each centroid
                    float minDist = FLT_MAX;
                    int minIdx = 0;
                    for(int k=0; k<nLABCentroid; ++k)
                    {
                        float ld=(l-pLABCentroids[3*k]);
                        float ad=(a-pLABCentroids[3*k+1]);
                        float bd=(b-pLABCentroids[3*k+2]);

                        float dist =ld*ld + ad*ad + bd*bd;
                        if(dist < minDist){
                            minDist = dist;
                            minIdx = k;
                        }
                    }
                    pFea[minIdx] += 1.f/cell_sizeQ;
                }
            }
            pFea += nLABCentroid;
        }
    }
    std::vector<cv::Mat> fv0;
    std::vector<cv::Mat> fv1;
    cv::split(feaHog, fv0);
    cv::split(feaLab, fv1);
    fv0.insert(fv0.end(), fv1.begin(), fv1.end());
    cv::merge(fv0, feaHog);
#endif
    return 0;
}

cv::Mat StapleTracker::getScaleFeaCF(int idx) {
    int tw = m_cfs[idx].target_size1[0]*m_cfs[idx].currentScaleFactor;
    int th = m_cfs[idx].target_size1[1]*m_cfs[idx].currentScaleFactor;

    cv::Mat image = cv::Mat(m_img.height, m_img.width,
                            CV_8UC3, m_img.data[0]);
    cv::Mat feaTmp, feas;
    for (int i=0; i<m_scale_num; i++)
    {
        int scale_tw = MAX_T(floor(tw*m_scale_factors[i]),m_scale_cell_size);
        int scale_th = MAX_T(floor(th*m_scale_factors[i]),m_scale_cell_size);
        cv::Rect roi;
        roi.width = scale_tw;
        roi.height = scale_th;
        roi.x = round(m_cfs[idx].pos[0] - roi.width/2.f);
        roi.y = round(m_cfs[idx].pos[1] - roi.height/2.f);
        if(roi.width<=0||roi.height<=0)
        {
            printf("roi.width or roi.height no more than 0 No.1 \n");
        }

         cv::Mat z;
         z = RectTools::subwindow(image, roi,
                                 cv::BORDER_REPLICATE);
         if (z.cols != m_cfs[idx].scale_norm_tz[0] || z.rows != m_cfs[idx].scale_norm_tz[1]) {

         if(z.cols>m_cfs[idx].scale_norm_tz[1]){   cv::resize(z, z, cv::Size(m_cfs[idx].scale_norm_tz[0], m_cfs[idx].scale_norm_tz[1]),0,0,1);}
            else{

             cv::resize(z, z, cv::Size(m_cfs[idx].scale_norm_tz[0],
                                       m_cfs[idx].scale_norm_tz[1]),0,0,3);
          }
        }

        getOneScaleFeaCF(z, feaTmp);
        feaTmp = feaTmp.reshape(1,1);
        feaTmp = feaTmp * m_scale_hann[i];
        if (0==i)
            feas = feaTmp.clone();
        else
            feas.push_back(feaTmp);
    }
    return feas;
}

int StapleTracker::trainScaleCF(int idx, float lr, bool isInit)
{
    cv::Mat fea = getScaleFeaCF(idx);
    fea = FFTTools::fftd1d(fea, 1);

    cv::Mat num = cv::Mat(cv::Size(fea.cols, fea.rows),
                          CV_32FC2);
    cv::Mat den = cv::Mat(cv::Size(m_scale_num, 1),
                          CV_32F);

    float *ptr1, *ptr2, *ptr3;
    ptr1 = (float *)(den.data);
    ptr2 = (float *)(fea.data);
    for (int i=0; i<m_scale_num; i++)
    {
        double val = 0;
        for (int j=0; j<fea.cols; j++)
        {
            val += ptr2[0]*ptr2[0] + ptr2[1]*ptr2[1];
            ptr2 +=2;
        }
        ptr1[i] = val;
    }

    ptr1 = (float *)(num.data);
    ptr2 = (float *)(fea.data);
    ptr3 = (float *)(m_scale_y.data);
    for (int i=0; i<num.rows; i++)
    {
        for (int j=0; j<num.cols; j++)
        {
            ptr1[0] = ptr2[0]*ptr3[0] + ptr2[1]*ptr3[1];
            ptr1[1] = ptr2[0]*ptr3[1] - ptr2[1]*ptr3[0];
            ptr1 += 2;
            ptr2 += 2;
        }
        ptr3 += 2;
    }

    if (isInit)
    {
        m_cfs[idx].num_scale = num.clone();
        m_cfs[idx].den_scale = den.clone();
    }
    else
    {
        m_cfs[idx].num_scale=(1-lr)*m_cfs[idx].num_scale + lr*num;
        m_cfs[idx].den_scale=(1-lr)*m_cfs[idx].den_scale + lr*den;
    }
    return 0;
}
