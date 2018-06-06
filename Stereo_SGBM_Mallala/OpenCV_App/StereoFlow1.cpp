//#include "stdafx.h"
#include <stdio.h>
#include "cv.h"
#include "highgui.h"
#include <math.h>
#include "Utils.h"
#include <vector>
#include <algorithm>
#include "cvaux.h"


static const double pi = 3.14159265358979323846;

static void VecSmooth(std::vector<float> & vecSm, const std::vector<float> &vec, int Len =0)
{
	// /*
	if(Len <=0) Len = 20;
	float * filter = new float[Len];
	float Tlen =0;
	for(int i = 0; i < Len; i++) {
		filter[i] = sin(pi*(i+0.5f)/Len);
		Tlen += filter[i];
	}
	int last = (int) vec.size()-1;

	float avg, sum;
	for(int j = 0; j<= last; j++)
	{
		sum = 0;
		for(int i =0; i<Len; i++)
		{
			int pos = j+i-Len/2;
			if(pos <0) pos = -pos;
			//if(pos>last) pos = 2*last -pos;
			if(pos>last) {sum = avg*Tlen; break;}
			sum += filter[i]*vec[pos];
		}
			
		avg= sum/Tlen;
		
		vecSm.push_back(avg);;
	}
	delete filter;
}


static CvPoint2D32f GetMedian(CvPoint2D32f *f1, CvPoint2D32f *f2, int numFeatures, int height, int width)
{
	std::vector<float> vx, vy;
	for(int i =0; i< numFeatures; i++)
	{
		if (f1[i].x < width*0.1 || f1[i].x > width*0.9) continue;
		if (f2[i].x < width*0.1 || f2[i].x > width*0.9) continue;
		if (f1[i].y < height*0.1 || f1[i].y > height*0.9) continue;
		if (f2[i].y < height*0.1 || f2[i].y > height*0.9) continue;
		vx.push_back(f2[i].x-f1[i].x); 
		vy.push_back(f2[i].y-f1[i].y);
	}
	sort( vx.begin( ), vx.end( ) );
	sort( vy.begin( ), vy.end( ) );
	size_t mid = vx.size()/2;
	CvPoint2D32f p; 
	p.x =vx[mid]; p.y =vy[mid];
	return p;
}

static void RemoveBoundary(CvPoint2D32f *f1, CvPoint2D32f *f2, int & numFeatures, int height, int width, char * found = NULL)
{
	std::vector<int> vi;
	int boundary = min(11, 0.1*min(height, width));
	for(int i =0; i< numFeatures; i++)
	{
		if (found!=NULL && found[i]==0) continue;
		if (f1[i].x < boundary || f1[i].x > width-boundary) continue;
		if (f2[i].x < boundary || f2[i].x > width-boundary) continue;
		if (f1[i].y < boundary || f1[i].y > height-boundary) continue;
		if (f2[i].y < boundary || f2[i].y > height-boundary) continue;
		vi.push_back(i); 
	}

	numFeatures = vi.size();
	for(int i =0; i< numFeatures; i++)
	{
		f1[i] = f1[vi[i]];
		f2[i] = f2[vi[i]];
	}
}


static float GetMoveRotaion(CvPoint2D32f *f1, CvPoint2D32f *f2, int numFeatures, int height, int width, CvPoint2D32f & m, char * found = NULL)
{
	// /*
	CvPoint2D32f * F1 = new CvPoint2D32f[numFeatures];
	CvPoint2D32f * F2 = new CvPoint2D32f[numFeatures];
	for(int i = 0; i< numFeatures; i++)
	{
		F1[i].x = f1[i].x - width/2;
		F1[i].y = f1[i].y - height/2;
		F2[i].x = f2[i].x - width/2;
		F2[i].y = f2[i].y - height/2;
	}

	CvMat* H;
	double error_lim = 3.0;
	do
	{
		H = ransac_xform(F1, F2, numFeatures, lsq_homog, 4, 0.01,
			homog_xfer_err, error_lim);
		error_lim *= 2;
	} while (!H);
		
	float C,S, Ang;
	C = H->data.fl[0] + H->data.fl[4];
	S = H->data.fl[3] - H->data.fl[1];
	Ang = atan2(S/C,1);

	m.x = H->data.fl[2];
	m.y = H->data.fl[5];

	delete F1;
	delete F2;
	cvReleaseMat(&H);
	return Ang;
	// */

	 /*
	std::vector<float> vx, vy, smvx, smvy;
	int boundary = min(11, 0.1*min(height, width));
	for(int i =0; i< numFeatures; i++)
	{
		vx.push_back(f2[i].x-f1[i].x); 
		vy.push_back(f2[i].y-f1[i].y);
	}
	sort( vx.begin( ), vx.end( ) );
	sort( vy.begin( ), vy.end( ) );
	VecSmooth(smvx, vx, 5);
	VecSmooth(smvy, vy, 5);

	double xDiff =smvx[smvx.size()-1]-smvx[0];
	int s = smvx.size ()/5;
	int pos = 0;
	for(int i =0; i< smvx.size()-s; i++)
	{
		if(smvx[s+i]-smvx[i] < xDiff)
		{
			xDiff = smvx[s+i]-smvx[i];
			pos = i+s/2;
		}
	}
	m.x = smvx[pos];

	double yDiff =smvy[smvy.size()-1]-smvy[0];
	pos = 0;
	for(int i =0; i< smvy.size ()-s; i++)
	{
		if(smvy[s+i]-smvy[i] < yDiff)
		{
			yDiff = smvy[s+i]-smvy[i];
			pos = i+s/2;
		}
	}
	m.y = smvy[pos];

	// find rotation
	
	// Get good matching points
	std::vector<int> vi;
	float err_limt = 1.0;
	do{
		for(int i=0; i<numFeatures; i++)
		{
			//if(abs(f2[i].x-f1[i].x-m.x)<5*xDiff && abs(f2[i].y-f1[i].y-m.y)<5*yDiff)
			if(abs(f2[i].x-f1[i].x-m.x)<err_limt && abs(f2[i].y-f1[i].y-m.y)<err_limt)
				vi.push_back(i);
		}

		s = vi.size();
		if(s<20) {
			vi.clear();
			err_limt *=2;
		}
	} while(s<20);
	// Estimate rotation
	CvMat* XTrans = cvCreateMat(3,3,CV_32FC1);
	CvMat* tmp1 = cvCreateMat(3,3,CV_32FC1);
	CvMat* tmp2 = cvCreateMat(3,3,CV_32FC1);
	float Ang=0;
	// estimation and refining
		CvMat* M1 = cvCreateMat(3,s,CV_32FC1);
		float *data1 = M1->data.fl;
		CvMat* M2 = cvCreateMat(3,s,CV_32FC1);
		float *data2 = M2->data.fl;
	
		// [x'] = [ cos(a) -sin(a) ] [x]   [Tx]
		// [y'] = [ sin(a)  cos(a) ] [y] + [Ty]
		// is equivalent to
		// [x'] = [ cos(a) -sin(a) Tx/w] [x]
		// |y'| = | sin(a)  cos(a) Ty/w| |y|
		// [ w] = [ 0        0      1]   [w]
		float w =100;
		for(int i =0; i<s; i++)
		{
			data1[i] = f1[vi[i]].x-width/2;
			data1[i+s] = f1[vi[i]].y - height/2;
			data1[i+s+s] = w;
			data2[i] = f2[vi[i]].x - width/2;
			data2[i+s] = f2[vi[i]].y - height/2;
			data2[i+s+s] = w;
		}
	
		CvMat* M3 = cvCreateMat(s,3,CV_32FC1);
		cvTranspose(M1, M3);
	
		cvMatMul(M1, M3, tmp2);
		cvInvert(tmp2, tmp1);
		cvMatMul(M2, M3, tmp2);
		cvMatMul(tmp2, tmp1, XTrans);
		
		cvReleaseMat(&M1);
		cvReleaseMat(&M2);
		cvReleaseMat(&M3);
	
		float C,S;
		C = XTrans->data.fl[0] + XTrans->data.fl[4];
		S = XTrans->data.fl[3] - XTrans->data.fl[1];
		Ang = atan2(S/C,1);

		m.x = XTrans->data.fl[2]*w;
		m.y = XTrans->data.fl[5]*w;

	cvReleaseMat(&tmp1);
	cvReleaseMat(&tmp2);
	cvReleaseMat(&XTrans);
	return Ang;
	 // */
}

static void SmoothRotation(std::vector<float> & vecSmRotation, const std::vector<float> &vecRotation)
{
	// /*
	int Len = 101;
	float * filter = new float[Len];
	float tlen =0;
	for(int i = 0; i < Len; i++)
	{
		filter[i] = sin(pi*(i+0.5f)/(Len+1));
		tlen += filter[i];
	}
	int last = (int) vecRotation.size()-1;

	float avg, sum;
	for(int j = 0; j<= last; j++)
	{
		sum = 0;
		for(int i =0; i<Len; i++)
		{
			int pos = j+i-Len/2;
			if(pos <0) pos = -pos;
			//if(pos>last) pos = 2*last -pos;
			if(pos>last) {sum = avg*tlen; break;}
			sum += filter[i]*vecRotation[pos];
		}
			
		avg= sum/tlen;
		
		vecSmRotation.push_back(avg);;
	}
	delete filter;
	// */
	
	 /* 
	int fLen = 100;
	int last = (int) vecRotation.size()-1;
	fLen = min(last,fLen);

	// initialisation
	float sum = 0, avg = 0, deviation;
	for(int i =0; i<fLen; i++)
		sum += vecRotation[i];
	avg = sum/fLen;

	// smoothing
	float a = 0.98f;
	for(int j = 0; j<= last; j++)
	{
		avg = a*avg + (1-a)*vecRotation[j];
		vecSmRotation.push_back(avg);
	}
	// */
}

static void GetShifts(std::vector<CvPoint2D32f> & vecShifts, const std::vector<CvPoint2D32f> &vecMotion)
{
	// /*
	int Len = 100;
	int last = (int) vecMotion.size()-1;

	float * filter = new float[Len];
	float tlen =0;
	for(int i = 0; i < Len; i++){
		//filter[i] = 1-cos(2*pi*(i+0.5f)/Len);
		filter[i] = sin(pi*(i+0.5f)/Len);
		tlen += filter[i];
	}

	// initialisation
	CvPoint2D32f p;
	double sumX = 0;
	double sumY = 0;

	for(int j = 0; j<= last; j++)
	{
		sumX = sumY = 0;
		for(int i =0; i<Len; i++)
		{
			int pos = j+i-Len/2;
			if(pos <0) pos = -pos;
			//if(pos>last) pos = 2*last -pos;
			if(pos>last) {sumX = p.x*tlen; sumY=p.y*tlen; break;}
			sumX += filter[i]*vecMotion[pos].x;
			sumY += filter[i]*vecMotion[pos].y;
		}
		p.x = sumX/tlen;
		p.y = sumY/tlen;
		vecShifts.push_back(p);;
	}
	delete filter;
	// */

	 /*
	int Len = 101;
	int last = (int) vecMotion.size()-1;
	Len = min(Len,last);
	float sumX = 0;
	float sumY = 0;
	for(int i =0; i<Len; i++)
	{
		sumX += vecMotion[i].x;
		sumY += vecMotion[i].y;
	}

	// get average
	CvPoint2D32f p;
	p.x = sumX/Len;
	p.y = sumY/Len;
	float a = 0.98f;
	for(int j = 0; j<= last; j++)
	{
		//p.x = sumX/Len;
		//p.y = sumY/Len;
		p.x = a*p.x + (1-a)*vecMotion[j].x;
		p.y = a*p.y + (1-a)*vecMotion[j].y;
		vecShifts.push_back(p);
		//sumX += (vecMotion[min(last, j+Len+1)].x - vecMotion[j].x);
		//sumY += (vecMotion[min(last, j+Len+1)].y - vecMotion[j].y);
	} // */
}

static void Shift(CvPoint2D32f * p, const IplImage *frame, IplImage *shift, float ang=0)
{
	int height	 = frame->height;
	int width	 = frame->width;
	int step	 = frame->widthStep;  
	int channels = frame->nChannels;
	uchar *dataF = (uchar *)frame->imageData;
	uchar *dataS = (uchar *)shift->imageData;

	float fx, fy, x, y;
	if(ang == 0) {
		for(int i =0; i< height; i++)
		{
			y = i - p->y;
			for(int j=0; j< width; j++)
			{
				x = j - p->x;
				if(x < 0 || (int)(x+1) >= width || y < 0 || (int)(y+1) >= height)
				{
					for(int k=0;k<channels;k++)
					{
						dataS[i*step + j*channels+k] = 0;
					}
					continue;
				}
				for(int k=0;k<channels;k++)
				{
					int ix = (int) x, iy= (int)y;
					if(iy+1<height && ix+1<width) 
					dataS[i*step + j*channels+k] = (int)
					( (dataF[iy*step + (ix+1)*channels+k]*(x-ix) + dataF[iy*step + ix*channels+k]*(ix+1-x))*(1+iy-y)
					+ (dataF[(iy+1)*step + (ix+1)*channels+k]*(x-ix) + dataF[(iy+1)*step + ix*channels+k]*(ix+1-x))*(y-iy) );
					else 
						dataS[i*step + j*channels+k] = dataF[iy*step + ix*channels+k];
				}
			}
		}
	}
	else //ang !=0
	{
		float C = cos(ang), S = -sin(ang);
		for(int i =0; i< height; i++)
		{
			fy = i - p->y;
			for(int j=0; j< width; j++)
			{
				fx = j - p->x;
				x = C*(fx-width/2) - S*(fy-height/2)+width/2;
				y = S*(fx-width/2) + C*(fy-height/2) + height/2;
				if(x < 0 || (int)(x+1) >= width || y < 0 || (int)(y+1) >= height)
				{
					for(int k=0;k<channels;k++)
					{
						dataS[i*step + j*channels+k] = 0;
					}
					continue;
				}
				for(int k=0;k<channels;k++)
				{
					int ix = (int) x, iy= (int)y;
					if(iy+1<height && ix+1<width) 
					dataS[i*step + j*channels+k] = (int)
					( (dataF[iy*step + (ix+1)*channels+k]*(x-ix) + dataF[iy*step + ix*channels+k]*(ix+1-x))*(1+iy-y)
					+ (dataF[(iy+1)*step + (ix+1)*channels+k]*(x-ix) + dataF[(iy+1)*step + ix*channels+k]*(ix+1-x))*(y-iy) );
					else 
						dataS[i*step + j*channels+k] = dataF[iy*step + ix*channels+k];
				}
			}
		}
	}

}

static CvRect GetMaxRect(IplImage* img1, IplImage* img2)
{
	CvPoint points[4];
	int w = img1->width; 
	int h = img1->height;
	int c = img1->nChannels;
	uchar* data1 = (uchar*) img1->imageData;
	uchar* data2 = (uchar*) img2->imageData;
	IplImage* tmp = cvCreateImage(cvSize(w,h), IPL_DEPTH_8U, 1);
	uchar* dataT = (uchar*) tmp->imageData;

	for(int i= 0; i< h; i++)
	{	
		int pos = i*img1->widthStep;
		int tpos = i*tmp->widthStep;
		for(int j= 0; j< w; j++)
		{	uchar p = 0;
			for(int k=0;k<c;k++)
			{
				if(data1[pos+j*c+k] >20 && data2[pos+j*c+k] >20)
				{
					p = 255;
					break;
				}
			}
			dataT[tpos+j]= p;
		}
	}

	for(int i=0; i<4; i++) {points[i].x = -1; points[i].y = -1;}

	for(int i= 0; i< h; i++)
	{	
		int pos = i*tmp->widthStep;
		for(int j= 0; j< w; j++)
		{
			if(dataT[pos+j]!=0) 
			{
				if(j==0) {
					if(i==0 || dataT[pos-tmp->widthStep+1] ==0)
					{
						if(points[0].x <0 || (w-j)*(h-i)>(w-points[0].x)*(h-points[0].y))
						{
							points[0].x =j; points[0].y = i;
						}
					}
				}
				else if(i==0 && j>0 && dataT[pos+tmp->widthStep+j-1] ==0)
				{
					if(points[0].x <0 || (w-j)*(h-i)>(w-points[0].x)*(h-points[0].y))
					{
						points[0].x =j; points[0].y = i;
					}
				}
				else if(i>0 && i<h-1 && j>0 && j<w-1 && dataT[pos+tmp->widthStep+j-1] ==0 && dataT[pos-tmp->widthStep+j+1] ==0)
				{
					if(points[0].x <0 || (w-j)*(h-i)>(w-points[0].x)*(h-points[0].y))
					{
						points[0].x =j; points[0].y = i;
					}
				}

				if(j==0) {
					if(i==h-1 || dataT[pos+tmp->widthStep+1] ==0)
						if(points[1].x <0 || (w-j)*i>(w-points[1].x)*points[1].y)
						{
							points[1].x =j; points[1].y = i;
						}
				}
				else if(i==h-1 && j>0 && dataT[pos-tmp->widthStep+j-1] ==0)
				{
					if(points[1].x <0 || (w-j)*i>(w-points[1].x)*points[1].y)
					{
						points[1].x =j; points[1].y = i;
					}
				}
				else if(i>0 && i<h-1 && j>0 &&  j<w-1 && dataT[pos+tmp->widthStep+j+1] ==0 && dataT[pos-tmp->widthStep+j-1] ==0)
				{
					if(points[1].x <0 || (w-j)*i>(w-points[1].x)*points[1].y)
					{
						points[1].x =j; points[1].y = i;
					}
				}

				break;
			}
		}

		for(int j= w-1; j>=0; j--)
		{
			if(dataT[pos+j]!=0) 
			{
				if(j ==w-1) {
					if(i==0 || dataT[pos-tmp->widthStep-1] ==0)
						if(points[2].x <0 || j*(h-i)>points[2].x*(h-points[2].y))
						{
							points[2].x =j; points[2].y = i;
						}
				}
				else if(i==0 && j< w-1 && dataT[pos+tmp->widthStep+j+1] ==0)
				{
					if(points[2].x <0 || j*(h-i)>points[2].x*(h-points[2].y))
					{
						points[2].x =j; points[2].y = i;
					}
				}
				else if(i>0 && i<h-1 && j>0 && j<w-1 && dataT[pos+tmp->widthStep+j+1] ==0 && dataT[pos-tmp->widthStep+j-1] ==0)
				{
					if(points[2].x <0 || j*(h-i)>points[2].x*(h-points[2].y))
					{
						points[2].x =j; points[2].y = i;
					}
				}

				if(j==w-1) {
					if(i==h-1 || dataT[pos+tmp->widthStep+1] ==0)
					if(points[3].x <0 || j*i> points[3].x*points[3].y) 
					{
						points[3].x =j; points[3].y = i;
					}
				}
				else if(i==h-1 && j<w-1 && dataT[pos-tmp->widthStep+j+1] ==0)
				{
					if(points[3].x <0 || j*i> points[3].x*points[3].y) 
					{
						points[3].x =j; points[3].y = i;
					}
				}
				else if(i>0 && i<h-1 && j>0 &&  j<w-1 && dataT[pos+tmp->widthStep+j-1] ==0 && dataT[pos-tmp->widthStep+j+1] ==0)
				{
					if(points[3].x <0 || j*i> points[3].x*points[3].y) 
					{
						points[3].x =j; points[3].y = i;
					}
				}

				break;
			}
		}
	}
	
	cvReleaseImage(&tmp);

	points[0].x = max(points[0].x, points[1].x);
	points[0].y = max(points[0].y, points[2].y);
	points[3].x = min(points[2].x, points[3].x);
	points[3].y = min(points[1].y, points[3].y);

	return cvRect(points[0].x, points[0].y, points[3].x-points[0].x, points[3].y-points[0].y);
}

static void caiGetSubImage(IplImage* img, IplImage* subImg, CvRect r)
{
	int c = subImg->nChannels;
	int w = subImg->width;
	int h = subImg->height;

	uchar* data0 = (uchar*) img->imageData;
	uchar* data = (uchar*) subImg->imageData;
	
	for(int i= 0; i< h; i++)
	{	
		int pos = i*subImg->widthStep;
		int pos0 = (i+r.y)*img->widthStep;
		for(int j= 0; j< w; j++)
		{	
			for(int k=0;k<c;k++)
			{
				data[pos+j*c+k] = data0[pos0+(j+r.x)*c+k];
			}
		}
	}
}

static void Norm(IplImage *frame1_1C, IplImage *frame2_1C, IplImage * pyramid1, IplImage * pyramid2, CvSize imSize)
{
	//  /*
	// pyramids are temporarily used for normalisation.
	cvSmooth(frame1_1C, pyramid1, CV_BLUR, 201,201,1,1);
	cvSmooth(frame2_1C, pyramid2, CV_BLUR, 201,201,1,1);
	// normalisation
	uchar * tp1 = (uchar*) pyramid1->imageData;
	uchar * tp2 = (uchar*) pyramid2->imageData;
	uchar * tp0 = (uchar*) frame2_1C->imageData;
	for(int i= 0; i< imSize.height; i++)
	{	
		int pos = i*frame1_1C->widthStep;
		for(int j= 0; j< imSize.width; j++)
		{	
			int pixel1 = tp1[pos+j]; int pixel2 = tp2[pos+j]; 
			int pixel0 = tp0[pos+j];
			pixel0 = pixel0 - (pixel2 - pixel1);
			tp0[pos+j] = (uchar) min(255,max(0, pixel0));
		}
	}
	// */

	 /*
	uchar * tp1 = (uchar*) frame1_1C->imageData;
	uchar * tp2 = (uchar*) frame2_1C->imageData;
	double sum =0, sumW =0, w;
	for(int i= 0; i< imSize.height; i++)
	{	int pos = i*frame1_1C->widthStep;
		for(int j= 0; j< imSize.width; j++)
		{	
			double pixel1 = tp1[pos+j]; double pixel2 = tp2[pos+j]; 
			if(pixel1>1 && pixel2>1) {
				w = min(1,(pixel1+pixel2)/128);
				sum += (pixel2-pixel1)*w;
				sumW += w;
			}
		}
	}
	double mean = sum/sumW;
	for(int i= 0; i< imSize.height; i++)
	{	int pos = i*frame1_1C->widthStep;
		for(int j= 0; j< imSize.width; j++)
		{	
			double pixel2 = tp2[pos+j]; 
			if(pixel2>1) {
				w = 1;//min(1,pixel2/64);
				pixel2 -= (mean*w);
				tp2[pos+j] = (uchar) min(255,max(0, pixel2));
			}
		}
	}
	// */
}


static void gridFlow(IplImage *img1, IplImage *img2, CvPoint2D32f *f1, CvPoint2D32f *f2, 
				CvPoint2D32f *f21, int & numFeatures, IplImage * pyramid1, IplImage * pyramid2, 
				CvSize flow_Win, char * flow_found, float * flow_error, float tolerance,
				IplImage *dMap, int minD=0, int maxD = 100)
{
	CvTermCriteria optical_flow_termination_criteria
			= cvTermCriteria( CV_TERMCRIT_ITER | CV_TERMCRIT_EPS, 100, .5 );

	float scale = 250.0f/(maxD-minD);
	int level = 3;
	cvCalcOpticalFlowPyrLK(img1, img2, pyramid1, pyramid2, f1, f2, numFeatures, 
		flow_Win, level, flow_found, flow_error, optical_flow_termination_criteria, 0 );
	
	cvCalcOpticalFlowPyrLK(img2, img1, pyramid1, pyramid2, f2, f21, numFeatures, 
		flow_Win, level, flow_found, flow_error, optical_flow_termination_criteria, 0 );

	int featurePos = 0;
	for(int i = 0; i < numFeatures; i++)
	{
		if(abs(f21[i].x-f1[i].x)>2*tolerance || abs(f21[i].y-f1[i].y) >tolerance
			|| abs(f1[i].y-f2[i].y) >tolerance)
		{
			// the flow of the point need to be re-estimated.
			f1[featurePos] = f1[i];
			f2[featurePos] = f2[i];
			featurePos++;
			continue;
		}
		
		// take the estimate of the flow
		int pos = (int) (f1[i].y)*(dMap->widthStep)+ (int)(f1[i].x);
		float tmp = min(maxD-minD, max(0,f2[i].x-f1[i].x-minD));
		dMap->imageData[pos] = (uchar) (tmp*scale);
	}

	numFeatures = featurePos;
}


/* This is just an inline that allocates images.  I did this to reduce clutter in the
 * actual computer vision algorithmic code.  Basically it allocates the requested image
 * unless that image is already non-NULL.  It always leaves a non-NULL image as-is even
 * if that image's size, depth, and/or channels are different than the request.
 */
inline static void allocateOnDemand( IplImage **img, CvSize size, int depth, int channels )
{
	if ( *img != NULL )	return;

	*img = cvCreateImage( size, depth, channels );
	if ( *img == NULL )
	{
		fprintf(stderr, "Error: Couldn't allocate image.  Out of memory?\n");
		exit(-1);
	}
}

int main(int argc, char** argv)
{
	char* outfile1 = "imageA.jpg";
	char* outfile2 = "imageB.jpg";
	
	IplImage * imageA = NULL, * imageB = NULL;
    if( (imageA = cvLoadImage( argv[1], 1)) == 0 )
        return -1;
	if( (imageB = cvLoadImage( argv[2], 1)) == 0 )
        return -1;

	
	CvRect r = GetMaxRect(imageA, imageB);
	/* Create a windows called "Optical Flow" for visualizing the output.
	 * Have the window automatically change its size to match the output.
	 */
	cvNamedWindow("EpipolarA", CV_WINDOW_AUTOSIZE);
	//cvNamedWindow("EpipolarB", CV_WINDOW_AUTOSIZE);

	
	if(imageA->width != imageB->width || imageA->height != imageB->height)
	{
		fprintf(stderr, "Error: image sizes are not same.\n");
		return -1;
	}
	
	CvSize imSize = cvSize(r.width, r.height);

		static IplImage *frame = NULL, *frame1 = NULL, *frame1_1C = NULL, *frame2 = NULL, *frame2_1C = NULL, *eig_image = NULL, *temp_image = NULL, *pyramid1 = NULL, *pyramid2 = NULL;
		static IplImage *smallRect = NULL, *eig_s = NULL, *tmp_s = NULL, * dMap = NULL; 
		
		allocateOnDemand( &frame1, imSize, IPL_DEPTH_8U, 3);
		caiGetSubImage(imageA, frame1, r);
		allocateOnDemand( &frame1_1C, imSize, IPL_DEPTH_8U, 1 );
		cvConvertImage(frame1, frame1_1C);
		cvSaveImage("Left.jpg", frame1);
		allocateOnDemand( &frame2, imSize, IPL_DEPTH_8U, 3);
		caiGetSubImage(imageB, frame2, r);
		allocateOnDemand( &frame2_1C, imSize, IPL_DEPTH_8U, 1 );
		cvConvertImage(frame2, frame2_1C);
		cvSaveImage("Right.jpg", frame2);

		int step = 5;
		int minD = -25;
		int maxD = 60;
		float tolerance = 1.0;
		int winSize = 3;
		int start = step/2;
		int totalNum = (imSize.height/step) * (imSize.width/step);

		// This is some workspace for the algorithm.
		// (The algorithm actually carves the image into pyramids of different resolutions.)
		allocateOnDemand( &pyramid1, imSize, IPL_DEPTH_8U, 1 );
		allocateOnDemand( &pyramid2, imSize, IPL_DEPTH_8U, 1 );
		
		 /*
		IplImage* corners = cvCreateImage(imSize, IPL_DEPTH_32F, 1);
		cvSmooth( frame1_1C, pyramid1, CV_GAUSSIAN, 7);
		cvPreCornerDetect( pyramid1, corners, 7 );
		cvAbs(corners, corners);
		cvNormalize( corners, pyramid1, 0, 255, CV_MINMAX );
		cvAdd(pyramid1, frame1_1C, frame1_1C);
		//cvCanny( pyramid1, frame1_1C, 100, 255, 7 );
		//cvCopyImage(pyramid1, frame1_1C);
		cvSmooth( frame2_1C, pyramid1, CV_GAUSSIAN, 7);
		cvPreCornerDetect( pyramid1, corners, 7 );
		cvAbs(corners, corners);
		cvNormalize( corners, pyramid1, 0, 255, CV_MINMAX );
		cvAdd(pyramid1, frame2_1C, frame2_1C);
		//cvCanny( pyramid1, frame2_1C, 100, 255, 7 );
		//cvCopyImage(pyramid1, frame2_1C);
		//cvMax(frame1_1C, frame2_1C, pyramid2);
		cvSaveImage("Edge.jpg", frame1_1C);
		cvReleaseImage( &corners);
		// */


		CvPoint2D32f* frame1_features = new CvPoint2D32f[totalNum];
		
		int it = 0;
		for(int i= step/2; i<= imSize.height-step+start; i += step)
		{
			for(int j= step/2; j<= imSize.width-step+start; j += step)
			{	if((uchar)frame1_1C->imageData[i*frame1_1C->widthStep+j] ==0) continue;
				if((uchar)frame2_1C->imageData[i*frame2_1C->widthStep+j] ==0) continue;
				frame1_features[it].x = j;
				frame1_features[it].y = i;
				it++;
			}
		}
		totalNum = it;
		int numfeatures = totalNum;
		/* Shi and Tomasi Feature Tracking! */
		 /*
		//Preparation: Allocate the necessary storage.
		allocateOnDemand( &eig_image, imSize, IPL_DEPTH_32F, 1 );
		allocateOnDemand( &temp_image, imSize, IPL_DEPTH_32F, 1 );
		cvGoodFeaturesToTrack(frame1_1C, eig_image, temp_image, frame1_features, &numfeatures, .01, .01, NULL);
		// */ 
		
		
		/* Pyramidal Lucas Kanade Optical Flow! */
		/* This array will contain the locations of the points from frame 1 in frame 2. */
		CvPoint2D32f* frame2_features = new CvPoint2D32f[totalNum];
		CvPoint2D32f* frame21_features = new CvPoint2D32f[totalNum];

		char* optical_flow_found_feature = new char[totalNum];
		float* optical_flow_feature_error = new float[totalNum];




		//cvCopyImage(frame1_1C, pyramid1);
		//cvEqualizeHist( pyramid1, frame1_1C);
		//cvCopyImage(frame2_1C, pyramid1);
		//cvEqualizeHist( pyramid1, frame2_1C);
		Norm(frame1_1C, frame2_1C, pyramid1, pyramid2, imSize);

		allocateOnDemand( &dMap, imSize, IPL_DEPTH_8U, 1 );
		for(int loop = 0; loop <4; loop++)
		{
			gridFlow(frame1_1C, frame2_1C, frame1_features, frame2_features, 
				frame21_features, numfeatures, pyramid1, pyramid2, 
				cvSize(winSize,winSize), optical_flow_found_feature, optical_flow_feature_error, 
				tolerance, dMap, minD, maxD);
			tolerance += 0.5;
			winSize = winSize*1.8;
			winSize = (winSize/2)*2+1;
		}

	// post processing
	float scale = 250.0f/(maxD-minD);
	for(int i = 0; i < numfeatures; i++)
	{
		int pos = (int) (frame1_features[i].y)*(dMap->widthStep)+ (int)(frame1_features[i].x);
		float tmp = min(maxD-minD, max(0,frame2_features[i].x-frame1_features[i].x-minD));
		dMap->imageData[pos] = (uchar) (tmp*scale);
	}
	
	if(step > 1)
	{
		int start = step/2;
		uchar* iData = (uchar*) (dMap->imageData);
		for(int j= start; j< imSize.height-step; j++)
		{	
			int jj = ((j-start)/step)*step + start;
			float j1 = j - jj;
			float j2 = step - j1;
			int pos = dMap->widthStep*j;
			for(int i= start; i< imSize.width-step; i++)
			{
				int ii = ((i-start)/step)*step + start;
				//if(j==jj && i==ii) continue;
				float i1 = i - ii;
				float i2 = step - i1;
				float tmp = (iData[jj*dMap->widthStep+ii]*i2*j2 
					+ iData[jj*dMap->widthStep+ii+step]*i1*j2
					+ iData[(jj+step)*dMap->widthStep+ii]*i2*j1
					+ iData[(jj+step)*dMap->widthStep+ii+step]*i1*j1)/(step*step);
				dMap->imageData[pos+i] = (uchar) (tmp);
			}
		}
	}

		//IplImage * vx = NULL, *vy = NULL;
		//allocateOnDemand( &vx, imSize, IPL_DEPTH_32F, 1 );
		//allocateOnDemand( &vy, imSize, IPL_DEPTH_32F, 1 );
		//cvCalcOpticalFlowHS(frame1_1C, frame2_1C, 0, vx, vy, 20,
        //                  optical_flow_termination_criteria );
		//cvConvertScale(vx, frame2_1C,20);
		//cvSaveImage("Testt.jpg", frame2_1C);

		

	cvShowImage("EpipolarA", dMap);
	//cvShowImage("EpipolarB", frame2);
	cvSaveImage("TestD.jpg", dMap);
	//cvSaveImage("TestB.jpg", frame2);
		/*
	CvMat* points1, *points2;
	CvMat* status, * fundamental_matrix;

	points1 = cvCreateMat(1,numfeatures,CV_32FC2);
	points2 = cvCreateMat(1,numfeatures,CV_32FC2);
	status = cvCreateMat(1,numfeatures,CV_8UC1);

	// Fill the points here ... 
	for( int i = 0; i < numfeatures; i++ )
	{
		points1->data.fl[2*i] = frame1_features[i].x;
		points1->data.fl[2*i+1] = frame1_features[i].y;
		points2->data.fl[2*i] = frame2_features[i].x;
		points2->data.fl[2*i+1] = frame2_features[i].y;
		//cvSetReal2D(points2,0,i,frame2_features[i].x);
		//cvSetReal2D(points2,1,i,frame2_features[i].y);
	}

	fundamental_matrix = cvCreateMat(3,3,CV_32FC1);
	int fm_count = cvFindFundamentalMat( points1,points2,fundamental_matrix, CV_FM_RANSAC,1,0.999,status );
	//int fm_count = cvFindFundamentalMat( points1,points2,fundamental_matrix, CV_FM_LMEDS,0.5,0.999,status );
	if( fm_count > 0 )
		printf("Fundamental matrix was found\n");
	else
	{
		printf("Fundamental matrix was not found\n");
		return -1;
	}
	
	int matchedCount =0;
	for( int i = 0; i < numfeatures; i++)
	{
		if(status->data.s[i] >0) matchedCount++;
		if(matchedCount>=20) break;
	}
	CvMat* pointsc1 = cvCreateMat(1,matchedCount,CV_32FC2);
	CvMat* pointsc2 = cvCreateMat(1,matchedCount,CV_32FC2);
	for( int i = 0, j =0; i < numfeatures, j < matchedCount; i++)
	{
		if(status->data.s[i] <1) continue;
		pointsc1->data.fl[2*j] = frame1_features[i].x;
		pointsc1->data.fl[2*j+1] = frame1_features[i].y;
		pointsc2->data.fl[2*j] = frame2_features[i].x;
		pointsc2->data.fl[2*j+1] = frame2_features[i].y;
		j++;
	}
	CvMat* H1 = cvCreateMat(3,3,CV_32FC1); 
	CvMat* H2 = cvCreateMat(3,3,CV_32FC1);
	cvStereoRectifyUncalibrated(pointsc1,pointsc2,fundamental_matrix, imSize, H1, H2, 1);
	cvWarpPerspective( imageA, frame1, H1);
	cvWarpPerspective( imageB, frame2, H2);

	cvShowImage("EpipolarA", frame1);
	cvShowImage("EpipolarB", frame2);
	//cvSaveImage("RectifiedA.jpg", frame1);
	//cvSaveImage("RectifiedB.jpg", frame2);
	cvConvertImage(frame1, frame1_1C);
	cvConvertImage(frame2, frame2_1C);
	// */
	
	// /*
	// stereo matching
	CvMat* disp = cvCreateMat( imSize.height, imSize.width, CV_16S );
    CvMat* vdisp = cvCreateMat( imSize.height, imSize.width, CV_8U );
		//Setup for finding stereo corrrespondences
        CvStereoBMState *BMState = cvCreateStereoBMState();
        assert(BMState != 0);
        BMState->preFilterSize=41;
        BMState->preFilterCap=31;
        BMState->SADWindowSize=41;
        BMState->minDisparity= minD;
        BMState->numberOfDisparities=128;
        BMState->textureThreshold=10;
        BMState->uniquenessRatio=15;
    
	cvFindStereoCorrespondenceBM( frame2_1C, frame1_1C, disp,BMState);
    cvNormalize( disp, vdisp, 0, 256, CV_MINMAX );
    cvNamedWindow( "disparity" );
    cvShowImage( "disparity", vdisp );
	cvSaveImage("disparity.jpg", vdisp);
	// */

	 /* gragh cut
	CvMat* disparity_left = cvCreateMat( imSize.height, imSize.width, CV_16S );
	CvMat* disparity_right = cvCreateMat( imSize.height, imSize.width, CV_16S );
	CvStereoGCState* state = cvCreateStereoGCState( 16, 4 );
	cvFindStereoCorrespondenceGC( frame2_1C, frame1_1C, disparity_left, disparity_right, state, 0 );
	cvReleaseStereoGCState( &state );
	IplImage * disparity_left_visual = NULL;
	allocateOnDemand( &disparity_left_visual, imSize, IPL_DEPTH_8U, 1 );
	cvNormalize( disparity_left, disparity_left_visual, 0, 256, CV_MINMAX );
	cvConvertScale( disparity_left, disparity_left_visual, 10 );
	cvNamedWindow( "disparity" );
    cvShowImage( "disparity", disparity_left_visual );
	cvSave( "disparity.jpg", disparity_left_visual );
	// */
	cvWaitKey();
			delete frame1_features;
			delete frame2_features;
			delete frame21_features;
			delete optical_flow_found_feature;
			delete optical_flow_feature_error;
			cvReleaseImage(&frame1);
			cvReleaseImage(&frame1_1C);
			cvReleaseImage(&frame2);
			cvReleaseImage(&frame2_1C);
			cvReleaseImage(&eig_image);
			cvReleaseImage(&temp_image);
			cvReleaseImage(&pyramid1);
			cvReleaseImage(&pyramid2);
			cvReleaseImage(&smallRect);
			cvReleaseImage(&eig_s);
			cvReleaseImage(&tmp_s);

	cvDestroyWindow("Optical Flow");
	//cvReleaseImage(&rectified);
	cvReleaseImage(&frame);
}