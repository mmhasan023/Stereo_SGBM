#include "opencv2/calib3d/calib3d.hpp"
#include "opencv2/imgproc/imgproc.hpp"
#include "opencv2/imgcodecs.hpp"
#include "opencv2/highgui/highgui.hpp"
#include "opencv2/core/utility.hpp"
#include <iostream>"
#include <stdio.h>
#include <vector>
#include<fstream>

using namespace cv;
using namespace std;

#ifndef Max
#define Max(a,b)            (((a) > (b)) ? (a) : (b))
#endif

#ifndef Min
#define Min(a,b)            (((a) < (b)) ? (a) : (b))
#endif

typedef struct {
	Point2f points[4];
} Coordinates;

vector<string> split(const string& str, const string& delim)
{
	vector<string> tokens;
	size_t prev = 0, pos = 0;
	do
	{
		pos = str.find(delim, prev);
		if (pos == string::npos) pos = str.length();
		string token = str.substr(prev, pos - prev);
		if (!token.empty()) tokens.push_back(token);
		prev = pos + delim.length();
	} while (pos < str.length() && prev < str.length());
	return tokens;
}

void readMaskCoordinates(string fileName, vector<string> & sDate, vector<vector<Coordinates>> & vCoordinates)
{
	ifstream myfile(fileName);
	string line0, line1, sPlot = "_2016";
	vector<string> vStr;
	int pos;

	if (myfile.is_open())
	{
		while (myfile.good())
		{
			getline(myfile, line0);
			pos = line0.find(sPlot);
			if (pos < 0) continue;
			else {
				std::string token;
				stringstream ss(line0);

				getline(ss, token, ','); // remove the first token
				while (getline(ss, token, ',')) {
					if (token.length() <= 3) continue;
					string str = token;
					//transform(str.begin(), str.end(), str.begin(), ::toupper);
					if (str[0] == '0') str = str.substr(1);
					sDate.push_back(str);
				}

				vCoordinates.resize(sDate.size());
				getline(myfile, line0); // skip a line

				for (int i = 0; i < 60; i++)
				{
					getline(myfile, line0);
					vector<string> tokens = split(line0, ",");
					for (int j = 0; j < sDate.size(); j++)
					{
						Coordinates c; 
						c.points[0].x = atof(tokens[4 * j + 1].c_str()); c.points[0].y = 0;
						c.points[1].x = atof(tokens[4 * j + 2].c_str()); c.points[1].y = 3456;
						c.points[2].x = atof(tokens[4 * j + 3].c_str()); c.points[2].y = 3456;
						c.points[3].x = atof(tokens[4 * j + 4].c_str()); c.points[3].y = 0;
						vCoordinates[j].push_back(c);
					}
				}
				break;
			}
		}

		myfile.close();
	}
	else
		cout << "Can't open file" << endl;

}



void readHeight(string fileName, vector<string> & sDate, vector<vector<int>> & vHeight)
{
	ifstream myfile(fileName);
	string line0, line1, sPlot = "Plot";
	vector<string> vStr;
	int pos;

	if (myfile.is_open())
	{

		while (myfile.good())
		{
			getline(myfile, line0);
			pos = line0.find(sPlot);
			if (pos < 0) continue;
			else {
				std::string token;
				stringstream ss(line0);

				getline(ss, token, ','); // remove the first token
				while (getline(ss, token, ',')) {
					if (token.length() <= 3) continue;
					string str = token;
					//transform(str.begin(), str.end(), str.begin(), ::toupper);
					if (str[0] == '0') str = str.substr(1);
					str = str + "_2016";
					sDate.push_back(str);
				}

				vHeight.resize(sDate.size());

				for (int i = 0; i < 60; i++)
				{
					getline(myfile, line0);
					stringstream sl(line0);

					getline(sl, token, ','); // remove the first token
					int j = 0;
					while (getline(sl, token, ',')) {
						if (token =="" || token =="  ") continue;
						string str = token;
						int num = atoi(str.c_str());
						vHeight[j].push_back(num);
						j++;
					}

				}
				break;
			}
		}
	
		myfile.close();
	}
	else
		cout << "Can't open file" << endl;

}

static void print_help()
{
    printf("\nDemo stereo matching converting L and R images into disparity and point clouds\n");
    printf("\nUsage: stereo_match <left_image> <right_image> [--algorithm=bm|sgbm|sgbm3way] [--blocksize=<block_size>]\n"
           "[--max-disparity=<max_disparity>] [--scale=scale_factor>] [-i=<intrinsic_filename>] [-e=<extrinsic_filename>]\n"
           "[--no-display] [-o=<disparity_image>] [-p=<point_cloud_file>]\n");
}

inline bool exists(const std::string& name) {
	struct stat buffer;
	return (stat(name.c_str(), &buffer) == 0);
}


void histgram_Stereo(Mat & stereoImg, vector<unsigned int> & hist, Rect & mask, int ground = 20)
{
	hist.resize(256);
	unsigned int total = 0, coverage =0;
	for (int i = 0; i < 256; i++)
	{
		hist[i] = 0;
	}

	for (int x = mask.x; x < mask.x + mask.width; x++)
	{
		for (int y = mask.y; y < mask.y + mask.height; y++)
		{
			uchar c = stereoImg.at<uchar>(y, x);
			if (c > ground / 2 && c< 250) {
				total++;
				hist[c]++;
			}
			if (c > ground*1.5) coverage++;
		}
	}

	unsigned int threshold, sum=0;
	coverage *= 0.1;
	for (int i = 250; i > 0; i--)
	{
		sum += hist[i];
		if (sum >= coverage) 
		{
			threshold = i; break;
		}
	}
	threshold = (threshold - ground)*.2 + ground;
	threshold = Max(threshold, ground*1.5);
	
	coverage = 0;
	for (int i = 0; i < 256; i++) { hist[i] = 0; }

	for (int x = mask.x; x < mask.x + mask.width; x++)
	{
		for (int y = mask.y; y < mask.y + mask.height; y++)
		{
			uchar c = stereoImg.at<uchar>(y, x);
			if (c > threshold && c< 250) {
					hist[c]++;
			}
		}
	}

	ofstream Hist_File("Hist_Stereo.csv");         //Opening file to print info to
	Hist_File << "Coverage;" << 100.0* coverage / total << endl << endl;          //Headings for file
	for (int i = 0; i <= 250; i++) {
		Hist_File << hist[i] << endl;
	}
	Hist_File.close();
}


void writeRectificationMatrix(string filename, Mat & H1, Mat & H2, Mat Homo = Mat::eye(3, 3, CV_64F))
{
	FileStorage fs(filename, FileStorage::WRITE);

	fs << "H1" << H1;
	fs << "H2" << H2;
	fs << "Homo" << Homo;
	fs.release();
	return;
}

int readRectificationMatrix(string filename, Mat & H1, Mat & H2, Mat & Homo)
{
	FileStorage fs(filename, FileStorage::READ);

	fs["H1"] >> H1;
	fs["H2"] >> H2;
	fs["Homo"] >> Homo;
	fs.release();
	if (H1.cols <3 || H2.cols <3) return 0;
	return 1;
}

void trim(string& s)
{
	int pos = -1;
	pos = s.find_first_not_of(" \t");
	if (pos >= 0)
		s = s.substr(pos);
	else s = "";
	pos = s.find_last_not_of(" \t");
	if (pos >= 0)
		s = s.substr(0, pos + 1);
	else s = "";
}

void distortionCorrectionGeneric_Div(IplImage *src, IplImage *dst, std::vector<double> & lensCoef, int kNum, double cx, double cy)
{
	int height = src->height;
	int width = src->width;
	int step = src->widthStep;
	int channels = src->nChannels;
	uchar *dataS = (uchar *)src->imageData;
	uchar *dataD = (uchar *)dst->imageData;

	double x, y, r1, r2, dxx, dxy, dyy, Dr, Dr1, temp;
	double minD = -DBL_MAX;
	int ix, iy;

	// xu-xdc = (xd-xdc)/[1+B2*rd^2 + B3*rd^3...] --> xu = x/D(r);
	// yu-ydc = (yd-ydc)/[1+B2*rd^2 + B3*rd^3...] --> yu = y/D(r);
	//
	// dxu/dx = [D(r)-x^2 \sum k_i*(i+1)*r^(i-1)]/D(r)^2	// dxx;  k = 1, ...
	// dxu/dy = -xy*[\sum k_i*(i+1)*r^(i-1)]/D(r)^2		// dxy
	// dyu/dx = -xy*[\sum k_i*(i+1)*r^(i-1)]/D(r)^2		// dxy
	// dyu/dy = [D(r)-y^2 \sum k_i*(i+1)*r^(i-1)]/D(r)^2	// dyy

	cvZero(dst);
	CvPoint2D32f * posd = new CvPoint2D32f[width*height];
	for (int i = 0; i<width*height; i++)
	{
		posd[i].x = minD; posd[i].y = minD;
	}

	for (int i = 0; i< height; i++)
	{
		for (int j = 0; j< width; j++)
		{
			// dxx = [Dr-x^2*Dr1]/Dr^2
			// dxy = -x*y*Dr1/Dr^2
			// dyy = [Dr-y^2*Dr1]/Dr^2
			r2 = (i - cy)*(i - cy) + (j - cx)*(j - cx);
			r1 = sqrt(r2);
			Dr = 1; Dr1 = 0; temp = 1; // temp is a temp variable
			dxx = dxy = dyy = 0; // fr'(xd), fr'(yd)
			for (int k = 0; k<kNum; k++)
			{
				Dr += lensCoef[k] * r2*temp;
				Dr1 += lensCoef[k] * temp*(k + 2);
				temp = temp*r1;
			}
			dxx = (Dr - (j - cx)*(j - cx)*Dr1) / (Dr*Dr); // dxx  actually is dxu/dx
			dxy = -(j - cx)*(i - cy)*Dr1 / (Dr*Dr); // dxy  actually is dxu/dy = dyu/dx
			dyy = (Dr - (i - cy)*(i - cy)*Dr1) / (Dr*Dr); // dyy  actually is dyu/dy
			x = (j - cx) / Dr + cx; // x is xu
			y = (i - cy) / Dr + cy; // y is yu
			ix = (int)(x + 0.5); iy = (int)(y + 0.5);
			if (ix >= 0 && ix < width && iy >= 0 && iy < height)
			{
				// dx = [dxu * dyy- dyu* dxy]/(dxx*dyy-dxy^2)
				// dy = [dyu * dxx - dxu *dyx]/(dxx*dyy-dxy^2)
				temp = dxx*dyy - dxy*dxy;
				posd[ix + iy*width].x = j + ((ix - x)*dyy - (iy - y)*dxy) / temp;
				posd[ix + iy*width].y = i + ((iy - y)*dxx - (ix - x)*dxy) / temp;
			}
		}
	}

	// fill gaps
	for (;;)
	{
		bool loop = false;
		for (int i = 0; i< height; i++)
		{
			for (int j = 0; j< width; j++)
			{
				long pos = i*width + j;
				if (posd[pos].x>-100) continue; // not a gap
				if (j>0 && j<width - 1 && posd[pos - 1].x>-100 && posd[pos + 1].x>-100)
				{
					posd[pos].x = (posd[pos - 1].x + posd[pos + 1].x) / 2.0;
					posd[pos].y = (posd[pos - 1].y + posd[pos + 1].y) / 2.0;
					if (posd[pos].x >-100) loop = true;
				}
				else if (i>0 && i<height - 1 && posd[pos - width].x>-100 && posd[pos + width].x>-100)
				{
					posd[pos].x = (posd[pos - width].x + posd[pos + width].x) / 2.0;
					posd[pos].y = (posd[pos - width].y + posd[pos + width].y) / 2.0;
					if (posd[pos].x >-100) loop = true;
				}
				else if (j>0 && i>0 && posd[pos - 1 - width].x>-100 && j<width - 1 && i<height - 1 && posd[pos + 1 + width].x>-100)
				{
					posd[pos].x = (posd[pos - 1 - width].x + posd[pos + 1 + width].x) / 2.0;
					posd[pos].y = (posd[pos - 1 - width].y + posd[pos + 1 + width].y) / 2.0;
					if (posd[pos].x >-100) loop = true;
				}
				else if (i>0 && j<width - 1 && posd[pos - width + 1].x>-100 && i<height - 1 && j>0 && posd[pos + width].x>-100)
				{
					posd[pos].x = (posd[pos - width + 1].x + posd[pos + width - 1].x) / 2.0;
					posd[pos].y = (posd[pos - width + 1].y + posd[pos + width - 1].y) / 2.0;
					if (posd[pos].x >-100) loop = true;
				}
			}
		}

		// use dx and dy to fill gaps
		for (int i = 0; i< height; i++)
		{
			for (int j = 0; j< width; j++)
			{
				long pos = i*width + j;
				if (posd[pos].x>minD) continue; // not a gap
				if (j>0 && posd[pos - 1].x >= 0)
				{
					r2 = (posd[pos - 1].y - cy)*(posd[pos - 1].y - cy) + (posd[pos - 1].x - cx)*(posd[pos - 1].x - cx);
					r1 = sqrt(r2);
					Dr = 1; Dr1 = 0; temp = 1; // temp is a temp variable
					dxx = dxy = dyy = 0; // fr'(xd), fr'(yd)
					for (int k = 0; k<kNum; k++)
					{
						Dr += lensCoef[k] * r2*temp;
						Dr1 += lensCoef[k] * temp*(k + 2);
						temp = temp*r1;
					}
					dxx = (Dr - (posd[pos - 1].x - cx)*(posd[pos - 1].x - cx)*Dr1) / (Dr*Dr); // dxx  actually is dxu/dx
					dxy = -(posd[pos - 1].x - cx)*(posd[pos - 1].y - cy)*Dr1 / (Dr*Dr); // dxy  actually is dxu/dy = dyu/dx
					dyy = (Dr - (posd[pos - 1].y - cy)*(posd[pos - 1].y - cy)*Dr1) / (Dr*Dr); // dyy  actually is dyu/dy
					temp = dxx*dyy - dxy*dxy;
					posd[pos].x = posd[pos - 1].x + 1.0*dyy / temp;
					posd[pos].y = posd[pos - 1].y - 1.0*dxy / temp;
					loop = true;
				}
				else if (j<width - 1 && posd[pos + 1].x >= 0)
				{
					r2 = (posd[pos + 1].y - cy)*(posd[pos + 1].y - cy) + (posd[pos + 1].x - cx)*(posd[pos + 1].x - cx);
					r1 = sqrt(r2);
					Dr = 1; Dr1 = 0; temp = 1; // temp is a temp variable
					dxx = dxy = dyy = 0; // fr'(xd), fr'(yd)
					for (int k = 0; k<kNum; k++)
					{
						Dr += lensCoef[k] * r2*temp;
						Dr1 += lensCoef[k] * temp*(k + 2);
						temp = temp*r1;
					}
					dxx = (Dr - (posd[pos + 1].x - cx)*(posd[pos + 1].x - cx)*Dr1) / (Dr*Dr); // dxx  actually is dxu/dx
					dxy = -(posd[pos + 1].x - cx)*(posd[pos + 1].y - cy)*Dr1 / (Dr*Dr); // dxy  actually is dxu/dy = dyu/dx
					dyy = (Dr - (posd[pos + 1].y - cy)*(posd[pos + 1].y - cy)*Dr1) / (Dr*Dr); // dyy  actually is dyu/dy
					temp = dxx*dyy - dxy*dxy;
					posd[pos].x = posd[pos + 1].x - 1.0*dyy / temp;
					posd[pos].y = posd[pos + 1].y + 1.0*dxy / temp;
					loop = true;
				}
				else if (i>0 && posd[pos - width].x >= 0)
				{
					r2 = (posd[pos - width].y - cy)*(posd[pos - width].y - cy) + (posd[pos - width].x - cx)*(posd[pos - width].x - cx);
					r1 = sqrt(r2);
					Dr = 1; Dr1 = 0; temp = 1; // temp is a temp variable
					dxx = dxy = dyy = 0; // fr'(xd), fr'(yd)
					for (int k = 0; k<kNum; k++)
					{
						Dr += lensCoef[k] * r2*temp;
						Dr1 += lensCoef[k] * temp*(k + 2);
						temp = temp*r1;
					}
					dxx = (Dr - (posd[pos - width].x - cx)*(posd[pos - width].x - cx)*Dr1) / (Dr*Dr); // dxx  actually is dxu/dx
					dxy = -(posd[pos - width].x - cx)*(posd[pos - width].y - cy)*Dr1 / (Dr*Dr); // dxy  actually is dxu/dy = dyu/dx
					dyy = (Dr - (posd[pos - width].y - cy)*(posd[pos - width].y - cy)*Dr1) / (Dr*Dr); // dyy  actually is dyu/dy
					temp = dxx*dyy - dxy*dxy;

					posd[pos].x = posd[pos - width].x - 1.0*dxy / temp;
					posd[pos].y = posd[pos - width].y + 1.0*dxx / temp;
					loop = true;
				}
				else if (i<height - 1 && posd[pos + width].x >= 0)
				{
					r2 = (posd[pos + width].y - cy)*(posd[pos + width].y - cy) + (posd[pos + width].x - cx)*(posd[pos + width].x - cx);
					r1 = sqrt(r2);
					Dr = 1; Dr1 = 0; temp = 1; // temp is a temp variable
					dxx = dxy = dyy = 0; // fr'(xd), fr'(yd)
					for (int k = 0; k<kNum; k++)
					{
						Dr += lensCoef[k] * r2*temp;
						Dr1 += lensCoef[k] * temp*(k + 2);
						temp = temp*r1;
					}
					dxx = (Dr - (posd[pos + width].x - cx)*(posd[pos + width].x - cx)*Dr1) / (Dr*Dr); // dxx  actually is dxu/dx
					dxy = -(posd[pos + width].x - cx)*(posd[pos + width].y - cy)*Dr1 / (Dr*Dr); // dxy  actually is dxu/dy = dyu/dx
					dyy = (Dr - (posd[pos + width].y - cy)*(posd[pos + width].y - cy)*Dr1) / (Dr*Dr); // dyy  actually is dyu/dy
					temp = dxx*dyy - dxy*dxy;

					posd[pos].x = posd[pos + width].x + 1.0*dxy / temp;
					posd[pos].y = posd[pos + width].y - 1.0*dxx / temp;
					loop = true;
				}
			}
		}
		if (!loop) break;
	} // for(;;)

	  // 4 sides
	for (int i = 0; i< height; i++)
	{
		long pos = i*width;
		if (posd[pos].x <-100) posd[pos] = posd[pos + 1];
		if (posd[pos + width - 1].x <-100) posd[pos + width - 1] = posd[pos + width - 2];
	}
	for (int j = 0; j< width; j++)
	{
		if (posd[j].x <-100) posd[j] = posd[j + width];
		long pos = height*width - width + j;
		if (posd[pos].x<-100) posd[pos] = posd[pos - width];
	}

	// lens distortion correction
	int ix1, iy1;
	float d[4];
	for (int i = 0; i< height; i++)
	{
		for (int j = 0; j< width; j++)
		{
			ix = floor(posd[i*width + j].x);
			iy = floor(posd[i*width + j].y);
			ix1 = ix + 1; iy1 = iy + 1;
			d[0] = posd[i*width + j].x - ix;  d[1] = ix1 - posd[i*width + j].x;
			d[2] = posd[i*width + j].y - iy;  d[3] = iy1 - posd[i*width + j].y;

			if (ix >= 0 && iy >= 0 && iy1<height && ix1<width)
			{
				for (int k = 0; k<channels; k++)
					dataD[i*step + j*channels + k] = (uchar)
					((dataS[iy*step + ix1*channels + k] * d[0] + dataS[iy*step + ix*channels + k] * d[1])*d[3]
						+ (dataS[iy1*step + ix1*channels + k] * d[0] + dataS[iy1*step + ix*channels + k] * d[1])*d[2]);
			}
			else {
				ix = posd[i*width + j].x + 0.5f;
				iy = posd[i*width + j].y + 0.5f;
				if (ix<0) ix = 0; if (ix >= width) ix = width - 1;
				if (iy<0) iy = 0; if (iy >= height) iy = height - 1;
				for (int k = 0; k<channels; k++)
					dataD[i*step + j*channels + k] = dataS[iy*step + ix*channels + k];
			}
		}
	}
	delete[] posd;
}

void distortionCorrection_Div(Mat & src, Mat & dst, std::vector<double> & lensCoef, int kNum, double cx, double cy)
{
	IplImage *imSrc = cvCloneImage(&(IplImage)src); IplImage * imDst = cvCloneImage(imSrc);
	vector<double> betaVec;
	betaVec.push_back(lensCoef[0]);
	for (int i = 1; i< lensCoef.size(); i++)
	{
		betaVec.push_back(0);
		betaVec.push_back(lensCoef[i]);
	}
	distortionCorrectionGeneric_Div(imSrc, imDst, betaVec, kNum, cx, cy);
	dst = cvarrToMat(imDst, true);

	cvReleaseImage(&imSrc);
	cvReleaseImage(&imDst);
}

// estimate dtm when plant is green
void DTM(const Mat & img, const Mat & dMap, Mat & dtm, int threshold = 10)
{
	dtm = dMap & (dMap <= threshold);
	Mat tmp = dtm.clone();
	Vec3b p;
	for (int x = 0; x < img.cols; x++)
	{
		for (int y = 0; y < img.rows; y++)
		{
			p = img.at<Vec3b>(y, x); // bgr
			if (p[1] >= p[0] && p[1] >= p[2]) { // green pixel
				dtm.at<uchar>(y, x) = 0;
			}
		}
	}
	imwrite("Test1.png", dtm);
	int size = 15;

	uchar u, u1;
	bool bloop;
	for (int loop = 0; loop < 6; loop++)
	{
		bloop = false;
		for (int x = 0; x < dMap.cols; x++)
		{
			for (int y = 0; y < dMap.rows; y++)
			{
				int minv = 255;
				if (dtm.at<uchar>(y, x) > 1) continue;
				bloop = true;
				for (int k = -size; k <= size; k++)
				{
					if (y + k < 0 || y + k >= dMap.rows) continue;

					u = dtm.at<uchar>(y + k, x);
					if (u > 0 && minv > u) minv = u;
				}
				if (minv < 255) { dtm.at<uchar>(y, x) = minv; }
			}
		}

		for (int y = 0; y < dMap.rows; y++)
		{
			for (int x = 0; x < dMap.cols; x++)
			{
				if (dtm.at<uchar>(y, x) > 1) continue;
				bloop = true;
				int minv = 255;
				for (int k = -size; k <= size; k++)
				{
					if (x + k < 0 || x + k >= dMap.cols) continue;

					u = dtm.at<uchar>(y, x + k);
					if (u > 0 && minv > u) minv = u;
				}
				if (minv < 255) { dtm.at<uchar>(y, x) = minv; }
			}
		}
		if (bloop)continue;
		break;
	}
	Mat element = getStructuringElement(MORPH_RECT, Size(2 * 2 + 1, 2 * 2 + 1), Point(2, 2));
	dilate(dtm, dtm, element);
	imwrite("Test2.png", dtm);
	erode(dtm, dtm, element);
	element = getStructuringElement(MORPH_RECT, Size(2 * (size / 2) + 1, 2 * (size / 2) + 1), Point(size / 2, size / 2));

	erode(dtm, tmp, element);
	dilate(tmp, dtm, element);
	element = getStructuringElement(MORPH_RECT, Size(2 * size + 1, 2 * size + 1), Point(size, size));
	dilate(dtm, tmp, element);
	dilate(tmp, dtm, element);
	erode(dtm, dtm, element);//dilate(tmp, dtm, element);

							 // restore the final detail
							 /*
							 for (int x = 0; x < dMap.cols; x++)
							 {
							 for (int y = 0; y < dMap.rows; y++)
							 {
							 u = dtm.at<uchar>(y, x);
							 u1 = dMap.at<uchar>(y, x);
							 if (u1 <= 1) continue;
							 if ( u1 < u ) dtm.at<uchar>(y, x) = u1;
							 }
							 }
							 // */
	//imwrite("Test.png", dtm);
	imwrite("dtm.png", dtm);
}

void DTM(const Mat & dMap, Mat & dtm, int threshold = 10)
{
	dtm = dMap & (dMap <= threshold);
	Mat tmp = dtm.clone();
	imwrite("Test1.png", dtm);
	int size = 30;

	uchar u, u1;
	bool bloop;
	for (int loop = 0; loop < 6; loop++)
	{
		bloop = false;
		for (int x = 0; x < dMap.cols; x++)
		{
			for (int y = 0; y < dMap.rows; y++)
			{
				int minv = 255;
				if (dtm.at<uchar>(y, x) > 1) continue;
				bloop = true;
				for (int k = -size; k <= size; k++)
				{
					if (y + k < 0 || y + k >= dMap.rows) continue;

					u = dtm.at<uchar>(y + k, x);
					if (u > 0 && minv > u) minv = u;
				}
				if (minv < 255) { dtm.at<uchar>(y, x) = minv; }
			}
		}

		for (int y = 0; y < dMap.rows; y++)
		{
			for (int x = 0; x < dMap.cols; x++)
			{
				if (dtm.at<uchar>(y, x) > 1) continue;
				bloop = true;
				int minv = 255;
				for (int k = -size; k <= size; k++)
				{
					if (x + k < 0 || x + k >= dMap.cols) continue;
					
					u = dtm.at<uchar>(y, x + k);
					if (u > 0 && minv > u) minv = u;
				}
				if (minv < 255) { dtm.at<uchar>(y, x) = minv; }
			}
		}
		if (bloop )continue;
		break;
	}
	Mat element = getStructuringElement(MORPH_RECT, Size(2 * 2 + 1, 2 * 2 + 1), Point(2, 2));
	dilate(dtm, dtm, element);
	imwrite("Test2.png", dtm);
	erode(dtm, dtm, element);
	element = getStructuringElement(MORPH_RECT, Size(2 * (size/2) + 1, 2 *(size/2) + 1), Point(size/2, size/2));

	erode(dtm, tmp, element); 
	dilate(tmp, dtm, element);
	element = getStructuringElement(MORPH_RECT, Size(2 * size + 1, 2 * size + 1), Point(size, size));
	dilate(dtm, tmp, element);
	dilate(tmp, dtm, element);
	erode(dtm, dtm, element);//dilate(tmp, dtm, element);

	// restore the final detail
	 /* 
	for (int x = 0; x < dMap.cols; x++)
	{
		for (int y = 0; y < dMap.rows; y++)
		{
			u = dtm.at<uchar>(y, x);
			u1 = dMap.at<uchar>(y, x);
			if (u1 <= 1) continue;
			if ( u1 < u ) dtm.at<uchar>(y, x) = u1; 
		}
	}
	// */
	imwrite("Test.png", dtm);

}

int detectShift(const Mat & img, const Mat & dMap, Mat & bar)
{
	Mat Mask = Mat(img.rows, img.cols, CV_8UC1);
	Vec3b p;// bgr
	for (int x =0; x < Mask.cols; x++)
	{
		for (int y =0; y < Mask.rows; y++)
		{
			p = img.at<Vec3b>(y, x);
			if (p[0] >= p[1] && p[0] >= p[2]) {
				Mask.at<uchar>(y, x) = 255;
			}
			else Mask.at<uchar>(y, x) = 0;
		}
	}

	Mat grad_x;
	int  erosion_size = 3;
	Mat element = getStructuringElement(MORPH_ELLIPSE,
		Size(2 * erosion_size + 1, 2 * erosion_size + 1),
		Point(erosion_size, erosion_size));
	erode(Mask, grad_x, element);
	dilate(grad_x, Mask, element);

	//imwrite("TestM.jpg", Mask);

	erosion_size = 1;
	element = getStructuringElement(MORPH_ELLIPSE,
		Size(2 * erosion_size + 1, 2 * erosion_size + 1),
		Point(erosion_size, erosion_size));
	dilate(Mask, bar, element);
	//imwrite("TestM1.jpg",bar);

	vector<unsigned int> hist1, hist2; hist1.resize(256); hist2.resize(256);
	unsigned int total = 0, coverage = 0;
	for (int i = 0; i < 256; i++)
	{
		hist1[i] = 0; hist2[i] = 0;
	}

	// histogram
	for (int x = 0; x < Mask.cols; x++)
	{
		for (int y = 0; y < Mask.rows; y++)
		{
			uchar u = Mask.at<uchar>(y, x);
			if (u == 0) continue;
			u= dMap.at<uchar>(y, x);
			if (u == 0) continue;
			if (x <= Mask.cols/2) 
				hist1[u]++;
			else 
				hist2[u]++;
		}
	}

	int ref1, ref2, w1=hist1[1], w2=hist2[1];
	for (int i = 2; i < 255; i++)
	{
		if (w1 < hist1[i]) {
			ref1 = i; w1 = hist1[i]; 
		}
		if (w2 < hist2[i]) {
			ref2 = i; w2 = hist2[i];
		}
	}

	if(w1/2 < w2)
		return ref2;
	else return ref1;
}


int heightFromDmap(string & dirctory, double & beta, vector<string> & v3Files, vector<string> & v4Files,
	const vector<int> & vHeight, const vector<Coordinates> vCoordinates) 
{
	float dFilter[] = { 0.15f, 0.2f, 0.3f, 0.2f, 0.15f };
	vector<double> betaVec;
	betaVec.push_back(beta);
	string dir3, dir4, fname3, fname4;
	string cName3, cName4, dName, dName1, dName2;
	int pos = 0, pos1 = 0;
	Mat imageA, imageB, IA, dMap;
	dir3 = dirctory + "\\3" + "\\";
	float scale = 0.35f;
	Coordinates corrd;
	string fileExt1 = ".JPG", fileExt2=".jpg";

	std::ofstream outHeightfile;
	outHeightfile.open("Height.csv", ios::app);
	for (int i = 0; i < dirctory.length(); i++)
	{
		char c = dirctory[i];
		if (c >= '1' && c <= '9')
		{
			pos = i; break;
		}
	}
	outHeightfile << "\n" << dirctory.substr(pos) << ",";
	cout << "\n" << dirctory.substr(pos) << "\n";

	vector<Point2f> objPoints, scenePoints;
	int totalNumImage = Min(v3Files.size(), v4Files.size());
	int xmin, xmax, ymin, ymax;
	for (int i = 0; i < totalNumImage; i++)
	{
		corrd = vCoordinates[i / 3];
		xmin = Max(corrd.points[0].x, corrd.points[1].x) *scale;
		xmax = Min(corrd.points[2].x, corrd.points[3].x) *scale;
		ymin = 100; ymax = corrd.points[1].y * scale - 100;
		cout << v3Files[i] << ", ";
		fname3 = dir3 + v3Files[i];

		pos = v3Files[i].find(fileExt1); if (pos < 0) pos = v3Files[i].find(fileExt2);
		cName3 = v3Files[i].substr(0, pos);

		dName = dirctory + "\\" + cName3 + "_dmap.jpg";
		dName1 = dirctory + "\\" + cName3 + "5Norm_dmap.jpg";
		dName2 = dirctory + "\\" + cName3 + "5Norm_dmap.csv";

		if (exists(fname3) && exists(dName))
		{
			IA = imread(fname3, 1);
			dMap = imread(dName, 0);
		}
		else
			continue;

		CvSize imSize = IA.size();

		double cx = IA.cols / 2.0;
		double cy = IA.rows / 2.0;
		distortionCorrection_Div(IA, imageA, betaVec, betaVec.size(), cx, cy);
		// /*
		Mat H1, H2, Homo;
		string xmlName = dirctory + "\\RectificationMatrix.xml";
		if (exists(xmlName))
			readRectificationMatrix(xmlName, H1, H2, Homo);
		else return -1;

		warpPerspective(imageA, IA, H1, imSize);
		// */

		if (scale != 1.f)
		{
			Mat temp1;
			int method = scale < 1 ? INTER_AREA : INTER_CUBIC;
			resize(IA, temp1, Size(), scale, scale, method);
			IA = temp1;
		}

		Size img_size = IA.size();

		Mat bar;
		int ref = detectShift(IA, dMap, bar);
		int diff = 120 * scale - ref;

		// histogram
		vector<unsigned int> hist; hist.resize(256);
		unsigned int total = 0, coverage = 0;
		for (int i = 0; i < 256; i++)
		{
			hist[i] = 0;
		}
		
		// /*
		for (int x = 0; x < IA.cols; x++)
		{
			for (int y = 0; y < IA.rows; y++)
			{
				uchar u = dMap.at<uchar>(y, x);
				int tmpD = diff + (int)(u);
				if (tmpD < 0) tmpD = 0;
				if (tmpD > 255) tmpD = 255;
				dMap.at<uchar>(y, x) = (uchar)(tmpD);
				if (x <= xmin || x >= xmax || y <= ymin || y >= ymax) 
					continue;

				u = (uchar)(tmpD);
				if (u == 0) continue;
				hist[u]++;
			}
		}
		// */
		imwrite(dName1, dMap);

		//histogram filetring
		pos = 0;
		for (int i = 1; i <= 10; i++)
		{
			if (i < 6) hist[i] = 0;
			if (hist[i] > hist[pos]) pos = i;
		}
		pos1 = 0;
		for (int i = 1; i <= 10; i++)
		{
			if (i!= pos && hist[i] > hist[pos1]) pos1 = i;
		}
		hist[pos] = hist[pos1];

		vector<float> hist1; hist1.resize(256);
		for (int i = 0; i < 256; i++)
		{
			hist1[i] = 0;
			for (int j = 0; j <= 4; j++)
			{
				if (i + j - 2 < 0 || i + j - 2 > 255) continue;
				hist1[i] += dFilter[j] * hist[i + j - 2];
			}
		}
		for (int i = 0; i < 256; i++)
		{
			hist[i] = 0;
			float t1 = 0;
			for (int j = 0; j <= 4; j++)
			{
				if (i + j - 2 < 0 || i + j - 2 > 255) continue;
				t1 += dFilter[j] * hist1[i + j - 2];
			}
			hist[i] = (int)t1;
		}

		std::ofstream outHistogram;
		outHistogram.open(dName2);
		for (int i = 0; i < 256; i++)
		{
			float t1 = i*2.23 * 190 / (5184 * scale);
			float t2 = t1 * 190;
			t1 = (int)(10 * t2 / (39.88 + t1));
			t1 /= 10;  // one decimal
			outHistogram << t1 << "," << hist[i] << "\n";
		}
		outHistogram.close();


		 /*
		for (int x = 0; x < IA.cols; x++)
		{
			for (int y = 0; y < IA.rows; y++)
			{
				if (x > xmin + 50 && x < xmax - 50) continue;
				if (bar.at<uchar>(y, x) >0) continue;
				uchar u = dMap.at<uchar>(y, x);
				int tmpD = diff + (int)(u);
				if (tmpD < 0) tmpD = 0;
				if (tmpD > 255) tmpD = 255;
				u = (uchar)(tmpD);
				if (u == 0) continue;
					hist[u]++;
			}
		}
		// */

		// average height
		float aveDisparity = 0, sum=0, sumw = 0;
		for (int i = 10; i < 250; i++)
		{
			sum += hist[i] * i;
			sumw += hist[i];
		}
		aveDisparity = sum / sumw;
		//cout << aveDispaity << "\n";

		// re-estimate 
		sum = 0; sumw = 0;
		int thre = Min(250, 2 * aveDisparity);
		//int thre = 220;
		for (int i = 10; i <thre; i++)
		{
			sumw += hist[i];
		}
		sumw = sumw*0.02; // 2 %
		
		// estimate disparity
		sum = 0;
		int dis = 0;
		for (int i = thre-1; i >0; i--)
		{
			sum += hist[i];
			if (sum >= sumw) { dis = i; break; }
		}

		float t1 = dis*2.23 * 190 / (5184 * scale);
		float t2 = t1 * 190;
		float height = (int) (100*t2/(39.88 + t1));
		height /= 100;  // two decimal

		if(ref >=180 * scale || ref < 100 * scale) outHeightfile << height << "?,"; // not reliable
		else outHeightfile << height << ",";
		
	}
	outHeightfile.close();
	return 0;
}


int Stereo_SGBM(string & dirctory, double & beta, vector<string> & v3Files, vector<string> & v4Files )
{
	cout << "\n" << dirctory << ":" << "\n";
	string fileExt1 = ".JPG", fileExt2 = ".jpg";
	enum { STEREO_BM = 0, STEREO_SGBM = 1, STEREO_HH = 2, STEREO_VAR = 3, STEREO_3WAY = 4 };
	int alg = STEREO_SGBM;
	int SADWindowSize, numberOfDisparities, minDisparities;
	float scale;
	//algorithm = sgbm --blocksize = 7 --min - disparity = -20 --max - disparity = 256 --scale = 1.0 - o = disp.png
	Ptr<StereoSGBM> sgbm = StereoSGBM::create(0, 16, 3);
	minDisparities = -20;
	numberOfDisparities = 256;
	SADWindowSize = 5;
	scale = 0.35f;

	if (alg < 0)
	{
		printf("Command-line parameter error: Unknown stereo algorithm\n\n");
		print_help();
		return -1;
	}
	if (numberOfDisparities < 1 || numberOfDisparities % 16 != 0)
	{
		printf("Command-line parameter error: The max disparity (--maxdisparity=<...>) must be a positive integer divisible by 16\n");
		print_help();
		return -1;
	}
	if (scale < 0)
	{
		printf("Command-line parameter error: The scale factor (--scale=<...>) must be a positive floating-point number\n");
		return -1;
	}
	if (SADWindowSize < 1 || SADWindowSize % 2 != 1)
	{
		printf("Command-line parameter error: The block size (--blocksize=<...>) must be a positive odd number\n");
		return -1;
	}

	vector<double> betaVec;
	betaVec.push_back(beta);
	string dir3, dir4, fname3, fname4;
	string cName3, cName4, dName;
	int pos=0;
	Mat imageA, imageB, IA, IB;
	dir3 = dirctory + "\\3" + "\\";
	dir4 = dirctory + "\\4" + "\\";

	CvSize ori_imSize;
	vector<Point2f> objPoints, scenePoints;
	int totalNumImage = Min(v3Files.size(), v4Files.size());
	for (int i = 0; i < totalNumImage; i++)
	{
		cout << v3Files[i] << ", ";
		fname3 = dir3 + v3Files[i];
		fname4 = dir4 + v4Files[i];
		pos = v3Files[i].find(fileExt1); if (pos < 0) pos = v3Files[i].find(fileExt2);
		cName3 = v3Files[i].substr(0, pos);

		dName = dirctory+ "\\"+ cName3 + "_dmap.jpg";

		cName3 = dir3 + cName3 + "_C.jpg";
		
		pos = v4Files[i].find(fileExt1); if (pos < 0) pos = v4Files[i].find(fileExt2);
		cName4 = v4Files[i].substr(0, pos);
		cName4 = dir4 + cName4 + "_C.png";

		if (exists(fname3) && exists(fname4))
		{
			IA = imread(fname3, 1);
			IB = imread(fname4, 1);
		}
		else
			continue;

		CvSize imSize = IA.size();
		ori_imSize = IA.size();
		if (imSize.width != IB.cols || imSize.height != IB.rows)
		{
			fprintf(stderr, "Error: image sizes are not same.\n");
			return -1;
		}

		double cx = IA.cols / 2.0;
		double cy = IA.rows / 2.0;
		distortionCorrection_Div(IA, imageA, betaVec, betaVec.size(), cx, cy);
		distortionCorrection_Div(IB, imageB, betaVec, betaVec.size(), cx, cy);
		if (imageA.cols == 0 || imageB.cols == 0)   return -1;


		if (!imageA.data || !imageB.data)
		{
			std::cout << " --(!) Error reading images " << std::endl; return -1;
		}

		Mat H1, H2, Homo;
		string xmlName = dirctory + "\\RectificationMatrix.xml";
		if (exists(xmlName) )
			readRectificationMatrix(xmlName, H1, H2, Homo);
		else return -1;

		warpPerspective(imageA, IA, H1, imSize);
		warpPerspective(imageB, IB, H2, imSize);

		imwrite(cName3, IA);
		//imwrite(cName4, IB);

		int color_mode = alg == STEREO_BM ? 0 : -1;
		Mat trans = Mat::eye(2, 3, CV_64F); trans.at<double>(0, 2) = minDisparities;
		Mat IC;
		warpAffine(IB, IC, trans, IB.size());

		if (scale != 1.f)
		{
			Mat temp1, temp2;
			int method = scale < 1 ? INTER_AREA : INTER_CUBIC;
			resize(IA, temp1, Size(), scale, scale, method);
			IA = temp1;
			resize(IC, temp2, Size(), scale, scale, method);
			IC = temp2;
		}

		Size img_size = IA.size();

		numberOfDisparities = numberOfDisparities > 0 ? numberOfDisparities : ((img_size.width / 8) + 15) & -16;

		sgbm->setPreFilterCap(63);
		int sgbmWinSize = SADWindowSize > 0 ? SADWindowSize : 3;
		sgbm->setBlockSize(sgbmWinSize);

		int cn = IA.channels();

		sgbm->setP1(8 * cn*sgbmWinSize*sgbmWinSize);
		sgbm->setP2(32 * cn*sgbmWinSize*sgbmWinSize);
		sgbm->setMinDisparity(0);
		sgbm->setNumDisparities(numberOfDisparities);
		sgbm->setUniquenessRatio(10);
		sgbm->setSpeckleWindowSize(100);
		sgbm->setSpeckleRange(32);
		sgbm->setDisp12MaxDiff(1);

		sgbm->setMode(StereoSGBM::MODE_SGBM);

		Mat disp, disp8;

		//int64 t = getTickCount();
		sgbm->compute(IA, IC, disp);
		//t = getTickCount() - t;
		//printf("Time elapsed: %fms\n", t * 1000 / getTickFrequency());

		if (alg != STEREO_VAR)
			disp.convertTo(disp8, CV_8U, 256 / (numberOfDisparities*16.));
		else
			disp.convertTo(disp8, CV_8U);

		imwrite(dName, disp8);
	}

	return 1;

}



void Stereo_SGBM(vector<string> & vBasDir, vector<double> & vBeta, vector<vector<string>> & v3ImFiles, vector<vector<string>> & v4ImFiles)
{
	for (int i = 0; i < vBasDir.size(); i++)
	{
		Stereo_SGBM(vBasDir[i], vBeta[i], v3ImFiles[i], v4ImFiles[i]);
	}
}

void heightFromDMap(vector<string> & vBasDir, vector<double> & vBeta, vector<vector<string>> & v3ImFiles, vector<vector<string>> & v4ImFiles, 
	vector<string> & sDateH, vector<vector<int>> & vHeight, vector<string> & sDateC, vector<vector<Coordinates>> & vCoordinates)
{
	for (int i = 0; i < vBasDir.size(); i++)
	{
		int hIdx = -1, cIdx=-1;
		for (int h = 0; h < sDateH.size(); h++)
		{
			int pos = vBasDir[i].find(sDateH[h]);
			if (pos >= 0) { hIdx = h; break; }
		}

		for (int c = 0; c < sDateC.size(); c++)
		{
			int pos = vBasDir[i].find(sDateC[c]);
			if (pos >= 0) { cIdx = c; break; }
		}
		if (hIdx < 0 || cIdx < 0) continue;
		heightFromDmap(vBasDir[i], vBeta[i], v3ImFiles[i], v4ImFiles[i], vHeight[hIdx], vCoordinates[cIdx]);
	}
}

int main(int argc, char** argv)
{
	//Mat dMap = imread("disp50.png", 0);
	//Mat dtm;
	//Mat org = imread("RectifiedA.jpg", 1);
	//DTM(org, dMap, dtm, 50);
	//DTM(dMap, dtm, 60);

	std::ofstream outHeightfile;
	outHeightfile.open("Height.csv");
	outHeightfile << ",";
	for (int i = 1; i <= 60; i++)
	{
		outHeightfile << i<<"a," << i << "b," << i << "c,";
	}
	outHeightfile.close();
	
	vector<string> sDateH;
	vector<vector<int>> vHeight;
	vector<string> sDateC;
	vector<vector<Coordinates>> vCoordinates;
	readHeight("plant_heights.csv",  sDateH,  vHeight);
	readMaskCoordinates("mask_coordinates.csv", sDateC, vCoordinates);

	char* infile = "file.txt"; // dir *.jpg /s/w

	int pos = 0;
	vector<string> vBasDir;
	vector<vector<string>> v3ImFiles;
	vector<vector<string>> v4ImFiles;
	vector<string> vTmpFiles;

	std::ifstream readfile(infile);
	string baseDir, line0, tmpStr, line1;
	string tokenExt = ".JPG", tokenDir = "Directory of", tokenR = "RIGHT", tokenL = "LEFT";
	string tokenCam1 = "3", tokenCam2 = "4", tokenSetup = "setup";
	bool bLeft = false;

	int lengthD = tokenDir.length();
	int camIdx = 0;
	if (readfile.is_open())
	{
		while (readfile.good())
		{
			getline(readfile, line0);
			if (line0.length() < 5) continue;
			//baseDir = line0;
			pos = line0.find(tokenDir);
			if (pos > 0)
			{
				if (vTmpFiles.size()>0) {
					if (camIdx == 3) v3ImFiles.push_back(vTmpFiles);
					if (camIdx == 4) v4ImFiles.push_back(vTmpFiles);
					vTmpFiles.clear();
				}

				line0 = line0.substr(lengthD + 1, line0.length() - lengthD);
				trim(line0);
				if (line0[line0.length() - 1] == '3') camIdx = 3;
				else if (line0[line0.length() - 1] == '4') camIdx = 4;
				else {
					camIdx = 0; continue;
				}

				if (baseDir.length() == 0)
				{
					baseDir = line0.substr(0, line0.length() - 2);
					vBasDir.push_back(baseDir);
				}
				else
				{
					pos = line0.find(baseDir);
					if (pos < 0)
					{
						baseDir = line0.substr(0, line0.length() - 2);
						vBasDir.push_back(baseDir);
					}
				}
				continue;
			}

			if (camIdx == 0) continue;

			// split a line
			std::string token;
			stringstream ss(line0);

			while (getline(ss, token, ' ')) {
				if (token.length() < 4) continue;
				string str = token;
				transform(str.begin(), str.end(), str.begin(), ::toupper);
				pos = str.find(tokenExt);
				if (pos >= 0) vTmpFiles.push_back(token);
			}
		}

	}

	vector<double> vBeta;
	for (int i = 0; i < vBasDir.size(); i++)
		vBeta.push_back(-4.613100000000000e-009);
	//vBeta[4] = vBeta[9] = -1.100000000000e-009;
	//vBeta[23] = -2.6067100000000000e-009;

	Stereo_SGBM(vBasDir, vBeta, v3ImFiles, v4ImFiles);
	
	heightFromDMap(vBasDir, vBeta, v3ImFiles, v4ImFiles, sDateH, vHeight, sDateC, vCoordinates);

	/*
	string fileNameCSV = infile;
	fileNameCSV = fileNameCSV.substr(0, fileNameCSV.length() - 4) + ".csv";

	vector<string> fileJpgList;
	if (readfile.is_open())
	{
	while (readfile.good())
	{
	getline(readfile, line0);
	trim(line0);
	if (line0 == "") continue;
	fileJpgList.push_back(line0);
	}
	}
	readfile.close();

	vector<double> betaVec(0);
	readBetavec("BetaVec.xml", betaVec);

	ofstream outputfile;
	outputfile.open(fileNameCSV);
	for (int i = 0; i < fileJpgList.size(); i++)
	{
	//double cov = coverag(betaVec, fileJpgList[i], baseDir, mask);
	//cov = cov * 100;
	//cout << i << "  ";
	//if (cov >3.0) outputfile << fileJpgList[i] << "," << cov << endl;
	}
	outputfile.close();
	*/

	return 0;
}

