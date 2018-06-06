
#include "utils.h"

#include <gsl/gsl_randist.h>

#include <time.h>

//// Local Function Prototypes ///
int calc_min_inliers( int, int, double, double );
//static __inline struct feature* get_match( struct feature*, int );
//int get_matched_features( struct feature*, int, int, struct feature*** );
static std::vector<int> draw_ransac_sample(int, int, gsl_rng* );
void extract_corresp_pts( std::vector<int> vSample, CvPoint2D32f *f1, CvPoint2D32f *f2, 
						 CvPoint2D32f**, CvPoint2D32f** );
std::vector<int> find_consensus( CvPoint2D32f *f1, CvPoint2D32f *f2, int n,
				   CvMat* M, ransac_err_fn err_fn, double err_tol);
static __inline void release_mem( CvPoint2D32f*, CvPoint2D32f*);

// */


//
//Calculates a best-fit image transform from image feature correspondences
//using RANSAC.
//
//For more information refer to:
//
//Fischler, M. A. and Bolles, R. C.  Random sample consensus: a paradigm for
//model fitting with applications to image analysis and automated cartography.
//<EM>Communications of the ACM, 24</EM>, 6 (1981), pp. 381--395.
//
//@param features an array of features; only features with a non-NULL match
//	of type mtype are used in homography computation
//@param n number of features in feat
//@param mtype determines which of each feature's match fields to use
//	for model computation; should be one of FEATURE_FWD_MATCH,
//	FEATURE_BCK_MATCH, or FEATURE_MDL_MATCH; if this is FEATURE_MDL_MATCH,
//	correspondences are assumed to be between a feature's img_pt field
//	and its match's mdl_pt field, otherwise correspondences are assumed to
//	be between the the feature's img_pt field and its match's img_pt field
//@param xform_fn pointer to the function used to compute the desired
//	transformation from feature correspondences
//@param m minimum number of correspondences necessary to instantiate the
//	model computed by xform_fn
//@param p_badxform desired probability that the final transformation
//	returned by RANSAC is corrupted by outliers (i.e. the probability that
//	no samples of all inliers were drawn)
//@param err_fn pointer to the function used to compute a measure of error
//	between putative correspondences and a computed model
//@param err_tol correspondences within this distance of a computed model are
//	considered as inliers
//@param inliers if not NULL, output as an array of pointers to the final
//	set of inliers
//@param n_in if not NULL and \a inliers is not NULL, output as the final
//	number of inliers
//
//@return Returns a transformation matrix computed using RANSAC or NULL
//	on error or if an acceptable transform could not be computed.
//ransac_xform (f1, f2, featNum, lsq_homog, 4, 0.01,
//			homog_xfer_err, 3.0, NULL, NULL );
CvMat* ransac_xform( CvPoint2D32f *f1, CvPoint2D32f *f2, int n,
					ransac_xform_fn xform_fn, int m, double p_badxform,
					ransac_err_fn err_fn, double err_tol)
{
	//struct feature** matched, ** sample, ** consensus, ** consensus_max = NULL;
	//struct ransac_data* rdata;
	std::vector<int> vSample, consensus, consensus_max;
	CvPoint2D32f* pts, * mpts;
	CvMat* M = NULL;
	int mtype = 1; //FEATURE_FWD_MATCH;
	
	gsl_rng* rng;
	double p, in_frac = RANSAC_INLIER_FRAC_EST;
	int i, in_min, in_max = 0, k = 0;

	if( n < m )
	{
		fprintf( stderr, "Warning: not enough matches to compute xform, %s" \
			" line %d\n", __FILE__, __LINE__ );
		return M;
	}

	// initialize random number generator 
	rng = gsl_rng_alloc( gsl_rng_mt19937 );
	gsl_rng_set( rng, time(NULL) );

	in_min = calc_min_inliers( n, m, RANSAC_PROB_BAD_SUPP, p_badxform );
	p = pow( 1.0 - pow( in_frac, m ), k );
	i = 0;
	while( p > p_badxform )
	{
		vSample = draw_ransac_sample(n, m, rng );
		extract_corresp_pts(vSample, f1, f2, &pts, &mpts );
		M = xform_fn( pts, mpts, m );

		if( M ) {
			consensus = find_consensus( f1, f2, n, M, err_fn, err_tol);
			if( consensus.size() > in_max )
			{
				consensus_max = consensus;
				in_max = consensus.size();
				in_frac = (double)in_max / n;
			}
		
			cvReleaseMat( &M );
		}

		release_mem( pts, mpts);
		p = pow( 1.0 - pow( in_frac, m ), ++k );
	}

	// calculate final transform based on best consensus set //
	if( in_max >= in_min )
	{
		extract_corresp_pts( consensus_max, f1, f2, &pts, &mpts );
		M = xform_fn( pts, mpts, in_max );
		consensus = find_consensus( f1, f2, n, M, err_fn, err_tol);
		cvReleaseMat( &M );
		release_mem( pts, mpts);
		extract_corresp_pts( consensus, f1, f2, &pts, &mpts );
		M = xform_fn( pts, mpts, consensus.size() );
		
		release_mem( pts, mpts);
	}

	gsl_rng_free( rng );

	return M;
}
// */


//
//Calculates a least-squares planar homography from point correspondeces.
//
//@param pts array of points
//@param mpts array of corresponding points; each pts[i], i=0..n-1, corresponds
//	to mpts[i]
//@param n number of points in both pts and mpts; must be at least 4
//
//@return Returns the 3 x 3 least-squares planar homography matrix that
//	transforms points in pts to their corresponding points in mpts or NULL if
//	fewer than 4 correspondences were provided
//
CvMat* lsq_homog( CvPoint2D32f* pts, CvPoint2D32f* mpts, int n )
{
	CvMat* H, * A, * B, X;
	float x[9];
	int i;

	if( n < 4 )
	{
		fprintf( stderr, "Warning: too few points in lsq_homog(), %s line %d\n",
			__FILE__, __LINE__ );
		return NULL;
	}

	// set up matrices so we can unstack homography into X; AX = B //
	A = cvCreateMat( 2*n, 8, CV_32FC1 );
	B = cvCreateMat( 2*n, 1, CV_32FC1 );
	X = cvMat( 8, 1, CV_32FC1, x );
	H = cvCreateMat(3, 3, CV_32FC1);
	cvZero( A );
	for( i = 0; i < n; i++ )
	{
		cvmSet( A, i, 0, pts[i].x );
		cvmSet( A, i+n, 3, pts[i].x );
		cvmSet( A, i, 1, pts[i].y );
		cvmSet( A, i+n, 4, pts[i].y );
		cvmSet( A, i, 2, 1.0 );
		cvmSet( A, i+n, 5, 1.0 );
		cvmSet( A, i, 6, -pts[i].x * mpts[i].x );
		cvmSet( A, i, 7, -pts[i].y * mpts[i].x );
		cvmSet( A, i+n, 6, -pts[i].x * mpts[i].y );
		cvmSet( A, i+n, 7, -pts[i].y * mpts[i].y );
		cvmSet( B, i, 0, mpts[i].x );
		cvmSet( B, i+n, 0, mpts[i].y );
	}
	cvSolve( A, B, &X, CV_SVD );
	x[8] = 1.0;
	X = cvMat( 3, 3, CV_32FC1, x );
	cvConvert( &X, H );

	cvReleaseMat( &A );
	cvReleaseMat( &B );
	return H;
}
// */


//
//Calculates the transfer error between a point and its correspondence for
//a given homography, i.e. for a point x, it's correspondence x', and
//homography H, computes d(x', Hx)^2.
//
//@param pt a point
//@param mpt pt's correspondence
//@param H a homography matrix
//
//@return Returns the transfer error between pt and mpt given H
//
double homog_xfer_err( CvPoint2D32f pt, CvPoint2D32f mpt, CvMat* H )
{
	CvPoint2D32f xpt = persp_xform_pt( pt, H );

	return sqrt( (xpt.x-mpt.x)*(xpt.x-mpt.x)+(xpt.y-mpt.y)*(xpt.y-mpt.y) );
}



//
//Performs a perspective transformation on a single point.  That is, for a
//point (x, y) and a 3 x 3 matrix T this function returns the point
//(u, v), where
//
//[x' y' w']^T = T * [x y 1]^T,
//
//and
//
//(u, v) = (x'/w', y'/w').
//
//Note that affine transforms are a subset of perspective transforms.
//
//@param pt a 2D point
//@param T a perspective transformation matrix
//
//@return Returns the point (u, v) as above.
//
CvPoint2D32f persp_xform_pt( CvPoint2D32f pt, CvMat* T )
{
	CvMat XY, UV;
	float xy[3] = { pt.x, pt.y, 1.0 }, uv[3] = { 0 };

	cvInitMatHeader( &XY, 3, 1, CV_32FC1, xy, CV_AUTOSTEP );
	cvInitMatHeader( &UV, 3, 1, CV_32FC1, uv, CV_AUTOSTEP );
	cvMatMul( T, &XY, &UV );
	CvPoint2D32f rslt = cvPoint2D32f( uv[0] / uv[2], uv[1] / uv[2] );

	return rslt;
}



/*
//Returns a feature's match according to a specified match type
//
//@param feat feature
//@param mtype match type, one of FEATURE_FWD_MATCH, FEATURE_BCK_MATCH, or
//FEATURE_MDL_MATCH
//
//@return Returns feat's match corresponding to mtype or NULL for bad mtype
//
static __inline struct feature* get_match( struct feature* feat, int mtype )
{
	if( mtype == FEATURE_MDL_MATCH )
		return feat->mdl_match;
	if( mtype == FEATURE_BCK_MATCH )
		return feat->bck_match;
	if( mtype == FEATURE_FWD_MATCH )
		return feat->fwd_match;
	return NULL;
}



//
//Finds all features with a match of a specified type and stores pointers
//to them in an array.  Additionally initializes each matched feature's
//feature_data field with a ransac_data structure.
//
//@param features array of features
//@param n number of features in features
//@param mtype match type, one of FEATURE_{FWD,BCK,MDL}_MATCH
//@param matched output as an array of pointers to features with a match of
//the specified type
//
//@return Returns the number of features output in matched.
//
int get_matched_features( struct feature* features, int n, int mtype,
struct feature*** matched )
{
	struct feature** _matched;
	struct ransac_data* rdata;
	int i, m = 0;

	_matched = calloc( n, sizeof( struct feature* ) );
	for( i = 0; i < n; i++ )
		if( get_match( features + i, mtype ) )
		{
			rdata = malloc( sizeof( struct ransac_data ) );
			memset( rdata, 0, sizeof( struct ransac_data ) );
			rdata->orig_feat_data = features[i].feature_data;
			_matched[m] = features + i;
			_matched[m]->feature_data = rdata;
			m++;
		}
		*matched = _matched;
		return m;
}

// */

//
//Calculates the minimum number of inliers as a function of the number of
//putative correspondences.  Based on equation (7) in
//
//Chum, O. and Matas, J.  Matching with PROSAC -- Progressive Sample Consensus.
//In <EM>Conference on Computer Vision and Pattern Recognition (CVPR)</EM>,
//(2005), pp. 220--226.
//
//@param n number of putative correspondences
//@param m min number of correspondences to compute the model in question
//@param p_badsupp prob. that a bad model is supported by a correspondence
//@param p_badxform desired prob. that the final transformation returned is bad
//
//@return Returns the minimum number of inliers required to guarantee, based
//	on p_badsupp, that the probability that the final transformation returned
//	by RANSAC is less than p_badxform
//
int calc_min_inliers( int n, int m, double p_badsupp, double p_badxform )
{
	double pi, sum;
	int i, j;

	for( j = m+1; j <= n; j++ )
	{
		sum = 0;
		for( i = j; i <= n; i++ )
		{
			pi = ( i - m ) * log( p_badsupp ) + ( n - i + m ) * log( 1.0 - p_badsupp ) +
				gsl_sf_lnchoose( n - m, i - m );
			sum += exp( pi );
		}
		if( sum < p_badxform )
			break;
	}
	return j;
}



//
//Draws a RANSAC sample from a set of features.
//
//@param features array of pointers to features from which to sample
//@param n number of features in features
//@param m size of the sample
//@param rng random number generator used to sample
//
//@return Returns an array of pointers to the sampled features; the sampled
//	field of each sampled feature's ransac_data is set to 1
//
struct std::vector<int> draw_ransac_sample(int n, int m, gsl_rng* rng )
{
	std::vector<int> vSample;
	int i, x;

	for( i = 0; i < m; i++ )
	{
		bool found = false;
		do
		{
			found = false;
			x = gsl_rng_uniform_int( rng, n );
			for(int k =0; k<vSample.size(); k++)
			{
				if(vSample[k] == x) {found = true; break;}
			}
			if (!found) vSample.push_back(x);
		}
		while(found);
	}

	return vSample;
}




//Extrancs raw point correspondence locations from a set of features
//
//@param features array of features from which to extract points and match
//	points; each of these is assumed to have a match of type mtype
//@param n number of features
//@param mtype match type; if FEATURE_MDL_MATCH correspondences are assumed
//	to be between each feature's img_pt field and it's match's mdl_pt field,
//	otherwise, correspondences are assumed to be between img_pt and img_pt
//@param pts output as an array of raw point locations from features
//@param mpts output as an array of raw point locations from features' matches
//
void extract_corresp_pts( std::vector<int> vSample, CvPoint2D32f *f1, CvPoint2D32f *f2,
						 CvPoint2D32f** pts, CvPoint2D32f** mpts )
{
	CvPoint2D32f* _pts, * _mpts;
	int n = vSample.size();

	_pts = (CvPoint2D32f*) calloc( n, sizeof( CvPoint2D32f ) );
	_mpts = (CvPoint2D32f*) calloc( n, sizeof( CvPoint2D32f ) );

	for( int i = 0; i < n; i++ )
	{
		_pts[i] = f1[vSample[i]];
		_mpts[i] = f2[vSample[i]];
	}

	*pts = _pts;
	*mpts = _mpts;
}



//
//For a given model and error function, finds a consensus from a set of
//feature correspondences.
//
//@param features set of pointers to features; every feature is assumed to
//	have a match of type mtype
//@param n number of features in features
//@param mtype determines the match field of each feature against which to
//	measure error; if this is FEATURE_MDL_MATCH, correspondences are assumed
//	to be between the feature's img_pt field and the match's mdl_pt field;
//	otherwise matches are assumed to be between img_pt and img_pt
//@param M model for which a consensus set is being found
//@param err_fn error function used to measure distance from M
//@param err_tol correspondences within this distance of M are added to the
//	consensus set
//@param consensus output as an array of pointers to features in the
//	consensus set
//
//@return Returns the number of points in the consensus set
//
std::vector<int> find_consensus( CvPoint2D32f *f1, CvPoint2D32f *f2, int n,
				   CvMat* M, ransac_err_fn err_fn, double err_tol)
{
	CvPoint2D32f pt, mpt;
	double err;
	int i, in = 0;

	std::vector<int> consensus;

	for( i = 0; i < n; i++ )
	{
		pt = f1[i]; mpt = f2[i];
		err = err_fn( pt, mpt, M );
		if( err <= err_tol )
			consensus.push_back(i);
	}

	return consensus;
}



//
//Releases memory and reduces code size above
//
//@param pts1 an array of points
//@param pts2 an array of points
//@param features an array of pointers to features; can be NULL
//
static __inline void release_mem( CvPoint2D32f* pts1, CvPoint2D32f* pts2)
{
	free( pts1 );
	free( pts2 );
}


static CvMat* CaiXTrans(CvPoint2D32f *f1, CvPoint2D32f *f2, int featNum, IplImage* img1, IplImage* img2)
{
		CvMat* H;
		H = ransac_xform( f1, f2, featNum, lsq_homog, 4, 0.01,
			homog_xfer_err, 3.0);
		if( H )
		{
			IplImage* xformed;
			xformed = cvCreateImage( cvGetSize( img2 ), IPL_DEPTH_8U, 3 );
			cvWarpPerspective( img1, xformed, H, 
				CV_INTER_LINEAR + CV_WARP_FILL_OUTLIERS,
				cvScalarAll( 0 ) );
			cvNamedWindow( "Xformed", 1 );
			/*cvSmooth(xformed, tmp, CV_GAUSSIAN, 9, 9, 0, 0);
			cvConvert(tmp, xformed);
			cvSmooth(img2, tmp, CV_GAUSSIAN, 9, 9, 0,0);
			cvConvert(tmp, img2);
			cvSub(xformed, img2, img3, 0);
			cvSaveImage("out0.jpg", img3);
			cvErode(img3, img3,0,2);
			cvShowImage( "Matches", img3 );
			cvSub(img2, xformed, tmp, 0);
			cvErode(tmp, tmp,0,2);
			cvShowImage( "Xformed", tmp );
			cvSaveImage("Xform.jpg", xformed);
			cvSaveImage("out1.jpg", tmp);
			
			cvShowImage( "Xformed", xformed );
			cvSaveImage("Xform.jpg", xformed);
			cvScale(img2, tmp, 0.5,0);
			cvScale(xformed, img3, 0.5,0);
			cvAdd(img3, tmp, xformed, 0);
			cvSaveImage("out1.jpg", xformed);

			cvWaitKey( 0 );
			*/
			cvReleaseImage( &xformed );
			cvReleaseMat( &H );
		}
	
	// 

	//cvReleaseImage( &img3 );
	//cvReleaseImage( &img1 );
	//cvReleaseImage( &img2 );
	//cvReleaseImage( &tmp );
	//kdtree_release( kd_root );
	//free( feat1 );
	//free( feat2 );
}