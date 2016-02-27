#pragma once
#include "atlimage.h"
#include "ImageMosaic.h"
#include "atltypes.h"
#include <list>
using namespace std;
typedef struct
{
	unsigned char M[37];
	float dist;
	short i;
	short j;  //角点坐标
}_CORNER;

typedef struct
{
	float D0;   //从待配准到基准的对应点间的欧式距离的平方
	float D1;   //从基准到待配准的对应点间的欧式距离的平方
	bool inliers; // 内点为1，外点为0
	short i;
	short j;
	short m;
	short n;
}_TEMP;
typedef struct
{
	float data;
	bool flag;
	short i;
	short j;
	short m;
	short n;
}_PAIR;
typedef struct
{
	bool flag;
	short i;
	short j;
	short m;
	short n;
}_PAIR2;
typedef struct
{
	float data;
	unsigned char M[37];
	short i;
	short j;  //角点坐标
}_DOT;
class CHarrisCorners
{
public:
	CHarrisCorners(void);
	~CHarrisCorners(void);
	void cornersDetecting(float* &Ii, float* &Ji,short height1,short width1,short height2, short width2, bool*& corner1, bool*& corner2);
	void harrisCornerDetecting(short height, short width, float*& J , bool*& corner ,bool*& IsEdge);
	void harrisCornerDetecting(short height, short width, float*& J , bool*& _corner,bool*& IsEdge,short minx,short maxx,short miny,short maxy);

	bool cornersMatching(float*& Ii,float*& Ji,bool*& corner1, bool*& corner2);
	bool cornersMatching(float*& Ii,float*& Ji,bool*& corner1, bool*& corner2, short m_nHeight1,short m_nWidth1,short m_nHeight2, short m_nWidth2,bool m_bTwice, bool m_bOnce, bool m_bOrigin);

	void coarseMatching(float*& Ig, float*& Jg, short height1, short width1, short height2, short width2, bool*& corner1, bool*& corner2,float dbMaxR);
	void coarseMatching(float*& Ig, float*& Jg, short height1, short width1, short height2, short width2, bool*& corner1, bool*& corner2,
							_PAIR2*& RansacPoints,int num, float dbMax, int nRadius,int nWvltNum);

	bool ransacMatching(_TEMP*& _temp, int n,int &maxN, short H1, short W1, short H2, short W2);
	bool ransacMatching(_TEMP*& _temp, int n1, int maxN,_PAIR2* pair, int n2,  list<_PAIR2> &listofRansac, short H1, short W1, short H2, short W2);
	bool judgeFourPoints(_PAIR2* temp, short H1, short W1, short H2, short W2, short n);
	void getProjectionPara(_PAIR2* temp, float* pDbProjectionPara0,bool flag);
	void calMatProduct(float* pDbSrc1, float* pDbSrc2, float* pDbDest,short y, short x, short z);
	BOOL calInvMatrix(float* pDbSrc, short nLen);
	void L_MOptimize(_PAIR2* temp, int temp_N,bool flag); 
	void imageFusion(CImage &pNewImage1,CImage &pNewImage2,CImage &m_newImage);
	void getAftProDim(short *left,short *right,short *top,short *bottom,short *newHeight,short *newWidth,short height1, short width1, short height2, short width2);
	short nHeight1;
	short nWidth1;
	short nHeight2;
	short nWidth2;

};
