#include "StdAfx.h"
#include "HarrisCorners.h"
#include "function.h"
#include "stdlib.h"
#include "time.h"
#include <list>
using namespace std;

float H[9];
float InvH[9];
_PAIR2* RansacPoints3;
_PAIR2* RansacPoints2;
_PAIR2* RansacPoints1;
_PAIR2* RansacPoints0;
int m_nNum3;
int m_nNum2;
int m_nNum1;
int m_nNum0;

list<_PAIR> m_listTribl;
list<_PAIR2> m_listTwice;
list<_PAIR2> m_listOnce;
list<_PAIR2> m_listOrigin;
CHarrisCorners::CHarrisCorners(void)
{
}

CHarrisCorners::~CHarrisCorners(void)
{
}
void CHarrisCorners::cornersDetecting(float* &Ii, float* &Ji,short height1,short width1,short height2, short width2, bool*& corner1, bool*& corner2)
{
	nHeight1 = height1;
	nWidth1  = width1;
	nHeight2 = height2;
	nWidth2  = width2;
	/******************************************对基准图像的角点检测*****************************************/
	bool* IsEdge1 = new bool[nHeight1 * nWidth1];
	memset(IsEdge1,0,nHeight1 * nWidth1 * sizeof(bool));
	prescreening(nHeight1 ,nWidth1 ,Ii ,IsEdge1);
	harrisCornerDetecting(nHeight1,nWidth1,Ii,corner1,IsEdge1);
	delete[]IsEdge1;
	IsEdge1 = NULL;
	/*******************************************对待配准图像的角点检测*****************************************/
	bool* IsEdge2 = new bool[nHeight2 * nWidth2];
	memset(IsEdge2,0,nHeight2 * nWidth2 * sizeof(bool));
	prescreening(nHeight2 ,nWidth2 ,Ji ,IsEdge2);
	harrisCornerDetecting(nHeight2,nWidth2,Ji,corner2 ,IsEdge2);
	delete[]IsEdge2;
	IsEdge2 = NULL;
}

void CHarrisCorners::harrisCornerDetecting(short height, short width, float*& J , bool*& _corner,bool*& IsEdge)
{

	#define     cim(ROW,COL)      cim[width * (COL) + (ROW)]
	#define _corner(ROW,COL)  _corner[width * (COL) + (ROW)]
 

	//--------------------------------------------------------------------------
	//                     第一步：利用差分算子对图像进行滤波
	//--------------------------------------------------------------------------
	
	//定义水平方向差分算子并求Jx
	float dx[25] = {-2,-1,0,1,2,
					-2,-1,0,1,2,
					-2,-1,0,1,2,
					-2,-1,0,1,2,
					-2,-1,0,1,2};

	//定义垂直方向差分算子并求Jy
	float dy[25] = {-2,-2,-2,-2,-2,
					-1,-1,-1,-1,-1,
					 0, 0, 0, 0, 0,
					 1, 1, 1, 1, 1,
					 2, 2, 2, 2, 2};


	float *cim = new float[width * height];
	memset(cim,0,width*height*sizeof(float));

	//使用5×5的高斯模板
	//计算模板参数
	float sigma = 0.8f;
	short gausswidth = 5;
	float h[5 * 5];
	for(short i = 0;i < gausswidth;i++)
	{
		for(short j = 0;j < gausswidth;j++)
		{
			h[i * gausswidth + j] = exp(-1 * ((i - short(gausswidth / 2)) * (i - short(gausswidth / 2)) + (j - short(gausswidth / 2)) * (j - short(gausswidth / 2))) / (2 * sigma * sigma));
		}
	}
			//归一化：使模板参数之和为1
	float total = 0;
	for(short i = 0;i < gausswidth * gausswidth;i++)
	{
			total += h[i];
	}
	for(short i = 0;i < gausswidth;i++)
	{
		for(short j = 0;j < gausswidth;j++)
		{
			h[i * gausswidth + j] /= total;
		}
	}
	//计算cim,角点量
	for(short j = 5; j < height - 5; j++)
	{
		for(short i = 5; i < width - 5; i++)
		{
			if(IsEdge[j * width + i] == 1)
			{
				float Jx2[25];
				float Jy2[25];
				float Jxy[25];
				for(int n = 0;n < 5; n++)
				{
					for(int m = 0; m < 5; m++)
					{
						float Jx = mbys(J,width,height,i+m-2,j+n-2,dx,5,5);
						float Jy = mbys(J,width,height,i+m-2,j+n-2,dy,5,5);
						Jx2[n*5+m] = Jx * Jx;
						Jy2[n*5+m] = Jy * Jy;
						Jxy[n*5+m] = Jx * Jy;
					}
				}
				float _Jx2 = mbys(Jx2,5,5,h,gausswidth,gausswidth);
				float _Jy2 = mbys(Jy2,5,5,h,gausswidth,gausswidth);
				float _Jxy = mbys(Jxy,5,5,h,gausswidth,gausswidth);
				float dd = ((_Jx2 * _Jy2 - _Jxy * _Jxy) / (_Jx2 + _Jy2 + 0.000001f)) * (2 * atan(sigma) / PI);
				cim(i,j) += dd;
			}

		}
	}

	float SUM[4] = {0,0,0,0};
	float EV[4] = {0,0,0,0};
	int sn = 0;

	for(short j = 0; j < height / 2; j++)
	{
		for(short i = 0; i < width / 2; i++)
		{
			if(cim(i,j) > 0&& IsEdge[j * width + i] == 1)
			{
				sn++;
				SUM[0] += cim(i,j);
			}
		}
	}
	EV[0] = SUM[0] / sn;
	for(short j = 0; j < height / 2; j++)
	{
		for(short i = 0; i < width / 2; i++)
		{
			if(cim(i,j) > EV[0] && IsEdge[j * width + i] == 1)
			{
				_corner(i,j) = 1;
			}
			else
			{
				_corner(i,j) = 0;
			}
		}
	}
	sn = 0;
	for(short j = 0; j < height / 2; j++)
	{
		for(short i = width / 2; i < width; i++)
		{
			if(cim(i,j) > 0 && IsEdge[j * width + i] == 1)
			{
				sn++;
				SUM[1] += cim(i,j);
			}
		}
	}
	EV[1] = SUM[1] / sn;
	for(short j = 0; j < height / 2; j++)
	{
		for(short i = width / 2; i < width; i++)
		{
			if(cim(i,j) > EV[1] && IsEdge[j * width + i] == 1)
			{
				_corner(i,j) = 1;
			}
			else
			{
				_corner(i,j) = 0;
			}
		}
	}
	sn = 0;
	for(short j = height / 2; j < height; j++)
	{
		for(short i = 0; i < width / 2; i++)
		{
			if(cim(i,j) > 0 && IsEdge[j * width + i] == 1)
			{
				sn++;
				SUM[2] += cim(i,j);
			}
		}
	}
	EV[2] = SUM[2] / sn;
	for(short j = height / 2; j < height; j++)
	{
		for(short i = 0; i < width / 2; i++)
		{
			if(cim(i,j) > EV[2] && IsEdge[j * width + i] == 1)
			{
				_corner(i,j) = 1;
			}
			else
			{
				_corner(i,j) = 0;
			}
		}
	}
	sn = 0;
	for(short j = height / 2; j < height; j++)
	{
		for(short i = width / 2; i < width; i++)
		{
			if(cim(i,j) > 0 && IsEdge[j * width + i] == 1)
			{
				sn++;
				SUM[3] += cim(i,j);
			}
		}
	}
	EV[3] = SUM[3] / sn;
	for(short j = height / 2; j < height; j++)
	{
		for(short i = width / 2; i < width; i++)
		{
			if(cim(i,j) > EV[3] && IsEdge[j * width + i] == 1)
			{
				_corner(i,j) = 1;
			}
			else
			{
				_corner(i,j) = 0;
			}
		}
	}
	//进行局部非极大值抑制以获得角点
	for(short j = 5; j < height - 5; j++)
	{
		for(short i = 5; i < width - 5; i++)
		{
			short size = 7;
			float max = 0.0f;
			if(j>short(size/2) && j<height-short(size/2) && i>short(size/2) && i<width-short(size/2)&&_corner(i,j) == 1)
			{
				for(short n=0;n<size;n++)
				{
					for(short m=0;m<size;m++)
					{						
						if(cim(i+m-short(size/2),j+n-short(size/2))>max)
						{
							max=cim(i+m-short(size/2),j+n-short(size/2));
						}
					}
				}
				for(short n=0;n<size;n++)
				{
					for(short m=0;m<size;m++)
					{						
						if(cim(i+m-short(size/2),j+n-short(size/2)) < max)
						{
							_corner(i+m-short(size/2),j+n-short(size/2)) = 0;
						}
					}
				}

			}
		}
	}
	delete[]cim;
	cim = NULL;	
}
//金字塔中间层图像角点检测
void CHarrisCorners::harrisCornerDetecting(short height, short width, float*& J , bool*& _corner,bool*& IsEdge,short minx,short maxx,short miny,short maxy)
{

	#define     cim(ROW,COL)      cim[width * (COL) + (ROW)]
	#define _corner(ROW,COL)  _corner[width * (COL) + (ROW)]
 

	//--------------------------------------------------------------------------
	//                     第一步：利用差分算子对图像进行滤波
	//--------------------------------------------------------------------------
	
	//定义水平方向差分算子并求Jx
	float dx[25] = {-2,-1,0,1,2,
					-2,-1,0,1,2,
					-2,-1,0,1,2,
					-2,-1,0,1,2,
					-2,-1,0,1,2};

	//定义垂直方向差分算子并求Jy
	float dy[25] = {-2,-2,-2,-2,-2,
					-1,-1,-1,-1,-1,
					 0, 0, 0, 0, 0,
					 1, 1, 1, 1, 1,
					 2, 2, 2, 2, 2};


	float *cim = new float[width * height];
	memset(cim,0,width*height*sizeof(float));
	
	//使用5×5的高斯模板
	//计算模板参数
	float sigma = float(0.8);
	short gausswidth = 5;
	float h[5 * 5];
	for(short i = 0;i < gausswidth;i++)
	{
		for(short j = 0;j < gausswidth;j++)
		{
			h[i * gausswidth + j] = exp(-1 * ((i - short(gausswidth / 2)) * (i - short(gausswidth / 2)) + (j - short(gausswidth / 2)) * (j - short(gausswidth / 2))) / (2 * sigma * sigma));
		}
	}
	//归一化：使模板参数之和为1
	float total = 0;
	for(short i = 0;i < gausswidth * gausswidth;i++)
	{
			total += h[i];
	}
	for(short i = 0;i < gausswidth;i++)
	{
		for(short j = 0;j < gausswidth;j++)
		{
			h[i * gausswidth + j] /= total;
		}
	}
	//计算cim,角点量
	#define max1(a,b)            (((a) > (b)) ? (a) : (b))
	#define min1(a,b)            (((a) < (b)) ? (a) : (b))
	for(short j = max1(miny-40,5); j < min1(maxy+40,height-5); j++)
	{
		for(short i = max1(minx-40,5); i < min1(maxx+40,width-5); i++)
		{
			if(IsEdge[j * width + i] == 1)
			{
				float Jx2[25];
				float Jy2[25];
				float Jxy[25];
				for(int n = 0;n < 5; n++)
				{
					for(int m = 0; m < 5; m++)
					{
						float Jx = mbys(J,width,height,i+m-2,j+n-2,dx,5,5);
						float Jy = mbys(J,width,height,i+m-2,j+n-2,dy,5,5);
						Jx2[n*5+m] = Jx * Jx;
						Jy2[n*5+m] = Jy * Jy;
						Jxy[n*5+m] = Jx * Jy;
					}
				}
				float _Jx2 = mbys(Jx2,5,5,h,gausswidth,gausswidth);
				float _Jy2 = mbys(Jy2,5,5,h,gausswidth,gausswidth);
				float _Jxy = mbys(Jxy,5,5,h,gausswidth,gausswidth);
				float dd = ((_Jx2 * _Jy2 - _Jxy * _Jxy) / (_Jx2 + _Jy2 + 0.000001f)) * (2 * atan(sigma) / PI);
				cim(i,j) += dd;
			}
		}
	}


	float SUM[4] = {0,0,0,0};
	float EV[4] = {0,0,0,0};
	int sn = 0;

	for(short j = 0; j < height / 2; j++)
	{
		for(short i = 0; i < width / 2; i++)
		{
			if(cim(i,j) > 0&& IsEdge[j * width + i] == 1)
			{
				sn++;
				SUM[0] += cim(i,j);
			}
		}
	}
	EV[0] = SUM[0] / sn;
	for(short j = 0; j < height / 2; j++)
	{
		for(short i = 0; i < width / 2; i++)
		{
			if(cim(i,j) > EV[0] && IsEdge[j * width + i] == 1)
			{
				_corner(i,j) = 1;
			}
			else
			{
				_corner(i,j) = 0;
			}
		}
	}
	sn = 0;
	for(short j = 0; j < height / 2; j++)
	{
		for(short i = width / 2; i < width; i++)
		{
			if(cim(i,j) > 0 && IsEdge[j * width + i] == 1)
			{
				sn++;
				SUM[1] += cim(i,j);
			}
		}
	}
	EV[1] = SUM[1] / sn;
	for(short j = 0; j < height / 2; j++)
	{
		for(short i = width / 2; i < width; i++)
		{
			if(cim(i,j) > EV[1] && IsEdge[j * width + i] == 1)
			{
				_corner(i,j) = 1;
			}
			else
			{
				_corner(i,j) = 0;
			}
		}
	}
	sn = 0;
	for(short j = height / 2; j < height; j++)
	{
		for(short i = 0; i < width / 2; i++)
		{
			if(cim(i,j) > 0 && IsEdge[j * width + i] == 1)
			{
				sn++;
				SUM[2] += cim(i,j);
			}
		}
	}
	EV[2] = SUM[2] / sn;
	for(short j = height / 2; j < height; j++)
	{
		for(short i = 0; i < width / 2; i++)
		{
			if(cim(i,j) > EV[2] && IsEdge[j * width + i] == 1)
			{
				_corner(i,j) = 1;
			}
			else
			{
				_corner(i,j) = 0;
			}
		}
	}
	sn = 0;
	for(short j = height / 2; j < height; j++)
	{
		for(short i = width / 2; i < width; i++)
		{
			if(cim(i,j) > 0 && IsEdge[j * width + i] == 1)
			{
				sn++;
				SUM[3] += cim(i,j);
			}
		}
	}
	EV[3] = SUM[3] / sn;
	for(short j = height / 2; j < height; j++)
	{
		for(short i = width / 2; i < width; i++)
		{
			if(cim(i,j) > EV[3] && IsEdge[j * width + i] == 1)
			{
				_corner(i,j) = 1;
			}
			else
			{
				_corner(i,j) = 0;
			}
		}
	}
	//进行局部非极大值抑制以获得角点
	for(short j = 5; j < height - 5; j++)
	{
		for(short i = 5; i < width - 5; i++)
		{
			short size = 7;
			float max = 0.0f;
			if(j>short(size/2) && j<height-short(size/2) && i>short(size/2) && i<width-short(size/2)&&_corner(i,j) == 1)
			{
				for(short n=0;n<size;n++)
				{
					for(short m=0;m<size;m++)
					{						
						if(cim(i+m-short(size/2),j+n-short(size/2))>max)
						{
							max=cim(i+m-short(size/2),j+n-short(size/2));
						}
					}
				}
				for(short n=0;n<size;n++)
				{
					for(short m=0;m<size;m++)
					{						
						if(cim(i+m-short(size/2),j+n-short(size/2)) < max)
						{
							_corner(i+m-short(size/2),j+n-short(size/2)) = 0;
						}
					}
				}

			}
		}
	}
	delete[]cim;
	cim = NULL;	
}
//不需要进行金字塔匹配时的原始图像匹配，或者是金字塔匹配时底层图像的匹配
bool CHarrisCorners::cornersMatching(float*& Ii,float*& Ji,bool*& corner1, bool*& corner2)
{

	//进行粗匹配
	coarseMatching(Ii,Ji,nHeight1,nWidth1,nHeight2,nWidth2,corner1,corner2,0.8f);
	int MatchingPoints_N = -1;
	_TEMP* MatchingPoints = new _TEMP[m_listTribl.size()];
	memset(MatchingPoints,0,m_listTribl.size() * sizeof(_TEMP));
	list<_PAIR>::iterator iter;
	for(iter = m_listTribl.begin();iter != m_listTribl.end();iter++)
	{
		MatchingPoints_N++;
		MatchingPoints[MatchingPoints_N].inliers = 0;
		MatchingPoints[MatchingPoints_N].i = iter->i;
		MatchingPoints[MatchingPoints_N].j = iter->j;
		MatchingPoints[MatchingPoints_N].m = iter->m;
		MatchingPoints[MatchingPoints_N].n = iter->n;
	}
	m_listTribl.clear();
	delete[]Ii;
	delete[]Ji;
	Ii = NULL;
	Ji = NULL;
	//进行细匹配
	int max_N = 0;
	bool flag = ransacMatching(MatchingPoints,MatchingPoints_N,max_N,nHeight1,nWidth1,nHeight2,nWidth2);//求H
	if(flag)
	{
		delete[]MatchingPoints;
		MatchingPoints = NULL;
		return 1;
	}
	int pair1_N = 0;
	for(int a = 0; a <= MatchingPoints_N; a++)
	{
		if(MatchingPoints[a].inliers == 1)
		{
			pair1_N++;
		}
	}
	_PAIR2 *pair1 = new _PAIR2[pair1_N];
	memset(pair1,0,pair1_N * sizeof(_PAIR2));//保存细匹配点对
	pair1_N = -1;
	for(int a = 0; a <= MatchingPoints_N; a++)
	{
		if(MatchingPoints[a].inliers == 1)
		{
			pair1_N++;
			pair1[pair1_N].i = MatchingPoints[a].i;
			pair1[pair1_N].j = MatchingPoints[a].j;
			pair1[pair1_N].m = MatchingPoints[a].m;
			pair1[pair1_N].n = MatchingPoints[a].n;
		}
	}
	list<_PAIR2> listofRansac;
	flag = ransacMatching(MatchingPoints,MatchingPoints_N,max_N,pair1,pair1_N,listofRansac,nHeight1,nWidth1,nHeight2,nWidth2);
	if(flag)
	{
		listofRansac.clear();
		delete[]pair1;
		pair1 = NULL;
		delete[]MatchingPoints;
		MatchingPoints = NULL;
		return 1;
	}
	delete[]pair1;
	pair1 = NULL;
	delete[]MatchingPoints;
	MatchingPoints = NULL;	
	_PAIR2* RansacPoints = new _PAIR2[listofRansac.size()];
	memset(RansacPoints,0,listofRansac.size() * sizeof(_PAIR2));//保存细匹配点对
	int m_nNum = listofRansac.size();
	int pair2_N = -1;
	for( list<_PAIR2>::iterator iter = listofRansac.begin();iter != listofRansac.end();iter++)
	{
		pair2_N++;
		RansacPoints[pair2_N].i = iter->i;
		RansacPoints[pair2_N].j = iter->j;
		RansacPoints[pair2_N].m = iter->m;
		RansacPoints[pair2_N].n = iter->n;
	}
	listofRansac.clear();
	//L-M优化变化矩阵H
	L_MOptimize(RansacPoints,m_nNum,0);
	//L-M优化变化矩阵InvH
	L_MOptimize(RansacPoints,m_nNum,1);
	delete[]RansacPoints;
	RansacPoints = NULL;
	return 0;
}
//金字塔匹配时，底层以上图像的匹配
bool CHarrisCorners::cornersMatching(float*& Ii,float*& Ji,bool*& corner1, bool*& corner2,short m_nHeight1,short m_nWidth1,short m_nHeight2, short m_nWidth2, 
									 bool m_bTwice, bool m_bOnce, bool m_bOrigin)
{
	nHeight1 = m_nHeight1;
	nWidth1 = m_nWidth1;
	nHeight2 = m_nHeight2;
	nWidth2 = m_nWidth2;
	if(m_bTwice == 1)
	{
		//粗匹配
		coarseMatching(Ii,Ji,nHeight1,nWidth1,nHeight2,nWidth2,corner1,corner2,0.85f);
		int MatchingPoints_N = -1;
		_TEMP* MatchingPoints = new _TEMP[m_listTribl.size()];
		memset(MatchingPoints,0,m_listTribl.size() * sizeof(_TEMP));
		list<_PAIR>::iterator iter;
		for(iter = m_listTribl.begin();iter != m_listTribl.end();iter++)
		{
			MatchingPoints_N++;
			MatchingPoints[MatchingPoints_N].i = iter->i;
			MatchingPoints[MatchingPoints_N].j = iter->j;
			MatchingPoints[MatchingPoints_N].m = iter->m;
			MatchingPoints[MatchingPoints_N].n = iter->n;
		}
		m_listTribl.clear();
		delete[]Ii;
		delete[]Ji;
		Ii = NULL;
		Ji = NULL;
		//进行细匹配
		int max_N = 0;
		bool flag = ransacMatching(MatchingPoints,MatchingPoints_N,max_N,nHeight1,nWidth1,nHeight2,nWidth2);//求H
		if(flag)
		{
			delete[]MatchingPoints;
			MatchingPoints = NULL;
			return 1;
		}
		int pair1_N = 0;
		for(int a = 0; a <= MatchingPoints_N; a++)
		{
			if(MatchingPoints[a].inliers == 1)
			{
				pair1_N++;
			}
		}
		_PAIR2 *pair1 = new _PAIR2[pair1_N];
		memset(pair1,0,pair1_N * sizeof(_PAIR2));//保存细匹配点对
		pair1_N = -1;
		for(int a = 0; a <= MatchingPoints_N; a++)
		{
			if(MatchingPoints[a].inliers == 1)
			{
				pair1_N++;
				pair1[pair1_N].i = MatchingPoints[a].i;
				pair1[pair1_N].j = MatchingPoints[a].j;
				pair1[pair1_N].m = MatchingPoints[a].m;
				pair1[pair1_N].n = MatchingPoints[a].n;
			}
		}
		list<_PAIR2> listofRansac;
		flag = ransacMatching(MatchingPoints,MatchingPoints_N,max_N,pair1,pair1_N,listofRansac,nHeight1,nWidth1,nHeight2,nWidth2);
		if(flag)
		{
			listofRansac.clear();
			delete[]pair1;
			pair1 = NULL;
			delete[]MatchingPoints;
			MatchingPoints = NULL;
			return 1;
		}
		delete[]pair1;
		pair1 = NULL;
		delete[]MatchingPoints;
		MatchingPoints = NULL;	
		RansacPoints2 = new _PAIR2[listofRansac.size()];
		memset(RansacPoints2,0,listofRansac.size() * sizeof(_PAIR2));//保存细匹配点对
		m_nNum2 = listofRansac.size();
		int pair2_N = -1;
		for( list<_PAIR2>::iterator iter = listofRansac.begin();iter != listofRansac.end();iter++)
		{
			pair2_N++;
			RansacPoints2[pair2_N].i = iter->i;
			RansacPoints2[pair2_N].j = iter->j;
			RansacPoints2[pair2_N].m = iter->m;
			RansacPoints2[pair2_N].n = iter->n;
		}
		listofRansac.clear();
	}
	if(m_bOnce == 1)
	{
		//粗匹配
		coarseMatching(Ii,Ji,nHeight1,nWidth1,nHeight2,nWidth2,corner1,corner2,RansacPoints2,m_nNum2,0.8f,10,1);
		delete[]RansacPoints2;
		RansacPoints2 = NULL;
		int MatchingPoints_N = -1;
		m_nNum1 = m_listOnce.size();
		RansacPoints1 = new _PAIR2[m_listOnce.size()];
		memset(RansacPoints1,0,m_listOnce.size() * sizeof(_PAIR2));//保存细匹配点对
		list<_PAIR2>::iterator iter;
		for(iter = m_listOnce.begin();iter != m_listOnce.end();iter++)
		{
			MatchingPoints_N++;
			RansacPoints1[MatchingPoints_N].i = iter->i;
			RansacPoints1[MatchingPoints_N].j = iter->j;
			RansacPoints1[MatchingPoints_N].m = iter->m;
			RansacPoints1[MatchingPoints_N].n = iter->n;
		}
		m_listOnce.clear();
		delete[]Ii;
		delete[]Ji;
		Ii = NULL;
		Ji = NULL;
	}
	if(m_bOrigin == 1)
	{
		//粗匹配
		coarseMatching(Ii,Ji,nHeight1,nWidth1,nHeight2,nWidth2,corner1,corner2,RansacPoints1,m_nNum1,0.8f,7,0);
		delete[]RansacPoints1;
		RansacPoints1 = NULL;
		int MatchingPoints_N = -1;
		_TEMP* MatchingPoints = new _TEMP[m_listOrigin.size()];
		memset(MatchingPoints,0,m_listOrigin.size() * sizeof(_TEMP));
		list<_PAIR2>::iterator iter;
		for(iter = m_listOrigin.begin();iter != m_listOrigin.end();)
		{
			int i = iter->i;
			int j = iter->j;
			int m = iter->m;
			int n = iter->n;
			float R1 = Corres(Ii,Ji,nWidth1,nWidth2,i,j,m,n);
			float max = R1;
			bool IsMax = 0;
			for(int b = j-2;b <= j+2;b++)
			{
				for(int a  = i-2;a <= i+2;a++)
				{
					for(int d = n-2;d <= n+2;d++)
					{
						for(int c  = m-2;c <= m+2;c++)
						{
							float R = Corres(Ii,Ji,nWidth1,nWidth2,a,b,c,d);
							if(R > max)
							{
								IsMax = 1;
								max = R;
							}
						}
					}
				}
			}
			if(!IsMax)
			{
				iter++;
			}
			else
			{
				for(int b = j-2;b <= j+2;b++)
				{
					for(int a  = i-2;a <= i+2;a++)
					{
						for(int d = n-2;d <= n+2;d++)
						{
							for(int c  = m-2;c <= m+2;c++)
							{
								float R = Corres(Ii,Ji,nWidth1,nWidth2,a,b,c,d);
								if(R == max)
								{
									_PAIR2 pair;
									pair.flag = 0;
									pair.i = a;
									pair.j = b;
									pair.m = c;
									pair.n = d;
									m_listOrigin.push_back(pair);
								}
							}
						}
					}
				}
				iter = m_listOrigin.erase(iter);
			}

		}
		for(iter = m_listOrigin.begin();iter != m_listOrigin.end();iter++)
		{
			MatchingPoints_N++;
			MatchingPoints[MatchingPoints_N].i = iter->i;
			MatchingPoints[MatchingPoints_N].j = iter->j;
			MatchingPoints[MatchingPoints_N].m = iter->m;
			MatchingPoints[MatchingPoints_N].n = iter->n;
		}
		m_listOrigin.clear();
		delete[]Ii;
		delete[]Ji;
		Ii = NULL;
		Ji = NULL;
		//进行细匹配
		int max_N = 0;
		bool flag = ransacMatching(MatchingPoints,MatchingPoints_N,max_N,nHeight1,nWidth1,nHeight2,nWidth2);//求H
		if(flag)
		{
			delete[]MatchingPoints;
			MatchingPoints = NULL;
			return 1;
		}
		int pair1_N = 0;
		for(int a = 0; a <= MatchingPoints_N; a++)
		{
			if(MatchingPoints[a].inliers == 1)
			{
				pair1_N++;
			}
		}
		_PAIR2 *pair1 = new _PAIR2[pair1_N];
		memset(pair1,0,pair1_N * sizeof(_PAIR2));//保存细匹配点对
		pair1_N = -1;
		for(int a = 0; a <= MatchingPoints_N; a++)
		{
			if(MatchingPoints[a].inliers == 1)
			{
				pair1_N++;
				pair1[pair1_N].i = MatchingPoints[a].i;
				pair1[pair1_N].j = MatchingPoints[a].j;
				pair1[pair1_N].m = MatchingPoints[a].m;
				pair1[pair1_N].n = MatchingPoints[a].n;
			}
		}
		list<_PAIR2> listofRansac;
		flag = ransacMatching(MatchingPoints,MatchingPoints_N,max_N,pair1,pair1_N,listofRansac,nHeight1,nWidth1,nHeight2,nWidth2);
		if(flag)
		{
			listofRansac.clear();
			delete[]pair1;
			pair1 = NULL;
			delete[]MatchingPoints;
			MatchingPoints = NULL;
			return 1;
		}
		delete[]pair1;
		pair1 = NULL;
		delete[]MatchingPoints;
		MatchingPoints = NULL;	
		RansacPoints0 = new _PAIR2[listofRansac.size()];
		memset(RansacPoints0,0,listofRansac.size() * sizeof(_PAIR2));//保存细匹配点对
		m_nNum0 = listofRansac.size();
		int pair2_N = -1;
		for( list<_PAIR2>::iterator iter = listofRansac.begin();iter != listofRansac.end();iter++)
		{
			pair2_N++;
			RansacPoints0[pair2_N].i = iter->i;
			RansacPoints0[pair2_N].j = iter->j;
			RansacPoints0[pair2_N].m = iter->m;
			RansacPoints0[pair2_N].n = iter->n;
		}
		listofRansac.clear();
		//L-M优化变化矩阵H
		L_MOptimize(RansacPoints0,m_nNum0,0);
		//L-M优化变化矩阵InvH
		L_MOptimize(RansacPoints0,m_nNum0,1);
		delete[]RansacPoints0;
		RansacPoints0 = NULL;
	}
	return 0;
}
//不需要进行金字塔匹配时的原始图像角点粗匹配，或者是金字塔匹配时底层图像的角点粗匹配
void CHarrisCorners::coarseMatching(float*& Ig, float*& Jg, short height1, short width1, short height2, short width2, bool*& corner1, bool*& corner2, float dbMaxR)
{
	#define Jg(ROW,COL)    Jg[width2 * (COL) + (ROW)]
	#define Ig(ROW,COL)    Ig[width1 * (COL) + (ROW)]
	#define corner1(ROW,COL)    corner1[width1 * (COL) + (ROW)]
	#define corner2(ROW,COL)    corner2[width2 * (COL) + (ROW)]
	list<_CORNER> listofcorner1;
	list<_CORNER>::iterator iterc1;
	list<_CORNER> listofcorner2;
	list<_CORNER>::iterator iterc2;
	for(short j = 7;j < height1 - 7;j++)
	{												
		for(short i = 7;i < width1 - 7;i++)
		{
			if(corner1(i,j) == 1)
			{
				_CORNER corner;
				corner.i = i;
				corner.j = j;
				GetM(Ig,width1,i,j,corner.M);
				listofcorner1.push_back(corner);
			}
		}
	}
	for(short j = 7;j < height2 - 7;j++)
	{												
		for(short i = 7;i < width2 - 7;i++)
		{
			if(corner2(i,j) == 1)
			{
				_CORNER corner;
				corner.i = i;
				corner.j = j;
				GetM(Jg,width2,i,j,corner.M);
				listofcorner2.push_back(corner);
			}
		}
	}
	delete[]corner2;
	corner2 = NULL;
	delete[]corner1;
	corner1 = NULL;
	list<_CORNER> list1;
	list<_CORNER>::iterator iter1;
	for(iterc1 = listofcorner1.begin(); iterc1 != listofcorner1.end();iterc1++)
	{	
		list<_DOT> dot;
		int i = iterc1->i;
		int j = iterc1->j;
		float max = 0.8f;

		for(iterc2 = listofcorner2.begin(); iterc2 != listofcorner2.end();iterc2++)
		{										
			int m = iterc2->i;
			int n = iterc2->j;
			float R = corresX(iterc1->M,iterc2->M,37);
			if(R > dbMaxR)
			{
				_DOT temp;
				temp.data = R;
				temp.i = m;
				temp.j = n;
				memcpy(temp.M,iterc2->M,37 * sizeof(unsigned char));
				if(R >= max)
				{
					max = R;
					dot.push_back(temp);
				}	
			}							
		}
		list<_DOT>::iterator iter;
		for(iter = dot.begin();iter != dot.end();)
		{
			if(iter->data != max)
			{
				iter = dot.erase(iter);
			}
			else
			{
				iter++;
			}
		}
		if(list1.size() == 0)
		{
			for(iter = dot.begin();iter != dot.end();iter++)
			{
				_CORNER corner;
				corner.i = iter->i;
				corner.j = iter->j;
				memcpy(corner.M,iter->M,37 * sizeof(unsigned char));
				list1.push_back(corner);
			}
		}
		else
		{
			for(iter = dot.begin();iter != dot.end();iter++)
			{
				_CORNER corner;
				corner.i = iter->i;
				corner.j = iter->j;
				memcpy(corner.M,iter->M,37 * sizeof(unsigned char));
				bool flag = 0;
				for(iter1 = list1.begin();iter1 != list1.end();iter1++)
				{
					if(iter->i == iter1->i && iter->j == iter1->j)
					{
						flag = 1;
						break;
					}
				}
				if(flag == 0)
				{
					list1.push_back(corner);
				}
			}
		}
		dot.clear();

	}

	list<_PAIR>::iterator iter2;
	for(iter1 = list1.begin();iter1 != list1.end();iter1++)
	{
		short m = iter1->i;
		short n = iter1->j;
		float max = 0.8f;
		list<_DOT> dot;
		for(iterc1 = listofcorner1.begin(); iterc1 != listofcorner1.end();iterc1++)
		{	
			int i = iterc1->i;
			int j = iterc1->j;
			float R = corresX(iterc1->M,iter1->M,37);
			if(R > dbMaxR)
			{
				_DOT temp;
				temp.data = R;
				temp.i = i;
				temp.j = j;
				memcpy(temp.M,iterc1->M,37 * sizeof(unsigned char));
				if(R >= max)
				{
					max = R;
					dot.push_back(temp);
				}	
			}
		}
		list<_DOT>::iterator iter;
		for(iter = dot.begin();iter != dot.end();)
		{
			if(iter->data != max)
			{
				iter = dot.erase(iter);
			}
			else
			{
				iter++;
			}
		}

		if(m_listTribl.size() == 0)
		{
			for(iter = dot.begin();iter != dot.end();iter++)
			{
				_PAIR pair;
				pair.data = iter->data;
				pair.flag = 1;
				pair.i = iter->i;
				pair.j = iter->j;
				pair.m = m;
				pair.n = n;
				m_listTribl.push_back(pair);
			}
		}
		else
		{
			for(iter = dot.begin();iter != dot.end();iter++)
			{
				_PAIR pair;
				pair.data = iter->data;
				pair.flag = 1;
				pair.i = iter->i;
				pair.j = iter->j;
				pair.m = m;
				pair.n = n;
				bool flag = 0;
				for(iter2 = m_listTribl.begin();iter2 != m_listTribl.end();iter2++)
				{
					if(iter->i == iter2->i && iter->j == iter2->j && iter->data < iter2->data)
					{
						flag = 1;
						break;
					}
					if(iter->i == iter2->i && iter->j == iter2->j && iter->data > iter2->data)
					{
						iter2->flag = 0;
					}
				}
				if(flag == 0)
				{
					m_listTribl.push_back(pair);
				}
			}
		}
		dot.clear();

	}
	list1.clear();
	listofcorner1.clear();
	listofcorner2.clear();
	for(iter2 = m_listTribl.begin();iter2 != m_listTribl.end();)
	{
		if(iter2->flag == 0)
		{
			iter2 = m_listTribl.erase(iter2);
		}
		else
		{
			iter2++;
		}
	}
	int a = 0;
}
//金字塔匹配时，底层以上图像的角点粗匹配
void CHarrisCorners::coarseMatching(float*& Ig, float*& Jg, short height1, short width1, short height2, short width2, 
										bool*& _corner1, bool*& _corner2, _PAIR2*& RansacPoints,int num, float dbMax, int nRadius,int nWvltNum)
{
	#define Jg(ROW,COL)    Jg[width2 * (COL) + (ROW)]
	#define Ig(ROW,COL)    Ig[width1 * (COL) + (ROW)]
	#define _corner1(ROW,COL)    _corner1[width1 * (COL) + (ROW)]
	#define _corner2(ROW,COL)    _corner2[width2 * (COL) + (ROW)]
	short minx1 = 2 * RansacPoints[0].i;
	short miny1 = 2 * RansacPoints[0].j;
	short minx2 = 2 * RansacPoints[0].m;
	short miny2 = 2 * RansacPoints[0].n;

	short maxx1 = minx1;
	short maxy1 = miny1;
	short maxx2 = minx2;
	short maxy2 = miny2;
	for(int x = 0;x < num;x++)
	{
		short i = 2 * RansacPoints[x].i;
		short j = 2 * RansacPoints[x].j;
		short m = 2 * RansacPoints[x].m;
		short n = 2 * RansacPoints[x].n;
		if(i < minx1)
		{
			minx1 = i;
		}
		if(i > maxx1)
		{
			maxx1 = i;
		}
		if(j < miny1)
		{
			miny1 = j;
		}
		if(j > maxy1)
		{
			maxy1 = j;
		}
		if(m < minx2)
		{
			minx2 = m;
		}
		if(m > maxx2)
		{
			maxx2 = m;
		}
		if(n < miny2)
		{
			miny2 = n;
		}
		if(n > maxy2)
		{
			maxy2 = n;
		}
	}
	bool *IsEdge1 = new bool[height1*width1];
	memset(IsEdge1,0,height1*width1*sizeof(bool));
	prescreening(height1 ,width1 ,Ig ,IsEdge1);
	harrisCornerDetecting(height1,width1,Ig ,_corner1,IsEdge1,minx1,maxx1,miny1,maxy1);
	delete[]IsEdge1;
	IsEdge1 = NULL;

	bool *IsEdge2 = new bool[height2*width2];
	memset(IsEdge2,0,height2*width2*sizeof(bool));
	prescreening(height2 ,width2 ,Jg ,IsEdge2);
	harrisCornerDetecting(height2,width2,Jg ,_corner2,IsEdge2,minx2,maxx2,miny2,maxy2);
	delete[]IsEdge2;
	IsEdge2 = NULL;

	for(int x = 0;x < num;x++)
	{
		short i = 2 * RansacPoints[x].i;
		short j = 2 * RansacPoints[x].j;
		short m = 2 * RansacPoints[x].m;
		short n = 2 * RansacPoints[x].n;
		if(i >= nRadius && i < width1 - nRadius && j >= nRadius && j < height1 - nRadius && m >= nRadius && m < width2 - nRadius && n >= nRadius && n < height2 - nRadius)
		{

			list<_PAIR2> list0;
			for(short b = j - nRadius;b <= j + nRadius;b++)
			{
				for(short a = i - nRadius;a <= i + nRadius;a++)
				{
					if(_corner1(a,b) == 1)
					{
						list<_DOT> dot;
						float max = 0.7f;
						for(short d = n - nRadius;d <= n + nRadius;d++)
						{										//遍历待配准的图像，从左上角开始
							for(short c = m - nRadius;c <= m + nRadius;c++)
							{
								if(_corner2(c,d) == 1)
								{
									float R = Corres(Ig,Jg,width1,width2,a,b,c,d);
									if(R > dbMax)
									{
										_DOT temp;
										temp.data = R;
										temp.i = c;
										temp.j = d;
										if(R >= max)
										{
											max = R;
											dot.push_back(temp);
										}	
									}
								}
							}							
						}
						list<_DOT>::iterator iter;
						for(iter = dot.begin();iter != dot.end();)
						{
							if(iter->data != max)
							{
								iter = dot.erase(iter);
							}
							else
							{
								iter++;
							}
						}
						for(iter = dot.begin();iter != dot.end();iter++)
						{
								_PAIR2 pair;
								pair.flag = 0;
								pair.i = a;
								pair.j = b;
								pair.m = iter->i;
								pair.n = iter->j;
								list0.push_back(pair);
						}
						dot.clear();
					}
				}
			}

			list<_PAIR2> list1;
			for(short b = n - nRadius;b <= n + nRadius;b++)
			{
				for(short a = m - nRadius;a <= m + nRadius;a++)
				{
					if(_corner2(a,b) == 1)
					{
						list<_DOT> dot;
						float max = 0.7f;
						for(short d = j - nRadius;d <= j + nRadius;d++)
						{										//遍历待配准的图像，从左上角开始
							for(short c = i - nRadius;c <= i + nRadius;c++)
							{
								if(_corner1(c,d) == 1)
								{
									float R = Corres(Jg,Ig,width2,width1,a,b,c,d);
									if(R > dbMax)
									{
										_DOT temp;
										temp.data = R;
										temp.i = c;
										temp.j = d;
										if(R >= max)
										{
											max = R;
											dot.push_back(temp);
										}	
									}
								}
							}							
						}
						list<_DOT>::iterator iter;
						for(iter = dot.begin();iter != dot.end();)
						{
							if(iter->data != max)
							{
								iter = dot.erase(iter);
							}
							else
							{
								iter++;
							}
						}
						for(iter = dot.begin();iter != dot.end();iter++)
						{
								_PAIR2 pair;
								pair.flag = 0;
								pair.m = a;
								pair.n = b;
								pair.i = iter->i;
								pair.j = iter->j;
								list1.push_back(pair);
						}
						dot.clear();
					}
				}
			}
			for(list<_PAIR2>::iterator iter0 = list0.begin();iter0 != list0.end();)
			{
				bool f = 0;
				for(list<_PAIR2>::iterator iter1 = list1.begin();iter1 != list1.end();)
				{
					if(iter0->i==iter1->i && iter0->j==iter1->j && iter0->m==iter1->m && iter0->n==iter1->n)
					{
						iter1 = list1.erase(iter1);
						f = 1;
						break;
					}
					else
					{
						iter1++;
					}
				}
				if(f == 0)
				{
					iter0 = list0.erase(iter0);
				}
				else
				{
					_PAIR2 pair;
					pair.flag = 0;
					pair.m = iter0->m;
					pair.n = iter0->n;
					pair.i = iter0->i;
					pair.j = iter0->j;
				    if(nWvltNum == 2)
					{
						m_listTwice.push_back(pair);
					}
				    if(nWvltNum == 1)
					{
						m_listOnce.push_back(pair);
					}
				    if(nWvltNum == 0)
					{
						m_listOrigin.push_back(pair);
					}
					iter0++;
				}
			}
			list0.clear();
			list1.clear();

		}
	}
	delete[]_corner1;
	_corner1 = NULL;
	delete[]_corner2;
	_corner2 = NULL;

}
//第一次Ransac匹配
bool CHarrisCorners::ransacMatching(_TEMP*& _temp, int n,int &maxN, short H1, short W1, short H2, short W2)//_temp表示粗匹配点对，n是粗匹配点对数，
																											//maxN表示第一次Ransac匹配获得的内点数，H、W是匹配图像的高和宽
{
	if(n <= 4)
	{
		AfxMessageBox(_T("匹配错误！"));
		return 1;
	}

	//获得采样次数N
	short N = 0;
	float p = 0.95f;//粗匹配点为内点的概率
	float epsilon = 0.8f;//任何一对匹配点不是内点的概率
	N = short(log(1 - p) / log(1 - (1 - epsilon) * (1 - epsilon) * (1 - epsilon) * (1 - epsilon)) + 0.5f);
	//N = 2500;
	//short nume = int((float(n)) * 0.2f);
	short nume = 4;
	short* S = new short[N];        //每次采样处理后获得的内点个数
	memset(S,0,(N) * sizeof(short));
	float* sigma2 = new float[N];     //每次采样获得的个对应点间的欧式距离的方差
	memset(sigma2,0,N * sizeof(float));
	short* NUM = new short[N * 4];        //每次采样处理后获得的内点个数
	memset(NUM,0,(N * 4) * sizeof(short));

	srand((unsigned)time(NULL));
	short* max = new short[N];
	memset(max,0,N * sizeof(short));
	short k = 0;

	for(short i = 0;i < N;i++)
	{
		//随机选取四对匹配点对
flag0:	short num = -1;
		short num0 = 0;
		short num1 = 0;
		short num2 = 0;
		short num3 = 0;
		short g = 0; 
		float pDbProjectionPara0[8]; //保存变换矩阵系数
		float H0[9]; //保存变换矩阵（从待配准到基准）
		float pDbProjectionPara1[8]; //保存变换矩阵系数
		float InvH1[9]; //保存变换矩阵（从ji准到daipei准）
		_PAIR2 temp[4];
		float sum1 = 0.0f;
		float sum2 = 0.0f;
		num0 = rand() % (n + 1);
		temp[0].i = _temp[num0].i;
		temp[0].j = _temp[num0].j;
		temp[0].m = _temp[num0].m;
		temp[0].n = _temp[num0].n;
		g = rand() % (n + 1);
		while(g == num0)
		{
			g = rand() % (n + 1);
		}
		num1 = g;
		temp[1].i = _temp[num1].i;
		temp[1].j = _temp[num1].j;
		temp[1].m = _temp[num1].m;
		temp[1].n = _temp[num1].n;
		g = rand() % (n + 1);
		while(g == num0 || g == num1)
		{
			g = rand() % (n + 1);
		}
		num2 = g;
		temp[2].i = _temp[num2].i;
		temp[2].j = _temp[num2].j;
		temp[2].m = _temp[num2].m;
		temp[2].n = _temp[num2].n;
		g = rand() % (n + 1);
		while(g == num0 || g == num1 || g == num2)
		{
			g = rand() % (n + 1);
		}
		num3 = g;
		temp[3].i = _temp[num3].i;
		temp[3].j = _temp[num3].j;
		temp[3].m = _temp[num3].m;
		temp[3].n = _temp[num3].n;

		NUM[3 * i + 0] = num0;
		NUM[3 * i + 1] = num1;
		NUM[3 * i + 2] = num2;
		NUM[3 * i + 3] = num3;
		bool IsContinue = 0;
		if(i > 0)
		{
			for(int a = 0;a < i;a++)
			{
				if(  (NUM[3 * i + 0] != NUM[3 * a + 0] && NUM[3 * i + 0] != NUM[3 * a + 1] && NUM[3 * i + 0] != NUM[3 * a + 2] && NUM[3 * i + 0] != NUM[3 * a + 3]) || 
				     (NUM[3 * i + 1] != NUM[3 * a + 0] && NUM[3 * i + 1] != NUM[3 * a + 1] && NUM[3 * i + 1] != NUM[3 * a + 2] && NUM[3 * i + 1] != NUM[3 * a + 3]) || 
				     (NUM[3 * i + 2] != NUM[3 * a + 0] && NUM[3 * i + 2] != NUM[3 * a + 1] && NUM[3 * i + 2] != NUM[3 * a + 2] && NUM[3 * i + 2] != NUM[3 * a + 3]) || 
				     (NUM[3 * i + 3] != NUM[3 * a + 0] && NUM[3 * i + 3] != NUM[3 * a + 1] && NUM[3 * i + 3] != NUM[3 * a + 2] && NUM[3 * i + 3] != NUM[3 * a + 3])  
				   )
				{
				}
				else
				{
					IsContinue = 1;
					break;
				}
			}
		}
		if(IsContinue)
		{
			goto flag0;
		}
		if(judgeFourPoints(temp,H1,W1,H2,W2,4) == 0)
		{
			goto flag0;
		}
		IsContinue = 0;
		for(short a = 0;a < 3;a++)
		{
			for(short b = a + 1;b < 4;b++)
			{
				short x1 = temp[a].i;
				short y1 = temp[a].j;
				short w1 = temp[a].m;
				short z1 = temp[a].n;
				short x2 = temp[b].i;
				short y2 = temp[b].j;
				short w2 = temp[b].m;
				short z2 = temp[b].n;
				float t1 = sqrt(float(((x1 - w1) - (x2 - w2)) * ((x1 - w1) - (x2 - w2)) + ((y1 - z1) - (y2 - z2)) * ((y1 - z1) - (y2 - z2))));
				float t2 = sqrt(float(((x1 - w1) + (x2 - w2)) * ((x1 - w1) + (x2 - w2)) + ((y1 - z1) + (y2 - z2)) * ((y1 - z1) + (y2 - z2))));
				float g = 2 * (t1 / t2);
				if(g > 2)
				{
					IsContinue = 1;
					break;
				}
			}
		}
		if(IsContinue)
		{
			goto flag0;
		}
		// 根据前面所得四个粗匹配点计算从待配准图象到基准图象的仿射变换系数
		getProjectionPara(temp,pDbProjectionPara0,0);
		// 根据前面所得四个粗匹配点计算从ji准图象到daipei准图象的仿射变换系数
		getProjectionPara(temp,pDbProjectionPara1,1);
		for(short a = 0;a < 9;a++)
		{
			if(a == 8)
			{
				 H0[8] = 1;
				 InvH1[8] = 1;
				 break;
			}
			H0[a] = pDbProjectionPara0[a];
			InvH1[a] = pDbProjectionPara1[a];
		}
		//计算相应的欧式距离平方
		for(int a = 0;a <= n;a++)
		{
			float x = H0[0] * (_temp[a].m) + H0[1] * _temp[a].n + H0[2];
			float y = H0[3] * (_temp[a].m) + H0[4] * _temp[a].n + H0[5];
			float w = H0[6] * (_temp[a].m) + H0[7] * _temp[a].n + 1;
			x = x / w;
			y = y / w;
			_temp[a].D0 = sqrt((_temp[a].i - x) * (_temp[a].i - x) + (_temp[a].j - y) * (_temp[a].j - y));
			float x1 = 0;
			float y1 = 0;
			float w1 = 0;
			float u1 = 0;
			float v1 = 0;
			x1 = InvH1[0] * (_temp[a].i) + InvH1[1] * _temp[a].j + InvH1[2];
			y1 = InvH1[3] * (_temp[a].i) + InvH1[4] * _temp[a].j + InvH1[5];
			w1 = InvH1[6] * (_temp[a].i) + InvH1[7] * _temp[a].j + 1;
			x1 = x1 / w1;
			y1 = y1 / w1;
			_temp[a].D1 = sqrt((_temp[a].m - x1) * (_temp[a].m - x1) + (_temp[a].n - y1) * (_temp[a].n - y1));
		}

		float t = 1.5f;
		if(i == 0)
		{
			for(int a = 0;a <= n;a++)	//统计每次采样获得的内点个数
			{
				if(_temp[a].D0 < t && _temp[a].D1 < t)
				{
					sum1 += (_temp[a].D0 + _temp[a].D1)/2; 
					S[i]++;
				}
			}
			if(S[i] < nume)
			{
				S[i] = 0;
				goto flag0;
			}
			float u = sum1 / S[i];
			for(int a = 0;a <= n;a++)	
			{
				if(_temp[a].D0 < t && _temp[a].D1 < t)
				{
					sum2 += ((_temp[a].D0 + _temp[a].D1)/2 - u) * ((_temp[a].D0 + _temp[a].D1)/2 - u);
					_temp[a].inliers = 1;
				}
			}
			sigma2[i] = sum2 / S[i];
			k = i;
			max[k] = S[i];
			maxN = S[i];
			memcpy(H,H0,9 * sizeof(float));
			memcpy(InvH,InvH1,9 * sizeof(float));
			continue;
					
		}
		if(i >= 1)
		{
			for(int a = 0;a <= n;a++)	//统计每次采样获得的内点个数
			{
				if(_temp[a].D0 < t && _temp[a].D1 < t)
				{
					S[i]++;
					sum1 += (_temp[a].D0 + _temp[a].D1)/2;
				}
			}
			if(S[i] < nume)
			{
				S[i] = 0;
				goto flag0;
			}
			float u = sum1 / S[i];
			for(int a = 0;a <= n;a++)	
			{
				if(_temp[a].D0 < t && _temp[a].D1 < t)
				{
					sum2 += (_temp[a].D0 - u) * (_temp[a].D0 - u); 
					sum2 += ((_temp[a].D0 + _temp[a].D1)/2 - u) * ((_temp[a].D0 + _temp[a].D1)/2 - u);
				}
			}
			sigma2[i] = sum2 / S[i];
			if(S[i] > max[k])
			{
				memcpy(H,H0,9 * sizeof(float));
				memcpy(InvH,InvH1,9 * sizeof(float));
				k = i;
				max[k] = S[i];
				maxN = S[i];
				for(int a = 0;a <= n;a++)
				{
					if(_temp[a].D0 < t && _temp[a].D1 < t && u < 1.0f)
					{
						_temp[a].inliers = 1;
					}
				}
				continue;
			}
			if(S[i] == max[k])
			{
				if(sigma2[i] < sigma2[k])
				{
					k = i;
					max[k] = S[i];
					maxN = S[i];
					memcpy(H,H0,9 * sizeof(float));
					memcpy(InvH,InvH1,9 * sizeof(float));
					for(int a = 0;a <= n;a++)
					{
						if(_temp[a].D0 < t && _temp[a].D1 < t && u < 1.0f)
						{
							_temp[a].inliers = 1;
						}
					}
				}
			}
		}
	}

	if(max[k] <= 4)
	{
		AfxMessageBox(_T("匹配错误！"));
		delete[]S;
		delete[]sigma2;
		delete[]NUM;
		delete[]max;
		S = NULL;
		sigma2 = NULL;
		NUM = NULL;
		max = NULL;
		return 1;
	}
	delete[]S;
	delete[]sigma2;
	delete[]NUM;
	delete[]max;
	S = NULL;
	sigma2 = NULL;
	NUM = NULL;
	max = NULL;
	return 0;
}

bool CHarrisCorners::ransacMatching(_TEMP*& _temp, int n1,int max_N, _PAIR2* pair, int n2,  list<_PAIR2> &listofRansac, short H1, short W1, short H2, short W2)
{
	if(n2 <= 4)
	{
		AfxMessageBox(_T("匹配错误！"));
		return 1;
	}
	//获得采样次数N
	short N = 0;
	float p = 0.99f;//粗匹配点为内点的概率
	float epsilon = 0.7f;//任何一对匹配点不是内点的概率
	N = 2*short(log(1 - p) / log(1 - (1 - epsilon) * (1 - epsilon) * (1 - epsilon) * (1 - epsilon)) + 0.5f);
	int N0 = (n2+1)*(n2)*(n2-1)*(n2-2)/24;
	short nume = max_N;
	short* S = new short[N];        //每次采样处理后获得的内点个数
	memset(S,0,(N) * sizeof(short));
	float* sigma2 = new float[N];     //每次采样获得的个对应点间的欧式距离的方差
	memset(sigma2,0,N * sizeof(float));
	short* NUM = new short[N * 4];        //每次采样处理后获得的内点个数
	memset(NUM,0,(N * 4) * sizeof(short));

	srand((unsigned)time(NULL));
	short* max = new short[N];
	memset(max,0,N * sizeof(short));
	short k = 0;
	float fAngle = 2;
	for(short i = 0;i < N;i++)
	{
		//随机选取四对匹配点对
flag0:	short num = -1;
		short num0 = 0;
		short num1 = 0;
		short num2 = 0;
		short num3 = 0;

		short g = 0; 
		float pDbProjectionPara0[8]; //保存变换矩阵系数
		float H0[9]; //保存变换矩阵（从待配准到基准）
		float pDbProjectionPara1[8]; //保存变换矩阵系数
		float InvH1[9]; //保存变换矩阵（从ji准到daipei准）
		_PAIR2 temp[4];
		float sum1 = 0.0f;
		float sum2 = 0.0f;
		num0 = rand() % (n2 + 1);
		temp[0].flag = 0;
		temp[0].i = pair[num0].i;
		temp[0].j = pair[num0].j;
		temp[0].m = pair[num0].m;
		temp[0].n = pair[num0].n;
		g = rand() % (n2 + 1);
		while(g == num0)
		{
			g = rand() % (n2 + 1);
		}
		num1 = g;
		temp[1].flag = 0;
		temp[1].i = pair[num1].i;
		temp[1].j = pair[num1].j;
		temp[1].m = pair[num1].m;
		temp[1].n = pair[num1].n;
		g = rand() % (n2 + 1);
		while(g == num0 || g == num1)
		{
			g = rand() % (n2 + 1);
		}
		num2 = g;
		temp[2].flag = 0;
		temp[2].i = pair[num2].i;
		temp[2].j = pair[num2].j;
		temp[2].m = pair[num2].m;
		temp[2].n = pair[num2].n;
		g = rand() % (n2 + 1);
		while(g == num0 || g == num1 || g == num2)
		{
			g = rand() % (n2 + 1);
		}
		num3 = g;
		temp[3].flag = 0;
		temp[3].i = pair[num3].i;
		temp[3].j = pair[num3].j;
		temp[3].m = pair[num3].m;
		temp[3].n = pair[num3].n;


		NUM[3 * i + 0] = num0;
		NUM[3 * i + 1] = num1;
		NUM[3 * i + 2] = num2;
		NUM[3 * i + 3] = num3;
		bool IsContinue = 0;
		if(i > 0)
		{
			for(int a = 0;a < i;a++)
			{
				if(  (NUM[3 * i + 0] != NUM[3 * a + 0] && NUM[3 * i + 0] != NUM[3 * a + 1] && NUM[3 * i + 0] != NUM[3 * a + 2] && NUM[3 * i + 0] != NUM[3 * a + 3]) || 
				     (NUM[3 * i + 1] != NUM[3 * a + 0] && NUM[3 * i + 1] != NUM[3 * a + 1] && NUM[3 * i + 1] != NUM[3 * a + 2] && NUM[3 * i + 1] != NUM[3 * a + 3]) || 
				     (NUM[3 * i + 2] != NUM[3 * a + 0] && NUM[3 * i + 2] != NUM[3 * a + 1] && NUM[3 * i + 2] != NUM[3 * a + 2] && NUM[3 * i + 2] != NUM[3 * a + 3]) || 
				     (NUM[3 * i + 3] != NUM[3 * a + 0] && NUM[3 * i + 3] != NUM[3 * a + 1] && NUM[3 * i + 3] != NUM[3 * a + 2] && NUM[3 * i + 3] != NUM[3 * a + 3])  
				   )
				{
				}
				else
				{
					IsContinue = 1;
					break;
				}
			}
		}
		if(IsContinue)
		{
			goto flag0;
		}
		if(judgeFourPoints(temp,H1,W1,H2,W2,4) == 0)
		{
			goto flag0;
		}
		IsContinue = 0;
		for(short a = 0;a < 3;a++)
		{
			for(short b = a + 1;b < 4;b++)
			{
				short x1 = temp[a].i;
				short y1 = temp[a].j;
				short w1 = temp[a].m;
				short z1 = temp[a].n;
				short x2 = temp[b].i;
				short y2 = temp[b].j;
				short w2 = temp[b].m;
				short z2 = temp[b].n;
				float t1 = sqrt(float(((x1 - w1) - (x2 - w2)) * ((x1 - w1) - (x2 - w2)) + ((y1 - z1) - (y2 - z2)) * ((y1 - z1) - (y2 - z2))));
				float t2 = sqrt(float(((x1 - w1) + (x2 - w2)) * ((x1 - w1) + (x2 - w2)) + ((y1 - z1) + (y2 - z2)) * ((y1 - z1) + (y2 - z2))));
				float g = 2 * (t1 / t2);
				if(g > 2)
				{
					IsContinue = 1;
					break;
				}
			}
		}
		if(IsContinue)
		{
			goto flag0;
		}
		// 根据前面所得四个粗匹配点计算从待配准图象到基准图象的仿射变换系数
		getProjectionPara(temp,pDbProjectionPara0,0);
		// 根据前面所得四个粗匹配点计算从ji准图象到daipei准图象的仿射变换系数
		getProjectionPara(temp,pDbProjectionPara1,1);
		for(short a = 0;a < 9;a++)
		{
			if(a == 8)
			{
				 H0[8] = 1;
				 InvH1[8] = 1;
				 break;
			}
			H0[a] = pDbProjectionPara0[a];
			InvH1[a] = pDbProjectionPara1[a];
		}
		//计算相应的欧式距离平方
		for(int a = 0;a <= n1;a++)
		{
			float x = H0[0] * (_temp[a].m) + H0[1] * _temp[a].n + H0[2];
			float y = H0[3] * (_temp[a].m) + H0[4] * _temp[a].n + H0[5];
			float w = H0[6] * (_temp[a].m) + H0[7] * _temp[a].n + 1;
			x = x / w;
			y = y / w;
			_temp[a].D0 = sqrt((_temp[a].i - x) * (_temp[a].i - x) + (_temp[a].j - y) * (_temp[a].j - y));
			float x1 = 0;
			float y1 = 0;
			float w1 = 0;
			float u1 = 0;
			float v1 = 0;
			x1 = InvH1[0] * (_temp[a].i) + InvH1[1] * _temp[a].j + InvH1[2];
			y1 = InvH1[3] * (_temp[a].i) + InvH1[4] * _temp[a].j + InvH1[5];
			w1 = InvH1[6] * (_temp[a].i) + InvH1[7] * _temp[a].j + 1;
			x1 = x1 / w1;
			y1 = y1 / w1;
			_temp[a].D1 = sqrt((_temp[a].m - x1) * (_temp[a].m - x1) + (_temp[a].n - y1) * (_temp[a].n - y1));
		}

		float t = 1.5f;
		if(i == 0)
		{
			for(int a = 0;a <= n1;a++)	//统计每次采样获得的内点个数
			{
				if(_temp[a].D0 < t && _temp[a].D1 < t)
				{
					sum1 += (_temp[a].D0 + _temp[a].D1)/2; 
					S[i]++;
					_temp[a].inliers = 1;
				}
				else
				{
					_temp[a].inliers = 0;
				}
			}
			if(S[i] < nume)
			{
				S[i] = 0;
				goto flag0;
			}
			float u = sum1 / S[i];
			for(int a = 0;a <= n1;a++)	
			{
				if(_temp[a].D0 < t && _temp[a].D1 < t)
				{ 
					sum2 += ((_temp[a].D0 + _temp[a].D1)/2 - u) * ((_temp[a].D0 + _temp[a].D1)/2 - u);
					_PAIR2 pair;
					pair.i = _temp[a].i;
					pair.j = _temp[a].j;
					pair.m = _temp[a].m;
					pair.n = _temp[a].n;
					listofRansac.push_back(pair);
				}
			}
			sigma2[i] = sum2 / S[i];
			k = i;
			max[k] = S[i];
			memcpy(H,H0,9 * sizeof(float));
			memcpy(InvH,InvH1,9 * sizeof(float));
			continue;
					
		}
		if(i >= 1)
		{
			for(int a = 0;a <= n1;a++)	//统计每次采样获得的内点个数
			{
				if(_temp[a].D0 < t && _temp[a].D1 < t)
				{
					S[i]++;
					sum1 += (_temp[a].D0 + _temp[a].D1)/2;
				}
			}
			float u = sum1 / S[i];
			for(int a = 0;a <= n1;a++)	
			{
				if(_temp[a].D0 < t && _temp[a].D1 < t)
				{
					sum2 += (_temp[a].D0 - u) * (_temp[a].D0 - u); 
					sum2 += ((_temp[a].D0 + _temp[a].D1)/2 - u) * ((_temp[a].D0 + _temp[a].D1)/2 - u);
				}
			}
			sigma2[i] = sum2 / S[i];
			if(S[i] > max[k])
			{
				memcpy(H,H0,9 * sizeof(float));
				memcpy(InvH,InvH1,9 * sizeof(float));
				k = i;
				max[k] = S[i];
				listofRansac.clear();
				for(int a = 0;a <= n1;a++)
				{
					if(_temp[a].D0 < t && _temp[a].D1 < t && u < 1.0f)
					{
						_PAIR2 pair;
						pair.i = _temp[a].i;
						pair.j = _temp[a].j;
						pair.m = _temp[a].m;
						pair.n = _temp[a].n;
						listofRansac.push_back(pair);
					}

				}
				continue;
			}
			if(S[i] == max[k])
			{
				if(sigma2[i] < sigma2[k])
				{
					k = i;
					max[k] = S[i];
					memcpy(H,H0,9 * sizeof(float));
					memcpy(InvH,InvH1,9 * sizeof(float));
					listofRansac.clear();
					for(int a = 0;a <= n1;a++)
					{
						if(_temp[a].D0 < t && _temp[a].D1 < t && u < 1.0f)
						{
							_PAIR2 pair;
							pair.i = _temp[a].i;
							pair.j = _temp[a].j;
							pair.m = _temp[a].m;
							pair.n = _temp[a].n;
							listofRansac.push_back(pair);
						}
					}
					int ff = 0;
				}
			}
		}
	}
	if(max[k] <= 4)
	{
		AfxMessageBox(_T("匹配错误！"));
		delete[]S;
		delete[]sigma2;
		delete[]NUM;
		delete[]max;
		S = NULL;
		sigma2 = NULL;
		NUM = NULL;
		max = NULL;
		return 1;
	}
	delete[]S;
	delete[]sigma2;
	delete[]NUM;
	delete[]max;
	S = NULL;
	sigma2 = NULL;
	NUM = NULL;
	max = NULL;
	return 0;
}

bool CHarrisCorners::judgeFourPoints(_PAIR2* temp, short H1, short W1, short H2, short W2,  short n)
{
	if(     (fabs((float)temp[0].m - temp[1].m) <= 2 && fabs((float)temp[0].m - temp[2].m) <= 2 && fabs((float)temp[1].m - temp[2].m) <= 2) 
		 || (fabs((float)temp[0].m - temp[1].m) <= 2 && fabs((float)temp[0].m - temp[3].m) <= 2 && fabs((float)temp[1].m - temp[3].m) <= 2)      //若取得的四对匹配点中有三个点是共线的则重新随机取点；
		 || (fabs((float)temp[0].m - temp[2].m) <= 2 && fabs((float)temp[0].m - temp[3].m) <= 2 && fabs((float)temp[2].m - temp[3].m) <= 2)  
		 || (fabs((float)temp[2].m - temp[1].m) <= 2 && fabs((float)temp[1].m - temp[3].m) <= 2 && fabs((float)temp[2].m - temp[3].m) <= 2) 
		 || (fabs((float)temp[0].n - temp[1].n) <= 2 && fabs((float)temp[0].n - temp[2].n) <= 2 && fabs((float)temp[1].n - temp[2].n) <= 2) 
		 || (fabs((float)temp[0].n - temp[1].n) <= 2 && fabs((float)temp[0].n - temp[3].n) <= 2 && fabs((float)temp[1].n - temp[3].n) <= 2) 
		 || (fabs((float)temp[0].n - temp[2].n) <= 2 && fabs((float)temp[0].n - temp[3].n) <= 2 && fabs((float)temp[2].n - temp[3].n) <= 2)  
		 || (fabs((float)temp[2].n - temp[1].n) <= 2 && fabs((float)temp[1].n - temp[3].n) <= 2 && fabs((float)temp[2].n - temp[3].n) <= 2)   )
	{
			return 0;
	}
	float d1 = sqrt(float((temp[0].i - temp[1].i) * (temp[0].i - temp[1].i) + (temp[0].j - temp[1].j) * (temp[0].j - temp[1].j)));
	float d2 = sqrt(float((temp[2].i - temp[3].i) * (temp[2].i - temp[3].i) + (temp[2].j - temp[3].j) * (temp[2].j - temp[3].j)));
	float c1 = sqrt(float((temp[0].m - temp[1].m) * (temp[0].m - temp[1].m) + (temp[0].n - temp[1].n) * (temp[0].n - temp[1].n)));
	float c2 = sqrt(float((temp[2].m - temp[3].m) * (temp[2].m - temp[3].m) + (temp[2].n - temp[3].n) * (temp[2].j - temp[3].n)));
	float k = (d1 / d2) * (c2 / c1);
	if(k > 1.5f || k < 0.5f)
	{
		return 0;
	}

	_PAIR2 t;
	for(short i = 0;i < 3;i++)
	{
		for(short j = 0;j < 3 - i;j++)
		{
			if(temp[j].j > temp[j + 1].j)
			{
				t = temp[j];
				temp[j] = temp[j + 1];
				temp[j + 1] = t;
			}
		}
	}
	if(temp[0].i > temp[1].i)
	{
		t = temp[0];
		temp[0] = temp[1];
		temp[1] = t;
	}
	if(temp[2].i > temp[3].i)
	{
		t = temp[2];
		temp[2] = temp[3];
		temp[3] = t;
	}
	float l1 = sqrt(float((temp[0].i - temp[1].i) * (temp[0].i - temp[1].i) + (temp[0].j - temp[1].j) * (temp[0].j - temp[1].j)));
	float l2 = sqrt(float((temp[0].i - temp[2].i) * (temp[0].i - temp[2].i) + (temp[0].j - temp[2].j) * (temp[0].j - temp[2].j)));
	float l3 = sqrt(float((temp[2].i - temp[1].i) * (temp[2].i - temp[1].i) + (temp[2].j - temp[1].j) * (temp[2].j - temp[1].j)));
	float r = (l1 + l2 + l3) / 2;
	float s1 = sqrt(r * (r - l1) * (r - l2) * (r - l3));
	float l4 = sqrt(float((temp[3].i - temp[2].i) * (temp[3].i - temp[2].i) + (temp[3].j - temp[2].j) * (temp[3].j - temp[2].j)));
	float l5 = sqrt(float((temp[3].i - temp[1].i) * (temp[3].i - temp[1].i) + (temp[3].j - temp[1].j) * (temp[3].j - temp[1].j)));
	r = (l5 + l4 + l3) / 2;
	float s2 = sqrt(r * (r - l5) * (r - l4) * (r - l3));
	float s = s1 + s2;
	if(s < 10 || s > 0.75 * H1 * W1)
	{
		return 0;
	}

	for(short i = 0;i < 3;i++)
	{
		for(short j = 0;j < 3 - i;j++)
		{
			if(temp[j].n > temp[j + 1].n)
			{
				t = temp[j];
				temp[j] = temp[j + 1];
				temp[j + 1] = t;
			}
		}
	}
	if(temp[0].m > temp[1].m)
	{
		t = temp[0];
		temp[0] = temp[1];
		temp[1] = t;
	}

	if(temp[2].m > temp[3].m)
	{
		t = temp[2];
		temp[2] = temp[3];
		temp[3] = t;
	}
	l1 = sqrt(float((temp[0].m - temp[1].m) * (temp[0].m - temp[1].m) + (temp[0].n - temp[1].n) * (temp[0].n - temp[1].n)));
	l2 = sqrt(float((temp[0].m - temp[2].m) * (temp[0].m - temp[2].m) + (temp[0].n - temp[2].n) * (temp[0].n - temp[2].n)));
	l3 = sqrt(float((temp[2].m - temp[1].m) * (temp[2].m - temp[1].m) + (temp[2].n - temp[1].n) * (temp[2].n - temp[1].n)));
	r = (l1 + l2 + l3) / 2;
	s1 = sqrt(r * (r - l1) * (r - l2) * (r - l3));
	l4 = sqrt(float((temp[3].m - temp[2].m) * (temp[3].m - temp[2].m) + (temp[3].n - temp[2].n) * (temp[3].n - temp[2].n)));
	l5 = sqrt(float((temp[3].m - temp[1].m) * (temp[3].m - temp[1].m) + (temp[3].n - temp[1].n) * (temp[3].n - temp[1].n)));
	r = (l5 + l4 + l3) / 2;
	s2 = sqrt(r * (r - l5) * (r - l4) * (r - l3));
	s = s1 + s2;
	if(s < 10 || s > 0.75 * H2 * W2)
	{
		return 0;
	}
	return 1;
}
// 该函数根据得到的4对配准的特征点，计算投影变换系数。
void CHarrisCorners::getProjectionPara(_PAIR2* temp, float* pDbProjectionPara,bool flag)
{
	// pDbBMatrix中存放的是基准图象中特征点的坐标，
	float *pDbBMatrix;
	pDbBMatrix = new float[8];

	// pDbMatrix
	// 大小为8*8
	float *pDbMatrix;
	pDbMatrix = new float[8*8];

	// pDbInvMatrix为临时变量，存放的是pDbMatrix的逆
	// 大小为8*8
	float *pDbInvMatrix;
	pDbInvMatrix = new float[8*8];

	// 给矩阵赋值
	if(flag == 0)//计算H0时
	{
		for(short count = 0; count < 4; count++)
		{
			pDbBMatrix[2 * count] = float(temp[count].i);
			pDbBMatrix[2 * count + 1] = float(temp[count].j);	

			pDbMatrix[16 * count + 0] = float(temp[count].m);
			pDbMatrix[16 * count + 1] = float(temp[count].n);
			pDbMatrix[16 * count + 2] = 1;
			pDbMatrix[16 * count + 3] = 0;
			pDbMatrix[16 * count + 4] = 0;
			pDbMatrix[16 * count + 5] = 0;
			pDbMatrix[16 * count + 6] = -1 * (float((temp[count].m) * temp[count].i));
			pDbMatrix[16 * count + 7] = -1 * (float(temp[count].n * temp[count].i));
			pDbMatrix[8 * (2 * count + 1) + 0] =  0;
			pDbMatrix[8 * (2 * count + 1) + 1] =  0;
			pDbMatrix[8 * (2 * count + 1) + 2] =  0;
			pDbMatrix[8 * (2 * count + 1) + 3] =  float(temp[count].m);
			pDbMatrix[8 * (2 * count + 1) + 4] =  float(temp[count].n);
			pDbMatrix[8 * (2 * count + 1) + 5] =  1;
			pDbMatrix[8 * (2 * count + 1) + 6] =  -1 * (float((temp[count].m) * temp[count].j));
			pDbMatrix[8 * (2 * count + 1) + 7] =  -1 * (float(temp[count].n * temp[count].j));

			pDbInvMatrix[16 * count + 0] = float(temp[count].m);
			pDbInvMatrix[16 * count + 1] = float(temp[count].n);
			pDbInvMatrix[16 * count + 2] = 1;
			pDbInvMatrix[16 * count + 3] = 0;
			pDbInvMatrix[16 * count + 4] = 0;
			pDbInvMatrix[16 * count + 5] = 0;
			pDbInvMatrix[16 * count + 6] = -1 * (float((temp[count].m) * temp[count].i));
			pDbInvMatrix[16 * count + 7] = -1 * (float(temp[count].n * temp[count].i));
			pDbInvMatrix[8 * (2 * count + 1) + 0] =  0;
			pDbInvMatrix[8 * (2 * count + 1) + 1] =  0;
			pDbInvMatrix[8 * (2 * count + 1) + 2] =  0;
			pDbInvMatrix[8 * (2 * count + 1) + 3] =  float(temp[count].m);
			pDbInvMatrix[8 * (2 * count + 1) + 4] =  float(temp[count].n);
			pDbInvMatrix[8 * (2 * count + 1) + 5] =  1;
			pDbInvMatrix[8 * (2 * count + 1) + 6] =  -1 * (float((temp[count].m) * temp[count].j));
			pDbInvMatrix[8 * (2 * count + 1) + 7] =  -1 * (float(temp[count].n * temp[count].j));
		}
	}
	if(flag == 1)//计算H1时
	{
		for(short count = 0; count < 4; count++)
		{
			pDbBMatrix[2 * count] = float(temp[count].m);
			pDbBMatrix[2 * count + 1] = float(temp[count].n);	

			pDbMatrix[16 * count + 0] = float(temp[count].i);
			pDbMatrix[16 * count + 1] = float(temp[count].j);
			pDbMatrix[16 * count + 2] = 1;
			pDbMatrix[16 * count + 3] = 0;
			pDbMatrix[16 * count + 4] = 0;
			pDbMatrix[16 * count + 5] = 0;
			pDbMatrix[16 * count + 6] = -1 * (float((temp[count].m) * temp[count].i));
			pDbMatrix[16 * count + 7] = -1 * (float(temp[count].j * (temp[count].m)));
			pDbMatrix[8 * (2 * count + 1) + 0] =  0;
			pDbMatrix[8 * (2 * count + 1) + 1] =  0;
			pDbMatrix[8 * (2 * count + 1) + 2] =  0;
			pDbMatrix[8 * (2 * count + 1) + 3] =  float(temp[count].i);
			pDbMatrix[8 * (2 * count + 1) + 4] =  float(temp[count].j);
			pDbMatrix[8 * (2 * count + 1) + 5] =  1;
			pDbMatrix[8 * (2 * count + 1) + 6] =  -1 * (float(temp[count].i * temp[count].n));
			pDbMatrix[8 * (2 * count + 1) + 7] =  -1 * (float(temp[count].n * temp[count].j));

			pDbInvMatrix[16 * count + 0] = float(temp[count].i);
			pDbInvMatrix[16 * count + 1] = float(temp[count].j);
			pDbInvMatrix[16 * count + 2] = 1;
			pDbInvMatrix[16 * count + 3] = 0;
			pDbInvMatrix[16 * count + 4] = 0;
			pDbInvMatrix[16 * count + 5] = 0;
			pDbInvMatrix[16 * count + 6] = -1 * (float((temp[count].m) * temp[count].i));
			pDbInvMatrix[16 * count + 7] = -1 * (float(temp[count].j * (temp[count].m)));
			pDbInvMatrix[8 * (2 * count + 1) + 0] =  0;
			pDbInvMatrix[8 * (2 * count + 1) + 1] =  0;
			pDbInvMatrix[8 * (2 * count + 1) + 2] =  0;
			pDbInvMatrix[8 * (2 * count + 1) + 3] =  float(temp[count].i);
			pDbInvMatrix[8 * (2 * count + 1) + 4] =  float(temp[count].j);
			pDbInvMatrix[8 * (2 * count + 1) + 5] =  1;
			pDbInvMatrix[8 * (2 * count + 1) + 6] =  -1 * (float(temp[count].i * temp[count].n));
			pDbInvMatrix[8 * (2 * count + 1) + 7] =  -1 * (float(temp[count].n * temp[count].j));
		}
	}
	// 计算pDbInvMatrix的逆矩阵,存放在pDbInvMatrix
	calInvMatrix(pDbInvMatrix,8);
	// 计算投影变换系数
	calMatProduct(pDbInvMatrix,pDbBMatrix,pDbProjectionPara,8,1,8);
	// 释放内存
	delete[]pDbBMatrix;
	delete[]pDbMatrix;
	delete[]pDbInvMatrix;
	pDbBMatrix   = NULL;
	pDbMatrix    = NULL;
	pDbInvMatrix = NULL;
}
// 该函数计算两个矩阵的相乘，然后将相乘的结果存放在pDbDest中。
void CHarrisCorners::calMatProduct(float* pDbSrc1, float* pDbSrc2, float* pDbDest,short y, short x, short z)
{

	for(short i=0;i<y;i++)
	{
		for(short j=0;j<x;j++)
		{
			pDbDest[i*x + j] = 0;
			for(int m=0;m<z;m++)
			{
				pDbDest[i*x + j] += pDbSrc1[i*z + m]*pDbSrc2[m*x + j];
			}
		}
	}
}
// 该函数计算矩阵pDbSrc的逆矩阵，其中pDbSrc的大小为nLen*nLen
BOOL CHarrisCorners::calInvMatrix(float* pDbSrc, short nLen)//calInvMatrix(pDbInvMatrix, 8);
{
	short *is,*js,i,j,k;
	float d,p;
	is = new short[nLen];
	memset(is,0,nLen * sizeof(short));
	js = new short[nLen];
	memset(js,0,nLen * sizeof(short));
	for(k=0;k<nLen;k++)
	{
		d=0.0f;
		for(i=k;i<nLen;i++)
		{
			for(j=k;j<nLen;j++)
			{
				p=fabs(pDbSrc[i*nLen + j]);
				if(p>d)
				{
					d     = p; 
					is[k] = i;
					js[k] = j;
				}
			}
		}
		if(d+1.0==1.0)
		{
			delete[]is;
			delete[]js;
			is = NULL;
			js = NULL;
			return FALSE;
		}
		if(is[k] != k)
		{
			for(j=0;j<nLen;j++)
			{
				p = pDbSrc[k*nLen + j];
				pDbSrc[k*nLen + j] = pDbSrc[(is[k]*nLen) + j];
				pDbSrc[(is[k])*nLen + j] = p;
			}
		}
		if(js[k] != k)
		{
			for(i=0; i<nLen; i++)
			{
				p = pDbSrc[i*nLen + k];
				pDbSrc[i*nLen + k] = pDbSrc[i*nLen + (js[k])];
				pDbSrc[i*nLen + (js[k])] = p;
			}
		}
		pDbSrc[k*nLen + k]=1.0f / pDbSrc[k*nLen + k];
		for(j=0; j<nLen; j++)
		{
			if(j != k)
			{
				pDbSrc[k*nLen + j]*=pDbSrc[k*nLen + k];
			}
		}
		for(i=0; i<nLen; i++)
		{
			if(i != k)
			{
				for(j=0; j<nLen; j++)
				{
					if(j!=k)
					{
						pDbSrc[i*nLen + j] -= pDbSrc[i*nLen + k]*pDbSrc[k*nLen + j];
					}
				}
			}
		}
		for(i=0; i<nLen; i++)
		{
			if(i != k)
			{
				pDbSrc[i*nLen + k] *= -pDbSrc[k*nLen + k];
			}
		}
	}
	for(k=nLen-1; k>=0; k--)
	{
		if(js[k] != k)
		{
			for(j=0; j<nLen; j++)
			{
				p = pDbSrc[k*nLen + j];
				pDbSrc[k*nLen + j] = pDbSrc[(js[k])*nLen + j];
				pDbSrc[(js[k])*nLen + j] = p;
			}
		}
		if(is[k] != k)
		{
			for(i=0; i<nLen; i++)
			{
				p = pDbSrc[i*nLen + k];
				pDbSrc[i*nLen + k] = pDbSrc[i*nLen +(is[k])];
				pDbSrc[i*nLen + (is[k])] = p;
			}
		}
	}
	delete[]is;
	delete[]js;
	is = NULL;
	js = NULL;
	return TRUE;	
}
//L-M优化
void CHarrisCorners::L_MOptimize(_PAIR2* temp, int temp_N,bool flag) 
{
	//最大迭代次数取为10
	float epsilon = 0.7f * (temp_N);
	float u = 0.01f;

	float F[11];
	float m0[11];
	float m1[11];
	float m2[11];
	float m3[11];
	float m4[11];
	float m5[11];;
	float m6[11];
	float m7[11];
	if(flag == 0)
	{
		for(int k = 0;k < 11;k++)
		{
			if( k == 0)
			{
				m0[k] = H[0] * 10000;
				m1[k] = H[1] * 10000;
				m2[k] = H[2];
				m3[k] = H[3] * 10000;
				m4[k] = H[4] * 10000;
				m5[k] = H[5];
				m6[k] = H[6] * 10000;
				m7[k] = H[7] * 10000;
			}

			for(int a = 0;a < temp_N;a++)
			{
				float x = (m0[k] * (temp[a].m) / 10000 + m1[k] * temp[a].n / 10000 + m2[k]) / (m6[k] * (temp[a].m) / 10000 + m7[k] * temp[a].n / 10000 + 1);
				float y = (m3[k] * (temp[a].m) / 10000 + m4[k] * temp[a].n / 10000 + m5[k]) / (m6[k] * (temp[a].m) / 10000 + m7[k] * temp[a].n / 10000 + 1);
				float D = sqrt((x - temp[a].i) * (x - temp[a].i) + (y - temp[a].j) * (y - temp[a].j));
				F[k] += D; 
			}

			if(F[k] > epsilon)
			{
				if(k != 0)
				{
					if(F[k] < F[k - 1])
					{
						u = u / 10;
						if(k > 9)
						{
							H[0] = m0[k] / 10000;
							H[1] = m1[k] / 10000;
							H[2] = m2[k];
							H[3] = m3[k] / 10000;
							H[4] = m4[k] / 10000;
							H[5] = m5[k] ;
							H[6] = m6[k] / 10000;
							H[7] = m7[k] / 10000;
							break;
						}
						else
						{
							float* J = new float[8 * (temp_N)];
							memset(J,0,8 * (temp_N) * sizeof(float));
							float* JT = new float[(temp_N) * 8];
							memset(JT,0,8 * (temp_N) * sizeof(float));
							float* f = new float[temp_N];
							memset(f,0,(temp_N) * sizeof(float));
							for(int b = 0;b < temp_N;b++)
							{
								float errorx = (m0[k] * ((temp[b].m)/10000) + m1[k] * (temp[b].n/10000) + m2[k]) / (m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1) - temp[b].i;
								float errory = (m3[k] * ((temp[b].m)/10000) + m4[k] * (temp[b].n/10000) + m5[k]) / (m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1) - temp[b].j;
								float error = sqrt(sqrt(errorx * errorx + errory * errory));
								float deverror = 1 / (4 * error * error * error);
								J[8 * b + 0] = 2 * deverror * errorx * ((temp[b].m)/10000) / (m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1);
								J[8 * b + 1] = 2 * deverror * errorx * (temp[b].n/10000) / (m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1);
								J[8 * b + 2] = 2 * deverror * errorx / (m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1);
								J[8 * b + 3] = 2 * deverror * errory * ((temp[b].m)/10000) / (m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1);
								J[8 * b + 4] = 2 * deverror * errory * (temp[b].n/10000) / (m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1);
								J[8 * b + 5] = 2 * deverror * errory / (m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1);
								J[8 * b + 6] = 2 * deverror * ((temp[b].m)/10000) * (-1 * errorx * (m0[k] * ((temp[b].m)/10000) + m1[k] * (temp[b].n/10000) + m2[k]) - errory * (m3[k] * ((temp[b].m)/10000) + m4[k] * (temp[b].n/10000) + m5[k]))
												   / ((m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1) * (m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1));
								J[8 * b + 7] = 2 * deverror * (temp[b].n/10000) * (-1 * errorx * (m0[k] * ((temp[b].m)/10000) + m1[k] * (temp[b].n/10000) + m2[k]) - errory * (m3[k] * ((temp[b].m)/10000) + m4[k] * (temp[b].n/10000) + m5[k]))
												   / ((m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1) * (m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1));

								f[b] = sqrt(sqrt(errorx * errorx + errory * errory));

							}
							for(int b = 0;b < temp_N;b++)
							{
								for(int a = 0;a < 8;a++)
								{
									JT[8 * a + b] = J[8 * b + a];
								}
							}
							//保存JT * J
							float* pDbTemp1 = new float[8 * 8];
							memset(pDbTemp1,0,8 * 8 * sizeof(float));
							//保存u * I
							float* E =new float[8 *8];
							memset(E,0,8 * 8 * sizeof(float));
							//保存JT * J + E的逆
							float* pDbInvTemp = new float[8 * 8];
							memset(J,0,8 * 8 * sizeof(float));
							for(int b = 0;b < 8;b++)
							{
								for(int a = 0;a < 8;a++)
								{
									E[8 * b + a] = 0;
									if(a == b)
									{
										E[8 * b + a] = u;
									}
								}
							}
							// 计算JT * J,存放在pDbTemp1
							calMatProduct(JT,J,pDbTemp1,8,8,temp_N);
							//计算JT * J + E,存放在pDbTemp1
							for(int b = 0;b < 8;b++)
							{
								for(int a = 0;a < 8;a++)
								{
									pDbTemp1[8 * b + a] = pDbTemp1[8 * b + a] + E[8 * b + a];
								}
							}
							memcpy(pDbInvTemp,pDbTemp1,64 * sizeof(float));
							// 计算pDbTemp1的逆矩阵,存放在pDbInvTemp
							calInvMatrix(pDbInvTemp,8);
							// 计算pDbInvTemp * JT,存放在pDbTemp2
							float* pDbTemp2 = new float[8 * (temp_N)];
							memset(pDbTemp2,0,8 * (temp_N) * sizeof(float));
							calMatProduct(pDbInvTemp,JT,pDbTemp2,8,temp_N,8);
							// 计算pDbTemp2 * f,存放在pDbTemp3
							float* pDbTemp3 = new float[8 * 1];
							memset(pDbTemp3,0,8 * sizeof(float));
							calMatProduct(pDbTemp2,f,pDbTemp3,8,1,temp_N);
							m0[k+1] = m0[k] - pDbTemp3[0];
							m1[k+1] = m1[k] - pDbTemp3[1];
							m2[k+1] = m2[k] - pDbTemp3[2];
							m3[k+1] = m3[k] - pDbTemp3[3];
							m4[k+1] = m4[k] - pDbTemp3[4];
							m5[k+1] = m5[k] - pDbTemp3[5];
							m6[k+1] = m6[k] - pDbTemp3[6];
							m7[k+1] = m7[k] - pDbTemp3[7];
							delete[]J;
							delete[]JT;
							delete[]f;
							delete[]pDbTemp1;
							delete[]pDbInvTemp;
							delete[]pDbTemp2;
							delete[]pDbTemp3;
							J  = NULL;
							JT = NULL;
							f  = NULL;
							pDbTemp1 = NULL;
							pDbTemp2 = NULL;
							pDbTemp3 = NULL;
							pDbInvTemp = NULL;
							continue;
						}
					}
					else  //F[k] >= F[k - 1]的情况
					{
						m0[k] = m0[k-1];
						m1[k] = m1[k-1];
						m2[k] = m2[k-1];
						m3[k] = m3[k-1];
						m4[k] = m4[k-1];
						m5[k] = m5[k-1];
						m6[k] = m6[k-1];
						m7[k] = m7[k-1];
						u = 10 * u;
					}

				}
				else // k = 0的情况
				{
					u = u / 10;
					float* J = new float[8 * (temp_N)];
					memset(J,0,8 * (temp_N) * sizeof(float));
					float* JT = new float[(temp_N) * 8];
					memset(JT,0,8 * (temp_N) * sizeof(float));
					float* f = new float[temp_N];
					memset(f,0,(temp_N) * sizeof(float));
					for(int b = 0;b < temp_N;b++)
					{
						float errorx = (m0[k] * ((temp[b].m)/10000) + m1[k] * (temp[b].n/10000) + m2[k]) / (m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1) - temp[b].i;
						float errory = (m3[k] * ((temp[b].m)/10000) + m4[k] * (temp[b].n/10000) + m5[k]) / (m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1) - temp[b].j;
						float error = sqrt(sqrt(errorx * errorx + errory * errory));
						float deverror = 1 / (4 * error * error * error);
						J[8 * b + 0] = 2 * deverror * errorx * ((temp[b].m)/10000) / (m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1);
						J[8 * b + 1] = 2 * deverror * errorx * (temp[b].n/10000) / (m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1);
						J[8 * b + 2] = 2 * deverror * errorx / (m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1);
						J[8 * b + 3] = 2 * deverror * errory * ((temp[b].m)/10000) / (m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1);
						J[8 * b + 4] = 2 * deverror * errory * (temp[b].n/10000) / (m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1);
						J[8 * b + 5] = 2 * deverror * errory / (m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1);
						J[8 * b + 6] = 2 * deverror * ((temp[b].m)/10000) * (-1 * errorx * (m0[k] * ((temp[b].m)/10000) + m1[k] * (temp[b].n/10000) + m2[k]) - errory * (m3[k] * ((temp[b].m)/10000) + m4[k] * (temp[b].n/10000) + m5[k]))
											   / ((m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1) * (m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1));
						J[8 * b + 7] = 2 * deverror * (temp[b].n/10000) * (-1 * errorx * (m0[k] * ((temp[b].m)/10000) + m1[k] * (temp[b].n/10000) + m2[k]) - errory * (m3[k] * ((temp[b].m)/10000) + m4[k] * (temp[b].n/10000) + m5[k]))
												   / ((m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1) * (m6[k] * ((temp[b].m)/10000) + m7[k] * (temp[b].n/10000) + 1));

						f[b] = sqrt(sqrt(errorx * errorx + errory * errory));
					}

					for(int b = 0;b < temp_N;b++)
					{
						for(int a = 0;a < 8;a++)
						{
							JT[8 * a + b] = J[8 * b + a];
						}
					}
					//保存JT * J
					float* pDbTemp1 = new float[8 * 8];
					memset(pDbTemp1,0,8 * 8 * sizeof(float));
					//保存u * I
					float* E =new float[8 *8];
					memset(E,0,8 * 8 * sizeof(float));
					//保存JT * J + E的逆
					float* pDbInvTemp = new float[8 * 8];
					memset(J,0,8 * 8 * sizeof(float));
					for(int b = 0;b < 8;b++)
					{
						for(int a = 0;a < 8;a++)
						{
							E[8 * b + a] = 0;
							if(a == b)
							{
								E[8 * b + a] = u;
							}
						}
					}
					// 计算JT * J,存放在pDbTemp1
					calMatProduct(JT,J,pDbTemp1,8,8,temp_N);
					//计算JT * J + E,存放在pDbTemp1
					for(int b = 0;b < 8;b++)
					{
						for(int a = 0;a < 8;a++)
						{
							pDbTemp1[8 * b + a] = pDbTemp1[8 * b + a] + E[8 * b + a];
						}
					}
					memcpy(pDbInvTemp,pDbTemp1,64 * sizeof(float));
					// 计算pDbTemp1的逆矩阵,存放在pDbInvTemp
					calInvMatrix(pDbInvTemp,8);
					// 计算pDbInvTemp * JT,存放在pDbTemp2
					float* pDbTemp2 = new float[8 * (temp_N)];
					memset(pDbTemp2,0,8 * (temp_N) * sizeof(float));
					calMatProduct(pDbInvTemp,JT,pDbTemp2,8,temp_N,8);
					// 计算pDbTemp2 * f,存放在pDbTemp3
					float* pDbTemp3 = new float[8 * 1];
					memset(pDbTemp3,0,8 * sizeof(float));
					calMatProduct(pDbTemp2,f,pDbTemp3,8,1,temp_N);
					m0[k+1] = m0[k] - pDbTemp3[0];
					m1[k+1] = m1[k] - pDbTemp3[1];
					m2[k+1] = m2[k] - pDbTemp3[2];
					m3[k+1] = m3[k] - pDbTemp3[3];
					m4[k+1] = m4[k] - pDbTemp3[4];
					m5[k+1] = m5[k] - pDbTemp3[5];
					m6[k+1] = m6[k] - pDbTemp3[6];
					m7[k+1] = m7[k] - pDbTemp3[7];
					delete[]J;
					delete[]JT;
					delete[]f;
					delete[]pDbTemp1;
					delete[]pDbInvTemp;
					delete[]pDbTemp2;
					delete[]pDbTemp3;
					J = NULL;
					JT = NULL;
					f = NULL;
					pDbTemp1 = NULL;
					pDbTemp2 = NULL;
					pDbTemp3 = NULL;
					pDbInvTemp = NULL;
					continue;
				}
			}
			else // F[k] <= epsilon的情况
			{
				H[0] = m0[k] / 10000;
				H[1] = m1[k] / 10000;
				H[2] = m2[k];
				H[3] = m3[k] / 10000;
				H[4] = m4[k] / 10000;
				H[5] = m5[k];
				H[6] = m6[k] / 10000;
				H[7] = m7[k] / 10000;
				break;
			}
		}
	}
	if(flag == 1)
	{
		for(int k = 0;k < 11;k++)
		{
			if( k == 0)
			{
				m0[k] = InvH[0] * 10000;
				m1[k] = InvH[1] * 10000;
				m2[k] = InvH[2];
				m3[k] = InvH[3] * 10000;
				m4[k] = InvH[4] * 10000;
				m5[k] = InvH[5];
				m6[k] = InvH[6] * 10000;
				m7[k] = InvH[7] * 10000;
			}

			for(int a = 0;a <= temp_N;a++)
			{
				float x = (m0[k] * (temp[a].i) / 10000 + m1[k] * temp[a].j / 10000 + m2[k]) / (m6[k] * (temp[a].i) / 10000 + m7[k] * temp[a].j / 10000 + 1);
				float y = (m3[k] * (temp[a].i) / 10000 + m4[k] * temp[a].j / 10000 + m5[k]) / (m6[k] * (temp[a].i) / 10000 + m7[k] * temp[a].j / 10000 + 1);
				float D = sqrt((x - (temp[a].m)) * (x - (temp[a].m)) + (y - temp[a].n) * (y - temp[a].n));
				F[k] += D; 
			}

			if(F[k] > epsilon)
			{
				if(k != 0)
				{
					if(F[k] < F[k - 1])
					{
						u = u / 10;
						if(k > 9)
						{
							InvH[0] = m0[k] / 10000;
							InvH[1] = m1[k] / 10000;
							InvH[2] = m2[k];
							InvH[3] = m3[k] / 10000;
							InvH[4] = m4[k] / 10000;
							InvH[5] = m5[k] ;
							InvH[6] = m6[k] / 10000;
							InvH[7] = m7[k] / 10000;
							break;
						}
						else
						{
							float* J = new float[8 * (temp_N)];
							memset(J,0,8 * (temp_N) * sizeof(float));
							float* JT = new float[(temp_N) * 8];
							memset(JT,0,8 * (temp_N) * sizeof(float));
							float* f = new float[temp_N];
							memset(f,0,(temp_N) * sizeof(float));
							for(int b = 0;b < temp_N;b++)
							{
								float errorx = (m0[k] * ((temp[b].i)/10000) + m1[k] * (temp[b].j/10000) + m2[k]) / (m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1) - (temp[b].m);
								float errory = (m3[k] * ((temp[b].i)/10000) + m4[k] * (temp[b].j/10000) + m5[k]) / (m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1) - temp[b].n;
								float error = sqrt(sqrt(errorx * errorx + errory * errory));
								float deverror = 1 / (4 * error * error * error);
								J[8 * b + 0] = 2 * deverror * errorx * ((temp[b].i)/10000) / (m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1);
								J[8 * b + 1] = 2 * deverror * errorx * (temp[b].j/10000) / (m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1);
								J[8 * b + 2] = 2 * deverror * errorx / (m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1);
								J[8 * b + 3] = 2 * deverror * errory * ((temp[b].i)/10000) / (m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1);
								J[8 * b + 4] = 2 * deverror * errory * (temp[b].j/10000) / (m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1);
								J[8 * b + 5] = 2 * deverror * errory / (m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1);
								J[8 * b + 6] = 2 * deverror * ((temp[b].i)/10000) * (-1 * errorx * (m0[k] * ((temp[b].i)/10000) + m1[k] * (temp[b].j/10000) + m2[k]) - errory * (m3[k] * ((temp[b].i)/10000) + m4[k] * (temp[b].j/10000) + m5[k]))
												   / ((m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1) * (m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1));
								J[8 * b + 7] = 2 * deverror * (temp[b].j/10000) * (-1 * errorx * (m0[k] * ((temp[b].i)/10000) + m1[k] * (temp[b].j/10000) + m2[k]) - errory * (m3[k] * ((temp[b].i)/10000) + m4[k] * (temp[b].j/10000) + m5[k]))
												   / ((m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1) * (m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1));

								f[b] = sqrt(sqrt(errorx * errorx + errory * errory));

							}
							for(int b = 0;b < temp_N;b++)
							{
								for(int a = 0;a < 8;a++)
								{
									JT[8 * a + b] = J[8 * b + a];
								}
							}
							//保存JT * J
							float* pDbTemp1 = new float[8 * 8];
							memset(pDbTemp1,0,8 * 8 * sizeof(float));
							//保存u * I
							float* E =new float[8 *8];
							memset(E,0,8 * 8 * sizeof(float));
							//保存JT * J + E的逆
							float* pDbInvTemp = new float[8 * 8];
							memset(J,0,8 * 8 * sizeof(float));
							for(int b = 0;b < 8;b++)
							{
								for(int a = 0;a < 8;a++)
								{
									E[8 * b + a] = 0;
									if(a == b)
									{
										E[8 * b + a] = u;
									}
								}
							}
							// 计算JT * J,存放在pDbTemp1
							calMatProduct(JT,J,pDbTemp1,8,8,temp_N);
							//计算JT * J + E,存放在pDbTemp1
							for(int b = 0;b < 8;b++)
							{
								for(int a = 0;a < 8;a++)
								{
									pDbTemp1[8 * b + a] = pDbTemp1[8 * b + a] + E[8 * b + a];
								}
							}
							memcpy(pDbInvTemp,pDbTemp1,64 * sizeof(float));
							// 计算pDbTemp1的逆矩阵,存放在pDbInvTemp
							calInvMatrix(pDbInvTemp,8);
							// 计算pDbInvTemp * JT,存放在pDbTemp2
							float* pDbTemp2 = new float[8 * (temp_N)];
							memset(pDbTemp2,0,8 * (temp_N) * sizeof(float));
							calMatProduct(pDbInvTemp,JT,pDbTemp2,8,temp_N,8);
							// 计算pDbTemp2 * f,存放在pDbTemp3
							float* pDbTemp3 = new float[8 * 1];
							memset(pDbTemp3,0,8 * sizeof(float));
							calMatProduct(pDbTemp2,f,pDbTemp3,8,1,temp_N);
							m0[k+1] = m0[k] - pDbTemp3[0];
							m1[k+1] = m1[k] - pDbTemp3[1];
							m2[k+1] = m2[k] - pDbTemp3[2];
							m3[k+1] = m3[k] - pDbTemp3[3];
							m4[k+1] = m4[k] - pDbTemp3[4];
							m5[k+1] = m5[k] - pDbTemp3[5];
							m6[k+1] = m6[k] - pDbTemp3[6];
							m7[k+1] = m7[k] - pDbTemp3[7];
							delete[]J;
							delete[]JT;
							delete[]f;
							delete[]pDbTemp1;
							delete[]pDbInvTemp;
							delete[]pDbTemp2;
							delete[]pDbTemp3;
							J  = NULL;
							JT = NULL;
							f  = NULL;
							pDbTemp1 = NULL;
							pDbTemp2 = NULL;
							pDbTemp3 = NULL;
							pDbInvTemp = NULL;
							continue;
						}
					}
					else  //F[k] >= F[k - 1]的情况
					{
						m0[k] = m0[k-1];
						m1[k] = m1[k-1];
						m2[k] = m2[k-1];
						m3[k] = m3[k-1];
						m4[k] = m4[k-1];
						m5[k] = m5[k-1];
						m6[k] = m6[k-1];
						m7[k] = m7[k-1];
						u = 10 * u;
					}

				}
				else // k = 0的情况
				{
					u = u / 10;
					float* J = new float[8 * (temp_N)];
					memset(J,0,8 * (temp_N) * sizeof(float));
					float* JT = new float[(temp_N) * 8];
					memset(JT,0,8 * (temp_N) * sizeof(float));
					float* f = new float[temp_N];
					memset(f,0,(temp_N) * sizeof(float));
					for(int b = 0;b < temp_N;b++)
					{
						float errorx = (m0[k] * ((temp[b].i)/10000) + m1[k] * (temp[b].j/10000) + m2[k]) / (m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1) - (temp[b].m);
						float errory = (m3[k] * ((temp[b].i)/10000) + m4[k] * (temp[b].j/10000) + m5[k]) / (m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1) - temp[b].n;
						float error = sqrt(sqrt(errorx * errorx + errory * errory));
						float deverror = 1 / (4 * error * error * error);
						J[8 * b + 0] = 2 * deverror * errorx * ((temp[b].i)/10000) / (m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1);
						J[8 * b + 1] = 2 * deverror * errorx * (temp[b].j/10000) / (m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1);
						J[8 * b + 2] = 2 * deverror * errorx / (m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1);
						J[8 * b + 3] = 2 * deverror * errory * ((temp[b].i)/10000) / (m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1);
						J[8 * b + 4] = 2 * deverror * errory * (temp[b].j/10000) / (m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1);
						J[8 * b + 5] = 2 * deverror * errory / (m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1);
						J[8 * b + 6] = 2 * deverror * ((temp[b].i)/10000) * (-1 * errorx * (m0[k] * ((temp[b].i)/10000) + m1[k] * (temp[b].j/10000) + m2[k]) - errory * (m3[k] * ((temp[b].i)/10000) + m4[k] * (temp[b].j/10000) + m5[k]))
											   / ((m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1) * (m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1));
						J[8 * b + 7] = 2 * deverror * (temp[b].j/10000) * (-1 * errorx * (m0[k] * ((temp[b].i)/10000) + m1[k] * (temp[b].j/10000) + m2[k]) - errory * (m3[k] * ((temp[b].i)/10000) + m4[k] * (temp[b].j/10000) + m5[k]))
												   / ((m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1) * (m6[k] * ((temp[b].i)/10000) + m7[k] * (temp[b].j/10000) + 1));

						f[b] = sqrt(sqrt(errorx * errorx + errory * errory));
					}

					for(int b = 0;b < temp_N;b++)
					{
						for(int a = 0;a < 8;a++)
						{
							JT[8 * a + b] = J[8 * b + a];
						}
					}
					//保存JT * J
					float* pDbTemp1 = new float[8 * 8];
					memset(pDbTemp1,0,8 * 8 * sizeof(float));
					//保存u * I
					float* E =new float[8 *8];
					memset(E,0,8 * 8 * sizeof(float));
					//保存JT * J + E的逆
					float* pDbInvTemp = new float[8 * 8];
					memset(J,0,8 * 8 * sizeof(float));
					for(int b = 0;b < 8;b++)
					{
						for(int a = 0;a < 8;a++)
						{
							E[8 * b + a] = 0;
							if(a == b)
							{
								E[8 * b + a] = u;
							}
						}
					}
					// 计算JT * J,存放在pDbTemp1
					calMatProduct(JT,J,pDbTemp1,8,8,temp_N);
					//计算JT * J + E,存放在pDbTemp1
					for(int b = 0;b < 8;b++)
					{
						for(int a = 0;a < 8;a++)
						{
							pDbTemp1[8 * b + a] = pDbTemp1[8 * b + a] + E[8 * b + a];
						}
					}
					memcpy(pDbInvTemp,pDbTemp1,64 * sizeof(float));
					// 计算pDbTemp1的逆矩阵,存放在pDbInvTemp
					calInvMatrix(pDbInvTemp,8);
					// 计算pDbInvTemp * JT,存放在pDbTemp2
					float* pDbTemp2 = new float[8 * (temp_N)];
					memset(pDbTemp2,0,8 * (temp_N) * sizeof(float));
					calMatProduct(pDbInvTemp,JT,pDbTemp2,8,temp_N,8);
					// 计算pDbTemp2 * f,存放在pDbTemp3
					float* pDbTemp3 = new float[8 * 1];
					memset(pDbTemp3,0,8 * sizeof(float));
					calMatProduct(pDbTemp2,f,pDbTemp3,8,1,temp_N);
					m0[k+1] = m0[k] - pDbTemp3[0];
					m1[k+1] = m1[k] - pDbTemp3[1];
					m2[k+1] = m2[k] - pDbTemp3[2];
					m3[k+1] = m3[k] - pDbTemp3[3];
					m4[k+1] = m4[k] - pDbTemp3[4];
					m5[k+1] = m5[k] - pDbTemp3[5];
					m6[k+1] = m6[k] - pDbTemp3[6];
					m7[k+1] = m7[k] - pDbTemp3[7];
					delete[]J;
					delete[]JT;
					delete[]f;
					delete[]pDbTemp1;
					delete[]pDbInvTemp;
					delete[]pDbTemp2;
					delete[]pDbTemp3;
					J = NULL;
					JT = NULL;
					f = NULL;
					pDbTemp1 = NULL;
					pDbTemp2 = NULL;
					pDbTemp3 = NULL;
					pDbInvTemp = NULL;
					continue;
				}
			}
			else // F[k] <= epsilon的情况
			{
				InvH[0] = m0[k] / 10000;
				InvH[1] = m1[k] / 10000;
				InvH[2] = m2[k];
				InvH[3] = m3[k] / 10000;
				InvH[4] = m4[k] / 10000;
				InvH[5] = m5[k];
				InvH[6] = m6[k] / 10000;
				InvH[7] = m7[k] / 10000;
				break;
			}
		}
	}
}
/******************************************************************************************************/
/*************************************************图像融合拼接*****************************************/
/******************************************************************************************************/

void Trans(float &s , float u)
{
	float aa = fabs(u);
	if(aa < 1)
	{
		s = 1 - 2 * aa * aa + aa * aa * aa; 
	}
	else if(aa <= 2 && aa >= 1)
	{
		s = 4 - 8 * aa + 5 * aa * aa - aa * aa * aa;
	}
	else
	{
		s = 0;
	}
}
void GetS(float* s , float v ,float b)
{
	for(int i = 0;i < 8;i++)
	{
		if(i == 0)
		{
				float u1 = v + 1;
				float t1;
				Trans(t1,u1);
				s[0] = t1;
		}
		if(i == 1)
		{
				float u2 = v;
				float t2;
				Trans(t2,u2); 
				s[1] = t2;
		}
		if(i == 2)
		{
				float u3 = 1 - v;
				float t3;
				Trans(t3,u3);
				s[2] = t3;
		}
		if(i == 3)
		{
				float u4 = 2 - v;
				float t4;
				Trans(t4,u4);
				s[3] = t4;
		}
		if(i == 4)
		{
				float u5 = 1 + b;
				float t5;
				Trans(t5,u5);
				s[4] = t5;
		}
		if(i == 5)
		{
				float u6 = b;
				float t6;
				Trans(t6,u6);
				s[5] = t6;
		}
		if(i == 6)
		{
				float u7 = 1 - b;
				float t7;
				Trans(t7,u7);
				s[6] = t7;
		}
		if(i == 7)
		{
				float u8 = 2 - b;
				float t8;
				Trans(t8,u8);
				s[7] = t8;
		}
	}
}

void CHarrisCorners::imageFusion(CImage& pNewImage1,CImage& pNewImage2,CImage &m_newImage)
{
//	// TODO: Add your control notification handler code here
	BYTE* pArray1 = (BYTE*)pNewImage1.GetBits(); //获得数据
	int pitch1 = pNewImage1.GetPitch();
	int bitCount1 = pNewImage1.GetBPP() / 8;  //每像素位数
	BYTE* pArray2 = (BYTE*)pNewImage2.GetBits(); //获得数据
	int pitch2 = pNewImage2.GetPitch();
	int bitCount2 = pNewImage2.GetBPP() / 8;  //每像素位数
	//新图像的高宽
	short newHeight = 0;
	short newWidth = 0;
	//新图像left,right ,top,bottom在投影变换坐标系内的坐标值
	short left = 0;
	short right = 0;
	short top = 0;
	short bottom = 0;
	// 计算图象经过投影变换后新图像的高宽和新图像left、right 、top、bottom在投影变换坐标系内的坐标值
	getAftProDim(&left,&right,&top,&bottom,&newHeight,&newWidth,nHeight1,nWidth1,nHeight2,nWidth2);

	// 对新图象进行赋值
	m_newImage.Create(newWidth,newHeight,24);
	BYTE* nArray = (BYTE *) m_newImage.GetBits();
	int npitch = m_newImage.GetPitch();
	int bitCount = m_newImage.GetBPP() / 8;  //每像素位数

	int dis = 30;

	//计算补偿值
	float dif_R = 0;
	float dif_G = 0;
	float dif_B = 0;
	int num = 0;
	float dif_Y = 0;
	for(short j = 0;j < newHeight;j++)
	{
		for(short i = 0;i < newWidth;i++)
		{
			if(i < nWidth1 - left && i >= -1 * left && j < nHeight1 - top && j >= -1 * top)
			{
				float a = H[6] * (i + left) - H[0];
				float b = H[7] * (i + left) - H[1];
				float c = H[2] - (i + left);
				float d = H[6] * (j + top) - H[3];
				float e = H[7] * (j + top) - H[4];
				float f = H[5] - (j + top);
				float tx = (b * f - c * e) / (b * d - a * e);
				float ty = (c * d - a * f) / (b * d - a * e);
				if(tx >= 0 && tx < nWidth2 && ty >= 0 && ty < nHeight2)
				{
					num++;
					float u1 = tx - int(tx);
					float u2 = ty - int(ty);
					int y = int(ty);
					int x = int(tx);
					if(x < 1 || x >= nWidth2 - 2 || y < 1 || y >= nHeight2 - 2)
					{
						continue;
					}
					float s[8];
					GetS(s ,u1,u2);
					float R[4],G[4],B[4];
					for(int m = 0;m < 4;m++)
					{
						R[m] = (float)(*(pArray2 + pitch2 * y + (x + m - 1) * bitCount2));
						G[m] = (float)(*(pArray2 + pitch2 * y + (x + m - 1) * bitCount2 + 1));
						B[m] = (float)(*(pArray2 + pitch2 * y + (x + m - 1) * bitCount2 + 2));

					}
					float g1[3];
					g1[0] = s[0] * R[0] + s[1] * R[1] + s[2] * R[2] + s[3] * R[3];
					g1[1] = s[0] * G[0] + s[1] * G[1] + s[2] * G[2] + s[3] * G[3];
					g1[2] = s[0] * B[0] + s[1] * B[1] + s[2] * B[2] + s[3] * B[3];

					for(int m = 0;m < 4;m++)
					{
						R[m] = (float)(*(pArray2 + pitch2 * (y - 1) + (x + m - 1) * bitCount2));
						G[m] = (float)(*(pArray2 + pitch2 * (y - 1) + (x + m - 1) * bitCount2 + 1));
						B[m] = (float)(*(pArray2 + pitch2 * (y - 1) + (x + m - 1) * bitCount2 + 2));

					}
					float g2[3];
					g2[0] = s[0] * R[0] + s[1] * R[1] + s[2] * R[2] + s[3] * R[3];
					g2[1] = s[0] * G[0] + s[1] * G[1] + s[2] * G[2] + s[3] * G[3];
					g2[2] = s[0] * B[0] + s[1] * B[1] + s[2] * B[2] + s[3] * B[3];

					for(int m = 0;m < 4;m++)
					{
						R[m] = (float)(*(pArray2 + pitch2 * (y + 1) + (x + m - 1) * bitCount2));
						G[m] = (float)(*(pArray2 + pitch2 * (y + 1) + (x + m - 1) * bitCount2 + 1));
						B[m] = (float)(*(pArray2 + pitch2 * (y + 1) + (x + m - 1) * bitCount2 + 2));

					}
					float g3[3];
					g3[0] = s[0] * R[0] + s[1] * R[1] + s[2] * R[2] + s[3] * R[3];
					g3[1] = s[0] * G[0] + s[1] * G[1] + s[2] * G[2] + s[3] * G[3];
					g3[2] = s[0] * B[0] + s[1] * B[1] + s[2] * B[2] + s[3] * B[3];

					for(int m = 0;m < 4;m++)
					{
						R[m] = (float)(*(pArray2 + pitch2 * (y + 2) + (x + m - 1) * bitCount2));
						G[m] = (float)(*(pArray2 + pitch2 * (y + 2) + (x + m - 1) * bitCount2 + 1));
						B[m] = (float)(*(pArray2 + pitch2 * (y + 2) + (x + m - 1) * bitCount2 + 2));

					}
					float g4[3];
					g4[0] = s[0] * R[0] + s[1] * R[1] + s[2] * R[2] + s[3] * R[3];
					g4[1] = s[0] * G[0] + s[1] * G[1] + s[2] * G[2] + s[3] * G[3];
					g4[2] = s[0] * B[0] + s[1] * B[1] + s[2] * B[2] + s[3] * B[3];

					float t1 = s[4] * g2[0] + s[5] * g1[0] + s[6] * g3[0] + s[7] * g4[0];
					if(t1 > 255)
					{
						t1 = 255.0f;
					}
					if(t1 < 0)
					{
						t1 = 0;
					}
					float t2 = s[4] * g2[1] + s[5] * g1[1] + s[6] * g3[1] + s[7] * g4[1];
					if(t2 > 255)
					{
						t2 = 255.0f;
					}
					if(t2 < 0)
					{
						t2 = 0;
					}
					float t3 = s[4] * g2[2] + s[5] * g1[2] + s[6] * g3[2] + s[7] * g4[2];
					if(t3 > 255)
					{
						t3 = 255.0f;
					}
					if(t3 < 0)
					{
						t3 = 0;
					}
					float y2 = 0.299f*t1 + 0.587f*t2+0.114f*t3;

					float R1 = (float)(*(pArray1 + pitch1 * (j + top) + (i + left) * bitCount1));
					float G1 = (float)(*(pArray1 + pitch1 * (j + top) + (i + left) * bitCount1 + 1));
					float B1 = (float)(*(pArray1 + pitch1 * (j + top) + (i + left) * bitCount1 + 2));
					float y1 = 0.299f*R1 + 0.587f*G1+0.114f*B1;
					dif_R += t1 - R1;
					dif_G += t2 - G1;
					dif_B += t3 - B1;
					dif_Y += y2 - y1;
				}
			}
		}
	}
	dif_R = dif_R / num;
	dif_G = dif_G / num;
	dif_B = dif_B / num;
	//float dif = (dif_R + dif_G + dif_B)/3;
	float dif = dif_Y / num;

	//不在基准图像内，但在待配准图像的投影区域
	for(short j = 0;j < newHeight;j++)
	{
		for(short i = 0;i < newWidth;i++)
		{
			if(i >= nWidth1 - left || i < -1 * left || j >= nHeight1 - top || j < -1 * top)
			{				
				float a = H[6] * (i + left) - H[0];
				float b = H[7] * (i + left) - H[1];
				float c = H[2] - (i + left);
				float d = H[6] * (j + top) - H[3];
				float e = H[7] * (j + top) - H[4];
				float f = H[5] - (j + top);
				float tx = (b * f - c * e) / (b * d - a * e);
				float ty = (c * d - a * f) / (b * d - a * e);
				if(tx >= 0 && tx < nWidth2 && ty >= 0 && ty < nHeight2)
				{
					float u1 = tx - int(tx);
					float u2 = ty - int(ty);
					int y = int(ty);
					int x = int(tx);
					if(x < 1 || x >= nWidth2 - 2 || y < 1 || y >= nHeight2 - 2)
					{
						continue;
					}
					float s[8];
					GetS(s ,u1,u2);
					float R[4],G[4],B[4];
					for(int m = 0;m < 4;m++)
					{
						R[m] = (float)(*(pArray2 + pitch2 * y + (x + m - 1) * bitCount2));
						G[m] = (float)(*(pArray2 + pitch2 * y + (x + m - 1) * bitCount2 + 1));
						B[m] = (float)(*(pArray2 + pitch2 * y + (x + m - 1) * bitCount2 + 2));

					}
					float g1[3];
					g1[0] = s[0] * R[0] + s[1] * R[1] + s[2] * R[2] + s[3] * R[3];
					g1[1] = s[0] * G[0] + s[1] * G[1] + s[2] * G[2] + s[3] * G[3];
					g1[2] = s[0] * B[0] + s[1] * B[1] + s[2] * B[2] + s[3] * B[3];

					for(int m = 0;m < 4;m++)
					{
						R[m] = (float)(*(pArray2 + pitch2 * (y - 1) + (x + m - 1) * bitCount2));
						G[m] = (float)(*(pArray2 + pitch2 * (y - 1) + (x + m - 1) * bitCount2 + 1));
						B[m] = (float)(*(pArray2 + pitch2 * (y - 1) + (x + m - 1) * bitCount2 + 2));

					}
					float g2[3];
					g2[0] = s[0] * R[0] + s[1] * R[1] + s[2] * R[2] + s[3] * R[3];
					g2[1] = s[0] * G[0] + s[1] * G[1] + s[2] * G[2] + s[3] * G[3];
					g2[2] = s[0] * B[0] + s[1] * B[1] + s[2] * B[2] + s[3] * B[3];

					for(int m = 0;m < 4;m++)
					{
						R[m] = (float)(*(pArray2 + pitch2 * (y + 1) + (x + m - 1) * bitCount2));
						G[m] = (float)(*(pArray2 + pitch2 * (y + 1) + (x + m - 1) * bitCount2 + 1));
						B[m] = (float)(*(pArray2 + pitch2 * (y + 1) + (x + m - 1) * bitCount2 + 2));

					}
					float g3[3];
					g3[0] = s[0] * R[0] + s[1] * R[1] + s[2] * R[2] + s[3] * R[3];
					g3[1] = s[0] * G[0] + s[1] * G[1] + s[2] * G[2] + s[3] * G[3];
					g3[2] = s[0] * B[0] + s[1] * B[1] + s[2] * B[2] + s[3] * B[3];

					for(int m = 0;m < 4;m++)
					{
						R[m] = (float)(*(pArray2 + pitch2 * (y + 2) + (x + m - 1) * bitCount2));
						G[m] = (float)(*(pArray2 + pitch2 * (y + 2) + (x + m - 1) * bitCount2 + 1));
						B[m] = (float)(*(pArray2 + pitch2 * (y + 2) + (x + m - 1) * bitCount2 + 2));

					}
					float g4[3];
					g4[0] = s[0] * R[0] + s[1] * R[1] + s[2] * R[2] + s[3] * R[3];
					g4[1] = s[0] * G[0] + s[1] * G[1] + s[2] * G[2] + s[3] * G[3];
					g4[2] = s[0] * B[0] + s[1] * B[1] + s[2] * B[2] + s[3] * B[3];

					float t1 = s[4] * g2[0] + s[5] * g1[0] + s[6] * g3[0] + s[7] * g4[0];
					if(t1 > 255)
					{
						t1 = 255.0f;
					}
					if(t1 < 0)
					{
						t1 = 0;
					}
					float t2 = s[4] * g2[1] + s[5] * g1[1] + s[6] * g3[1] + s[7] * g4[1];
					if(t2 > 255)
					{
						t2 = 255.0f;
					}
					if(t2 < 0)
					{
						t2 = 0;
					}
					float t3 = s[4] * g2[2] + s[5] * g1[2] + s[6] * g3[2] + s[7] * g4[2];
					if(t3 > 255)
					{
						t3 = 255.0f;
					}
					if(t3 < 0)
					{
						t3 = 0;
					}
					*(nArray + npitch * j + i * bitCount)     = (unsigned char)(t1);
					*(nArray + npitch * j + i * bitCount + 1) = (unsigned char)(t2);
					*(nArray + npitch * j + i * bitCount + 2) = (unsigned char)(t3);
				}
			}
		}
	}
	//在基准图像内的点

	for(short j = 0;j < newHeight;j++)
	{
		for(short i = 0;i < newWidth;i++)
		{
			if(i < nWidth1 - left && i >= -1 * left && j < nHeight1 - top && j >= -1 * top)
			{	
				float R = (float)(*(pArray1 + pitch1 * (j + top) + (i + left) * bitCount1)) /*+ dif*/;
				float G = (float)(*(pArray1 + pitch1 * (j + top) + (i + left) * bitCount1 + 1))/* + dif*/;
				float B = (float)(*(pArray1 + pitch1 * (j + top) + (i + left) * bitCount1 + 2)) /*+ dif*/;
				float y = 0.299f*R + 0.587f*R + 0.114f*B;
				float u = -0.147f*R - 0.289f*G + 0.436f*B;
				float v = 0.615f*R - 0.515f*G - 0.1f*B;
				y = y + dif;
				R = y+1.14f*v;
				G = y - 0.39f*u - 0.58f*v;
				B = y + 2.03f*u;
				R = R > 255 ? 255 : (R <0 ? 0 : R);
				G = G > 255 ? 255 : (G <0 ? 0 : G);
				B = B > 255 ? 255 : (B <0 ? 0 : B);
				*(nArray + npitch * j + i * bitCount)     = (unsigned char)(R);
				*(nArray + npitch * j + i * bitCount + 1) = (unsigned char)(G);
				*(nArray + npitch * j + i * bitCount + 2) = (unsigned char)(B);
			}
		}
	}




	//在待配准图像的投影区域,也在基准图像边缘的过渡区域

	for(short j = 0;j < newHeight;j++)
	{
		for(short i = 0;i < newWidth;i++)
		{
			if(i < nWidth1 - left && i >= 0 - left && j < nHeight1 - top && j >= 0 - top)
			{
				if(i <= nWidth1 - dis - left && i >= dis - left && j <= nHeight1 - dis - top && j >= dis - top)
				{}
				else 
				{		
					float a = H[6] * (i + left) - H[0];
					float b = H[7] * (i + left) - H[1];
					float c = H[2] - (i + left);
					float d = H[6] * (j + top) - H[3];
					float e = H[7] * (j + top) - H[4];
					float f = H[5] - (j + top);
					float tx = (b * f - c * e) / (b * d - a * e);
					float ty = (c * d - a * f) / (b * d - a * e);
					if(tx >= 0 && tx < nWidth2 && ty >= 0 && ty < nHeight2)
					{
						float u1 = tx - int(tx);
						float u2 = ty - int(ty);
						short y = int(ty);
						short x = int(tx);
						if(x < 1 || x >= nWidth2 - 2 || y < 1 || y >= nHeight2 - 2)
						{
							continue;
						}
						float R = (float)(*(pArray2 + pitch2 * y + (x) * bitCount2));
						float G = (float)(*(pArray2 + pitch2 * y + (x) * bitCount2 + 1));
						float B = (float)(*(pArray2 + pitch2 * y + (x) * bitCount2 + 2));
						float t1 = u2 * (u1 * R + (1 - u1) * R)
							+ (1 - u2) * (u1 * R +(1 - u1) * R);
						float t2 = u2 * (u1 * G + (1 - u1) * G)
							+ (1 - u2) * (u1 * G +(1 - u1) * G);
						float t3 = u2 * (u1 * B + (1 - u1) * B)
							+ (1 - u2) * (u1 * B +(1 - u1) * B);
						float w = 0;

						float d1 = float(nWidth1 - 1 - (i + left));
						float d2 = float(i + left - 0);
						float d3 = float(nHeight1 - 1 - (j + top));
						float d4 = float(j + top);
						float min = 30;
						if(d1 < min)
						{
							min = d1;
						}
						if(d2 < min)
						{
							min = d2;
						}
						if(d3 < min)
						{
							min = d3;
						}
						if(d4 < min)
						{
							min = d4;
						}
						w = (30 - min) / 30;

						R = (float)(*(pArray1 + pitch1 * (j + top) + (i + left) * bitCount1))/* + dif*/;
						G = (float)(*(pArray1 + pitch1 * (j + top) + (i + left) * bitCount1 + 1))/* + dif*/;
						B = (float)(*(pArray1 + pitch1 * (j + top) + (i + left) * bitCount1 + 2)) /*+ dif*/;
						float Y = 0.299f*R + 0.587f*R + 0.114f*B;
						float U = -0.147f*R - 0.289f*G + 0.436f*B;
						float V = 0.615f*R - 0.515f*G - 0.1f*B;
						Y = Y + dif;
						R = Y+1.14f*V;
						G = Y - 0.39f*U - 0.58f*V;
						B = Y + 2.03f*U;
						float H = w * t1 + (1 - w) * R;
						if(H > 255.0f)
						{
							H = 255.0f;
						}
						if(H < 0)
						{
							H = 0;
						}
						float S = w * t2 + (1 - w) * G;
						if(S > 255.0f)
						{
							S = 255.0f;
						}
						if(S < 0)
						{
							S = 0;
						}
						float I = w * t3 + (1 - w) * B;
						if(I > 255.0f)
						{
							I = 255.0f;
						}
						if(I < 0)
						{
							I = 0;
						}
						//HSItoRGB(H,S,I,&R,&G,&B);
						*(nArray + npitch * j + i * bitCount)     = (unsigned char)(H);
						*(nArray + npitch * j + i * bitCount + 1) = (unsigned char)(S);
						*(nArray + npitch * j + i * bitCount + 2) = (unsigned char)(I);
					}
				}
			}
		}
	}


	//不在基准图像区域内且不在待配准图像的投影区域
	for(short j = 0;j < newHeight;j++)
	{
		for(short i = 0;i < newWidth;i++)
		{
			if(i >= nWidth1 - left || i < -1 * left || j >= nHeight1 - top || j < -1 * top)
			{
				float a = H[6] * (i + left) - H[0];
				float b = H[7] * (i + left) - H[1];
				float c = H[2] - (i + left);
				float d = H[6] * (j + top) - H[3];
				float e = H[7] * (j + top) - H[4];
				float f = H[5] - (j + top);
				float tx = (b * f - c * e) / (b * d - a * e);
				float ty = (c * d - a * f) / (b * d - a * e);
				int x = int(tx + 0.5f);
				int y = int(ty + 0.5f);
				if(tx <= 1 || tx >= nWidth2 - 2 || ty <= 1 || ty >= nHeight2 - 2)
				{				
					*(nArray + npitch * j + i * bitCount)     = 0;
					*(nArray + npitch * j + i * bitCount + 1) = 0;
					*(nArray + npitch * j + i * bitCount + 2) = 0;
				}
			}
		}
	}
}

// 计算图象经过投影变换后的尺寸大小
void CHarrisCorners::getAftProDim(short *left,short *right,short *top,short *bottom,short *newHeight,short *newWidth,short height1, short width1, short height2, short width2)
{
	//算出新图像left,right ,top,bottom在投影变换坐标系内的坐标值
	short left1 = 0;
	short right1 = width1 - 1;
	short top1 = 0;
	short bottom1 = height1 - 1;
	//待配准图像左上角点坐标(m1,n1)，并计算在投影变换后的坐标(_m1,_n1)
	float m1 = 0.0f;
	float n1 = 0.0f;
	float _m1 = H[0] * m1 + H[1] * n1 + H[2];
	float _n1 = H[3] * m1 + H[4] * n1 + H[5];
	float  w1 = H[6] * m1 + H[7] * n1 + 1;
	_m1 = _m1 / w1;
	_n1 = _n1 / w1;
	if(_m1 < left1)
	{
		left1 = (int)_m1;
	}
	if(_m1 > right1)
	{
		right1 = (int)_m1 + 1;
	}
	if(_n1 < top1)
	{
		top1 = (int)_n1;
	}
	if(_n1 > bottom1)
	{
		bottom1 = (int)_n1 + 1;
	}
	//待配准图像右上角点坐标(m2,n2)，并计算在投影变换后的坐标(_m2,_n2)
	float m2 = float((width2 - 1));
	float n2 = 0.0f;
	float _m2 = H[0] * m2 + H[1] * n2 + H[2];
	float _n2 = H[3] * m2 + H[4] * n2 + H[5];
	float  w2 = H[6] * m2 + H[7] * n2 + 1;
	_m2 = _m2 / w2;
	_n2 = _n2 / w2;
	if(_m2 < left1)
	{
		left1 = (int)_m2;
	}
	if(_m2 > right1)
	{
		right1 = (int)_m2 + 1;
	}
	if(_n2 < top1)
	{
		top1 = (int)_n2;
	}
	if(_n2 > bottom1)
	{
		bottom1 = (int)_n2 + 1;
	}
	//待配准图像左下角点坐标(m3,n3)，并计算在投影变换后的坐标(_m3,_n3)
	float m3 = 0.0f;
	float n3 = float(height2 - 1);
	float _m3 = H[0] * m3 + H[1] * n3 + H[2];
	float _n3 = H[3] * m3 + H[4] * n3 + H[5];
	float  w3 = H[6] * m3 + H[7] * n3 + 1;
	_m3 = _m3 / w3;
	_n3 = _n3 / w3;
	if(_m3 < left1)
	{
		left1 = (int)_m3;
	}
	if(_m3 > right1)
	{
		right1 = (int)_m3 + 1;
	}
	if(_n3 < top1)
	{
		top1 = (int)_n3;
	}
	if(_n3 > bottom1)
	{
		bottom1 = (int)_n3 + 1;
	}
	//待配准图像右下角点坐标(m4,n4)，并计算在投影变换后的坐标(_m4,_n4)
	float m4 = float((width2 - 1));
	float n4 = float(height2 - 1);
	float _m4 = H[0] * m4 + H[1] * n4 + H[2];
	float _n4 = H[3] * m4 + H[4] * n4 + H[5];
	float  w4 = H[6] * m4 + H[7] * n4 + 1;
	_m4 = _m4 / w4;
	_n4 = _n4 / w4;
	if(_m4 < left1)
	{
		left1 = (int)_m4;
	}
	if(_m4 > right1)
	{
		right1 = (int)_m4 + 1;
	}
	if(_n4 < top1)
	{
		top1 = (int)_n4;
	}
	if(_n4 > bottom1)
	{
		bottom1 = (int)_n4 + 1;
	}
	*newHeight = bottom1 - top1 + 1;
	*newWidth = right1- left1 + 1;
	*left = left1;
	*right = right1;
	*top = top1;
	*bottom = bottom1;
}