#include "StdAfx.h"
#include "DWT.h"

CDWT::CDWT(void)
: m_bTwice(false)
, m_bOnce(false)
{
}

CDWT::~CDWT(void)
{
}

void CDWT::WvltTransTwice(CImage &pImage,float* &I)
{
	//获取图像的属性参数
	BYTE* pArray = (BYTE*)pImage.GetBits(); //获得数据
	int pitch = pImage.GetPitch();
	int bitCount = pImage.GetBPP() / 8;  //每像素位数
	int height = pImage.GetHeight();
	int width = pImage.GetWidth();

	//图像矩阵坐标与像素数据	
	int x,y;
	unsigned char tempR, tempG, tempB;
	float fTempBufforDisp;
	short MaxPixVal,MinPixVal,Diff;
	short **spOriginData, **spTransData0, **spTransData1;
	//分配图像小波变换所用的数据空间
	spOriginData = new short* [height];
	spTransData0 = new short* [height];
	spTransData1 = new short* [height];

	for(int i = 0; i < height; i ++)
	{
		spOriginData[i] = new short [width];
		spTransData0[i] = new short [width];
		spTransData1[i] = new short [width];
	}

	//从设备缓存中获取原始图像数据
	for(y = 0; y < height; y++)
	{
		for( x = 0; x < width; x++)
		{
			tempR = char((*(pArray + pitch * y + bitCount * x)));
			tempG = char((*(pArray + pitch * y + bitCount * x + 1)));
			tempB = char((*(pArray + pitch * y + bitCount * x + 2)));
			spOriginData[y][x] = short(0.299f * tempR + 0.587f * tempG + 0.114f * tempB);
		}
	}
	//完成图像的两次小波变换
	DWT_TwoLayers(spOriginData,spTransData0,spTransData1,height,height/2,width,width/2,2,1.414f);
	MaxPixVal=spTransData1[0][0];
	MinPixVal=spTransData1[0][0];
	//计算得到图像小波系数的极大值与极小值
	for( y=0; y<height/4; y++)
	{
		for( x=0; x<width/4; x++)
		{
			if(MaxPixVal<spTransData1[y][x])
			{
				MaxPixVal=spTransData1[y][x];
			}
			if(MinPixVal>spTransData1[y][x])
			{
				MinPixVal=spTransData1[y][x];
			}
		}
	}
	//计算获得小波系数的极值差
	Diff=MaxPixVal-MinPixVal;
	//小波系数经过处理后，放入显示缓存中
	for(y=0; y<height / 4; y++)
	{
		for(x=0; x<width / 4; x++)
		{
			//因为小波变换后的小波系数有可能超过255甚至更多，那么就将
			//小波系数的范围映射到0~255区间内，以后出现类似的处理，目的都是一样的
			fTempBufforDisp=spTransData1[y][x];
			fTempBufforDisp-=MinPixVal;
			fTempBufforDisp*=255;
			fTempBufforDisp/=Diff;
			I[y*(width/4) + x] = fTempBufforDisp;
		}
	}

	//删除临时的数据空间
	for(int a = 0; a < height; a ++)
	{
		delete[]spOriginData[a];
		spOriginData[a] = NULL;
		delete[]spTransData0[a];
		spTransData0[a] = NULL;
		delete[]spTransData1[a];
		spTransData1[a] = NULL;
	}
	delete[]spOriginData;
	spOriginData = NULL;
	delete[]spTransData0;
	spTransData0 = NULL;
	delete[]spTransData1;
	spTransData1 = NULL;
}
void CDWT::DWT_TwoLayers(short**& spOriginData, short**& spTransData0, short**& spTransData1, int nHeight, int nHeight_H, int nWidth, int nWidth_H, int layer, float fRadius)
{
	//利用循环完成两次小波变换
	for(int i=1; i<=layer; i++)
	{
		DWT_Once(spOriginData,spTransData0,spTransData1,nHeight,nHeight_H,nWidth,nWidth_H,i,fRadius);
		nHeight   = nHeight/2;		
		nWidth    = nWidth/2;
		nHeight_H = nHeight/2;	
		nWidth_H  = nWidth/2;
	}

}
void CDWT::WvltTransOnce(CImage &pImage,float* &I)
{
	//获取图像的属性参数 
	BYTE* pArray = (BYTE*)pImage.GetBits(); //获得数据
	int pitch = pImage.GetPitch();
	int bitCount = pImage.GetBPP() / 8;  //每像素位数
	int height = pImage.GetHeight();
	int width = pImage.GetWidth();

	//图像矩阵坐标与像素数据	
	int x,y;
	unsigned char tempR, tempG, tempB;
	float fTempBufforDisp;
	short MaxPixVal,MinPixVal,Diff;
	short **spOriginData, **spTransData0, **spTransData1;
	//分配图像小波变换所用的数据空间
	spOriginData = new short* [height];
	spTransData0 = new short* [height];
	spTransData1 = new short* [height];

	for(int i = 0; i < height; i ++)
	{
		spOriginData[i] = new short [width];
		spTransData0[i] = new short [width];
		spTransData1[i] = new short [width];
	}

	//从设备缓存中获取原始图像数据
	for(y = 0; y < height; y++)
	{
		for( x = 0; x < width; x++)
		{
			tempR = char((*(pArray + pitch * y + bitCount * x)));
			tempG = char((*(pArray + pitch * y + bitCount * x + 1)));
			tempB = char((*(pArray + pitch * y + bitCount * x + 2)));
			spOriginData[y][x] = short(0.299f * tempR + 0.587f * tempG + 0.114f * tempB);
		}
	}
	//完成一次图像小波变换
	DWT_Once(spOriginData,spTransData0,spTransData1,height,height/2,width,width/2,1,1.414f);
	MaxPixVal = spTransData1[0][0];
	MinPixVal = spTransData1[0][0];
	for( y = 0; y < height/2; y++)
	{
		for( x = 0; x < width/2; x++)
		{
			if(MaxPixVal < spTransData1[y][x])
			{
				MaxPixVal = spTransData1[y][x];
			}
			if(MinPixVal > spTransData1[y][x])
			{
				MinPixVal = spTransData1[y][x];
			}
		}
	}
	Diff = MaxPixVal - MinPixVal;
	for(y = 0; y < height / 2; y++)
	{
		for(x = 0; x < width / 2; x++)
		{
			//因为小波变换后的小波系数有可能超过255甚至更多，那么就将
			//小波系数的范围映射到0~255区间内，以后出现类似的处理，目的都是一样的
			fTempBufforDisp = spTransData1[y][x];
			fTempBufforDisp -= MinPixVal;
			fTempBufforDisp *= 255;
			fTempBufforDisp /= Diff;
			I[y*(width/2) + x] = fTempBufforDisp;
		}
	}
	//删除临时的数据空间
	for(int a = 0; a < height; a ++)
	{
		delete[]spOriginData[a];
		spOriginData[a] = NULL;
		delete[]spTransData0[a];
		spTransData0[a] = NULL;
		delete[]spTransData1[a];
		spTransData1[a] = NULL;
	}
	delete[]spOriginData;
	spOriginData = NULL;
	delete[]spTransData0;
	spTransData0 = NULL;
	delete[]spTransData1;
	spTransData1 = NULL;
}
void CDWT::DWT_Once(short**& spOriginData, short**& spTransData0, short**& spTransData1, int nHeight, int nHeight_H, int nWidth, int nWidth_H, int layer, float fRadius)
{
	int Trans_W,				//图像扫描线控制：横坐标
		Trans_H,				//图像扫描线控制：纵坐标
		Trans_M,				//图像矩阵的横坐标
		Trans_N;				//图像矩阵的纵坐标
	short Trans_Coeff0;			//小波变换系数
    signed short Trans_Coeff1;
	fRadius = 1.414f;				//变换滤波系数
	//本模块完成变换系数的赋值采样
	//行变换,第一次（layer=1时）时nHeight即为原始图像的高度值
    for(Trans_H = 0; Trans_H < nHeight; Trans_H++)            
	{
		if(layer == 1)
		{
			 //layer=1时，nWidth_H为原始图像宽度值的一半
			for(Trans_N = 0; Trans_N < nWidth_H; Trans_N++)          
			{
				Trans_W = Trans_N<<1;
	            if (fRadius == 2)
				{
					spTransData0[Trans_H][Trans_N] = (spOriginData[Trans_H][Trans_W]);
                    spTransData0[Trans_H][nWidth_H + Trans_N] = (spOriginData[Trans_H][Trans_W + 1]);
				}
	            else
				{
                    spTransData0[Trans_H][Trans_N] = (spOriginData[Trans_H][Trans_W]- 128);		
                    spTransData0[Trans_H][nWidth_H + Trans_N] = (spOriginData[Trans_H][Trans_W + 1] - 128);	
				}
	   		}
		}
		//若变换层数大于1,则仅采样低频的小波系数
		if(layer > 1)
		{
			for(Trans_N = 0; Trans_N < nWidth_H; Trans_N++)
			{
				Trans_W = Trans_N<<1;
				spTransData0[Trans_H][Trans_N] = spTransData1[Trans_H][Trans_W];
				spTransData0[Trans_H][nWidth_H + Trans_N] = spTransData1[Trans_H][Trans_W + 1];
			}
		}
	}
	for(Trans_H = 0; Trans_H < nHeight; Trans_H++)
	{
		for(Trans_N = 0; Trans_N < nWidth_H - 1; Trans_N++)
		{
			//奇偶数值和的一半
			Trans_Coeff1 = ((spTransData0[Trans_H][Trans_N] + spTransData0[Trans_H][Trans_N + 1])>>1);	
			//逻辑非操作后数值加1
			Trans_Coeff1 = ~Trans_Coeff1 + 1;	
			//系数预测
			spTransData0[Trans_H][nWidth_H + Trans_N] = spTransData0[Trans_H][nWidth_H + Trans_N] + Trans_Coeff1;	
		}
		//完成一个偶系数的边界处理
		Trans_Coeff1 = ((spTransData0[Trans_H][nWidth_H - 1] + spTransData0[Trans_H][nWidth_H - 2])>>1);
		Trans_Coeff1 = ~Trans_Coeff1 + 1;
		spTransData0[Trans_H][nWidth - 1] = spTransData0[Trans_H][nWidth - 1] + Trans_Coeff1;
		//完成一个奇系数的边界处理
		Trans_Coeff0 = ((spTransData0[Trans_H][nWidth_H] + spTransData0[Trans_H][nWidth_H + 1])>>2);
		spTransData0[Trans_H][0] = spTransData0[Trans_H][0] + Trans_Coeff0;
		//提升，整数到整数的变换
		for(Trans_N = 1; Trans_N < nWidth_H; Trans_N++)
		{
			Trans_Coeff0 = ((spTransData0[Trans_H][nWidth_H + Trans_N] + spTransData0[Trans_H][nWidth_H + Trans_N - 1])>>2);
			spTransData0[Trans_H][Trans_N] = spTransData0[Trans_H][Trans_N] + Trans_Coeff0;
		}

	}//水平方向的变换结束
	//竖直方向的变换开始，数据源未水平变换后的小波系数
	for(Trans_M = 0; Trans_M < nHeight; Trans_M++)
	{
		for(Trans_N = 0; Trans_N < nWidth_H; Trans_N++)
		{
			spTransData0[Trans_M][Trans_N] = short(spTransData0[Trans_M][Trans_N] * fRadius);
			spTransData0[Trans_M][Trans_N + nWidth_H] = short(spTransData0[Trans_M][Trans_N + nWidth_H] / fRadius);
		}
	}
	//行提升后的数据在spTransData0中，spTransData0中的数据自然奇偶有序
	for(Trans_N = 0; Trans_N < nWidth_H; Trans_N++)
	{
		//列变换
		for(Trans_M = 0; Trans_M < nHeight_H; Trans_M++)
		{
			Trans_H = Trans_M<<1;
			//频带LL部分
			spTransData1[Trans_M][Trans_N] = spTransData0[Trans_H][Trans_N];
			//频带HL部分
			spTransData1[nHeight_H + Trans_M][Trans_N] = spTransData0[Trans_H + 1][Trans_N];
			//频带LH部分
			spTransData1[Trans_M][nWidth_H + Trans_N] = spTransData0[Trans_H][nWidth_H + Trans_N];	
			//频带HH部分
			spTransData1[nHeight_H + Trans_M][nWidth_H + Trans_N] = spTransData0[Trans_H + 1][nWidth_H + Trans_N];
		}
		//第一次提升奇数坐标系数
		for(Trans_M = 0; Trans_M < nHeight_H - 1; Trans_M++)
		{
			//竖直方向的变换 
			Trans_Coeff1 = ((spTransData1[Trans_M][Trans_N]+spTransData1[Trans_M+1][Trans_N])>>1);
			Trans_Coeff1=~Trans_Coeff1+1;
			spTransData1[nHeight_H+Trans_M][Trans_N] = spTransData1[nHeight_H+Trans_M][Trans_N]+Trans_Coeff1;
			Trans_Coeff1 = ((spTransData1[Trans_M][nWidth_H+Trans_N]+spTransData1[Trans_M+1][nWidth_H+Trans_N])>>1);
			Trans_Coeff1=~Trans_Coeff1+1;
			spTransData1[nHeight_H+Trans_M][nWidth_H+Trans_N] = spTransData1[nHeight_H+Trans_M][nWidth_H+Trans_N]+Trans_Coeff1;
		}
		Trans_Coeff1 = ((spTransData1[nHeight_H-1][Trans_N]+spTransData1[nHeight_H-2][Trans_N])>>1);
		Trans_Coeff1=~Trans_Coeff1+1;
		spTransData1[nHeight-1][Trans_N] = spTransData1[nHeight-1][Trans_N]+Trans_Coeff1;
		Trans_Coeff1 = ((spTransData1[nHeight_H-1][nWidth_H+Trans_N]+spTransData1[nHeight_H-2][nWidth_H+Trans_N])>>1);
		Trans_Coeff1=~Trans_Coeff1+1;
		//边界处理
		spTransData1[nHeight-1][nWidth_H+Trans_N] = spTransData1[nHeight-1][nWidth_H+Trans_N]+Trans_Coeff1;

		Trans_Coeff0 = ((spTransData1[nHeight_H][Trans_N]+spTransData1[nHeight_H+1][Trans_N])>>2);
		spTransData1[0][Trans_N] = spTransData1[0][Trans_N]+Trans_Coeff0;
		Trans_Coeff0 = ((spTransData1[nHeight_H][nWidth_H+Trans_N]+spTransData1[nHeight_H+1][nWidth_H+Trans_N])>>2);
		//边界处理
		spTransData1[0][nWidth_H+Trans_N] = spTransData1[0][nWidth_H+Trans_N]+Trans_Coeff0;
		//第一次提升偶数坐标系数
		for(Trans_M = 1; Trans_M < nHeight_H; Trans_M++)
		{
			Trans_Coeff0 = ((spTransData1[nHeight_H+Trans_M][Trans_N]+spTransData1[nHeight_H+Trans_M-1][Trans_N])>>2);
			spTransData1[Trans_M][Trans_N] = spTransData1[Trans_M][Trans_N]+Trans_Coeff0;
			Trans_Coeff0 = ((spTransData1[nHeight_H+Trans_M][nWidth_H+Trans_N]+spTransData1[nHeight_H+Trans_M-1][nWidth_H+Trans_N])>>2);
			spTransData1[Trans_M][nWidth_H+Trans_N] = spTransData1[Trans_M][nWidth_H+Trans_N]+Trans_Coeff0;
		}
	}
	//存放小波系数，LL频带的系数进行幅值增强处理，其它高频频带的系数则削弱其幅值
	for(Trans_N = 0; Trans_N < nWidth; Trans_N++)
	{
		for(Trans_M = 0; Trans_M < nHeight_H; Trans_M++)
		{
			spTransData1[Trans_M][Trans_N] = short(spTransData1[Trans_M][Trans_N] * fRadius);
			spTransData1[Trans_M+nHeight_H][Trans_N] = short(spTransData1[Trans_M+nHeight_H][Trans_N] / fRadius);
		}
	}
}
