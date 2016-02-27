#pragma once
#include "atlimage.h"
//class CImageMosaicDlg;
class CDWT
{
public:
	CDWT(void);
	virtual ~CDWT(void);
	bool m_bTwice;
	bool m_bOnce;
	void WvltTransTwice(CImage &pImage,float*& I);
	void DWT_TwoLayers(short**& spOriginData, short**& spTransData0, short**& spTransData1,
						int nHeight, int nHeight_H, int nWidth, int nWidth_H, int layer, float fRadius);
	void WvltTransOnce(CImage &pImage,float*& I);
	void DWT_Once(short**& spOriginData, short**& spTransData0, short**& spTransData1, 
				   int nHeight, int nHeight_H, int nWidth, int nWidth_H, int layer, float fRadius);
};
