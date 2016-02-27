// ImageMosaicDlg.h : header file
//

#pragma once
#include "afxwin.h"
#include "atlimage.h"
#include "DWT.h"
#include "HarrisCorners.h"
// CImageMosaicDlg dialog

class CImageMosaicDlg : public CDialog
{
// Construction
public:
	CImageMosaicDlg(CWnd* pParent = NULL);	// standard constructor
	 ~CImageMosaicDlg();

// Dialog Data
	enum { IDD = IDD_IMAGEMOSAIC_DIALOG };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV support


// Implementation
protected:
	HICON m_hIcon;

	// Generated message map functions
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void Loadimage();
	afx_msg void StartImageMosaic();
	bool Imagemosaic(CImage& baseImage,CImage& stitchingImage,CImage& resultImage);
	void OnClose(void);

public:
	//list<ImageGroup> srcImageGroup;
	short m_nHeight1;
	short m_nWidth1;
	short m_nHeight2;
	short m_nWidth2;
	short m_nImageNum;

	list<CString> pathNameList;
private:
	CHarrisCorners harris; 
public:
	bool wvltTransTwiceProcess(CImage &baseImage,CImage &stitchingImage);
	bool wvltTransOnceProcess(CImage &baseImage,CImage &stitchingImage);
	bool originImageProcess(CImage &baseImage,CImage &stitchingImage,CImage &resultImage);
	afx_msg void Cancel();
	void imageCopy(CImage &destImg,CImage &srcImg);
	bool imageGroupMosaic(list<CString> &_list,CImage& resultImage);
	void getYUV(CImage &image,float*& P);
};
