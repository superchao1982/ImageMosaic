#pragma once
#include "afxwin.h"
#include "atlimage.h" 
// AftRegDlg dialog

class AftRegDlg : public CDialog
{
	DECLARE_DYNAMIC(AftRegDlg)
public:
	//CImage* p;
	AftRegDlg();
	AftRegDlg(CWnd* pParent,CImage& image);   // standard constructor
	virtual ~AftRegDlg();	

// Dialog Data
	enum { IDD = IDD_DIALOG_NewImage };

protected:
	HICON m_hIcon;
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support
	virtual BOOL OnInitDialog();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void Imagesave();
	void OnClose(void);
	void PostNcDestroy(void);
	afx_msg void Cancel();
	CImage m_AftRegImage;
	void imageCopy(CImage &destImg,CImage &srcImg);
	//CImage* m_AftRegImage_1;
protected:
	void OnPaint(void);
public:
	CStatic imageShow;
	void CalImageLocation2(void);
	bool m_bCalImgLoc2;
};
