// ImageMosaicDlg.cpp : implementation file
//

#include "stdafx.h"
#include "ImageMosaic.h"
#include "ImageMosaicDlg.h"
#include "DWT.h"
#include "AftRegDlg.h"
#include "direct.h"
#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// CAboutDlg dialog used for App About

class CAboutDlg : public CDialog
{
public:
	CAboutDlg();

// Dialog Data
	enum { IDD = IDD_ABOUTBOX };

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

// Implementation
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialog(CAboutDlg::IDD)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialog)
END_MESSAGE_MAP()


// CImageMosaicDlg dialog




CImageMosaicDlg::CImageMosaicDlg(CWnd* pParent /*=NULL*/)
	: CDialog(CImageMosaicDlg::IDD, pParent)
	, m_nHeight1(0)
	, m_nWidth1(0)
	, m_nHeight2(0)
	, m_nWidth2(0)
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}
CImageMosaicDlg::~CImageMosaicDlg()
{
	
}
void CImageMosaicDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialog::DoDataExchange(pDX);
	//DDX_Control(pDX, IDC_ImageShow1, m_PictureShow1);
	//DDX_Control(pDX, IDC_ImageShow2, m_PictureShow2);
}

BEGIN_MESSAGE_MAP(CImageMosaicDlg, CDialog)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	//}}AFX_MSG_MAP
	ON_BN_CLICKED(IDC_LoadImage1, &CImageMosaicDlg::Loadimage)
	ON_BN_CLICKED(IDC_ImageMosaic, &CImageMosaicDlg::StartImageMosaic)
	ON_BN_CLICKED(IDCANCEL, &CImageMosaicDlg::Cancel)
END_MESSAGE_MAP()


// CImageMosaicDlg message handlers

BOOL CImageMosaicDlg::OnInitDialog()
{
	CDialog::OnInitDialog();

	// Add "About..." menu item to system menu.

	// IDM_ABOUTBOX must be in the system command range.
	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != NULL)
	{
		CString strAboutMenu;
		strAboutMenu.LoadString(IDS_ABOUTBOX);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// Set the icon for this dialog.  The framework does this automatically
	//  when the application's main window is not a dialog
	SetIcon(m_hIcon, TRUE);			// Set big icon
	SetIcon(m_hIcon, FALSE);		// Set small icon
	//CalImageLocation();
	ShowWindow(SW_NORMAL);
	AfxGetMainWnd()->CenterWindow();
	// TODO: Add extra initialization here

	return TRUE;  // return TRUE  unless you set the focus to a control
}

void CImageMosaicDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal();
	}
	else
	{
		CDialog::OnSysCommand(nID, lParam);
	}
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CImageMosaicDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // device context for painting

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// Center icon in client rectangle
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// Draw the icon
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialog::OnPaint();
	}
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CImageMosaicDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}
/********************************************************************************************/
/*******************************************加载图像*****************************************/
/********************************************************************************************/
void CImageMosaicDlg::Loadimage()
{
	// TODO: Add your control notification handler code here
	//获取CImage支持的图像文件的过滤字符串
	CFileDialog   dlgFile(TRUE,NULL,NULL,OFN_FILEMUSTEXIST|OFN_ALLOWMULTISELECT, _T("All   Files   (*.*)|*.*||")); 
	char   path[_MAX_DIR]; 
	_getcwd(path,   _MAX_DIR); 
	dlgFile.m_ofn.lpstrInitialDir   = (LPCWSTR)path; 

	if(dlgFile.DoModal()   !=   IDOK)  
	{
		return; 
	}

	//operate   the   files   selected   one   by   one   
	POSITION   pos= dlgFile.GetStartPosition(); 
	int a = -1;
	while   (pos   !=   NULL) 
	{ 
		a++;
		CString filename = dlgFile.GetNextPathName(pos);
		pathNameList.push_back(filename);
	} 
	if(a < 1)
	{
		AfxMessageBox(_T("请继续加载图像！"));
		return;
	}
	m_nImageNum = a + 1;
}
void CImageMosaicDlg::imageCopy(CImage &destImg,CImage &srcImg)		//图的复制 
{
	int width = srcImg.GetWidth();
	int pitch = srcImg.GetPitch();
	int height = srcImg.GetHeight();
	int bytesPerPixel = srcImg.GetBPP() / 8;
	destImg.Create(width, height, bytesPerPixel * 8, 0);
	int _pitch = destImg.GetPitch();

	BYTE *pDestD = (BYTE *)destImg.GetBits();
	BYTE *pSrcD = (BYTE *)srcImg.GetBits();

	for (int y=0; y<height; y++)
	{
		memcpy(pDestD + y * _pitch, pSrcD + y * pitch, abs(pitch));
	}

	if (srcImg.GetBPP() <= 8)
	{
		RGBQUAD  D_pal[256];
		srcImg.GetColorTable(0, 256, D_pal);
		destImg.SetColorTable(0, 256, D_pal);
	}
}
void CImageMosaicDlg::StartImageMosaic()
{
	float t = (float)GetTickCount();
	if(m_nImageNum == 2)//拼接图像幅数为2
	{
		list<CString>::iterator iter1 = pathNameList.begin();
		CImage stitchingImage;
		stitchingImage.Load((*iter1));

		iter1 = pathNameList.erase(iter1);
		CImage baseImage;
		baseImage.Load((*iter1));
		iter1 = pathNameList.erase(iter1);
		CImage resultImage;
		if(Imagemosaic(baseImage,stitchingImage,resultImage))
		{
			pathNameList.clear();
			return;				
		}
		pathNameList.clear();
		t = ((float)GetTickCount() - t) / 1000;
		CString str;
		str.Format(_T("%f"),t);
		AfxMessageBox(str);
		 //显示合并后的图象
		AftRegDlg* pDlg;
		pDlg = new AftRegDlg(NULL,resultImage);
		pDlg->Create(IDD_DIALOG_NewImage);
		pDlg->ShowWindow(SW_NORMAL);
		pDlg->CenterWindow();
		return;
	}
	else if(m_nImageNum == 3)//拼接图像幅数为3
	{
		list<CString>::iterator iter1 = pathNameList.begin();
		CImage stitchingImage;
		stitchingImage.Load((*iter1));
		iter1 = pathNameList.erase(iter1);
		CImage baseImage;
		baseImage.Load((*iter1));
		iter1 = pathNameList.erase(iter1);
		CImage transitImage;
		if(Imagemosaic(baseImage,stitchingImage,transitImage))
		{
			pathNameList.clear();
			return;				
		}
		stitchingImage.Destroy();
		stitchingImage.Load((*iter1));

		CImage resultImage;
		if(Imagemosaic(transitImage,stitchingImage,resultImage))
		{
			pathNameList.clear();
			return;				
		}
		pathNameList.clear();
		t = ((float)GetTickCount() - t) / 1000;
		CString str;
		str.Format(_T("%f"),t);
		AfxMessageBox(str);
		 //显示合并后的图象
		AftRegDlg* pDlg;
		pDlg = new AftRegDlg(NULL,resultImage);
		pDlg->Create(IDD_DIALOG_NewImage);
		pDlg->ShowWindow(SW_NORMAL);
		pDlg->CenterWindow();
		return;
	}
	else//拼接图像幅数为>=4
	{
		//首先将图像分两组
		int n = pathNameList.size();
		list<CString> list1;
		list<CString> list2;
		int a =0;
		for(list<CString>::iterator iter = pathNameList.begin(); iter != pathNameList.end();)
		{
			a++;
			if(a <= n / 2)
			{
			list1.push_back(*iter);
			iter = pathNameList.erase(iter);
			}
			else
			{
				list2.push_back(*iter);
				iter = pathNameList.erase(iter);
			}
		}
		pathNameList.clear();
		list2.reverse();
		//拼接第一组图像
		CImage resultImage1;
		if(imageGroupMosaic(list1,resultImage1))
		{
			return;
		}
		//拼接第二组图像
		CImage resultImage2;
		if(imageGroupMosaic(list2,resultImage2))
		{
			return;
		}
		//将两组图像的拼接图拼接起来得到最终图像
		CImage resultImage;
		if(Imagemosaic(resultImage2,resultImage1,resultImage))
		{
			return;
		}
		else
		{
			t = ((float)GetTickCount() - t) / 1000;
			CString str;
			str.Format(_T("%f"),t);
			AfxMessageBox(str);
			 //显示合并后的图象
			AftRegDlg* pDlg;
			pDlg = new AftRegDlg(NULL,resultImage);
			pDlg->Create(IDD_DIALOG_NewImage);
			pDlg->ShowWindow(SW_NORMAL);
			pDlg->CenterWindow();
		}
	}
}
//分组图像拼接
bool CImageMosaicDlg::imageGroupMosaic(list<CString> &_list,CImage& resultImage)
{
	bool m_pflag = 0;
	m_nImageNum = _list.size();
	list<CString>::iterator iter1 = _list.begin();
	CImage stitchingImage;
	stitchingImage.Load((*iter1));
	iter1 = _list.erase(iter1);
	CImage baseImage;
	baseImage.Load((*iter1));
	iter1 = _list.erase(iter1);

	if(m_nImageNum == 2)
	{
		m_pflag = 0;
		if(Imagemosaic(baseImage,stitchingImage,resultImage))
		{
			_list.clear();
			return 1;				
		}
		_list.clear();
		return 0;
	}
	else
	{
		CImage transitImage;
		if(Imagemosaic(baseImage,stitchingImage,transitImage))
		{
			_list.clear();
			return 1;
		}
		
		for(;iter1 != _list.end();)
		{
			int a = _list.size();
			m_nImageNum = a + 1;
			baseImage.Destroy();
			baseImage.Load((*iter1));
			iter1 = _list.erase(iter1);
			if(m_nImageNum == 2)
			{
				if(Imagemosaic(baseImage,transitImage,resultImage))
				{
					_list.clear();
					return 1;					
				}
				_list.clear();
				return 0;
			}
			else
			{
				CImage transitImage2;
				stitchingImage.Destroy();
				imageCopy(stitchingImage, transitImage);
				transitImage.Destroy();
				if(Imagemosaic(baseImage,stitchingImage,transitImage))
				{
					_list.clear();
					return 0;
				}
			}
		}
	}
	return 0;
}
//两幅图像拼接
bool CImageMosaicDlg::Imagemosaic(CImage& baseImage,CImage& stitchingImage,CImage& resultImage)
{

	m_nHeight1 = baseImage.GetHeight();
	m_nWidth1 = baseImage.GetWidth();
	m_nHeight2 = stitchingImage.GetHeight();
	m_nWidth2 = stitchingImage.GetWidth();
	if(m_nHeight1 <= 700 || m_nWidth1 <= 700 || m_nHeight2 <= 700 || m_nWidth2 <= 700)
	{
		CHarrisCorners harris0;
		bool* pCorner1 = new bool[m_nHeight1 * m_nWidth1];
		memset(pCorner1,0,m_nHeight1 * m_nWidth1 * sizeof(bool));
		bool* pCorner2 = new bool[m_nHeight2 * m_nWidth2]; 
		memset(pCorner2,0,m_nHeight2 * m_nWidth2 * sizeof(bool));
		float* pGrayValue1 = new float[m_nHeight1 * m_nWidth1];
		float* pGrayValue2 = new float[m_nHeight2 * m_nWidth2];
		getYUV(baseImage,pGrayValue1);
		getYUV(stitchingImage,pGrayValue2);
		harris0.cornersDetecting(pGrayValue1,pGrayValue2,m_nHeight1,m_nWidth1,m_nHeight2,m_nWidth2,pCorner1,pCorner2);
		if(harris0.cornersMatching(pGrayValue1,pGrayValue2,pCorner1,pCorner2))
		{
			return 1;
		}
		harris0.imageFusion(baseImage,stitchingImage,resultImage);

	}
	else 
	{
		if(wvltTransTwiceProcess(baseImage,stitchingImage))
		{
			return 1;
		}
		if(wvltTransOnceProcess(baseImage,stitchingImage))
		{
			return 1;
		}
		if(originImageProcess(baseImage,stitchingImage,resultImage))
		{
			return 1;
		}
	}
	return 0;
}

/******************************************只需要进行两次小波变换*****************************************/
bool CImageMosaicDlg::wvltTransTwiceProcess(CImage& baseImage,CImage& stitchingImage)
{
	//两次小波变换
	CDWT pTrans;
	CImage* pWvltTransTwiceImage1 = new CImage;
	pWvltTransTwiceImage1->Create(m_nWidth1 / 4,m_nHeight1 / 4,24);
	CImage* pWvltTransTwiceImage2 = new CImage;
	pWvltTransTwiceImage2->Create(m_nWidth2 / 4,m_nHeight2 / 4,24);

	float* pGrayValue1 = new float[(m_nWidth1/4)*(m_nHeight1/4)];
	pTrans.WvltTransTwice(baseImage,pGrayValue1);
	float* pGrayValue2 = new float[(m_nWidth2/4)*(m_nHeight2/4)];
	pTrans.WvltTransTwice(stitchingImage,pGrayValue2);
	//对两次小波变换图像进行Harris角点检测
	bool* pWvltTwiceCorner1 = new bool[pWvltTransTwiceImage1->GetHeight() * pWvltTransTwiceImage1->GetWidth()];
	memset(pWvltTwiceCorner1,0,pWvltTransTwiceImage1->GetHeight() * pWvltTransTwiceImage1->GetWidth() * sizeof(bool));
	bool* pWvltTwiceCorner2 = new bool[pWvltTransTwiceImage2->GetHeight() * pWvltTransTwiceImage2->GetWidth()];
	memset(pWvltTwiceCorner2,0,pWvltTransTwiceImage2->GetHeight() * pWvltTransTwiceImage2->GetWidth() * sizeof(bool));
	harris.cornersDetecting(pGrayValue1,pGrayValue2,m_nHeight1/4,m_nWidth1/4,m_nHeight2/4,m_nWidth2/4, pWvltTwiceCorner1, pWvltTwiceCorner2);

	//对两次小波变换图像进行角点匹配
	delete[]pWvltTransTwiceImage1;
	pWvltTransTwiceImage1 = NULL;
	delete[]pWvltTransTwiceImage2;
	pWvltTransTwiceImage2 = NULL;
	if(harris.cornersMatching(pGrayValue1,pGrayValue2,pWvltTwiceCorner1,pWvltTwiceCorner2,m_nHeight1/4,m_nWidth1/4,m_nHeight2/4,m_nWidth2/4,1,0,0))
	{
		return 1;
	}
	return 0;
}
bool CImageMosaicDlg::wvltTransOnceProcess(CImage &baseImage,CImage &stitchingImage)
{
	//一次小波变换
	CDWT pTrans;
	CImage* pWvltTransOnceImage1 = new CImage;
	pWvltTransOnceImage1->Create(m_nWidth1 / 2,m_nHeight1 / 2,24);
	CImage* pWvltTransOnceImage2 = new CImage;
	pWvltTransOnceImage2->Create(m_nWidth2 / 2,m_nHeight2 / 2,24);

	float* pGrayValue1 = new float[m_nHeight1/2 * m_nWidth1/2];
	pTrans.WvltTransOnce(baseImage,pGrayValue1);
	float* pGrayValue2 = new float[m_nHeight2/2 * m_nWidth2/2];
	pTrans.WvltTransOnce(stitchingImage,pGrayValue2);
	//对一次小波变换图像进行Harris角点检测
	bool* pWvltOnceCorner1 = new bool[m_nHeight1/2 * m_nWidth1/2];
	memset(pWvltOnceCorner1,0,m_nHeight1/2 * m_nWidth1/2 * sizeof(bool));
	bool* pWvltOnceCorner2 = new bool[m_nHeight2/2 * m_nWidth2/2];
	memset(pWvltOnceCorner2,0,m_nHeight2/2 * m_nWidth2/2 * sizeof(bool));


	//对一次小波变换图像进行角点匹配
	delete[]pWvltTransOnceImage1;
	pWvltTransOnceImage1 = NULL;
	delete[]pWvltTransOnceImage2;
	pWvltTransOnceImage2 = NULL;
	if(harris.cornersMatching(pGrayValue1,pGrayValue2,pWvltOnceCorner1,pWvltOnceCorner2,m_nHeight1/2,m_nWidth1/2,m_nHeight2/2,m_nWidth2/2 ,0,1,0))
	{
		return 1;
	}
	//harris.imageFusion(pWvltTransOnceImage1,pWvltTransOnceImage2);
	return 0;
}
bool CImageMosaicDlg::originImageProcess(CImage& baseImage,CImage& stitchingImage,CImage& resultImage)
{
	//对原图像进行Harris角点检测
	bool* pOriginCorner1 = new bool[m_nHeight1 *  m_nWidth1];
	memset(pOriginCorner1,0,m_nHeight1 * m_nWidth1 * sizeof(bool));
	bool* pOriginCorner2 = new bool[m_nHeight2 * m_nWidth2];
	memset(pOriginCorner2,0,m_nHeight2 * m_nWidth2 * sizeof(bool));

	float* pGrayValue1 = new float[m_nHeight1 * m_nWidth1];
	getYUV(baseImage,pGrayValue1);//getYUV(CImage* &pNewImage1,float*& P)
	float* pGrayValue2 = new float[m_nHeight2 * m_nWidth2];
	getYUV(stitchingImage,pGrayValue2);

	//对原图像进行角点匹配
	if(harris.cornersMatching(pGrayValue1,pGrayValue2,pOriginCorner1,pOriginCorner2,m_nHeight1,m_nWidth1,m_nHeight2, m_nWidth2,0,0,1))
	{
		return 1;
	}
	harris.imageFusion(baseImage,stitchingImage,resultImage);
	return 0;
}
void CImageMosaicDlg::getYUV(CImage &image,float*& P)
{
	float R = 0;
	float G = 0;
	float B = 0;
	BYTE* pArray1 = (BYTE*)image.GetBits(); //获得数据
	int pitch1 = image.GetPitch();
	int bitCount1 = image.GetBPP() / 8;  //每像素位数
	int nWidth1 = image.GetWidth();
	int nHeight1 = image.GetHeight();

	for (short y = 0; y < nHeight1; y++) 
	{
		for (short x = 0; x < nWidth1; x++) 
		{	
			R = (float)(*(pArray1 + pitch1 * y + x * bitCount1));
			G = (float)(*(pArray1 + pitch1 * y + x * bitCount1 + 1));
			B = (float)(*(pArray1 + pitch1 * y + x * bitCount1 + 2));
			P[y * nWidth1 + x] = 0.299f * R + 0.587f * G + 0.114f * B;
		}
	}
}
void CImageMosaicDlg::Cancel()
{
	// TODO: Add your control notification handler code here
	OnCancel();
}