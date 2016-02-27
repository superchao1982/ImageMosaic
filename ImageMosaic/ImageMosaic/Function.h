#include "math.h"
#include "atlimage.h"
#define PI 3.1415926f
#define P1(ROW,COL) P1[width1 * (COL) + (ROW)]
#define P2(ROW,COL) P2[width2 * (COL) + (ROW)]
float mbys(float*& im,int imW,int imH,int i, int j, float *tp,int tpW,int tpH) 
{
	#define     im(ROW,COL)     im[imW * (COL) + (ROW)]
	#define     tp(ROW,COL)     tp[tpW * (COL) + (ROW)]
	float a = 0;
	for(int n = 0;n < tpH;n++)
	{
		for(int m = 0;m < tpW;m++)
		{
				float b = im(i + m - int(tpW/2),j + n - int(tpH/2));
				float c = tp(m,n);
				a += im(i + m - int(tpW/2),j + n - int(tpH/2)) * tp(m,n);
		}
	}
	return a;
}
float mbys(float* im,int imW,int imH, float *tp,int tpW,int tpH) 
{
	float a = 0;
	for(int n = 0;n < tpH;n++)
	{
		for(int m = 0;m < tpW;m++)
		{
			float b = im[n*tpW + m];
			float c = tp[n*tpW + m];
			a += b * c;
		}
	}
	return a;
}

void prescreening(short height ,short width ,float* &I ,bool* &IsEdge)
{
	#define      num(ROW,COL)       num[width * (COL) + (ROW)]	
	char *num = new char[width * height];
	memset(num,0,width*height*sizeof(char));
	for(int j = 1; j < height-1; j++)
	{
		for(int i = 1; i < width-1; i++)
		{	
			int n = -1;
			for(int b = j-1; b <= j+1; b++)
			{
				for(int a = i-1; a <= i+1; a++)
				{	
					if(abs(I[b * width + a]-I[j * width + i])<=25)
					{
						n++;
					}
				}
			}
			num(i,j) = n;
			if(n == 8 || n == 0 || n==7/* || n == 1*/)
			{
				IsEdge[j * width + i] = 0;
			}
			else
			{
				IsEdge[j * width + i] = 1;
			}
		}
	}
	for(int j = 2; j < height-2; j++)
	{
		for(int i = 2; i < width-2; i++)
		{	
			if(IsEdge[j * width + i] == 1)
			{
				int n = 0;
				for(int b = j - 1;b <= j + 1;b++)
				{
					for(int a = i - 1;a <= i + 1;a++)
					{
						if(a != i&& b != j)
						{
							if(num(a,b) < num(i,j) && (num(a,b) != 0 || num(a,b) != 1))
							{
								IsEdge[j * width + i] = 0;
							}
							else
							{
								//if(num(a,b) >= 5 && num(i,j) == 5)
								//{
								//	n++;
								//}
							}
						}
					}
				}
				//if(n == 8&&((num(i-1,j) == 5&&num(i+1,j) == 5)||(num(i-1,j-1) == 5&&num(i+1,j+1) == 5)||(num(i-1,j+1) == 5&&num(i+1,j-1) == 5)))
				//{
				//	IsEdge[j * width + i] = 0;
				//}
			}
		}
	}
	delete[]num;
	num = 0;
}
float corresX(unsigned char P1[],unsigned char P2[],short nCount)
{
	float sum1 = 0.0f;
	float sum2 = 0.0f;
	float sum3 = 0.0f;
	float sum4 = 0.0f;
	float sum5 = 0.0f;
	float u1 = 0.0f;
	float u2 = 0.0f;
	for(short k = 0;k < nCount;k++)
	{
		sum1 += (float)P1[k];
		sum2 += (float)P2[k];
	}
	u1 = sum1 / nCount;
	u2 = sum2 / nCount;
	for(short k = 0;k < nCount;k++)
	{
		sum3 += (P1[k]-u1) * (P2[k]-u2);
		sum4 += (P1[k]-u1) * (P1[k]-u1);
		sum5 += (P2[k]-u2) * (P2[k]-u2);
	}
	float R = fabs(sum3 / (sqrt(sum4) * sqrt(sum5)));
	return R;

}
/********************************************************************************************/
/*****************************************D7相关系数计算*************************************/
/********************************************************************************************/
float Corres(float *&P1,float *&P2,short width1,short width2,short i,short j,short m,short n)
{
	unsigned char pBuffer1[37];
	unsigned char pBuffer2[37];
	short c = 0;
	short d = 0;
	short e = 0;
	short f = 0;
	short r = 0;
	short g = -1;
	for(f = j - 3,d = n - 3;f <= j + 3,d <= n + 3;f++,d++)
	{
		if((f == j - 3 && d == n - 3) || (f == j + 3 && d == n + 3))
		{
			for(e = i - 1,c = m - 1;e <= i + 1,c <= m + 1;e++,c++)
			{
				g++;
				if(P1(e,f) > 255)
				{
					P1(e,f) = 255;
				}
				if(P1(e,f) < 0)
				{
					P1(e,f) = 0;
				}
				pBuffer1[g]   = (unsigned char)P1(e,f);
				if(P2(c,d) > 255)
				{
					P2(c,d) = 255;
				}
				if(P2(c,d) < 0)
				{
					P2(c,d) = 0;
				}
				pBuffer2[g]   = (unsigned char)P2(c,d);
			}
		}
		else if((f == j - 2 && d == n - 2) || (f == j + 2 && d == n + 2))
		{
			for(e = i - 2,c = m - 2;e <= i + 2,c <= m + 2;e++,c++)
			{
				g++;
				if(P1(e,f) > 255)
				{
					P1(e,f) = 255;
				}
				if(P1(e,f) < 0)
				{
					P1(e,f) = 0;
				}
				pBuffer1[g]   = (unsigned char)P1(e,f);
				if(P2(c,d) > 255)
				{
					P2(c,d) = 255;
				}
				if(P2(c,d) < 0)
				{
					P2(c,d) = 0;
				}
				pBuffer2[g]   = (unsigned char)P2(c,d);
			}
		}
		else
		{
			for(e = i - 3,c = m - 3;e <= i + 3,c <= m + 3;e++,c++)
			{
				g++;
				if(P1(e,f) > 255)
				{
					P1(e,f) = 255;
				}
				if(P1(e,f) < 0)
				{
					P1(e,f) = 0;
				}
				pBuffer1[g]   = (unsigned char)P1(e,f);
				if(P2(c,d) > 255)
				{
					P2(c,d) = 255;
				}
				if(P2(c,d) < 0)
				{
					P2(c,d) = 0;
				}
				pBuffer2[g]   = (unsigned char)P2(c,d);
			}
		}
	}

	float R = corresX(pBuffer1,pBuffer2,37);
	return R;
}
void GetM(float *&P1,short width1,short i,short j,unsigned char M[])
{
	short e = 0;
	short f = 0;
	short r = 0;
	short g = -1;
	for(f = j - 3;f <= j + 3;f++)
	{
		if((f == j - 3) || (f == j + 3))
		{
			for(e = i - 1;e <= i + 1;e++)
			{
				g++;
				if(P1(e,f) > 255)
				{
					P1(e,f) = 255;
				}
				if(P1(e,f) < 0)
				{
					P1(e,f) = 0;
				}
				M[g]   = (unsigned char)P1(e,f);
			}
		}
		else if((f == j - 2) || (f == j + 2))
		{
			for(e = i - 2;e <= i + 2;e++)
			{
				g++;
				if(P1(e,f) > 255)
				{
					P1(e,f) = 255;
				}
				if(P1(e,f) < 0)
				{
					P1(e,f) = 0;
				}
				M[g]   = (unsigned char)P1(e,f);
			}
		}
		else
		{
			for(e = i - 3;e <= i + 3;e++)
			{
				g++;
				if(P1(e,f) > 255)
				{
					P1(e,f) = 255;
				}
				if(P1(e,f) < 0)
				{
					P1(e,f) = 0;
				}
				M[g]   = (unsigned char)P1(e,f);
			}
		}
	}
}
float max3v(float r,float g,float b)
{
	float max = r;
	if(g > max)
	{
		max = g;
	}
	if(b > max)
	{
		max = b;
	}
	return max;
}
float min3v(float r,float g,float b)
{
	float min = r;
	if(g < min)
	{
		min = g;
	}
	if(b < min)
	{
		min = b;
	}
	return min;
}
void RGBtoYUV(float R,float G,float B,float &H,float &S,float &L)
{
	float h = 0;
	float s = 0;
	float l = 0;
	float r = R / 255.0f;
	float g = G / 255.0f;
	float b = B / 255.0f;
	float max = max3v(r,g,b);
	float min = min3v(r,g,b);
	if(max == min)
	{
		h = 0;
	}
	else if(max == r && g >= b)
	{
		h = 60.0f * (g - b) / (max - min);
	}
	else if(max = r && g < b)
	{
		h = 60.0f * (g - b)/(max - min) + 360.0f;
	}
	else if(max == g)
	{
		h = 60.0f * (b - r)/(max - min) + 120.0f;
	}
	else if(max == b)
	{
		h = 60.0f * (r - g) / (max - min) + 240.0f; 
	}
	l = (max + min) / 2.0f;
	if(l == 0 || max == min)
	{
		s = 0;
	}
	else if(0 < l && l <= 0.5f)
	{
		s = (max - min) / (max + min);
	}
	else if(l > 0.5f)
	{
		s = (max - min) / (2 - (max + min));
	}
	H = (h > 360)?360:((h<0)?0:h);
	S = ((s>1)?1:((s<0)?0:s)) * 100;
	L = ((l>1)?1:((l<0)?0:l)) * 100;
}
void YUVtoRGB(float H,float S,float L,float &R,float &G,float &B)
{
	float h = H;
	float s = S / 100.0f;
	float l = L / 100.0f;
	float r = 0;
	float g = 0;
	float b = 0;
	if(S == 0)
	{
		r = g = b = l * 255.0f;
	}
	else
	{
		float q = (l < 0.5f)?(l*(1.0f+s)):(l+s-(l*s));
		float p = 2.0f * l - q;
		float Hk = h / 360.0f;
		float T[3];
		T[0] = Hk + 0.3333333f;
		T[1] = Hk;
		T[2] = Hk - 0.3333333f;
		for(int i = 0; i < 3;i++)
		{
			if(T[i] < 0)
			{
				T[i] += 1.0f;
			}
			if(T[i] > 1)
			{
				T[i] -= 1.0f;
			}
			if((T[i]*6) < 1)
			{
				T[i] = p +((q - p) * 6.0f * T[i]);
			}
			else if((T[i]*2.0f) < 1)
			{
				T[i] = q;
			}
			else if((T[i]*3.0f) < 2)
			{
				T[i] = p + (q - p) * ((2.0f/3.0f) - T[i]) * 6.0f;
			}
			else
			{
				T[i] = p;
			}
		}
		r = T[0] * 255.0f;
		g = T[1] * 255.0f;
		b = T[2] * 255.0f;
	}
	R = r > 255 ? 255 : (r <0 ? 0 : r);
	G = g > 255 ? 255 : (g <0 ? 0 : g);
	B = b > 255 ? 255 : (b <0 ? 0 : b);
}