#include "img_convert.h"

// #define RGB2YUV(R,G,B,Y,U,V) \
//     (((Y)=(((R)*19595+(G)*38470+(B)*7471+32767)>>16)), \
//     ((U)=(((-R)*11059-(G)*21709+(B)*32768+32767)>>16)),\
//     ((V)=(((R)*32768-(G)*27439-(B)*5329+32767)>>16)))

#define RGB2YCbCr(R, G, B, Y, Cb, Cr)         \
    (((Y)=(((R)*19595+(G)*38470+(B)*7471+32767)>>16)),\
    ((Cb)=(((-R)*11059-(G)*21709+(B)*32768+32767)>>16)+128),\
    ((Cr)=(((R)*32768-(G)*27439-(B)*5329+32767)>>16)+128))

#define RGB2GRAY(R,G,B,Y) \
    ((Y)=(((R)*19595+(G)*38470+(B)*7471+32767)>>16))


#define YUV2RGB(R,G,B,Y,U,V)                             \
    (((R) = 256 * Y + 358 * (V - 128) >> 8,(R) = sakTrunc(R, 0, 255)),\
	((G) = 256 * Y - 183 * (V - 128) - 87 * (U - 128) >> 8, (G) = sakTrunc(G, 0, 255)),\
	((B) = 256 * Y + 454 * (U - 128) >> 8, (B) = sakTrunc(B, 0, 255)))


#define sakTrunc(V,I,A)	\
	((V)<(I)?(I):((V)>(A)?(A):(V)))


static int BGRBGR2YYUUVV(Image_T &src, Image_T &dst)
{
	int width = src.width;
	int height = src.height;
	unsigned char *pBGR = src.data[0];
	unsigned char *pY = dst.data[0];
	unsigned char *pU = dst.data[1];
	unsigned char *pV = dst.data[2];
	for (int h = 0; h < height; h++)
	{
		for (int w = 0; w < width; w++)
		{
			RGB2YCbCr(pBGR[3 * w + 2], pBGR[3 * w + 1], pBGR[3 * w],
				pY[w], pU[w], pV[w]);
		}
		pBGR += src.pitch[0];
		pY += dst.pitch[0];
		pU += dst.pitch[1];
		pV += dst.pitch[2];
	}
	return 0;
}

static int BGRBGR2GRAY(Image_T &src, Image_T &dst)
{
	int width = src.width;
	int height = src.height;
	unsigned char *pBGR = src.data[0];
	unsigned char *pY = dst.data[0];
	for (int h = 0; h < height; h++)
	{
		for (int w = 0; w < width; w++)
		{
			RGB2GRAY(pBGR[3 * w + 2], pBGR[3 * w + 1], pBGR[3 * w],
				pY[w]);
		}
		pBGR += src.pitch[0];
		pY += dst.pitch[0];
	}
	return 0;
}

static int RGBARGBA2YYUUVV(Image_T &src, Image_T &dst)
{
	int width = src.width;
	int height = src.height;
	unsigned char *pRGBA = src.data[0];
	unsigned char *pY = dst.data[0];
	unsigned char *pU = dst.data[1];
	unsigned char *pV = dst.data[2];
	for (int h = 0; h < height; h++)
	{
		for (int w = 0; w < width; w++)
		{
			RGB2YCbCr(pRGBA[4 * w], pRGBA[4 * w + 1], pRGBA[4 * w + 2],
				pY[w], pU[w], pV[w]);
		}
		pRGBA += src.pitch[0];
		pY += dst.pitch[0];
		pU += dst.pitch[1];
		pV += dst.pitch[2];
	}
	return 0;
}

static int RGBARGBA2GRAY(Image_T &src, Image_T &dst)
{
	int width = src.width;
	int height = src.height;
	unsigned char *pRGBA = src.data[0];
	unsigned char *pY = dst.data[0];
	for (int h = 0; h < height; h++)
	{
		for (int w = 0; w < width; w++)
		{
			RGB2GRAY(pRGBA[4 * w], pRGBA[4 * w + 1], pRGBA[4 * w + 2],
				pY[w]);
		}
		pRGBA += src.pitch[0];
		pY += dst.pitch[0];
	}
	return 0;
}


static int YYUUVV2BGRBGR(Image_T &src, Image_T &dst)
{
	int width = src.width;
	int height = src.height;
	unsigned char *pBGR = dst.data[0];
	unsigned char *pY = src.data[0];
	unsigned char *pU = src.data[1];
	unsigned char *pV = src.data[2];
	for (int h = 0; h < height; h++)
	{
		for (int w = 0; w < width; w++)
		{
			YUV2RGB(pBGR[3 * w + 2], pBGR[3 * w + 1], pBGR[3 * w], \
				pY[w], pU[w], pV[w]);
		}
		pBGR += dst.pitch[0];
		pY += src.pitch[0];
		pU += src.pitch[1];
		pV += src.pitch[2];
	}
	return 0;
}


int img_convert(Image_T &src, Image_T &dst)
{
	if ((IMG_FMT_BGRBGR == src.format) &&
		(IMG_FMT_YYUUVV == dst.format))
		return BGRBGR2YYUUVV(src, dst);

	else if ((IMG_FMT_BGRBGR == src.format) &&
		(IMG_FMT_GRAY == dst.format))
		return BGRBGR2GRAY(src, dst);

	else if ((IMG_FMT_RGBARGBA == src.format) &&
		(IMG_FMT_YYUUVV == dst.format))
		return RGBARGBA2YYUUVV(src, dst);

	else if ((IMG_FMT_RGBARGBA == src.format) &&
		(IMG_FMT_GRAY == dst.format))
		return RGBARGBA2GRAY(src, dst);
	else if ((IMG_FMT_YYUUVV == src.format) &&
		(IMG_FMT_BGRBGR == dst.format))
		return YYUUVV2BGRBGR(src, dst);

	return -1;
}