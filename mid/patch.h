//*===* PATCH.H *==============================================*/
#include <stdio.h>
typedef int *INTPTR;  

class Patch
{
	int ok;
	int NX;
	int NY;
	int X0;
	int Y0;
	int DX;
	int DY;	int CX;
	int CY;
	int NX_patch;
	int NY_patch;
	INTPTR *Inds_map;
	int PatchFrameSize;
	int NorgpixInPatchPix;
	int Open()
	{
		int i;
		Inds_map = NULL;
		Inds_map = new INTPTR[PatchFrameSize];
		if(Inds_map == NULL) return 0;
		for( i = 0; i < PatchFrameSize; i++ ) Inds_map[i] = 0;
		for( i = 0; i < PatchFrameSize; i++ )
		{
			Inds_map[i] = new int[CX*CY];
			if( Inds_map[i] == NULL ) return 0;
		}
		MapIndices();
		return 1;
	};
	void Close()
	{
		int i;
		for( i = 0; i < PatchFrameSize; i++ )
		if( Inds_map[i] != 0 ) delete[] Inds_map[i];
		if( Inds_map != 0 ) delete[] Inds_map;
	};

	void MapIndices(void)
	{
		int ind_patch;
		int ind_org;
		int ind_patch_x;
		int ind_patch_y;
		int ind_x;
		int ind_y;
		int beg_ind_x;
		int end_ind_x;
		int beg_ind_y;
		int end_ind_y;
		int i_map;
		for( ind_patch = 0; ind_patch < PatchFrameSize; ind_patch++ )
		{
			ind_patch_x = ind_patch % NX_patch;
			ind_patch_y = ind_patch / NX_patch;
			beg_ind_x = X0 + CX*ind_patch_x;
			end_ind_x = beg_ind_x + CX;
			beg_ind_y = Y0 + CY*ind_patch_y;
			end_ind_y = beg_ind_y + CY;
			i_map = 0;
			for( ind_x = beg_ind_x; ind_x < end_ind_x; ind_x++ )
				for( ind_y = beg_ind_y; ind_y < end_ind_y; ind_y++ )
				{
					ind_org = ind_y * NX + ind_x;
					Inds_map[ind_patch][i_map++] = ind_org;
				}
		}
	};
	unsigned char getPatchPixel(int pix_ind, unsigned char *buff_in)
	{
		int i;
		int s = 0;
		for( i = 0; i < NorgpixInPatchPix; i++ ) 
            s += buff_in[Inds_map[pix_ind][i]]; 
		return (s/NorgpixInPatchPix);
	};
    double dgetPatchPixel(int pix_ind, double *dbuff_in)
	{
		int i;
		double ds = 0;
		for( i = 0; i < NorgpixInPatchPix; i++ ) 
            ds += dbuff_in[Inds_map[pix_ind][i]]; 
		return (ds/NorgpixInPatchPix);
	};

	public:
	void Convert(unsigned char *frame, unsigned char *patch)
	{
		int i;
		if(!ok) return;
		for( i = 0; i < PatchFrameSize; i++ )
		{
			patch[i] = getPatchPixel(i,frame);
		}
	};


	void dConvert(double *dframe, double *dpatch)
	{
		int i;
		if(!ok) return;
		for( i = 0; i < PatchFrameSize; i++ )
		{
			dpatch[i] = dgetPatchPixel(i,dframe);
		}
	};


	void Print(unsigned char *ptch)
	{
		int i, j;
		for(i = 0; i < NY_patch; i++)
		{
			for(j = 0; j < NX_patch; j++) 
                printf("%2d ",ptch[j + i*NX_patch]);
			printf("\n");
		}
	};
 

	Patch(int nx, int ny, int x0, int y0, int dx, int dy, int cx, int cy)
	{
		NX = nx;
		NY = ny;
		X0 = x0;
		Y0 = y0;
		DX = dx;
		DY = dy;
		CX = cx;
		CY = cy;
		NX_patch = DX/CX;
		NY_patch = DY/CY;
		PatchFrameSize = NX_patch * NY_patch;
		NorgpixInPatchPix = CX*CY;
		ok = Open();
	};

	~Patch()
	{
		Close();
	};
};

