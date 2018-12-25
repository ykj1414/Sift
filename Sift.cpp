#include "Sift.h"
#include <stdio.h>
#include <math.h>
#include <malloc.h>



float sigma_ary[LAYER_NUM + 3];
int mat_size[LAYER_NUM + 3];
float sigma15_ary[LAYER_NUM];
int mat15_size[LAYER_NUM];
float *cur_offset = (float*)calloc(3,sizeof(float));
float rotate_angle;

void buildGaussianPyramid(Mat GrayImage, Mat ***GaussianPyramid, int image_width, int image_height, int octave_num = 4, float sigma = 1.6);

void buildDogPyramid(Mat***GaussianPyramid, Mat ***DogPyramid, int octave_num = 4);

void findExtemePoint(Mat ***GaussianPyramid,Mat***DogPyramid,std::vector<SFKP>&siftpoints, int octave_num = 4, float diff_threshold = DIFF_THRESHOLD);

bool judgeLocalExtreme(Mat***DogPyramid, int *row, int *column, int noctave, int *layer, float diff_threshold = DIFF_THRESHOLD);

void findOrientation(Mat ***GaussianPyramid,std::vector<SFKP>&siftpoints, int image_height, int image_width,int noctave,int layer,int x,int y);

float ***siftDescriptor(float **Magnitude_matrix, float **Orientation_matrix, float orientation, int x, int y, int in_radius, int out_radius);

float* judgeRange(float degree);

float tri_linear(float center_x, float center_y, int point_x, int point_y, float distance_db);

void buildGaussianPyramid(Mat GrayImage, Mat *** GaussianPyramid, int image_width, int image_height, int octave_num, float sigma)
{
	sigma_ary[0] = sigma;// sqrt(sigma*sigma - SIGMA_INITIAL*SIGMA_INITIAL);
	float k = pow(2, 1.f / LAYER_NUM);
	if (cvCeil(6 * sigma_ary[0] + 1) %2 == 0)
	{
		mat_size[0] = cvCeil(6 * sigma_ary[0] + 1) + 1;
	}
	else
		mat_size[0] = cvCeil(6 * sigma_ary[0] + 1);
	for (int i = 1;i<LAYER_NUM+3;i++)
	{
		sigma_ary[i] = pow(k, i - 1)*sigma_ary[0] * sqrt(k*k - 1.f);
		if (cvCeil(6 * sigma_ary[i] + 1) % 2 == 0)
		{
			mat_size[i] = cvCeil(6 * sigma_ary[i] + 1) + 1;
		}
		else
			mat_size[i] = cvCeil(6 * sigma_ary[i] + 1);
	}
	for (int i = 0;i<LAYER_NUM;i++)
	{
		sigma15_ary[i] = sigma_ary[i + 1] * 1.5;
		if (cvCeil(6 * sigma15_ary[i] + 1) % 2 == 0)
		{
			mat15_size[i] = cvCeil(6 * sigma15_ary[i] + 1) + 1;
		}
		else
			mat15_size[i] = cvCeil(6 * sigma15_ary[i] + 1);
	}
	for (int i = 0;i<octave_num;i++)
	{
		if (i == 0)
		{
			Mat_<float>&dst = (Mat_<float>)***GaussianPyramid;
			const Mat_<float>&src = GrayImage;
			Mat_<float> d_src;
			resize(src, d_src, Size(src.cols << 1, src.rows << 1), 0, 0, INTER_LINEAR);
			float sig_negtive = sqrtf(std::max(sigma*sigma - SIGMA_INITIAL*SIGMA_INITIAL*4 , 0.8));
			GaussianBlur(d_src,dst, Size(), sig_negtive, sig_negtive);
			//resize(d_src, dst, Size(d_src.rows>>1, d_src.cols>>1), 0, 0, INTER_NEAREST);
			//for (int x = 0;x<dst.rows;x++)
			//{
			//	float *p = GaussianPyramid[0][0]->ptr<float>(x);
			//	for (int y = 0;y<dst.cols;y++)
			//	{
			//		*p++ = cvFloor(*p);
			//	}
			//}
		}
		else {
			const Mat_<float>& src = **(*(GaussianPyramid+i-1)+LAYER_NUM);
			Mat_<float> &dst = (Mat_<float>)**(*(GaussianPyramid+i)+0);
			resize(src, dst, Size(src.cols>>1,src.rows>>1),0,0,INTER_NEAREST);
		}			
		for (int j = 1;j<LAYER_NUM+3;j++)
		{
			const Mat_<float>&src = **(*(GaussianPyramid+i)+j-1);
			Mat_<float> &dst = (Mat_<float> )**(*(GaussianPyramid+i)+j);
			GaussianBlur(src,dst, Size(mat_size[j],mat_size[j]), sigma_ary[j],sigma_ary[j]);
			//for (int x = 0; x < dst.rows; x++)
			//{
			//	float *p = GaussianPyramid[i][j]->ptr<float>(x);
			//	for (int y = 0; y < dst.cols; y++)
			//	{
			//		*p++ = cvFloor(*p);
			//	}
			//}
		}
	}
}

void buildDogPyramid(Mat *** GaussianPyramid, Mat *** DogPyramid, int octave_num)
{
	for (int i = 0;i<octave_num;i++)
	{
		for (int j = 0;j<LAYER_NUM+2;j++)
		{
			const Mat_<float> &src1 = *GaussianPyramid[i][j];
			const Mat_<float> &src2 = *GaussianPyramid[i][j + 1];
			Mat_<float> &dst = (Mat_<float>)*DogPyramid[i][j];
			subtract(src2, src1, dst);
			/*for (int k = 0;k<DogPyramid[i][j]->rows;k++)
			{
				const float *sp1 = src1.ptr<float>(k);
				const float *sp2 = src2.ptr<float>(k);
				float *dp = dst.ptr<float>(k);
				for ( int  m = 0;m<DogPyramid[i][j]->cols;m++)
				{
					*dp++ = *sp1++ - *sp2++;
				}
			}*/
		}
	}
}

void findExtemePoint(Mat*** GaussianPyramid,Mat *** DogPyramid,std::vector<SFKP>&siftpoints, int octave_num, float diff_threshold)
{
	int sum = 0;
	for (int i = 0;i<octave_num;i++)
	{
		for (int j = 1;j<LAYER_NUM+1;j++)
		{
			const Mat_<float>&pre = *DogPyramid[i][j - 1];
			const Mat_<float> &cur = *DogPyramid[i][j];
			const Mat_<float> &next = *DogPyramid[i][j + 1];
			for (int x = 5;x<cur.rows-5;x++)
			{
				const float*preptr = pre.ptr<float>(x);
				const float*preptrup = pre.ptr<float>(x - 1);
				const float*preptrdown = pre.ptr<float>(x + 1);
				const float*curptr = cur.ptr<float>(x);
				const float*curptrup = cur.ptr<float>(x - 1);
				const float*curptrdown = cur.ptr<float>(x + 1);
				const float*nextptr = next.ptr<float>(x);
				const float*nextptrup = next.ptr<float>(x - 1);
				const float*nextptrdown = next.ptr<float>(x + 1);
				for (int y = 5; y < cur.cols - 5; y++)
				{
					if (abs(curptr[y])>=diff_threshold*255)
					{
						if ((curptr[y]>0 && curptr[y]>=curptr[y-1] && curptr[y]>=curptr[y+1]&&
							curptr[y]>=curptrup[y-1]&&curptr[y]>=curptrup[y]&&curptr[y]>=curptrup[y+1]&&
							curptr[y]>=curptrdown[y-1]&&curptr[y]>=curptrdown[y]&&curptr[y]>=curptrdown[y+1]&&
							curptr[y]>=preptrdown[y-1]&&curptr[y]>=preptrdown[y]&&curptr[y]>=preptrdown[y+1]&&
							curptr[y]>=preptr[y-1]&&curptr[y]>=preptr[y]&&curptr[y]>=preptr[y+1]&&
							curptr[y]>=preptrup[y-1]&&curptr[y]>=preptrup[y]&&curptr[y]>=preptrup[y+1]&&
							curptr[y]>=nextptr[y]&&curptr[y]>=nextptr[y-1]&&curptr[y]>=nextptr[y+1]&&
							curptr[y]>=nextptrup[y]&&curptr[y]>=nextptrup[y-1]&&curptr[y]>=nextptrup[y+1]&&
							curptr[y]>=nextptrdown[y]&&curptr[y]>=nextptrdown[y-1]&&curptr[y]>=nextptrdown[y+1]) || 
							(curptr[y] < 0 && curptr[y] <= curptr[y - 1] && curptr[y] <= curptr[y + 1] && 
								curptr[y] <= curptrup[y - 1] && curptr[y] <= curptrup[y] && curptr[y] <= curptrup[y + 1] && 
								curptr[y] <= curptrdown[y - 1] && curptr[y] <= curptrdown[y] && curptr[y] <= curptrdown[y + 1] && 
								curptr[y] <= preptrdown[y - 1] && curptr[y] <= preptrdown[y] && curptr[y] <= preptrdown[y + 1] && 
								curptr[y] <= preptr[y - 1] && curptr[y] <= preptr[y] && curptr[y] <= preptr[y + 1] && 
								curptr[y] <= preptrup[y - 1] && curptr[y] <= preptrup[y] && curptr[y] <= preptrup[y + 1] && 
								curptr[y] <= nextptr[y] && curptr[y] <= nextptr[y - 1] && curptr[y] <= nextptr[y + 1] && 
								curptr[y] <= nextptrup[y] && curptr[y] <= nextptrup[y - 1] && curptr[y] <= nextptrup[y + 1] && 
								curptr[y] <= nextptrdown[y] && curptr[y] <= nextptrdown[y - 1] && curptr[y] <= nextptrdown[y + 1]))
						{
							int x1 = x; int y1 = y; int r = j; int o = i;
							if (!judgeLocalExtreme(DogPyramid,&x1, &y1, o, &r, diff_threshold))
								continue;
							else
							{
								findOrientation(GaussianPyramid, siftpoints, cur.rows, cur.cols, o, r, x1, y1);
								sum++;
							}
						}
					}
				}
			}
		}
	}
	printf("%d\n", sum);
}

bool judgeLocalExtreme(Mat ***DogPyramid, int * row, int * column,int noctave, int * layer, float diff_threshold)
{
	int xborder = DogPyramid[noctave][*layer]->rows-1;
	int yborder = DogPyramid[noctave][*layer]->cols-1;
	int i = 0;
	while (i<5)
	{
		i++;
		float cur2 = (float)DogPyramid[noctave][*layer]->at<float>(*row, *column) * 2.f;
		float  dx = (float)(DogPyramid[noctave][*layer]->at<float>(*row + 1, *column) - DogPyramid[noctave][*layer]->at<float>(*row - 1, *column)) / (2.f * 255);
		float  dy = (float)(DogPyramid[noctave][*layer]->at<float>(*row, *column+1) - DogPyramid[noctave][*layer]->at<float>(*row, *column-1)) / (2.f * 255);
		float  dr = (float)(DogPyramid[noctave][*layer + 1]->at<float>(*row, *column) - DogPyramid[noctave][*layer-1]->at<float>(*row, *column)) /(2.f*255);
		float dxx = (float)(DogPyramid[noctave][*layer]->at<float>(*row + 1, *column) + DogPyramid[noctave][*layer]->at<float>(*row -1, *column) - cur2) / 255.f;
		float dyy = (float)(DogPyramid[noctave][*layer]->at<float>(*row, *column+1) + DogPyramid[noctave][*layer]->at<float>(*row, *column-1) - cur2) / 255.f;
		float drr = (float)(DogPyramid[noctave][*layer+1]->at<float>(*row, *column) + DogPyramid[noctave][*layer-1]->at<float>(*row , *column) - cur2) / 255.f;
		float dxy = (float)(DogPyramid[noctave][*layer]->at<float>(*row + 1, *column + 1) + DogPyramid[noctave][*layer]->at<float>(*row - 1, *column - 1) -
			DogPyramid[noctave][*layer]->at<float>(*row + 1, *column - 1) - DogPyramid[noctave][*layer]->at<float>(*row - 1, *column + 1))/(4.f*255);
		float dxr = (float)(DogPyramid[noctave][*layer + 1]->at<float>(*row + 1, *column) - DogPyramid[noctave][*layer + 1]->at<float>(*row - 1, *column) +
			DogPyramid[noctave][*layer - 1]->at<float>(*row + 1, *column) - DogPyramid[noctave][*layer - 1]->at<float>(*row - 1, *column)) / (4.f * 255);
		float dyr = (float)(DogPyramid[noctave][*layer + 1]->at<float>(*row, *column+1) - DogPyramid[noctave][*layer + 1]->at<float>(*row, *column-1) +
			DogPyramid[noctave][*layer - 1]->at<float>(*row, *column+1) - DogPyramid[noctave][*layer - 1]->at<float>(*row, *column-1)) / (4.f * 255);
		float det = dxx*dyy - dxy*dxy;
		float tr = dxx + dyy;
		Vec3f d_once(dy, dx, dr);
		Matx33f d_twice(dyy, dxy, dyr,dxy, dxx, dxr, dyr, dxr, drr);
		Vec3f X = d_twice.solve(d_once, DECOMP_LU);
		Vec3f negX(-X[0], -X[1], -X[2]);
		float delta = d_once.dot(negX);
		if (abs(X[0])<0.5&&abs(X[1])<0.5&&abs(X[2])<0.5)
		{
			if (abs((float)DogPyramid[noctave][*layer]->at<float>(*row, *column) / 255.f + delta / 2.f) < diff_threshold)
			{
				return false;
			}
			if (tr*tr / det < (EDGE_THRESHOLD + 1)*(EDGE_THRESHOLD + 1) / EDGE_THRESHOLD&&det>0)
			{
				float *f_offset = (float*)calloc(3,sizeof(float));
				f_offset[0] = -X[1];
				f_offset[1] = -X[0];
				f_offset[2] = -X[2];
				cur_offset = f_offset;
				return true;
			}
			else
				return false;
		}
		*row = *row + cvRound(-X[1]);
		*column = *column + cvRound(-X[0]);
		*layer = *layer + cvRound(-X[2]);
		if (*row<1 || *row>xborder - 1 || *column<1 || *column>yborder - 1 || *layer<1 || *layer>LAYER_NUM)
		{
			return false;
		}
	}
	return false;
}

void findOrientation(Mat ***GaussianPyramid, std::vector<SFKP> &siftpoints, int image_height,int image_width, int noctave, int layer, int x, int y)
{
	{
		int o = noctave;
		int r = layer-1;
		{
			float **Magnitude_Matrix = (float**)calloc(mat15_size[r], sizeof(float*));
			float **Orientation_Matrix = (float**)calloc(mat15_size[r], sizeof(float*));
			for (int i = 0; i < mat15_size[r]; i++)
			{
				Magnitude_Matrix[i] = (float*)calloc(mat15_size[r], sizeof(float));
				Orientation_Matrix[i] = (float*)calloc(mat15_size[r], sizeof(float));
			}
			int in_radius = (mat_size[r+ 1] - 1) >> 1;
			int out_radius = (mat15_size[r] - 1) >> 1;
			for (int x1 = x - out_radius; x1 < x + out_radius + 1; x1++)
			{
				if (x1 <= 0 || x1 >= image_height- 1)
					continue;
				int x_tran = x1 - x + out_radius;
				for (int y1 = y - out_radius; y1 < y + out_radius + 1; y1++)
				{
					if (y1 <= 0 || y1 >= image_width  - 1)
						continue;
					else
					{
						int y_tran = y1 - y + out_radius;
						float dx = (float)GaussianPyramid[o][r + 1]->at<float>(x1 + 1, y1) - (float)GaussianPyramid[o][r + 1]->at<float>(x1 - 1, y1);
						float dy = (float)GaussianPyramid[o][r + 1]->at<float>(x1, y1 + 1) - (float)GaussianPyramid[o][r + 1]->at<float>(x1, y1 - 1);
						float gweight = 1.f/(2.f*Pi*sigma15_ary[r]*sigma15_ary[r])*exp(-((x1 - x)*(x1 - x) + (y1 - y)*(y1 - y)) / (sigma15_ary[r] * sigma15_ary[r]* 2.f));
						Magnitude_Matrix[x_tran][y_tran] = sqrtf(dx*dx + dy*dy)*gweight;
						Orientation_Matrix[x_tran][y_tran] = fastAtan2(dx, dy);
						//printf("%f,%f\n", Magnitude_Matrix[x_tran][y_tran], Orientation_Matrix[x_tran][y_tran]);
					}
				}
			}
			float *histogram = (float*)calloc(36, sizeof(float));
			for (int i = 0; i < mat15_size[r] - 1; i++)
			{
				for (int j = 0; j < mat15_size[r] - 1; j++)
				{
					int index = Orientation_Matrix[i][j] / 10;
					histogram[index] = histogram[index] + Magnitude_Matrix[i][j];
				}
			}
			int total_max = 1;
			float *max = (float*)calloc(total_max, sizeof(float));
			int *ort_index = (int*)calloc(total_max, sizeof(int));
			for (int i = 0; i < 36; i++)
			{
				printf("%f\t", histogram[i]);
				if (max[total_max - 1] < histogram[i])
				{
					max[total_max - 1] = histogram[i];
					ort_index[total_max - 1] = i + 1;
				}
			}
			for (int i = 0; i < 36; i++)
			{
				if (histogram[i] > 0.8*max[0] && i != ort_index[0] - 1)
				{
					max = (float*)realloc(max, ++total_max * sizeof(float));
					ort_index = (int*)realloc(ort_index, total_max * sizeof(int));
					max[total_max - 1] = histogram[i];
					ort_index[total_max - 1] = i + 1;
				}
				//printf("%f\t", histogram[i]);
			}
			/*for (int m = 0; m < out_radius * 2 + 1; m++)
			{
				for (int n = 0; n < out_radius * 2 + 1; n++)
				{
					float gweight = 1.f / (2.f*Pi*sigma15_ary[r] * sigma15_ary[r])*exp(-((m - out_radius)*(m - out_radius) + (n - out_radius)*(n - out_radius)) / (sigma15_ary[r] * sigma15_ary[r] * 2.f));
					Magnitude_Matrix[m][n] = gweight*Magnitude_Matrix[m][n];
				}
			}*/
			for (int i = 0; i < total_max; i++)
			{
				SFKP *s1 = (SFKP*)malloc(sizeof(SFKP));
				s1->octave_num = o;
				s1->layer_num = r;
				s1->column = o>0?y<<(o-1):y>>1;
				s1->row = o>0?x<<(o-1):x>>1;
				s1->magnitude = max[i];
				s1->offset = cur_offset;
				float left_his = histogram[(ort_index[i] - 2) % 36];
				float right_his = histogram[(ort_index[i]) % 36];
				s1->orientation = 10.f*((float)ort_index[i]  -  (right_his - left_his) / (left_his + right_his - 2.f * histogram[ort_index[i] - 1]));
				s1->descripter = siftDescriptor(Magnitude_Matrix, Orientation_Matrix, s1->orientation, x, y, in_radius, out_radius);
				float p_sum = 0.f;
				for (int i = 0; i < 4; i++)
				{
					for (int j = 0; j < 4; j++)
					{
						for (int k = 0; k < 8; k++)
						{
							p_sum += s1->descripter[i][j][k];
						}
					}
				}
				for (int i = 0; i < 4; i++)
				{
					for (int j = 0; j < 4; j++)
					{
						for (int k = 0; k < 8; k++)
						{
							s1->descripter[i][j][k] /= p_sum;
						}
					}
				}
				siftpoints.push_back(*s1);
			}
			for (int i = 0; i < mat15_size[r]; i++)
			{
				free(Magnitude_Matrix[i]);
				free(Orientation_Matrix[i]);
			}
			free(Magnitude_Matrix);
			free(Orientation_Matrix);
			free(ort_index);
			free(max);
			free(histogram);
		}
	}
}

float *** siftDescriptor(float ** Magnitude_matrix, float ** Orientation_matrix, float orientation, int x, int y, int in_radius, int out_radius)
{
	int midy = out_radius;
	int midx = midy;
	int sx = midx - in_radius;
	int sy = sx;
	float ***result = (float***)malloc(sizeof(float**) * 4);
	for (int i = 0 ;i<4;i++)
	{
		result[i] = (float**)malloc(sizeof(float*) * 4);
		for (int j = 0;j<4;j++)
		{
			result[i][j] = (float*)calloc(8, sizeof(float));
		}
	}
	float **Magnitude_trans = (float**)calloc((2 * in_radius + 1), sizeof(float));
	float **Orientation_trans = (float**)calloc((2 * in_radius + 1), sizeof(float));
	for (int i =0;i<2*in_radius+1;i++)
	{
		Magnitude_trans[i] = (float*)calloc(2 * in_radius + 1, sizeof(float));
		Orientation_trans[i] = (float*)calloc(2 * in_radius + 1, sizeof(float));
		for (int j = 0;j<2*in_radius+1;j++)
		{
			int dx = sx + i - midx;
			int dy = sy + j - midy;
			double ag_sin = sin((360-orientation)*Pi / 180.f);
			double ag_cos = cos((360-orientation)*Pi / 180.f);
			int rowindex = midx + cvRound(dy*ag_sin + dx*ag_cos);
			int colindex = midy + cvRound(dy*ag_cos - dx*ag_sin);
			Magnitude_trans[i][j] = Magnitude_matrix[midx + cvRound(dy*ag_sin + dx*ag_cos)][midy + cvRound(dy*ag_cos - dx*ag_sin)];
			Orientation_trans[i][j] = Orientation_matrix[midx + cvRound(dy*ag_sin + dx*ag_cos)][midy + cvRound(dy*ag_cos - dx*ag_sin)];
		}
	}
	float db = in_radius / 2.f;
	for (int i = 1;i<=2*in_radius+1;i++)
	{
		for (int j = 1;j<=2*in_radius+1;j++)
		{
			float *index_weight = judgeRange(Orientation_trans[i-1][j-1]);
			int left_side = index_weight[0];
			int right_side = left_side % 8 + 1;
			float left_weight = index_weight[1];
			float right_weight = 1.f - left_weight;
			free(index_weight);
			int px = cvCeil(i / db - 0.5)+1;
			int py = cvCeil(j / db - 0.5)+1;
			if (px == 1 && py == 1)
			{
				float weight00 = tri_linear(db*(px - 1) + db / 2, db*(py - 1) + db / 2, i, j, db);
				result[0][0][left_side - 1] = weight00*left_weight*Magnitude_trans[i - 1][j - 1];
				result[0][0][right_side - 1] = weight00*right_weight*Magnitude_trans[i - 1][j - 1];
			}
			else if (px == 1 && py == 5)
			{
				float weight00 = tri_linear(db*(px - 1) + db / 2, db*(py - 2) + db / 2, i, j, db);
				result[0][3][left_side - 1] = weight00*left_weight*Magnitude_trans[i - 1][j - 1];
				result[0][3][right_side - 1] = weight00*right_weight*Magnitude_trans[i - 1][j - 1];
			}
			else if (px == 5 && py == 1)
			{
				float weight00 = tri_linear(db*(px - 2) + db / 2, db*(py - 1) + db / 2, i, j, db);
				result[3][0][left_side - 1] = weight00*left_weight*Magnitude_trans[i - 1][j - 1];
				result[3][0][right_side - 1] = weight00*right_weight*Magnitude_trans[i - 1][j - 1];
			}
			else if (px == 5 && py == 5)
			{
				float weight00 = tri_linear(db*(px - 2) + db / 2, db*(py - 2) + db / 2, i, j, db);
				result[3][3][left_side - 1] = weight00*left_weight*Magnitude_trans[i - 1][j - 1];
				result[3][3][right_side - 1] = weight00*right_weight*Magnitude_trans[i - 1][j - 1];
			}
			else if (px==1&&py>1&&py<5)
			{
				float weight01 = tri_linear(db*(px - 1) + db / 2, db*(py - 1) + db / 2, i, j, db);
				result[0][py - 1][left_side - 1] = weight01*left_weight*Magnitude_trans[i - 1][j - 1];
				result[0][py - 1][right_side - 1] = weight01*right_weight*Magnitude_trans[i - 1][j - 1];
				float weight00 = tri_linear(db*(px - 1) + db / 2, db*(py - 2) + db / 2, i, j, db);
				result[0][py - 2][left_side - 1] = weight00*left_weight*Magnitude_trans[i - 1][j - 1];
				result[0][py - 2][right_side - 1] = weight00*right_weight*Magnitude_trans[i - 1][j - 1];
			}
			else if (px==5 &&py>1&&py<5)
			{
				float weight01 = tri_linear(db*(px - 2) + db / 2, db*(py - 1) + db / 2, i, j, db);
				result[3][py - 1][left_side - 1] = weight01*left_weight*Magnitude_trans[i - 1][j - 1];
				result[3][py - 1][left_side - 1] = weight01*right_weight*Magnitude_trans[i - 1][j - 1];
				float weight00 = tri_linear(db*(px - 2) + db / 2, db*(py - 2) + db / 2, i, j, db);
				result[3][py - 2][left_side - 1] = weight00*left_weight*Magnitude_trans[i - 1][j - 1];
				result[3][py - 2][right_side - 1] = weight00*right_weight*Magnitude_trans[i - 1][j - 1];
			}
			else if (py==5&&px>1&&px<5)
			{
				float weight01 = tri_linear(db*(px - 1) + db / 2, db*(py - 2) + db / 2, i, j, db);
				result[px - 1][3][left_side - 1] = weight01*left_weight*Magnitude_trans[i - 1][j - 1];
				result[px - 1][3][right_side - 1] = weight01*right_weight*Magnitude_trans[i - 1][j - 1];
				float weight00 = tri_linear(db*(px - 2) + db / 2, db*(py - 2) + db / 2, i, j, db);
				result[px - 2][3][left_side - 1] = weight00*left_weight*Magnitude_trans[i - 1][j - 1];
				result[px - 2][3][right_side - 1] = weight00*right_weight*Magnitude_trans[i - 1][j - 1];
			}
			else if (py==1&&px>1&&px<5)
			{
				float weight01 = tri_linear(db*(px - 1) + db / 2, db*(py - 1) + db / 2, i, j, db);
				result[px - 1][0][left_side - 1] = weight01*left_weight*Magnitude_trans[i - 1][j - 1];
				result[px - 1][0][right_side - 1] = weight01*right_weight*Magnitude_trans[i - 1][j - 1];
				float weight00 = tri_linear(db*(px - 2) + db / 2, db*(py - 1) + db / 2, i, j, db);
				result[px - 2][0][left_side - 1] = weight00*left_weight*Magnitude_trans[i - 1][j - 1];
				result[px - 2][0][right_side - 1] = weight00*right_weight*Magnitude_trans[i - 1][j - 1];
			}
			else
			{
				float weight11 = tri_linear(db*(px - 1) + db / 2, db*(py - 1) + db / 2, i, j, db);
				result[px - 1][py - 1][left_side - 1] = weight11*left_weight*Magnitude_trans[i - 1][j - 1];
				result[px - 1][py - 1][right_side - 1] = weight11*right_weight*Magnitude_trans[i - 1][j - 1];
				float weight00 = tri_linear(db*(px - 2) + db / 2, db*(py - 2) + db / 2, i, j, db);
				result[px - 2][py - 2][left_side - 1] = weight11*left_weight*Magnitude_trans[i - 1][j - 1];
				result[px - 2][py - 2][right_side - 1] = weight11*right_weight*Magnitude_trans[i - 1][j - 1];
				float weight10 = tri_linear(db*(px - 1) + db / 2, db*(py - 2) + db / 2, i, j, db);
				result[px - 1][py - 2][left_side - 1] = weight11*left_weight*Magnitude_trans[i - 1][j - 1];
				result[px - 1][py - 2][right_side - 1] = weight11*right_weight*Magnitude_trans[i - 1][j - 1];
				float weight01 = tri_linear(db*(px - 2) + db / 2, db*(py - 1) + db / 2, i, j, db);
				result[px - 2][py - 1][left_side - 1] = weight11*left_weight*Magnitude_trans[i - 1][j - 1];
				result[px - 2][py - 1][right_side - 1] = weight11*right_weight*Magnitude_trans[i - 1][j - 1];
			}
		}
	}
	for (int i = 0; i < 2 * in_radius + 1; i++)
	{
		free(Magnitude_trans[i]);
		free(Orientation_trans[i]);
	}
	return result;
}


float* judgeRange(float degree) {
	float *result = (float*)malloc(sizeof(float) * 2);
	if (degree<22.5 && degree>=0)
	{
		result[0] = 8;
		result[1] = 1.f - (degree + 22.5) / 45.f;
	}
	else if (degree >= 337.5 &&degree <= 360)
	{
		result[0] = 8;
		result[1] = 1 - (degree - 337.5) / 45;
	}
	else if (degree >= 22.5 && degree < 67.5)
	{
		result[0] = 1;
		result[1] = 1 - (degree - 22.5) / 45;
	}
	else if (degree >= 67.5  && degree <112.5)
	{
		result[0] = 2;
		result[1] = 1 - (degree - 67.5) / 45;
	}
	else if (degree >= 112.5 && degree<157.5)
	{
		result[0] = 3;
		result[1] = 1 - (degree - 112.5) / 45;
	}
	else if (degree >= 157.5 && degree <202.5)
	{
		result[0] = 4;
		result[1] = 1 - (degree - 157.5) / 45;
	}
	else if (degree >= 202.5 &&degree < 247.5)
	{
		result[0] = 5;
		result[1] = 1 - (degree - 202.5) / 45;
	}
	else if (degree >= 247.5 && degree <292.5)
	{
		result[0] = 6;
		result[1] = 1 - (degree - 247.5) / 45;
	}
	else if (degree >= 292.5 && degree<337.5)
	{
		result[0] = 7;
		result[1] = 1 - (degree - 292.5) / 45;
	}
	return result;
}

float tri_linear(float center_x, float center_y, int point_x, int point_y, float distance_db)
{
	return (1.f - abs((float)(center_x - point_x) / distance_db))*(1.f - abs((float)(center_y - point_y) / distance_db));
}

void showMatchImage(Mat src1, std::vector<SFKP>siftpoints1, Mat src2, std::vector<SFKP>siftpoints2,Mat output, std::vector<MP> matchedpoints)
{
	int sum = 0;
	int op_row = src1.rows > src2.rows ? src2.rows : src1.rows;
	int op_col = src1.rows > src2.rows ? src1.cols : src2.cols;
	for (int i = 0; i < output.rows; i++)
	{
		uchar *op = output.ptr<uchar>(i);
		if (i >= op_row)
		{
			const uchar *sp = src1.rows > src2.rows ? src1.ptr<uchar>(i) : src2.ptr<uchar>(i);
			for (int j = 0; j < op_col*output.channels(); j++)
			{
				*op++ = *sp++;
			}
		}
		else
		{
			const uchar *sp1 = src1.rows > src2.rows ? src1.ptr<uchar>(i) : src2.ptr<uchar>(i);
			const uchar *sp2 = src1.rows > src2.rows ? src2.ptr<uchar>(i) : src1.ptr<uchar>(i);
			for (int j = 0; j < output.cols*output.channels(); j++)
			{
				if (j < op_col*output.channels())
				{
					*op++ = *sp1++;
				}
				else
					*op++ = *sp2++;
			}
		}
	}
	for (int i = 0; i < matchedpoints.size();)
	{
		int min_dist_index = i;
		int next_min_index = i;
		float rotate_angle = 0.f;
		int step = 1;
		if (std::abs(matchedpoints[i].rotation-rotate_angle)<=10)
		{
			Point pt1 = Point(siftpoints1[matchedpoints[min_dist_index].input_num].column, siftpoints1[matchedpoints[min_dist_index].input_num].row);
			Point pt2 = Point(siftpoints2[matchedpoints[min_dist_index].compare_num].column + src1.cols, siftpoints2[matchedpoints[min_dist_index].compare_num].row);
			circle(output, pt1, 5, Scalar((10 * i + i) % 255, i % 255, (i + 5 * i) % 255), 2, 10);
			circle(output, pt2, 5, Scalar((10 * i + i) % 255, i % 255, (i + 5 * i) % 255), 2, 10);
			line(output, pt1, pt2, Scalar((10 * i + i) % 255, i % 255, (i + 5 * i) % 255), 2);
			sum++;
		}
		//while (i<matchedpoints.size()-step && siftpoints1[matchedpoints[i].input_num].row == siftpoints1[matchedpoints[i + step].input_num].row&&
		//	siftpoints1[matchedpoints[i].input_num].column == siftpoints1[matchedpoints[i + step].input_num].column)
		//{
		//	next_min_index = matchedpoints[min_dist_index].eu_distance > matchedpoints[i + step].eu_distance ? min_dist_index : (i + step);
		//	min_dist_index = matchedpoints[min_dist_index].eu_distance > matchedpoints[i + step].eu_distance ? (i+step): min_dist_index;
		//	++step;
		//}
		//for (int j = 1;j<=step&&i+j<matchedpoints.size();++j)
		//{
		//	if (matchedpoints[next_min_index].eu_distance>matchedpoints[i+j].eu_distance&&
		//		matchedpoints[i+j].eu_distance!=matchedpoints[min_dist_index].eu_distance)
		//	{
		//		next_min_index = i + j;
		//		break;
		//	}
		//}
		//if (next_min_index == min_dist_index || matchedpoints[min_dist_index].eu_distance / matchedpoints[next_min_index].eu_distance <= 0.8)
		//{
		//	Point pt1 = Point(siftpoints1[matchedpoints[min_dist_index].input_num].column, siftpoints1[matchedpoints[min_dist_index].input_num].row);
		//	Point pt2 = Point(siftpoints2[matchedpoints[min_dist_index].compare_num].column + src1.cols, siftpoints2[matchedpoints[min_dist_index].compare_num].row);
		//	circle(output, pt1, 3, Scalar(10 * i % 255, i % 255, 5 * i % 255), 2, 10);
		//	circle(output, pt2, 3, Scalar(10 * i % 255, i % 255, 5 * i % 255), 2, 10);
		//	line(output, pt1, pt2, Scalar(10*i%255, i%255, 5*i%255),2);
		//	sum++;
		//}
		i += step;
	}
	printf("%d,%d",matchedpoints.size(),sum);
	imshow("siftÍ¼ÏñÆ¥Åä", output);
}

void buildSiftPrograme(Mat src, std::vector<SFKP> &siftpoints, float diff_threshold,float sigma)
{
	if (src.channels() == 3 || src.channels() == 4)
	{
		cvtColor(src, src, COLOR_RGB2GRAY);
	}
	int octave_num = cvRound(log(src.rows < src.cols ? src.rows : src.cols) / log(2) - 2)+1;
	//octave_num = octave_num > 4 ? 4 : octave_num;
	Mat ***GaussianPyramid;
	GaussianPyramid = (Mat***)malloc(sizeof(Mat**)*octave_num);
	for (int i = 0; i < octave_num; i++)
	{
		GaussianPyramid[i] = (Mat**)malloc(sizeof(Mat_<float>*)*(LAYER_NUM + 3));
		for (int j = 0; j < LAYER_NUM + 3; j++)
		{
			GaussianPyramid[i][j] = new Mat_<float>(i>0 ? src.rows >> (i - 1) : src.rows * 2, i>0 ? src.cols >> (i - 1) : src.cols * 2);
		}
	}
	buildGaussianPyramid(src, GaussianPyramid, src.cols, src.rows,octave_num,sigma);
	Mat ***DogPyramid = (Mat***)malloc(sizeof(Mat**)*octave_num);
	for (int i = 0; i < octave_num; i++)
	{
		DogPyramid[i] = (Mat**)malloc(sizeof(Mat_<float>*)*(LAYER_NUM + 2));
		for (int j = 0; j < LAYER_NUM + 2; j++)
		{
			DogPyramid[i][j] = new Mat_<float>(i>0?src.rows >> (i-1):src.rows*2, i>0?src.cols >> (i-1):src.cols*2);
		}
	}
	buildDogPyramid(GaussianPyramid, DogPyramid,octave_num);
	//for (int j = 0; j<src.cols; ++j)
	//{
	//	printf("%f,", DogPyramid[0][1]->at<float>(src.rows - 1, j));
	//}
	findExtemePoint(GaussianPyramid,DogPyramid,siftpoints,octave_num, diff_threshold);
	//findOrientation(GaussianPyramid,siftpoints,src.rows, src.cols);

	for (int i = 0; i < octave_num; i++)
	{
		for (int j = 0; j < LAYER_NUM + 2; j++)
		{
			if (j == LAYER_NUM + 1)
			{
				free(GaussianPyramid[i][LAYER_NUM+2]);
			}
			free(DogPyramid[i][j]);
		}
		free(DogPyramid[i]);
		free(GaussianPyramid[i]);
	}
	free(DogPyramid);
	free(GaussianPyramid);
}

void showSiftPoints(std::string title,Mat src, std::vector<SFKP>& siftpoints, int radius, Scalar sc, int linetype)
{
	std::ofstream f1("H:\\MatlabCodes\\siftDemoV4\\points.txt");
	for (int i = 0; i < siftpoints.size(); i++)
	{
		float x = siftpoints[i].row;
		x += siftpoints[i].offset[0] * (float)(siftpoints[i].octave_num>0?1<<(siftpoints[i].octave_num-1):0.5);
		float y = siftpoints[i].column;
		y += siftpoints[i].offset[1] * (float)(siftpoints[i].octave_num>0?1 << (siftpoints[i].octave_num-1):0.5);
		Point ptstart(y,x);
		//if ((int)y==7&&(int)x==6)
		//{
		//	//printf("%f", siftpoints[i].orientation);
		//	for (int m = 0; m < 4; m++)
		//	{
		//		for (int n = 0; n < 4; n++)
		//		{
		//			for (int k = 0; k < 8; k++)
		//			{
		//				printf("%f,", siftpoints[i].descripter[m][n][k]);
		//			}
		//		}
		//	}
		//	printf("\n");
		//}
		circle(src,ptstart,radius, sc,1,linetype);
		int r = siftpoints[i].layer_num + 1;
		float s_radius = mat_size[r] << siftpoints[i].octave_num;//siftpoints[i].octave_num>0?mat_size[r]<<(siftpoints[i].octave_num-1):mat_size[r]>>1;
		float dy = cos(siftpoints[i].orientation * Pi / 180.f );
		float dx = sin(siftpoints[i].orientation  * Pi / 180.f);
		//float dx = sin(siftpoints[i].orientation * 10 < 180 ? siftpoints[i].orientation * 10 * Pi / 180.f : siftpoints[i].orientation * 10 * Pi / 180.f - Pi);
		Point ptend(y + s_radius*dx, x + s_radius*dy);
		f1 <<y << "," << x <<","<< siftpoints[i].column + radius*cos(dy)<<","<< siftpoints[i].row + radius*sin(dx)<<std::endl;
		line(src, ptstart, ptend, sc);
	}
	f1.close();
	printf("%d", siftpoints.size());
	imshow(title, src);
}

void matchSiftPoints(std::vector<SFKP>& siftpoints1, std::vector<SFKP>& siftpoints2, std::vector<MP>& matchPoints,float dis_threshold)
{
	for (int o = 0; o<siftpoints1.size(); o++)
	{
		float min_dist = 0.8f;
		float sec_dist = 1.f;
		int cp_pos = -1;
		for (int r = 0; r<siftpoints2.size(); r++)
		{
			float t_distance = 0.f;
			for (int i = 0; i < 4; i++)
			{
				for (int j = 0; j < 4; j++)
				{
					for (int k = 0; k < 8; k++)
					{
						float p_distance = siftpoints1[o].descripter[i][j][k] - siftpoints2[r].descripter[i][j][k];
						t_distance += p_distance*p_distance;
						if (t_distance> dis_threshold)
							break;
					}
					if (t_distance> dis_threshold)
						break;
				}
				if (t_distance> dis_threshold)
					break;
			}
			//t_distance = sqrt(t_distance);
			if (t_distance < dis_threshold)
			{
				sec_dist = t_distance;
				if (min_dist>t_distance)
				{
					sec_dist = min_dist;
					min_dist = t_distance;
					cp_pos = r;
				}
			}
		}
		if (cp_pos>-1)
		{
			if (min_dist/sec_dist<=0.8 || min_dist == sec_dist)
			{
				MP *mp = (MP*)malloc(sizeof(MP));
				mp->compare_num = cp_pos;
				mp->input_num = o;
				mp->eu_distance = min_dist;
				mp->rotation = siftpoints1[o].orientation - siftpoints2[cp_pos].orientation;
				matchPoints.push_back(*mp);
				rotate_angle += mp->rotation;
			}
		}
	}
	rotate_angle /= matchPoints.size();
}


