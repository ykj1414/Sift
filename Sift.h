#pragma once
#include <opencv.hpp>
#include <vector>

#define Pi		3.1415926

#define SIGMA_INITIAL	0.5

#define EDGE_THRESHOLD 10

#define  DIFF_THRESHOLD 0.04

#define LAYER_NUM		3
using namespace cv;
using std::vector;
typedef struct siftkeypoints {
	unsigned char octave_num;
	unsigned char layer_num;
	float *offset;
	int row;
	int column;
	float magnitude;
	float orientation;
	float ***descripter;
}SFKP;

typedef struct matchpoints {
	int input_num;
	int compare_num;
	float eu_distance;
	float rotation;
}MP;

void buildSiftPrograme(Mat src, std::vector<SFKP>&siftpoints,float diff_threshold = DIFF_THRESHOLD,float sigma = 1.6);

void showSiftPoints(std::string title,Mat src, std::vector<SFKP>&siftpoints,int radius = 1,Scalar sc = Scalar(0,0,255),int linetype = 8);

void matchSiftPoints(std::vector<SFKP>&sifpoints1, std::vector<SFKP>&sifpoints2,std::vector<MP>&matchPoints,float dis_threshold = 0.01);

void showMatchImage(Mat src1, std::vector<SFKP>siftpoints1, Mat src2, std::vector<SFKP>siftpoints2, Mat output, std::vector<MP> matchedpoints);