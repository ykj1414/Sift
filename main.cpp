#include "Sift.h"
#include <iostream>
#include <opencv.hpp>
#include <vector>
#include <stdio.h>
#include <malloc.h>
using std::vector;
using std::cout;
using std::endl;
using std::distance;
using namespace cv;

int main() {
	//获取运行时间
	double t = (double)getTickCount();
	Mat	src1 = imread("IMG_5031.JPG");
	Mat src2 = imread("IMG_5030.JPG");
	if (src1.empty()||src2.empty())
	{
		printf("文件不存在或格式错误!!!");
		return -1;
	}
	std::vector<SFKP>siftpoints1;
	std::vector<SFKP>siftpoints2;
	std::vector<MP>mpoints;
	buildSiftPrograme(src1, siftpoints1, 0.03, 1.6);
	//showSiftPoints("scene", src1, siftpoints1,2);
	buildSiftPrograme(src2, siftpoints2, 0.03,1.6);
	//showSiftPoints("box",src2, siftpoints2, 3);
	Mat output = Mat(src1.rows >= src2.rows ? src1.rows : src2.rows, src1.cols + src2.cols, CV_8UC3);
	matchSiftPoints(siftpoints1, siftpoints2, mpoints, 0.03);
	showMatchImage(src1, siftpoints1, src2, siftpoints2, output, mpoints);
	t = ((double)getTickCount() - t) / cvGetTickFrequency() * 1e-6;  //单位为s
	cout << "the time is :" << t << endl;
	waitKey(0);
}