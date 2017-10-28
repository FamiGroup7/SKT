#pragma once
#include "QubeKe.h"

int QubeKe::count = 0;

Point QubeKe::calcCenter() {
	double xCenter, yCenter, zCenter;
	xCenter = nodes[1].x + nodes[0].x;
	yCenter = nodes[2].y + nodes[0].y;
	zCenter = nodes[4].z + nodes[0].z;
	return Point(xCenter / 2, yCenter / 2, zCenter / 2);
}

double QubeKe::calcVolume()
{
	double hx, hy, hz;
	hx = fabs(nodes[1].x - nodes[0].x);
	if (hx == 0)hx = 1;
	hy = fabs(nodes[2].y - nodes[0].y);
	if (hy == 0)hy = 1;
	hz = fabs(nodes[4].z - nodes[0].z);
	if (hz == 0)hz = 1;
	return hx*hy*hz;
}

void QubeKe::initQube(double xLeft, double yLeft, double zLeft, double xRight, double yRight, double zRight)
{
	nodes[0] = Point(xLeft, yLeft, zLeft);
	nodes[1] = Point(xRight, yLeft, zLeft);
	nodes[2] = Point(xLeft, yRight, zLeft);
	nodes[3] = Point(xRight, yRight, zLeft);

	nodes[4] = Point(xLeft, yLeft, zRight);
	nodes[5] = Point(xRight, yLeft, zRight);
	nodes[6] = Point(xLeft, yRight, zRight);
	nodes[7] = Point(xRight, yRight, zRight);
}

QubeKe::QubeKe()
{
	number = count++;
	nodes = new Point[8];
}

QubeKe::QubeKe(Point node1, Point node2)
{
	number = count++;
	nodes = new Point[8];
	double xLeft = node1.x;
	double yLeft = node1.y;
	double zLeft = node1.z;
	double xRight = node2.x;
	double yRight = node2.y;
	double zRight = node2.z;
	initQube(xLeft, yLeft, zLeft, xRight, yRight, zRight);
}

QubeKe::QubeKe(double xLeft, double yLeft, double zLeft, double xRight, double yRight, double zRight)
{
	number = count++;
	nodes = new Point[8];
	initQube(xLeft, yLeft, zLeft, xRight, yRight, zRight);
}

QubeKe::~QubeKe()
{
	delete[]nodes;
}

ostream & operator<<(ostream & os, const QubeKe & qube)
{
	os << "{" << endl;
	os << "\t" << "number = " << qube.number << endl;
	for (size_t i = 0; i < 8; i++)
	{
		os << "\t" << qube.nodes[i] << endl;
	}
	os << "}" << endl;
	return os;
}