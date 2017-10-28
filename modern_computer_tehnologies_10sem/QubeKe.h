#pragma once
#include "Point.h"

class QubeKe
{
public:
	QubeKe();
	QubeKe(Point node1, Point node2);
	QubeKe(double xLeft,double yLeft,double zLeft,double xRight,double yRight,double zRight);
	~QubeKe();
	Point calcCenter();
	double calcVolume();
	void initQube(double xLeft, double yLeft, double zLeft, double xRight, double yRight, double zRight);
	friend ostream& operator<<(ostream& os, const QubeKe& qube);
	
	int number;
	Point *nodes;
private:
	static int count;
};