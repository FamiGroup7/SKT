#pragma once

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <math.h>
#include <iomanip>

using namespace std;

class Point
{
public:
	Point();
	Point(double x, double y, double z);
	~Point();
	double x, y, z;
	friend ostream& operator<<(ostream& os, const Point& point);
	Point operator- (Point& left);
private:

};
