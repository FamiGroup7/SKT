#include "Point.h"

Point::Point()
{
	x = y = z = 0;
}

Point::Point(double newX, double newY, double newZ)
{
	x = newX;
	y = newY;
	z = newZ;
}

Point::~Point()
{
}

Point Point::operator-(Point & left)
{
	return Point(left.x - this->x, left.y - this->y, left.z - this->z);
}

ostream & operator<<(ostream & os, const Point & point)
{
	os << "{ x = " << point.x << "; y = " << point.y << "; z = " << point.z << "}";
	return os;
}
