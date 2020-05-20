#pragma once

#include <cassert>
#include <vector>
#include <algorithm>
#include <Eigen/Dense>
#include <iostream>
#include <iomanip>

class GeometryItem
{
public:
	GeometryItem();
	GeometryItem(const int index);
	~GeometryItem();

	void setIndex(const int index);
	const int getIndex() const;

private:
	int index = -1;
};

