/*
 * espPoints.h
 *
 *  Created on: Jun 25, 2018
 *      Author: cmyers
 */

#ifndef ESPPOINTS_H_
#define ESPPOINTS_H_

#define _USE_MATH_DEFINES

#include <vector>
#include <string>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <algorithm>

#include "atomicRadii.h"
#include "customStructs.hpp"



class espPoints {

	double maxX, maxY, maxZ, minX, minY, minZ;
	double head, gridSpace;
	double density;
	vecD intervals;
	int xPts, yPts, zPts;
	vecD radii;
	vecD xEsp, yEsp, zEsp;

	void formBoxDims(const std::vector<double>&, const std::vector<double>&,
			const std::vector<double>&);
	void assignRadii(const std::vector<int>&);

	void printESPFile(const std::string outFileLoc, double scale = 1.0);
public:
	espPoints();
	virtual ~espPoints();


	void genPoints(const std::vector<int>&, const vecD& xCoord,
			const vecD& yCoord, const vecD& zCoord, const std::string);

	void genPoints_MK(const std::vector<int>&, const vecD& xCoord,
			const vecD& yCoord, const vecD& zCoord, const std::string);
	
	void printPoints(std::string);

};

#endif /* ESPPOINTS_H_ */
