/*
 * espPoints.cpp
 *
 *  Created on: Jun 25, 2018
 *      Author: cmyers
 */

#include "espPoints.h"

using namespace std;

double atomicRadii::bohrRadius_ang = 0.529177249;

espPoints::espPoints() {
	// TODO Auto-generated constructor stub

	maxX = -DBL_MAX; maxY = -DBL_MAX; maxZ = -DBL_MAX;
	minX = DBL_MAX; minY = DBL_MAX; minZ = DBL_MAX;

	// ChElPG point parameters
	xPts = 0; yPts = 0; zPts = 0;
	head = 4.0;
	gridSpace = 0.3;

	// Merz-Kollman parameters
	//intervals = {1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0};
	intervals = {1.4, 1.6, 1.8, 2.0};
	//intervals = {1.2};
	density = 7.0;
}

espPoints::~espPoints() {
	// TODO Auto-generated destructor stub
}

void espPoints::formBoxDims(const vector<double>& xCoord,
		const vector<double>& yCoord, const vector<double>& zCoord)
{
	//double bohrToAng = 0.529177249;

	//find maximum and minimum positions + head-room
	for(int i = 0; i < (signed)xCoord.size(); i++)
	{
		maxX = max(xCoord[i] + head, maxX);
		maxY = max(yCoord[i] + head, maxY);
		maxZ = max(zCoord[i] + head, maxZ);

		minX = min(xCoord[i] - head, minX);
		minY = min(yCoord[i] - head, minY);
		minZ = min(zCoord[i] - head, minZ);
	}

	//number of grid points in each direction
	xPts = (int)ceil((maxX - minX)/gridSpace) + 1;
	yPts = (int)ceil((maxY - minY)/gridSpace) + 1;
	zPts = (int)ceil((maxZ - minZ)/gridSpace) + 1;

	//print out results to console
	printf(" Number of steps:    [ %4d %4d %4d ]\n", xPts, yPts, zPts);
	printf(" Box coordinates:    [ %6.2f  %6.2f ]\n", minX, maxX);
	printf("                     [ %6.2f  %6.2f ]\n", minY, maxY);
	printf("                     [ %6.2f  %6.2f ]\n", minZ, maxZ);
	printf(" Max num. of points:     %8d\n", xPts*yPts*zPts);
	printf(" Total volume:           %8.2f\n", (maxX - minX)*(maxY - minY)*(maxZ - minZ));
}

void espPoints::assignRadii(const vector<int> &atomIdx)
{
	for (int n = 0; n < (signed)atomIdx.size(); n++)
	{
		switch(atomIdx[n])
		{
			case 1:
				radii.push_back(1.45);
				break;
			case 6:
				radii.push_back(1.5);
				break;
			case 7:
				radii.push_back(2.0);
				break;
			case 8:
				radii.push_back(2.0);
				break;
			case 15:
				radii.push_back(2.5);
				break;
			case 16:
				radii.push_back(2.5);
				break;
		}
	}
}

void espPoints::printESPFile(const std::string outFileLoc, double scale)
{
	//print grid points to file
	ofstream outEsp;
	outEsp.open(outFileLoc.c_str(), std::ofstream::binary);
	outEsp.precision(8);
	outEsp << scientific;
	for(int i = 0; i < (signed)xEsp.size(); i++)
	{
		outEsp << setw(15) << xEsp[i] * scale << "  ";
		outEsp << setw(15) << yEsp[i] * scale << "  ";
		outEsp << setw(15) << zEsp[i] * scale << "\n";
	}
	outEsp.close();
}

void espPoints::genPoints_MK(const vector<int>& atomType, const vector<double>& xCoord,
		const vector<double>& yCoord, const vector<double>& zCoord, const string outFileLoc)
{
	double xDiff, yDiff, zDiff, norm;
	double xPoint = 0.0, yPoint = 0.0, zPoint = 0.0;
	double scale = 1.0, dTheta, dPhi, pRad, radius;
	double theta, phi, sinT, cosT;
	const double angToBohr = 1.8897259886;
	int numTheta, numPhi;
	int idx;
	int n_atoms = (int)xCoord.size();
	bool accept;
	vector<int> nbrList;

	xEsp.reserve(5000*xCoord.size());
	yEsp.reserve(5000*xCoord.size());
	zEsp.reserve(5000*xCoord.size());
	
	printf("Point generation will use an arc density of\n %.5f between points.\n", density);

	/*  Loop over vdW radii intervals */
	for(int i = 0; i < (int)intervals.size(); i++)
	//for(int i = 0; i < 1; i++)
	{
		scale = intervals[i];
		printf(" Generating points for vdW radii interval %.4f\n", scale);

		/*  Loop over all atoms */
		for(int n = 0; n < n_atoms; n++)
		{
			radius = atomicRadii::Bondi_radii(atomType[n])*scale;
			numTheta = floor(radius * M_PI * density);
			dTheta = M_PI/numTheta;

			nbrList.clear();

			/*  Search for possible intersecting neighbors */
			for(int m = 0; m < n_atoms; m++)
			{
				xDiff = xCoord[n] - xCoord[m];
				yDiff = yCoord[n] - yCoord[m];
				zDiff = zCoord[n] - zCoord[m];

				norm = sqrt(xDiff*xDiff + yDiff*yDiff + zDiff*zDiff);
				if(norm <= (radius + atomicRadii::Bondi_radii(atomType[n])) * scale && m != n)
					nbrList.push_back(m);
			}

			/*  Loop over theta angles */
			for(int t = 0; t < numTheta; t++)
			{
				theta = t * dTheta;
				sinT = sin(theta);
				cosT = cos(theta);

				pRad = radius * sinT;
				numPhi = floor(pRad * 2 * M_PI * density);
				dPhi = 2 * M_PI / numPhi;

				/*  Loop over phi angles */
				for(int p = 0; p < numPhi; p++)
				{
					phi = p * dPhi;

					xPoint = xCoord[n] + radius * (sinT * cos(phi));
					yPoint = yCoord[n] + radius * (sinT * sin(phi));
					zPoint = zCoord[n] + radius * cosT;

					/*  Loop over neighborlist to check for intersection */
					accept = true;
					for(int s = 0; s < (int)nbrList.size(); s++)
					{
						idx = nbrList[s];
						xDiff = xPoint - xCoord[idx];
						yDiff = yPoint - yCoord[idx];
						zDiff = zPoint - zCoord[idx];

						norm = sqrt(xDiff*xDiff + yDiff*yDiff + zDiff*zDiff);

						if(norm <= atomicRadii::Bondi_radii(atomType[idx])* scale)
						{
							accept = false;
							break;
						}
					}
					if(accept)
					{
						xEsp.push_back(xPoint);
						yEsp.push_back(yPoint);
						zEsp.push_back(zPoint);
					}
				}
			}
		}
	}
	printf(" Number of point used:  %10i\n", (int)xEsp.size());
	printf(" Number of points/atom: %10.0f\n", (double)xEsp.size() / ((double) n_atoms));
	printf(" Printing points to file...");
	printESPFile(outFileLoc, angToBohr);
	printf("done\n");
}

void espPoints::genPoints(const vector<int>& atomType, const vector<double>& xCoord,
		const vector<double>& yCoord, const vector<double>& zCoord, const string outFileLoc)
{
	formBoxDims(xCoord, yCoord, zCoord);
	assignRadii(atomType);

	int xIdx = 0, yIdx = 0, zIdx = 0, rem = 0, i = 0;
	const int lenCoords = xCoord.size();
	const int yzPts = yPts * zPts;
	const double angToBohr = 1.8897259886;

	double xDiff, yDiff, zDiff, norm;
	double xPoint = 0.0, yPoint = 0.0, zPoint = 0.0;

	bool accept1, accept2;
	//double minDist = std::numeric_limits<double>::max();

	//loop through all possible grid points
	for(int n = 0; n < xPts*yPts*zPts; n++)
	{
		xIdx = n / yzPts;
		rem = n % yzPts;
		yIdx = rem / zPts;
		zIdx = rem % zPts;

		xPoint = minX + xIdx*gridSpace;
		yPoint = minY + yIdx*gridSpace;
		zPoint = minZ + zIdx*gridSpace;
		accept1 = true;
		accept2 = false;
		//minDist = std::numeric_limits<double>::max();
		for(i = 0; i < lenCoords; i++)
		{
			xDiff = xCoord[i] - xPoint;
			yDiff = yCoord[i] - yPoint;
			zDiff = zCoord[i] - zPoint;

			norm = sqrt(xDiff*xDiff + yDiff*yDiff + zDiff*zDiff);

			//accept grid point if it is within vdW radii and the head-space
			//if((norm > radii[i]) && (norm < head))
			//if((norm < atomicRadii::bohrRadius_ang) || (norm > head))
			if((norm < 1.2*atomicRadii::Bondi_radii(atomType[i])))
			{
				accept1 = false;
				break;
			}
			if (norm < head)
				accept2 = true;
		}
		if (accept1 && accept2)
		{
			xEsp.push_back(xPoint);
			yEsp.push_back(yPoint);
			zEsp.push_back(zPoint);
		}
	}

	printf(" Number of point used: %8i\n", (int)xEsp.size());

	//print grid points to file
	printESPFile(outFileLoc, angToBohr);
}






