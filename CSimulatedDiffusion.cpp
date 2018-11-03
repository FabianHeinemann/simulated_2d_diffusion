#include <new>
#include <fstream>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/variate_generator.hpp>
#include <math.h>
#include "CSimulatedDiffusion.h"

int inline floor_int(double x);

CSimulatedDiffusion::CSimulatedDiffusion()
{
	nParticles		= 0;
	x				= 0;
	y				= 0;
	edges			= 0;
	nEdges			= 0;
	boxWidth		= 10*1e-6; 
	boxHeight		= boxWidth;
	n				= 0;
	nUpdate			= 100;
	deltaHashmap	= boxWidth;
	vector<Edge*>	edgesInABox;
	edgesInABox.clear();
	spatialHashmap.push_back(edgesInABox);
	Tmax			= 0.5;
	pjump			= 0.1;

	xTraj			= 0;
	yTraj			= 0;

	spotx1			= boxWidth/2 - boxWidth/10;
	spotx2			= boxWidth/2;
	spotx3			= boxWidth/2 + boxWidth/10;
	spoty1			= boxHeight/2 - boxHeight/10;
	spoty2			= boxHeight/2;
	spoty3			= boxHeight/2 + boxWidth/10;
}


CSimulatedDiffusion::~CSimulatedDiffusion()
{
	if (x)
	{
		delete[] x;
		x = 0;
	}
	if (y)
	{
		delete[] y;
		y = 0;
	}
	if (edges)
	{
		delete[] edges;
		edges = 0;
	}
	if (xTraj)
	{
		delete[] xTraj;
	}
	if (yTraj)
	{
		delete[] yTraj;
	}
}

void CSimulatedDiffusion::start()
{
	unsigned int		m, l;
	double				dx, dy;
	float				sigma;

	// Detection spots
	double*				F1;
	double*				F2;
	double*				F3;
	double*				F4;
	double*				F5;
	double*				F6;
	double*				F7;
	double*				F8;
	double*				F9;

	vector<Collision>	collisionList;
	bool				reflected;

	// Set simulation parameters	
	const double			r0	= 250e-9;
	const double			r02	= 3.125e-14;

	nParticles			=	1000;
	D					=	10e-12;
	dt					=	10*1e-6;
	#ifndef SPEED_TEST
	nSteps				=	floor_int(Tmax / dt);	
	#else
	nSteps				=	100000;
	#endif
	boxWidth			=	10*1e-6; 
	boxHeight			=	10*1e-6; 
	sigma				=	float(sqrt(2*D*dt));

	// Reserve memory for trajectory
	xTraj = new double[floor_int(nSteps/nUpdate)];
	yTraj = new double[floor_int(nSteps/nUpdate)];

	// Reserve memory for fluorescence traces	
	F1 = new double[nSteps];	
	F2 = new double[nSteps];
	F3 = new double[nSteps];
	F4 = new double[nSteps];
	F5 = new double[nSteps];
	F6 = new double[nSteps];
	F7 = new double[nSteps];
	F8 = new double[nSteps];
	F9 = new double[nSteps];
	
	for (n = 0; n < nSteps; n++)
	{
		F1[n] = 0;
		F2[n] = 0;
		F3[n] = 0;
		F4[n] = 0;
		F5[n] = 0;
		F6[n] = 0;
		F7[n] = 0;
		F8[n] = 0;
		F9[n] = 0;
	}

	// Random number generator seed
	// Mersenne twister algorithm
	boost::mt11213b mtGenerator(unsigned int (time(0)));	
	boost::uniform_01<boost::mt11213b> getRandomNumber01(mtGenerator);
	
	// Initialize particle positions
	x = new double[nParticles];
	y = new double[nParticles];

	for (n = 0; n < nParticles; n++)
	{
		x[n] = getRandomNumber01()*boxWidth;		
		y[n] = getRandomNumber01()*boxHeight;
		//x[n] = boxWidth/2;
		//y[n] = boxHeight/2;
	}
	// Uncomment next two lines for real simulation, only for trajectory starting from center
	x[0] = boxWidth/2.0;
	y[0] = boxHeight/2.0;

	#ifdef GUI_VERSION
	frame->panel->Refresh();
	#endif
	
	// Normally distributed random numbers for simulation
	boost::normal_distribution<double> distribution(0.0,sigma);
	boost::variate_generator<boost::mt11213b&, boost::normal_distribution<double>> getPos(mtGenerator, distribution); 

	#ifdef SPEED_TEST
	clock_t startTime = clock();
	#endif

	// Simulation loop
	for (n = 0; n < nSteps; n++)
	{
		// Loop over all particles
		for (m = 0; m < nParticles; m++)
		{
			// Update vector for positions
			dx = getPos();			
			dy = getPos();

			// Check for collisions
			while(isColliding(x[m], y[m], dx, dy, &collisionList))
			{
				reflected = false;
				// The collisions are sorted according to distance
				// Check them one by one
				for (l = 0; l < collisionList.size(); l++)
				{
					// Try to jump
					if (getRandomNumber01() >= pjump)
					{
						// Calculate refection 
						// Variables x[m], y[m], dx, dy are updated by the function
						reflect(x[m], y[m], dx, dy, &collisionList[l]);
						reflected = true;
						break;
					}
				}	

				// All possible reflections failed
				// Leave collision loop and perform the step
				if (reflected == false)
				{
					break;
				}
			}

			// Update position
			x[m] = x[m] + dx;
			y[m] = y[m] + dy;			

			// Periodic bounray conditions (customized modulo operation)
			if (x[m] < 0 || x[m] > boxWidth)
			{
				x[m] = x[m] - boxWidth * floor_int(x[m] / boxWidth);
			}
			if (y[m] < 0 || y[m] > boxHeight)
			{
				y[m] = y[m] - boxHeight * floor_int(y[m] / boxHeight);
			}

			// Acquire fluorescence signals
			F1[n] += exp(-((x[m]-spotx1)*(x[m]-spotx1)+(y[m]-spoty1)*(y[m]-spoty1))/r02);			
			F2[n] += exp(-((x[m]-spotx2)*(x[m]-spotx2)+(y[m]-spoty1)*(y[m]-spoty1))/r02);
			F3[n] += exp(-((x[m]-spotx3)*(x[m]-spotx3)+(y[m]-spoty1)*(y[m]-spoty1))/r02);
			F4[n] += exp(-((x[m]-spotx1)*(x[m]-spotx1)+(y[m]-spoty2)*(y[m]-spoty2))/r02);
			F5[n] += exp(-((x[m]-spotx2)*(x[m]-spotx2)+(y[m]-spoty2)*(y[m]-spoty2))/r02);
			F6[n] += exp(-((x[m]-spotx3)*(x[m]-spotx3)+(y[m]-spoty2)*(y[m]-spoty2))/r02);
			F7[n] += exp(-((x[m]-spotx1)*(x[m]-spotx1)+(y[m]-spoty3)*(y[m]-spoty3))/r02);
			F8[n] += exp(-((x[m]-spotx2)*(x[m]-spotx2)+(y[m]-spoty3)*(y[m]-spoty3))/r02);
			F9[n] += exp(-((x[m]-spotx3)*(x[m]-spotx3)+(y[m]-spoty3)*(y[m]-spoty3))/r02);
		}
			
		#ifdef CONSOLE_VERSION
		if (n % nUpdate == 0)
		{			
			cout << (n*100.0/nSteps) << "% (" << dt*n << " s of " << Tmax << " s).      " << "\r";
			cout << flush;
		}
		#endif

	}	

	// Calculate autocorrelation
	vector<double> G1, G2, G3, G4, G5, G6, G7, G8, G9, Gav, tau;
	vector<vector<double>*> Glist;
	vector<string> header;

	#ifdef CONSOLE_VERSION
	cout << "Correlating.                    ";
	#endif
	multipleTauCorrelate(F1, nSteps, dt, G1, tau);
	Glist.push_back(&G1);
	header.push_back("G1");
	multipleTauCorrelate(F2, nSteps, dt, G2, tau);
	Glist.push_back(&G2);
	header.push_back("G2");
	multipleTauCorrelate(F3, nSteps, dt, G3, tau);
	Glist.push_back(&G3);
	header.push_back("G3");
	multipleTauCorrelate(F4, nSteps, dt, G4, tau);
	Glist.push_back(&G4);
	header.push_back("G4");
	multipleTauCorrelate(F5, nSteps, dt, G5, tau);
	Glist.push_back(&G5);
	header.push_back("G5");
	multipleTauCorrelate(F6, nSteps, dt, G6, tau);
	Glist.push_back(&G6);
	header.push_back("G6");
	multipleTauCorrelate(F7, nSteps, dt, G7, tau);
	Glist.push_back(&G7);
	header.push_back("G7");
	multipleTauCorrelate(F8, nSteps, dt, G8, tau);
	Glist.push_back(&G8);
	header.push_back("G8");
	multipleTauCorrelate(F9, nSteps, dt, G9, tau);
	Glist.push_back(&G9);
	header.push_back("G9");	

	// Calculate average autocorrelation
	calculateAverage(&Glist, &Gav);
	#ifdef CONSOLE_VERSION
	cout << " put to List ";
	#endif
	Glist.push_back(&Gav);
	header.push_back("Gav");

	#ifdef CONSOLE_VERSION
	cout << "Saving " << outputFileStr << endl;
	#endif
	saveAutocorrelations(outputFileStr, &Glist, &tau, &header);

	// Clean up	
	delete[] F1;
	delete[] F2;
	delete[] F3;
	delete[] F4;
	delete[] F5;
	delete[] F6;
	delete[] F7;
	delete[] F8;
	delete[] F9;
}

void CSimulatedDiffusion::buildHashmap()
{
	unsigned int n, m;
	unsigned int nx = int(boxWidth/deltaHashmap);
	unsigned int ny = int(boxHeight/deltaHashmap);
	double x1, y1, x2, y2;
	vector<Edge*> edgesInABox;

	// Reserve space in spatialHashmap
	spatialHashmap.clear();
	spatialHashmap.reserve(nx*ny);

	// Loop over all boxes of the hashmap and fill them with all contained / intersecting edges	
	for (n = 0; n < nx*ny; n++)
	{
		// Generate coordinate of box
		x1 = n%nx*deltaHashmap;
		y1 = n/ny*deltaHashmap;
		x2 = (n%nx+1)*deltaHashmap;
		y2 = (n/ny+1)*deltaHashmap;

		// Test all edges, if they fall in the box		
		edgesInABox.clear();
		for (m = 0; m < nEdges; m++)
		{		
			if (edgeInBox(&edges[m], x1, y1, x2, y2))
			{
				edgesInABox.push_back(&edges[m]);
			}
		}
		// Put edge in spatial hashmap
		spatialHashmap.push_back(edgesInABox);		
	}
}

bool CSimulatedDiffusion::edgeInBox(Edge* edge, double x1, double y1, double x2, double y2)
{
	bool returnValue = false;
	Edge boxEdge1,boxEdge2,boxEdge3,boxEdge4;	

	// Check if edge is inside box
	if (x1 < edge->x1 && x2 > edge->x1 && x1 < edge->x2 && x2 > edge->x2 &&
		y1 < edge->y1 && y2 > edge->y1 && y1 < edge->y2 && y2 > edge->y2)
	{
		returnValue = true;
	}

	// Check intersection with all four edges
	// Top
	boxEdge1.x1 = x1;
	boxEdge1.y1 = y1;
	boxEdge1.x2 = x2;
	boxEdge1.y2 = y1;
	if (edgesIntersect(edge, &boxEdge1))
	{
		returnValue = true;
	}

	// Right
	boxEdge2.x1 = x2;
	boxEdge2.y1 = y1;
	boxEdge2.x2 = x2;
	boxEdge2.y2 = y2;
	if (edgesIntersect(edge, &boxEdge2))
	{
		returnValue = true;
	}
	// Bottom
	boxEdge3.x1 = x2;
	boxEdge3.y1 = y2;
	boxEdge3.x2 = x1;
	boxEdge3.y2 = y2;
	if (edgesIntersect(edge, &boxEdge3))
	{
		returnValue = true;
	}

	// Left
	boxEdge4.x1 = x1;
	boxEdge4.y1 = y2;
	boxEdge4.x2 = x1;
	boxEdge4.y2 = y1;
	if (edgesIntersect(edge, &boxEdge4))
	{
		returnValue = true;
	}

	return returnValue;
}

bool CSimulatedDiffusion::edgesIntersect(Edge* edge1, Edge* edge2)
{
	// http://paulbourke.net/geometry/lineline2d/

	double ua, ub;
	double denominator;

	denominator = ((edge2->y2-edge2->y1)*(edge1->x2-edge1->x1) - (edge2->x2-edge2->x1)*(edge1->y2-edge1->y1));

	// Parallel lines
	if (denominator == 0)
	{		
		if ((edge2->x2-edge2->x1)*(edge1->y1-edge2->y1) - (edge2->y2-edge2->y1)*(edge1->x1-edge2->x1) == 0 &&
			(edge2->x2-edge1->x1)*(edge1->y1-edge2->y1) - (edge1->y2-edge1->y1)*(edge1->x1-edge2->x1) == 0)
		{
			// Coincident edges
			return true;
		}
		return false;
	}
	else
	{
		ua = ((edge2->x2-edge2->x1)*(edge1->y1-edge2->y1) - (edge2->y2-edge2->y1)*(edge1->x1-edge2->x1)) / denominator;
		ub = ((edge1->x2-edge1->x1)*(edge1->y1-edge2->y1) - (edge1->y2-edge1->y1)*(edge1->x1-edge2->x1)) / denominator;

		if (0 <= ua && ua <= 1 && 0 <= ub && ub <= 1)
		{
			return true;
		}
	}

	return false;
}

bool CSimulatedDiffusion::edgesIntersect(double x1, double y1, double x2, double y2, Edge* edge2, double& xIntersect, double& yIntersect)
{
	// http://paulbourke.net/geometry/lineline2d/

	double ua, ub;
	double denominator;
	
	denominator = ((edge2->y2-edge2->y1)*(x2-x1) - (edge2->x2-edge2->x1)*(y2-y1));

	// Parallel lines
	if (denominator == 0)
	{		
		if ((edge2->x2-edge2->x1)*(y1-edge2->y1) - (edge2->y2-edge2->y1)*(x1-edge2->x1) == 0 &&
			(edge2->x2-x1)*(edge2->y2-y1) - (y2-y1)*(x1-edge2->x1) == 0)
		{
			// Coincident edges
			xIntersect = x1;
			yIntersect = y1;
			return true;
		}
		return false;
	}
	else
	{
		ua = ((edge2->x2-edge2->x1)*(y1-edge2->y1) - (edge2->y2-edge2->y1)*(x1-edge2->x1))/denominator;
		
		if (0 <= ua && ua <= 1)
		{
			ub = ((x2-x1)*(y1-edge2->y1) - (y2-y1)*(x1-edge2->x1)) / denominator;
			if (0 <= ub && ub <= 1)
			{
				xIntersect = x1 + ua*(x2-x1);
				yIntersect = y1 + ua*(y2-y1);
				return true;			
			}
		}
	}

	return false;
}

bool inline CSimulatedDiffusion::isColliding(double x, double y, double dx, double dy, vector<Collision>* collisionList)
{	
	int							index0, index1;
	
	static vector<int>			hashmapIndexes;
	static unsigned int nx		= unsigned int(boxWidth/deltaHashmap);

	// Determine indexes of (x, y) and (x+dy, y+dy)
	index0 = getHashmapIndex(x, y);
	index1 = getHashmapIndex(x+dx, y+dy);

	// Check how many boxes are affected
	if (index0 == index1)
	{
		// Add to list, if not emty
		if (spatialHashmap[index0].size() != 0)
		{
			// Clear list of hashmap indexes
			hashmapIndexes.clear();
			hashmapIndexes.push_back(index0);
		}
		else
		{
			return false;
		}
	}
	else
	{		
		// Horizontal or vertical movement, the two boxes are sufficient
		// Other cases are treated here
		if (!(abs(index0-index1) == 1 || abs(index0-index1) == nx))
		{
			hashmapIndexes.clear();
			getHashmapIndexList(x, y, x+dx, y+dy, &hashmapIndexes);
			//*textCtrl << int(hashmapIndexes.size());
		}
		else
		{
			hashmapIndexes.clear();
			hashmapIndexes.push_back(index0);
			hashmapIndexes.push_back(index1);
		}
	}

	unsigned int				n, m; 
	double						xCollision;
	double						yCollision;
	static Collision			collision;
	static Edge*				edge;

	// Clear list of collisions
	collisionList->clear();

	// Check all boxes for collisions
	for (n = 0; n < hashmapIndexes.size(); n++)
	{
		for (m = 0; m < spatialHashmap[hashmapIndexes[n]].size(); m++)
		{					
			if (edgesIntersect(x, y, x+dx, y+dy, spatialHashmap[hashmapIndexes[n]][m], xCollision, yCollision))
			{
				// Test, if element is already in collisionList
				// This can happen since neighbouring hashmap boxes can contain the same edges
				edge = spatialHashmap[hashmapIndexes[n]][m];				

				if (!containsEdge(collisionList, edge))
				{
					collision.edge		= edge;
					collision.x			= xCollision;
					collision.y			= yCollision;
					collision.distance	= ((x-xCollision)*(x-xCollision)+(y-yCollision)*(y-yCollision));
					if(collision.distance > 0)
					{
						collisionList->push_back(collision);			
					}
				}
			}
		}
	}

	if (collisionList->size())
	{
		if (collisionList->size() > 1)
		{
			sort(collisionList->begin(), collisionList->end());
		}
		return true;
	}

	return false;
}

bool CSimulatedDiffusion::containsEdge(vector<Collision>* collisionList, Edge* edge)
{
	if (collisionList->size() == 0)
	{
		return false;
	}

	static unsigned int n;

	for (n = 0; n < collisionList->size(); n++)
	{
		if (collisionList->at(n).edge == edge)
		{
			return true;
		}
	}

	return false;
}

void inline CSimulatedDiffusion::getHashmapIndexList(double x1, double y1, double x2, double y2, vector<int>* indexes)
{
	double xTemp, yTemp;
	int nx, n;

	// Bring arguments in order, so that y1 < y2 and x1 < x2
	if (x1 > x2)
	{
		xTemp = x1;
		x1 = x2;
		x2 = xTemp;
	}
	if (y1 > y2)
	{
		yTemp = y1;
		y1 = y2;
		y2 = yTemp;
	}

	yTemp = y1;
	while (getHashmapIndex(x1, yTemp) <= getHashmapIndex(x1, y2) && yTemp <= boxWidth)
	{	
		nx = getHashmapIndex(x2, yTemp);
		for (n = getHashmapIndex(x1, yTemp); n <= nx; n++)
		{
			indexes->push_back(n);
		}

		yTemp += deltaHashmap;
	}
}


// Fast floor, see: http://ldesoras.free.fr/doc/articles/rounding_en.pdf
int inline floor_int(double x)
{
	const float round_towards_m_i = -0.5f;
	int i;
	__asm
	{
		fld x
		fadd st, st (0)
		fadd round_towards_m_i
		fistp i
		sar i, 1
	}
	return (i);
}

int inline CSimulatedDiffusion::getHashmapIndex(double x, double y)
{
	static unsigned int nx		= unsigned int(boxWidth/deltaHashmap);

	if (x < 0 || x > boxWidth)
	{
		x = x - boxWidth * floor(x / boxWidth);
	}
	if (y < 0 || y > boxHeight)
	{
		y = y - boxHeight * floor(y / boxHeight);
	}

	// Fast replacement for <math.h> floor, which is extremely slow
	return floor_int(x/deltaHashmap) + floor_int(y/deltaHashmap)*nx; 
}

void inline CSimulatedDiffusion::reflect(double& x, double& y, double& dx, double& dy, Collision* collision)
{
	double nx, ny;
	double wx, wy;
	double dot;
	static double epsilon = 1e-6;
	static double oneMinusEpsilon = 1-epsilon;


	// Normal vector of edge
	nx	= collision->edge->nx;
	ny	= collision->edge->ny;

	collision->edge->nx = collision->edge->nx;
	collision->edge->ny = collision->edge->ny;

	wx = x + dx - collision->x;
	wy = y + dy - collision->y;

	// Dot product
	dot = wx*nx + wy*ny;
	
	// New update vector
	// Remaining part of old dx, dy after reflection	
	dx = (wx - 2*nx*dot)*oneMinusEpsilon;
	dy = (wy - 2*ny*dot)*oneMinusEpsilon;

	// Move already a very small step (fraction epsilon) in the new direction
	// Otherwise the paricle would collide again with the same edge
	x = collision->x + dx*epsilon;
	y = collision->y + dy*epsilon;
}


bool CSimulatedDiffusion::loadFile(string filename)
{
	string					line, str;
	fstream					filestr;
	bool					edgeData = false;
	unsigned int			edgeCount	= 0;
	
	double					l;		

	filestr.open(filename.c_str(), fstream::in);

	if (!filestr.good())
	{
		#ifdef CONSOLE_VERSION
		cout << "Voronoi mesh file " << filename <<" not found!\n";
		#endif
		return false;
	}

	deltaHashmap			= 250e-9;

	// Count edges
	while (filestr.good())
	{
		getline(filestr,line);
		if (strcmp(line.c_str(), "EDGES END") == 0)
		{
			edgeData = false;
			break;
		}
		if (edgeData)
		{
			edgeCount++;
		}
		if (strcmp(line.c_str(), "EDGES START") == 0)
		{
			edgeData = true;
		}
	}

	// Allocate memory and read edges
	if (edgeCount)
	{
		edges	= new Edge[edgeCount];
		nEdges	= edgeCount;
	}

	// Read edges
	edgeCount = 0;
	filestr.seekg(ios_base::beg);
	while (filestr.good())
	{
		getline(filestr,line);
		if (strcmp(line.c_str(), "EDGES END") == 0)
		{
			edgeData = false;
			break;
		}
		if (edgeData)
		{
			str = line.substr(0, line.find("\t"));
			edges[edgeCount].x1= atof(str.c_str());
			edges[edgeCount].x1 *= 1e-9;

			line = line.substr(line.find("\t")+1, line.length());
			str = line.substr(0, line.find("\t"));
			edges[edgeCount].y1 = atof(str.c_str());
			edges[edgeCount].y1 *= 1e-9;

			line = line.substr(line.find("\t")+1, line.length());
			str = line.substr(0, line.find("\t"));
			edges[edgeCount].x2 = atof(str.c_str());
			edges[edgeCount].x2 *= 1e-9;

			line = line.substr(line.find("\t")+1, line.length());
			str = line.substr(0, line.find("\t"));
			edges[edgeCount].y2 = atof(str.c_str());
			edges[edgeCount].y2 *= 1e-9;

			// Normal vector
			edges[edgeCount].nx = edges[edgeCount].y2-edges[edgeCount].y1;
			edges[edgeCount].ny = edges[edgeCount].x1-edges[edgeCount].x2;

			l = sqrt(edges[edgeCount].nx*edges[edgeCount].nx+edges[edgeCount].ny*edges[edgeCount].ny);
			edges[edgeCount].nx /= l;
			edges[edgeCount].ny /= l;

			edgeCount++;
		}
		if (strcmp(line.c_str(), "EDGES START") == 0)
		{
			edgeData = true;
		}
	}
	buildHashmap();

	#ifdef CONSOLE_VERSION
	cout << "Voronoi mesh file " << filename <<" loaded.\n";
	#endif

	filestr.close();

	return true;
}

// Calculate autocorrelation of input using the multiple tau algorithm
// Result is returned in G and tau vectors
void CSimulatedDiffusion::multipleTauCorrelate(double* countrateArray, unsigned int nPoints, double deltaT, vector<double>& G, vector<double>& tau)
{
	unsigned int j, m, n, step;
	double Fav, Gtmp;
	vector<double> countrate;
	
	G.clear();
	tau.clear();

	try
	{
		countrate.reserve(nPoints);
	}
	catch (exception& e)
	{
		cout << "Memory allocation for multiple tau failed: " << e.what() << endl;
	}
	// Copy array to the countrate vector
	for (n = 0; n < nPoints; n++)
	{
		countrate.push_back(countrateArray[n]);
	}

	// 
	m = 16;
	step = 0;

	// Normalization factor (denominator of G(tau))
	Fav = mean(&countrate);

	// Calculate autocorrelation of first m bins
	for (n = 0; n < m; n++)
	{		
		Gtmp = 0;

		for (j = 0; j + n + 1< countrate.size(); j++)
		{
			Gtmp += countrate[j]*countrate[j+n+1];
		}
		Gtmp = Gtmp/((nPoints-n-1)*Fav*Fav) - 1;

		G.push_back(Gtmp);
		tau.push_back(deltaT*(n+1));
	}

	while(countrate.size() >= m)
	{
		step += 1;

		// Bin countrate
		bin(countrate);

		Fav = mean(&countrate);

		for (n = 0; n < m/2; n++)
		{		
			Gtmp = 0;

			for (j = 0; j + n + 1 + m/2 < countrate.size(); j++)
			{
				Gtmp += countrate[j]*countrate[j+m/2+n+1];
			}
			Gtmp = Gtmp/((countrate.size()-n-1-m/2)*Fav*Fav) - 1;
			G.push_back(Gtmp);

			tau.push_back(deltaT*(n+1+m/2)*pow(2.0, int(step)));
		}
	}
}

// Calculates the average of the given array
double CSimulatedDiffusion::mean(vector<double>* countrate)
{
	unsigned int n;
	double result = 0;

	for (n= 0; n < countrate->size(); n++)
	{
		result = result + countrate->at(n);
	}

	return result/countrate->size();
}

// Bin the countrate by combining two bins
void CSimulatedDiffusion::bin(vector<double>& countrate)
{
	vector<double> countrateBinned;	

	unsigned int nPoints = floor_int(countrate.size()/2);
	countrateBinned.reserve(nPoints);

	for (n = 0; n < nPoints; n++)
	{
		countrateBinned.push_back((countrate[2*n]+countrate[2*n+1])/2);
	}
	
	countrate = countrateBinned;
}

// Saves the autocorrelation data
void CSimulatedDiffusion::saveAutocorrelations(string filename, vector<vector<double>*>* G, vector<double>* tau, vector<string>* header)
{
	ofstream file;
	file.open(filename.c_str());
	unsigned int n, m;

	if (file.is_open() && file.good())
	{
		// Write header
		file << "tau";
		for (n = 0; n < header->size(); n++)
		{
			file << "\t" << header->at(n);
		}
		file << "\n";

		// Write data
		for (n = 0; n < tau->size(); n++)
		{
			file << tau->at(n);

			for (m = 0; m < G->size(); m++)
			{
				file << "\t" << G->at(m)->at(n);
			}
			file << "\n";
		}

		file.close();
	}
}

// Calculate averaged autocorrelation data
void CSimulatedDiffusion::calculateAverage(vector<vector<double>*>* G, vector<double>* Gav)
{
	unsigned int n, m;
	double Gtmp;
	Gav->clear();
	
	// Loop over all times
	for (n = 0; n < G->at(0)->size(); n++)
	{
		Gtmp = 0;
		// Loop over the individual curves
		for (m = 0; m < G->size(); m++)
		{
			// Sum up time n of all curves
			Gtmp += G->at(m)->at(n);
		}

		// Divide by number of curves
		Gtmp /= G->size();	
		Gav->push_back(Gtmp);
	}
}

// Calculate filament density (Filament length / Area, 1/µm).
double CSimulatedDiffusion::getDensity()
{
	unsigned int n;
	double length = 0;

	for (n = 0; n < nEdges; n++)
	{
		length = length + sqrt((edges[n].x1-edges[n].x2)*(edges[n].x1-edges[n].x2)+(edges[n].y1-edges[n].y2)*(edges[n].y1-edges[n].y2));
	}

	return length*1e6/(10*10);
}