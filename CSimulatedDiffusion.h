#pragma once

#include <vector>
#include <string>

//#define SPEED_TEST
#define CONSOLE_VERSION

#ifdef SPEED_TEST
	#include <time.h>
#endif

using namespace std;

struct Edge
{
	// Coordinates
	double x1;
	double y1;
	double x2;
	double y2;

	// Normal vector
	double nx;
	double ny;
};

class Collision
{
public:
	Edge*	edge;
	double	distance;
	double  x;
	double  y;
	bool operator<(Collision& collision){return distance < collision.distance;}
};

class CSimulatedDiffusion
{
public:
	CSimulatedDiffusion();
	~CSimulatedDiffusion();

	// Start the simulation
	void			start();

	// Load a voronoi mesh
	bool			loadFile(string filename);

#ifdef GUI_VERSION
	void			setFrame(CFrame* frame);
#endif

	string			outputFileStr;

	unsigned int	nParticles;
	unsigned int	n;
	unsigned int	nUpdate;

#ifdef GUI_VERSION
	CFrame*			frame;
#endif

	// Vector containing particle positions
	double*			x;
	double*			y;

	// Jump probability
	double			pjump;

	// Array containg the voronoi edges
	Edge*			edges;
	unsigned int	nEdges;

	// Diffusion coefficient
	double			D;
	
	// Stepsize of simulation
	double			dt;
	
	// Total simulation time
	double			Tmax;
	
	// Total number of simulation steps
	unsigned int	nSteps;
	
	// Size of simulation box
	double			boxWidth;
	double			boxHeight;

	// Positions of detection spots
	double			spotx1, spotx2, spotx3;
	double			spoty1, spoty2, spoty3;

	// Storage for trajectorys
	double*			xTraj;
	double*			yTraj;

	// Spatial hashmap
	double			deltaHashmap;	
	vector<vector<Edge*>>	spatialHashmap;	

	// Initialize spatial hashmap
	void			buildHashmap();

	// Calculate autocorrelation of input using the multiple tau algorithm
	// Result is returned in G and tau vectors
	void			multipleTauCorrelate(double* countrateArray, unsigned int nPoints, double deltaT, vector<double>& G, vector<double>& tau);

	// Saves the autocorrelation data
	void			saveAutocorrelations(string filename, vector<vector<double>*>* G, vector<double>* tau, vector<string>* header);

	// Calculate average correlation
	void			calculateAverage(vector<vector<double>*>* G, vector<double>* Gav);

	// Calculates the average of the given array
	double			mean(vector<double>* countrate);

	// Bin the countrate vector by combining two bins
	void			bin(vector<double>& countrate);

	// Determine if edge is contained in a given rectangle
	bool			edgeInBox(Edge* edge, double x1, double y1, double x2, double y2);
	
	// Determine if edges intersect
	bool			edgesIntersect(Edge* edge1, Edge* edge2);
	bool			edgesIntersect(double x1, double y1, double x2, double y2, Edge* edge2, double& xIntersect, double& yIntersect);

	// Determine if given vector intersects with mesh. Generates collision vector.
	bool			isColliding(double x, double y, double dx, double dy, vector<Collision>* collisionList);

	// Reflect the particle according to the given collision. Updates x, y, dx, dy
	void			reflect(double& x, double& y, double& dx, double& dy, Collision* collision);

	// Get a list of all hashmap Indexes which are inside the rectange x1, y1, x2, y2
	void			getHashmapIndexList(double x1, double y1, double x2, double y2, vector<int>* indexes);
	
	// Get the hashmap index of x and y
	int				getHashmapIndex(double x, double y);
	
	// Test if collisionList contains a collision with edge
	bool			containsEdge(vector<Collision>* collisionList, Edge* edge);

	// Calculate filament density (Filament length / Area, 1/µm).
	double			getDensity();
};

