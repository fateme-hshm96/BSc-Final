#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <string>
#include <mpi.h>

//#include "gnuplot-iostream.h"

using namespace std;

struct Point {
	int id; 		bool flag;
	double x; 		double y;
	double kDist;	double dist;
	vector<Point> neighborhood;
};

bool compareTwoPoints(Point a, Point b){
	return a.dist < b.dist;
}

void sortDataset(vector<Point> *dataset) {
	sort(dataset->begin(), dataset->end(), compareTwoPoints);
}

double manhattanDist(Point a, Point b) {
	return abs(a.x - b.x) + abs(a.y - b.y);
}

double euclideanDist(Point a, Point b) {
	return sqrt(pow(a.x - b.x, 2) + pow(a.y - b.y, 2));
}

double kDistanse(vector<Point> *dataset, Point *p, int K) {
	/// compute-intensive step of the LOF algorithm ///

	if (p->flag)	/// k-distance for this point has been already calculated
		return p->kDist;

	double dist;

	for (auto& o : *dataset) {	/// find the distnace of each point in dataset with point p
		dist = manhattanDist(o, *p);
		//dist = euclideanDist(o, *p);
		o.dist = dist;
	}

	sortDataset(dataset);

	double kDist = dataset->at(K - 1).dist;
	p->kDist = kDist;

	p->neighborhood.clear();
	for (int i = 0; i < K && i < dataset->size(); i++) {
		p->neighborhood.push_back((dataset->at(i)));
	}

	int ii = K;
	while (dataset->at(ii).dist == kDist && ii < dataset->size()) {
		p->neighborhood.push_back(dataset->at(ii));
		ii++;
	}

	//printf("p(%d) neighbors:\n\n", p->id);
	//for (auto o : p->neighborhood) {
		//printf("%d\n", o.id);
	//}

	p->flag = true;

	return kDist;
}

double reachDist(Point p, Point o) {
	//return max(o.kDist, euclideanDist(p, o));

	//printf("reachDist: pid: %d , oid: %d - okdist: %f , opdist: %f\n",
		//p.id, o.id, o.kDist, manhattanDist(p, o));

	return max(o.kDist, manhattanDist(p, o));
}

double localReachabilityDensity(vector<Point> *dataset, Point *p, int K) {

	double sigmaRD = 0;

	double rd;
	for (auto& o : p->neighborhood) {
		rd = reachDist(*p, o);
		sigmaRD += rd;
	}

	printf("LRD: sigmaRD: %f, size: %d, sigmaRD/size: %f, 1/(sigmaRD/size): %f\n",
		sigmaRD, p->neighborhood.size(), (sigmaRD / p->neighborhood.size()),
		1.0 / (sigmaRD / p->neighborhood.size()));

	return 1.0 / (sigmaRD / p->neighborhood.size());
}

double LOF(vector<Point> *dataset, Point *p, int K) {
	//// for objects deep inside a cluster, their LOFs are close to 1 ////

	kDistanse(dataset, p, K);			/// k-distance for p

	for (auto& n1 : p->neighborhood) {	/// k-distance for p's neighbors
		kDistanse(dataset, &n1, K);
	}

	for (auto& n2 : p->neighborhood) {	/// k-distance for neighbors of p's neighbors
		for (auto& n1 : n2.neighborhood) {
			kDistanse(dataset, &n1, K);
		}
	}

	
	double lrdP = localReachabilityDensity(dataset, p, K);

	double sigmaLRDneigh = 0;
	for (auto& o : p->neighborhood) {
		sigmaLRDneigh += localReachabilityDensity(dataset, &o, K);
	}

	printf("\nLOF: sigmaLRDneigh: %f\n", sigmaLRDneigh);
	printf("\nLOF: lrdP: %f\n", lrdP);
	printf("LOF: p.neigh.size: %f\n", (double)(p->neighborhood.size()));

	return  ((sigmaLRDneigh / lrdP) / (double)(p->neighborhood.size()));
}

void LOFBounds(vector<Point> *dataset, Point *p, int K) {

	for (auto& n3 : p->neighborhood) {	/// k-distance for neighbors of neighbors of p's neighbors
		for (auto& n2 : n3.neighborhood) {
			for (auto& n1 : n2.neighborhood) {
				kDistanse(dataset, &n1, K);
			}
		}
	}

	double directMin = 10000, directMax = 0.0001, rd;
	double indirectMin = 10000, indirectMax = 0.0001;

	for (auto& n1 : p->neighborhood) {
		rd = reachDist(*p, n1);

		if (rd < directMin)
			directMin = rd;

		if (rd > directMax)
			directMax = rd;
	}

	for (auto& n1 : p->neighborhood) {
		for (auto& n2 : n1.neighborhood) {
			rd = reachDist(n1, n2);

			if (rd < indirectMin)
				indirectMin = rd;

			if (rd > indirectMax)
				indirectMax = rd;
		}
	}
	printf("\nmin: %f, max: %f\n", directMin / indirectMax, directMax / indirectMin);
}

void readCSV(vector<Point> *dataset, char* file_name) {
	//// read .csv file ////
	ifstream  data;
	data.open(file_name, ios::in);

	vector<vector<string> > dataList;
	string line;
	int id = 0;
	while (getline(data, line)) {

		string token = line.substr(0, line.find(','));
		line.erase(0, line.find(',') + 1);

		double x = atof(token.c_str());		double y = atof(line.c_str());

		Point p;
		p.x = x * 100;		p.y = y * 100;
		p.id = id++;		p.flag = false;
		(*dataset).push_back(p);
	}

}

int main(int argc, char *argv[]) {

	/// Initialize the MPI environment ///
	MPI_Init(&argc, &argv);

	/// Get the number of processes ///
	int world_size;
	MPI_Comm_size(MPI_COMM_WORLD, &world_size);

	/// Get the rank of the process ///
	int world_rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

	int K = 5;
	vector<Point> dataset;
	Point p;
	p.x = 0;		p.y = 0;
	p.flag = false;	p.id = 121212;

	if (world_rank == 0) {
		readCSV(&dataset, "F:\\Fatemeh\\UNI\\testMPI\\testMPI\\test.csv");

		printf("dataset size: %d\n", dataset.size());

		//p.x = 1 * dataset.at(2).x + 10;	p.y = -101 * dataset.at(2).y + 10;
		//p.flag = false;						p.id = 1188;
	}

	//double lof = LOF(&dataset, &p, K);
	//printf("\nLOF: %f\n", lof);

	//LOFBounds(&dataset, &p, K);
	

	printf("Hello world from rank %d out of %d processors\n", world_rank, world_size);

	//system("pause");

	int token;
	if (world_rank != 0) {
		MPI_Recv(&token, 1, MPI_INT, world_rank - 1, 0,
			MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("Process %d received token %d from process %d\n",
			world_rank, token, world_rank - 1);
	}
	else {
		// Set the token's value if you are process 0
		token = -1;
	}
	MPI_Send(&token, 1, MPI_INT, (world_rank + 1) % world_size,
		0, MPI_COMM_WORLD);

	// Now process 0 can receive from the last process.
	if (world_rank == 0) {
		MPI_Recv(&token, 1, MPI_INT, world_size - 1, 0,
			MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		printf("Process %d received token %d from process %d\n",
			world_rank, token, world_size - 1);
	}

	/// Finalize the MPI environment. ///
	MPI_Finalize();

	return 0;
}