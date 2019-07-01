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
	int neighborhood_size;
	int neighborhood[20];
};

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
		p.neighborhood_size = 0;
		(*dataset).push_back(p);
	}

}

bool compareTwoPoints(Point a, Point b) {
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
	
	if (p->flag) {	/// k-distance for this point has been already calculated
		//printf("k-distance for this point(%d) has been already calculated\n", p->id);
		return p->kDist;
	}

	p->flag = true;

	double dist;

	for (auto& o : *dataset) {	/// find the distnace of each point in dataset with point p
		dist = manhattanDist(o, *p);
		//dist = euclideanDist(o, *p);
		if (dist != 0)
			o.dist = dist;
		else
			o.dist = 100000000000;
	}
	sortDataset(dataset);

	double kDist = dataset->at(K - 1).dist;
	p->kDist = kDist;

	for (int i = 0; i < K && i < dataset->size(); i++) {
		p->neighborhood[i] = (dataset->at(K-1-i)).id;
		p->neighborhood_size++;
	}

	int i = K;
	while (dataset->at(i).dist == kDist && i < dataset->size()) {
		p->neighborhood[i] = dataset->at(i).id;
		p->neighborhood_size++;
		i++;
	}

	return kDist;
}

double reachDist(Point p, Point o) {
	//return max(o.kDist, euclideanDist(p, o));

	//printf("reachDist: pid: %d , oid: %d - okdist: %f , opdist: %f\n",
	//p.id, o.id, o.kDist, manhattanDist(p, o));

	return max(o.kDist, manhattanDist(p, o));
}

double localReachabilityDensity(vector<Point> *dataset, Point *p, int K) {
//	printf("localReachabilityDensity for p id %d\n", p->id);

	double sigmaRD = 0;
	int i = 0;
	double rd;

	for (int i = 0; i < p->neighborhood_size; i++) {
		//printf("localReachabilityDensity:\tp id: %d - nei id: %d\n", p->id, dataset->at(p->neighborhood[i]).id);
		rd = reachDist(*p, (dataset->at(p->neighborhood[i])));
		sigmaRD += rd;
	}
	
//	printf("\tLRD: sigmaRD: %f,   size: %d,   sigmaRD/size: %f,   1/(sigmaRD/size): %f\n",
//		sigmaRD, p->neighborhood_size, (sigmaRD / p->neighborhood_size),
//		1.0 / (sigmaRD / p->neighborhood_size));

	return 1.0 / (sigmaRD / p->neighborhood_size);
}

double LOF(vector<Point> *dataset, vector<Point> *dataset_copy, Point *p, int K) {
	//// for objects deep inside a cluster, their LOFs are close to 1 ////

	kDistanse(dataset_copy, p, K);						/// k-distance for p

	printf("p(%d): ", p->id);
	for (int i = 0; i < p->neighborhood_size; i++) {	/// k-distance for p's neighbors
		printf("%d\t", p->neighborhood[i]);
		kDistanse(dataset_copy, &(dataset->at(p->neighborhood[i])), K);
	}
	printf("\n");

	for (int i = 0; i < p->neighborhood_size; i++) {	/// k-distance for neighbors of p's neighbors
		for (int j = 0; j < dataset->at(p->neighborhood[i]).neighborhood_size; j++) {
			kDistanse(dataset_copy, &(dataset->at(dataset->at(p->neighborhood[i]).neighborhood[j])), K);
		}
	}

	double lrdP = localReachabilityDensity(dataset, p, K);

	int i = 0;
	double sigmaLRDneigh = 0;

	for (int i = 0; i < p->neighborhood_size; i++) {
		sigmaLRDneigh += localReachabilityDensity(dataset, &(dataset->at(p->neighborhood[i])), K);
	}

//	printf("\nLOF: sigmaLRDneigh: %f\n", sigmaLRDneigh);
//	printf("\nLOF: lrdP: %f\n", lrdP);
//	printf("LOF: p.neigh.size: %f\n", (double)(p->neighborhood_size));

	return  ((sigmaLRDneigh / lrdP) / (double)(p->neighborhood_size));
}

void LOFBounds(vector<Point> *dataset, vector<Point> *dataset_copy, Point *p, int K) {

	for (int i = 0; i < p->neighborhood_size; i++) {
		for (int j = 0; j < dataset->at(p->neighborhood[i]).neighborhood_size; j++) {
			for (int k = 0; k < dataset->at(dataset->at(p->neighborhood[i]).neighborhood[j]).neighborhood_size; k++) {
				kDistanse(dataset_copy,
					&(dataset->at((dataset->at(dataset->at(p->neighborhood[i]).neighborhood[j])).neighborhood[k])), K);
			}
		}
	}

	double directMin = 10000, directMax = 0.0001, rd;
	double indirectMin = 10000, indirectMax = 0.0001;

	for (int i = 0; i < p->neighborhood_size; i++) {
		rd = reachDist(*p, dataset->at(p->neighborhood[i]));

		if (rd < directMin)
			directMin = rd;

		if (rd > directMax)
			directMax = rd;
	}

	for (int i = 0; i < p->neighborhood_size; i++) {
		for (int j = 0; j < dataset->at(p->neighborhood[i]).neighborhood_size; j++) {
			rd = reachDist(dataset->at(p->neighborhood[i]), dataset->at(dataset->at(p->neighborhood[i]).neighborhood[j]));

			if (rd < indirectMin)
				indirectMin = rd;

			if (rd > indirectMax)
				indirectMax = rd;
		}
	}
	printf("\nmin: %f, max: %f\n", directMin / indirectMax, directMax / indirectMin);
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

	MPI_Status stat;
	

	int K = 5;
	vector<Point> dataset, dataset_copy;

//	readCSV(&dataset, "F:\\Fatemeh\\UNI\\testMPI\\testMPI\\test.csv");
//	printf("dataset size: %d\n", dataset.size());

	int dataset_size;
	if (world_rank == 0) {
		readCSV(&dataset, "F:\\Fatemeh\\UNI\\testMPI\\testMPI\\test.csv");
		printf("dataset size: %d\n", dataset.size());
		dataset_size = dataset.size();
	}

	MPI_Bcast(&dataset_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (world_rank != 0) {
		dataset.resize(dataset_size);
	}

	printf("rankkkk: %d - dataset_size: %d , dataset size: %d\n", world_rank, dataset_size, dataset.size());

	MPI_Bcast((void*)dataset.data(), dataset_size * sizeof(Point), MPI_BYTE, 0, MPI_COMM_WORLD);

	MPI_Barrier(MPI_COMM_WORLD);

	dataset_copy = dataset;

	/*
	if (world_rank == 0) {
		MPI_Send((void*)dataset.data(), dataset_size * sizeof(Point), MPI_BYTE, 1, 0, MPI_COMM_WORLD);
		printf("sent");
	}
	else if (world_rank == 1) {
		printf("about to rcv");
		MPI_Recv((void*)dataset.data(), dataset_size * sizeof(Point), MPI_BYTE, 0, 0, MPI_COMM_WORLD, &stat);
	}
	*/

	Point p;

	p.x = 10.55;
	p.y = -19.61;
	p.flag = false;
	p.id = 1188;
	p.neighborhood_size = 0;

	double lof = LOF(&dataset, &dataset_copy, &p, K);
	printf("\nrank: %d ------ LOF: %f\n", world_rank, lof);

	LOFBounds(&dataset, &dataset_copy, &p, K);


	printf("Hello world from rank %d out of %d processors\n", world_rank, world_size);

		system("pause");


	/// Finalize the MPI environment. ///
	MPI_Finalize();

	return 0;
}