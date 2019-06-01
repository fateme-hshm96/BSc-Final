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
	if (p->flag)
		return p->kDist;

	double dist;

	for (auto& o : *dataset) {	// find the distnace of each point in dataset with point p
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
	for (auto o : p->neighborhood) {
		//printf("%d\n", o.id);
	}

	p->flag = true;

	return kDist;
}

double reachDist(Point p, Point o) {
	//return max(o.kDist, euclideanDist(p, o));
	printf("reachDist: pid: %d , oid: %d - okdist: %f , opdist: %f\n", p.id, o.id, o.kDist, manhattanDist(p, o));
	return max(o.kDist, manhattanDist(p, o));
}

double localReachabilityDensity(vector<Point> *dataset, Point *p, int K) {

	double sigmaRD = 0;

	double rd;
	for (auto& o : p->neighborhood) {
		rd = reachDist(*p, o);
		sigmaRD += rd;
	}

	printf("LRD: sigmaRD: %f, size: %d, sigmaRD/size: %f, 1/(sigmaRD/size): %f\n", sigmaRD, p->neighborhood.size(),
		(sigmaRD / p->neighborhood.size()), 1.0 / (sigmaRD / p->neighborhood.size()));
	return 1.0 / (sigmaRD / p->neighborhood.size());
}

double LOF(vector<Point> *dataset, Point *p, int K) {
	//// for objects deep inside a cluster, their LOFs are close to 1 
	kDistanse(dataset, p, K);

	for (auto& n1 : p->neighborhood) {
		kDistanse(dataset, &n1, K);
	}
	for (auto& n1 : p->neighborhood) {
		for (auto& n2 : n1.neighborhood) {
			kDistanse(dataset, &n2, K);
		}
	}
	double lrdP = localReachabilityDensity(dataset, p, K);

	double sigmaLRDneigh = 0;
	for (auto& o : p->neighborhood) {
		kDistanse(dataset, &o, K);
		sigmaLRDneigh += localReachabilityDensity(dataset, &o, K);
	}

	printf("\nLOF: sigmaLRDneigh: %f\n", sigmaLRDneigh);
	printf("\nLOF: lrdP: %f\n", lrdP);
	printf("LOF: p.neigh.size: %f\n", (double)(p->neighborhood.size()));

	return  ((sigmaLRDneigh / lrdP) / (double)(p->neighborhood.size()));
}

int main(int argc, char *argv[]) {

	int par_rank, par_size;

	//// Init MPI 
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &par_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &par_size);

	int K = 6;
	vector<Point> dataset;

	//// read .csv file
	ifstream  data;
	data.open("test.csv", ios::in);

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
		dataset.push_back(p);

		Point p1;
		p1.x = x * 100 + 50;		p1.y = y * 100 + 50;
		p1.id = id++;		p1.flag = false;
		dataset.push_back(p1);

		Point p2;
		p2.x = x * 200 + 40;		p2.y = y * 200 + 40;
		p2.id = id++;		p2.flag = false;
		dataset.push_back(p2);

		Point p3;
		p3.x = x * 200;		p3.y = y * 200;
		p3.id = id++;		p3.flag = false;
		dataset.push_back(p3);
	}

	printf("size: %d\n", dataset.size());

	Point p;
	p.x = 1 * dataset.at(2).x + 10;	p.y = 1 * dataset.at(2).y + 10;
	p.flag = false;					p.id = 1188;


	//double kDistP = kDistanse(&dataset, &p, K);
	//printf("kDistP: %f\n", kDistP);

	//double LRD = localReachabilityDensity(&dataset, &p, K);
	//printf("LRD: %f\n", LRD);

	double lof = LOF(&dataset, &p, K);
	printf("LOF: %f\n", lof);

	system("pause");

	//// Finalize MPI 
	MPI_Finalize();

	return 0;
}