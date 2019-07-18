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

vector<Point>* sub_KNN(vector<Point>* sub_dataset, Point* p, int K, int num_elements) {
	vector<Point> result;
	
	double dist;
	for (int i = 0; i < num_elements; i++) {
		dist = manhattanDist(sub_dataset->at(i), *p);
		//dist = euclideanDist(sub_dataset->at(i), *p);
		if (sub_dataset->at(i).id != p->id)
			sub_dataset->at(i).dist = dist;
		else
			sub_dataset->at(i).dist = 100000000000;
	}
	sortDataset(sub_dataset);

	for (int i = 0; i < K; i++)
		result.push_back(sub_dataset->at(i));

	return &result;
}

double kDistanse(vector<Point> *dataset, Point *p, int K) {
	/// compute-intensive step of the LOF algorithm ///
	
	if (p->flag) {	/// k-distance for this point has been already calculated
	//	printf("k-distance for this point(%d) has been already calculated\n", p->id);
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

double reachDist(vector<Point> *dataset_copy, Point p, Point o, int K) {
	//double kDist = kDistanse(dataset_copy, &o, K);
	
	//return max(o.kDist, euclideanDist(p, o));

	//printf("reachDist: pid: %d , oid: %d - okdist: %f , opdist: %f\n",
	//p.id, o.id, o.kDist, manhattanDist(p, o));

	return max(o.kDist, manhattanDist(p, o));
}

double localReachabilityDensity(vector<Point> *dataset, vector<Point> *dataset_copy, Point *p, int K) {
//	printf("localReachabilityDensity for p id %d\n", p->id);

	double sigmaRD = 0;
	int i = 0;
	double rd;

	for (int i = 0; i < p->neighborhood_size; i++) {
		//printf("localReachabilityDensity:\tp id: %d - nei id: %d\n", p->id, dataset->at(p->neighborhood[i]).id);
		rd = reachDist(dataset_copy, *p, (dataset->at(p->neighborhood[i])), K);
		sigmaRD += rd;
	}
	
	printf("\tLRD: id: %d -  sigmaRD: %f,   size: %d,   sigmaRD/size: %f,   1/(sigmaRD/size): %f\n",
		p->id, sigmaRD, p->neighborhood_size, (sigmaRD / p->neighborhood_size),
		1.0 / (sigmaRD / p->neighborhood_size));

	return 1.0 / (sigmaRD / p->neighborhood_size);
}


void parallel_kDistance(vector<Point> *dataset_copy,
	Point *p, int K, int rank, int world_size) {
	vector<Point> result, final_res;

	for (auto& q : *dataset_copy) {
		double dist = manhattanDist(q, *p);
		//dist = euclideanDist(q, *p);
		if (q.id != p->id)
			q.dist = dist;
		else
			q.dist = 100000000000;
	}

	sortDataset(dataset_copy);

	int i;
	for (i = 0; i < K + 2 && i < dataset_copy->size(); i++) {
		result.push_back(dataset_copy->at(i));
	}

	double lastDist = result.back().dist;
	/*
	while (dataset_copy->at(i).dist == lastDist && i < dataset_copy->size()) {
	result.push_back(dataset_copy->at(i));
	i++;
	}
	*/
	int res_size = result.size() * sizeof(Point);
	//printf("rank: %d, res_size: %d\n\n", rank, result.size());

	int *recvcounts = NULL;
	int * offset = NULL;
	if (rank == 0) {
		recvcounts = (int*)malloc(world_size * sizeof(int));
		offset = (int*)malloc(world_size * sizeof(int));
	}
	MPI_Gather(&res_size, 1, MPI_INT, recvcounts, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		offset[0] = 0;
		for (int i = 1; i < world_size; i++) {
			offset[i] = offset[i - 1] + recvcounts[i - 1] + 1;
			//	printf("offset %d: %d --- %d\n", i, offset[i - 1], recvcounts[i - 1]);
		}
		int final_res_size = (offset[world_size - 1] + recvcounts[world_size - 1]) / sizeof(Point);
		//printf("offset i-1: %d, recvciont i-1: %d, final res size: %d\n", offset[world_size - 1],
		//	recvcounts[world_size - 1], final_res_size);
		final_res.resize(final_res_size);
	}

	MPI_Gather((void*)result.data(), result.size() * sizeof(Point), MPI_BYTE,
		(void*)final_res.data(), result.size() * sizeof(Point), MPI_BYTE, 0, MPI_COMM_WORLD);

	//	MPI_Gatherv((void*)result.data(), result.size() * sizeof(Point), MPI_BYTE,
	//	(void*)final_res.data(), recvcounts, offset, MPI_BYTE, 0, MPI_COMM_WORLD);

	if (rank == 0) {
		sortDataset(&final_res);
		i = 0;
		//	printf("parallel NEIGHBOR ID");
		for (i = 0; i < K && i < final_res.size(); i++) {
			p->neighborhood[i] = final_res.at(i).id;
			//	printf(" %d - ", p->neighborhood[i]);
			p->neighborhood_size++;
		}

		lastDist = final_res.at(i - 1).dist;
		while (dataset_copy->at(i).dist == lastDist && i < dataset_copy->size()) {
			p->neighborhood[i] = final_res.at(i).id;
			//	printf(" %d - ", p->neighborhood[i]);
			p->neighborhood_size++;
			i++;
		}
		printf("\n");

		p->dist = lastDist;
		if (rank == 0)
			//	printf("nei_size: %d   -  id: %d _________lastDist: %f\n", p->neighborhood_size, p->id, lastDist);
			p->flag = true;
	}
	MPI_Barrier(MPI_COMM_WORLD);
}

void all_kDistance(vector<Point> *dataset, Point *p, int K, int rank, int world_size) {
	vector<Point> sub_dataset;

	int sub_dataset_size;
	if (rank == 0) {
		sub_dataset_size = dataset->size() / world_size;
	}
	MPI_Bcast(&sub_dataset_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

	sub_dataset.resize(sub_dataset_size);

	MPI_Scatter((void*)dataset->data(), sub_dataset_size * sizeof(Point), MPI_BYTE, (void*)sub_dataset.data(),
		sub_dataset_size * sizeof(Point), MPI_BYTE, 0, MPI_COMM_WORLD);

	//printf("rank: %d, sub: %f - %f\n", rank, sub_dataset.at(0).x, sub_dataset.at(0).y);

	vector<Point> sub_copy;
	sub_copy = sub_dataset;

	parallel_kDistance(&sub_copy, p, K, rank, world_size);

	int neighborhood_size;
	if (rank == 0)
		neighborhood_size = p->neighborhood_size;
	MPI_Bcast(&neighborhood_size, 1, MPI_INT, 0, MPI_COMM_WORLD);

	for (int i = 0; i < neighborhood_size; i++) {	/// k-distance for p's neighbors
		Point q;
		if (rank == 0) {
			q.flag = false;
			q.neighborhood_size = 0;
			q.id = (dataset->at(p->neighborhood[i])).id;
			q.x = (dataset->at(p->neighborhood[i])).x;
			q.y = (dataset->at(p->neighborhood[i])).y;
			//printf("qx: %f, qy: %f\n", q.x, q.y);
		}
		MPI_Bcast(&q, sizeof(Point), MPI_BYTE, 0, MPI_COMM_WORLD);

		parallel_kDistance(&sub_copy, &q, K, rank, world_size);

		if (rank == 0) {
			(dataset->at(p->neighborhood[i])).neighborhood_size = q.neighborhood_size;
			for (int j = 0; j < q.neighborhood_size; j++) {
				(dataset->at(p->neighborhood[i])).neighborhood[j] = q.neighborhood[j];
			}
		}
	}
	printf("\n");


	for (int i = 0; i < neighborhood_size; i++) {	/// k-distance for neighbors of p's neighbors
		int neighborhood_size2;
		if (rank == 0) {
			neighborhood_size2 = dataset->at(p->neighborhood[i]).neighborhood_size;
			printf("\t%d : neighborhood_size2 : %d\n", p->neighborhood[i], neighborhood_size2);
		}
		MPI_Bcast(&neighborhood_size2, 1, MPI_INT, 0, MPI_COMM_WORLD);

		for (int j = 0; j < neighborhood_size2; j++) {
			Point q;
			if (rank == 0) {
				q.flag = false;
				q.neighborhood_size = 0;
				q.id = (dataset->at(dataset->at(p->neighborhood[i]).neighborhood[j])).id;
				q.x = (dataset->at(dataset->at(p->neighborhood[i]).neighborhood[j])).x;
				q.y = (dataset->at(dataset->at(p->neighborhood[i]).neighborhood[j])).y;
			}
			MPI_Bcast(&q, sizeof(Point), MPI_BYTE, 0, MPI_COMM_WORLD);

			parallel_kDistance(&sub_copy, &q, K, rank, world_size);

			if (rank == 0) {
				(dataset->at(dataset->at(p->neighborhood[i]).neighborhood[j])).neighborhood_size = q.neighborhood_size;
				for (int k = 0; k < q.neighborhood_size; k++) {
					(dataset->at(dataset->at(p->neighborhood[i]).neighborhood[j])).neighborhood[k] = q.neighborhood[k];
				}
			}

		}
	}
	printf("%d are we done?\n", rank);
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
		rd = reachDist(dataset_copy,*p, dataset->at(p->neighborhood[i]), K);

		if (rd < directMin)
			directMin = rd;

		if (rd > directMax)
			directMax = rd;
	}

	for (int i = 0; i < p->neighborhood_size; i++) {
		for (int j = 0; j < dataset->at(p->neighborhood[i]).neighborhood_size; j++) {
			rd = reachDist(dataset_copy, dataset->at(p->neighborhood[i]), 
				dataset->at(dataset->at(p->neighborhood[i]).neighborhood[j]), K);

			if (rd < indirectMin)
				indirectMin = rd;

			if (rd > indirectMax)
				indirectMax = rd;
		}
	}
	printf("\nmin: %f, max: %f\n", directMin / indirectMax, directMax / indirectMin);
}

double LOF(vector<Point> *dataset, vector<Point> *dataset_copy, Point *p, int K) {
	//// for objects deep inside a cluster, their LOFs are close to 1 ////

	/*
	kDistanse(dataset_copy, p, K);						/// k-distance for p

	
	printf("p(%d): ", p->id);
	for (int i = 0; i < p->neighborhood_size; i++) {	/// k-distance for p's neighbors
		printf("%d\t", p->neighborhood[i]);
		kDistanse(dataset_copy, &(dataset->at(p->neighborhood[i])), K);
	}
	printf("\n");
	*/

	double lrdP = localReachabilityDensity(dataset, dataset_copy, p, K);

	int i = 0;
	double sigmaLRDneigh = 0;

	for (int i = 0; i < p->neighborhood_size; i++) {
		sigmaLRDneigh += localReachabilityDensity(dataset, dataset_copy, &(dataset->at(p->neighborhood[i])), K);
	}

//	printf("\nLOF: sigmaLRDneigh: %f\n", sigmaLRDneigh);
//	printf("\nLOF: lrdP: %f\n", lrdP);
//	printf("LOF: p.neigh.size: %f\n", (double)(p->neighborhood_size));

	return  ((sigmaLRDneigh / lrdP) / (double)(p->neighborhood_size));
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

	//	double in_x = atof(argv[1]);
	//	double in_y = atof(argv[2]);

	int K = 5;
	vector<Point> dataset, copy;

	if (world_rank == 0) {
		readCSV(&dataset, "F:\\Fatemeh\\UNI\\testMPI\\testMPI\\test.csv");
		printf("dataset size: %d\n", dataset.size());
	}

	Point p;
	// p.x = in_x;	p.y = in_y;
	p.x = 80.55;	p.y = -129.61;
	p.flag = false;
	p.id = 99999;
	p.neighborhood_size = 0;

	all_kDistance(&dataset, &p, K, world_rank, world_size);
	
	if (world_rank == 0) {
		copy = dataset;
		double lof = LOF(&dataset, &copy, &p, K);
		printf("%f", lof);
	}

	//LOFBounds(&dataset, &copy, &p, K);

	//system("pause");

	/// Finalize the MPI environment. ///
	MPI_Finalize();

	return 0;
}