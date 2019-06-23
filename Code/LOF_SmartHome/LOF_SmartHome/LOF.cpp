#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <string>
#include <mpi.h>

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
	for (auto o : p->neighborhood) {
		//printf("%d\n", o.id);
	}

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

void readFile() {
	FILE *file;
	file = fopen("test.csv", "r");

	if (file == NULL) {
		printf("no file to open \n");
		return;
	}
	char c = fgetc(file);
	printf("%c\n", c);
}

void readCSV(vector<Point> *dataset) {
	//// read .csv dataset file ////
	float x, y;

	FILE *fp = fopen("test.csv", "r");
	if (fp) {
		char linee[100];
		int i = 0;
		while (fgets(linee, 100, fp)) {
			sscanf(linee, "%f,%f", x, y);
			printf("%s\n", linee);
		}
	}

	ifstream  data;
	data.open("test.csv", ios::in);

	vector<vector<string>> dataList;
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

		Point p1;
		p1.x = x * 100 + 50;		p1.y = y * 100 + 50;
		p1.id = id++;		p1.flag = false;
		(*dataset).push_back(p1);

		Point p2;
		p2.x = x * 200 + 40;		p2.y = y * 200 + 40;
		p2.id = id++;		p2.flag = false;
		(*dataset).push_back(p2);

		Point p3;
		p3.x = x * 200;		p3.y = y * 200;
		p3.id = id++;		p3.flag = false;
		(*dataset).push_back(p3);
		printf("%f - %f\n", p.x, p.y);
	}
	printf("done reading file\n");
}

const char* getfield(char* line, int num)
{
	const char* tok;
	for (tok = strtok(line, ";");
		tok && *tok;
		tok = strtok(NULL, ";\n"))
	{
		if (!--num)
			return tok;
	}
	return NULL;
}


int main(int argc, char *argv[]) {

	FILE* stream = fopen("input", "r");

	char line[1024];
	while (fgets(line, 1024, stream))
	{
		char* tmp = strdup(line);
		printf("Field 3 would be %s\n", getfield(tmp, 3));
		// NOTE strtok clobbers tmp
		free(tmp);
	}
	
	system("pause");
	return 0;
}