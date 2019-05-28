#include <iostream>
#include <vector>
#include <algorithm>
#include <iomanip>
#include <fstream>
#include <string>

using namespace std;

struct Point {
	double x;
	double y;
	double kDist;
	double dist;
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

	return kDist;
}


int main() {
	int K = 5;
	vector<Point> dataset;

	//// read .csv file
	ifstream  data;
	data.open("test.csv", ios::in);

	vector<vector<string> > dataList;
	string line;
	while (getline(data, line)) {

		string token = line.substr(0, line.find(','));
		line.erase(0, line.find(',') + 1);

		double x = atof(token.c_str());
		double y = atof(line.c_str());

		Point p;
		p.x = x;
		p.y = y;
		dataset.push_back(p);
	}

	printf("size: %d\n", dataset.size());

	Point p;
	p.x = 3;
	p.y = 4;

	double kDistP = kDistanse(&dataset, &p, K);
	printf("kDistP: %f\n", kDistP);

	system("pause");
}
