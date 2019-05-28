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
		p.dist = x * y + 10; //just to test sort function
		dataset.push_back(p);
	}

	sortDataset(&dataset);

	for (Point p : dataset) {
		printf("%f\n", p.dist);
	}
	printf("size: %d\n", dataset.size());
	system("pause");
}
