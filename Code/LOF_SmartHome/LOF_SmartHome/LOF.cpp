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
	system("pause");
}
