#include <iostream>

using namespace std;

#include "polynomial.h"

int main() {
	cout.precision(1000);
	string poly;
	cin >> poly;
	Polynomial p(poly);
	cout << p << endl;
	vector<mpf_class> rs = allRoots(p, -1000000000.0, 1000000000.0);
	for (int i = 0; i < rs.size(); i++) {
		cout << rs[i] << " ";
	}
	cout << endl;
}
