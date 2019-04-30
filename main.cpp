#include <iostream>
#include <regex>

using namespace std;

#include "polynomial.h"

int main() {
	string poly;
	cin >> poly;

	Polynomial p = parsePolynomial(poly);

	cout << p << endl;
	
	cout.precision(100);
	vector<mpf_class> rs = allRoots(p, -1000000000.0, 1000000000.0);
	for (int i = 0; i < rs.size(); i++) {
		cout << rs[i] << " ";
	}
	cout << endl;
}
