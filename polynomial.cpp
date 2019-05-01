#include "polynomial.h"
#include <sstream>
#include <cctype>
#include <iostream>
#include <regex>

static const mpf_class PRECISION = 0.000000000000000000001;

int Monomial::getDegree() const {
	return degree;
}

std::string Monomial::getVar() const {
	return var;
}

mpf_class Monomial::getCof() const {
	return cof;
}

Monomial Monomial::derivative() const {
	if (degree == 0) {
		return Monomial(0, "", 0);
	}
	if (degree == 1) {
		return Monomial(cof, "", 0);
	}
	return Monomial(cof*degree, var, degree-1);	
}

mpf_class Monomial::eval(mpf_class x) const {
	if (var == "") {
		return cof;
	} else {
		mpf_class a, b = 7;
		mpf_pow_ui(a.get_mpf_t(), x.get_mpf_t(), degree);
		return cof*a;
	}
}

mpz_class Polynomial::round(const mpf_class& f) const {
	mpz_class floor(f);
	mpz_class ceil = floor;
	if (f > 0) {
		ceil += 1;
		if (f-floor < ceil-f)
			return floor;
		return ceil;
	} else {
		ceil -= 1;
		if (f-floor > ceil-f)
			return floor;
		return ceil;
	}
}

int Polynomial::getDegree() const {
	return degree;
}

void Polynomial::addTerm(Monomial m) {
	terms.push_back(m);
	if (m.getDegree() > degree)
		degree = m.getDegree();
}

Polynomial Polynomial::derivative() const {
	Polynomial deriv;
	for (int i = 0; i < terms.size(); i++) {
		Monomial d = terms[i].derivative();
		if (d.getCof() != 0) {
			deriv.addTerm(d);
		}
	}	
	return deriv;
}
mpf_class Polynomial::eval(mpf_class x) const {
	mpf_class result = 0.0;
	for (int i = 0; i < terms.size(); i++) {
		result += terms[i].eval(x);
	}
	return result;
}

mpf_class Polynomial::root(mpf_class lower, mpf_class upper) const {
	if (eval(lower)*eval(upper) > 0) //No root
		return upper+1;

	if (eval(lower) == 0)
		return lower;
	if (eval(upper) == 0)
		return upper;

	bool rising = eval(lower) < eval(upper) ? true : false;
	mpf_class r = (upper+lower)/2.0;
	mpf_class result = eval(r);
	mpf_class prev = 0.0;

	while (diff(result, prev) > PRECISION) {
		if (rising) {
			if (result < 0) {
				lower = r;
			} else {
				upper = r;
			}
		} else {
			if (result < 0) {
				upper = r;
			} else {
				lower = r;
			}
		}
		r = (upper+lower)/2.0;
		prev = result;
		result = eval(r);
	}

	mpz_class integer = round(r);
	if (eval(integer) == 0) {
		return integer;
	}
	return r;
}

mpf_class Polynomial::diff(const mpf_class& a, const mpf_class& b) const {
	mpf_class temp1, temp2;
	mpf_abs(temp1.get_mpf_t(), a.get_mpf_t());
	mpf_abs(temp2.get_mpf_t(), b.get_mpf_t());
	temp1 -= temp2;
	mpf_abs(temp2.get_mpf_t(), temp1.get_mpf_t());
	return temp2;
}

std::ostream& operator<<(std::ostream& os, const Polynomial& p) {
	for (int i = 0; i < p.terms.size(); i++) {
		os << p.terms[i].getCof() << p.terms[i].getVar();
		if (p.terms[i].getDegree() > 1) {
			os << '^' << p.terms[i].getDegree();
		}
		os << " + ";
	}
	return os;
}

std::vector<mpf_class> allRoots(const Polynomial& p, mpf_class lower, mpf_class upper) {
	std::vector<mpf_class> roots, points;
	if (p.getDegree() > 1) {
		points = allRoots(p.derivative(), lower, upper);
		mpf_class r;

		if (points.size() == 0) {
			roots.push_back(p.root(lower, upper));
			return roots;
		}

		//First
	       	r = p.root(lower, points[0]);
		if (r != points[0]+1) { // Root exists
			roots.push_back(r);
		}
		
		// From 0 to i
		for (int i = 1; i < points.size(); i++) {
			r = p.root(points[i-1], points[i]);
			if (r != points[i]+1) {
				roots.push_back(r);
			}
		}
		
		// Last
		r = p.root(points[points.size()-1], upper);
		if (r != upper+1) {
			roots.push_back(r);
		}
	} else {
		mpf_class r = p.root(lower, upper);
		if (r != upper+1)
			roots.push_back(p.root(lower, upper));
	}
	return roots;
}

Polynomial parsePolynomial(std::string s) {
	Polynomial out;
	std::regex splite("([+-]?[^-+]+)");
	std::regex termparse("([-+][0-9]+\\.*[[0-9]*)([a-z])\\^*([0-9]+)*|([-+][0-9]+\\.*[[0-9]*)|([+-])([a-z])\\^*([0-9]+)*");
	std::sregex_token_iterator rend;
	std::sregex_token_iterator m(s.begin(), s.end(), splite);
	for (; m != rend; m++) {
		std::smatch matches;
		std::string str = *m;
		std::regex_search(str, matches, termparse);

		int start = 1;
		while (matches.str(start) == "") {start++;}

		mpf_class cof;
		if (matches.str(start) == "-")
			cof = -1;
		else if (matches.str(start) == "+")
			cof = 1;
		else
			cof = matches.str(start)[0] != '+' ? matches.str(start) : matches.str(start).substr(1);
		std::string var = matches.str(start+1);
		int deg = (matches.str(start+2) == "" ? 1 : stoi(matches.str(start+2)));
		out.addTerm(Monomial(cof, var, deg));
	}
	return out;
}

