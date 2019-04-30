#include <gmpxx.h>
#include <string>
#include <vector>
#include <ostream>

class Monomial {
	private:
		mpf_class cof;
		std::string var;
		int degree;
	public:
		Monomial(mpf_class c = 0, std::string v = "", int deg = 0) : cof(c),  var(v), degree(deg) {};
		int getDegree() const;
		std::string getVar() const;
		mpf_class getCof() const;
		Monomial derivative() const; 
		mpf_class eval(mpf_class x) const;
};

class Polynomial {
	private:
		int degree;
		std::vector<Monomial> terms;
		mpz_class round(const mpf_class&) const;
		mpf_class diff(const mpf_class& a, const mpf_class& b) const;
	public:
		Polynomial() : degree(0) {};
		int getDegree() const;
		void addTerm(Monomial m);
		Polynomial derivative() const;
		mpf_class eval(mpf_class x) const;
		mpf_class root(mpf_class lower, mpf_class upper) const;

	private:
		friend std::ostream& operator<<(std::ostream& os, const Polynomial& p);

};

std::vector<mpf_class> allRoots(const Polynomial& p, mpf_class lower, mpf_class upper);
Polynomial parsePolynomial(std::string);
