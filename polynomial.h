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
		Monomial(std::string term);
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
		mpz_class round(const mpf_class&);
		std::vector<std::string> split(std::string s);
		mpf_class diff(const mpf_class& a, const mpf_class& b);
	public:
		Polynomial(std::string poly);
		int getDegree() const;
		void addTerm(Monomial m);
		Polynomial derivative() const;
		mpf_class eval(mpf_class x) const;
		mpf_class root(mpf_class lower, mpf_class upper);

	private:
		friend std::ostream& operator<<(std::ostream& os, const Polynomial& p);

};

std::vector<mpf_class> allRoots(Polynomial p, mpf_class lower, mpf_class upper);
