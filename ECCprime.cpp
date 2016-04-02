/*#include <NTL/ZZ_p.h>
#include <NTL/pair.h>

using namespace NTL;
using namespace std;

typedef Pair<int, int> ECCpoint;

struct ECCpoint {
	ECCpoint(int a, int b) : a(a), b(b) {}
	int a, b;
};

inline int square(const int & x) {
	return x * x;
}

int p;

int M(int x, int y) {
	return MulMod(x, y, p);
}

int A(int x, int y) {
	return AddMod(x, y, p);
}

int S(int x, int y) {
	return SubMod(x, y, p);
}

int I(int x) {
	return InvMod(x, p);
}

int D(int x, int y) {
	return M(x, I(y));
}

class ECCprime {
public:
	ECCprime(int A, int B, int P) : a(A), b(B), p(P) {}
	const ECCpoint & add(const ECCpoint& P1, const ECCpoint& P2) {
		int X = S(S(M(D(S(P2.b, P1.b), S(P2.a, P1.a)), D(S(P2.b, P1.b), S(P2.a, P1.a))), P1.a), P2.a);
		return ECCpoint(X, S(M(D((S(P2.b, P1.b)), S(P2.a, P1.a)), S(P1.a, X)), P1.b));
	}
	const ECCpoint & doubling(const ECCpoint& P) {
		int tmp = D((A(M(3, M(P.a, P.a)), a)), M(2, P.b));
		int X = S(M(tmp, tmp), M(2, P.a));
		return ECCpoint(X, S(M(tmp, S(P.a, X)), P.b));
	}
	const ECCpoint & negate(const ECCpoint& P) {
		return ECCpoint(P.a, -P.b);
	}
	ECCpoint mult(ZZ& c, const ECCpoint P) {
		while (!(bit(c, 1))) {
			P = doubling(P);
			RightShift(c, 1);
		}
		ECCpoint Q = P;
		P = doubling(P);
		RightShift(c, 1);
		while (c != 0) {
			if (!(bit(c, 1))) {
				Q = add(Q, P);
			}
			P = doubling(P);
			RightShift(c, 1);
		}
		return Q;
	}

private:
	int a, b;
	int p;
};

int main() {
	cin >> p;
	int a, b;
	int x1, y1, x2, y2;
	cin >> a >> b;
	ECCprime E(a, b, p);
	cin >> x1 >> y1;
	//cin >> x2 >> y2;
	ECCpoint P(x1, y1);
	ECCpoint res = E.doubling(P);
	cout << res.a << " " << res.b << endl;
	system("pause");
	return 0;
}*/

#include <NTL/ZZ_pE.h>
#include <NTL/pair.h>
#include <NTL/ZZ_p.h>

using namespace NTL;
using namespace std;

typedef Pair<ZZ_p, ZZ_p> ECCpoint;
/*
struct ECCpoint {
	ECCpoint(ZZ_pE a, ZZ_pE b) : a(a), b(b) {}
	ZZ_pE a, b;
};*/

inline ZZ_p square(const ZZ_p & x) {
	return x * x;
}

class ECCprime {
public:
	ECCprime(ZZ_p A, ZZ_p B) : a(A), b(B) {	
	}
	ECCpoint add(const ECCpoint& P1, const ECCpoint& P2) {
		ZZ_p X = square((P2.b - P1.b) / (P2.a - P1.a)) - P1.a - P2.a;
		ZZ_p Y = ((P2.b - P1.b) / (P2.a - P1.a)) * (P1.a - X) - P1.b;
		return ECCpoint(X, Y);
	}
	ECCpoint doubling(const ECCpoint& P) {
		ZZ_p tmp = (3 * square(P.a) + a) / (2 * P.b);
		ZZ_p X = square(tmp) - 2 * P.a;
		return ECCpoint(X, tmp * (P.a - X) - P.b);
	}
	ECCpoint negate(const ECCpoint& P) {
		return ECCpoint(P.a, -P.b);
	}
	ECCpoint mult(ZZ c, ECCpoint P) {
		while (!(bit(c, 0))) {
			P = doubling(P);
			c >>= 1;
		}
		ECCpoint Q = P;
		P = doubling(P);
		c >>= 1;
		while (c != 0) {
			if (bit(c, 0)) {
				Q = add(Q, P);
			}
			P = doubling(P);
			c >>= 1;
		}
		return Q;
	}

private:
	ZZ_p a, b;
};

#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEX.h>

#define NUM_BRANCHES 32

class triple {
	private:
		ZZ Aj,Bj,Rj;
	public:
		triple(ZZ Aj,ZZ Bj,ZZ Rj) : Aj(Aj), Bj(Bj), Rj(Rj) {}
		~triple() {
			delete Aj;
			delete Bj;
			delete Rj;
		}
		ZZ getAj() {
			return Aj;
		}
		ZZ getBj() {
			return Bj;
		}
		ZZ getRj() {
			return Rj;
		}
};

ZZ partition_function(ECCpoint P) {
	return (P.a % NUM_BRANCHES)+1;
}

ZZ pollard_attack(ECCpoint P, ECCpoint Q, ZZ order_p) {
	ZZ Aj,Bj;
	triple triples_array[] = new triple[NUM_BRANCHES];
	for(int j=0; j<NUM_BRANCHES; j++) {
		Aj=RandomBnd(order_p);
		Bj=RandomBnd(order_p);
		triples[i]=new triple(Aj,Bj,add(mult(Aj,P),mult(Bj,Q));
	}
	ZZ ctag, dtag, ctags, dtags;
	ECCpoint Xtag, Xtags;
	ctags = RandomBnd(order_p);
	dtags = RandomBnd(order_p);
	Xtags = add(mult(ctags,P),mult(dtags,Q));
	Xtag=Xtags;
	ZZ index;
	do {
		index = partition_function(Xtag);
		Xtag = add(Xtag, (triples_array[index]).getRj);
		ctag = (ctag + (triples_array[index]).getAj)%order_p;
		dtag = (dtag + (triples_array[index]).getBj)%order_p;
		for(int i = 0; i<=1; i++) {
			index = partition_function(Xtags);
			Xtags = add(Xtags, (triples_array[index]).getRj);
			ctags = (ctags + (triples_array[index]).getAj)%order_p;
			dtags = (dtags + (triples_array[index]).getBj)%order_p;
		}
	} while (Xtag!=Xtags);
	if(dtag == dtags) {
		cerr << "Error" << endl;
		return 0;
	}	else {
		return ((ctag-ctags)/(dtag-dtags))%order_p;
	}	
}

int main() {
   ZZ p;
   cin >> p;
   ZZ_p::init(p);
	ZZ_p a, b;
	ZZ_p x1, y1, x2, y2;
	cin >> a >> b;
	ECCprime E(a, b);
	cin >> x1 >> y1;
	cin >> x2 >> y2;
	ECCpoint P1(x1, y1), P2(x2, y2);
	/*ECCpoint res = E.add(P1, P2);
	cout << res.a << " " << res.b << endl;
	ZZ c;
	cin >> c;
	ECCpoint res2 = E.mult(c, P1);
	cout << res2.a << " " << res2.b << endl;*/
	ZZ c;
	cin >> c;
	while (P1 != E.mult(c, P2)) {
		c++;
		cout << c << endl;
	}
	cout << c << endl;
	system("pause");
	return 0;
}