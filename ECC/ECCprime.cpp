#include <NTL/ZZ_pE.h>
#include <NTL/pair.h>
#include <NTL/ZZ_p.h>

using namespace NTL;
using namespace std;

#include <fstream>

inline ZZ_p square(const ZZ_p & x) {
	return x * x;
}

struct ECCpoint {
	ZZ_p x, y, z;
	ECCpoint(ZZ_p X, ZZ_p Y, ZZ_p Z) : x(X), y(Y), z(Z) {
	}
	ECCpoint(long X, long Y, long Z) : x(conv<ZZ_p>(X)), y(conv<ZZ_p>(Y)), z(conv<ZZ_p>(Z)) {
	}
	ECCpoint() {}
	bool operator==(const ECCpoint& P) const {
		if (z == 0) {
			return P.z == 0;
		}
		else if (P.z == 0) {
			return false;
		}
		ZZ_p z2 = square(z), Pz2 = square(P.z);
		return ((x * Pz2) == (P.x * z2)) && ((y * Pz2 * P.z) == (P.y * z * z2));
	}
	bool operator!=(const ECCpoint& P) const {
		return !(*this == P);
	}
};

ostream& operator<<(ostream& cout, const ECCpoint& P) {
	if (P.z == 0) {
		cout << "infinity" << endl;
	} else {
		cout << "(" << P.x / (P.z * P.z) << ", " << P.y / (P.z * P.z * P.z) << ")";
	}
	return cout;
}

//typedef Pair<ZZ_p, ZZ_p> ECCpoint;
/*
struct ECCpoint {
	ECCpoint(ZZ_pE a, ZZ_pE b) : a(a), b(b) {}
	ZZ_pE a, b;
};*/

class triple;

#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEX.h>

#define NUM_BRANCHES 16
#define BRANCH_MASK 0xf

class triple {
	private:
		ZZ Aj,Bj;
		ECCpoint Rj;
	public:
		triple() {}
		triple(ZZ Aj,ZZ Bj,ECCpoint Rj) : Aj(Aj), Bj(Bj), Rj(Rj) {}
		~triple() {
		}
		ZZ getAj() {
			return Aj;
		}
		ZZ getBj() {
			return Bj;
		}
		ECCpoint getRj() {
			return Rj;
		}
};


int inline partition_function(ECCpoint P) {
	//return (rep(P.x) % NUM_BRANCHES);
	if (P.z == 0) {
		return (rep(P.x) % BRANCH_MASK);
	}
	return (rep(P.x / (P.z * P.z))) % BRANCH_MASK;
}

class ECCprime {
public:
	ECCpoint inf;
	ECCprime(ZZ_p A, ZZ_p B) : a(A), b(B), inf(1, 1, 0) {	
	}
	ECCpoint doubling(const ECCpoint& P) {
		//cout << "*@" << endl;
		//cout << P << endl;
		if ((P.y == 0) || (P.z == 0)) {
			return inf;
		}
		ZZ_p y2 = square(P.y);
		ZZ_p S = 4*P.x*y2;
		ZZ_p M = 3*square(P.x) + a*square(square(P.z));
		ZZ_p X = square(M) - 2*S;
		ZZ_p Y = M*(S - X) - 8*square(y2);
		ZZ_p Z = 2*P.y*P.z;
		//cout << ECCpoint(X, Y, Z) << endl;
		return ECCpoint(X, Y, Z);
		/*ZZ_p tmp = (3 * square(P.x) + a) / (2 * P.y);
		ZZ_p X = square(tmp) - 2 * P.x;
		return ECCpoint(X, tmp * (P.x - X) - P.y, conv<ZZ_p>(1));*/
	}
	ECCpoint add(const ECCpoint& P1, const ECCpoint& P2) {
		//cout << "+" << endl;
		//cout << P1 << " " << P2 << endl;

		if (P1 == inf) {
			return P2;
		}
		if (P2 == inf) {
			return P1;
		}

		ZZ_p z12 = square(P1.z), z22 = square(P2.z);
		//cout << "z12= " << z12 << " z22= " << z22 << endl;
		ZZ_p U1 = P1.x*z22;
		ZZ_p U2 = P2.x*z12;
		//cout << "U1= " << U1 << " U2= " << U2 << endl;
		ZZ_p S1 = P1.y*z22*P2.z;
		ZZ_p S2 = P2.y*z12*P1.z;
		if (U1 == U2) {
			if (S1 != S2) {
				return inf;
			} else {
				return doubling(P1);
			}
		}
		ZZ_p H = U2 - U1;
		ZZ_p R = S2 - S1;
		ZZ_p h2 = square(H);
		ZZ_p h3 = h2 * H;
		ZZ_p tmp = U1*h2;
		ZZ_p X3 = square(R) - h3 - 2*tmp;
		ZZ_p Y3 = R*(tmp - X3) - S1*h3;
		ZZ_p Z3 = H*P1.z*P2.z;
		//cout << ECCpoint(X3, Y3, Z3) << endl;
		return ECCpoint(X3, Y3, Z3);
/*
		ZZ_p X = square((P2.y - P1.y) / (P2.x - P1.x)) - P1.x - P2.x;
		ZZ_p Y = ((P2.y - P1.y) / (P2.x - P1.x)) * (P1.x - X) - P1.y;
		return ECCpoint(X, Y, conv<ZZ_p>(1));*/
	}
	ECCpoint negate(const ECCpoint& P) {
		return ECCpoint(P.x, -P.y, P.z);
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
	/*
	229
1 44
5 116
155 166
*/

//#define DEBUG

	ZZ pollard_attack(ECCpoint P, ECCpoint Q, ZZ order_p, long mySeed) {
		ZZ Aj,Bj; 
		long as[NUM_BRANCHES] = {79, 206, 87, 219}; 
		long bs[NUM_BRANCHES] = {163, 19, 109, 68};
		triple triples_array[NUM_BRANCHES];

		SetSeed(conv<ZZ>(123));

		for(int j=0; j<NUM_BRANCHES; j++) {
#ifdef DEBUG
			Aj=conv<ZZ>(as[j]); //RandomBnd(order_p);
			Bj=conv<ZZ>(bs[j]); //RandomBnd(order_p);
#else
			Aj=RandomBnd(order_p);
			Bj=RandomBnd(order_p);
#endif
			triples_array[j] = triple(Aj,Bj,add(mult(Aj,P),mult(Bj,Q)));
		}
		ZZ ctag, dtag, ctags, dtags;
		ECCpoint Xtag, Xtags;
#ifdef DEBUG
		ctag = ctags = conv<ZZ>(54);//RandomBnd(order_p);
		dtag = dtags = conv<ZZ>(175);// RandomBnd(order_p);
#else
		ctag = ctags = RandomBnd(order_p);
		dtag = dtags = RandomBnd(order_p);
#endif
		Xtag = Xtags = add(mult(ctags,P),mult(dtags,Q));
		cout << "c'= " << ctag << " " << "d'= " << dtag << endl;
		cout << "c''= " << ctags << " " << "d''= " << dtags << endl;
		cout << "x' = " << Xtag << " " << "x''= " << Xtags << endl;

		int index;
		long long ccc = 1;
		do {
			index = partition_function(Xtag);
			Xtag = add(Xtag, triples_array[index].getRj());
			ctag = (ctag + (triples_array[index]).getAj())%order_p;
			dtag = (dtag + (triples_array[index]).getBj())%order_p;
			for(int i = 0; i<=1; i++) {
				index = partition_function(Xtags);
				/*cout << "index= " << index << " " << "R= " << (triples_array[index]).getRj() << endl;
				cout << "R.z= " << (triples_array[index]).getRj().z << endl;
				cout << "x''= " << Xtags << "x''.z= " << Xtags.z << endl;*/
				Xtags = add(Xtags, (triples_array[index]).getRj());
				ctags = (ctags + (triples_array[index]).getAj())%order_p;
				dtags = (dtags + (triples_array[index]).getBj())%order_p;
			}
			ccc++;
			/*cout << "c'= " << ctag << " " << "d'= " << dtag << endl;
			cout << "c''= " << ctags << " " << "d''= " << dtags << endl;
			cout << "x' = " << Xtag << " " << "x''= " << Xtags << endl;*/
		} while (Xtag!=Xtags);
		if(dtag == dtags) {
			cerr << "Error" << endl;
			return conv<ZZ>(0);
		}	else {
			cout << "Iterations: " << ccc << endl;
			cout << ctag << " , " << ctags << " | " << dtag << " , " << dtags << endl;
			SubMod(ctag, ctag, ctags, conv<ZZ>(order_p));
			SubMod(dtags, dtags, dtag, conv<ZZ>(order_p));
			InvMod(dtags, dtags, conv<ZZ>(order_p));
			MulMod(ctag, ctag, dtags, conv<ZZ>(order_p));
			return ctag;
		}	
	}
private:
	ZZ_p a, b;
};

int main() {
   ZZ p;
#ifdef DEBUG
		p = conv<ZZ>(229);
		cout << p << endl;
#else
   cin >> p;
#endif
   ZZ_p::init(p);
	ZZ_p a, b;
	ZZ_p x1, y1, x2, y2;
#ifdef DEBUG
	a = conv<ZZ_p>(1);
	b = conv<ZZ_p>(44);
#else
	cin >> a >> b;
#endif
	ECCprime E(a, b);
#ifdef DEBUG
	x1 = conv<ZZ_p>(5);
	y1 = conv<ZZ_p>(116);
	x2 = conv<ZZ_p>(155);
	y2 = conv<ZZ_p>(166);
#else
	cin >> x1 >> y1;
	cin >> x2 >> y2;
#endif
	ZZ n;
#ifdef DEBUG
	n = conv<ZZ>(239);
#else
	cin >> n;
#endif
	ECCpoint P(x1, y1, conv<ZZ_p>(1)), Q(x2, y2, conv<ZZ_p>(1));
	/*while (true) {
		ZZ ctag, dtag, ctags, dtags;
		cin >> ctag >> dtag >> ctags >> dtags;
		ECCpoint xtag = E.add(E.mult(ctag,P),E.mult(dtag,Q));
		cout << xtag << endl;
		ECCpoint xtags = E.add(E.mult(ctags, P), E.mult(dtags, Q));
		cout << xtags << endl;
	}*/
	/*ECCpoint res = E.add(P1, P2);
	cout << res.a << " " << res.b << endl;
	ZZ c;
	cin >> c;
	ECCpoint res2 = E.mult(c, P1);
	cout << res2.a << " " << res2.b << endl;*/
	cout << "Starting Attack  " <<endl;
	long mySeed;
	cin >> mySeed;
	ZZ l = E.pollard_attack(P, Q, n, mySeed);
	cout << "Private key = " << l << endl;
	cout << "Attack Ended " << endl;
	cout << "Q = " << Q << endl;
	cout << "lP = " << E.mult(l, P) << endl;
	system("pause");
	return 0;
}