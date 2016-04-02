#include <NTL/ZZ_pE.h>
#include <NTL/pair.h>
#include <NTL/ZZ_p.h>
#include <vector>

using namespace std;

using namespace NTL;
using namespace std;

#include <fstream>
#include <windows.h>

ZZ p;

inline ZZ_p square(const ZZ_p & x) {
	return x * x;
}

ZZ R, R1;
long n;
vector<ZZ_p> powers;

ZZ MontMul(const ZZ& a, const ZZ& b, const ZZ& p) {
	ZZ res = conv<ZZ>(0);

	for (long  i=0; i<n; i++) {
		if (bit(a, i)) {
			res += b;
		}
		if (bit(res, 0)) {
			res += p;
		}
		res /= 2;
		
		if (res > p) {
			cout << "res: " << res << endl;
		}
	}
	return res;
}

const ZZ_p MEXP(const ZZ_p& a) {
	ZZ p = ZZ_p::modulus();
	ZZ X = rep(a*conv<ZZ_p>(R));
	ZZ e = p-2;
	ZZ A = R;
	for (long  i=0; i<n; i++) {
		A = MontMul(A, A, p);
		if (bit(e,i) == 1) {
			A = MontMul(X, A, p);
		}
		cout << A << endl;
	}
	return conv<ZZ_p>(A*R1);
}


const ZZ_p& MINV(const ZZ_p& a) {
	ZZ u = ZZ_p::modulus(), v = conv<ZZ>(a), r = conv<ZZ>(0), s = conv<ZZ>(1);
	ZZ p = ZZ_p::modulus();

	long long k=0;
	while (v != 0) {
		if (!bit(u, 0)) {
			u >>= 1;
			s <<= 1;
		} else if (!bit(v, 0)) {
			v >>= 1;
			r <<= 1;
		} else if (u > v) {
			u = (u - v) >>= 1;
			r = r + s;
			s <<= 1;
		} else {
			v = (v - u) >>= 1;
			s = s + r;
			r <<= 1;
		}
		++k;
	}

	ZZ_p ans = conv<ZZ_p>(r);
	ans *= powers[k-n];
	//cout << ans << endl;
	
	return ans;
}

void MDIV(const ZZ_p& a, const ZZ_p& b, ZZ_p& res) {
	ZZ u = ZZ_p::modulus(), v = conv<ZZ>(b), r = conv<ZZ>(0), s = conv<ZZ>(a);
	ZZ p = ZZ_p::modulus();

	long long k=0;
	while (v != 0) {
		if (!bit(u, 0)) {
			u >>= 1;
			s <<= 1;
		} else if (!bit(v, 0)) {
			v >>= 1;
			r <<= 1;
		} else if (u > v) {
			u = (u - v) >>= 1;
			r = r + s;
			s <<= 1;
		} else {
			v = (v - u) >>= 1;
			s = s + r;
			r <<= 1;
		}
		++k;
	}

	res = conv<ZZ_p>(r)*powers[k-n];
}

typedef Pair<ZZ_p, ZZ_p> Affine;

struct Jacobian {
	ZZ_p x, y, z, z2;
	Jacobian(const ZZ_p& X, const ZZ_p& Y, const ZZ_p& Z) : x(X), y(Y), z(Z) {
		z2 = square(z);
	}
	Jacobian(const Affine& P) : x(P.a), y(P.b), z(conv<ZZ_p>(1)), z2(conv<ZZ_p>(1)) {
	}
	Jacobian(long X, long Y, long Z) : x(conv<ZZ_p>(X)), y(conv<ZZ_p>(Y)), z(conv<ZZ_p>(Z)) {
	}
	Jacobian() {}
	bool operator==(const Jacobian& P) const {
		/*
		if (z == 0) {
			return P.z == 0;
		}
		else if (P.z == 0) {
			return false;
		}*/
		return ((x * P.z2) == (P.x * z2)) && ((y * P.z2 * P.z) == (P.y * z * z2));
	}
	bool operator!=(const Jacobian& P) const {
		// can we earn something by implementing this one instead?
		return !(*this == P);
	}
};


ostream& operator<<(ostream& cout, const Jacobian& P) {
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

#ifdef DEBUG
#define NUM_BRANCHES 4
#else
#define NUM_BRANCHES 32
#endif

#ifdef HALF_BRANCH_TEST
#define BRANCH_MASK 0xf
#else
#define BRANCH_MASK (NUM_BRANCHES-1)
#endif

class triple {
private:
	ZZ Aj,Bj;
	Affine Rj;
public:
	triple() {}
	triple(ZZ Aj,ZZ Bj,Affine Rj) : Aj(Aj), Bj(Bj), Rj(Rj) {}
	~triple() {
	}
	const ZZ& getAj() {
		return Aj;
	}
	const ZZ& getBj() {
		return Bj;
	}
	const Affine& getRj() {
		return Rj;
	}
};

//#define H_AB

int inline partition_function(const Affine& P) {
	//return (rep(P.x) % NUM_BRANCHES);
	/*
	ZZ_p a = P.x;

	if (P.z != 0) {
		a *= inv(P.z2);
	}
	
	return conv<int>(rep(a)) & BRANCH_MASK;
	*/
#ifdef H_AB
	return conv<int>(rep(P.a*P.b)) & BRANCH_MASK;
#else
	return conv<int>(rep(P.a)) & BRANCH_MASK;
#endif
}

#define __WORDSIZE 64

class optimize { 
private:
	ZZ n,k,m,p; // n is the number of bits of given number, m equals some integer i times w (wordsize) s.t. m=i*w >= n.
	ZZ r; // the number which we are looking for is inverse.
public:
	optimize(const ZZ& pp) : p(pp) {
		n = NumBits(p);
		m=__WORDSIZE;
		while(m<n) {
			m+=m; 
		}
		// cout<< "passed the c'tor" << endl;
	}
	void AlmMonInv(const ZZ& a) {
		//cout << "n= " << endl;
		//cout << "m= " << endl;
		ZZ u,v;
		ZZ s;
		r=0; s=1;
		u = p;
		v=a;
		k = 0;
		//cout << "starting the loop" << endl;
		//cout << "u = " << u << " v = " << v << " s = " << s << " r = " << r << " k = " << k << endl;
		while(v>0) {
			if(u % 2 == 0) {
				u /=2;	s *= 2;
			} else if(v%2 == 0) {
				v /=2; r*=2;
			} else if(u>v) {
				u=u-v;	u /= 2; r=r+s; s*=2;
			} else if(v>=u) {
				v=v-u;	v /=2; s=s+r; r*=2;
			}
			k++;
		//	cout << "u = " << u << " v = " << v << " r = " << r << " s = " << s << " k = " << k << endl;
		}
	//	cout << "ended the loop" << endl;
		//cout << " k " << k << " r " << r << " p " << p << endl;

		if(r>=p) {
			r = r-p;
		}
		//cout << "r1 = " << r << endl;
		r = p-r;
		//cout << "phase1 ended r= " << r << endl;
	}
	ZZ_p Inv() {
		if(k>m) {
			//cout << " k > m" << endl;
			ZZ powcond = PowerMod(conv<ZZ>(2), -m, p);
			r = MulMod(r,powcond,p);
			k = k-m;
		}
		ZZ powpow = PowerMod(conv<ZZ>(2),-k,p);
		r = MulMod(r,powpow,p);
		ZZ_p $ = conv<ZZ_p>(r);
		//cout << "phase2 ended r = " << $ << endl;
		return $;
	}
};

int inline partition_function(const Jacobian& P) {
	if (P.z == 0) {
		return conv<int>(P.x) & BRANCH_MASK;
	}
	 //*power(P.z2, e);
	//optimize I(ZZ_p::modulus());
	//I.AlmMonInv(conv<ZZ>(P.z2));
	//a *= I.Inv();
		
	//ZZ_p x = MEXP(P.z2);
	//ZZ_p y = inv(P.z2);

	//cout << "x " << x << "y " << y << " R^R1" << R*R1 << endl;
	ZZ_p a = MEXP(P.z2);
	ZZ_p y = inv(P.z2);
	a *= P.x;

	//a *= MINV(P.z2); 
	//a /= P.z2;
	//MDIV(P.x, P.z2, a); 
	return conv<int>(a) & BRANCH_MASK;
}

class ECCprime {
public:
	Jacobian inf;
	ECCprime(ZZ_p A, ZZ_p B,ZZ p) : a(A), b(B), p(p), inf(1, 1, 0) {	
	}

	Jacobian doubling(const Jacobian& P) const {
		if ((P.y == 0) || (P.z == 0)) {
			cout << "inf" << endl;
			return inf;
		}
		ZZ_p y2 = square(P.y);
		ZZ_p S = 4*P.x*y2;
		ZZ_p M = 3*square(P.x) + a*square(P.z2);
		ZZ_p X = square(M) - 2*S;
		return Jacobian(X, M*(S - X) - 8*square(y2), 2*P.y*P.z);
	}

	Jacobian add(const Jacobian& P1, const Affine& P2) {
		if (P1 == inf) {
			return Jacobian(P2.a, P2.b, conv<ZZ_p>(1));
		}
		ZZ_p U2 = P2.a*P1.z2;
		ZZ_p S2 = P2.b*P1.z2*P1.z;
		if (P1.x == U2) {
			if (P1.y != S2) {
				// how many times?
				cout << "inf" << endl;
				return inf;
			} else {
				return doubling(P1);
			}
		}
		U2 -= P1.x;
		S2 -= P1.y;
		ZZ_p h2 = square(U2);
		ZZ_p h3 = h2 * U2;
		ZZ_p tmp = P1.x*h2;
		ZZ_p X3 = square(S2) - h3 - 2*tmp;

		return Jacobian(X3, S2*(tmp - X3) - P1.y*h3, U2*P1.z);
	}

	void add_i(Jacobian& P1, const Affine& P2) {
		if (P1 == inf) {
			P1.x = P2.a;
			P1.y = P2.b;
			P1.z = conv<ZZ_p>(1);
			return;
		}
		ZZ_p U2 = P2.a*P1.z2;
		ZZ_p S2 = P2.b*P1.z2*P1.z;
		if (P1.x == U2) {
			if (P1.y != S2) {
				P1 = inf;
			} else {
				P1 = doubling(P1);
			}
			return;
		}
		ZZ_p H = U2 - P1.x;
		ZZ_p R = S2 - P1.y;
		ZZ_p h2 = square(H);
		ZZ_p h3 = h2 * H;
		ZZ_p tmp = P1.x*h2;
		P1.x = square(R) - h3 - 2*tmp;
		P1.y = R*(tmp - P1.x) - P1.y*h3;
		P1.z *= H;
	}

	Jacobian add(const Jacobian& P1, const Jacobian& P2) const {
		if (P1 == inf) {
			return P2;
		}
		if (P2 == inf) {
			return P1;
		}
		ZZ_p U1 = P1.x*P2.z2;
		ZZ_p U2 = P2.x*P1.z2;
		ZZ_p S1 = P1.y*P2.z2*P2.z;
		ZZ_p S2 = P2.y*P1.z2*P1.z;
		if (U1 == U2) {
			if (S1 != S2) {
				cout << "inf" << endl;
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
		return Jacobian(X3, R*(tmp - X3) - S1*h3, H*P1.z*P2.z);
	}
	Jacobian negate(const Jacobian& P) const {
		return Jacobian(P.x, -P.y, P.z);
	}
	Jacobian mult(ZZ c, Jacobian P) const {
		while (!(bit(c, 0))) {
			P = doubling(P);
			c >>= 1;
		}
		Jacobian Q = P;
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
    Affine add(const Affine& P1, const Affine& P2) {
		if (P2.a == P1.a) {
			return Affine(conv<ZZ_p>(0),conv<ZZ_p>(0));
		}
		ZZ_p tmp = (P2.b - P1.b) / (P2.a - P1.a);
		ZZ_p X = square(tmp) - P1.a - P2.a;
		return Affine(X, tmp * (P1.a - X) - P1.b);
	}
	Affine doubling(const Affine& P) {
		ZZ_p tmp = (3 * square(P.a) + a) / (2 * P.b);
		ZZ_p X = square(tmp) - 2 * P.a;
		return Affine(X, tmp * (P.a - X) - P.b);
	}
	Affine negate(const Affine& P) {
		return Affine(P.a, -P.b);
	}
	Affine mult(ZZ c, const Affine& P1) {
		if (c == 0) {
			return Affine(conv<ZZ_p>(0),conv<ZZ_p>(0));
		}
		Affine P = P1;
		while (!(bit(c, 0))) {
			P = doubling(P);
			c >>= 1;
		}
		Affine Q = P;
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

#define INPUT_TYPE Affine
#define SKIP_INF 1

	ZZ pollardAttack(INPUT_TYPE P, INPUT_TYPE Q, ZZ order_p) {
		ZZ ctag, dtag, ctags, dtags, ctagss, dtagss;
		INPUT_TYPE Xtag, Xtags, Xtagss;
		//ZZ e = ZZ_p::modulus()-2;
		
#ifdef DEBUG
		ctag = ctags = conv<ZZ>(54);//RandomBnd(order_p);
		dtag = dtags = conv<ZZ>(175);// RandomBnd(order_p);
#else
		ctag = ctags = RandomBnd(order_p);
		dtag = dtags = RandomBnd(order_p);
#endif
		Xtag = Xtags = add(mult(ctags,P),mult(dtags,Q));
		cout << "c'= " << ctag << " d'= " << dtag << " x' = " << Xtag << endl;
		/*
		cout << "c''= " << ctags << " " << "d''= " << dtags << endl;
		cout << "x' = " << Xtag << " " << "x''= " << Xtags << endl;
		*/
		DWORD startTime = ::GetTickCount();

		int index;
		long long ccc = 1;
		long long iii = 1;

#ifdef FLOYD_CYCLE
		do {
			index = partition_function(Xtag);
			Xtag = add(Xtag, triples_array[index].getRj());
			ctag = (ctag + (triples_array[index]).getAj())%order_p;
			dtag = (dtag + (triples_array[index]).getBj())%order_p;
			ccc++;
			for(int i = 0; i<=1; i++) {
				index = partition_function(Xtags);
				/*cout << "index= " << index << " " << "R= " << (triples_array[index]).getRj() << endl;
				cout << "R.z= " << (triples_array[index]).getRj().z << endl;
				cout << "x''= " << Xtags << "x''.z= " << Xtags.z << endl;*/
				Xtags = add(Xtags, (triples_array[index]).getRj());
				ctags = (ctags + (triples_array[index]).getAj())%order_p;
				dtags = (dtags + (triples_array[index]).getBj())%order_p;
				ccc++;
			}

			/*cout << "c'= " << ctag << " " << "d'= " << dtag << endl;
			cout << "c''= " << ctags << " " << "d''= " << dtags << endl;
			cout << "x' = " << Xtag << " " << "x''= " << Xtags << endl;*/
		} while (Xtag!=Xtags);
#else // BRENT_CYCLE
		long long dist = 2, lamda = 2, half_way = 1;

		index = partition_function(Xtag);
#ifdef SKIP_INF
		while (triples_array[index].getRj().a == Xtag.a) {
			cout << "avoiding inf" << endl;
			index = (index+1) & BRANCH_MASK;
		}
#endif
		Xtagss = add(Xtag, (triples_array[index]).getRj());
		ctagss = (ctag + (triples_array[index]).getAj())%order_p;
		dtagss = (dtag + (triples_array[index]).getBj())%order_p;

		index = partition_function(Xtagss);

#ifdef SKIP_INF
		while (triples_array[index].getRj().a == Xtagss.a) {
			cout << "avoiding inf" << endl;
			index = (index+1) & BRANCH_MASK;
		}
#endif
		Xtags = add(Xtagss, (triples_array[index]).getRj());
		ctags = (ctagss + (triples_array[index]).getAj())%order_p;
		dtags = (dtagss + (triples_array[index]).getBj())%order_p;

		while (Xtag!=Xtags) {
			if (Xtags == Xtagss) {
				Xtag = Xtagss;
				ctag = ctagss;
				dtag = dtagss;
				cout << "half-way hit after " << lamda << " steps out of " << dist << endl;
				break;
			}

			if (dist == lamda) {
				Xtag = Xtags;
				ctag = ctags;
				dtag = dtags;
				half_way = dist;
				dist *= 2;
				lamda = 0;
			} else if (lamda == half_way) {
				Xtagss = Xtags;
				ctagss = ctags;
				dtagss = dtags;
			} 

			index = partition_function(Xtags);
			//cout << "index: " << index << endl;
#ifdef SKIP_INF
			while (triples_array[index].getRj().a == Xtags.a) {
				cout << "avoiding inf" << endl;
				index = (index+1) & BRANCH_MASK;
			}
#endif
			// we can calc histogram on c/d, and compute the result only at the end
			// 32 multiplications, 32 additions
			// but the hist is also inc (addition)...
			Xtags = add(Xtags, (triples_array[index]).getRj());
			ctags += triples_array[index].getAj();
			ctags %= order_p;
			dtags += triples_array[index].getBj();
			dtags %= order_p;
			//dtags = (dtags + (triples_array[index]).getBj())%order_p;

			//ctags = (ctags + (triples_array[index]).getAj())%order_p;
			//dtags = (dtags + (triples_array[index]).getBj())%order_p;

			lamda++;
			ccc++;
			iii++;
			if (iii == 1000000000) {
				cout << "iteration count: " << ccc << endl;
				cout << "Xtag: " << Xtag << endl;
				cout << "Xtags: " << Xtags << endl;
				cout << "Xtagss: " << Xtagss << endl;
				cout << lamda << " " << dist << endl;
				iii = 0;
			}
		};
		cout << "exiting after " << lamda << " steps out of " << dist << endl;

#endif

		DWORD endTime = ::GetTickCount();

		if(dtag == dtags) {
			cerr << "Error" << endl;
			return conv<ZZ>(0);
		}	else {
			cout << "Iterations: " << ccc << endl;
			cout << "ExecutionTime: " << double(endTime - startTime)/1000.0 << " sec " << endl;

			cout << ctag << " , " << ctags << " | " << dtag << " , " << dtags << endl;

			SubMod(ctag, ctag, ctags, conv<ZZ>(order_p));
			SubMod(dtags, dtags, dtag, conv<ZZ>(order_p));
			InvMod(dtags, dtags, conv<ZZ>(order_p));
			MulMod(ctag, ctag, dtags, conv<ZZ>(order_p));
			return ctag;
		}	
	}

	void initHashFunction(Affine P, Affine Q, ZZ order_p) {
	
#ifdef DEBUG
		long as[NUM_BRANCHES] = {79, 206, 87, 219}; 
		long bs[NUM_BRANCHES] = {163, 19, 109, 68};
#endif

		ZZ Aj,Bj; 
		Affine Rj;
		for(int j=0; j<NUM_BRANCHES; j++) {
#ifdef DEBUG
			Aj=conv<ZZ>(as[j]); //RandomBnd(order_p);
			Bj=conv<ZZ>(bs[j]); //RandomBnd(order_p);
#else
			//do {
				Aj=RandomBnd(order_p);
				Bj=RandomBnd(order_p);
				Rj = add(mult(Aj,P),mult(Bj,Q));
				/*
				cout << "A[" << j << "] = " << Aj << endl; 
				cout << "B[" << j << "] = " << Bj << endl; 
				cout << "R[" << j << "] = " << Rj << endl; 
				*/
			//} while (Rj.z == 0);
#endif
			triples_array[j] = triple(Aj,Bj,Rj); //Affine(Rj.x / Rj.z2, Rj.y / (Rj.z2 * Rj.z)));
		}
	}

private:
	ZZ_p a, b;
	ZZ p;
	triple triples_array[NUM_BRANCHES];
};

#include <vector>

using namespace std;

ZZ chineseRemainder(vector<ZZ> a, vector<ZZ> b, long size) {
	ZZ N = conv<ZZ>(1);
	vector<ZZ> Ni(size), Di(size);
	ZZ x = conv<ZZ>(0);
	for(long i = 0; i < size; i++) {
		N*=b[i];
	}
	for(long i = 0; i < size; i++) {
		Ni[i] = N/b[i];
		ZZ temp = Ni[i];
		rem(temp, temp, b[i]);
		Di[i] = InvMod(temp,b[i]);
	}
	for(long i = 0; i < size; i++) {
		x = x + (a[i]*Ni[i]*Di[i]);
	}
	x = MulMod(x,1,N);
	return x;
}

ZZ tryOrder(const vector<ZZ>& factors, int n, int k, ZZ current, const ECCprime& E, const Jacobian& P) {
	if (k == n) {
		if ((E.mult(current, P)).z == 0) {
			return current;
		} else {
			return conv<ZZ>(0);
		}
	} else {
		ZZ op1 = tryOrder(factors, n, k + 1, current, E, P);
		if (op1 != 0) {
			return op1;
		}
		return tryOrder(factors, n, k + 1, current * factors[k], E, P);
	}
}

ZZ tryOrder(const vector<ZZ>& factors, int n, const ECCprime& E, const Jacobian& P) {
	return tryOrder(factors, n, 0, conv<ZZ>(1), E, P);
}

void factor(ZZ n, vector<ZZ>& factors) {
	ZZ p = conv<ZZ>(2);
	while ((p * p) <= n) {
		while ((n % p) == 0) {
			n /= p;
			factors.push_back(p);
		}
		++p;
	}
	if (n > 1) {
		factors.push_back(n);
	}
}

ZZ findOrder(const ZZ& p, const ECCprime& E, const Jacobian& P) {
	for (int k = 1; k <= 6; ++k) {
		vector<ZZ> factors;
		factor(power(p, k) - 1, factors);
		ZZ order = tryOrder(factors, factors.size(), E, P);
		if (order != 0) {
			return order;
		}
	}
	return conv<ZZ>(0); // error
}

void findLP(const ZZ& n, vector<ZZ>& ls, vector<ZZ>& ps, ECCprime E, const Affine& P, const Affine& Q) {
	vector<ZZ> factors;
	factor(n, factors);
	int f = factors.size();
	ZZ c1, zi, c2;
	Affine c2p, Pol, P0;
	for (int i = 0; i < f; ++i) {
		if ((!i) || factors[i] != factors[i - 1]) {
			c1 = n / factors[i];
			Pol = Q;
			c2p = P;
			c2 = conv<ZZ>(1);
			P0 = E.mult(c1, P);
			Affine Q0 = E.mult(c1, Pol);
			/*cout << P0 << " " << Q0 << endl;
			system("pause");*/
			if (P0 == Q0) {
				zi = conv<ZZ>(1);
			} else if ((Q0.a == 0) && (Q0.b == 0)) {
				zi = conv<ZZ>(0);
			} else {
				E.initHashFunction(P0, Q0, factors[i]);
				zi = E.pollardAttack(P0, Q0, factors[i]);
			}
			cout << zi << endl;
			ls.push_back(zi);
			ps.push_back(factors[i]);
		} else {
			ps[ps.size() - 1] *= factors[i];
			c1 /= factors[i];
			c2 *= factors[i];
			Pol = E.add(Pol, E.negate(E.mult(zi, c2p)));
			Affine Qi = E.mult(c1, Pol);
			if (P0 == Qi) {
				zi = conv<ZZ>(1);
			} else if ((Qi.a == 0) && (Qi.b == 0)) {
				zi = conv<ZZ>(0);
			} else {
				E.initHashFunction(P0, Qi, factors[i]);
				zi = E.pollardAttack(P0, Qi, factors[i]);
			}
			cout << zi << endl;
			c2p = E.mult(factors[i], c2p);
			ls[ls.size() - 1] += zi * c2;
		}
	}
}

int main() {
	while (1) {
#ifdef DEBUG
		p = conv<ZZ>(229);
		cout << p << endl;
#else
	    cin >> p;
		if (p == 0) {
			break;
		}
		SetSeed(p);
#endif
/*
		ZZ powers = p;

		for (int k=0; k<6; k++) {
			powers *= p;
			cout << powers << endl;
		}
		return 0;
		*/
	    ZZ_p::init(p);

		ZZ_p a, b;
		ZZ_p x1, y1, x2, y2;
#ifdef DEBUG
		a = conv<ZZ_p>(1);
		b = conv<ZZ_p>(44);
#else
		cin >> a >> b;
#endif
		
#ifdef DEBUG
		x1 = conv<ZZ_p>(5);
		y1 = conv<ZZ_p>(116);
		x2 = conv<ZZ_p>(155);
		y2 = conv<ZZ_p>(166);
#else
		cin >> x1 >> y1;
		cin >> x2 >> y2;
#endif
		if ((a == 0) && (b == 0)) {
			ZZ_p c1 = y1 * y1 - x1 * x1 * x1;
			ZZ_p c2 = y2 * y2 - x2 * x2 * x2;
			a = (c1 - c2) / (x1 - x2);
			b = c1 - a * x1;
			cout << a << " " << b << endl;
		}

		ECCprime E(a, b, p);

		ZZ n;

#ifdef DEBUG
		n = conv<ZZ>(239);
#else
		cin >> n;
#endif
		
		Affine P(x1, y1 /*, conv<ZZ_p>(1)*/), Q(x2, y2/*, conv<ZZ_p>(1)*/);
		
		if (n == 0) {
			n = findOrder(p, E, P);
			cout << n << endl;
			vector<ZZ> ls;
			vector<ZZ> ps;
			findLP(n, ls, ps, E, P, Q);
			int ss = ls.size();
			for (int i = 0; i < ss; ++i) {
				cout << ls[i] << " " << ps[i] << endl;
			}
			cout << chineseRemainder(ls, ps, ps.size()) << endl;
			continue;
		}
	
		SetSeed(p);
		E.initHashFunction(P, Q, n);			
		cout << "Starting Attack  " << endl;
		//SetSeed(mySeed);
		ZZ l = E.pollardAttack(P, Q, n);

		cout << "Private key = " << l << endl;
		cout << "Attack Ended " << endl;
		cout << "Q = " << Q << endl;
		cout << "lP = " << E.mult(l, P) << endl;
		//mySeed += 456;
		//mySeed *= 6364136223846793005;
		//mySeed += 1442695040888963407;
	}
	system("pause");
	return 0;
}