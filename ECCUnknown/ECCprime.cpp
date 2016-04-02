#include <NTL/ZZ_pE.h>
#include <NTL/pair.h>
#include <NTL/ZZ_p.h>

using namespace NTL;
using namespace std;

#include <fstream>
#include <windows.h>

ZZ p;

inline ZZ_p square(const ZZ_p & x) {
	return x * x;
}

struct Jacobian {
	ZZ_p x, y, z, z2;
	Jacobian(const ZZ_p& X, const ZZ_p& Y, const ZZ_p& Z) : x(X), y(Y), z(Z) {
		z2 = square(z);
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

typedef Pair<ZZ_p, ZZ_p> Affine;


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

#define __WORDSIZE 64

class optimize;

class optimize { 
private:
	ZZ n,k,m,p; // n is the number of bits of given number, m equals some integer i times w (wordsize) s.t. m=i*w >= n.
	ZZ a,r; // the number which we are looking for is inverse.
public:
	optimize(ZZ_p aa, ZZ pp) : a(conv<ZZ>(aa)), p(pp) {
		n = NumBits(p);
		m=__WORDSIZE;
		while(m<n) {
			m+=m; 
		}
		// cout<< "passed the c'tor" << endl;
	}
	void AlmMonInv() {
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
		if(r>=p) {
			
			r = r-p;
		}
		//cout << "r1 = " << r << endl;
		r = p-r;
		//cout << "phase1 ended r= " << r << endl;
	}
	ZZ_p Inv() {
		if(k>m) {
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

class triple;

#include <NTL/ZZ_pXFactoring.h>
#include <NTL/ZZ_pEX.h>

#ifdef DEBUG
#define NUM_BRANCHES 4
#else
#define NUM_BRANCHES 32
#endif

#define BRANCH_MASK (NUM_BRANCHES-1)

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

int inline partition_function(const Jacobian& P) {
	//return (rep(P.x) % NUM_BRANCHES);
	ZZ_p a = P.x;
	if (P.z != 0) {
		optimize I(P.z2, p);
		I.AlmMonInv();
		a *= I.Inv();
	}

	return conv<int>(rep(a)) & BRANCH_MASK;
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
		ZZ_p H = U2 - P1.x;
		ZZ_p R = S2 - P1.y;
		ZZ_p h2 = square(H);
		ZZ_p h3 = h2 * H;
		ZZ_p tmp = P1.x*h2;
		ZZ_p X3 = square(R) - h3 - 2*tmp;

		return Jacobian(X3, R*(tmp - X3) - P1.y*h3, H*P1.z);
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
	Jacobian negate(const Jacobian& P) {
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
		ZZ_p X = square((P2.b - P1.b) / (P2.a - P1.a)) - P1.a - P2.a;
		return Affine(X, ((P2.b - P1.b) / (P2.a - P1.a)) * (P1.a - X) - P1.b);
	}
	Affine doubling(const Affine& P) {
		ZZ_p tmp = (3 * square(P.a) + a) / (2 * P.b);
		ZZ_p X = square(tmp) - 2 * P.a;
		return Affine(X, tmp * (P.a - X) - P.b);
	}
	Affine negate(const Affine& P) {
		return Affine(P.a, -P.b);
	}
	Affine mult(ZZ c, Affine P) {
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

	ZZ pollard_attack(Jacobian P, Jacobian Q, ZZ order_p) {
		ZZ ctag, dtag, ctags, dtags;
		Jacobian Xtag, Xtags;

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
		/*
		 power = lam = 1
    tortoise = x0
    hare = f(x0)  # f(x0) is the element/node next to x0.
    while tortoise != hare:
        if power == lam:  # time to start a new power of two?
            tortoise = hare
            power *= 2
            lam = 0
        hare = f(hare)
        lam += 1
		*/
		long long power = 1, lamda = 1;
		index = partition_function(Xtags);
		Xtags = add(Xtags, (triples_array[index]).getRj());
		ctags = (ctags + (triples_array[index]).getAj())%order_p;
		dtags = (dtags + (triples_array[index]).getBj())%order_p;

		while (Xtag!=Xtags) {
			if (power == lamda) {
				Xtag = Xtags;
				ctag = ctags;
				dtag = dtags;
				power *= 2;
				lamda = 0;
			}

			index = partition_function(Xtags);
			/*cout << "index= " << index << " " << "R= " << (triples_array[index]).getRj() << endl;
			cout << "R.z= " << (triples_array[index]).getRj().z << endl;
			cout << "x''= " << Xtags << "x''.z= " << Xtags.z << endl;*/
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
			if (iii == 100000000) {
				cout << "iteration count: " << ccc << endl;
				iii = 0;
			}
		};
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

	void init_iteration_function(Jacobian P, Jacobian Q, ZZ order_p) {
	
#ifdef DEBUG
		long as[NUM_BRANCHES] = {79, 206, 87, 219}; 
		long bs[NUM_BRANCHES] = {163, 19, 109, 68};
#endif

		ZZ Aj,Bj; 
		Jacobian Rj;
		for(int j=0; j<NUM_BRANCHES; j++) {
#ifdef DEBUG
			Aj=conv<ZZ>(as[j]); //RandomBnd(order_p);
			Bj=conv<ZZ>(bs[j]); //RandomBnd(order_p);
#else
			do {
				Aj=RandomBnd(order_p);
				Bj=RandomBnd(order_p);
				Rj = add(mult(Aj,P),mult(Bj,Q));
			} while (Rj.z == 0);
#endif
			triples_array[j] = triple(Aj,Bj,Affine(Rj.x / Rj.z2, Rj.y / (Rj.z2 * Rj.z)));
		}
	}

private:
	ZZ_p a, b;
	ZZ p;
	triple triples_array[NUM_BRANCHES];
};

#include <vector>

using namespace std;

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

#ifdef DEBUG
		x1 = conv<ZZ_p>(5);
		y1 = conv<ZZ_p>(116);
		x2 = conv<ZZ_p>(155);
		y2 = conv<ZZ_p>(166);
#else
		cin >> x1 >> y1;
		cin >> x2 >> y2;
#endif
		ECCprime E(a, b, p);
#ifdef MON
		// examlpe
		optimize A(conv<ZZ_p>(54), p); // initilliaze a variable of the class optimize with the number which you are looking for is inverse - a
		// and the group order - p
		A.AlmMonInv(); // run those two methods on the variable you initillized and then the second one will return the inverse.
		ZZ_p inverse = A.Inv();
		// end of example
		cout << inverse << endl;
		cout << inverse*conv<ZZ_p>(54) << endl;
#endif
		if ((a==0) && (b==0)) {
			ZZ_p c1 = y1 * y1 - x1 * x1 * x1;
			ZZ_p c2 = y2 * y2 - x2 * x2 * x2;
			a = (c1 - c2) / (x1 - x2);
			b = c1 - a * x1;
			cout << a << endl << b << endl;
		}
		ZZ n;
		Jacobian P(x1, y1, conv<ZZ_p>(1)), Q(x2, y2, conv<ZZ_p>(1));

#ifdef DEBUG
		n = conv<ZZ>(239);
#else
		cin >> n;
		if (n == 0) {
			n = findOrder(p, E, P);
			cout << n << endl;
			system("pause");
			return 0;
		}
#endif
		//Jacobian P(x1, y1, conv<ZZ_p>(1)), Q(x2, y2, conv<ZZ_p>(1));
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

	//	ZZ mySeed = conv<ZZ>(0x44e458c2);
		SetSeed(p);
		E.init_iteration_function(P, Q, n);

		//ZZ mySeed = conv<ZZ>(123);
		for (int i=0; i<10; i++) {
			cout << "Starting Attack  " <<endl;
			//SetSeed(mySeed);
			//E.init_iteration_function(P, Q, n);
			ZZ l = E.pollard_attack(P, Q, n);
			cout << "Private key = " << l << endl;
			cout << "Attack Ended " << endl;
			cout << "Q = " << Q << endl;
			cout << "lP = " << E.mult(l, P) << endl;
			//mySeed += 456;
			//mySeed *= 6364136223846793005;
			//mySeed += 1442695040888963407;
		}
	}
	system("pause");
	return 0;
}