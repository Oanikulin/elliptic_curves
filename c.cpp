#include <bits/stdc++.h>
using namespace std;

template <typename T>
struct Polynomial {
    vector<T> polynom;
    void normal() {
        while (polynom.size() && polynom.back() == T(0))
            polynom.pop_back();
    }
    Polynomial( vector<T>& k) {
        for (auto val : k)
            polynom.push_back(val);
        normal();
    }
    Polynomial(T k = T()) {
        polynom.push_back(k);
        normal();
    }
    template<typename Iterator,
            typename ValueType = typename std::iterator_traits<Iterator>::value_type>
    Polynomial(Iterator beg, Iterator en) {
        for (; beg != en; beg = next(beg)) {
            polynom.push_back(*beg);
        }
        normal();
    }
    T operator[] (size_t it)  {
        if (it >= polynom.size())
            return T(0);
        return polynom[it];
    }
    size_t size()  {
        return polynom.size();
    }

    int Degree()  {
        return static_cast<int>(polynom.size()) - 1;
    }
    T back()  {
        return polynom.back();
    }
    T operator() (T pt) {
        T ans = T(0);
        if (polynom.size())
            ans += polynom[0];
        T pow = pt;
        for (size_t i = 1; i < polynom.size(); ++i) {
            ans += polynom[i] * pow;
            pow *= pt;
        }
        return ans;
    }
    Polynomial<T> operator+ (Polynomial<T> sec) {
        vector<T> tmp(std::max(size(), sec.size()));
        for (size_t i = 0; i < std::max(size(), sec.size()); ++i) {
            if (i < size())
                tmp[i] += polynom[i];
            if (i < sec.size())
                tmp[i] += sec[i];
        }
        return Polynomial<T>(tmp);
    }
    Polynomial<T> operator- (Polynomial<T> sec) {
        vector<T> tmp(std::max(size(), sec.size()));
        for (size_t i = 0; i < std::max(size(), sec.size()); ++i) {
            if (i < size())
                tmp[i] = tmp[i] + polynom[i];
            if (i < sec.size())
                tmp[i] = tmp[i] - sec[i];
        }
        return Polynomial<T>(tmp);
    }
    Polynomial<T> operator* (Polynomial<T> sec) {
        vector<T> tmp(size() + sec.size());
        for (size_t i = 0; i < size(); ++i)
            for (size_t j = 0; j < sec.size(); ++j)
                tmp[i + j] = tmp[i + j] + polynom[i] * sec[j];
        return Polynomial(tmp);
    }
    Polynomial<T>& operator *= (Polynomial& other) {
        (*this) = (*this) * other;
        return (*this);
    }
    Polynomial<T>& operator += (Polynomial& other) {
        (*this) = (*this) + other;
        return (*this);
    }
    Polynomial<T>& operator -= (Polynomial& other) {
        (*this) = (*this) - other;
        return (*this);
    }
    Polynomial<T>& operator += (T other) {
        (*this) = (*this) + Polynomial(other);
        return (*this);
    }
    Polynomial<T> operator -= (T other) {
        (*this) = (*this) - Polynomial(other);
        return (*this);
    }
    Polynomial<T>& operator *= (T other) {
        (*this) = (*this) * Polynomial(other);
        return (*this);
    }
    void swap(Polynomial& other) {
        std::swap(polynom, other.polynom);
    }
};

template<typename T>
bool operator== ( Polynomial<T>& fr,  Polynomial<T>& sec) {
    if (fr.size() != sec.size())
        return 0;
    for (size_t i = 0; i < fr.size(); ++i)
        if (fr[i] != sec[i])
            return 0;
    return 1;
}

template<typename T>
bool operator== ( Polynomial<T>& fr, T sec) {
    return (fr == Polynomial<T>(sec));
}

template<typename T>
bool operator== (T sec,  Polynomial<T>& fr) {
    return (fr == Polynomial<T>(sec));
}

template<typename T>
bool operator!= ( Polynomial<T>& fr,  Polynomial<T>& sec) {
    return !(fr == sec);
}

template<typename T>
bool operator!= ( Polynomial<T>& fr, T sec) {
    return !(fr == Polynomial<T>(sec));
}

template<typename T>
bool operator!= (T sec,  Polynomial<T>& fr) {
    return !(fr == Polynomial<T>(sec));
}

template<typename T>
Polynomial<T> operator+ ( Polynomial<T>& fr, T sec) {
    return fr + Polynomial<T>(sec);
}

template<typename T>
Polynomial<T> operator+ (T sec,  Polynomial<T>& fr) {
    return fr + Polynomial<T>(sec);
}

template<typename T>
Polynomial<T> operator- ( Polynomial<T>& fr, T sec) {
    return fr - Polynomial<T>(sec);
}

template<typename T>
Polynomial<T> operator- (T sec,  Polynomial<T>& fr) {
    return Polynomial<T>(sec) - fr;
}

template<typename T>
Polynomial<T> operator* ( Polynomial<T>& fr, T sec) {
    return fr * Polynomial<T>(sec);
}

template<typename T>
Polynomial<T> operator* (T sec,  Polynomial<T>& fr) {
    return Polynomial<T>(sec) * fr;
}

template<typename T>
std::pair<Polynomial<T>, Polynomial<T>> get_div(Polynomial<T> fir,
                                                Polynomial<T> sec) {
    if (fir.size() < sec.size())
        return {Polynomial<T>(T(0)), fir};
    vector <T> coef(fir.size());
    Polynomial<T> tmp = fir;
    for (int i = fir.Degree(); i >= sec.Degree(); --i) {
        if (tmp[i] == T(0))
            continue;
        vector<T> pw(tmp.size());
        pw[tmp.Degree() - sec.size() + 1] = tmp[i] / sec[sec.size() - 1];
        coef[tmp.Degree() - sec.size() + 1] = tmp[i] / sec[sec.size() - 1];
        tmp = tmp - (sec * Polynomial<T>(pw));
    }
    return {Polynomial<T>(coef), tmp};
}

template<typename T>
Polynomial<T> operator/(Polynomial<T> first, Polynomial<T> sec) {
    return get_div(first, sec).first;
}

template<typename T>
Polynomial<T> operator%(Polynomial<T> first,  Polynomial<T> sec) {
    return get_div(first, sec).second;
}

const int maxn = 1e4 + 14, lg = 15;

/*
  ######################################################################
  #######################   THE   BIG   INT   ##########################
*/
const int base = 100000000;
const int base_digits = 8;
struct bigint {
	vector<long long> a;
	int sign;
	/*<arpa>*/
	int size(){
		if(a.empty())return 0;
		int ans=(a.size()-1)*base_digits;
		int ca=a.back();
		while(ca)
			ans++,ca/=10;
		return ans;
	}
	/*</arpa>*/
	bigint() :
		sign(1) {
	}

	bigint(long long v) {
		*this = v;
	}

	void operator=(const bigint &v) {
		sign = v.sign;
		a = v.a;
	}

	void operator=(long long v) {
		sign = 1;
		a.clear();
		if (v < 0)
			sign = -1, v = -v;
		for (; v > 0; v = v / base)
			a.push_back(v % base);
	}

	bigint operator+(const bigint &v) const {
		if (sign == v.sign) {
			bigint res = v;

			for (int i = 0, carry = 0; i < (int) max(a.size(), v.a.size()) || carry; ++i) {
				if (i == (int) res.a.size())
					res.a.push_back(0);
				res.a[i] += carry + (i < (int) a.size() ? a[i] : 0);
				carry = res.a[i] >= base;
				if (carry)
					res.a[i] -= base;
			}
			return res;
		}
		return *this - (-v);
	}

    bigint abs() const {
		bigint res = *this;
		res.sign *= res.sign;
		return res;
	}


	bigint operator-(const bigint &v) const {
		if (sign == v.sign) {
			if (abs() >= v.abs()) {
				bigint res = *this;
				for (int i = 0, carry = 0; i < (int) v.a.size() || carry; ++i) {
					res.a[i] -= carry + (i < (int) v.a.size() ? v.a[i] : 0);
					carry = res.a[i] < 0;
					if (carry)
						res.a[i] += base;
				}
				res.trim();
				return res;
			}
			return -(v - *this);
		}
		return *this + (-v);
	}

	void operator*=(int v) {
		if (v < 0)
			sign = -sign, v = -v;
		for (int i = 0, carry = 0; i < (int) a.size() || carry; ++i) {
			if (i == (int) a.size())
				a.push_back(0);
			long long cur = a[i] * (long long) v + carry;
			carry = (int) (cur / base);
			a[i] = (int) (cur % base);
			//asm("divl %%ecx" : "=a"(carry), "=d"(a[i]) : "A"(cur), "c"(base));
		}
		trim();
	}

	bigint operator*(int v) const {
		bigint res = *this;
		res *= v;
		return res;
	}

	void operator*=(long long v) {
		if (v < 0)
			sign = -sign, v = -v;
		if(v > base){
			*this = *this * (v / base) * base + *this * (v % base);
			return ;
		}
		for (int i = 0, carry = 0; i < (int) a.size() || carry; ++i) {
			if (i == (int) a.size())
				a.push_back(0);
			long long cur = a[i] * (long long) v + carry;
			carry = (int) (cur / base);
			a[i] = (int) (cur % base);
			//asm("divl %%ecx" : "=a"(carry), "=d"(a[i]) : "A"(cur), "c"(base));
		}
		trim();
	}

	bigint operator*(long long v) const {
		bigint res = *this;
		res *= v;
		return res;
	}

	friend pair<bigint, bigint> divmod(const bigint &a1, const bigint &b1) {
		int norm = base / (b1.a.back() + 1);
		bigint a = a1.abs() * norm;
		bigint b = b1.abs() * norm;
		bigint q, r;
		q.a.resize(a.a.size());

		for (int i = a.a.size() - 1; i >= 0; i--) {
			r *= base;
			r += a.a[i];
			int s1 = r.a.size() <= b.a.size() ? 0 : r.a[b.a.size()];
			int s2 = r.a.size() <= b.a.size() - 1 ? 0 : r.a[b.a.size() - 1];
			int d = ((long long) base * s1 + s2) / b.a.back();
			r -= b * d;
			while (r < 0)
				r += b, --d;
			q.a[i] = d;
		}

		q.sign = a1.sign * b1.sign;
		r.sign = a1.sign;
		q.trim();
		r.trim();
		return make_pair(q, r / norm);
	}

	bigint operator/(const bigint &v) const {
		return divmod(*this, v).first;
	}

	bigint operator%(const bigint &v) const {
		return divmod(*this, v).second;
	}

	void operator/=(int v) {
		if (v < 0)
			sign = -sign, v = -v;
		for (long long i = (int) a.size() - 1, rem = 0; i >= 0; --i) {
			long long cur = a[i] + rem * (long long) base;
			a[i] = (int) (cur / v);
			rem = (int) (cur % v);
		}
		trim();
	}

    void operator/=(unsigned long long v) {
		if (v < 0)
			sign = -sign, v = -v;
		for (long long i = (int) a.size() - 1, rem = 0; i >= 0; --i) {
			unsigned long long cur = a[i] + rem * (unsigned long long) base;
			a[i] = (long long) (cur / v);
			rem = (unsigned long long) (cur % v);
		}
		trim();
	}

	bigint operator/(int v) const {
		bigint res = *this;
		res /= v;
		return res;
	}


	bigint operator-() const {
		bigint res = *this;
		res.sign = -sign;
		return res;
	}

	int operator%(int v) const {
		if (v < 0)
			v = -v;
		int m = 0;
		for (int i = a.size() - 1; i >= 0; --i)
			m = (a[i] + m * (long long) base) % v;
		return m * sign;
	}

    unsigned long long operator%(unsigned long long v) const {
		unsigned long long m = 0;
		for (int i = a.size() - 1; i >= 0; --i)
			m = (a[i] + m * (long long) base) % v;
		return m * sign;
	}

	void operator+=(const bigint &v) {
		*this = *this + v;
	}
	void operator-=(const bigint &v) {
		*this = *this - v;
	}
	void operator/=(const bigint &v) {
		*this = *this / v;
	}

	bool operator<(const bigint &v) const {
		if (sign != v.sign)
			return sign < v.sign;
		if (a.size() != v.a.size())
			return a.size() * sign < v.a.size() * v.sign;
		for (int i = a.size() - 1; i >= 0; i--)
			if (a[i] != v.a[i])
				return a[i] * sign < v.a[i] * sign;
		return false;
	}

	bool operator>(const bigint &v) const {
		return v < *this;
	}
	bool operator<=(const bigint &v) const {
		return !(v < *this);
	}
	bool operator>=(const bigint &v) const {
		return !(*this < v);
	}
	bool operator==(const bigint &v) const {
		return !(*this < v) && !(v < *this);
	}
	bool operator!=(const bigint &v) const {
		return *this < v || v < *this;
	}

	void trim() {
		while (!a.empty() && !a.back())
			a.pop_back();
		if (a.empty())
			sign = 1;
	}

	bool isZero() const {
		return a.empty() || (a.size() == 1 && !a[0]);
	}

};



template<typename T>
int sgn(T a) {
    if (a > 0)
        return 1;
    if (a < 0)
        return -1;
    return 0;
}

unsigned long long MOD;

template<typename T, typename T2>
T rem_pow(T val, T2 deg) {
    T res(1);
    while (deg) {
        if (deg & 1) {
            --deg;
            res = res * val;
        }
        val = val * val;
        deg /= 2;
    }
    return res;
}

struct Remainder{
    unsigned long long val, p;

    Remainder() {
        val = 0;
        p = MOD;
    }

    Remainder(long long a) {
        val = a % MOD;
        p = MOD;
    }

    Remainder(long long a, unsigned int b) {
        val = a % b;
        p = MOD;
    }

    Remainder operator * (Remainder other) {
        if (p != other.p)
            throw "kek";
        return Remainder((val * 1ULL * other.val) % p, p);
    }

    Remainder operator * (long long other) {
        return Remainder((val * 1ULL * other) % p, p);
    }

    Remainder operator + (Remainder other) {
        return Remainder((val * 1ULL + other.val) % p, p);
    }

    Remainder operator - (Remainder other) {
        return Remainder((val * 1ULL - other.val + p) % p, p);
    }

    Remainder operator - (long long other) {
        return Remainder((val * 1ULL - other + p) % p, p);
    }

    Remainder operator / (Remainder other) {
        //cout << "division " << p << " "<< val << " " << other.val << " " << rem_pow(other, p - 2).val << endl;
        return (*this) * rem_pow(other, p - 2);
    }

    Remainder operator % (Remainder other) {
        return Remainder(val - ((*this)/other).val * other.val + p, p);
    }

    bool operator == (Remainder other) {
        return val == other.val;
    }

    bool operator != (Remainder other) {
        return val != other.val;
    }
};

template<typename T>
std::ostream& operator<<(std::ostream& out,  Polynomial<T> m) {
    m.normal();
    int deg = m.Degree();
    if (deg == -1) {
        return out;
    }
    for (int i = 0; i <= m.Degree(); ++i) {
        out << m[i].val << " ";
    }
    return out;
}

int char_to_number(char symbol) {
  if (symbol >= 48 && symbol <= 57)
    return symbol - 48;
  if (symbol >= 65 && symbol <= 90)
    return symbol - 55;
  if (symbol >= 97 && symbol <= 122)
    return symbol - 61;
  if (symbol == 32)
    return 62;
  if (symbol == 46)
    return 63;
  assert(false);
}

vector<pair<Polynomial<Remainder>, Polynomial<Remainder>>> ans;

template<typename T, typename T2>
T pow(T val, T2 deg, T md) {
    T res(1);
    while (deg) {
        if (deg & 1) {
            --deg;
            res = res * val % md;
        }
        val = val * val % md;
        deg /= 2;
    }
    return res;
}



template<typename T>
void encrypt(T num, T p, T g, T k, long long mod) {
    int b = rand();
    ans.push_back({num * pow(k, b,p) % p , pow(g, b,p)});
}

vector<Remainder> parse_string(string &s, long long p) {
    vector<Remainder> tmp;
    stringstream ss;
    ss << s;
    long long cur;
    while (ss >> cur) {
        if (cur < 0)
            cur += p;
        tmp.push_back(Remainder(cur, p));
    }
    return tmp;
}

int main() {
    srand(time(0));
    unsigned long long p;
    cin >> p;
    MOD = p;
    string s;
    getline(cin, s);
    getline(cin, s);
    vector<Remainder> fr = parse_string(s, p);
    getline(cin, s);
    vector<Remainder> sec = parse_string(s, p);
    getline(cin, s);
    vector<Remainder> thr = parse_string(s, p);
    Polynomial<Remainder> f(fr), g(sec), k(thr);
    getline(cin, s);
    reverse(s.begin(), s.end());
    bigint a;
    int deg = f.Degree();
    for (char c : s) {
        a = a * 64 + char_to_number(c);
    }
    int cnt = 0;
    vector<Remainder> pol;
    while (!a.isZero()) {
        unsigned int cur = a % (unsigned long long)(p);
        a /= p;
        ++cnt;
        pol.push_back(Remainder(cur, p));
        if (cnt == deg) {
            encrypt(Polynomial<Remainder>(pol), f, g, k, p);
            cnt = 0;
            pol.resize(0);
        }
    }

    if (cnt) {
        encrypt(Polynomial<Remainder>(pol), f, g, k, p);
    }
    for (auto k : ans) {
        cout << k.second << endl << k.first << endl;
    }
}

