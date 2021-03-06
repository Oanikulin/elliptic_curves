#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
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

struct Remainder {
    unsigned long long val, p;

    Remainder(unsigned long long a) {
        val = 1;
        p = a;
    }

    Remainder(unsigned long long a, unsigned long long b) {
        val = a;
        p = b;
    }

    Remainder operator * (Remainder other) {
        if (p != other.p)
            throw "kek";
        return Remainder((val * 1ULL * other.val) % p, p);
    }

    Remainder operator * (unsigned long long other) {
        return Remainder((val * 1ULL * other) % p, p);
    }
};

template<typename T, typename T2>
T pow(T val, T2 deg) {
    T res(val.p);
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

char number_to_char(int num) {
    if (num < 10)
        return '0' + num;
    if (num < 36)
        return 'A' + num - 10;
    if (num < 36 + 26)
        return 'a' + num - 36;
    if (num == 62)
        return 32;
    if (num == 63)
        return 46;
    assert(false);
}

vector<pair<unsigned int, unsigned int>> ans;

template<typename T>
void encrypt(T num, T p, T g, T k) {
    ans.push_back({num * 1ULL * k % p, g});
}

template <typename T>
T invert(T a, unsigned int mod) {
    return pow(a, mod - 2);
}

template<typename T>
T decrypt(T a, T b, unsigned int p, unsigned int k) {
    return a * invert(pow(b, k), p);
}

signed main() {
    unsigned long long p, g, k;
    cin >> p >> g >> k;
    string s;
    getline(cin, s);
    getline(cin, s);
    reverse(s.begin(), s.end());
    bigint a(0);
    for (char c : s) {
        a = a * (long long)64 + (unsigned long long)char_to_number(c);
    }
    while (!a.isZero()) {
        unsigned long long cur = a % (unsigned long long)(p);
        encrypt(cur, p, g, k);
        a /= (unsigned long long)p;
    }
    for (auto k : ans) {
        cout << k.second << " " << k.first << endl;
    }
}
