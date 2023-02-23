#include "int2048.h"
#include <complex>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <string>

namespace sjtu {
constexpr long long BASE = 1000000000ll;
const double pi = acos(-1);

typedef std::complex<double> cp;

std::string removezero(std::string s) {
  // 去除字符串前导 0
  std::string ret;
  int siz = s.size();
  int cur = 0;
  if (s[0] == '-') {
    ret.push_back('-');
    ++cur;
  }
  for (; cur < siz && s[cur] == '0'; ++cur) {
  }
  if (cur == siz) {
    ret = "0";
    return ret;
  }
  ret.append(s.substr(cur, siz - cur));
  if (ret.size() == 2 && ret[0] == '-' && ret[1] == '0') {
    ret = "0";
  }
  return ret;
}

int2048::int2048() {
  vec.clear();
  sig = false;
  vec.push_back(0);
}
int2048::int2048(long long x) {
  vec.clear();
  if (x == 0) {
    sig = false;
    vec.push_back(0);
    return;
  }
  if (x < 0) {
    sig = true;
    x = -x;
  } else {
    sig = false;
  }
  while (x > 0) {
    vec.push_back(x % BASE);
    x /= BASE;
  }
}
int2048::int2048(const std::string &s) {
  vec.clear();
  std::string tmp = removezero(s);
  int len = tmp.size();
  if (len == 0) {
    sig = false;
    vec.push_back(0);
    return;
  }
  if (tmp[0] == '-') {
    sig = true;
    tmp = tmp.substr(1, len - 1);
    len = tmp.size();
  } else {
    sig = false;
  }
  int i = len - 1;
  for (; i >= 9; i -= 9) {
    long long now = 0;
    for (int j = i - 8; j <= i; ++j) {
      now = now * 10ll + tmp[j] - '0';
    }
    vec.push_back(now);
  }
  vec.push_back(0);
  int cur = vec.size() - 1;
  for (int j = 0; j <= i; ++j) {
    vec[cur] = vec[cur] * 10ll + tmp[j] - '0';
  }
}
int2048::int2048(const int2048 &x) {
  sig = x.sig;
  vec.clear();
  vec.resize(x.vec.size());
  for (long unsigned int i = 0; i < x.vec.size(); ++i) {
    vec[i] = x.vec[i];
  }
}
int2048::int2048(const int2048 &&x) {
  if (&x == this) return;
  sig = x.sig;
  vec.clear();
  vec.resize(x.vec.size());
  for (long unsigned int i = 0; i < x.vec.size(); ++i) {
    vec[i] = x.vec[i];
  }
}

void cut(const int2048 &x, int2048 &y, int low, int high) {
  // 从 low 到 high 截断
  y.vec.clear();
  y.sig = x.sig;
  if (low > high) {
    y.vec.push_back(0);
    return;
  }
  y.vec.resize(high - low + 1);
  for (int i = low; i <= high; ++i) {
    y.vec[i - low] = x.vec[i];
  }
}

bool unsigned_cmp(const int2048 &x, const int2048 &y) {
  // 大于等于号（绝对值），用于无符号减法
  if (x.vec.size() != y.vec.size()) {
    return (x.vec.size() > y.vec.size());
  }
  int siz = x.vec.size();
  for (int i = siz - 1; i >= 0; --i) {
    if (x.vec[i] != y.vec[i]) {
      return (x.vec[i] > y.vec[i]);
    }
  }
  return true;  // x == y
}

void int2048::read(const std::string &s) {
  vec.clear();
  sig = false;
  std::string tmp = removezero(s);
  int len = tmp.size();
  if (len == 0) {
    sig = false;
    vec.push_back(0);
    return;
  }
  if (tmp[0] == '-') {
    sig = true;
    tmp = tmp.substr(1, len - 1);
    len = tmp.size();
  } else {
    sig = false;
  }
  int i = len - 1;
  for (; i >= 9; i -= 9) {
    long long now = 0;
    for (int j = i - 8; j <= i; ++j) {
      now = now * 10ll + tmp[j] - '0';
    }
    vec.push_back(now);
  }
  vec.push_back(0);
  int cur = vec.size() - 1;
  for (int j = 0; j <= i; ++j) {
    vec[cur] = vec[cur] * 10ll + tmp[j] - '0';
  }
}
void int2048::print() {
  std::cout << (*this);
}

int2048 &int2048::unsigned_add(const int2048 &y) {
  // 无符号加法
  if (vec.size() < y.vec.size()) {
    vec.resize(y.vec.size());
  }
  int siz = vec.size(), carry = 0;
  for (int i = 0; i < siz; ++i) {
    if (i >= (int)y.vec.size()) {
      vec[i] = vec[i] + carry;
      if (vec[i] > BASE) {
        vec[i] -= BASE;
        carry = 1;
      } else {
        carry = 0;
      }
      continue;
    }
    vec[i] = vec[i] + y.vec[i] + carry;
    if (vec[i] >= BASE) {
      vec[i] -= BASE;
      carry = 1;
    } else {
      carry = 0;
    }
  }
  if (carry) {
    vec.push_back(1ll);
  }
  return *this;
}

int2048 unsigned_add(int2048 x, const int2048 &y) {
  if (x.vec.size() < y.vec.size()) {
    x.vec.resize(y.vec.size());
  }
  int siz = x.vec.size(), carry = 0;
  for (int i = 0; i < siz; ++i) {
    if (i >= (int)y.vec.size()) {
      x.vec[i] = x.vec[i] + carry;
      if (x.vec[i] > BASE) {
        x.vec[i] -= BASE;
        carry = 1;
      } else {
        carry = 0;
      }
      continue;
    }
    x.vec[i] = x.vec[i] + y.vec[i] + carry;
    if (x.vec[i] >= BASE) {
      x.vec[i] -= BASE;
      carry = 1;
    } else {
      carry = 0;
    }
  }
  if (carry) {
    x.vec.push_back(1ll);
  }
  return x;
}

int2048 &int2048::unsigned_minus(const int2048 &y) {
  // 无符号减法，要求大数减小数
  int siz = vec.size(), carry = 0;
  for (int i = 0; i < siz; ++i) {
    if (i >= (int)y.vec.size()) {
      vec[i] = vec[i] - carry;
      if (vec[i] < 0) {
        // carry
        carry = 1;
        vec[i] += BASE;
      } else {
        carry = 0;
      }
      continue;
    }
    vec[i] = vec[i] - y.vec[i] - carry;
    if (vec[i] < 0) {
      // carry
      carry = 1;
      vec[i] += BASE;
    } else {
      carry = 0;
    }
  }
  for (int i = siz - 1; i >= 0; --i) {
    if (vec[i]) {
      break;
    }
    vec.pop_back();
  }
  if (vec.empty()) {
    vec.push_back(0ll);
  }
  return *this;
}
int2048 unsigned_minus(int2048 x, const int2048 &y) {
  int siz = x.vec.size(), carry = 0;
  for (int i = 0; i < siz; ++i) {
    if (i >= (int)y.vec.size()) {
      x.vec[i] = x.vec[i] - carry;
      if (x.vec[i] < 0) {
        // carry
        carry = 1;
        x.vec[i] += BASE;
      } else {
        carry = 0;
      }
      continue;
    }
    x.vec[i] = x.vec[i] - y.vec[i] - carry;
    if (x.vec[i] < 0) {
      // carry
      carry = 1;
      x.vec[i] += BASE;
    } else {
      carry = 0;
    }
  }
  for (int i = siz - 1; i >= 0; --i) {
    if (x.vec[i]) {
      break;
    }
    x.vec.pop_back();
  }
  if (x.vec.empty()) {
    x.vec.push_back(0ll);
  }
  return x;
}

int2048 &int2048::cmp_minus(const int2048 &y) {
  // 真正的无符号减法，|x| - |y|
  bool flag = unsigned_cmp(*this, y);
  if (flag) {
    this->unsigned_minus(y);
    sig = false;
  } else {
    int2048 ret = sjtu::unsigned_minus(y, *this);
    vec.clear();
    vec.resize(ret.vec.size());
    for (unsigned i = 0; i < ret.vec.size(); ++i) {
      vec[i] = ret.vec[i];
    }
    sig = true;
  }
  return (*this);
}
int2048 cmp_minus(int2048 x, const int2048 &y) {
  bool flag = unsigned_cmp(x, y);
  if (flag) {
    x.unsigned_minus(y);
    x.sig = false;
    return x;
  } else {
    x = sjtu::unsigned_minus(y, x);
    x.sig = true;
    return x;
  }
}

int2048 &int2048::add(const int2048 &y) {
  // 有符号加法
  if (sig) {
    if (y.sig) {
      this->unsigned_add(y);
      this->sig = true;
    } else {
      // -|x| + |y|
      // |y| - |x| = -(|x| - |y|)
      this->cmp_minus(y);
      if (vec.size() == 1 && !vec[0]) {
        sig = false;
      } else {
        this->sig ^= 1;
      }
    }
  } else {
    if (y.sig) {
      // |x| - |y|
      this->cmp_minus(y);
    } else {
      this->unsigned_add(y);
    }
  }
  return (*this);
}
int2048 add(int2048 x, const int2048 &y) {
  if (x.sig) {
    if (y.sig) {
      x.unsigned_add(y);
      x.sig = true;
    } else {
      // -|x| + |y|
      // |y| - |x|
      x.cmp_minus(y);
      if (x.vec.size() == 1 && !x.vec[0]) {
        x.sig = false;
      } else {
        x.sig ^= 1;
      }
    }
  } else {
    if (y.sig) {
      // |x| - |y|
      x.cmp_minus(y);
    } else {
      x.unsigned_add(y);
    }
  }
  return x;
}

int2048 &int2048::minus(const int2048 &y) {
  // 有符号减法
  if (sig) {
    if (y.sig) {
      // -|x| + |y| = |y| - |x|
      this->cmp_minus(y);
      if (vec.size() == 1 && !vec[0]) {
        sig = false;
      } else {
        this->sig ^= 1;
      }
    } else {
      // -|x| - |y|
      // -(|x| + |y|)
      this->unsigned_add(y);
      sig = true;
    }
  } else {
    if (y.sig) {
      // |x| + |y|
      this->unsigned_add(y);
    } else {
      // |x| - |y|
      this->cmp_minus(y);
    }
  }
  return (*this);
}
int2048 minus(int2048 x, const int2048 &y) {
  if (x.sig) {
    if (y.sig) {
      // -|x| + |y| = |y| - |x|
      x.cmp_minus(y);
      if (x.vec.size() == 1 && !x.vec[0]) {
        x.sig = false;
      } else {
        x.sig ^= 1;
      }
    } else {
      // -|x| - |y|
      // -(|x| + |y|)
      x.unsigned_add(y);
      x.sig = true;
    }
  } else {
    if (y.sig) {
      // |x| + |y|
      x.unsigned_add(y);
    } else {
      // |x| - |y|
      x.cmp_minus(y);
    }
  }
  return x;
}

// FFT
int n;
std::vector<cp> a, b, omg, inv;
std::vector<int> ans;
void fft(std::vector<cp> &x, const std::vector<cp> &y) {
  int root = 0;
  while ((1 << root) < n) {
    ++root;
  }
  // 把 i 放到最后的位置上
  for (int i = 0; i < n; ++i) {
    int t = 0;
    for (int j = 0; j < root; ++j) {
      if ((i >> j) & 1) {
        t |= (1 << (root - j - 1));
      }
    } // t 为 i 翻转之后的结果
    if (i < t) {
      std::swap(x[i], x[t]);
    }
  }

  for (int j = 2; j <= n; j <<= 1) {
    int m = j >> 1; // 卡一半位置，分治
    for (int k = 0; k < n; k += j) { // 处理此时独立的各区间
      for (int i = 0; i < m; ++i) { // 扫左半边
        // w(cur(n) = j, i) == w(j * (n / j), i * n / j);
        cp t = y[n / j * i] * x[i + k + m];
        x[i + k + m] = x[i + k] - t;
        x[i + k] = x[i + k] + t;
      }
    }
  }
}
int2048& int2048::mul(const int2048 &y) {
  // 有符号乘法
  bool flag = !(sig == y.sig);
  if (vec.size() == 1 && !vec[0]) {
    sig = false;
    return (*this);
  }
  if (y.vec.size() == 1 && !y.vec[0]) {
    sig = false;
    vec.clear();
    vec.resize(1, 0);
    return (*this);
  } // 特判 0

  // 拆位
  a.clear();
  b.clear();
  long long cur;
  for (unsigned long i = 0; i < vec.size() - 1; ++i) {
    // 算九位进去
    cur = vec[i];
    for (int j = 0; j < 9; ++j) {
      a.push_back({(double)(cur % 10ll), 0});
      cur /= 10ll;
    }
  }
  cur = vec[vec.size() - 1];
  while (cur > 0ll) {
    a.push_back({(double)(cur % 10ll), 0});
    cur /= 10ll;
  }
  for (unsigned long i = 0; i < y.vec.size() - 1; ++i) {
    // 算九位进去
    cur = y.vec[i];
    for (int j = 0; j < 9; ++j) {
      b.push_back({(double)(cur % 10ll), 0});
      cur /= 10ll;
    }
  }
  cur = y.vec[y.vec.size() - 1];
  while (cur > 0ll) {
    b.push_back({(double)(cur % 10ll), 0});
    cur /= 10ll;
  }
  int lena = a.size(), lenb = b.size();
  n = 1;
  while (n < lena + lenb) {
      n <<= 1;
  }

  // 处理单位根
  a.resize(n);
  b.resize(n);
  omg.resize(n, 0);
  inv.resize(n, 0);
  for (int i = 0; i < n; ++i) {
    omg[i] = cp(cos(2 * pi * i / n), sin(2 * pi * i / n));
    inv[i] = std::conj(omg[i]);
  }

  fft(a, omg);
  fft(b, omg);
  for (int i = 0; i < n; ++i) {
    a[i] *= b[i];
  }
  // 把 a 换回系数表示
  fft(a, inv);
  ans.clear();
  ans.resize(n, 0);
  for (int i = 0; i < n; ++i) {
    ans[i] += floor(a[i].real() / n + 0.5);
    ans[i + 1] += ans[i] / 10;
    ans[i] %= 10;
  }
  int m;
  if (!ans[lena + lenb - 1]) {
    m = lena + lenb - 2;
  } else {
    m = lena + lenb - 1;
  }
  std::string ret;
  ret.resize(m + 1, 0);
  for (int i = m; i >= 0; --i) {
    ret[m - i] = ans[i] + '0';
  }
  this->read(ret);
  sig = flag;
  return (*this);
}

int2048& int2048::div_mul(const long long &y) {
  // 高精除法专用乘法，高精 * 单精，O(n)
  if (!y) {
    vec.clear();
    vec.push_back(0);
    sig = false;
    return (*this);
  }
  long long tmp_y = y;
  if (y < 0) {
    tmp_y = -y;
    sig ^= 1;
  }
  long long carry = 0;
  for (unsigned long i = 0; i < vec.size(); ++i) {
    vec[i] = vec[i] * tmp_y + carry;
    carry = vec[i] / BASE;
    vec[i] %= BASE;
  }
  if (carry) {
    vec.push_back(carry);
  }
  return (*this);
}
int2048 div_mul(int2048 x, const long long &y) {
  x.div_mul(y);
  return x;
}
int2048& int2048::div(const int2048 &y) {
  // 有符号除法
  // 从高位开始算起
  if (y.vec.size() == 1 && y.vec[0] == 0) {
    std::cout << "Error: Divided by zero." << std::endl;
    exit(1);
  }
  bool flag = (sig != y.sig);
  sig = false;
  int2048 tmp_y = y;
  tmp_y.sig = false;
  if (tmp_y > (*this)) {
    // 答案是 1
    sig = false;
    vec.clear();
    vec.push_back(0);
    return (*this);
  }

  int lenx = vec.size(), leny = y.vec.size();
  long long y_high = y.vec[y.vec.size() - 1], x_high = 0ll;
  std::vector<long long> ans;
  int2048 nowx;
  cut(*this, nowx, lenx - leny + 1, lenx - 1);
  nowx.div_mul(BASE);
  ans.reserve(lenx - leny + 1);
  for (int i = lenx - leny; i >= 0; --i) {
    nowx += (int2048)vec[i];
    // 试商
    if (nowx < tmp_y) {
      if (!ans.empty()) {
        ans.push_back(0);
      }
      nowx.div_mul(BASE);
      continue;
    }
    // 高位/高位试商
    x_high = nowx.vec[nowx.vec.size() - 1];
    if (x_high < y_high || nowx.vec.size() > y.vec.size()) {
      if (nowx.vec.size() >= 2) {
        x_high = x_high * BASE + nowx.vec[nowx.vec.size() - 2];
      }
      else {
        x_high *= BASE;
      }
    }
    
    long long l = x_high / (y_high + 1), r = (x_high + 1) / y_high;
    // <= x 的最大数
    while (l < r) {
      long long mid = (l + r + 1ll) >> 1ll;
      // 判定是否可行
      if (sjtu::div_mul(tmp_y, mid) <= nowx) {
        l = mid;
      } else {
        r = mid - 1;
      }
    }
    ans.push_back(l);
    nowx.unsigned_minus(sjtu::div_mul(tmp_y, l));
    nowx.div_mul(BASE);
  }
  vec.resize(ans.size());
  for (unsigned i = 0; i < ans.size(); ++i) {
    vec[i] = ans[ans.size() - i - 1];
  }
  sig = flag;
  return (*this);
}

int2048 &int2048::operator=(const int2048 &y) {
  sig = y.sig;
  vec.clear();
  vec.resize(y.vec.size());
  for (unsigned i = 0; i < y.vec.size(); ++i) {
    vec[i] = y.vec[i];
  }
  return (*this);
}
int2048 &int2048::operator=(const int2048 &&y) {
  sig = y.sig;
  vec = y.vec;
  vec.clear();
  vec.resize(y.vec.size());
  for (unsigned i = 0; i < y.vec.size(); ++i) {
    vec[i] = y.vec[i];
  }
  return (*this);
}

std::istream &operator>>(std::istream &in, int2048 &x) {
  // std::ios::sync_with_stdio(false);
  // std::cin.tie(nullptr);
  std::string s;
  in >> s;
  x.read(s);
  return in;
}
std::ostream &operator<<(std::ostream &os, const int2048 &x) {
  // std::ios::sync_with_stdio(false);
  // std::cout.tie(nullptr);
  if (x.sig) {
    os << '-';
  }
  int siz = x.vec.size();
  os << x.vec[siz - 1];
  for (int i = siz - 2; i >= 0; --i) {
    os << std::setw(9) << std::setfill('0') << x.vec[i];
  }
  return os;
}

int2048& int2048::operator+=(const int2048 &y) {
  this->add(y);
  return (*this);
}
int2048 operator+(int2048 x, const int2048 &y) {
  return add(x, y);
}
int2048& int2048::operator-=(const int2048 &y) {
  this->minus(y);
  return (*this);
}
int2048 operator-(int2048 x, const int2048 &y) {
  return minus(x, y);
}
int2048& int2048::operator*=(const int2048 &y) {
  this->mul(y);
  return (*this);
}
int2048 operator*(int2048 x, const int2048 &y) {
  x.mul(y);
  return x;
}
int2048& int2048::operator/=(const int2048 &y) {
  this->div(y);
  return (*this);
}
int2048 operator/(int2048 x, const int2048 &y) {
  x.div(y);
  return x;
}


bool operator==(const int2048 &x, const int2048 &y) {
  if (x.sig != y.sig) {
    return false;
  }
  if (x.vec.size() != y.vec.size()) {
    return false;
  }
  for (unsigned long i = 0; i < x.vec.size(); ++i) {
    if (x.vec[i] != y.vec[i]) {
      return false;
    }
  }
  return true;
}
bool operator!=(const int2048 &x, const int2048 &y) {
  return !(x == y);
}
bool operator<(const int2048 &x, const int2048 &y) {
  if (x.sig != y.sig) {
    if (x.sig == true) {
      return true;
    } else {
      return false;
    }
  }
  // 正负性相同
  bool sig = x.sig;
  if (x.vec.size() != y.vec.size()) {
    if (x.vec.size() < y.vec.size()) {
      return (sig ^ 1);
    } else {
      return sig;
    }
  }
  for (int i = x.vec.size() - 1; i >= 0; --i) {
    if (x.vec[i] != y.vec[i]) {
      if (x.vec[i] < y.vec[i]) {
        return (sig ^ 1);
      } else {
        return sig;
      }
    }
  }
  // 二者相等
  return false;
}
bool operator>(const int2048 &x, const int2048 &y) {
  if (x.sig != y.sig) {
    if (x.sig == true) {
      return false;
    } else {
      return true;
    }
  }
  // 正负性相同
  bool sig = x.sig;
  if (x.vec.size() != y.vec.size()) {
    if (x.vec.size() < y.vec.size()) {
      return sig;
    } else {
      return (sig ^ 1);
    }
  }
  for (int i = x.vec.size() - 1; i >= 0; --i) {
    if (x.vec[i] != y.vec[i]) {
      if (x.vec[i] < y.vec[i]) {
        return sig;
      } else {
        return (sig ^ 1);
      }
    }
  }
  // 二者相等
  return false;
}
bool operator<=(const int2048 &x, const int2048 &y) {
  if (x.sig != y.sig) {
    if (x.sig == true) {
      return true;
    } else {
      return false;
    }
  }
  // 正负性相同
  bool sig = x.sig;
  if (x.vec.size() != y.vec.size()) {
    if (x.vec.size() < y.vec.size()) {
      return (sig ^ 1);
    } else {
      return sig;
    }
  }
  for (int i = x.vec.size() - 1; i >= 0; --i) {
    if (x.vec[i] != y.vec[i]) {
      if (x.vec[i] < y.vec[i]) {
        return (sig ^ 1);
      } else {
        return sig;
      }
    }
  }
  // 二者相等
  return true;
}
bool operator>=(const int2048 &x, const int2048 &y) {
  if (x.sig != y.sig) {
    if (x.sig == true) {
      return false;
    } else {
      return true;
    }
  }
  // 正负性相同
  bool sig = x.sig;
  if (x.vec.size() != y.vec.size()) {
    if (x.vec.size() < y.vec.size()) {
      return sig;
    } else {
      return (sig ^ 1);
    }
  }
  for (int i = x.vec.size() - 1; i >= 0; --i) {
    if (x.vec[i] != y.vec[i]) {
      if (x.vec[i] < y.vec[i]) {
        return sig;
      } else {
        return (sig ^ 1);
      }
    }
  }
  // 二者相等
  return true;
}
double int_to_double(const int2048 &x) {
  double ret = 0;
  int siz = x.vec.size();
  for (int i = siz - 1; i >= 0; --i) {
    ret = ret * (double)BASE + x.vec[i];
  }

  if (x.sig) {
    ret = -ret;
  }
  return ret;
}

std::string int_to_string(const int2048 &x) {
  std::string str;
  int siz = x.vec.size();
  for (int i = 0; i < siz - 1; ++i) {
    int cur = x.vec[i];
    for (int j = 0; j < 9; ++j) {
      str.push_back(cur % 10 + '0');
      cur /= 10;
    }
  }
  // 最高位直接处理到 0
  int cur = x.vec[siz - 1];
  while (cur) {
    str.push_back(cur % 10 + '0');
    cur /= 10;
  }
  if (x.sig) {
    str.push_back('-');
  }
  std::reverse(str.begin(), str.end());
  return str;
}

long long to_ll(const int2048 &x) {
  long long ret = 0;
  int siz = x.vec.size();
  for (int i = siz - 1; i >= 0; --i) {
    ret = ret * BASE + x.vec[i];
  }
  if (x.sig) {
    ret = -ret;
  }
  return ret;
}

}  // namespace sjtu

int main() {
  using namespace sjtu;
  int2048 a, b, c, d;
  std::cin >> a >> b;

  std::cout << a / b << std::endl;
  return 0;
}

