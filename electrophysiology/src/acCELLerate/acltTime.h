/*
 * File: acltTime.h
 *
 * Institute of Biomedical Engineering, 
 * Karlsruhe Institute of Technology (KIT)
 * https://www.ibt.kit.edu
 * 
 * Repository: https://github.com/KIT-IBT/CardioMechanics
 *
 * License: GPL-3.0 (See accompanying file LICENSE or visit https://www.gnu.org/licenses/gpl-3.0.html)
 *
 */


#ifndef ACLT_TIME_H
#define ACLT_TIME_H

#include <string>
#include <cinttypes>
#include <iostream>
#include <regex>
#include <stdexcept>
#include <cmath>
#include <limits>

template<int P, typename T = int64_t>
struct POW {
  constexpr static inline T Get(T val) {
    return val * POW<P - 1, T>::Get(val);
  }
};

template<typename T>
struct POW<0, T> {
  constexpr static inline T Get(T) {return 1;}
};

template<int DECIMALS, class T>
class acltTimeTemplate {
 public:
  acltTimeTemplate() : value(0) {}

  acltTimeTemplate(const std::string &v) {this->parse(v);}

  acltTimeTemplate(const char *v) {if (v) this->parse(std::string(v)); else this->value = 0; }

  explicit acltTimeTemplate(double v) : value(v * POW<DECIMALS>::Get(10)) { /* TODO check value maybe? */}

  std::string toStr(int precision = -1) const {return this->toStr(0, precision);}

  std::string toStr(int, int) const;

  inline bool operator==(const acltTimeTemplate<DECIMALS, T> &rhs) const {return this->value == rhs.value;}

  inline bool operator!=(const acltTimeTemplate<DECIMALS, T> &rhs) const {return !(*this == rhs);}

  inline bool operator<(const acltTimeTemplate<DECIMALS, T> &rhs) const {return this->value < rhs.value;}

  inline bool operator>(const acltTimeTemplate<DECIMALS, T> &rhs) const {return rhs < *this;}

  inline bool operator<=(const acltTimeTemplate<DECIMALS, T> &rhs) const {return !(rhs < *this);}

  inline bool operator>=(const acltTimeTemplate<DECIMALS, T> &rhs) const {return !(*this < rhs);}

  inline acltTimeTemplate<DECIMALS, T> operator-() {
    acltTimeTemplate<DECIMALS, T> a = *this; a.value = -this->value; return a;
  }

  inline acltTimeTemplate<DECIMALS, T> & operator+=(const acltTimeTemplate<DECIMALS, T> &rhs) {
    this->value += rhs.value; return *this;
  }

  inline acltTimeTemplate<DECIMALS, T> & operator-=(const acltTimeTemplate<DECIMALS, T> &rhs) {
    this->value -= rhs.value; return *this;
  }

  inline acltTimeTemplate<DECIMALS, T> & operator*=(const acltTimeTemplate<DECIMALS, T> &rhs) {
    this->value *= rhs.value; return *this;
  }

  inline acltTimeTemplate<DECIMALS, T> & operator/=(const acltTimeTemplate<DECIMALS, T> &rhs) {
    this->value /= rhs.value; return *this;
  }

  inline acltTimeTemplate<DECIMALS, T> & operator%=(const acltTimeTemplate<DECIMALS, T> &rhs) {
    this->value %= rhs.value; return *this;
  }

  inline acltTimeTemplate<DECIMALS, T> & operator*=(const T &rhs) {this->value *= rhs; return *this;}

  inline acltTimeTemplate<DECIMALS, T> & operator/=(const T &rhs) {this->value /= rhs; return *this;}

  explicit inline operator double() {
    return value * (1. / POW<DECIMALS>::Get(10));
  }

  explicit inline operator T() {
    return value / POW<DECIMALS>::Get(10);
  }

  static constexpr T max() {
    return std::numeric_limits<T>::max() / POW<DECIMALS>::Get(10);
  }

  static constexpr T min() {
    return std::numeric_limits<T>::min() / POW<DECIMALS>::Get(10);
  }

 private:
  void parse(const std::string &);
  T value;
};  // class acltTimeTemplate

template<int DECIMALS, class T>
inline std::ostream & operator<<(std::ostream &lhs, acltTimeTemplate<DECIMALS, T> const &rhs) {
  return lhs << rhs.toStr(lhs.width(), lhs.precision());
}

template<int DECIMALS, class T>
inline std::istream & operator>>(std::istream &lhs, acltTimeTemplate<DECIMALS, T> &rhs) {
  std::string buf;

  lhs >> buf;
  rhs = buf;
  return lhs;
}

template<int DECIMALS, class T>
inline acltTimeTemplate<DECIMALS, T> operator+(acltTimeTemplate<DECIMALS, T> lhs,
                                               const acltTimeTemplate<DECIMALS, T> &rhs) {
  lhs += rhs;
  return lhs;
}

template<int DECIMALS, class T>
inline acltTimeTemplate<DECIMALS, T> operator-(acltTimeTemplate<DECIMALS, T> lhs,
                                               const acltTimeTemplate<DECIMALS, T> &rhs) {
  lhs -= rhs;
  return lhs;
}

template<int DECIMALS, class T>
inline acltTimeTemplate<DECIMALS, T> operator*(acltTimeTemplate<DECIMALS, T> lhs,
                                               const acltTimeTemplate<DECIMALS, T> &rhs) {
  lhs *= rhs;
  return lhs;
}

template<int DECIMALS, class T>
inline acltTimeTemplate<DECIMALS, T> operator/(acltTimeTemplate<DECIMALS, T> lhs,
                                               const acltTimeTemplate<DECIMALS, T> &rhs) {
  lhs /= rhs;
  return lhs;
}

template<int DECIMALS, class T>
inline acltTimeTemplate<DECIMALS, T> operator%(acltTimeTemplate<DECIMALS, T>        lhs,
                                               const acltTimeTemplate<DECIMALS, T> &rhs) {
  lhs %= rhs;
  return lhs;
}

template<int DECIMALS, class T>
inline acltTimeTemplate<DECIMALS, T> operator*(acltTimeTemplate<DECIMALS, T> lhs, const T &rhs) {
  lhs *= rhs;
  return lhs;
}

template<int DECIMALS, class T>
inline acltTimeTemplate<DECIMALS, T> operator/(acltTimeTemplate<DECIMALS, T> lhs, const T &rhs) {
  lhs /= rhs;
  return lhs;
}

template<int DECIMALS, class T>
std::string acltTimeTemplate<DECIMALS, T>::toStr(int leading, int precision) const {
  char buf[64];
  int  len = 0;

  int fw = (leading < 0 ? leading - DECIMALS : leading + DECIMALS);

  if ((this->value == 0) && (leading == 0)) {
    fw++;
  }

  // TODO : automatically determine format specifier or use to_string()
  long long v = value;
  snprintf(buf, sizeof buf, "%0*lld%n", fw, v, &len);

  std::string left(buf, len - DECIMALS);
  if (left.empty()) {
    left = "0";
  }
  std::string rest(buf + len - DECIMALS, precision < 0 ? DECIMALS : precision);
  if (precision < 0) {
    size_t notzero = rest.find_last_not_of("0");
    if (notzero == std::string::npos) {
      return left;
    } else {
      return left + "." + rest.substr(0, notzero+1);
    }
  } else if (precision == 0) {
    return left;
  } else {
    return left + "." + rest;
  }
}  // >::toStr

template<int DECIMALS, class T>
void acltTimeTemplate<DECIMALS, T>::parse(std::string const &str) {
  // Parsing strings of the format:
  // -1.1e-1.1

  this->value = 0;
  char *e = strdup(str.c_str());
  char *s = e;

  int exponent = DECIMALS;

  while (std::isspace(*s)) ++s;

  // TODO make sure T is signed if s is negative number...
  // Maybe also make sure T is large enough...
  this->value = std::strtoll(s, &s, 10);

  if (*s == '.') {
    while (*++s >= '0' && *s <= '9') {
      if (exponent) {
        this->value = this->value * 10 + *s - '0';
        --exponent;
      } else {
        long remainder = std::strtol(s, &s, 10);
        if (remainder) {
          throw std::runtime_error("The value " + str + " can not be represented exactly with this data type.");
        }
      }
    }
  }

  if ((*s == 'e') || (*s == 'E')) {
    exponent += strtol(++s, &s, 10);
    if (*s == '.') {
      throw std::runtime_error("No non-integer exponents possible at this time.");
    }
  }

  for (int i = 0; i < exponent; ++i) {
    // TODO check fpr overflow (e.g. <= max()/(10^exponent) or so, which is a constexpr...
    this->value *= 10;
  }

  free(e);
}  // >::parse

#if 0

// TODO: regex parse is elegant but too slow (think acltConditions file...)

template<int DECIMALS, class T>
void acltTimeTemplate<DECIMALS, T>::parse(std::string const &str) {
  std::regex re("([+-]?[0-9]*)(\\.([0-9]*))?([eE]([+-]?[0-9]*)(\\.([0-9]*))?)?", std::regex::extended);
  std::smatch match;

  if (!std::regex_match(str, match, re)) {
    throw std::runtime_error("Unable to parse the number: " + str);
  }

  // 1: integer part
  if (match[1].length()) {
    this->value = std::stoll(match[1].str()) * POW<DECIMALS>::Get(10);
  } else if (match[2].length() && !match[3].length()) {
    throw std::runtime_error("Unable to parse, decimal point without numbers: " + str);
  }

  // 2: decimal portion
  if (match[3].length()) {
    T val = std::stoll(match[3].str().substr(0, DECIMALS));
    for (int i = match[3].length(); i < DECIMALS; ++i) {
      val *= 10;
    }
    if (this->value < 0)
      this->value -= val;
    else
      this->value += val;
  }

  // exponential part
  if (match[4].length()) {
    bool neg_exp = false;
    if (match[5].length()) {
      int exp = std::stoi(match[5].str());
      if (exp > 0) {
        for (int i = 0; i < exp; ++i) {
          value *= 10;
        }
      } else {
        neg_exp = true;
        for (int i = 0; i < -exp; ++i) {
          value /= 10;
        }
      }
    } else if (!match[7].length()) {
      throw std::runtime_error("Unable to parse exponent: " + str);
    }

    // only lossy operation: decimal exponents...
    if (match[7].length()) {
      double exp = std::stod(match[6].str());
      if (neg_exp) {
        value /= std::pow(10, exp);
      } else {
        value *= std::pow(10, exp);
      }
    }
  }
}  // >::parse

#endif  // if 0

typedef acltTimeTemplate<10, int64_t> acltTime;
inline acltTime operator"" _s(const char *s) {return acltTime(std::string(s));}


#if TEST

# define TEST(a) do {                                                                                          \
    try {                                                                                                      \
      if (!(a)) {                                                                                              \
        std::cerr << "Test failed: " << # a << std::endl;                                                      \
      }                                                                                                        \
    } catch (std::exception &e) {                                                                              \
      std::cerr << "Test failed with exception: " << # a << " (Exception: '" << e.what() << "')" << std::endl; \
    } catch (...) {                                                                                            \
      std::cerr << "Test failed with an unknown exception: " << # a << std::endl;                              \
    }                                                                                                          \
} while (0)
# define THROWS(a) do {                                                             \
    try {                                                                           \
      (a); std::cerr << "Test sould have thrown but did not: " << # a << std::endl; \
    } catch (...) {                                                                 \
    }                                                                               \
} while (0)

int main(int argc, char *argv[]) {
  typedef acltTimeTemplate<10, int64_t> act;


  TEST(act("0").toStr() == "0");
  TEST(act("0e0").toStr() == "0");
  TEST(act("0.").toStr() == "0");
  TEST(act("0.0").toStr() == "0");
  TEST(act("0e1").toStr() == "0");
  TEST(act("0e1.1").toStr() == "0");
  TEST(act("0e1.").toStr() == "0");
  TEST(act("0.e1.").toStr() == "0");
  TEST(act("0.e.0").toStr() == "0");
  TEST(act("10").toStr(5, 0) == "00010");
  TEST(act("1e1").toStr(5, 0) == "00010");
  TEST(act("1e01").toStr(5, 0) == "00010");
  TEST(act("1e1.").toStr(5, 0) == "00010");
  TEST(act("1.e1.").toStr(5, 0) == "00010");

  THROWS(act(".") );
  THROWS(act("1e.") );
  THROWS(act(".e") );
  THROWS(act(".e1") );
  THROWS(act(".e1.") );
  THROWS(act(".e1.0") );
  THROWS(act(".e.") );
  THROWS(act(".1e.") );
  THROWS(act("1.e.") );
  THROWS(act("1.e-.") );

  act a("1"), b("1"), c("2");
  TEST(a == b);
  TEST(-a != b);
  TEST(-a == -b);
  TEST(a+b == b+a);
  TEST(a-b == b-a);
  TEST(a+b == c);
  TEST(a == 1_s);
  TEST(act("1e-4") == 0.0001_s);

  TEST( (double)1_s == 1.);
  TEST( (double)2e-1_s == .2);
  TEST(2e-5_s == 0.00002_s);

  TEST(.001_s == 0.001_s);
  TEST(act(".001") == 0.001_s);
  TEST(act(".001") == act(0.001) );

  return 0;
}  // main

#endif  // if TEST

#endif  // ifndef ACLT_TIME_H
