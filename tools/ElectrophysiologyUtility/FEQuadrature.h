/*
 * File: FEQuadrature.h
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


#ifndef FEQUADRATURE_H
#define FEQUADRATURE_H

#include <array>
#include <stdexcept>
#include <memory>

class FETetCentroid;
class FETet4Point;
class FETet8Point;

class FETetQuadrature {
 public:
  virtual int                           GetDegree() const                     = 0;
  virtual int                           GetNumberOfPoints() const             = 0;
  virtual double                        GetWeight(int i) const                = 0;
  virtual std::array<double, 4> const & GetParametricCoordinates(int i) const = 0;
  virtual double const                 *GetWeights() const                    = 0;

  static std::unique_ptr<FETetQuadrature> GetDefaultRuleByDegree(int deg);
  static std::unique_ptr<FETetQuadrature> GetDefaultRuleByPoints(int deg);
};

template<int N, int D>
class FETetQuadratureBaseTemplate : public FETetQuadrature {
 public:
  virtual int GetNumberOfPoints() const final {return N;}

  virtual int GetDegree() const {return D;}

  virtual double GetWeight(int i) const {return w[i];}

  virtual std::array<double, 4> const & GetParametricCoordinates(int i) const {return z[i];}

  virtual double const *GetWeights() const {return w.data();}

 protected:
  FETetQuadratureBaseTemplate(std::array<double, N> const &&w_, std::array<std::array<double, 4>, N> const &&z_) {
    w = w_; z = z_;
  }

  std::array<double, N> w;
  std::array<std::array<double, 4>, N> z;

 private:
  FETetQuadratureBaseTemplate() {w = {0}; z = {0};}
};


class FETetCentroid : public FETetQuadratureBaseTemplate<1, 1> {
 public:
  FETetCentroid() : FETetQuadratureBaseTemplate({ {1}
                                                }, {{std::array<double, 4>{{.25, .25, .25, .25}
                                                     }
                                                    }
                                                }) {}
};


class FETet4Point : public FETetQuadratureBaseTemplate<4, 2> {
 public:
  FETet4Point() : FETetQuadratureBaseTemplate({ {.25, .25, .25, .25}
                                              }, {{
                                                    std::array<double, 4>{{1.-3*z, z, z, z}
                                                    }, {{z, 1.-3*z, z, z}
                                                    },
                                                    {{z, z, 1.-3*z, z}
                                                    }, {{z, z, z, 1.-3*z}
                                                    }
                                                  }
                                              }) {}

 protected:
  /* See wuelfers16, Appendix B.2.4 */
  static constexpr double z = .13819660112501051518;
};


class FETet8Point : public FETetQuadratureBaseTemplate<8, 3> {
 protected:
  /* See wuelfers16, Appendix B.2.5 */
  static const double w[2];
  static const double g[2];

 public:
  FETet8Point() : FETetQuadratureBaseTemplate({ {w[0], w[0], w[0], w[0], w[1], w[1], w[1], w[1]}
                                              },
                                              {{std::array<double, 4>{{1.-3*g[0], g[0], g[0], g[0]}
                                                }, {{g[0], 1.-3*g[0], g[0], g[0]}
                                                },
                                                {{g[0], g[0], 1.-3*g[0], g[0]}
                                                }, {{g[0], g[0], g[0], 1.-3*g[0]}
                                                },
                                                {{1.-3*g[1], g[1], g[1], g[1]}
                                                }, {{g[1], 1.-3*g[1], g[1], g[1]}
                                                },
                                                {{g[1], g[1], 1.-3*g[1], g[1]}
                                                }, {{g[1], g[1], g[1], 1.-3*g[1]}
                                                }
                                               }
                                              }) {}
};

const double FETet8Point::w[2] = {.13852796651186214232, .11147203348813785768};
const double FETet8Point::g[2] = {.32805469671142664733, .10695227393293068277};


std::unique_ptr<FETetQuadrature> FETetQuadrature::GetDefaultRuleByDegree(int deg) {
  switch (deg) {
    case 1:
      return std::unique_ptr<FETetQuadrature>(new FETetCentroid());

    case 2:
      return std::unique_ptr<FETetQuadrature>(new FETet4Point());

    case 3:
      return std::unique_ptr<FETetQuadrature>(new FETet8Point());

    default:
      throw std::runtime_error("Not implemented.");
  }
}

std::unique_ptr<FETetQuadrature> FETetQuadrature::GetDefaultRuleByPoints(int np) {
  switch (np) {
    case 1:
      return std::unique_ptr<FETetQuadrature>(new FETetCentroid());

    case 4:
      return std::unique_ptr<FETetQuadrature>(new FETet4Point());

    case 8:
      return std::unique_ptr<FETetQuadrature>(new FETet8Point());

    default:
      throw std::runtime_error("Not implemented.");
  }
}

#endif  // ifndef FEQUADRATURE_H
