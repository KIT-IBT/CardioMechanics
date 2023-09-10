// created by Robin Moss

// based on the available Matlab code from cemrg.co.uk

#ifndef LAND17PARAMETERS_H
#define LAND17PARAMETERS_H

/// Heart failure specific alterations from Bollen et al. 2017. Implemented by Albert Dasi
#define HFM

#include <ParameterLoader.h>

namespace NS_Land17Parameters {
enum varType {
  VT_perm50 = vtFirst,
  VT_TRPN_n,
  VT_koff,
  VT_dr,
  VT_wfrac,
  VT_A_eff,
  VT_ktm_unblock,
  VT_beta_1,
  VT_beta_0,
  VT_gamma,
  VT_gamma_wu,
  VT_phi,
  VT_nperm,
  VT_ca50,
  VT_Tref,
  VT_nu,
  VT_mu,
  VT_xi,
  VT_a,
  VT_k,
  VT_eta_l,
  VT_eta_s,

  VT_k_ws,
  VT_k_uw,
  VT_cdw,
  VT_cds,
  VT_k_wu,
  VT_k_su,
  VT_A,
  VT_XSSS,
  VT_XWSS,
  VT_ktm_block,

  VT_XS_init,
  VT_XW_init,
  VT_TRPN_init,
  VT_TmBlocked_init,
  VT_ZETAS_init,
  VT_ZETAW_init,
  VT_Cd_init,
  #ifdef HFM
  VT_HFM_Multiplier,
  #endif // ifdef HFM
  vtLast
};
}  // namespace NS_Land17Parameters

using namespace NS_Land17Parameters;

class Land17Parameters : public vbNewForceParameters {
 public:
  Land17Parameters(const char *);
  ~Land17Parameters() {}

  void Init(const char *);
  void Calculate();
  void PrintParameters();
};

#endif  // ifndef LAND17PARAMETERS_H
