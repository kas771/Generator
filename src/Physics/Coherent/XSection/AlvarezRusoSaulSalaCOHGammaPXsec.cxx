//____________________________________________________________________________
/*
 Copyright (c) 2003-2019, The GENIE Collaboration
 For the full text of the license visit http://copyright.genie-mc.org
 or see $GENIE/LICENSE

 Author: Marco Roda
         University of Liverpool 
         August 23, 2019

 For the class documentation see the corresponding header file.

*/
//____________________________________________________________________________

#include <iostream>

#include <TMath.h>

#include "Framework/Algorithm/AlgConfigPool.h"
#include "Physics/XSectionIntegration/XSecIntegratorI.h"
#include "Framework/Conventions/GBuild.h"
#include "Framework/Conventions/Constants.h"
#include "Framework/Conventions/Units.h"
#include "Framework/Conventions/RefFrame.h"
#include "Framework/Messenger/Messenger.h"
#include "Framework/ParticleData/PDGCodes.h"
#include "Framework/ParticleData/PDGUtils.h"
#include "Framework/Utils/HadXSUtils.h"
#include "Framework/Utils/KineUtils.h"

#include "Physics/Coherent/XSection/AlvarezRusoSaulSalaCOHGammaPXsec.h"


using namespace genie;
using namespace genie::constants;
using namespace genie::utils;

//____________________________________________________________________________
AlvarezRusoSaulSalaCOHGammaPXsec::AlvarezRusoSaulSalaCOHGammaPXsec() :
XSecAlgorithmI("genie::AlvarezRusoSaulSalaCOHGammaPXsec")
{
  ;
}
//____________________________________________________________________________
AlvarezRusoSaulSalaCOHGammaPXsec::AlvarezRusoSaulSalaCOHGammaPXsec(string config) :
XSecAlgorithmI("genie::AlvarezRusoSaulSalaCOHGammaPXsec", config)
{
  ;
}
//____________________________________________________________________________
AlvarezRusoSaulSalaCOHGammaPXsec::~AlvarezRusoSaulSalaCOHGammaPXsec()
{
  ;
}
//____________________________________________________________________________
double AlvarezRusoSaulSalaCOHGammaPXsec::XSec( const Interaction * interaction, 
					       KinePhaseSpace_t kps) const
{

  if(! this -> ValidProcess    (interaction) ) return 0.;
  if(! this -> ValidKinematics (interaction) ) return 0.;

  const Kinematics &   kinematics = interaction -> Kine();
  const InitialState & init_state = interaction -> InitState();
  
  int A       = init_state.Tgt().A(); // mass number
  int Z       = init_state.Tgt().Z(); // atomic number
  
  std::string mode;
  if ( init_state.ProbePdg() > 0) {
    mode = "nu";
  } else {
    mode = "nubar";
  }
  
  std::string nucleus;
  if(Z == 6 && A == 12) nucleus = "12C";
  else if(Z == 18 && A == 40) nucleus = "40Ar";
  else {
    LOG( "AlvarezRusoSaulSala",pERROR) << "WARNING: Z = " << Z 
				       <<" and A = " << A 
				       << ";  wrong nucleus for NCgamma cross section. " ;
    return -1. ; 
  }

  Diff_Cross_Section dxsec(mode, nucleus) ;

  double E_nu = init_state.ProbeE(kRfLab); // neutrino energy
  
  const TLorentzVector p4_lep = kinematics.FSLeptonP4();
  const TLorentzVector p4_gamma = kinematics.HadSystP4();
  double E_lep = p4_lep.E();

  double xsec = dxsec.getDiffCrossSection( E_nu, 
					   E_lep, p4_lep.Theta(), 
					   p4_gamma.Theta(), p4_gamma.Phi());

  // This cross section is the following
  //
  //             d4 sigma
  // ------------------------------
  //  dE_g dTheta_l dTheta_g dPhi_g
  //
  
  // double xsec = fMultidiff->DXSec(E_nu, E_lep, p4_lep.Theta(), p4_lep.Phi(), p4_pi.Theta(), p4_pi.Phi());
  // xsec = xsec * 1E-38 * units::cm2;
  
  if (kps != kPSEgTlTgPgfE ) {
    xsec *= utils::kinematics::Jacobian(interaction, kPSEgTlTgPgfE, kps );
  }

  return xsec ;
}
//____________________________________________________________________________
double AlvarezRusoSaulSalaCOHGammaPXsec::Integral(const Interaction * interaction) const
{
  double xsec = fXSecIntegrator->Integrate(this,interaction);
  return xsec;
}
//____________________________________________________________________________
bool AlvarezRusoSaulSalaCOHGammaPXsec::ValidProcess(const Interaction * interaction) const
{
  if(interaction->TestBit(kISkipProcessChk)) return true;

  const InitialState & init_state = interaction->InitState();
  const ProcessInfo &  proc_info  = interaction->ProcInfo();
  const Target &       target     = init_state.Tgt();
  const XclsTag &      xcls       = interaction -> ExclTag() ;

  int nu = init_state.ProbePdg();

  if ( ! proc_info.IsCoherent() )  return false;
  if ( ! proc_info.IsWeakNC() )    return false;

  if ( target.A() <= 1 )          return false;  
  if ( target.HitNucIsSet() )     return false;
  
  int A = target.A() ;
  int Z = target.Z() ;
  
  // only 12C and 40Ar are implemented at the moment 
  if ( ( A!=12 || Z!= 6 ) && ( A!=40 || Z!=18 ) ) return false ;
  
  if (!pdg::IsNeutrino(nu) && !pdg::IsAntiNeutrino(nu)) return false;
  
  return true;
}
//____________________________________________________________________________
void AlvarezRusoSaulSalaCOHGammaPXsec::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AlvarezRusoSaulSalaCOHGammaPXsec::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfig();
}
//____________________________________________________________________________
void AlvarezRusoSaulSalaCOHGammaPXsec::LoadConfig(void)
{
  
  //-- load the differential cross section integrator
  fXSecIntegrator =
      dynamic_cast<const XSecIntegratorI *> (this->SubAlg("XSec-Integrator"));
  assert(fXSecIntegrator);

}
//____________________________________________________________________________

