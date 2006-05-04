//____________________________________________________________________________
/*
 Copyright (c) 2003-2006, GENIE Neutrino MC Generator Collaboration
 All rights reserved.
 For the licensing terms see $GENIE/USER_LICENSE.

 Author: Costas Andreopoulos <C.V.Andreopoulos@rl.ac.uk>
         CCLRC, Rutherford Appleton Laboratory - July 13, 2005

 For the class documentation see the corresponding header file.

 Important revisions after version 2.0.0 :

*/
//____________________________________________________________________________

#include "Conventions/Controls.h"
#include "Conventions/KinePhaseSpace.h"
#include "EVGCore/EVGThreadException.h"
#include "EVGModules/IMDKinematicsGenerator.h"
#include "GHEP/GHepRecord.h"
#include "GHEP/GHepFlags.h"
#include "Messenger/Messenger.h"
#include "Numerical/RandomGen.h"
#include "Utils/MathUtils.h"

using namespace genie;
using namespace genie::controls;

//___________________________________________________________________________
IMDKinematicsGenerator::IMDKinematicsGenerator() :
KineGeneratorWithCache("genie::IMDKinematicsGenerator")
{

}
//___________________________________________________________________________
IMDKinematicsGenerator::IMDKinematicsGenerator(string config) :
KineGeneratorWithCache("genie::IMDKinematicsGenerator", config)
{

}
//___________________________________________________________________________
IMDKinematicsGenerator::~IMDKinematicsGenerator()
{

}
//___________________________________________________________________________
void IMDKinematicsGenerator::ProcessEventRecord(GHepRecord * evrec) const
{
// Selects kinematic variables using the 'Rejection' method and adds them to
// the event record's summary

  Interaction * interaction = evrec->GetInteraction();

  //-- Get the random number generators
  RandomGen * rnd = RandomGen::Instance();

  //-- For the subsequent kinematic selection with the rejection method:
  //   Calculate the max differential cross section or retrieve it from the
  //   cache. Throw an exception and quit the evg thread if a non-positive
  //   value is found.
  double xsec_max = this->MaxXSec(evrec);

  //------ Try to select a valid inelastisity y
  register unsigned int iter = 0;
  const double e = 1E-6;

  double ymin = 0.0 + e; // the xsec algorithm would internally compute the
  double ymax = 1.0 - e; // kinematically allowd range, and return 0 if outside
  double dy   = ymax-ymin;

  while(1) {

     iter++;
     if(iter > kRjMaxIterations) {
        LOG("IMDKinematics", pWARN)
              << "*** Could not select a valid y after "
                                              << iter << " iterations";
        evrec->EventFlags()->SetBitNumber(kNoValidKinematics, true);
        genie::exceptions::EVGThreadException exception;
        exception.SetReason("Couldn't select kinematics");
        exception.SwitchOnFastForward();
        throw exception;
     }

     double y = ymin + dy * rnd->Random1().Rndm();
     interaction->GetKinematicsPtr()->Sety(y);

     LOG("IMDKinematics", pINFO) << "Trying: y = " << y;
     double xsec = fXSecModel->XSec(interaction, kPSyfE);
     double t    = xsec_max * rnd->Random1().Rndm();

     LOG("IMDKinematics", pINFO)
           << "xsec: (computed) = " << xsec << ", (generated) = " << t;
     assert(xsec < xsec_max);

     if(t < xsec) {
        // kinematical selection done.
        LOG("IMDKinematics", pINFO) << "Selected: y = " << y;

        // set the cross section for the selected kinematics
        evrec->SetDiffXSec(xsec);

        // lock selected kinematics & clear running values
        interaction->GetKinematicsPtr()->Sety(y, true);
        interaction->GetKinematicsPtr()->ClearRunningValues();

        return;
     }
  }// iterations
}
//___________________________________________________________________________
double IMDKinematicsGenerator::ComputeMaxXSec(
                                       const Interaction * interaction) const
{
// Computes the maximum differential cross section in the requested phase
// space. This method overloads KineGeneratorWithCache::ComputeMaxXSec
// method and the value is cached at a circular cache branch for retrieval
// during subsequent event generation.
// The computed max differential cross section does not need to be the exact
// maximum. The number used in the rejection method will be scaled up by a
// safety factor. But it needs to be fast - do not use a very small y step.

  const int N = 20;
  double max_xsec = -1.0;

  const double ymin = 0.;
  const double ymax = 1.;
  const double dy   = (ymax-ymin)/(N-1);

  for(int i=0; i<N; i++) {

     double y = ymin + i * dy;
     interaction->GetKinematicsPtr()->Sety(y);

     double xsec = fXSecModel->XSec(interaction, kPSyfE);
     max_xsec = TMath::Max(xsec, max_xsec);
  }//y

  // Apply safety factor, since value retrieved from the cache might
  // correspond to a slightly different energy.
  max_xsec *= fSafetyFactor;

  SLOG("IMDKinematics", pDEBUG) << interaction->AsString();
  SLOG("IMDKinematics", pDEBUG) << "Max xsec in phase space = " << max_xsec;
  SLOG("IMDKinematics", pDEBUG) << "Computed using alg = " << *fXSecModel;

  return max_xsec;
}
//___________________________________________________________________________
double IMDKinematicsGenerator::Energy(const Interaction * interaction) const
{
// Override the base class Energy() method to cache the max xsec for the
// neutrino energy in the LAB rather than in the hit nucleon rest frame.

  const InitialState & init_state = interaction->GetInitialState();
  double E = init_state.GetProbeE(kRfLab);
  return E;
}
//___________________________________________________________________________
void IMDKinematicsGenerator::Configure(const Registry & config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void IMDKinematicsGenerator::Configure(string config)
{
  Algorithm::Configure(config);
  this->LoadConfigData();
  this->LoadSubAlg();
}
//____________________________________________________________________________
void IMDKinematicsGenerator::LoadSubAlg(void)
{
  fXSecModel = dynamic_cast<const XSecAlgorithmI *> (
                            this->SubAlg("xsec-alg-name", "xsec-param-set"));
  assert(fXSecModel);
}
//____________________________________________________________________________
void IMDKinematicsGenerator::LoadConfigData(void)
{
  fSafetyFactor = fConfig->GetDoubleDef("max-xsec-safety-factor", 1.25);
  fEMin         = fConfig->GetDoubleDef("min-energy-cached",     -1.00);
}
//____________________________________________________________________________

