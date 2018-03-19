
/**
   \class testClassicSVfit4tau testClassicSVfit4tau.cc "TauAnalysis/ClassicSVfit/bin/testClassicSVfit4tau.cc"
   \brief Basic example of the use of the standalone version of the "classic" SVfit algorithm, customized for the hh->4tau analysis
*/

#include "TauAnalysis/ClassicSVfit4tau/interface/ClassicSVfit4tau.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit4tau/interface/svFitHistogramAdapter4tau.h"
//#include "TauAnalysis/SVfitTF/interface/HadTauTFCrystalBall2.h"

#include "TH1F.h"

using namespace classic_svFit;

int main(int argc, char* argv[])
{
  /*
     This is a single event for testing purposes.
  */

  // define MET
  double measuredMETx = 13.1779;
  double measuredMETy = 54.9923;

  // define MET covariance
  TMatrixD covMET(2, 2);
  covMET[0][0] = 829.442;
  covMET[1][0] = 161.376;
  covMET[0][1] = 161.376;
  covMET[1][1] = 574.186;

  // define lepton four vectors
  std::vector<MeasuredTauLepton> measuredTauLeptons;
  // decay products of first Higgs bosons
  measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToMuDecay,   21.4323, -0.65693,   0.800028, 0.10566));     // tau -> muon decay (Pt, eta, phi, mass)
  measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay,  37.126,   0.174097,  2.16421,  0.762756, 1)); // tau -> 1prong1pi0 hadronic decay (Pt, eta, phi, mass)
  // decay products of second Higgs bosons
  measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToMuDecay,   30.3555, -0.611901, -1.10589,  0.10566));     // tau -> muon decay (Pt, eta, phi, mass)
  measuredTauLeptons.push_back(MeasuredTauLepton(MeasuredTauLepton::kTauToHadDecay, 196.753,  -0.276988, -2.68803,  0.85196,  1)); // tau -> 1prong1pi0 hadronic decay (Pt, eta, phi, mass)
  /*
     tauDecayModes:  0 one-prong without neutral pions
                     1 one-prong with neutral pions
                    10 three-prong without neutral pions
  */

  int verbosity = 1;
  ClassicSVfit4tau svFitAlgo(verbosity);
#ifdef USE_SVFITTF
  //HadTauTFCrystalBall2* hadTauTF = new HadTauTFCrystalBall2();
  //svFitAlgo.setHadTauTF(hadTauTF);
  //svFitAlgo.enableHadTauTF();
#endif

  // run with mass constraint for each tau pair to match measured mass (125.06 GeV) of SM-like Higgs boson 
  double massContraint = 125.06;
  //svFitAlgo.addLogM_fixed(false);
  svFitAlgo.addLogM_fixed(true, 12.);
  //svFitAlgo.addLogM_dynamic(true, "(m/1000.)*15.");
  //svFitAlgo.setMaxObjFunctionCalls(100000); // CV: default is 100000 evaluations of integrand per event
  svFitAlgo.setLikelihoodFileName("testClassicSVfit.root");
  svFitAlgo.setDiTau1MassConstraint(massContraint);
  svFitAlgo.setDiTau2MassConstraint(massContraint);

  svFitAlgo.integrate(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
  bool isValidSolution_1stRun = svFitAlgo.isValidSolution();
  double mass_1stRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo.getHistogramAdapter())->getMass();
  double massErr_1stRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo.getHistogramAdapter())->getMassErr();
  double transverseMass_1stRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo.getHistogramAdapter())->getTransverseMass();
  double transverseMassErr_1stRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo.getHistogramAdapter())->getTransverseMassErr();

  if ( isValidSolution_1stRun ) {
    std::cout << "found valid solution: mass = " << mass_1stRun << " +/- " << massErr_1stRun << " (expected value = 395.003 +/- 12.8448),"
              << " transverse mass = " << transverseMass_1stRun << " +/- " << transverseMassErr_1stRun << " (expected value = 384.283 +/- 11.4043)" << std::endl;
  } else {
    std::cout << "sorry, failed to find valid solution !!" << std::endl;
  }
  if (std::abs((mass_1stRun - 395.003) / 395.003) > 0.001) return 1;
  if (std::abs((massErr_1stRun - 12.8448) / 12.8448) > 0.001) return 1;
  if (std::abs((transverseMass_1stRun - 384.283) / 384.283) > 0.001) return 1;
  if (std::abs((transverseMassErr_1stRun - 11.4043) / 11.4043) > 0.001) return 1;
 
  // re-run without mass constraint for each tau pair;
  // this mode will set the mass of the second tau pair to match the mass of the first tau pair,
  // while the mass of the first tau pair is allowed to freely vary within the fit
  std::cout << "\n\nTesting integration with di tau mass constraint set to " << massContraint << std::endl;
  svFitAlgo.setLikelihoodFileName("testClassicSVfit_withMassContraint.root");
  svFitAlgo.setDiTau1MassConstraint(-1.);
  svFitAlgo.setDiTau2MassConstraint(-1.);
  svFitAlgo.integrate(measuredTauLeptons, measuredMETx, measuredMETy, covMET);
  bool isValidSolution_2ndRun = svFitAlgo.isValidSolution();
  double mass_2ndRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo.getHistogramAdapter())->getMass();
  double massErr_2ndRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo.getHistogramAdapter())->getMassErr();
  double transverseMass_2ndRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo.getHistogramAdapter())->getTransverseMass();
  double transverseMassErr_2ndRun = static_cast<HistogramAdapterDiHiggs*>(svFitAlgo.getHistogramAdapter())->getTransverseMassErr();

  if ( isValidSolution_2ndRun ) {
    std::cout << "found valid solution: mass = " << mass_2ndRun << " +/- " << massErr_2ndRun << " (expected value = 436.01 +/- 90.3875),"
              << " transverse mass = " << transverseMass_2ndRun << " +/- " << transverseMassErr_2ndRun << " (expected value = 424.177 +/- 85.6547)" << std::endl;
  } else {
    std::cout << "sorry, failed to find valid solution !!" << std::endl;
  }
  if (std::abs((mass_2ndRun - 436.01) / 436.01) > 0.001) return 1;
  if (std::abs((massErr_2ndRun - 90.3875) / 90.3875) > 0.001) return 1;
  if (std::abs((transverseMass_2ndRun - 424.177) / 424.177) > 0.001) return 1;
  if (std::abs((transverseMassErr_2ndRun - 85.6547) / 85.6547) > 0.001) return 1;

  return 0;
}
