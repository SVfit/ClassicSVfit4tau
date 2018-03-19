#ifndef TauAnalysis_ClassicSVfit4tau_ClassicSVfitIntegrand4tau_h
#define TauAnalysis_ClassicSVfit4tau_ClassicSVfitIntegrand4tau_h

#include "TauAnalysis/ClassicSVfit/interface/ClassicSVfitIntegrandBase.h"
#include "TauAnalysis/ClassicSVfit/interface/MeasuredTauLepton.h"
#include "TauAnalysis/ClassicSVfit/interface/FittedTauLepton.h"
#include "TauAnalysis/ClassicSVfit4tau/interface/svFitHistogramAdapter4tau.h" // HistogramAdapterDiHiggs

#include <Math/Functor.h>
#include <TMatrixD.h>

namespace classic_svFit
{
  class ClassicSVfitIntegrand4tau : public ClassicSVfitIntegrandBase
  {
   public:
    ClassicSVfitIntegrand4tau(int);
    ~ClassicSVfitIntegrand4tau();

    void setDiTau1MassConstraint(double diTauMass);
    void setDiTau2MassConstraint(double diTauMass);

    /// set pointer to histograms used to keep track of pT, eta, phi, mass and transverse mass of di-Higgs system
    /// during Markov Chain integration
    void setHistogramAdapter(HistogramAdapterDiHiggs* histogramAdapter);

    /// set momenta of visible tau decay products
    void setLeptonInputs(const std::vector<classic_svFit::MeasuredTauLepton>&);

    /// evaluate Phase Space part of the integrand for given value of integration variables x
    double EvalPS(const double* x) const;

    /// evaluate the iComponent of the full integrand for given value of integration variables q.
    /// q is given in standarised range [0,1] for each dimension.
    double Eval(const double* q, unsigned int iComponent=0) const;

    /// static pointer to this (needed for interfacing the likelihood function calls to Markov Chain integration)
    static const ClassicSVfitIntegrand4tau* gSVfitIntegrand;

   protected:
    /// momenta of visible tau decay products and of reconstructed tau leptons
    MeasuredTauLepton measuredTauLepton1_;    
    mutable FittedTauLepton fittedTauLepton1_;
    bool leg1isLeptonicTauDecay_;
    bool leg1isHadronicTauDecay_;
    bool leg1isPrompt_;
    MeasuredTauLepton measuredTauLepton2_;  
    mutable FittedTauLepton fittedTauLepton2_;
    bool leg2isLeptonicTauDecay_;
    bool leg2isHadronicTauDecay_;
    bool leg2isPrompt_;
    MeasuredTauLepton measuredTauLepton3_;    
    mutable FittedTauLepton fittedTauLepton3_;
    bool leg3isLeptonicTauDecay_;
    bool leg3isHadronicTauDecay_;
    bool leg3isPrompt_;
    MeasuredTauLepton measuredTauLepton4_;  
    mutable FittedTauLepton fittedTauLepton4_;
    bool leg4isLeptonicTauDecay_;
    bool leg4isHadronicTauDecay_;
    bool leg4isPrompt_;

    /// visible masses of first and second tau pair
    mutable double ditau1_mVis_measured_; 
    mutable double ditau1_mVis2_measured_;
    mutable double ditau2_mVis_measured_; 
    mutable double ditau2_mVis2_measured_;

    double diTau1MassConstraint_;
    double diTau1MassConstraint2_;
    double diTau2MassConstraint_;
    double diTau2MassConstraint2_;

    HistogramAdapterDiHiggs* histogramAdapter_;
  };
}

#endif
