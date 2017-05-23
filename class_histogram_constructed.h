#ifndef __class_histogram_constructed_h__
#define __class_histogram_constructed_h__ 1

#include <vector>
#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <iomanip>
#include <stdlib.h>
#include <time.h>

#include "TROOT.h"
#include "TApplication.h"
#include "TTree.h"
#include "TH1D.h"
#include "TH1.h"
#include "TFile.h"
#include "TStyle.h"
#include "TString.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TLatex.h"
#include "TGraphErrors.h"
#include "TPaveLabel.h"
#include "TAxis.h"
#include "TFile.h"
#include "TSpectrum.h"
#include "TGraph.h" 

/// \file class_histogram_constructed.h
/// \brief Histogram_constructed

using namespace std;

/// \class histogram_constructed
/// \brief Class containing the root histograms
class Histogram_constructed
{
 protected :
  TString filename; ///< Root filename without .root 
  TFile *file; ///< Root file
  TH1D* histo[4]; ///< Array containing the histograms of each detector
 
 public :
  /// \brief Default constructor
  Histogram_constructed();
  /// \brief Surdefine Constructor
  Histogram_constructed(TString name);
  /// \brief Default destructor
  ~Histogram_constructed();

  /// \brief Get histograms
  /// \return An array with the four histograms
  vector<TH1D*> Get_histo() const;

  /// \brief operator += for Histo
  /// \brief use to make the sum of the histo
  /// \brief add the bincontent to the object
  /// \param hist: the Histo to add
  void operator+=(const Histogram_constructed & hist);

  /// \brief Get content of bin
  /// \param bin 
  /// \return a vector containing the content of bin for the four detectors
  vector<Int_t> Get_bin(Int_t bin) const;
  
  /// \brief Get background of histograms
  /// \param detector : number of the detector
  /// \param bin nb_bin :the peak start at bin and is nb_bin wide 
  /// \param bin_g bin_d : the background is evaluated between bin_g and bin_d
  /// \return The background of bin contents 
  Int_t Background(Int_t detector, Int_t bin, Int_t nb_bin, Int_t bin_g, Int_t bin_d, Double_t & Err) const;

  /// \brief Get integral of histograms
  /// \param detector : number of the detector
  /// \param bin nb_bin :the peak start at bin and is nb_bin wide 
  /// \param bin_g bin_d : the background is evaluated between bin_g and bin_d
  /// \return The sum of bin contents 
  Int_t Integrate(Int_t detector, Int_t bin, Int_t nb_bin, Int_t bin_g, Int_t bin_d, Double_t & Err) const;

  /// \brief Get the peak positions of histograms
  /// \return An array containing the peak positions for the four detectors
  vector<vector<Int_t> > Get_peakpos() const;
  
  /// \brief Compute the temperature of the sample with the 60Co peak at 1173keV
  /// \param &warm : reference of a warm Histogram
  /// \param &temp &Errtemp : the temperature and error calculated
  void temp1(const Histogram_constructed & warm, Int_t nbfileWarm, Double_t &temp, Double_t &ErrTemp) const;

  /// \brief Compute the temperature of the sample with the 60Co peak at 1332keV
  /// \param &warm : reference of a warm Histogram
  /// \param &temp &Errtemp : the temperature and error calculated
  void temp2(const Histogram_constructed & warm, Int_t nbfileWarm, Double_t &temp, Double_t &ErrTemp) const;
  
  /// \brief Make the calibration of the detector bin to keV
  /// \return a vector with the parameter of the calibration for the four detector
  vector<vector<Double_t> > Calibrate() const;

  /// \brief Make the inverse calibration of the detector keV to bin
  /// \return a vector with the parameter of the inverse calibration for the four detector
    vector<vector<Double_t> > CalibInv() const;
};

#endif


