#ifndef __class_histogram_h__
#define __class_histogram_h__ 1

#include "class_histogram_constructed.h"

/// \file class_histogram.h
/// \brief Histogram

using namespace std;

/// \class histogram
/// \brief Class containing the root histograms
class Histogram : public Histogram_constructed
{
 private :
  Double_t time; ///< Time of acquisition

 public :
  /// \brief Default constructor
  Histogram();
  /// \brief Surdefine Constructor
  Histogram(TString name);
  /// \brief Default destructor
  ~Histogram();
  
  /// \brief Get the run time
  /// \return runtime
  Double_t Get_time() const;

  /// \brief Get the dead time for each detector
  /// \return a vector containing the dead times of each detectors
  vector<Double_t> Get_deadtime() const;
  
  /// \brief Get the relative activity compare to a first histogram, tau is the decay constant
  /// \return The relative activity
  Double_t Get_activity(const Histogram &hist, Double_t lambda) const;

  /// \brief Compute the number of second from the 01/01/1970
  /// /return this number
  time_t Get_instant() const;
};

#endif


