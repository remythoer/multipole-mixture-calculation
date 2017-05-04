#include "class_histogram.h"

using namespace std;

Histogram::Histogram() : Histogram_constructed(), time(0.)
{}

Histogram::Histogram(TString name) 
{
  filename = name;
  file = new TFile(name+".root");
  histo[0] = (TH1D*)file->Get("Detector_0");
  histo[1] = (TH1D*)file->Get("Detector_1");
  histo[2] = (TH1D*)file->Get("Detector_2");
  histo[3] = (TH1D*)file->Get("Detector_3");
  ifstream mpafile;
  mpafile.open(name+".mpa");
  TString realtime;
  for (int j=0; j< 52; j++)
    {
      mpafile >> realtime;
    }
  realtime=realtime.Remove(0,9);
  time=atof(realtime);
}

Histogram::~Histogram()
{}

vector<Double_t> Histogram::Get_deadtime() const
{
  TString S_deadtime;
  Double_t D_deadtime;
  vector<Double_t> V_deadtime;
  ifstream mpafile;
  mpafile.open(filename+".mpa");
  for (int j=0; j<54; j++) {mpafile>>S_deadtime;}
  S_deadtime=S_deadtime.Remove(0,9);
  D_deadtime=atof(S_deadtime);
  V_deadtime.push_back(time-D_deadtime); /// \brief deadtime for detector 0
  for(int j=0; j<32; j++){mpafile>>S_deadtime;}
  S_deadtime=S_deadtime.Remove(0,9);
  D_deadtime=atof(S_deadtime);
  V_deadtime.push_back(time-D_deadtime); /// \brief deadtime for detector 1
  for(int j=0; j<32; j++){mpafile>>S_deadtime;}
  S_deadtime=S_deadtime.Remove(0,9);
  D_deadtime=atof(S_deadtime);
  V_deadtime.push_back(time-D_deadtime); /// \brief deadtime for detector 2
  for(int j=0; j<32; j++){mpafile>>S_deadtime;}
  S_deadtime=S_deadtime.Remove(0,9);
  D_deadtime=atof(S_deadtime);
  V_deadtime.push_back(time-D_deadtime); /// \brief deadtime for detector 3
  mpafile.close();
  return V_deadtime;
}

Double_t Histogram::Get_time() const
{
  return time;
}


