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

Double_t Histogram::Get_activity(const Histogram &hist, Double_t lambda) const
{
  time_t first=hist.Get_instant();
  time_t second=this->Get_instant();
  double delta_t=difftime(second, first);
  
  return exp(-delta_t*lambda);  
}

time_t Histogram::Get_instant() const
{
  time_t temps;
  struct tm date;
  TString instant;
  TString temp_instant;
  ifstream mpafile;
  double month;
  double day;
  double year;
  double hour;
  double minute;
  double second;
  mpafile.open(filename+".mpa");

  for(Int_t i=0; i<47; i++){mpafile>>instant;}
  temp_instant=instant;
  temp_instant.Remove(2,8);
  month=atof(temp_instant);
  temp_instant=instant;
  temp_instant.Remove(0,3);
  temp_instant.Remove(2,5);
  day=atof(temp_instant);
  temp_instant=instant;
  temp_instant.Remove(0,6);
  year=atof(temp_instant);

  mpafile>>instant;
  temp_instant=instant;
  temp_instant.Remove(2,6);
  hour=atof(temp_instant);
  temp_instant=instant;
  temp_instant.Remove(0,3);
  temp_instant.Remove(2,3);
  minute=atof(temp_instant);
  temp_instant=instant;
  temp_instant.Remove(0,6);
  second=atof(temp_instant);

  date.tm_year=year-1900;
  date.tm_mon=month-1;
  date.tm_mday=day;
  date.tm_hour=hour;
  date.tm_min=minute;
  date.tm_sec=second;
  
  temps=mktime(&date);
    
  return temps;
}
