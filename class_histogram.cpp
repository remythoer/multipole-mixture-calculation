#include "class_histogram.h"

using namespace std;

Histogram::Histogram() : filename(' '), file(NULL), time(0.)
{
  for(Int_t i=0;i<4;i++) histo[i]=NULL;
}

Histogram::Histogram(TString name) : filename(name)
{
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
{
  file->Close();
  delete file;
}

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

vector<TH1D*> Histogram::Get_histo() const
{
  vector<TH1D*> histograms;
  for(int i=0;i<4;i++)
    {
      histograms.push_back(histo[i]);
    }  
return histograms;
}

vector<Int_t> Histogram::Get_bin(Int_t bin) const
{
  vector<Int_t> content;
  for(Int_t i=0;i<4;i++) content.push_back(histo[i]->GetBinContent(bin));/// \brief loop on the four detectors
  return content;  
}

Int_t Histogram::Integrate(Int_t detector, Int_t bin, Int_t nb_bin, Int_t bin_g, Int_t bin_d) const
{
  Int_t integral(0);
  Int_t bin_max=bin+nb_bin;
  for(Int_t j=bin;j<bin_max;j++)
    {
       integral=integral+histo[detector]->GetBinContent(j);
    }
  integral=integral-this->Background(detector,bin,nb_bin,bin_g,bin_d); /// \brief Background substraction 
  return integral;
}

Int_t Histogram::Background(Int_t detector, Int_t bin, Int_t nb_bin, Int_t bin_g, Int_t bin_d) const
{
  Int_t background;
  Int_t bin_max=bin+nb_bin;
  Float_t ya = histo[detector]->GetBinContent(bin_g);
  Float_t yb = histo[detector]->GetBinContent(bin_d);
  Float_t xa = histo[detector]->GetBinCenter(bin_g);   
  Float_t xb = histo[detector]->GetBinCenter(bin_d);
  Float_t a = (yb-ya)/(xb-xa);
  Float_t b = yb-a*xb;
  Float_t e_bin_a = histo[detector]->GetBinCenter(bin);  
  Float_t e_bin_b = histo[detector]->GetBinCenter(bin_max);  
  background=a/2.*(e_bin_b*e_bin_b-e_bin_a*e_bin_a)+b*(e_bin_b-e_bin_a);
  return background;
}

Double_t Histogram::Get_time() const
{
  return time;
}

vector<vector<Int_t> > Histogram::Get_peakpos() const
{
  vector<vector<Int_t> > positions;
  for(Int_t i=0;i<4;i++)
    {
     TSpectrum *s = new TSpectrum(20);
     TH1D *d = new TH1D(*histo[i]);
     d->GetXaxis()->SetLimits(1450,1550);
     Int_t nfound = s->Search(histo[i],2,"",0.10);   
     Float_t *xpeaks = s->GetPositionX();
     vector<Int_t> bin;
     for(Int_t p=0;p<nfound;p++)
       { 
	 Float_t xp = xpeaks[p];
	 bin.push_back(histo[i]->GetXaxis()->FindBin(xp));
       } 
     positions.push_back(bin);
     delete s;
     delete d;
    }
  return positions;
}

vector<vector<Double_t> > Histogram::Calibrate() const
{
  vector<vector<Double_t> > cal;
  Int_t peak[3];
  Int_t energy[3];
  energy[0]=511;
  energy[1]=847;
  energy[2]=1238;
  for(Int_t i=0;i<4;i++)
    {
      TSpectrum *s1 = new TSpectrum(20);
      Int_t nfound = s1->Search(histo[i],2,"",0.10);   
      Float_t *xpeaks1 = s1->GetPositionX();
      for(Int_t j=0;j<4;j++)
      	{
      	  if(550 < histo[i]->GetXaxis()->FindBin(xpeaks1[j]) && histo[i]->GetXaxis()->FindBin(xpeaks1[j]) < 650) peak[0]=histo[i]->GetXaxis()->FindBin(xpeaks1[j]);
	  if(950 < histo[i]->GetXaxis()->FindBin(xpeaks1[j]) && histo[i]->GetXaxis()->FindBin(xpeaks1[j]) < 1050) peak[1]=histo[i]->GetXaxis()->FindBin(xpeaks1[j]);
      	  if(1400 < histo[i]->GetXaxis()->FindBin(xpeaks1[j]) && histo[i]->GetXaxis()->FindBin(xpeaks1[j]) < 1600) peak[2]=histo[i]->GetXaxis()->FindBin(xpeaks1[j]);
      	}
      delete s1;
      vector<Double_t> param;
      TGraph *calib = new TGraph(3,peak,energy);
      TF1 *f1=new TF1("f","[0]*x*x+[1]*x+[2]");
      calib->Fit(f1,"QN");
      Double_t par[3];
      f1->GetParameters(&par[0]);
      delete calib;
      delete f1;
      param.push_back(par[0]);
      param.push_back(par[1]);
      param.push_back(par[2]);
      cal.push_back(param);
    }
  return cal;
}

vector<vector<Double_t> > Histogram::CalibInv() const
{
  vector<vector<Double_t> > cal;
  Int_t peak[3];
  Int_t energy[3];
  energy[0]=511;
  energy[1]=847;
  energy[2]=1238;
  for(Int_t i=0;i<4;i++)
    {
      TSpectrum *s1 = new TSpectrum(20);
      Int_t nfound = s1->Search(histo[i],2,"",0.10);   
      Float_t *xpeaks1 = s1->GetPositionX();
      for(Int_t j=0;j<4;j++)
      	{
      	  if(550 < histo[i]->GetXaxis()->FindBin(xpeaks1[j]) && histo[i]->GetXaxis()->FindBin(xpeaks1[j]) < 650) peak[0]=histo[i]->GetXaxis()->FindBin(xpeaks1[j]);
	  if(950 < histo[i]->GetXaxis()->FindBin(xpeaks1[j]) && histo[i]->GetXaxis()->FindBin(xpeaks1[j]) < 1050) peak[1]=histo[i]->GetXaxis()->FindBin(xpeaks1[j]);
      	  if(1400 < histo[i]->GetXaxis()->FindBin(xpeaks1[j]) && histo[i]->GetXaxis()->FindBin(xpeaks1[j]) < 1600) peak[2]=histo[i]->GetXaxis()->FindBin(xpeaks1[j]);
      	}
      delete s1;
      vector<Double_t> param;
      TGraph *calib = new TGraph(3,energy,peak);
      TF1 *f1=new TF1("f","[0]*x*x+[1]*x+[2]");
      calib->Fit(f1,"QN");
      Double_t par[3];
      f1->GetParameters(&par[0]);
      delete calib;
      delete f1;
      param.push_back(par[0]);
      param.push_back(par[1]);
      param.push_back(par[2]);
      cal.push_back(param);
    }
  return cal;
}

void Histogram::temp1(const Histogram &warm, Int_t nbfileWarm, Double_t &temp, Double_t &Errtemp) const
{
  Double_t WW=(((double)this->Integrate(0,1400,17,1390,1427)+(double)this->Integrate(2,1341,13,1331,1364))/((double)this->Integrate(1,1381,10,1371,1401)+(double)this->Integrate(3,1393,10,1383,1413)))*(((double)warm.Integrate(1,1381,10,1371,1401)+(double)warm.Integrate(3,1393,10,1383,1413))/((double)warm.Integrate(0,1400,17,1390,1427)+(double)warm.Integrate(2,1341,13,1331,1364)));
  Double_t ErrWW=sqrt(WW*WW*(1./((double)this->Integrate(0,1400,17,1390,1427)+(double)this->Integrate(2,1341,13,1331,1364))+1./nbfileWarm/((double)warm.Integrate(1,1381,10,1371,1401)+(double)warm.Integrate(3,1393,10,1383,1413))+1./((double)this->Integrate(1,1381,10,1371,1401)+(double)this->Integrate(3,1393,10,1383,1413))+1./nbfileWarm/((double)warm.Integrate(0,1400,17,1390,1427)+(double)warm.Integrate(2,1341,13,1331,1364))));
  double As = TMath::Abs((WW-1.)*100.);
  double ErrAs=ErrWW*100.;
  temp=exp(4.59648-0.481448*log(As)-0.0868518*log(As)*log(As)+0.102367*log(As)*log(As)*log(As)-0.0591571*log(As)*log(As)*log(As)*log(As)+0.0154936*log(As)*log(As)*log(As)*log(As)*log(As)-0.00158925*log(As)*log(As)*log(As)*log(As)*log(As)*log(As));
  Errtemp=ErrAs*(47.7282*exp(log(As)*log(As)*(-0.0868518+0.102367*log(As)-0.0591571*log(As)*log(As)+0.0154936*log(As)*log(As)*log(As)-0.00158925*log(As)*log(As)*log(As)*log(As))))/pow(As,1.48145);
}

void Histogram::temp2(const Histogram &warm, Int_t nbfileWarm, Double_t &temp, Double_t &Errtemp) const
{
  Double_t WW=(((double)this->Integrate(0,1591,17,1581,1618)+(double)this->Integrate(2,1521,15,1511,1546))/((double)this->Integrate(1,1567,11,1557,1588)+(double)this->Integrate(3,1582,10,1572,1602)))*(((double)warm.Integrate(1,1567,11,1557,1588)+(double)warm.Integrate(3,1582,10,1572,1602))/((double)warm.Integrate(0,1591,17,1581,1618)+(double)warm.Integrate(2,1521,15,1511,1546)));
  Double_t ErrWW=sqrt(WW*WW*(1./((double)this->Integrate(0,1591,17,1581,1618)+(double)this->Integrate(2,1521,15,1511,1546))+1./nbfileWarm/((double)warm.Integrate(1,1567,11,1557,1588)+(double)warm.Integrate(3,1582,10,1572,1602))+1./((double)this->Integrate(1,1567,11,1557,1588)+(double)this->Integrate(3,1582,10,1572,1602))+1./nbfileWarm/((double)warm.Integrate(0,1591,17,1581,1618)+(double)warm.Integrate(2,1521,15,1511,1546))));
double As=TMath::Abs((WW-1.)*100.);
double ErrAs=ErrWW*100.;
temp=exp(4.59736-0.480986*log(TMath::Abs(As))-0.0879057*log(TMath::Abs(As))*log(TMath::Abs(As))+0.103411*log(TMath::Abs(As))*log(TMath::Abs(As))*log(TMath::Abs(As))-0.0596485*log(TMath::Abs(As))*log(TMath::Abs(As))*log(TMath::Abs(As))*log(TMath::Abs(As))+0.0156033*log(TMath::Abs(As))*log(TMath::Abs(As))*log(TMath::Abs(As))*log(TMath::Abs(As))*log(TMath::Abs(As))-0.00159848*log(TMath::Abs(As))*log(TMath::Abs(As))*log(TMath::Abs(As))*log(TMath::Abs(As))*log(TMath::Abs(As))*log(TMath::Abs(As)));
 Errtemp=ErrAs*(1./(pow(As,1.48099)))*exp(log(As)*log(As)*(-0.0879057 + 0.103411*log(As)-0.0596485*log(As)*log(As)+0.0156033*log(As)*log(As)*log(As)-0.00159848*log(As)*log(As)*log(As)*log(As)))*(-47.7244-17.4444*log(As)+30.7819*log(As)*log(As)-23.6738*log(As)*log(As)*log(As)+7.74095*log(As)*log(As)*log(As)*log(As)-0.951627*log(As)*log(As)*log(As)*log(As)*log(As));
}
