#include "class_histogram.h"
#include<dirent.h>
#include<sstream>
#ifndef WIN32
#include<sys/types.h>
#endif

using namespace std;

int main(int argc, char *argv[]) 
{
  //TApplication app("appplication",NULL,NULL); 
  Int_t j1(0), j2(0); /// \brief count the number of rootfiles
  ifstream HistoInputs;
  vector<string> peakref;
  vector<vector<Double_t> > valeur_integrale;
  vector<vector<Double_t> > valeur_integrale_cold;
  Double_t time1[4]={0};
  Double_t time2[4]={0};
  Double_t temp(0), Errtemp(0);
  Histogram warm("./ForRemy/RunForAnalysis/Warm_night_2015-04-21--62067");
  
  struct dirent* openfile = NULL;
  DIR* dir = NULL;
  dir = opendir("./ForRemy/RunForAnalysis/Warm"); /// \brief open directory
  if(dir==NULL){cout<<"error, directory not open"<<endl;}
  else
    { /// \brief if directory is open
      while ((openfile = readdir(dir)) != NULL) /// \brief loop on all files
	{
	  string name = openfile->d_name;    
	  if (name.compare(".") != 0 && name.compare("..") != 0)
	    { 
	      string::iterator it;
	      for(it=name.begin();it!=name.end();it++){if(*it=='.') break;}
	      it++;
	      if(*it!='r')continue;
	      name="./ForRemy/RunForAnalysis/"+name;
	      Histogram histo(name.erase(name.size()-5,5));
	      
	      Int_t i(0); /// \brief count the number of peaks
	      string ligne;
	      HistoInputs.open("./HistoForAnalysis.txt");
	      if(HistoInputs)
		{ 
		  while(getline(HistoInputs, ligne))
		    {
		      istringstream iss (ligne);
		      Int_t detector, bin_g, bin_d, min_bg, max_bg;
		      string ref;
		      iss >> detector >> bin_g >> bin_d >> min_bg >> max_bg >> ref; 
		      if(j1==0 && detector==0)
			{
			  peakref.push_back(ref);
			  valeur_integrale.push_back(vector<Double_t>(4,0));
			}
		      Int_t Integrale = histo.Integrate(detector,bin_g,bin_d-bin_g,bin_g-max_bg,bin_d+max_bg);
		      Int_t k=i-detector*25;
		      valeur_integrale[k][detector]=valeur_integrale[k][detector]+Integrale;
		      
		      
		      i++;
		    }
		}
	      else {cout<<"Error file "<<name<<" not open"<<endl;}
	      HistoInputs.close();
	      for(Int_t j=0; j<3; j++){time1[j]=time1[j]+histo.Get_time()-histo.Get_deadtime()[j];}	      	      
	      cout<<name<<" "<<j1<<" Time="<<histo.Get_time()<<endl;
	      cout<<"calibration detector 3 a="<<histo.Calibrate()[3][0]<<" b="<<histo.Calibrate()[3][1]<<" c="<<histo.Calibrate()[3][2]<<endl;
	      j1++;
	    }
	}
    }
  closedir(dir);

dir = opendir("./ForRemy/RunForAnalysis/Cold"); /// \brief open directory
  if(dir==NULL){cout<<"error, directory not open"<<endl;}
  else
    { /// \brief if directory is open
      while ((openfile = readdir(dir)) != NULL) /// \brief loop on all files
	{
	  string name = openfile->d_name;    
	  if (name.compare(".") != 0 && name.compare("..") != 0)
	    { 
	      string::iterator it;
	      for(it=name.begin();it!=name.end();it++){if(*it=='.') break;}
	      it++;
	      if(*it!='r')continue;
	      name="./ForRemy/RunForAnalysis/"+name;
	      Histogram histo(name.erase(name.size()-5,5));
	      
	      Int_t i(0); /// \brief count the number of peaks
	      string ligne;
	      HistoInputs.open("./HistoForAnalysis.txt");
	      if(HistoInputs)
		{
		  while(getline(HistoInputs, ligne))
		    {
		      istringstream iss (ligne);
		      Int_t detector, bin_g, bin_d, min_bg, max_bg;
		      string ref;
		      iss >> detector >> bin_g >> bin_d >> min_bg >> max_bg >> ref; 
		      if(j2==0 && detector==0)
			{
			  valeur_integrale_cold.push_back(vector<Double_t>(4,0));
			}
		      Int_t Integrale = histo.Integrate(detector,bin_g,bin_d-bin_g,bin_g-max_bg,bin_d+max_bg);
		      Int_t k=i-detector*25;
		      valeur_integrale_cold[k][detector]=valeur_integrale_cold[k][detector]+Integrale;
		      
		      i++;
		    }
		}
	      else {cout<<"Error file "<<name<<" not open"<<endl;}
	      HistoInputs.close();
	      for(Int_t j=0; j<3; j++){time2[j]=time2[j]+histo.Get_time()-histo.Get_deadtime()[j];}
	      histo.temp1(warm, 1, temp, Errtemp);
	      
	      cout<<name<<" "<<j2<<" T="<<temp<<"+/-"<<Errtemp<<" Time="<<histo.Get_time()<<endl;
	      cout<<"calibInv detector 3 a="<<histo.CalibInv()[3][0]<<" b="<<histo.CalibInv()[3][1]<<" c="<<histo.CalibInv()[3][2]<<endl;
	      j2++;
	    }
	}
    }
  closedir(dir);
  cout<<"out"<<endl;
  ofstream filewarm;
  filewarm.open("file_warm.dat");
  for(Int_t i=0; i<peakref.size(); i++)
    {
      for(Int_t j=0; j<3; j++){valeur_integrale[i][j]=valeur_integrale[i][j]/time1[j];}      
      filewarm<<peakref[i]<<" "<<valeur_integrale[i][0]<<" "<<valeur_integrale[i][1]<<" "<<valeur_integrale[i][2]<<" "<<valeur_integrale[i][3]<<endl;
    }
  filewarm.close();

  ofstream filecold;
  filecold.open("file_cold.dat");
  for(Int_t i=0; i<peakref.size(); i++)
    {
      for(Int_t j=0; j<3; j++){valeur_integrale_cold[i][j]=valeur_integrale_cold[i][j]/time2[j];}
      filecold<<peakref[i]<<" "<<valeur_integrale_cold[i][0]<<" "<<valeur_integrale_cold[i][1]<<" "<<valeur_integrale_cold[i][2]<<" "<<valeur_integrale_cold[i][3]<<endl;
    }
  filecold.close();
  
  cout<<"done"<<endl;
  //app.Run();
  return 0;
}
