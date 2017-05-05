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

  vector<vector<Double_t> > valeur_integrale_sum;
  vector<vector<Double_t> > valeur_integrale_cold_sum;

  Double_t time1[4]={0};
  Double_t time2[4]={0};
  Double_t temp(0), Errtemp(0);
  Histogram warm("./ForRemy/RunForAnalysis/Warm_night_2015-04-21--62067");
  Histogram_constructed sum_warm ("./sum_warm");
  Histogram_constructed sum_cold ("./sum_cold");

  
  struct dirent* openfile = NULL;
  DIR* dir = NULL;

  /// \brief Loop on the warm runs to make the sum

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

	      sum_warm += histo;
	      	     
	      for(Int_t j=0; j<4; j++){time1[j]=time1[j]+histo.Get_time()-histo.Get_deadtime()[j];} /// \brief deadtime of warm runs	      	      
	      cout<<name<<" "<<j1<<" Time="<<histo.Get_time()<<endl;
	      cout<<"calibration detector 3 a="<<histo.Calibrate()[3][0]<<" b="<<histo.Calibrate()[3][1]<<" c="<<histo.Calibrate()[3][2]<<endl;
	      j1++;
	    }
	}
    }
  closedir(dir);

/// \brief Loop on the cold runs to make the sum

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

	      sum_cold+=histo;
	      
	      for(Int_t j=0; j<4; j++){time2[j]=time2[j]+histo.Get_time()-histo.Get_deadtime()[j];}/// \brief deadtime of cold runs	      	      
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
  ofstream filecold;
  ofstream filew;
  
  /// \brief Loop on the peak and calculation of the integrale of each peak with the sum run

  HistoInputs.open("./HistoForAnalysis.txt");
  if(HistoInputs)
    {
       Int_t i(0);
       string ligne;
      while(getline(HistoInputs, ligne))
	{
	  istringstream iss (ligne);
	  Int_t detector, bin_g, bin_d, min_bg, max_bg;
	  string ref;
	  iss >> detector >> bin_g >> bin_d >> min_bg >> max_bg >> ref; 
	  if(detector==0)
	    {
	      peakref.push_back(ref);
	      valeur_integrale_cold_sum.push_back(vector<Double_t>(4,0));
	      valeur_integrale_sum.push_back(vector<Double_t>(4,0));
	    }
	  Int_t Integrale_cold = sum_cold.Integrate(detector,bin_g,bin_d-bin_g,bin_g-max_bg,bin_d+max_bg);
	  Int_t Integrale_warm = sum_warm.Integrate(detector,bin_g,bin_d-bin_g,bin_g-max_bg,bin_d+max_bg);
	  Int_t k=i-detector*25;
	  valeur_integrale_cold_sum[k][detector]=valeur_integrale_cold_sum[k][detector]+Integrale_cold;
	  valeur_integrale_sum[k][detector]=valeur_integrale_sum[k][detector]+Integrale_warm;
	  i++;
	}
    }
  else {cout<<"Error file not open"<<endl;}
  HistoInputs.close();

  /// \brief Output files, ncount warm, ncount cold and w=cold/warm

  filewarm.open("file_warm_sum.dat");
  filecold.open("file_cold_sum.dat");
  filew.open("file_w_sum.dat");
  for(Int_t i=0; i<peakref.size(); i++)
    {
      for(Int_t j=0; j<4; j++)
	{/// \brief time normalisation
	  valeur_integrale_sum[i][j]=valeur_integrale_sum[i][j]/time1[j];
	  valeur_integrale_cold_sum[i][j]=valeur_integrale_cold_sum[i][j]/time2[j];
	}      
      filewarm<<peakref[i]<<" "<<valeur_integrale_sum[i][0]<<" "<<valeur_integrale_sum[i][1]<<" "<<valeur_integrale_sum[i][2]<<" "<<valeur_integrale_sum[i][3]<<endl;
      filecold<<peakref[i]<<" "<<valeur_integrale_cold_sum[i][0]<<" "<<valeur_integrale_cold_sum[i][1]<<" "<<valeur_integrale_cold_sum[i][2]<<" "<<valeur_integrale_cold_sum[i][3]<<endl;
      filew<<peakref[i]<<" "<<valeur_integrale_cold_sum[i][0]/valeur_integrale_sum[i][0]<<" "<<valeur_integrale_cold_sum[i][1]/valeur_integrale_sum[i][1]<<" "<<valeur_integrale_cold_sum[i][2]/valeur_integrale_sum[i][2]<<" "<<valeur_integrale_cold_sum[i][3]/valeur_integrale_sum[i][3]<<endl;
    }
  filewarm.close();
  filecold.close();
  filew.close();


  cout<<"done"<<endl;
  //app.Run();
  return 0;
}
