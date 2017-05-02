#include "class_histogram.h"
#include<dirent.h>
#include<sstream>
#ifndef WIN32
#include<sys/types.h>
#endif

using namespace std;

int main(int argc, char *argv[]) 
{
  TApplication app("appplication",NULL,NULL); 
  Int_t j1(0), j2(0); /// \brief count the number of rootfiles
  ifstream HistoInputs;
  vector<string> peakref;
  vector<vector<Double_t> > valeur_integrale;
  vector<vector<Double_t> > valeur_integrale_cold;
  Double_t time1(0), time2(0);
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
	      HistoInputs.open("./ForRemy/HistoForAnalysis.txt");
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
		      Int_t k=i-detector*23;
		      valeur_integrale[k][detector]=valeur_integrale[k][detector]+Integrale;
		      
		      
		      i++;
		    }
		}
	      else {cout<<"Error file "<<name<<" not open"<<endl;}
	      HistoInputs.close();
	      time1=time1+histo.Get_time();
	      	      
	      cout<<name<<" "<<j1<<" Time="<<histo.Get_time()<<endl;
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
	      HistoInputs.open("./ForRemy/HistoForAnalysis.txt");
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
		      Int_t k=i-detector*23;
		      valeur_integrale_cold[k][detector]=valeur_integrale_cold[k][detector]+Integrale;
		      
		      i++;
		    }
		}
	      else {cout<<"Error file "<<name<<" not open"<<endl;}
	      HistoInputs.close();
	      time2=time2+histo.Get_time();
	      histo.temp1(warm, 1, temp, Errtemp);
	      
	      cout<<name<<" "<<j2<<" T="<<temp<<"+/-"<<Errtemp<<" Time="<<histo.Get_time()<<endl;
	      j2++;
	    }
	}
    }
  closedir(dir);
  cout<<"out"<<endl;
  cout<<"time1="<<time1<<" time2="<<time2<<endl;
  ofstream filewarm;
  filewarm.open("file_warm.dat");
  for(Int_t i=0; i<peakref.size(); i++)
    {
      valeur_integrale[i][0]=valeur_integrale[i][0]/time1;
      valeur_integrale[i][1]=valeur_integrale[i][1]/time1;
      valeur_integrale[i][2]=valeur_integrale[i][2]/time1;
      valeur_integrale[i][3]=valeur_integrale[i][3]/time1;
      
      filewarm<<peakref[i]<<" "<<valeur_integrale[i][0]<<" "<<valeur_integrale[i][1]<<" "<<valeur_integrale[i][2]<<" "<<valeur_integrale[i][3]<<endl;
    }
  filewarm.close();

  ofstream filecold;
  filecold.open("file_cold.dat");
  for(Int_t i=0; i<peakref.size(); i++)
    {
      valeur_integrale_cold[i][0]=valeur_integrale_cold[i][0]/time2;
      valeur_integrale_cold[i][1]=valeur_integrale_cold[i][1]/time2;
      valeur_integrale_cold[i][2]=valeur_integrale_cold[i][2]/time2;
      valeur_integrale_cold[i][3]=valeur_integrale_cold[i][3]/time2;

      filecold<<peakref[i]<<" "<<valeur_integrale_cold[i][0]<<" "<<valeur_integrale_cold[i][1]<<" "<<valeur_integrale_cold[i][2]<<" "<<valeur_integrale_cold[i][3]<<endl;
    }
  filecold.close();

  //app.Run();
  return 0;
}
