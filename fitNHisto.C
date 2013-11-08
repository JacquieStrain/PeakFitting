//To run: 
// >>.L fitNHisto.C
// >>fitNHisto()

#include "Fit/Fitter.h"
#include "Fit/BinData.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "Math/WrappedTF1.h"
#include "HFitInterface.h"

// example of fitting N histograms with the same function w/some common and some different parameters
// adapted from sample script I found on a ROOT forum

char *modGaus(int par_start, int par_a, int par_b, int par_g){
	char* function = new char[100];
	int amp=par_start;
	int mu=par_start+1;
	sprintf(function,"[%d]*exp(-0.5*pow(x-[%d],2.)/(pow([%d],2.)+[%d]*pow([%d],2.)+pow([%d]*[%d],2.)))",
					amp, mu, par_a, mu, par_b, mu, par_g);
	return function;
};

void fitNHisto() { 

//Double_t alpha = 0.38;
//Double_t beta = 0.03;
//Double_t gamma = 0.01;

//Need to learn how to use TObject to make this more flexible for n peaks
TF1 *f1 = new TF1("f1","gaus(0)",5,13);
f1->SetParameters(1, 9, 0.40); //"exact" sigma = 0.400749...
TF1 *f2 = new TF1("f2","gaus(0)",19,29);
f2->SetParameters(1, 24, 0.47); //"exact" sigma = 0.472864...
TF1 *f3 = new TF1("f3","gaus(0)",36,48);
f3->SetParameters(1, 42, 0.59); //"exact" sigma = 0.598832...
TF1 *f4 = new TF1("f4","gaus(0)",55,71);
f4->SetParameters(1, 63, 0.76); //"exact" sigma = 0.773305...

TH1F *h1 = new TH1F("h1","pk1",100,5,13);
h1->FillRandom("f1",2000);
TH1F *h2 = new TH1F("h2","pk2",100,19,29);
h2->FillRandom("f2",2000);
TH1F *h3 = new TH1F("h3","pk3",100,36,48);
h3->FillRandom("f3",2000);
TH1F *h4 = new TH1F("h4","pk4",100,55,71);
h4->FillRandom("f4",2000);

ROOT::Fit::BinData data; 
ROOT::Fit::FillData(data, h1); 
ROOT::Fit::FillData(data, h2);
ROOT::Fit::FillData(data, h3);
ROOT::Fit::FillData(data, h4); 
cout << "data size is " << data.Size() << endl;

int noPeaks = 4;
char function[1000];
char gaussian[100];
int aa=2, bb=3, gg=4;
for(int i=0; i<noPeaks; i++){
	if(i==0){
		sprintf(gaussian,"%s",modGaus(i,aa,bb,gg));
		strcpy(function,gaussian);
		}
	else{
		sprintf(gaussian,"+%s",modGaus(2*i+3,aa,bb,gg));
		strcat(function,gaussian);
		}
	}

TF1 *f0 = new TF1("f0",function,0,100);
f0->SetParameters(80, 9, 0.4, 0.03, 0.01, 80, 24, 80, 42, 80, 63);

ROOT::Math::WrappedTF1 wf(*f0);

ROOT::Fit::Fitter fitter;
fitter.SetFunction(wf);

fitter.Fit(data);
ROOT::Fit::FitResult result = fitter.Result();
result.Print(std::cout);

TCanvas *c1 = new TCanvas();
gStyle->SetOptStat(0);
c1->SetFillColor(0);
c1->SetFrameBorderMode(0);
int noDiv;
if(noPeaks%2 != 0){noDiv=(noPeaks+1)/2;}
else{noDiv=noPeaks/2;}
c1->Divide(noDiv,2);

f0->SetLineColor(6);
for(int i=0; i<noPeaks; i++){
	c1->cd(i+1);
	gPad->SetFillColor(0);
	gPad->SetFrameBorderMode(0);
	if(i==0){h1->Draw();} //Again, need to make this more flexible for n peaks with TObject(?)
	if(i==1){h2->Draw();}
	if(i==2){h3->Draw();}
	if(i==3){h3->Draw();}
	f0->Draw("Psame");
	}

c1->WaitPrimitive();
delete c1;
}

int main() { 

   fitNHisto(); 
}
