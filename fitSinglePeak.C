//To run:
//>>.L fitSinglePeak.C
//>>fitSinglePeak([name of file with histogram], [name of histogram], [start of range to fit], [end of range to fit])

//There are a few "TCanvas->WaitPrimitive()"s so you can digest the many outputs. 
//What that means is the program will hang until you click on the TCanvas and hit spacebar.

#include "TFile.h"
#include "Riostream.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TF1.h"
#include "TFitResult.h"

void fitSinglePeak(char filename[130], char histname[30], double start, double end){

TCanvas *c1 = new TCanvas();
gStyle->SetOptStat(0);
c1->SetFillColor(0);
c1->SetFrameBorderMode(0);

//Open the histogram and display the range that is about to be fitted.
TFile *f = new TFile(filename);
TH1D *spectrum = f->Get(histname);
gStyle->SetErrorX(0);
spectrum->GetXaxis()->SetRangeUser(start,end);
spectrum->SetMarkerStyle(20);
spectrum->Draw("E1");
c1->WaitPrimitive();
char function[10000];

//Fit with a quadratic bkgd + Gauss + step + tail
{
char quad_bkgd[50];
strcpy(quad_bkgd, "[0]+x*[1]+[2]*pow(x,2.)");
strcpy(function,quad_bkgd);
char gauss[100];
strcpy(gauss, "[3]*(1-[7])*exp(-0.5*pow((x-[4])/[5],2.))/([5]*sqrt(2.*TMath::Pi()))");
strcat(function,"+");
strcat(function,gauss);
char step[100];
strcpy(step, "TMath::Abs([6])*[3]*TMath::Erfc((x-[4])/([5]*sqrt(2.)))");
strcat(function,"+");
strcat(function,step);
char tail[400];
strcpy(tail,"[7]*[3]*exp((x-[4])/[8])*TMath::Erfc((x-[4])/(sqrt(2.0)*[5])+[5]/(sqrt(2.0)*[8]))/(2.*[8]*exp(-([5]/(sqrt(2.)*[8]))^2))");
strcat(function,"+");
strcat(function,tail);

//~*~*~*~*
/*
Initial guesses for fit. Change as needed.
If you are having problems getting the fit to converge, sometimes it helps to fix 
par 1 (slope of bkgd) and/or par 2 (x-squared term of bkgd) to zero.
Another thing that might help is either fixing par 7 (Htail) or par 8 (tau) to a reasonable value. 
(This might help because the 2 parameters are very correlated.)
*/
double slope_approx = (spectrum->GetBinContent(spectrum->GetXaxis()->FindBin(end)) - spectrum->GetBinContent(spectrum->GetXaxis()->FindBin(end)))/(end-start);
double yint_approx = 0.0;
double xterm_approx = 0.0;
double sigma_approx = (end-start)/(2.*10.);
double peak_approx = spectrum->GetXaxis()->GetBinCenter(spectrum->GetMaximumBin());
double area_approx = sigma_approx*sqrt(2.*TMath::Pi())*spectrum->GetBinContent(spectrum->GetMaximumBin());
double Hstep_approx = 1.*pow(10.,-4.);
double Htail_approx = 0.2;
double tau_approx = sigma_approx*1.5;
//~*~*~*~*

printf("*******************************************************************\n");
printf("Y-int approx = %7.5f \n", yint_approx);
printf("Slope approx = %7.5f \n", slope_approx);
printf("x^2 term approx = %7.5f \n", xterm_approx);
printf("Area approx = %7.5f \n", area_approx);
printf("Mu (gaus centroid) approx = %7.5f \n", peak_approx);
printf("Sigma approx = %7.5f \n", sigma_approx);
printf("Hstep approx= %7.7e \n", Hstep_approx);
printf("Htail approx = %7.5f \n", Htail_approx);
printf("Tau approx = %7.5f \n", tau_approx);
printf("*******************************************************************\n");

printf("\tFitting with a gaussian, step, and tail (plus a quadratic background)\n");
TF1 *fit = new TF1("fit",function,start,end);
fit->SetNpx(100000.0);
fit->SetParameter(0,yint_approx);
fit->SetParameter(1,slope_approx);
fit->SetParameter(2,xterm_approx);
fit->SetParameter(3,area_approx);
fit->SetParameter(4,peak_approx);
fit->SetParameter(5,sigma_approx);
fit->SetParameter(6,Hstep_approx);
fit->SetParameter(7,Htail_approx);
fit->SetParameter(8,tau_approx);
TFitResultPtr result = spectrum->Fit("fit","SRN"); //If fitting low stats, should add option "L" to fit.
TMatrixDSym cov = result->GetCovarianceMatrix();
result->Print("V");
fit->SetLineColor(2); //red
fit->Draw("Psame");
if(fit->GetNDF() != 0){
     printf("ChiSq: %.2f \n", fit->GetChisquare());
     printf("ChiSq/NDF: %.4f \n", fit->GetChisquare()/fit->GetNDF());
     printf("p value: %.4f \n", fit->GetProb());
     }
c1->WaitPrimitive();
}

//Find the centroid, sigma, area, and then plot the residuals.
{
TF1 *pk0 = new TF1("pk0", "[0]+[1]*x+[2]*pow(x,2.)",start,end);
TF1 *pk1 = new TF1("pk1", "[0]*(1.-[3])*exp(-0.5*((x-[1])/[2])^2)/([2]*sqrt(2.*TMath::Pi()))",start,end);
TF1 *pk2 = new TF1("pk2", "[0]*TMath::Abs([3])*TMath::Erfc((x-[1])/([2]*sqrt(2.0)))",start,end);
TF1 *pk3 = new TF1("pk3", "[0]*[4]*exp((x-[1])/[3])*TMath::Erfc((x-[1])/(sqrt(2.0)*[2])+[2]/(sqrt(2.0)*[3]))/(2.*[3]*exp(-([2]/(sqrt(2.)*[3]))^2))",start,end);

double MaxXFit, HalfMaxXFit, lowFWHM, highFWHM, FWHM, dFWHM;
double PeakCentroid, dPeakCentroid;

TH1D *residual = spectrum->Clone();
TH1D *sigmaResidual = spectrum->Clone();

//Assign variables
{
double Yint = fit->GetParameter(0);
double Slope = fit->GetParameter(1);
double X2term = fit->GetParameter(2);
double Area = fit->GetParameter(3);
double Mu = fit->GetParameter(4);
double Sigma = fit->GetParameter(5);
double N2 = fit->GetParameter(6);
double N3rel = fit->GetParameter(7);
double Tau = fit->GetParameter(8);
double dYint = fit->GetParError(0);
double dSlope = fit->GetParError(1);
double dX2term = fit->GetParError(2);
double dArea = fit->GetParError(3);
double dMu = fit->GetParError(4);
double dSigma = fit->GetParError(5);
double dN2 = fit->GetParError(6);
double dN3rel = fit->GetParError(7);
double dTau = fit->GetParError(8);
double GaussCentroid, TailCentroid, AreaGauss, AreaTail;
}

//pk4 is the gaussian+tail (i.e. just the signal; no linear bkgd or step)
//Use to find FWHM. Do not draw.
{
TF1 *pk4 = new TF1("pk4", "[0]*(1.-[4])*exp(-0.5*((x-[1])/[2])^2)/([2]*sqrt(2.*TMath::Pi())) + [0]*[4]*exp((x-[1])/[3])*TMath::Erfc((x-[1])/(sqrt(2.0)*[2])+[2]/(sqrt(2.0)*[3]))/(2.*[3]*exp(-([2]/(sqrt(2.)*[3]))^2))",start,end);
pk4->SetNpx(100000.0);
pk4->SetParameter(0,Area);
pk4->SetParameter(1,Mu);
pk4->SetParameter(2,Sigma);
pk4->SetParameter(3,Tau);
pk4->SetParameter(4,N3rel);
    
MaxXFit = pk4->GetMaximumX();
HalfMaxXFit = pk4->Eval(MaxXFit)/2.;
lowFWHM = pk4->GetX(HalfMaxXFit, MaxXFit-5.*Sigma, MaxXFit);
highFWHM = pk4->GetX(HalfMaxXFit, MaxXFit, MaxXFit+5.*Sigma);
FWHM = highFWHM - lowFWHM;
dFWHM = dSigma*FWHM/Sigma;
printf("FWHM of the FIT: %10.3f +/- %10.3f \n", FWHM, dFWHM);
printf("*******************************************************************\n");
}

//pk0 is the quadratic background.
//Draw in color: "salmon" (ROOT color 50)
{
pk0->SetParameter(0,Yint);
pk0->SetParameter(1,Slope);
pk0->SetParameter(2,X2term);
pk0->SetLineColor(50); //salmon
pk0->Draw("Psame");
}

//pk1 is the gaussian.
//Find GAUSSIAN centroid and area. Draw in color "green" (ROOT color 30)
{
pk1->SetParameter(0,Area);
pk1->SetParameter(1,Mu);
pk1->SetParameter(2,Sigma);
pk1->SetParameter(3,N3rel);
pk1->SetLineColor(30); //green
pk1->Draw("Psame");
GaussCentroid = Mu;
printf("GAUSSIAN centroid: %15.3f +/- %15.3f \n", GaussCentroid, dMu);
AreaGauss = 1.-N3rel;
printf("Relative Area of GAUSSIAN: %10.3f \n", AreaGauss);
printf("------------------------\n");
}

//pk2 is the step.
//Draw in color: "gray" (ROOT color 39)
{
pk2->SetParameter(0,Area);
pk2->SetParameter(3,N2);
pk2->SetParameter(1,Mu);
pk2->SetParameter(2,Sigma);
pk2->SetLineColor(39); //gray
pk2->Draw("Psame");         
}

//pk3 is the tail.
//Find TAIL centroid and area. Draw in color "brown" (ROOT color 28)
{
pk3->SetParameter(0,Area);
pk3->SetParameter(4,N3rel);
pk3->SetParameter(1,Mu);
pk3->SetParameter(2,Sigma);
pk3->SetParameter(3,Tau);
pk3->SetLineColor(28); //brown
pk3->Draw("Psame");
TailCentroid = GaussCentroid - Tau;
printf("TAIL centroid: %15.3f \n", TailCentroid);
AreaTail = N3rel;
printf("Relative Area of TAIL: %10.3f \n", AreaTail);
printf("------------------------\n");
}

//Find Peak Centroid.
{
printf("*******************************************************************\n");
PeakCentroid = GaussCentroid - (Tau*N3rel);
//partial derivatives of PeakCentroid wrt its variables
double dPCdM = 1.;
double dPCdT = -N3rel;
double dPCdN3rel = -Tau;
//uncertainty of PeakCentroid
dPeakCentroid = sqrt( pow(dPCdM*dMu,2.0) + pow(dPCdT*dTau,2.0) + pow(dPCdN3rel*dN3rel,2.0)
                    + 2.0*( cov[4][7]*dPCdM*dPCdN3rel + cov[4][8]*dPCdM*dPCdT
                              + cov[7][8]*dPCdN3rel*dPCdT) );
printf("Peak Centroid: %15.3f +/- %15.3f \n", PeakCentroid, dPeakCentroid);
printf("*******************************************************************\n"); 
}

c1->Update();
c1->WaitPrimitive();

//Calculate the variance
//Plot the spectrum, the peak centroid, the sigma, and the area
{
double Variance = pow(Sigma,2.) + pow(Tau*N3rel,2.);
double dVdSigma = 2.*Sigma;
double dVdTau = 2.*Tau*pow(N3rel,2.);
double dVdN3rel = 2.*N3rel*pow(Tau,2.);
double dVariance = sqrt( pow(dVdSigma*dSigma,2.) + pow(dVdTau*dTau,2.) + pow(dVdN3rel*dN3rel,2.)
						+ 2.*cov[5][8]*dVdSigma*dVdTau + 2.*cov[5][7]*dVdSigma*dVdN3rel + 2.*cov[7][8]*dVdTau*dVdN3rel);
double SqrtVar = sqrt(Variance);
double dSqrtVar = dVariance/(2.*SqrtVar);

TCanvas *c4 = new TCanvas();
gStyle->SetOptStat(0);
c4->SetFillColor(0);
c4->SetFrameBorderMode(0);

spectrum->Draw();

double n3sigma = PeakCentroid-(3.*SqrtVar);
double p3sigma = PeakCentroid+(3.*SqrtVar);
double n5sigma = PeakCentroid-(5.*SqrtVar);
double p5sigma = PeakCentroid+(5.*SqrtVar);

TH1D *signal5s = spectrum->Clone();
signal5s->GetXaxis()->SetRangeUser(n5sigma,p5sigma);
signal5s->SetFillColor(61);
signal5s->Draw("same");

TH1D *signal3s = spectrum->Clone();
signal3s->GetXaxis()->SetRangeUser(n3sigma,p3sigma);
signal3s->SetFillColor(65);
signal3s->Draw("same");

TH1D *bkgdMod = spectrum->Clone();
double bkgd_cont;
for(int i=0; i<bkgdMod->GetNbinsX(); i++){
	bkgd_cont = pk0->Eval(bkgdMod->GetBinCenter(i))+pk2->Eval(bkgdMod->GetBinCenter(i));
	bkgdMod->SetBinContent(i,bkgd_cont);
	}
bkgdMod->SetFillColor(13);
bkgdMod->SetFillStyle(3001);
bkgdMod->Draw("same");

fit->SetLineColor(1);
fit->Draw("Psame");

printf("sqrt(Variance): %10.3f +/- %10.3f \n", SqrtVar, dSqrtVar);
printf("10 sqrt(Variance) fit range = %15.1f to %15.1f \n", PeakCentroid-10.*SqrtVar, PeakCentroid+10.*SqrtVar);
double BkgdArea3 = (pk0->Integral(n3sigma,p3sigma) + pk2->Integral(n3sigma,p3sigma))/spectrum->GetBinWidth(spectrum->GetXaxis()->FindBin(PeakCentroid));
double BkgdArea5 = (pk0->Integral(n5sigma,p5sigma) + pk2->Integral(n5sigma,p5sigma))/spectrum->GetBinWidth(spectrum->GetXaxis()->FindBin(PeakCentroid));
double TotalArea3 = spectrum->Integral(spectrum->GetXaxis()->FindBin(n3sigma),spectrum->GetXaxis()->FindBin(p3sigma));
double TotalArea5 = spectrum->Integral(spectrum->GetXaxis()->FindBin(n5sigma),spectrum->GetXaxis()->FindBin(p5sigma));
double SignalArea3 = TotalArea3-BkgdArea3;
double SignalArea5 = TotalArea5-BkgdArea5;
printf("*******************************************************************\n");
printf("+/- 3 sigma \n");
printf("Bkgd Area: %10.3f \n", BkgdArea3);
printf("Total Area: %10.3f \n", TotalArea3);
printf("Signal Area: %10.3f +/- %10.3f \n",SignalArea3,pow(SignalArea3,0.5));
printf("*******************************************************************\n");

c4->Update();
c4->WaitPrimitive();
}

//Residual Plots
{
TCanvas *c2 = new TCanvas();
gStyle->SetOptStat(0);
c2->SetFillColor(0);
c2->SetFrameBorderMode(0);

TPad *c2b = new TPad("c2b", "c2b",0.01,0.01,0.99,0.49);
c2b->Draw();
c2b->cd();
c2b->SetFillColor(0);
c2b->SetFrameBorderMode(0);

c2->cd();
TPad *c2a = new TPad("c2a", "c2a",0.01,0.51,0.99,0.99);
c2a->Draw();
c2a->cd();
c2a->SetFillColor(0);
c2a->SetFrameBorderMode(0);

spectrum->Draw("E1");
fit->SetLineColor(4); //bright blue
fit->Draw("Psame");

//plot the residual: data-fit
{
c2b->cd();
gPad->SetGrid();
double aval;
for(int i=spectrum->GetXaxis()->FindBin(start); i<(spectrum->GetXaxis()->FindBin(end)+1); i++){
     aval = spectrum->GetXaxis()->GetBinCenter(i);
	if(spectrum->GetBinError(i)==0){
		residual->SetBinError(i, 1.);
		}
	else{
		residual->SetBinError(i, spectrum->GetBinError(i));
		}
     residual->SetBinContent(i, spectrum->GetBinContent(i)-fit->Eval(aval));
     aval = aval + spectrum->GetBinWidth(i);
     }
residual->GetYaxis() ->SetTitle("Residual (data - fit)");
residual->Draw("E0");
c2->Update();
c2->WaitPrimitive();

TCanvas *c3 = new TCanvas();
gStyle->SetOptStat(0);
c3->SetFillColor(0);
c3->SetFrameBorderMode(0);

TPad *c3b = new TPad("c3b", "c3b",0.01,0.01,0.99,0.49);
c3b->Draw();
c3b->cd();
c3b->SetFillColor(0);
c3b->SetFrameBorderMode(0);

c3->cd();
TPad *c3a = new TPad("c3a", "c3a",0.01,0.51,0.99,0.99);
c3a->Draw();
c3a->cd();
c3a->SetFillColor(0);
c3a->SetFrameBorderMode(0);

spectrum->Draw("E1");
fit->SetLineColor(4); //bright blue
fit->Draw("Psame");
}

//plot the residual in terms of sigma: data-fit/(data error)
{
c3b->cd();
gPad->SetGrid();
double sigmaE;
for(int i=spectrum->GetXaxis()->FindBin(start); i<(spectrum->GetXaxis()->FindBin(end)+1); i++){
     aval = spectrum->GetXaxis()->GetBinCenter(i);
     if(spectrum->GetBinError(i)==0){
	 	sigmaE = (spectrum->GetBinContent(i)-fit->Eval(aval))/1.;
		}
	 else{
	 	sigmaE = (spectrum->GetBinContent(i)-fit->Eval(aval))/spectrum->GetBinError(i);
	 	}
	 sigmaResidual->SetBinContent(i, sigmaE);
     sigmaResidual->SetBinError(i,1.0);
     aval = aval + spectrum->GetBinWidth(i);
     }
sigmaResidual->GetYaxis() ->SetTitle("Residual in terms of sigma (data - fit)/(data error)");
sigmaResidual->Draw("E0");
c3->Update();
c3->WaitPrimitive();
}
}
}
    
f->Close();
c1->Close();
c2->Close();
c3->Close();
c4->Close();
delete fit;
delete pk0;
delete pk1;
delete pk2;
delete pk3;
delete pk4;
}
