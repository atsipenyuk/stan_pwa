c = new TCanvas("Canvas","Plotting Canvas");
gStyle->SetPadTickY(1);
gStyle->SetOptStat(0);
t->Draw("y","","");
TH1F *h1 = (TH1F*)htemp->Clone("h1");
h1->SetTitle(";sigma;Number of samples");
h1->SetLabelSize(0.047, "X")
h1->SetLabelSize(0.047, "Y")
