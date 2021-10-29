// Ratio CT_{solid target}/CT_{liquid target} for each method and variable
void CanvasPartition(TCanvas *C,const Int_t Nx = 2,const Int_t Ny = 2,
                     Float_t lMargin = 0.15, Float_t rMargin = 0.05,
                     Float_t bMargin = 0.15, Float_t tMargin = 0.05);

void canvasCT_EX2(){
    TString target[] = {"D", "C", "Fe", "Pb"};
    TString method[] = {"SM", "CSMV"};
    TString var[] = {"Zh", "Pt2", "PhiPQ"};
    TString varxtitle[] = {"Z_{h}", "P_{t}^{2} [GeV^{2}]", "#phi_{PQ} [deg]"};

    // Open files
    std::vector<TFile*> input_Vec;

    int NtargetTS = sizeof(target)/sizeof(target[0]);
    int NmethodTS = sizeof(method)/sizeof(method[0]);
    int NvarTS = sizeof(var)/sizeof(var[0]);

    for (int imeth=0; imeth<NmethodTS; imeth++){
        for (int itarg=0; itarg<NtargetTS; itarg++){
            for (int ivar=0; ivar<NvarTS; ivar++){
                TString in_title = "ClosTest"+method[imeth]+"_"+target[itarg]+"1_EX_"+var[ivar]+".root";
                TFile *input = TFile::Open(in_title,"READ");
                input_Vec.push_back(input);
            }
        }
    }

    // Canvas options and variables
    Color_t color[3] = {kViolet+10, kRed, kBlack};
    float Nu_limits[] = {2.2, 3.2, 3.7, 4.2};
    float Q2_limits[] = {1.0, 1.3, 1.8, 4.1};

    gStyle->SetErrorX(0);
    gStyle->SetTitleStyle(0);
    gStyle->SetOptTitle(1);
    gStyle->SetMarkerStyle(7);

    // Fill canvases
    std::vector<TCanvas*> c_Vec;

    int Ncanvas = NmethodTS*NvarTS;
    for (int i=0; i<Ncanvas; i++){
        TCanvas *cratio = new TCanvas("cratio"+method[i/NvarTS]+"_"+var[i%NvarTS],"CT ratio "+method[i/NvarTS]+var[i%NvarTS],
                                      1200,600);

        // Number of PADS
        const Int_t Nx = 3;
        const Int_t Ny = 3;

        // Margins
        Float_t lMargin = 0.08;
        Float_t rMargin = 0.05;
        Float_t bMargin = 0.10;
        Float_t tMargin = 0.05;

        // Canvas setup
	    CanvasPartition(cratio,Nx,Ny,lMargin,rMargin,bMargin,tMargin);

        TPad *pad[Nx][Ny];

        for (int itarg=0; itarg<NtargetTS; itarg++){
            if (itarg==0) continue;
            for (int iNubin=0; iNubin<3; iNubin++){
                cratio->cd(0);

                pad[itarg-1][iNubin] = (TPad*) gROOT->FindObject(Form("pad_%i_%i",itarg-1,2-iNubin));
                pad[itarg-1][iNubin]->Draw();
                pad[itarg-1][iNubin]->cd();

                TString htitle = ";"+varxtitle[i%NvarTS]+";CT_{solid}/CT_{D}";

                TLegend *legend = new TLegend(0.01,0.05,0.70,0.30);
                legend->SetBorderSize(0);
                legend->SetNColumns(3);

                THStack *hratioStack = new THStack(Form("hratioStack"+method[i/NvarTS]+"_"+var[i%NvarTS]+"_"+target[itarg]+"%i",iNubin),htitle);
                double ratiomin = 0.5001; // 0.5001;
                double ratiomax = 1.4999; // 1.4999;
                hratioStack->SetMinimum(ratiomin);
                hratioStack->SetMaximum(ratiomax);
                for (int iQ2bin=0; iQ2bin<3; iQ2bin++){
                    int ratiosol = i%NvarTS + itarg*NvarTS + (i/NvarTS)*NtargetTS*NvarTS;
                    TH1F *hratio_solid = (TH1F*) input_Vec[ratiosol]->Get(Form("hratio%i",iNubin + iQ2bin*3));
                    int ratioD = i%NvarTS + 0*NvarTS + (i/NvarTS)*NtargetTS*NvarTS;
                    TH1F *hratio_liquid = (TH1F*) input_Vec[ratioD]->Get(Form("hratio%i",iNubin + iQ2bin*3));

                    hratio_solid->Divide(hratio_liquid);
                    // CTratio->Divide(hratio_solid,hratio_liquid);
                    hratio_solid->SetLineColor(color[iQ2bin]);
                    hratio_solid->SetMarkerColor(color[iQ2bin]);
                    hratioStack->Add(hratio_solid);
                    legend->AddEntry(hratio_solid,Form("%.1f < Q^{2} < %.1f GeV^{2}",Q2_limits[iQ2bin],Q2_limits[iQ2bin+1]),"lp");
                }
                hratioStack->Draw("NOSTACK");
                hratioStack->GetYaxis()->CenterTitle();
                hratioStack->GetYaxis()->SetTitleOffset(2.5);
                hratioStack->GetYaxis()->SetNdivisions(505);
                hratioStack->GetXaxis()->CenterTitle();
                hratioStack->GetXaxis()->SetTitleOffset(3.5);
                if (iNubin==0){
                    TText *targText;
                    if (var[i%NvarTS]=="Zh"||var[i%NvarTS]=="Pt2") targText = new TText(0.5,ratiomax+0.05,target[itarg]);
                    else if (var[i%NvarTS]=="PhiPQ") targText = new TText(0,ratiomax+0.05,target[itarg]);
                    targText->SetTextAlign(21);
                    targText->Draw();
                }
                if (itarg==3){
                    TLatex *nuText;
                    if (var[i%NvarTS]=="Zh"||var[i%NvarTS]=="Pt2") nuText = new TLatex(1.08,1.,Form("%.1f < #nu < %.1f [GeV]",
                                                                                       Nu_limits[iNubin],Nu_limits[iNubin+1]));
                    else if (var[i%NvarTS]=="PhiPQ") nuText = new TLatex(210.,1.,Form("%.1f < #nu < %.1f [GeV]",
                                                                         Nu_limits[iNubin],Nu_limits[iNubin+1]));
                    nuText->SetTextAlign(21);
                    nuText->SetTextAngle(90);
                    nuText->Draw();
                }
                if (itarg==3 && iNubin==0){
                    TText *titText;
                    if (var[i%NvarTS]=="Zh"||var[i%NvarTS]=="Pt2") titText = new TText(1.,ratiomax+0.05,method[i/NvarTS]+"_"+var[i%NvarTS]);
                    else if (var[i%NvarTS]=="PhiPQ") titText = new TText(180.,ratiomax+0.05,method[i/NvarTS]+"_"+var[i%NvarTS]);
                    titText->SetTextAlign(21);
                    titText->Draw();

                    legend->Draw();
                }
            }
        }
        cratio->cd();
        // TString par;
        // if (i==0) par = "(";
        // else if (i==3) par = ")";
        // else par = "";
        // cratio->Print("PlotRatio3.pdf"+par,"pdf");
        c_Vec.push_back(cratio);
    }

    // TString pdftitle = "PlotRatioOfCTs0020";
    // for (int ifin=0; ifin<Ncanvas; ifin++){
    //     TString par = "";
    //     if (ifin==0) par = "(";
    //     else if (ifin==Ncanvas-1) par = ")";
    //     c_Vec[(ifin/2)+(ifin%2)*NvarTS]->Print(pdftitle+".pdf"+par,"pdf");
    // }
    TString pdftitle = "PlotDoubleCTRatio3";
    for (int ifin=0; ifin<Ncanvas/2; ifin++){
        TString par = "";
        if (ifin==0) par = "(";
        else if (ifin==Ncanvas/2-1) par = ")";
        c_Vec[ifin]->Print(pdftitle+".pdf"+par,"pdf");
    }

}

//// FUNCTIONS

void CanvasPartition(TCanvas *C,const Int_t Nx,const Int_t Ny,
						Float_t lMargin, Float_t rMargin,
						Float_t bMargin, Float_t tMargin){
	if (!C) return;

	// Setup Pad layout:
	Float_t vSpacing = 0.0;
	Float_t vStep  = (1.- bMargin - tMargin - (Ny-1) * vSpacing) / Ny;

	Float_t hSpacing = 0.0;
	Float_t hStep  = (1.- lMargin - rMargin - (Nx-1) * hSpacing) / Nx;

	Float_t vposd,vposu,vmard,vmaru,vfactor;
	Float_t hposl,hposr,hmarl,hmarr,hfactor;

	for (Int_t i=0;i<Nx;i++) {
		if (i==0) {
			hposl = 0.0;
			hposr = lMargin + hStep;
			hfactor = hposr-hposl;
			hmarl = lMargin / hfactor;
			hmarr = 0.0;
		} else if (i == Nx-1) {
			hposl = hposr + hSpacing;
			hposr = hposl + hStep + rMargin;
			hfactor = hposr-hposl;
			hmarl = 0.0;
			hmarr = rMargin / (hposr-hposl);
		} else {
			hposl = hposr + hSpacing;
			hposr = hposl + hStep;
			hfactor = hposr-hposl;
			hmarl = 0.0;
			hmarr = 0.0;
		}

		for (Int_t j=0;j<Ny;j++) {
			if (j==0) {
				vposd = 0.0;
				vposu = bMargin + vStep;
				vfactor = vposu-vposd;
				vmard = bMargin / vfactor;
				vmaru = 0.0;
			} else if (j == Ny-1) {
				vposd = vposu + vSpacing;
				vposu = vposd + vStep + tMargin;
				vfactor = vposu-vposd;
				vmard = 0.0;
				vmaru = tMargin / (vposu-vposd);
			} else {
				vposd = vposu + vSpacing;
				vposu = vposd + vStep;
				vfactor = vposu-vposd;
				vmard = 0.0;
				vmaru = 0.0;
			}

			C->cd(0);

			char name[16];
			sprintf(name,"pad_%i_%i",i,j);
			TPad *pad = (TPad*) gROOT->FindObject(name);
			if (pad) delete pad;
			pad = new TPad(name,"",hposl,vposd,hposr,vposu);
			pad->SetLeftMargin(hmarl);
			pad->SetRightMargin(hmarr);
			pad->SetBottomMargin(vmard);
			pad->SetTopMargin(vmaru);

			pad->SetFrameBorderMode(0);
			pad->SetBorderMode(0);
			pad->SetBorderSize(0);

			pad->Draw();
		}
	}
}