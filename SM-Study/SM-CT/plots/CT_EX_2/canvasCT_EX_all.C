// Four targets, closure tests over all electron variable's range
void CanvasPartition(TCanvas *C,const Int_t Nx = 2,const Int_t Ny = 2,
                     Float_t lMargin = 0.15, Float_t rMargin = 0.05,
                     Float_t bMargin = 0.15, Float_t tMargin = 0.05);

void canvasCT_EX_all(){
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
    Color_t color[4] = {kRed, kMagenta, kBlue, kGreen};

    gStyle->SetErrorX(0);
    gStyle->SetTitleStyle(0);
    gStyle->SetOptTitle(1);
    gStyle->SetMarkerStyle(7);

    // Fill canvases
    std::vector<TCanvas*> c_Vec;

    int Ncanvas = NmethodTS*NvarTS;
    for (int i=0; i<Ncanvas; i++){
        TCanvas *cratio = new TCanvas("cratio"+method[i/NvarTS]+"_"+var[i%NvarTS],"CT ratio "+method[i/NvarTS]+var[i%NvarTS]);

        for (int itarg=0; itarg<NtargetTS; itarg++){
            cratio->cd(0);

            TString htitle = "CT ratio all ("+var[i%NvarTS]+"), "+method[i/NvarTS]+";"+varxtitle[i%NvarTS]+";Corr/True";

            int ratiobin = i%NvarTS + itarg*NvarTS + (i/NvarTS)*NtargetTS*NvarTS;
            TH1F *hratio_tmp = (TH1F*) input_Vec[ratiobin]->Get("hratio_all");
            hratio_tmp->SetTitle(htitle);

            double ratiomin = 0.0001; // 0.5001;
            double ratiomax = 1.9999; // 1.4999;
            hratio_tmp->SetMinimum(ratiomin);
            hratio_tmp->SetMaximum(ratiomax);

            hratio_tmp->GetYaxis()->CenterTitle();
            hratio_tmp->GetXaxis()->CenterTitle();

            hratio_tmp->SetLineColor(color[itarg]);
            hratio_tmp->SetMarkerColor(color[itarg]);
            hratio_tmp->Draw("SAME");
        }
        cratio->cd();
        c_Vec.push_back(cratio);
    }

    TString pdftitle = "PlotRatio_all0020";
    for (int ifin=0; ifin<Ncanvas; ifin++){
        TString par = "";
        if (ifin==0) par = "(";
        else if (ifin==Ncanvas-1) par = ")";
        c_Vec[(ifin/2)+(ifin%2)*NvarTS]->Print(pdftitle+".pdf"+par,"pdf");
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