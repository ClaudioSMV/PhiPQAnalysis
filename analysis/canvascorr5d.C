// Four targets acceptance comparison side by side
void CanvasPartition(TCanvas *C,const Int_t Nx = 2,const Int_t Ny = 2,
                     Float_t lMargin = 0.15, Float_t rMargin = 0.05,
                     Float_t bMargin = 0.15, Float_t tMargin = 0.05);

void canvascorr5d(TString target = "Fe", TString file_n = "*", TString binning_name = "", TString draw_opt = ""){

    if (binning_name==""){
		std::cerr << "[ERROR] write a correct extension for the acceptance/output file" << std::endl;
        return -1;
	}

    TString file_name;
    if (file_n!="*") file_name = "./Corr5d/Corr5d_"+target+file_n+"_"+binning_name+".root";
    else file_name = "./Corr5d/Corr5d_"+target+"A_"+binning_name+".root";
    TFile *input = TFile::Open(file_name,"READ");

    // Open files
    int NQ2 = 3;
    int NNu = 3;
    int Npads = NQ2*NNu;

    // Canvas options and variables
    Color_t color[3] = {kViolet+10, kRed, kBlack};
    float Q2_limits[] = {1.0, 1.3, 1.8, 4.1};
    float Nu_limits[] = {2.2, 3.2, 3.7, 4.2};

    gStyle->SetErrorX(0);
    gStyle->SetTitleStyle(0);
    gStyle->SetOptTitle(0);
    gStyle->SetMarkerStyle(7);

    // Fill canvases
    std::vector<TCanvas*> c_Vec;

    TString str_type[3] = {"Raw", "Corrected", "New corrected"};
    TString str_type_in[3] = {"raw", "corr", "corrnew"};
    TString str_var[3] = {"Zh", "Pt2", "PhiPQ"};
    // TString str_var_title[3] = {"Z_{h}", "Pt2", "PhiPQ"};

    // Number of PADS
    const Int_t Nx = NQ2;
    const Int_t Ny = NNu;

    // Margins
    Float_t lMargin = 0.06;
    Float_t rMargin = 0.12;
    Float_t bMargin = 0.06;
    Float_t tMargin = 0.06;

    Float_t vStep  = (1.- bMargin - tMargin) / Ny;
    Float_t hStep  = (1.- lMargin - rMargin) / Nx;

    for (int iType=0; iType<3; iType++){
        for (int iVar=0; iVar<3; iVar++){
            TCanvas *canv = new TCanvas(Form("canv%i%i",iType,iVar),str_type[iType]+" "+str_var[iVar]);
            CanvasPartition(canv,Nx,Ny,lMargin,rMargin,bMargin,tMargin);
            c_Vec.push_back(canv);
        }
    }

    for (int iType=0; iType<3; iType++){
        // Canvas setup
        TPad *pad[Nx][Ny];

        for (int iQ2=0; iQ2<NQ2; iQ2++){
            for (int iNu=0; iNu<NNu; iNu++){

                TString str_hname = Form("h"+str_type_in[iType]+"%i%i",iQ2,iNu);
                TH3D *h3d = (TH3D*) input->Get(str_hname);

                for (int iVar=0; iVar<3; iVar++){
                    c_Vec[iVar + 3*iType]->cd(0);
                    pad[iQ2][iNu] = (TPad*) gROOT->FindObject(Form("pad_%i_%i",iQ2,2-iNu));
                    pad[iQ2][iNu]->Draw();
                    pad[iQ2][iNu]->cd();

                    TH1D *hProj;
                    if (iVar==0) hProj = (TH1D*) h3d->ProjectionX(Form(str_type_in[iType]+"Zh%i%i",iQ2,iNu));
                    else if (iVar==1) hProj = (TH1D*) h3d->ProjectionY(Form(str_type_in[iType]+"Pt2%i%i",iQ2,iNu));
                    else if (iVar==2) hProj = (TH1D*) h3d->ProjectionZ(Form(str_type_in[iType]+"PhiPQ%i%i",iQ2,iNu));

                    // hProj->SetTitle(str_type[iType]);
                    hProj->GetYaxis()->SetTitle("Counts");
                    hProj->GetYaxis()->CenterTitle();
                    hProj->GetXaxis()->CenterTitle();
                    hProj->Draw(draw_opt);

                    if (iNu==0){
                        float x_crd = 0.5, y_crd = 1.-0.5/(1+hStep/tMargin);
                        if (iQ2==0) x_crd = (1+1/(1+hStep/lMargin))/2.;
                        else if (iQ2==2) x_crd = 0.5/(1 + rMargin/hStep);
                        TLatex *Q2_title = new TLatex(x_crd,y_crd,Form("%.1f < Q^{2} [GeV^{2}] < %.1f",Q2_limits[iQ2],Q2_limits[iQ2+1]));
                        // TString Q2_title = Form("%.1f < Q^{2} < %.1f",Q2_limits[iQ2],Q2_limits[iQ2+1]);
                        Q2_title->SetNDC();
                        Q2_title->SetTextSize(12);
                        Q2_title->SetTextAlign(22);
                        Q2_title->Draw();
                    }
                    if (iQ2==2){
                        float x_crd = 1 - (3./4)*(1/(1+hStep/rMargin)), y_crd = 0.5;
                        if (iNu==0) y_crd = 0.5/(1 + tMargin/vStep);
                        else if (iNu==2) y_crd = (1+1/(1+vStep/bMargin))/2.;
                        TLatex *Nu_title = new TLatex(x_crd,y_crd,Form("%.1f < #nu < %.1f",Nu_limits[iNu],Nu_limits[iNu+1]));
                        // TLatex *nuText = new TLatex(1.1,1.,Form("%.1f < #nu < %.1f [GeV]",Nu_limits[iNu],Nu_limits[iNu+1]));
                        // TString Nu_title = Form("%.1f < #nu < %.1f",Nu_limits[iNu],Nu_limits[iNu+1]);
                        Nu_title->SetNDC();
                        Nu_title->SetTextSize(12);
                        Nu_title->SetTextAlign(22);
                        Nu_title->SetTextAngle(90);
                        Nu_title->Draw();
                    }
                    if (iQ2==2 && iNu==0){
                        TText *titText = new TText(0.98,0.98,str_type[iType]);
                        titText->SetNDC();
                        titText->SetTextSize(12);
                        titText->SetTextAlign(33);
                        titText->Draw();
                    }
                }
                // TString htitle = ";"+varxtitle[i%2]+";Corr/True";
                // hratioStack->SetMinimum(0.5111);
                // hratioStack->SetMaximum(1.4999);
                // TString Q2_title = Form("%.1f < Q^{2} < %.1f",Q2_limits[iQ2],Q2_limits[iQ2+1]);
                // TString Nu_title = Form("%.1f < #nu < %.1f",Nu_limits[iNu],Nu_limits[iNu+1]);
                // hratioStack->GetYaxis()->CenterTitle();
                // hratioStack->GetXaxis()->CenterTitle();
            }
        }
    }

    TString pdftitle;
    if (file_n!="*") pdftitle = "P_"+target+file_n+"_"+binning_name;
    else pdftitle = "P_A"+target+"_"+binning_name;

    for (int ifin=0; ifin<9; ifin++){
        TString par = "";
        if (ifin==0) par = "(";
        else if (ifin==8) par = ")";
        c_Vec[ifin/3 + 3*(ifin%3)]->Print(pdftitle+".pdf"+par,"pdf");
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