// Four targets acceptance comparison side by side

void canvasCT(){
    TString target[] = {"D", "C", "Fe", "Pb"};
    TString method[] = {"SM", "CSMV"};

    std::vector<TFile*> input_Vec;

    for (int i=0; i<8; i++){
        TString in_title = "ClosTest"+method[i/4]+"_"+target[i%4]+"1.root";
        TFile *input = TFile::Open(in_title,"READ");
        input_Vec.push_back(input);
    }

    TCanvas *cratio = new TCanvas("cratio","Closure Test ratio");
    TCanvas *cratio_rb = new TCanvas("cratio_rb","Closure Test ratio");
    cratio->Divide(2,1);
    cratio_rb->Divide(2,1);

    std::vector<TLegend*> legend_Vec;

    THStack *hratioSM_stack = new THStack("hratioSM_stack","Closure Test ratio SM;#phi_{PQ} [deg];Corr/True");
    hratioSM_stack->SetMinimum(0.8);
    THStack *hratioCSMV_stack = new THStack("hratioCSMV_stack","Closure Test ratio CSMV;#phi_{PQ} [deg];Corr/True");
    hratioCSMV_stack->SetMinimum(0.8);
    THStack *hratioSM_rb_stack = new THStack("hratioSM_rb_stack","Closure Test ratio SM;#phi_{PQ} [deg];Corr/True");
    hratioSM_rb_stack->SetMinimum(0.95);
    THStack *hratioCSMV_rb_stack = new THStack("hratioCSMV_rb_stack","Closure Test ratio CSMV;#phi_{PQ} [deg];Corr/True");
    hratioCSMV_rb_stack->SetMinimum(0.95);

    Color_t color[4] = {kRed, kMagenta, kBlue, kGreen};

    for (int i=0; i<8; i++){
        if (i%4==0){
            TLegend *leg1 = new TLegend(0.1,0.1,0.4,0.3);
            TLegend *leg2 = new TLegend(0.1,0.1,0.4,0.3);
            leg1->SetBorderSize(0);
            leg2->SetBorderSize(0);
            leg1->SetNColumns(2);
            leg2->SetNColumns(2);
            legend_Vec.push_back(leg1);
            legend_Vec.push_back(leg2);
        }
        TH1F *hist = (TH1F*) input_Vec[i]->Get("hratio");
        hist->SetLineColor(color[i%4]);
        hist->SetMarkerColor(color[i%4]);

        legend_Vec[2*(i/4)]->AddEntry(hist,target[i%4]+"_"+method[i/4],"lp");

        if (i/4 == 0) hratioSM_stack->Add(hist);
        else hratioCSMV_stack->Add(hist);

        TH1F *hist_rb = (TH1F*) input_Vec[i]->Get("hratio_rebin");
        hist_rb->SetLineColor(color[i%4]);
        hist_rb->SetMarkerColor(color[i%4]);

        legend_Vec[2*(i/4)+1]->AddEntry(hist_rb,target[i%4]+"_"+method[i/4],"lp");

        if(i/4 == 0) hratioSM_rb_stack->Add(hist_rb);
        else hratioCSMV_rb_stack->Add(hist_rb);
    }

    cratio->cd(1);
    hratioSM_stack->Draw("NOSTACK");
    legend_Vec[0]->Draw();
    cratio->cd(2);
    hratioCSMV_stack->Draw("NOSTACK");
    legend_Vec[2]->Draw();
    cratio_rb->cd(1);
    hratioSM_rb_stack->Draw("NOSTACK");
    legend_Vec[1]->Draw();
    cratio_rb->cd(2);
    hratioCSMV_rb_stack->Draw("NOSTACK");
    legend_Vec[3]->Draw();

    cratio->Print("PlotRatio.pdf(","pdf");
    cratio_rb->Print("PlotRatio.pdf)","pdf");

}
