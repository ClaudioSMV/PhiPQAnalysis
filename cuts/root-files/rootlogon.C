{
  // rootlogon.C
  // Manuel Calderon de la Barca {// Add my own options here:
  TStyle *mcStyle = new TStyle("mcStyle", "Manuel's Root Styles");
  mcStyle->SetCanvasDefW(800);  // 800
  mcStyle->SetCanvasDefH(500);  // 600
  
  // use plain black on white colors
  Int_t icol=0; // WHITE
  mcStyle->SetFrameBorderMode(icol);
  mcStyle->SetFrameFillColor(icol);
  mcStyle->SetCanvasBorderMode(icol);
  mcStyle->SetCanvasColor(icol);
  mcStyle->SetPadBorderMode(icol);
  mcStyle->SetPadColor(icol);
  mcStyle->SetStatColor(icol);
  //mcStyle->SetFillColor(icol); // don't use: white fill color for *all* objects

  // set the paper & margin sizes
  mcStyle->SetPaperSize(20,26);

  // set margin sizes
  mcStyle->SetPadTopMargin(0.07);     // 0.05
  mcStyle->SetPadRightMargin(0.18);   // 0.05
  mcStyle->SetPadBottomMargin(0.1);  // 0.16
  mcStyle->SetPadLeftMargin(0.18);    // 0.16

  // set title offsets (for axis label)
  mcStyle->SetTitleXOffset(1);
  mcStyle->SetTitleYOffset(1.4);

  // use large fonts
  // Int_t font=72; // Helvetica italics
  Int_t font=42; // Helvetica
  Double_t tsize=0.04;
  Double_t Tsize=0.05;
  // Int_t font=43; // Helvetica
  // Double_t tsize=18;
  mcStyle->SetTextFont(font);

  mcStyle->SetTextSize(Tsize);
  mcStyle->SetLabelFont(font,"x");
  mcStyle->SetTitleFont(font,"x");
  mcStyle->SetLabelFont(font,"y");
  mcStyle->SetTitleFont(font,"y");
  mcStyle->SetLabelFont(font,"z");
  mcStyle->SetTitleFont(font,"z");
  
  mcStyle->SetLabelSize(tsize,"x");
  mcStyle->SetTitleSize(tsize,"x");
  mcStyle->SetLabelSize(tsize,"y");
  mcStyle->SetTitleSize(tsize,"y");
  mcStyle->SetLabelSize(tsize,"z");
  mcStyle->SetTitleSize(tsize,"z");

  // mcStyle->SetTitleStyle(0);

  // use bold lines and markers
  mcStyle->SetMarkerStyle(9); // 9
  mcStyle->SetMarkerSize(1);
  // mcStyle->SetHistLineWidth(2.);
  mcStyle->SetHistLineWidth(1.);
  mcStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars 
  //mcStyle->SetErrorX(0.001);
  // get rid of error bar caps
  mcStyle->SetEndErrorSize(0);

  // do not display any of the standard histogram decorations (default)
  mcStyle->SetOptTitle(); // 0 (default= without title)
  mcStyle->SetTitleX(0.1);
  mcStyle->SetTitleY(1);
  mcStyle->SetTitleW(0.8);
  // mcStyle->SetTitleH(0.5);
  //mcStyle->SetOptStat(1111);
  mcStyle->SetOptStat(0); // 0 (default)
  //mcStyle->SetOptFit(1111);
  mcStyle->SetOptFit(0); // 0 (default)

  // put tick marks on top and RHS of plots
  mcStyle->SetPadTickX(1);
  mcStyle->SetPadTickY(1);
  
  gROOT->SetStyle("mcStyle");
  gROOT->ForceStyle();
  std::cout << "Styles are Set!" << std::endl;

//  gInterpreter->AddIncludePath("-I/Users/lopez/temp/carolina/analysis/eic-smear/install/include/");

  return;
}
