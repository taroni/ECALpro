{
   gStyle->SetCanvasBorderMode(0);
   gStyle->SetPadBorderMode(0);
   gStyle->SetPadColor(0);
   gStyle->SetCanvasColor(0);
   gStyle->SetTitleColor(1);
   gStyle->SetStatColor(0);

   gStyle->SetLabelSize(0.03,"xyz"); // size of axis values

   // default canvas positioning
   gStyle->SetCanvasDefX(900);
   gStyle->SetCanvasDefY(20);
   gStyle->SetCanvasDefH(550);
   gStyle->SetCanvasDefW(540);

   gStyle->SetPadBottomMargin(0.15); // It was 0.1
   gStyle->SetPadTopMargin(0.15);
   gStyle->SetPadLeftMargin(0.15);
   gStyle->SetPadRightMargin(0.15);

   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);

   gStyle->SetFrameBorderMode(0);

   // US letter
   gStyle->SetPaperSize(20, 24);
   gStyle->SetOptStat(0);

   bool foundIt=true;
   TString FWLiteLib = "libFWCoreFWLite.so";

   const char *cmsbase=gSystem->Getenv("CMSSW_BASE");
   if (cmsbase==NULL) {
     cout << " CMSSW environment has not been setup -- "
   	  << " FWLite libraries will not be loaded\n" << endl;
     foundIt=false;
   } else {
     cout << " CMSSW environment has been setup \n" << endl;

     const char *search=gSystem->Getenv("LD_LIBRARY_PATH");
     string cms_path = search;
     
     const char* foundlib =gSystem->Which(search, FWLiteLib,(EAccessMode)0);
     
     if (! foundlib) {
       FWLiteLib = "libPhysicsToolsFWLite.so";
       foundlib =gSystem->Which(search, FWLiteLib, (EAccessMode)0);
       if (! foundlib) {
   	 cout << "Could not find any FWLite libraries to load " << endl;       
   	 foundIt=false;
       }
     }
   }
   if (foundIt){
     //cout << "Loading: " << FWLiteLib << endl;
     gSystem->Load(FWLiteLib);
     //AutoLibraryLoader::enable();
   }
    const Int_t NRGBs = 5;
    const Int_t NCont = 255;

    Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
    Double_t red[NRGBs]   = { 0.00, 0.00, 0.87, 1.00, 0.51 };
    Double_t green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
    Double_t blue[NRGBs]  = { 0.51, 1.00, 0.12, 0.00, 0.00 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
    gStyle->SetNumberContours(NCont);



}
