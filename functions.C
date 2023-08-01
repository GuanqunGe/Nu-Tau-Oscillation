#include <iostream>
#include <algorithm>

double get_entry(TTree * tree, int i, std::string field) {
    double v;
    tree->SetBranchAddress(field.c_str(), &v);
    tree->GetEvent(i);
    tree->ResetBranchAddresses();
    return v;
}





/*std::vector<double> pt_50(TTree * tree, std::string branch, std::vector<double> vect){
    
    std::vector<double> electron_pt;
    std::vector<double> proton_pt;
    std::vector<double> pT;
    
    for(int iEntry = 0; iEntry < vect.size(); iEntry++){
    
        tree->GetEntry(iEntry);
        
        double pt = 0; double px = 0; double py = 0; double electronpx = 0; double electronpy = 0; double p_pt = 0; double e_pt = 0;
        
        for(int k = 0; k < mctruth_daughters_pdg->size(); k++){
        
            if(mctruth_daughters_pdg->at(k) == 2212 && mctruth_daughters_status_code->at(k) == 1){
            
                px += mctruth_daughters_px->at(k);
                py += mctruth_daughters_py->at(k);
            }
            
            if(abs(mctruth_daughters_pdg->at(k)) == 11&& mctruth_daughters_status_code->at(k) == 1){
            
                electronpx = mctruth_daughters_px->at(k);
                electronpy = mctruth_daughters_py->at(k);
            }
            
            p_pt = sqrt(px*px + py*py);
            e_pt = sqrt(electronpx*electronpx + electronpy*electronpy);
            pt = sqrt((px+electronpx)*(px+electronpx) + (py+electronpy)*(py+electronpy));
            proton_pt.push_back(p_pt);
            electron_pt.push_back(e_pt);
            pT.push_back(pt);
        }
    }
    
     
} */

double mean(TTree * tree, std::string branch){

    double mean;
    for(int iEntry = 0; iEntry<branch.size(); iEntry++){
    
        mean+= get_entry(tree, iEntry, branch);
    }
    
    mean /= branch.size();
    return mean;
}

double std_dev(TTree * tree, std::string branch){

    double dev;
    double avg = mean(tree, branch);
    
    for(int iEntry = 0; iEntry < branch.size(); iEntry++){
    
        dev += (get_entry(tree, iEntry, branch)-avg)*(get_entry(tree, iEntry, branch)-avg);
    }
    
    dev /=branch.size();
    dev = sqrt(dev);
    return dev;
}

TGraph * scatter (TTree * tree1, TTree * tree2, std::string branch1, std::string branch2){

    int nentries = tree1->GetEntries();
    double s_var1[nentries];
    double s_var2[nentries];
    double b_var1[nentries];
    double b_var2[nentries];
    
    
    for(int iEntry = 0; iEntry < nentries; iEntry++){
        
        s_var1[iEntry] = get_entry(tree1, iEntry, branch1);
        s_var2[iEntry] = get_entry(tree1, iEntry, branch2);
        b_var1[iEntry] = get_entry(tree2, iEntry, branch1);
        b_var2[iEntry] = get_entry(tree2, iEntry, branch2);
    }
    TGraph * sig = new TGraph(nentries, s_var1, s_var2);
    TGraph * bkg = new TGraph(nentries, b_var1, b_var2);
    sig->SetTitle("scatter; branch1; branch2 ");
    sig->SetMarkerColor(kBlue);
    bkg->SetMarkerColor(kRed);
    sig->SetMarkerStyle(7);
    bkg->SetMarkerStyle(7);
    sig->Draw("ap");
    bkg->Draw("p");
    auto legend = new TLegend();
    legend->AddEntry(sig, "#nu_{#tau} Signal");
    legend->AddEntry(bkg, "#nu_{e} Background");
    legend->Draw();
    //scatter->Draw("ap");
    return sig;
}

std::vector<double> get_parameters (TH1D * hist, char * func){

    TF1 * g = (TF1*)hist->GetListOfFunctions()->FindObject(func);
    std::vector<double> params;
    params.push_back(g->GetParameter(0));
    params.push_back(g->GetParameter(1));
    
    return params;
} 

std::vector<double> get_significance (TH1D * background, TH1D * signal){
  
    //int numbins = background->GetXaxis()->GetNbins();
    int numbins = signal->FindLastBinAbove(0,1);
    std::vector<double> S;
    //double bins[numbins];
    double width = background->GetBinWidth(1);
    
    for(int i = 0; i<numbins; i++){
        double sig = signal->GetBinContent(i);
        double bkg = background->GetBinContent(i);
        double stat = sig/sqrt(bkg);
        
        if(stat > 1000000000){
        
            S.push_back(0);}
        else{
        
            S.push_back(stat);
        }
        //bins[i] = i * width;
    }
    
    return S;
}

std::vector<double> better_significance (TH1D * background, TH1D * signal){
  
    //int numbins = background->GetXaxis()->GetNbins();
    int numbins = signal->FindLastBinAbove(0,1);
    std::vector<double> S;
    //double bins[numbins];
    double width = background->GetBinWidth(1);
    
    for(int i = 0; i<numbins; i++){
        double sig = signal->GetBinContent(i);
        double bkg = background->GetBinContent(i);
        double nO = sig+bkg;
        double stat = sqrt(2*nO*log(1+sig/bkg)-2*sig);
        
        if(stat > 1000000000){
        
            S.push_back(0);}
        else{
        
            S.push_back(stat);
        }
        //bins[i] = i * width;
    }
    
    return S;
}

std::vector<double> get_bins (TH1D * hist){

    //int num = hist->GetXaxis()->GetNbins();
    int num = hist->FindLastBinAbove(0,1);
    std::vector<double> bins;
    double width = hist->GetBinWidth(1);
    
    for(int i = 0; i <num; i++){
    
        bins.push_back(i * width);
        //std::cout << i*width << " ";
    }
    
    return bins;
}

/* int get_cut_events(TTree * tree){

    int nentries = tree->GetEntries();
    
    std::vector<std::string> names;
    for (size_t i=0; i<tree->GetListOfBranches()->GetEntries(); i++) {
        names.push_back(tree->GetListOfBranches()->At(i)->GetName());
    }
    
    //order goes pt, philm, philh, phihm, protonpt, electronpt, electronE, events
    std::vector<std::vector<double> > table(names.size(), std::vector<double>(nentries));
    for(int i = 0; i<names.size(); i++){
        for(int j = 0; j< nentries; j++){
            double v = get_entry(tree, j, names[i]);
            table[i][j] = v;
            //std::cout << table[i * j] << " ";
        }
        //std::cout<<std::endl;
    }

    int cut_events = 0;
    int pt_cut = 0;
    for(size_t i = 0; i < nentries; i++){
        
        if(table[1][i] > 30 && table[3][i] > 40){
        
            //pt_cut = 0.45 + (180 - table[2][i]) * 0.000833;
            pt_cut = 180 - (0.9333 * table[i];
        
            if(table[6][i] < 2 && table[6][i] > 0.5 && table[2][i] < 140){
            
                cut_events++; 
            }
            
            if(table[6][i] >= 2 && table[6][i] < 2.75 && table[2][i] < 170){
            
                cut_events++;
            }
            
            if(table[6][i] >= 2.75){
            
                cut_events++;
            }
        }
    }
    
    return cut_events;
} */

int philm_phihm_cut(TTree * tree){

double nentries = tree->GetEntries();
double events = 0;

    std::vector<std::string> names;
    for (size_t i=0; i<tree->GetListOfBranches()->GetEntries(); i++) {
        names.push_back(tree->GetListOfBranches()->At(i)->GetName());
    }
    
    //order goes pt, philm, philh, phihm, protonpt, electronpt, electronE, events
    std::vector<std::vector<double> > table(names.size(), std::vector<double>(nentries));
    for(int i = 0; i<names.size(); i++){
        for(int j = 0; j< nentries; j++){
            double v = get_entry(tree, j, names[i]);
            table[i][j] = v;
            //std::cout << table[i * j] << " ";
        }
        //std::cout<<std::endl;
    }
    
    for(int k = 0; k < nentries; k++){
    
        if(table[1][k] > 35){
        
            if((205-0.98*table[1][k]) < table[3][k]){
            
                events++;
            }
        }
    }
    
    double ratio = events/nentries;
    std::cout << "total events: " << nentries << " cut events: " << events << " ratio: " << ratio << std::endl;
    
    int cuts = events;
    return cuts;

}

int pt_philh_cut(TTree * tree){

double nentries = tree->GetEntries();
double events = 0;

    std::vector<std::string> names;
    for (size_t i=0; i<tree->GetListOfBranches()->GetEntries(); i++) {
        names.push_back(tree->GetListOfBranches()->At(i)->GetName());
    }
    
    //order goes pt, philm, philh, phihm, protonpt, electronpt, electronE, events
    std::vector<std::vector<double> > table(names.size(), std::vector<double>(nentries));
    for(int i = 0; i<names.size(); i++){
        for(int j = 0; j< nentries; j++){
            double v = get_entry(tree, j, names[i]);
            table[i][j] = v;
            //std::cout << table[i * j] << " ";
        }
        //std::cout<<std::endl;
    }
    
    for(int k = 0; k< nentries; k++){
    
        if(table[0][k] > 0.5 && table[2][k] < 120){
            events++;
        }
    }

    double ratio = events/nentries;
    std::cout << "total events: " << nentries << " cut events: " << events << " ratio: " << ratio << std::endl;
    
    int cuts = events;
    return cuts;
    
}

TGraph * optimal_cut(TH1D * background, TH1D * signal){
   
   int bins = background->GetXaxis()->GetNbins();
   double width = background->GetBinWidth(1);
   std::cout << "bins: " << bins << std::endl;
   double bkg = 0;
   double sig = 0;
   double S[bins-2];
   double x[bins-2];
   double max = 0;
   double maxbin = 0;
   
   for(int i = 1; i < bins-1; i++){
   
       bkg = background->Integral(i,bins-1);
       if(bkg == 0){bkg = 1;}
       sig = signal->Integral(i, bins-1);
       S[i] = sig/sqrt(bkg);
       x[i] = i * width;
       if(S[i] > max){
           max = S[i]; 
           maxbin = x[i];
       }
   }
   
   std::cout<< "maxbin: " << maxbin << " max: " << max << std::endl; 
   
   TGraph * optimal_cut = new TGraph(bins-2, x, S);
   optimal_cut->SetTitle("Cuts;Corresponding Bins;Significance");
   optimal_cut->Draw("ac");
   return optimal_cut;
}

TGraph * backwards_cut(TH1D * background, TH1D * signal){
   
   int bins = background->GetXaxis()->GetNbins();
   double width = background->GetBinWidth(1);
   std::cout << "bins: " << bins << std::endl;
   double bkg = 0;
   double sig = 0;
   double S[bins-2];
   double x[bins-2];
   double max = 0;
   int maxbin = 0;
   
   for(int i = bins-1; i > 0; i--){
   
       bkg = background->Integral(1, i);
       if(bkg == 0){bkg = 1;}
       sig = signal->Integral(1, i);
       S[i] = sig/sqrt(bkg);
       x[i] = i * width;
       if(S[i] > max){max = S[i]; maxbin = x[i];}
   }
   
   std::cout<< "maxbin: " << maxbin << " max: " << max << std::endl; 
   
   TGraph * backwards_cut = new TGraph(bins-2, x, S);
   backwards_cut->GetYaxis()->SetRangeUser(0.,14.);
   backwards_cut->SetTitle("Reverse Cuts;Corresponding Bins; Significance");
   backwards_cut->Draw("ac");
   return backwards_cut;
}

double twovariable_cut(TH1D * sig_var1,TH1D * bkg_var1, TH1D * sig_var2, TH1D * bkg_var2){

int bins = sig_var1->GetXaxis()->GetNbins();
double width1 = sig_var1->GetBinWidth(1);
double width2 = sig_var2->GetBinWidth(1);
double var1_s = 0;
double var1_b = 0;
double var2_s = 0;
double var2_b = 0;
double max_1 = 0;
double max_2 = 0;
double var1 = 0;
double var2 = 0;
double x1 = 0;
double x2 = 0;
double maxbin_1 = 0;
double maxbin_2 = 0;

//std::vector< std::vector<double> > S;

for(int i = 1; i < bins-1; i++){

    var1_b = bkg_var1->Integral(i,bins);
    if(var1_b == 0){var1_b = 1;}
    var1_s = sig_var1->Integral(i, bins);
    var1 = var1_s/sqrt(var1_b);
    x1 = i * width1;
    if(var1 > max_1){
    
      max_1 = var1; 
      maxbin_1 = x1;
    }

    for(int j = bins-1; j > 0; j--){
    
        var2_b = bkg_var2->Integral(0, j);
        if(var2_b == 0){var2_b = 1;}
        var2_s = sig_var2->Integral(0, j);
        var2 = var2_s/sqrt(var2_b);
        x2 = j * width2;
        if(var2 > max_2){
            
            max_2 = var2;
            maxbin_2 = x2;
        }
        
        //std::cout << "pt bin: " << x1 << " phi_lh bin: " << x2 << " total significance: " << var1 + var2 << std::endl;
    }
}

std::cout << "max pt: " << maxbin_1 << " max phi_lm: " << maxbin_2 << " max significance: " << max_1 + max_2 <<std::endl;
return max_1 + max_2;

}

/*std::vector<float> get_correlation_matrix(std::string signal, std::string background) {
    int nentries = signal->GetEntries();
    
    TFile * f1 = new TFile(signal);
    TTree * sig = (TTree*)f1->Get("tT50");
    TH1D * t_pt_1p = (TH1D*)tau->Get("pt1p");
    TH1D * t_philm_1p = (TH1D*)tau->Get("philm_1p");
    TH1D * t_philh_1p = (TH1D*)tau->Get("philh_1p");
    TH1D * t_phihm_1p = (TH1D*)tau->Get("phihm_1p");
    TH1D * t_protonpt_1p = (TH1D*)tau->Get("protonpt_1p");
    TH1D * t_electronpt_1p = (TH1D*)tau->Get("electronpt_1p");
    TH1D * t_electronE_1p = (TH1D*)tau->Get("electronE_1p");
    
    TFile * f2 = new TFile(background);
    TTree * bkg = (TTree*)f2->Get("eT50");
    TH1D * e_pt_1p = (TH1D*)file->Get("pt1p");
    TH1D * e_philm_1p = (TH1D*)file->Get("philm_1p");
    TH1D * e_philh_1p = (TH1D*)file->Get("philh_1p");
    TH1D * e_phihm_1p = (TH1D*)file->Get("phihm_1p");
    TH1D * e_protonpt_1p = (TH1D*)file->Get("protonpt_1p");
    TH1D * e_electronpt_1p = (TH1D*)file->Get("electronpt_1p");
    TH1D * e_electronE_1p = (TH1D*)file->Get("electronE_1p");
  

    // get list of branch names
    std::vector<std::string> names;
    for (size_t i=0; i<tree->GetListOfBranches()->GetEntries(); i++) {
        names.push_back(tree->GetListOfBranches()->At(i)->GetName());
    }
    
    std::vector<double> minmax(names.size() * 2);
    for(int k = 0; k<names.size(); k++){
    
        if(names[k] == "phi_lm" || names[k] == "phi_lh" || names[k] == "phi_hm"){
        
            minmax[k * 0] = 0;
            minmax[k * 1] = 180;
        }
        
        if(names[k] == "proton_pt" || names[k] == "electron_pt" || names[k] == "pt"){
        
            minmax[k * 0] = 0;
            minmax[k * 1] = 2;
        }
        
        if(names[k] == "electron_E"){
        
            minmax[k * 0] = 0;
            minmax[k * 1] = 5;
        }
    }
    // convert the ntuple to a vector
    std::vector<double> table(names.size() * nentries);
    std::vector<double> binvalues(names.size() * nentries);
    for (size_t i=0; i<nentries; i++) {
        for (size_t j=0; j<names.size(); j++) {
            double v = get_entry(tree, i, names[j]);
            table[j + i * names.size()] = v;
        }
    }
    
    //find which bin each event would fall into
    for(int i = 0; i<names.size(); i++){
    
        TH1D * hist = new TH1D("hist", "hist", 200, minmax[i * 0], minmax[i * 1]);
        for(int j = 0; j<nentries; j++){
        
            double v = get_entry(tree, j, names[i]);
            table[i * j] = v;
            hist->Fill(v);
            binvalues[i * j] = hist->FindLastBinEntry(0,1);
            hist->Reset("ICES");
        }
    }
    
    std::vector<double> pt = better_significance(e_pt_1p, t_pt_1p);
    std::vector<double> philm = better_significance(e_philm_1p, t_philm_1p);
    std::vector<double> philh = better_significance(e_philh_1p, t_philh_1p);
    std::vector<double> phihm = better_significance(e_phihm_1p, t_phihm_1p);
    std::vector<double> protonpt = better_significance(e_protonpt_1p, t_prptonpt_1p);
    std::vector<double> electronpt = better_significance(e_electronpt_1p, t_electronpt_1p);
    std::vector<double> electronE = better_significance(e_electronE_1p, t_electronE_1p);
    std::vector<double> xbins02 = get_bins(t_pt_1p);
    std::vector<double> xbins05 = get_bins(t_electron_1p);
    std::vector<double> xbins0180 = get_bins(t_philm_1p);

    double* S_pt = &pt[1];
    double* S_philm = &philm[1];.q

    double* S_philh = &philh[1];
    double* S_phihm = &phihm[1];
    double* S_protonpt = &protonpt[1];
    double* S_electronpt = &electronpt[1];
    double* S_electronE = &electronE[1];
    double* bins02 = &xbins02[1];
    double* bins05 = &xbins05[1];
    double* bins0180 = &xbins0180[1];
    
    std::vector<double> matrix(names.size() * nentries)
    //find the significance for the given event and variable
    for(int i = 0; i< names.size(); i++){
    
        if(names[i] == "phi_lm" || names[i] == "phi_lh" || names[i] == "phi_hm"){
        
            for(int j = 0; j<nentries; j++){
        
                for(int k = 0; k < xbins0180.size(); k++){
            
                    if(binvalues[i * j] == bins0180[k]){
                    
                        if(i == 1){matrix[i * j] == S_philm[k]}
                        if(i == 2){matrix[i * j] == S_philh[k]}
                        if(i == 3){matrix[i * j] == S_phihm[k]}
                    
                    }
                }
            }
        }
        
        if(names[i] == "proton_pt" || names[i] == "electron_pt" || names[i] == "pt"){
        
            for(int j = 0; j<nentries; j++){
        
                for(int k = 0; k < xbins02.size(); k++){
            
                    if(binvalues[i * j] == bins02[k]){
                    
                        if(i == 0){matrix[i * j] == S_pt[k]}
                        if(i == 4){matrix[i * j] == S_protonpt[k]}
                        if(i == 5){matrix[i * j] == S_electronpt[k]}
                    }
                 }   
            }
        }
        
        if(names[i] == "electron_E"){
        
            for(int j = 0; j<nentries; j++){
        
                for(int k = 0; k < xbins05.size(); k++){
            
                    if(binvalues[i * j] == bins05[k]){
                      
                        marix[i * j] == S_electronE[k]
                    }
                }
            } 
        }        
    }
    
    
    
}  */



void functions(){

    TFile * f = new TFile("tree.root");
    TTree * t = (TTree*)f->Get("eT50");
    TTree * e1p = (TTree*)f->Get("e1p");
    TTree * e40 = (TTree*)f->Get("e1p40");
    TTree * e50f = (TTree*)f->Get("eT50_f");
    TTree * e1p_noc = (TTree*)f->Get("e1p_noc");
    TH1D * e_pt_1p = (TH1D*)f->Get("pt1p");
    TH1D * e_philm_1p = (TH1D*)f->Get("philm_1p");
    TH1D * e_philh_1p = (TH1D*)f->Get("philh_1p");
    TH1D * e_phihm_1p = (TH1D*)f->Get("phihm_1p");
    TH1D * e_protonpt_1p = (TH1D*)f->Get("protonpt_1p");
    TH1D * e_electronpt_1p = (TH1D*)f->Get("electronpt_1p");
    TH1D * e_electronE_1p = (TH1D*)f->Get("electronE_1p");
    TH1D * e_pt_1p40 = (TH1D*)f->Get("pt1p40");
    TH1D * e_philm_1p40 = (TH1D*)f->Get("philm_1p40");
    TH1D * e_philh_1p40 = (TH1D*)f->Get("philh_1p40");
    TH1D * e_phihm_1p40 = (TH1D*)f->Get("phihm_1p40");
    TH1D * e_protonpt_1p40 = (TH1D*)f->Get("protonpt_1p40");
    TH1D * e_electronpt_1p40 = (TH1D*)f->Get("electronpt_1p40");
    TH1D * e_electronE_1p40 = (TH1D*)f->Get("electronE_1p40");
    TH1D * e_pt_50 = (TH1D*)f->Get("pt_50");
    TH1D * e_philm_50 = (TH1D*)f->Get("philm_50");
    TH1D * e_philh_50 = (TH1D*)f->Get("philh_50");
    TH1D * e_phihm_50 = (TH1D*)f->Get("phihm_50");
    TH1D * e_protonpt_50 = (TH1D*)f->Get("protonpt_50");
    TH1D * e_electronpt_50 = (TH1D*)f->Get("electronpt_50");
    TH1D * e_electronE_50 = (TH1D*)f->Get("electronE_50");
    TH1D * e_pt_50f = (TH1D*)f->Get("pt_50f");
    TH1D * e_philm_50f = (TH1D*)f->Get("philm_50f");
    TH1D * e_philh_50f = (TH1D*)f->Get("philh_50f");
    TH1D * e_phihm_50f = (TH1D*)f->Get("phihm_50f");
    TH1D * e_protonpt_50f = (TH1D*)f->Get("protonpt_50f");
    TH1D * e_electronpt_50f = (TH1D*)f->Get("electronpt_50f");
    TH1D * e_electronE_50f = (TH1D*)f->Get("electronE_50f");
    TH1D * e_pt_1pnoc = (TH1D*)f->Get("pt_1pnoc");

    
    TFile * f2 = new TFile("tautree.root");
    TTree * t2 = (TTree*)f2->Get("tT50");
    TTree * t1p = (TTree*)f2->Get("t1p");
    TTree * t40 = (TTree*)f2->Get("t1p40");
    TTree * t50f = (TTree*)f2->Get("tT50_f");
    TTree * t1p_noc = (TTree*)f2->Get("t1p_noc");
    TH1D * t_pt_1p = (TH1D*)f2->Get("pt1p");
    TH1D * t_philm_1p = (TH1D*)f2->Get("philm_1p");
    TH1D * t_philh_1p = (TH1D*)f2->Get("philh_1p");
    TH1D * t_phihm_1p = (TH1D*)f2->Get("phihm_1p");
    TH1D * t_protonpt_1p = (TH1D*)f2->Get("protonpt_1p");
    TH1D * t_electronpt_1p = (TH1D*)f2->Get("electronpt_1p");
    TH1D * t_electronE_1p = (TH1D*)f2->Get("electronE_1p");
    TH1D * t_pt_1p40 = (TH1D*)f2->Get("pt1p40");
    TH1D * t_philm_1p40 = (TH1D*)f2->Get("philm_1p40");
    TH1D * t_philh_1p40 = (TH1D*)f2->Get("philh_1p40");
    TH1D * t_phihm_1p40 = (TH1D*)f2->Get("phihm_1p40");
    TH1D * t_protonpt_1p40 = (TH1D*)f2->Get("protonpt_1p40");
    TH1D * t_electronpt_1p40 = (TH1D*)f2->Get("electronpt_1p40");
    TH1D * t_electronE_1p40 = (TH1D*)f2->Get("electronE_1p40");
    TH1D * t_pt_50 = (TH1D*)f2->Get("pt_50");
    TH1D * t_philm_50 = (TH1D*)f2->Get("philm_50");
    TH1D * t_philh_50 = (TH1D*)f2->Get("philh_50");
    TH1D * t_phihm_50 = (TH1D*)f2->Get("phihm_50");
    TH1D * t_protonpt_50 = (TH1D*)f2->Get("protonpt_50");
    TH1D * t_electronpt_50 = (TH1D*)f2->Get("electronpt_50");
    TH1D * t_electronE_50 = (TH1D*)f2->Get("electronE_50");
    TH1D * t_pt_50f = (TH1D*)f2->Get("pt_50f");
    TH1D * t_philm_50f = (TH1D*)f2->Get("philm_50f");
    TH1D * t_philh_50f = (TH1D*)f2->Get("philh_50f");
    TH1D * t_phihm_50f = (TH1D*)f2->Get("phihm_50f");
    TH1D * t_protonpt_50f = (TH1D*)f2->Get("protonpt_50f");
    TH1D * t_electronpt_50f = (TH1D*)f2->Get("electronpt_50f");
    TH1D * t_electronE_50f = (TH1D*)f2->Get("electronE_50f");
    TH1D * t_pt_1pnoc = (TH1D*)f2->Get("pt_1pnoc");
    
    TFile * taudata = new TFile("/nevis/riverside/data/guanqun/NuTauOscillation/vertexed_MConly_Fullosc_NuTau_CC_20k.root");
    TTree * t_subrun = (TTree*)taudata->Get("run_subrun_tree");
    
    TFile * elecdata = new TFile("/nevis/riverside/data/guanqun/NuTauOscillation/vertexed_MConly_Intrinsic_NuE_CC_20k.root");
    TTree * e_subrun = (TTree*)elecdata->Get("run_subrun_tree");
    
    int t_cut = pt_philh_cut(t50f);
    int e_cut = pt_philh_cut(e50f);
    
    std::cout << "tau: " << t_cut/126 << " elec: " << e_cut/21 << " sig: " << (t_cut/126)/sqrt(e_cut/21) << std::endl;
    
    std::cout << "tau: " << t50f->GetEntries()/126 << " elec: " << e50f->GetEntries()/21 <<std::endl;
    
    std::cout << "tau: " << t2->GetEntries()/126 << " elec: " << t->GetEntries()/21 <<std::endl;
    
    
    TCanvas * a = new TCanvas("a", "pt vs. p_lm");
    scatter(t50f, e50f, "pt", "phi_lm");
    
    /*TPad * a1 = new TPad("a1", "a1", 0,0.5,1,1); a1->Draw();
    a1->cd();
    e_pt_50f->Scale(0.046875);
    t_pt_50f->Scale(0.007947);
    e_pt_50f->Draw("HIST");
    t_pt_50f->Draw("HIST, SAME");
    e_pt_50f->SetLineColor(kRed);
    t_pt_50f->SetLineColor(kBlue);
    e_pt_50f->SetFillStyle(3002);
    e_pt_50f->SetFillColor(kRed);
    t_pt_50f->SetFillColor(kBlue);
    t_pt_50f->SetFillStyle(3018);
    e_pt_50f->SetTitle("50 MeV Filtered Events: p^{T}");
    e_pt_50f->GetXaxis()->SetTitle("p^{T} (GeV)");
    e_pt_50f->GetYaxis()->SetTitle("POT Normalized # of Events");
    e_pt_50f->SetAxisRange(0,2 ,"X");
    
    auto l = new TLegend();
    l->AddEntry(t_pt_50f, "#nu_{#tau} Signal");
    l->AddEntry(e_pt_50f, "#nu_{e} Background");
    l->SetBorderSize(0);
    l->Draw(); 

    
    a->cd();
    TPad * a2 = new TPad("a2", "a2", 0,0,1,0.5); a2->Draw();
    a2->cd();
    //TGraph * sig1 = new TGraph(last_pt-1, x02, S1);
    //sig1->Draw("AC");
    //optimal_cut(e_pt_1p, t_pt_1p);
    //optimal_cut(e_pt_1p40, t_pt_1p40);
    //optimal_cut(e_pt_50, t_pt_50);
    optimal_cut(e_pt_50f, t_pt_50f); */
    
    TCanvas * b = new TCanvas("b", "pt vs. p_lh");
    scatter(t50f, e50f, "pt", "phi_lh");
 
    /*TPad * b1 = new TPad("b1", "b1", 0,0.5,1,1); b1->Draw();
    b1->cd();
    /*e_pt_1pnoc->Scale(0.046875);
    t_pt_1pnoc->Scale(0.007947);
    e_pt_1pnoc->Draw("HIST");
    t_pt_1pnoc->Draw("HIST, SAME");
    t_pt_1pnoc->SetLineColor(kRed);
    e_pt_1pnoc->SetTitle("One Proton Events Filtered: p^{T}");
    e_pt_1pnoc->GetXaxis()->SetTitle("p^{T} (GeV)");
    e_pt_1pnoc->GetYaxis()->SetTitle("Number of Events"); 
    
    e_protonpt_50f->Scale(0.046875);
    t_protonpt_50f->Scale(0.007947);
    e_protonpt_50f->Draw("HIST");
    t_protonpt_50f->Draw("SAME");
    t_protonpt_50f->SetLineColor(kRed);
    e_protonpt_50f->SetTitle("One Proton Events: Proton p^{T}");
    e_protonpt_50f->GetXaxis()->SetTitle("p^{T} (GeV)");
    e_protonpt_50f->GetYaxis()->SetTitle("Number of Events"); 
    
    b->cd();
    TPad * b2 = new TPad("b2", "b2",  0,0,1,0.5); b2->Draw();
    b2->cd();
    //TGraph * sig2 = new TGraph(last_ppt-1, x02, S2);
    //sig2->Draw("AC");
    //optimal_cut(e_protonpt_1p, t_protonpt_1p);
    //optimal_cut(e_protonpt_1p40, t_protonpt_1p40);
    optimal_cut(e_protonpt_50f, t_protonpt_50f); 
    //optimal_cut(e_pt_1pnoc, t_pt_1pnoc); */
    
    TCanvas * c = new TCanvas("c", "pt vs. p_hm");
    scatter(t50f, e50f, "pt", "phi_hm");

    /*TPad * c1 = new TPad("c1", "c1",0,0.5,1,1); c1->Draw();
    c1->cd(); 
    /*e_pt_1p->Scale(0.046875);
    t_pt_1p->Scale(0.007947);
    e_pt_1p->Draw("HIST");
    t_pt_1p->Draw("HIST, SAME");
    t_pt_1p->SetLineColor(kRed);
    e_pt_1p->SetTitle("One Proton Events: p^{T}");
    e_pt_1p->GetXaxis()->SetTitle("p^{T} (GeV)");
    e_pt_1p->GetYaxis()->SetTitle("Number of Events"); 
    
    
    e_electronpt_50f->Scale(0.046875);
    t_electronpt_50f->Scale(0.007947);
    e_electronpt_50f->Draw("HIST");
    t_electronpt_50f->Draw("SAME");
    t_electronpt_50f->SetLineColor(kRed);
    e_electronpt_50f->SetTitle("One Proton Events: Electron p^{T}");
    e_electronpt_50f->GetXaxis()->SetTitle("p^{T} (GeV)");
    e_electronpt_50f->GetYaxis()->SetTitle("Number of Events"); 
    
    c->cd();
    TPad * c2 = new TPad("c2", "c2", 0,0,1,0.5); c2->Draw();
    c2->cd();
    //TGraph * sig3 = new TGraph(last_ept-1, x02, S3);
    //sig3->Draw("AC");
    //optimal_cut(e_electronpt_1p, t_electronpt_1p);
    //optimal_cut(e_electronpt_1p40, t_electronpt_1p40);
    optimal_cut(e_electronpt_50f, t_electronpt_50f);
    //optimal_cut(e_pt_1p, t_pt_1p); */
    
    TCanvas * d = new TCanvas("d", "ept vs. p_lh");
    scatter(t50f, e50f, "electron_pt", "phi_lh");

    /*TPad * d1 = new TPad("d1", "d1", 0,0.5,1,1); d1->Draw();
    d1->cd();
    e_pt_50->Scale(0.046875);
    t_pt_50->Scale(0.007947);
    e_pt_50->Draw("HIST");
    t_pt_50->Draw("HIST, SAME");
    t_pt_50->SetLineColor(kRed);
    e_pt_50->SetTitle("50 MeV Filtered Events: p^{T}");
    e_pt_50->GetXaxis()->SetTitle("p^{T} (GeV)");
    e_pt_50->GetYaxis()->SetTitle("Number of Events");
    
    d->cd();
    TPad * d2 = new TPad("d2", "d2", 0,0,1,0.5); d2->Draw();
    d2->cd();
    //TGraph * sig4 = new TGraph(last_eE-1, x05, S4);
    //sig4->Draw("AC");
    //optimal_cut(e_electronE_1p, t_electronE_1p);
    //optimal_cut(e_electronE_1p40, t_electronE_1p40);
    //optimal_cut(e_electronE_50, t_electronE_50); 
    optimal_cut(e_pt_50, t_pt_50); */
    
    TCanvas * e = new TCanvas("e", "eE vs. p_lh");
    scatter(t50f, e50f, "electron_E", "phi_lh");
  
    /*TPad * e1 = new TPad("e1", "e1", 0,0.5,1,1); e1->Draw();
    e1->cd();
    e_philm_50f->Scale(0.046875);
    t_philm_50f->Scale(0.007947);
    e_philm_50f->Draw("HIST");
    t_philm_50f->Draw("SAME");
    t_philm_50f->SetLineColor(kRed);
    e_philm_50f->SetTitle("One Proton Events: #phi_{lm}");
    e_philm_50f->GetXaxis()->SetTitle("Angle (Degrees)");
    e_philm_50f->GetYaxis()->SetTitle("Number of Events");
    //e_philm_50f->SetAxisRange(0, 800, "Y");
    
    e->cd();
    TPad * e2 = new TPad("e2", "e2", 0,0,1,0.5); e2->Draw();
    e2->cd();
    //TGraph * sig5 = new TGraph(last_tlm-1, x0180, S5);
    //sig5->Draw("AC");
    //backwards_cut(e_philm_1p, t_philm_1p);
    //backwards_cut(e_philm_1p40, t_philm_1p40);
    //backwards_cut(e_philm_50, t_philm_50); 
    backwards_cut(e_philm_50f, t_philm_50f); */
    
    TCanvas * g = new TCanvas("g", "p_lm vs. p_hm");
    scatter(t50f, e50f, "phi_lm", "phi_hm");
  
    /*TPad * g1 = new TPad("g1", "g1", 0,0.5,1,1); g1->Draw();
    g1->cd();
    e_philh_50f->Scale(0.046875);
    t_philh_50f->Scale(0.007947);
    e_philh_50f->Draw("HIST");
    t_philh_50f->Draw("HIST, SAME");
    t_philh_50f->SetLineColor(kRed);
    e_philh_50f->GetXaxis()->SetTitle("Angle (Degrees)");
    e_philh_50f->GetYaxis()->SetTitle("Number of Events");
    e_philh_50f->SetTitle("One Proton Events: #phi_{lh}");
    
    g->cd();
    TPad * g2 = new TPad("g2", "g2",  0,0,1,0.5); g2->Draw();
    g2->cd();
    //TGraph * sig6 = new TGraph(last_tlh-1, x0180, S6);
    //sig6->Draw("AC");
    //backwards_cut(e_philh_1p, t_philh_1p);
    //backwards_cut(e_philh_1p40, t_philh_1p40);
    //backwards_cut(e_philh_50, t_philh_50);
    backwards_cut(e_philh_50f, t_philh_50f); */
    
    TCanvas * h = new TCanvas("h", "h");
  
    TPad * h1 = new TPad("h1", "h1",0,0.5,1,1); h1->Draw();
    h1->cd();
    e_phihm_50f->Scale(0.046875);
    t_phihm_50f->Scale(0.007947);
    e_phihm_50f->Draw("HIST");
    t_phihm_50f->Draw("SAME");
    t_phihm_50f->SetLineColor(kRed);
    e_phihm_50f->GetXaxis()->SetTitle("Angle (Degrees)");
    e_phihm_50f->GetYaxis()->SetTitle("Number of Events");
    //e_phihm_50f->SetAxisRange(0, 550, "Y");
    e_phihm_50f->SetTitle("One Proton Events: #phi_{hm}");
    
    h->cd();
    TPad * h2 = new TPad("h2", "h2",  0,0,1,0.5); h2->Draw();
    h2->cd();
    //TGraph * sig7 = new TGraph(last_thm-1, x0180, S7);
    //sig7->Draw("AC");
    //backwards_cut(e_phihm_1p, t_phihm_1p);
    //backwards_cut(e_phihm_1p40, t_phihm_1p40);
    //backwards_cut(e_phihm_50, t_phihm_50); 
    //backwards_cut(e_phihm_1pnoc, t_phihm_1pnoc);
    backwards_cut(e_phihm_50f, t_phihm_50f); 
    
    TCanvas * j = new TCanvas("j", "j");
    
    //scatter(t1p, e1p, "pt1p", "proton_pt");
    //scatter(t1p, "pt1p", "phi_lh")->Draw("p");
    //scatter(t1p, "pt1p", "phi_lh")->SetLineColor(kRed); */
   
    
   /* int t_cut_events40 = get_cut_events(t40);
    int e_cut_events40 = get_cut_events(e40);
    
    int t_cut_events = get_cut_events(t1p);
    int e_cut_events = get_cut_events(e1p);
    
    //std::cout<< "tau left: " << t_cut_events.size() << "  e left: " << e_cut_events.size() * 5.9 << std::endl;
    
std::cout << "tau events: " << t_cut_events/126 << "  e events: " << e_cut_events/21  << "  sig: " << (t_cut_events/126)/sqrt(e_cut_events/21) << std::endl;
std::cout << "tau events 40: " << t_cut_events40/126 << "  e events: " << e_cut_events40/21 << "  sig: " << (t_cut_events40/126)/sqrt(e_cut_events40/21) << std::endl;   

double give = twovariable_cut(t_pt_1p, e_pt_1p, t_philm_1p, e_philm_1p);
std::cout << give <<std::endl;

double giving = twovariable_cut(t_electronE_1p, t_electronE_1p, t_philh_1p, e_philh_1p);
std::cout << "eE vs. philh: " << giving <<std::endl;

std::cout << "just pt cut on 0.57: " << t_pt_1p->Integral(57,200)/sqrt(e_pt_1p->Integral(57,200));

std::cout << "eE cut on 0.025: " << "tau events - " << t_electronE_1p->Integral() << "  elec events: " << e_electronE_1p->Integral() << " significance - " << t_electronE_1p->Integral()/sqrt(e_electronE_1p->Integral()) <<std::endl;

std::cout << "philh cut on 154: " << "tau events - " << t_philh_1p->Integral()/126 << "  elec events: " << e_philh_1p->Integral(0, 154)/21 << " significance - " << (t_philh_1p->Integral()/126)/sqrt(e_philh_1p->Integral(0, 154)/21) << std::endl;

std::cout << "pt 50f cut on 0.65 GeV: " << " tau - " << t_pt_50f->Integral(65,200)/126 << "  elec - " << e_pt_50f->Integral(65,200)/21 << " significance - " << (t_pt_50f->Integral()/126)/sqrt(e_pt_50f->Integral(0, 154)/21) << std::endl;
*/

}

