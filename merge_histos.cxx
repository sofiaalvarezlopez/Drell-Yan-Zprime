#include <iostream>
#include <vector>
#include <string>
#include <TH1F.h>
#include <fstream>
using namespace std;

void merge_histos()
{
  
    //vector<std::string> names_files {"bkg1.root", "bkg2.root", "bkg3.root", "bkg4.root", "signal1.root", "signal2.root", "signal3.root", "signal4.root", "signal5.root", "signal6.root", "signal7.root", "signal8.root"};
    vector<std::string> names_files {"Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/DY+jets.root", "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/ttbar.root", 
    "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/ww.root", "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/wz.root", "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/zz.root", 
    "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/w+jets.root", "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/Zprime_tata_350.root", "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/Zprime_tata_1000.root", 
    "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/Zprime_tata_3000.root"};


    //vector<std::string> names_files_wo_ext {"$t\\bar{t}h$", "$t\\bar{t}t\\bar{t}$", "$t\\bar{t}b\\bar{b}/h$", "$WWZ$", "$m(Z')=350\\,GeV$", "$m(Z')=500\\,GeV$", "$m(Z')=750\\,GeV$", "$m(Z')=1000\\,GeV$", "$m(Z')=1500\\,GeV$", "$m(Z')=2000\\,GeV$", "$m(Z')=3000\\,GeV$", "$m(Z')=4000\\,GeV$"};
    vector<std::string> names_files_wo_ext {"DY+jets", "t\\bar{t}", "WW", "WZ", "ZZ", "W+jets", "$m(Z')=350\\,GeV$", "$m(Z')=1000\\,GeV$", "$m(Z')=3000\\,GeV$"};

    //vector<std::string> names {"t#bar{t}h", "t#bar{t}t#bar{t}", "t#bar{t}b#bar{b}/h", "WWZ", "t#bar{t}Z', m(Z')= 350 GeV", "t#bar{t}Z', m(Z')= 500 GeV", "t#bar{t}Z', m(Z')= 750 GeV", "t#bar{t}Z', m(Z')= 1000 GeV", "t#bar{t}Z', m(Z')= 1500 GeV", "t#bar{t}Z', m(Z')= 2000 GeV", "t#bar{t}Z', m(Z')= 3000 GeV", "t#bar{t}Z', m(Z')= 4000 GeV"};
    //vector<std::string> names {"DY+jets", "t#bar{t}", "WW", "WZ", "ZZ"};
    vector<std::string> names {"DY+jets", "t#bar{t}", "WW", "WZ", "ZZ", "w+jets", "m(Z')= 350 GeV", "m(Z')= 1000 GeV", "m(Z')= 3000 GeV"};

    //vector<std::string> plots {"ETA_muons", "PHI_muons", "PT_mu1", "PT_mu2", "DELTA_R muons", "M_REC muons", "Charge muons", "DELTA_PHI muons", "Cos_DELTA_PHI muons", "MET", "Cos_DELTA_PHI MET v. Muon lead", "Transverse_mass"};
    vector<std::string> plots {"Cos_DELTA_PHI taus", "Cos_DELTA_PHI MET tau1", "Cos_DELTA_PHI MET tau2", "DELTA_Pt taus", "DELTA_PHI taus", "DELTA_ETA taus", "DELTA_R taus", "PT_tau1", "PT_tau2", "ETA_taus", "PHI_taus", "P_tau1", "P_tau2"};

    //vector<std::string> plots_names {"#eta(#mu)", "#phi(#mu)", "p_{T}(#mu_{1})", "p_{T}(#mu_{2})", "#Delta R (#mu_{1}, #mu_{2})", "M(#mu_{1}, #mu_{2})", "Q(#mu_{1})Q(#mu_{2})", "#Delta #phi (#mu_{1}, #mu_{2})", "cos(#Delta #phi (#mu_{1}, #mu_{2}))", "MET", "cos(#Delta #phi (MET, #mu_{1}))", "M_{T}(#mu_{1}#mu_{2})"};
    vector<std::string> plots_names {"cos(#Delta #phi (#tau_{1}, #tau_{2}))", "cos(#Delta #phi (MET, #tau_{1}))", "cos(#Delta #phi (MET, #tau_{2}))", "#Delta p_{T} (#tau_{1}, #tau_{2})", "#Delta #phi (#tau_{1}, #tau_{2})", "#Delta #eta (#tau_{1}, #tau_{2})", 
    "#Delta R (#tau_{1}, #tau_{2})", "p_{T}(#tau_{1})", "p_{T}(#tau_{2})", "#eta(#tau)", "#phi(#tau)", "p(#tau_{1})", "p(#tau_{2})"};


    vector<std::string> x_labels {"cos(#Delta #phi (#tau_{1}, #tau_{2})) [a.u.]", "cos(#Delta #phi (MET, #tau_{1})) [a.u.]", "cos(#Delta #phi (MET, #tau_{2})) [a.u.]", "#Delta p_{T} (#tau_{1}, #tau_{2}) [GeV]",
    "#Delta #phi (#tau_{1}, #tau_{2}) [a.u.]", "#Delta #eta (#tau_{1}, #tau_{2}) [a.u.]", "#Delta R (#tau_{1}, #tau_{2}) [a.u]", "p_{T}(#tau_{1}) [GeV]", "p_{T}(#tau_{2}) [GeV]", "#eta(#tau) [GeV]", "#phi(#tau) [a.u.]", "p(#tau_{1}) [GeV]", "p(#tau_{2}) [GeV]"};


    //vector<int> colors {3, 7, 6, 5, 2, 4, 8, 9, 1, 43, 97, 38};
    //vector<int> colors {3, 7, 6, 4, 2};
    vector<int> colors {3, 7, 6, 5, 4, 8, 43, 97, 38};

    //vector<int> linestyles {1, 1, 1, 1, 10, 9, 8, 7, 6, 5, 4, 3};
    //vector<int> linestyles {1, 1, 1, 1, 1};
    vector<int> linestyles {1, 1, 1, 1, 1, 1, 10, 9, 8};

    TList *l = new TList();

    for(int i=0; i<plots.size(); i++)
    {
        THStack *hs = new THStack("hs", plots_names[i].c_str());
        TCanvas *c2 = new TCanvas(plots[i].c_str(),"Histos",1280,1024);  
        Double_t x_1,x_2;
        if (plots[i]=="M_dijet_partially" || plots[i]=="M_b_dijet_partially" || plots[i]=="M_b_dijet_fully" || plots[i]=="M_dijet_no" || plots[i]=="M_b_not_used_diff" || plots[i]=="M_b1b2"){
            x_1 = 0.15;
            x_2 = 0.35;
        }

        else if (plots[i]=="M_b_not_used_diff_before" || plots[i]=="M_b1b2_before"){
            x_1 = 0.45;
            x_2 = 0.65;
        }
        else{
            x_1 = 0.65;
            x_2 = 0.85;
        }

        auto legend = new TLegend(x_1,0.65,x_2,0.85);
    for (int j=0; j<names.size(); j++)
    {    
        TFile f(names_files[j].c_str());
        TH1F *h = (TH1F*)f.Get(plots[i].c_str());
        h->SetDirectory(0);

        if( (plots[i]=="N_Merged") || (plots[i]=="Eff_mu") || (plots[i]=="Eff_e") )
        {
          h->SetLineColor(colors[j]);
          h->SetLineStyle(linestyles[j]);
          h->SetLineWidth(2);   
        }

        else
        {
          if( (names_files[j] == "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/Zprime_tata_350.root") || (names_files[j] == "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/Zprime_tata_1000.root") || (names_files[j] == "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/Zprime_tata_3000.root") )
            {  
              h->SetLineColor(colors[j]);
              h->SetLineStyle(linestyles[j]);
              h->SetLineWidth(2);
            }

            else
            {
              h->SetFillColor(colors[j]);
              h->SetFillStyle(1001);
              h->SetLineColor(0);
            }
        }

        h->Scale(1.0/h->Integral());
    
        if( (plots[i]=="N_Merged") || (plots[i]=="Eff_mu") || (plots[i]=="Eff_e") )
        {
          legend->AddEntry(h,names[j].c_str(),"l");   
        }

        else
        {
          if( (names_files[j] == "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/Zprime_tata_350.root") || (names_files[j] == "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/Zprime_tata_1000.root") || (names_files[j] == "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/Zprime_tata_3000.root") || (names_files[j] == "signal6.root") || (names_files[j] == "signal7.root") || (names_files[j] == "signal8.root") )
            {
              legend->AddEntry(h,names[j].c_str(),"l");
            }

            else
            {
              legend->AddEntry(h,names[j].c_str(),"f");
            }
        }
     
        legend->SetBorderSize(0);
        hs->Add(h);
         
    }   

        hs->Draw("NOSTACK HIST");
        hs->GetXaxis()->SetTitle(x_labels[i].c_str());
        hs->GetYaxis()->SetTitle("a.u.");
        legend->Draw();
        l->Add(c2);
        std::string filename = plots[i] + ".png";
        c2->SaveAs(filename.c_str());
    }

    TFile* Output = new TFile("joined_tau.root", "RECREATE"); 
    l->Write();
    Output->Close();


       

    
}
