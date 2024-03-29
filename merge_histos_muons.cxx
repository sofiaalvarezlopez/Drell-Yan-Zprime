#include <iostream>
#include <vector>
#include <string>
#include <TH1F.h>
#include <fstream>
using namespace std;

void merge_histos_muons()
{
  
    //vector<std::string> names_files {"bkg1.root", "bkg2.root", "bkg3.root", "bkg4.root", "signal1.root", "signal2.root", "signal3.root", "signal4.root", "signal5.root", "signal6.root", "signal7.root", "signal8.root"};
    vector<std::string> names_files {"Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/muon_root_files/DY+jets.root", "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/muon_root_files/ttbar.root", "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/muon_root_files/ww.root", "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/muon_root_files/wz.root", "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/muon_root_files/zz.root"};

    vector<std::string> names_files_wo_ext {"DY+jets", "t\\bar{t}", "WW", "WZ", "ZZ"};

    vector<std::string> names {"DY+jets", "t#bar{t}", "WW", "WZ", "ZZ"};

    vector<std::string> plots {"ETA_muons", "PHI_muons", "PT_mu1", "PT_mu2", "DELTA_R muons", "M_REC muons", "Charge muons", "DELTA_PHI muons", "Cos_DELTA_PHI muons", "MET", "Cos_DELTA_PHI MET v. Muon lead", "Transverse_mass"};

    vector<std::string> plots_names {"#eta(#mu)", "#phi(#mu)", "p_{T}(#mu_{1})", "p_{T}(#mu_{2})", "#Delta R (#mu_{1}, #mu_{2})", "M(#mu_{1}, #mu_{2})", "Q(#mu_{1})Q(#mu_{2})", "#Delta #phi (#mu_{1}, #mu_{2})", "cos(#Delta #phi (#mu_{1}, #mu_{2}))", "MET", "cos(#Delta #phi (MET, #mu_{1}))", "M_{T}(#mu_{1}#mu_{2})"};

    vector<std::string> x_labels {"#eta(#mu) [a.u.]", "#phi(#mu) [a.u.]", "p_{T}(#mu_{1}) [GeV]", "p_{T}(#mu_{2}) [GeV]", "#Delta R (#mu_{1}, #mu_{2}) [a.u.]", "M(#mu_{1}, #mu_{2}) [GeV]", "Q(#mu_{1})Q(#mu_{2}) [e]", "#Delta #phi (#mu_{1}, #mu_{2}) [rad]", "cos(#Delta #phi (#mu_{1}, #mu_{2})) [a.u.]", "MET [GeV]", "cos(#Delta #phi (MET, #mu_{1})) [a.u.]", "M_{T}(#mu_{1}#mu_{2}) [GeV]"};


    //vector<int> colors {3, 7, 6, 5, 2, 4, 8, 9, 1, 43, 97, 38};
    vector<int> colors {6, 7, 2, 4, 3};

    //vector<int> linestyles {1, 1, 1, 1, 10, 9, 8, 7, 6, 5, 4, 3};
    vector<int> linestyles {1, 1, 1, 1, 1};

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
          if( (names_files[j] == "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/muon_root_files/DY+jets.root") || (names_files[j] == "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/muon_root_files/ttbar.root") || (names_files[j] == "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/muon_root_files/ww.root") || (names_files[j] == "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/muon_root_files/wz.root") || (names_files[j] == "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/muon_root_files/zz.root") || (names_files[j] == "signal6.root") || (names_files[j] == "signal7.root") || (names_files[j] == "signal8.root") )
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
    
        if( (plots[i]=="PT_mu1") || (plots[i]=="PT_mu2"))
        {
          legend->AddEntry(h,names[j].c_str(),"l");   
        }

        else
        {
          if( (names_files[j] == "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/muon_root_files/DY+jets.root") || (names_files[j] == "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/muon_root_files/ttbar.root") || (names_files[j] == "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/muon_root_files/ww.root") || (names_files[j] == "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/muon_root_files/wz.root") || (names_files[j] == "Documents/Universidad/Noveno Semestre/Proyecto Teórico:Computacional/pyroot/muon_root_files/zz.root") || (names_files[j] == "signal6.root") || (names_files[j] == "signal7.root") || (names_files[j] == "signal8.root") )
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
        hs->GetYaxis()->SetTitle("Num. Eventos (normalizados a la unidad) [a.u.]");
        legend->Draw();
        l->Add(c2);
        std::string filename = plots[i] + ".png";
        c2->SaveAs(filename.c_str());
    }

    TFile* Output = new TFile("joined_muons.root", "RECREATE"); 
    l->Write();
    Output->Close();


       

    
}
