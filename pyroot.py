import ROOT
from ROOT import TLorentzVector, TH1F

def PT(TLV):
    return TLV.Pt()

# Debemos buscar los dos muones de mayor p_T que satisfagan la condición de corte. 
# mu_list es una lista ordenada descendentemente según los p_T de los muones.
# Con este código, me estoy evitando comparaciones dobles y solo hago O(n). 
def cut(mu_list, pt_cut=20, delta_R=0.3):
  flag = True
  i, j = 0, 1
  mu1, mu2 = None, None
  while flag and i < len(mu_list):
    mu1, mu2 = mu_list[i], mu_list[j]
    if mu1.Pt() >= pt_cut and mu2.Pt() >= pt_cut and mu1.DeltaR(mu2) > 0.3:
      flag = False
    else:
      if j == len(mu_list) - 1 and i < len(mu_list): 
        i += 1
        j = i + 1
      else:
        j += 1
  return mu1, mu2



# Esto todavía no lo cambio, depende de la señal que me den.
signals = ["ttbarh", "ttbarbbar_noh"]
jobs = [1, 1]

c1 = ROOT.TCanvas("c1", "Titulo")
plot_PT_leptons = TH1F("PT_leptons", "PT_leptons", 100, 0.0, 1000.0)
plot_PT_mu1 = TH1F("PT_mu1", "PT_mu1", 100, 0.0, 1000.0)
plot_PT_mu2 = TH1F("PT_mu2", "PT_mu2", 100, 0.0, 1000.0)
plot_ETA_muons = TH1F("ETA_muons", "ETA_muons", 100, -8.0, 8.0)
plot_PHI_muons = TH1F("PHI_muons", "PHI_muons", 100, -8.0, 8.0)
plot_Delta_R_mu1_mu2 = TH1F("DELTA_R muons", "DELTA_R muons", 100, 0, 1.0)
plot_mRec_mu1_mu2 = TH1F("M_REC muons", "M_REC muons", 100, 0, 300.0)

for n_signal, signal in enumerate(signals):

  f = ROOT.TFile(signal + ".root", "recreate")

  for ind in range(1,jobs[n_signal]+1):
    directory= str("/disco3/SIMULACIONES/with_delphes/" + signal + "/" + signal + "_" + str(ind) + "/Events/run_01/tag_1_delphes_events.root")
    File = ROOT.TChain("Delphes;1")
    File.Add(directory)
    Number = File.GetEntries()

    print("Signal: " + signal + "_" + str(ind))

    for i in range(Number):
    	#print("Evento: " + str(i))
        Entry = File.GetEntry(i)

        jets = []
        electrons = []
        muons = []
        leptons = []
        METs = []

        EntryFromBranch_j = File.Jet.GetEntries()
        for j in range(EntryFromBranch_j):
		        #print("PT = {}, BTag = {}".format(PT,BTag) )
          jet = TLorentzVector()
          jet_PT, jet_Eta, jet_Phi, jet_M  = File.GetLeaf("Jet.PT").GetValue(j), File.GetLeaf("Jet.Eta").GetValue(j), File.GetLeaf("Jet.Phi").GetValue(j), File.GetLeaf("Jet.Mass").GetValue(j)
          jet.SetPtEtaPhiM(jet_PT, jet_Eta, jet_Phi, jet_M)
          jets.append(jet)

        EntryFromBranch_e = File.Electron.GetEntries()
        for j in range(EntryFromBranch_e):
          electron = TLorentzVector()
          electron_PT, electron_Eta, electron_Phi, electron_M  = File.GetLeaf("Electron.PT").GetValue(j), File.GetLeaf("Electron.Eta").GetValue(j), File.GetLeaf("Electron.Phi").GetValue(j), 520.998
          electron.SetPtEtaPhiM(electron_PT, electron_Eta, electron_Phi, electron_M)
          electrons.append(electron)

        EntryFromBranch_mu = File.Muon.GetEntries()
        for j in range(EntryFromBranch_mu):
          muon = TLorentzVector()
          muon_PT, muon_Eta, muon_Phi, muon_M  = File.GetLeaf("Muon.PT").GetValue(j), File.GetLeaf("Muon.Eta").GetValue(j), File.GetLeaf("Muon.Phi").GetValue(j), 105658.374
          muon.SetPtEtaPhiM(muon_PT, muon_Eta, muon_Phi, muon_M)
          muons.append(muon)

        EntryFromBranch_MET = File.MissingET.GetEntries()
        for j in range(EntryFromBranch_MET):
          MET = TLorentzVector()
          MET_PT, MET_Eta, MET_Phi, MET_M  = File.GetLeaf("MissingET.MET").GetValue(j), File.GetLeaf("MissingET.Eta").GetValue(j), File.GetLeaf("MissingET.Phi").GetValue(j), 0.0
          MET.SetPtEtaPhiM(MET_PT, MET_Eta, MET_Phi, MET_M)
          METs.append(MET)

        leptons = electrons + muons


        if (len(jets) >= 2 and len(leptons) != 0):

            jets.sort(reverse = True, key=PT)
            muons.sort(reverse = True, key=PT)
            j1, j2 = jets[0], jets[1]
            # Tomamos los dos muones de mayor p_T que satisfagan el corte. Cada uno es un TLorentzVector
            mu1, mu2 = cut(muons)
            # Hacemos la gráfica de p_T de los dos muones de mayor p_T
            plot_PT_mu1.Fill(mu1.Pt())
            plot_PT_mu2.Fill(mu2.Pt())
            # Hacemos la grafica de Delta_R de los dos muones con mayor pT
            plot_Delta_R_mu1_mu2.Fill(mu1.DeltaR(mu2))
            # Hacemos la gráfica de masa reconstruida de los dos muones de mayor pT
            plot_mRec_mu1_mu2.Fill(mu1.M() + mu2.M())

            for muon in muons:
                plot_ETA_muons.Fill(muon.Pt())
                plot_PHI_muons.Fill(muon.Pt())

            #MET
            for MET in METs:
                plot_MET.Fill(MET.Pt())




  plot_PT_leptons.Draw("HIST")
  plot_ETA_muons.Draw("HIST")
  plot_PHI_muons.Draw("HIST")
  plot_PT_mu1.Draw("HIST")
  plot_PT_mu2.Draw("HIST")
  plot_Delta_R_mu1_mu2.Draw("HIST")
  plot_mRec_mu1_mu2.Draw("HIST")




  plot_PT_leptons.Write()
  plot_ETA_muons.Write()
  plot_PHI_muons.Write()
  plot_PT_mu1.Write()
  plot_PT_mu2.Write()
  plot_Delta_R_mu1_mu2.Write()
  plot_mRec_mu1_mu2.Write()



  f.Close()

  plot_PT_leptons.Reset("ICESM")
  plot_ETA_muons.Reset("ICESM")
  plot_PHI_muons.Reset("ICESM")
  plot_PT_mu1.Reset("ICESM")
  plot_PT_mu2.Reset("ICESM")
  plot_Delta_R_mu1_mu2.Reset("ICESM")
  plot_mRec_mu1_mu2.Reset("ICESM")
  