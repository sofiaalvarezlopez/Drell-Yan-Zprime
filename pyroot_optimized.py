import ROOT
from ROOT import TLorentzVector, TH1F

def PT(TLV):
    return TLV.Pt()

# Debemos buscar los dos muones de mayor p_T que satisfagan la condición de corte. 
# mu_list es una lista ordenada descendentemente según los p_T de los muones.
# Con este código, me estoy evitando comparaciones dobles y solo hago O(n).
# Queremos dos muones de carga opuesta. Entonces, la multiplicación de sus cargas debe ser menor a 0.

def cut(mu_list, pt_cut=20, delta_R=0.3):
  i, j = 0, 1
  mu1, mu2 = None, None
  mu_cut = []
  cut = False
  while i < len(mu_list):
    mu1, mu2 = mu_list[i], mu_list[j]
    if mu1.Pt() >= pt_cut and mu2.Pt() >= pt_cut and mu1.DeltaR(mu2) > delta_R and mu1.charge()*mu2.charge() < 0:
      cut = True
      mu_cut.append((mu1, mu2))
    else:
      if j == len(mu_list) - 1 and i < len(mu_list): 
        i += 1
        j = i + 1
      else:
        j += 1
  return mu_cut, cut



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
        # Estos son los dos muones de mayor pT.
        mu1, mu2 = muons[0], muons[1]
        # Tomamos los muones que satisfagan el corte. Cada uno es un TLorentzVector.
        # Si es estrictamente la pareja de mayor pT, paso muons ordenada. Si no, puede estar desordenada y cogemos la primera
        # pareja que satisfaga la condición de corte. 
        muons_cut, cut = cut(muons)
        # Esto me da la pareja de muones de mayor p_T que además cumplen el Delta_R.
        leading_pair = muons_cut[0]
        # Cogemos los dos muones que satisfacen el corte
        mu1_cut, mu2_cut = leading_pair[0], leading_pair[1]

        # Hacemos la gráfica de p_T de los dos muones de mayor p_T
        plot_PT_mu1.Fill(mu1.Pt())
        plot_PT_mu2.Fill(mu2.Pt())
        # Hacemos la grafica de Delta_R de los dos muones con mayor pT y Delta_R > 0.3
        plot_Delta_R_mu1_mu2.Fill(mu1_cut.DeltaR(mu2_cut))
        # Sumamos los dos muones cut para hacer la gráfica de masa reconstruida. 
        mu_total = mu1 + mu2
        # Hacemos la gráfica de masa reconstruida de mu_total.
        plot_mRec_mu1_mu2.Fill(mu_total.M())

        for muon in muons:
            plot_ETA_muons.Fill(muon.Pt())
            plot_PHI_muons.Fill(muon.Pt())

        #MET
       # for MET in METs:
            #plot_MET.Fill(MET.Pt())




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
  
