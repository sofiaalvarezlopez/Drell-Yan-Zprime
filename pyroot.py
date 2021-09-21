import ROOT
import numpy as np
from ROOT import TLorentzVector, TH1F



def PT(TLV):
    return TLV.Pt()

def muon_filter(muon):
    return muon[1].Pt()


# Debemos buscar los muones de mayor p_T que satisfagan la condición de corte.
# mu_list es una lista ordenada descendentemente según los p_T de los muones.
# Retornamos todos los muones que satisfacen las condiciones.
# Note que la primera pareja sera la de mayor p_T y delta_R según el cut.
def cuts(mu_list, pt_cut=20, delta_R=0.3):
    global mu_cut
    global cut
    mu_cut = []
    cut = False
    for i, tuple_1 in enumerate(mu_list):
        for j, tuple_2 in enumerate(mu_list):
            if j > i:
                mu_1, mu_2 = tuple_1[1], tuple_2[1]
                if mu_1.Pt() >= pt_cut and mu_2.Pt() >= pt_cut and mu_1.DeltaR(mu_2) > delta_R:
                    cut = True
                    mu_cut.append((mu_1, mu_2))
    return mu_cut, cut


signals = ["w+jets", "ttbar", "ww", "wz", "zz", "DY+jets"]
ext = ["/root/data_docker/SIM_D1/SIMULACIONES/", "/root/data_docker/SIM_D2/SIMULACIONES/", "/root/data_docker/SIM_D3/SIMULACIONES/",
       "/root/data_docker/SIM_D3/SIMULACIONES/", "/root/data_docker/SIM_D3/SIMULACIONES/", "/root/data_docker/SIM_D3/SIMULACIONES/"]
sufijos = ["m_delphes_events.root", "m_delphes_events.root", "m_delphes_events.root", "m_delphes_events.root", "m_delphes_events.root",
           "tag_1_delphes_events.root"]
jobs = [1,1,1,1,1,1]

c1 = ROOT.TCanvas("c1", "Titulo")
plot_PT_mu1 = TH1F("PT_mu1", "PT_mu1", 100, 0.0, 1000.0)
plot_PT_mu2 = TH1F("PT_mu2", "PT_mu2", 100, 0.0, 1000.0)
plot_ETA_muons = TH1F("ETA_muons", "ETA_muons", 100, -8.0, 8.0)
plot_PHI_muons = TH1F("PHI_muons", "PHI_muons", 100, -8.0, 8.0)
plot_Delta_R_mu1_mu2 = TH1F("DELTA_R muons", "DELTA_R muons", 100, 0, 10.0)
plot_mRec_mu1_mu2 = TH1F("M_REC muons", "M_REC muons", 100, 0, 100.0)

plot_Delta_PHI_muons = TH1F("DELTA_PHI muons", "DELTA_PHI muons", 100, -8.0, 8.0)
plot_Cos_Delta_PHI_muons = TH1F("Cos_DELTA_PHI muons", "Cos_DELTA_PHI muons", 100, -1.0, 1.0)
plot_MET = TH1F("MET", "MET", 100, 0.0, 100.0)
plot_charge_muons = TH1F("Charge muons", "Charge muons", 2, -1.0, 1.0)
plot_cos_Delta_PHI_MET_muon_lead = TH1F("Cos_DELTA_PHI MET v. Muon lead", "Cos_DELTA_PHI MET v. Muon lead", 100, -1.0, 1.0)
plot_transverse_mass = TH1F("Transverse_mass", "Transverse_mass", 100, 0.0, 1000)

for n_signal, signal in enumerate(signals):

    f = ROOT.TFile(signal + ".root", "recreate")

    for ind in range(1, jobs[n_signal] + 1):
        directory = str(
            ext[n_signal] + signal + "/" + signal + "_" + str(ind) + "/Events/run_01/" + sufijos[n_signal])
        File = ROOT.TChain("Delphes;1")
        File.Add(directory)
        Number = File.GetEntries()

        print("Signal: " + signal + "_" + str(ind))

        for i in range(Number):
            # print("Evento: " + str(i))
            Entry = File.GetEntry(i)

            jets = []
            electrons = []
            muons = []
            leptons = []
            METs = []

            EntryFromBranch_j = File.Jet.GetEntries()
            for j in range(EntryFromBranch_j):
                jet = TLorentzVector()
                jet_PT, jet_Eta, jet_Phi, jet_M = File.GetLeaf("Jet.PT").GetValue(j), File.GetLeaf("Jet.Eta").GetValue(
                    j), File.GetLeaf("Jet.Phi").GetValue(j), File.GetLeaf("Jet.Mass").GetValue(j)
                jet.SetPtEtaPhiM(jet_PT, jet_Eta, jet_Phi, jet_M)
                jets.append(jet)

            EntryFromBranch_e = File.Electron.GetEntries()
            for j in range(EntryFromBranch_e):
                electron = TLorentzVector()
                electron_PT, electron_Eta, electron_Phi, electron_M = File.GetLeaf("Electron.PT").GetValue(
                    j), File.GetLeaf("Electron.Eta").GetValue(j), File.GetLeaf("Electron.Phi").GetValue(j), 0.000511
                electron.SetPtEtaPhiM(electron_PT, electron_Eta, electron_Phi, electron_M)
                electrons.append(electron)

            EntryFromBranch_mu = File.Muon.GetEntries()
            for j in range(EntryFromBranch_mu):
                muon = TLorentzVector()
                muon_PT, muon_Eta, muon_Phi, muon_M, muon_charge = File.GetLeaf("Muon.PT").GetValue(j), File.GetLeaf(
                    "Muon.Eta").GetValue(j), File.GetLeaf("Muon.Phi").GetValue(j), 0.1056583755, File.GetLeaf("Muon.Charge").GetValue(j)
                muon.SetPtEtaPhiM(muon_PT, muon_Eta, muon_Phi, muon_M)
                muons.append((muon_charge, muon))

            EntryFromBranch_MET = File.MissingET.GetEntries()
            for j in range(EntryFromBranch_MET):
                MET = TLorentzVector()
                MET_PT, MET_Eta, MET_Phi, MET_M = File.GetLeaf("MissingET.MET").GetValue(j), File.GetLeaf(
                    "MissingET.Eta").GetValue(j), File.GetLeaf("MissingET.Phi").GetValue(j), 0.0
                MET.SetPtEtaPhiM(MET_PT, MET_Eta, MET_Phi, MET_M)
                METs.append(MET)

            leptons = electrons + muons

            if (len(muons) >= 2):

                jets.sort(reverse=True, key=PT)
                muons.sort(reverse=True, key=muon_filter)
                # j1, j2 = jets[0], jets[1]
                # Estos son los dos muones de mayor pT.
                mu1, mu2 = muons[0][1], muons[1][1]
                # Tomamos los muones que satisfagan el corte. Cada uno es un TLorentzVector.
                # Si es estrictamente la pareja de mayor pT, paso muons ordenada. Si no, puede estar desordenada y cogemos la primera
                # pareja que satisfaga la condición
                muons_cut, cut = cuts(muons)
                # Esto me da la pareja de muones de mayor p_T que además cumplen el Delta_R.

                if cut:
                    leading_pair = muons_cut[0]
                    # Cogemos los dos muones que satisfacen el corte
                    mu1_cut, mu2_cut = leading_pair[0], leading_pair[1]

                ############################# HISTOS ###########################
                if cut:
                    # Hacemos la grafica de Delta_R de los dos muones con mayor pT y Delta_R > 0.3
                    plot_Delta_R_mu1_mu2.Fill(mu1_cut.DeltaR(mu2_cut))

                # Hacemos la gráfica de p_T de los dos muones de mayor p_T
                plot_PT_mu1.Fill(mu1.Pt())
                plot_PT_mu2.Fill(mu2.Pt())


                plot_Delta_PHI_muons.Fill(mu1.DeltaPhi(mu2))

                # Sumamos los dos muones para hacer la gráfica de masa reconstruida.
                mu_total = mu1 + mu2
                plot_mRec_mu1_mu2.Fill(mu_total.M())
                # Hacemos la gráfica de masa reconstruida de mu_total.
                for muon in muons:
                    plot_ETA_muons.Fill(muon[1].Eta())
                    plot_PHI_muons.Fill(muon[1].Phi())
                    plot_Cos_Delta_PHI_muons.Fill(np.cos(muon[1].Phi()))

                # Hacemos la gráfica de la carga de los muones
                charge = muons[0][0] * muons[1][0]
                if charge < 0: plot_charge_muons.Fill(-1)
                else: plot_charge_muons.Fill(1)

                # MET
                for MET in METs:
                    plot_MET.Fill(MET.Pt())

                plot_transverse_mass.Fill((METs[0] + mu1).Mt())


                # Cogemos el PHI del MET
                phi_met = METs[0].Phi()
                plot_cos_Delta_PHI_MET_muon_lead.Fill(np.cos(phi_met - mu1.Pt()))






    plot_ETA_muons.Draw("HIST")
    plot_PHI_muons.Draw("HIST")
    plot_PT_mu1.Draw("HIST")
    plot_PT_mu2.Draw("HIST")
    plot_Delta_R_mu1_mu2.Draw("HIST")
    plot_mRec_mu1_mu2.Draw("HIST")

    plot_Delta_PHI_muons.Draw("HIST")
    plot_Cos_Delta_PHI_muons.Draw("HIST")
    plot_MET.Draw("HIST")
    plot_charge_muons.Draw("HIST")
    plot_cos_Delta_PHI_MET_muon_lead.Draw("HIST")
    plot_transverse_mass.Draw("HIST")


    x_axis_muons_charge = plot_charge_muons.GetXaxis()
    x_axis_muons_charge.SetBinLabel(1, 'OS')
    x_axis_muons_charge.SetBinLabel(2, 'SS')
    plot_charge_muons.Draw("HIST")

    plot_ETA_muons.Write()
    plot_PHI_muons.Write()
    plot_PT_mu1.Write()
    plot_PT_mu2.Write()
    plot_Delta_R_mu1_mu2.Write()
    plot_mRec_mu1_mu2.Write()
    plot_charge_muons.Write()
    plot_Delta_PHI_muons.Write()
    plot_Cos_Delta_PHI_muons.Write()
    plot_MET.Write()
    plot_charge_muons.Write()
    plot_cos_Delta_PHI_MET_muon_lead.Write()
    plot_transverse_mass.Write()
    plot_charge_muons.Write()

    f.Close()

    plot_ETA_muons.Reset("ICESM")
    plot_PHI_muons.Reset("ICESM")
    plot_PT_mu1.Reset("ICESM")
    plot_PT_mu2.Reset("ICESM")
    plot_Delta_R_mu1_mu2.Reset("ICESM")
    plot_mRec_mu1_mu2.Reset("ICESM")
    plot_Delta_PHI_muons.Reset("ICESM")
    plot_Cos_Delta_PHI_muons.Reset("ICESM")
    plot_MET.Reset("ICESM")
    plot_charge_muons.Reset("ICESM")
    plot_cos_Delta_PHI_MET_muon_lead.Reset("ICESM")
    plot_transverse_mass.Reset("ICESM")
    plot_charge_muons.Reset("ICESM")