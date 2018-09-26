#!/bin/bash

if [ $# -lt 1 ] ; then
    echo "  "
    echo "  ./hMerging.sh <rootHistoFolder>"
   # echo "  For more information type on your shell:  <hadd --help>"
    echo "  "
    exit -1
fi

FOLDER="$1"

pushd $FOLDER

#mc
#rm -rf Histo_*_MC.root

#hadd -f -k Histo_avgTrackMultiplicity_MC.root        Histo_avgTrackMultiplicity_pthat*_MC.root 
#hadd -f -k Histo_dRmin_matching_MC.root              Histo_dRmin_matching_pthat*_MC.root 
#hadd -f -k Histo_flight3Dsignif_MC.root              Histo_flight3Dsignif_pthat*_MC.root 
#hadd -f -k Histo_jet_CSV_MC.root                     Histo_jet_CSV_pthat*_MC.root 
#hadd -f -k Histo_jet_JBP_MC.root                     Histo_jet_JBP_pthat*_MC.root 
#hadd -f -k Histo_jet_JP_MC.root                      Histo_jet_JP_pthat*_MC.root 
#hadd -f -k Histo_jet_TCHE_MC.root                    Histo_jet_TCHE_pthat*_MC.root 
#hadd -f -k Histo_jet_TCHP_MC.root                    Histo_jet_TCHP_pthat*_MC.root 
#hadd -f -k Histo_jet_eta_MC.root                     Histo_jet_eta_pthat*_MC.root 
#hadd -f -k Histo_jet_phi_MC.root                     Histo_jet_phi_pthat*_MC.root 
#hadd -f -k Histo_jet_pt_MC.root                      Histo_jet_pt_pthat*_MC.root 
#hadd -f -k Histo_massSV_MC.root                      Histo_massSV_pthat*_MC.root 
#hadd -f -k Histo_nrSV_MC.root                        Histo_nrSV_pthat*_MC.root 
#hadd -f -k Histo_seltrack_IP3D_MC.root               Histo_seltrack_IP3D_pthat*_MC.root 
#hadd -f -k Histo_seltrack_IP3Dsignif_MC.root         Histo_seltrack_IP3Dsignif_pthat*_MC.root 
#hadd -f -k Histo_tracks_Pt_MC.root                   Histo_tracks_Pt_pthat*_MC.root 
#hadd -f -k Histo_tracks_distanceToJetAxis_MC.root    Histo_tracks_distanceToJetAxis_pthat*_MC.root 
#hadd -f -k Histo_tracks_nrPixelHits_MC.root          Histo_tracks_nrPixelHits_pthat*_MC.root 

#data
rm -rf Histo_*_TotalData.root

hadd -f -k Histo_avgTrackMultiplicity_TotalData.root         Histo_avgTrackMultiplicity_Data_*.root 
hadd -f -k Histo_dRmin_matching_TotalData.root               Histo_dRmin_matching_Data_*.root 
hadd -f -k Histo_flight3Dsignif_TotalData.root               Histo_flight3Dsignif_Data_*.root 
hadd -f -k Histo_jet_CSV_TotalData.root                      Histo_jet_CSV_Data_*.root 
hadd -f -k Histo_jet_JBP_TotalData.root                      Histo_jet_JBP_Data_*.root 
hadd -f -k Histo_jet_JP_TotalData.root                       Histo_jet_JP_Data_*.root 
hadd -f -k Histo_jet_TCHE_TotalData.root                     Histo_jet_TCHE_Data_*.root 
hadd -f -k Histo_jet_TCHP_TotalData.root                     Histo_jet_TCHP_Data_*.root 
hadd -f -k Histo_jet_eta_TotalData.root                      Histo_jet_eta_Data_*.root 
hadd -f -k Histo_jet_phi_TotalData.root                      Histo_jet_phi_Data_*.root 
hadd -f -k Histo_jet_pt_TotalData.root                       Histo_jet_pt_Data_*.root 
hadd -f -k Histo_massSV_TotalData.root                       Histo_massSV_Data_*.root 
hadd -f -k Histo_nrSV_TotalData.root                         Histo_nrSV_Data_*.root 
hadd -f -k Histo_seltrack_IP3D_TotalData.root                Histo_seltrack_IP3D_Data_*.root
hadd -f -k Histo_seltrack_IP3Dsignif_TotalData.root          Histo_seltrack_IP3Dsignif_Data_*.root
hadd -f -k Histo_tracks_Pt_TotalData.root                    Histo_tracks_Pt_Data_*.root
hadd -f -k Histo_tracks_distanceToJetAxis_TotalData.root     Histo_tracks_distanceToJetAxis_Data_*.root
hadd -f -k Histo_tracks_nrPixelHits_TotalData.root           Histo_tracks_nrPixelHits_Data_*.root
 
popd
