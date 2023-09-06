import ROOT

def renameTH1(filename, oldName, newName):
    # Open the ROOT file in update mode
    file = ROOT.TFile(filename, "UPDATE")

    # Access the TH1 object by its name
    th1 = file.Get(oldName)

    # Check if the object exists
    if th1:
        # Rename the TH1 object
        th1.SetName(newName)
        th1.Write(newName,ROOT.TObject.kWriteDelete);

        # Update the directory in the file

        # Print a message to indicate success
        #print(f"Renamed {oldName} to {newName} in {filename}")
    #else:
        # Print an error message if the TH1 object doesn't exist
        #print(f"Error: {oldName} not found in {filename}")

    # Close the file
        print "ok"
        print th1.GetName()
    
    file.Save()
    file.Close()

def changeName():
    filename = "crab_Analysis_SingleMuon_Run2017_2018_CodeV73p3_v4_cutIndex3_rebinEta4_rebinIh4_rebinP2_rebinMass1_nPE200_25april_VR_SR_cp.root"
    oldName = "mass_obs_999ias100"
    newName = "data_obs"

    renameTH1(filename, oldName, newName)

changeName()
