#!/usr/bin/env python
# make_nice_pictures.py

# NAME
#    make_nice_pictures.py repacks graphs in the 'output.root' file
#
# SYNOPSIS
#    ./make_nice_pictures.py - creates a new file called
# 'nice_output.root' with a bit more useful ROOT canvas.

import numpy as np
import ROOT


# Read the input file & tree, create the output file & tree
f_in = ROOT.TFile("output.root", "read")
#t_in = ROOT.TTree("t", "t") # tree name, tree title
f_out = ROOT.TFile("nice_output.root", "recreate")
#t_out = ROOT.TTree("t", "t")


# Define a small function to get ROOT value extracting a bit easier.
# It returns a command (a string); if you execute this command,
# you get the current event value for the specified branch.
def GetValue(branch_name):
    return "f_in.t.GetLeaf(\"" + branch_name + "\").GetValue()"

# Create the histograms we shall fill
num_entries = int(f_in.t.GetEntries())
num_bins_2D = 180
num_bins_1D = 120
# Axis range for 2D histograms
x_min = -5
x_max = 5
y_min = -5
y_max = 5
# Axis range for 1D histograms
#z_min = 0.
#z_max = 2.
TH2D_names = ["theta_flat", "theta_f0_980",
              "theta_f0_600", "theta_f0_1370",
              "theta_f0_1500",
              "theta_f2_1270"]
TH2D_re_names = [name[:6] + "re_" + name[6:] for name in TH2D_names]
TH2D_im_names = [name[:6] + "im_" + name[6:] for name in TH2D_names]

TH1D_names = ["theta_background_flat_abs2",
              "theta_background_rho_770_abs2"]


# Define the input coordinates for our parameters
# - we plot them as lines
TH2D_coord = [-1.1789795785524197, 0.67794332606520935, 1.3694066410273278, 0.29107636714486301, 3.6949292785919234, -0.19364303809889216, 1.2136545544463624, -0.46587833440889037, 0.79127378037251639, -0.76412420750489707, -1.1437419735315568, -1.7612081926853904]
TH2D_re  = TH2D_coord[::2] # Separate real and imaginary coordinates
TH2D_im  = TH2D_coord[1::2]
TH1D_coord = [0.2, 0.2]

#x_min = [TH2D_re[i]*0.9 for i in range(len(TH2D_re))] 
#y_min = [TH2D_im[i]*0.9 for i in range(len(TH2D_im))] 
x_min = [-5 for i in range(len(TH2D_re))] 
y_min = [5 for i in range(len(TH2D_im))] 
z_min = [0, 0]

x_max = [-5 for i in range(len(TH2D_re))] 
y_max = [5 for i in range(len(TH2D_im))] 
z_max = [2, 2]

#y_min[1]=-5
#y_max[1]=5
x_min[1]=1.2
x_max[1]=2.0


TH2D_x_lines = [ROOT.TLine(TH2D_re[i], -5, TH2D_re[i], 5) for i in range(len(TH2D_re))]
TH2D_y_lines = [ROOT.TLine(-5, TH2D_im[i], 5, TH2D_im[i]) for i in range(len(TH2D_re))]

#TH2D_x_lines = [ROOT.TLine(TH2D_re[i], y_min[i], TH2D_re[i], y_max[i]) for i in range(len(TH2D_re))]
#TH2D_y_lines = [ROOT.TLine(x_min[i], TH2D_im[i], x_max[i], TH2D_im[i]) for i in range(len(TH2D_im))]
#

TH1D_lines = [ROOT.TLine(x, 0, x, 1000) for x in TH1D_coord]


# Reminder - initialization of the 2D histogram
# TH2D(const char* name, const char* title, 
#      Int_t nbinsx, Double_t xlow, Double_t xup, 
#      Int_t nbinsy, Double_t ylow, Double_t yup)

# Declare ROOT histograms for the tnames specified above
for i in range(len(TH2D_names)):
    # Create a histogram of type TH2D named hist_name
    exec(TH2D_names[i] + " = ROOT.TH2D(\"" + TH2D_names[i] + "\",\"" + 
         TH2D_names[i] + "\"," + str(num_bins_2D) + ",x_min[i],x_max[i],"
         + str(num_bins_2D) + ",y_min[i],y_max[i])\n")

# Declare reference histogram
#ref_hist = ROOT.TH2D("rho_770", "rho_770 (fixed parameter)", 1, -5, 5, 1, -5, 5)
#x_min, x_max, 1, y_min, y_max)


for i in range(len(TH1D_names)):
    # Create a branch of type TH1D named hist_name
    exec(TH1D_names[i] + " = ROOT.TH1D(\"" + TH1D_names[i] + "\",\"" + 
         TH1D_names[i] + "\"," + str(num_bins_1D) + ",z_min[i],z_max[i])\n")


# Loop over the number of samples
for i in range(num_entries):
    f_in.t.GetEntry(i) # Load the i-th event

    # Loop over histograms - fill each of them
    for k in range(len(TH2D_names)):
        exec(TH2D_names[k] + ".Fill(" + 
             GetValue(TH2D_re_names[k]) + "," +
             GetValue(TH2D_im_names[k]) + ")\n")

    for hist_k in TH1D_names:
        exec(hist_k + ".Fill(" + GetValue(hist_k) + ")\n")


# Create a canvas, plot the histograms
c1 = ROOT.TCanvas("c1","PWA parameter fitting")
# Calculate, how many pads do we need
num_pads_per_row = 2.0
TH2D_len = int(np.ceil((len(TH2D_names))/num_pads_per_row))
TH1D_len = int(np.ceil(len(TH1D_names)/num_pads_per_row))
c1.Divide(2,TH2D_len + TH1D_len,0,0)

# Plot complex parameters
for k in range(len(TH2D_names)):
    try:
        c1.cd(k + 1)
        exec(TH2D_names[k] + ".Draw(\"colz\")")

        TH2D_x_lines[k].SetLineColor(46)
        TH2D_y_lines[k].SetLineColor(46)

        TH2D_x_lines[k].Draw()
        TH2D_y_lines[k].Draw()
    except IndexError:
        pass

# One parameter was fixed as reference parameter - lets plot it as well
#c1.cd(len(TH2D_names)+1)
#ref_hist.Draw()
##reference_pad = c1.GetPad(len(TH2D_names)+1)
##reference_pad.SetName("Fixed amplitude - theta_rho_770 = 1.0 + 0.0j")
#fixed_x_line = ROOT.TLine(1.0,y_min,1.0,y_max)
#fixed_y_line = ROOT.TLine(x_min,1.0,x_max,1.0)
#fixed_x_line.SetLineColor(46)
#fixed_y_line.SetLineColor(46)
#fixed_x_line.Draw()
#fixed_y_line.Draw()
#

# Plot real-valued parameters
for k in range(len(TH1D_names)):
    try:
        c1.cd(TH2D_len * 2 + 1 + k)
        exec(TH1D_names[k] + ".Draw()")

        TH1D_lines[k].SetLineColor(46)
        TH1D_lines[k].Draw()

    except IndexError:
        pass


# The saves do not really work for some reason
## Save the data in f_out
#for name in TH2D_names:
#    exec("f_out.WriteObject(" + name + ", \"" + str(name) + "\")")
#for name in TH1D_names:
#    exec("f_out.WriteObject(" + name + ", \"" + str(name) + "\")")
#f_out.WriteObject(c1, "c1")
c1.SaveAs("nice_output_canvas.root")

f_out.Write()
f_in.Close()
f_out.Close()
         

    
