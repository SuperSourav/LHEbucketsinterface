import ROOT
from math import sqrt
from array import array

def efferrPoiss(n,d):
   err =  sqrt((1./n) + (1./d))
   err = (n/d)*err
   return err

def effdicts(d):
   k = d.keys()
   y = d.values()
   eff = [n/d for (d,n) in y]
   efferr = [efferrPoiss(n,d) for (d,n) in y]
   print y
   print eff
   print efferr
   return dict(zip(k, eff)), dict(zip(k, efferr))

def dictmaker(lines):
   d = {}
   tmp = [l.split() for l in lines]
   for t in tmp: d[int(t[0])] = [float(t[1]), float(t[2])]
   return d

def main(lines,chtype):
   d = dictmaker(lines[1:])
   print d
   eff, err = effdicts(d)
   print eff
   print err
   ROOT.gROOT.SetBatch(0)
   c = ROOT.TCanvas( 'c1', 'c1', 200, 10, 600, 400 )
   g = ROOT.TGraphErrors(len(err), array('f', eff.keys()), array('f', eff.values()), array('f', [0]*len(err)), array('f', err.values()))
   g.GetXaxis().SetTitle("B1wt")
   g.GetYaxis().SetTitle("buckets algo. efficiency")
   g.SetMarkerSize(1)
   g.SetMarkerStyle(8)
   g.SetMarkerColor(2)
   g.SetTitle(chtype)
   g.Draw("ap")
   c.Print("effvsB1wt%s.eps"%(chtype))
   return 

if __name__=="__main__":
   #chtype = "notallhad"
   chtype = "notallhad_correctB2masstarget"
   f = open("effB1wt%s.txt"%chtype, 'r')
   eb1wt = f.readlines()
   main(eb1wt, chtype)
