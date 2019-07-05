import ROOT
from math import sqrt
from array import array

def efferrPoiss(n,d):
   eff = (n/d)
   err = sqrt(eff*(1-eff)/d)
   #sqrt((1./n) + (1./d))
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
   for t in tmp: d[float(t[0])] = [float(t[1]), float(t[2])]
   return d

def main(lines,chtype):
   d = dictmaker(lines[1:])
   #print d
   eff, err = effdicts(d)
   #print eff
   #print err
   g = ROOT.TGraphErrors(len(err), array('f', eff.keys()), array('f', eff.values()), array('f', [0]*len(err)), array('f', err.values()))
   g.GetXaxis().SetTitle("B1wt")
   g.GetYaxis().SetTitle("buckets algo. efficiency")
   g.SetMarkerSize(1)
   g.SetMarkerStyle(8)
   #g.SetMarkerColor(2)
   g.SetTitle(chtype)
   return g

if __name__=="__main__":
   #chtype = "notallhad"
   smearwidth=[0.8, 1.0, 1.2]
   ROOT.gROOT.SetBatch(1)
   c = ROOT.TCanvas( 'c1', 'c1', 200, 10, 600, 400 )
   leg = ROOT.TLegend(0.65, 0.8, 0.90, 0.9)
   counter = 0
   chtype1 = "notallhad_correctB2masstarget_fakeleptonunsmeared"
   #chtype1 = "notallhad_correctB2masstarget"
   #ROOT.gPad.SetLogx(1)
   mg = ROOT.TMultiGraph()
   for sw in smearwidth:
      chtype = "notallhadsw%s_correctB2masstarget_fakeleptonunsmeared"%str(sw)
      #chtype = "notallhadsw%s_correctB2masstarget"%str(sw)
      f = open("effB1wt%s.txt"%chtype, 'r')
      eb1wt = f.readlines()
      g = main(eb1wt, chtype)
      g.SetMarkerColor(counter + 2)
      mg.Add(g, "p")
      leg.AddEntry(g, str(sw), "p")
      counter += 1
   mg.Draw("a")
   mg.GetXaxis().SetTitle("B1wt")
   mg.GetYaxis().SetTitle("buckets algo. efficiency")
   leg.Draw()
   c.Print("effvsB1wt%s.eps"%(chtype1))
   #mg.Draw("a")
   #c.Print("effvsB1wt%s_logx.eps"%(chtype1))
