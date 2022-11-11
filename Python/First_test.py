import ROOT
import math

c1 = ROOT.TCanvas('c1', 'mycanvas')

%jsroot on
c = ROOT.TCanvas()
f1 = ROOT.TF1("func1", "sin(x)", 0, 10)
f1.Draw()
c.Draw() # Necessary to make the graphics show!
c1.Draw()