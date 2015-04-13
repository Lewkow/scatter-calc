from matplotlib import pyplot as plt
import numpy as np
import matplotlib
import argparse

## Read file
## Return header
def reader(filename):
  i = 0
  data = {}
  header_lines = []
  content_end_lines = []
  with open(filename) as f:
    raw = f.readlines()  

  for r in raw:
    if "-" in r:
      header_lines.append(i)
      if len(header_lines) > 1:
        content_end_lines.append(i)
    i += 1

  for h in range(len(header_lines)):
    try:
      content = raw[header_lines[h]+1:content_end_lines[h]]
    except:
      content = raw[header_lines[h]+1:]
    x = []; y = [];
    for c in content:
      x.append(float(c.split(",")[0])) 
      y.append(float(c.split(",")[1]))
    data[raw[header_lines[h]]] = (x,y)
  return data

def header_parse(header):
  try:
    hs        = header.split("-")
    calc_type = hs[0].rstrip()
    energy    = "{:.3e}".format(float(hs[1]))
    proj      = hs[2].rstrip()
    targ      = hs[3].rstrip()
    pot_type  = hs[4].rstrip()
  except:
    calc_type = None; energy = None; proj = None; targ = None; pot_type = None;
  return (calc_type, energy, proj, targ, pot_type)

def plot_scatter(filename="./data/scatter_calc_out.dat"):
  data = reader(filename)
  for h in data.keys():
    ct, e, p, t, pot = header_parse(h)
    filename = "./plots/"+str(ct)+"_"+str(e)+"eV_"+str(p)+"-"+str(t)+"_"+str(pot)
    (x,y) = data[h]
    if ct == "dcs":
      plot_dcs(filename, h, x, y)

def plot_dcs(outfile, header, x, y, show_plot=False, save_plot=True):
  ct, e, p, t, pot = header_parse(header)
  fig = plt.figure()
  plt.semilogy(x, y, linewidth=3.0)
  plt.title(str(e)+" eV "+str(p)+"-"+str(t)+" "+pot)
  plt.xlabel("Scattering Angle [rad]")
  plt.ylabel("Differential Cross Section [a0^2]")
  if show_plot is True:
    plt.show()
  if save_plot is True:
    outf = outfile.split(".dat")[0]+".png"
    fig.savefig(outf, dpi=fig.dpi)


def main():
  plot_scatter()

if __name__ == "__main__":
  parser = argparse.ArgumentParser()
  parser.add_argument('-f', action='store', 
                            dest='filename', 
                            default='./data/scatter_calc_out.dat',
                            help='Filename to plot')

  results = parser.parse_args()
  if results.filename is not None:
    plot_scatter(results.filename)
  else:
    plot_scatter()

