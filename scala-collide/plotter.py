from matplotlib import pyplot as plt
import numpy as np
import matplotlib

def reader(filename):
  i = 0
  with open(filename) as f:
    raw = f.readlines()  
  header = raw[0].rstrip(); content = raw[1:];
  x = []; y = []; 
  for c in content:
    x.append(float(c.split(",")[0])) 
    y.append(float(c.split(",")[1]))
  return (header, x, y)

def plot_scatter(filename="./data/test_dcs.dat"):
  header, x, y = reader(filename)
  if header == "dcs":
    plot_dcs("./plots/test_dcs.png",x, y)

def plot_dcs(outfile, x, y):
  fig = plt.figure()
  plt.semilogy(x, y)
  plt.xlabel("Scattering Angle [rad]")
  plt.ylabel("Differential Cross Section [a0^2]")
  # plt.show()
  fig.savefig(outfile, dpi=fig.dpi)


def main():
  plot_scatter()

if __name__ == "__main__":
  main()