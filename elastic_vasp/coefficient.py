import os
import numpy as np
from glob import glob
import matplotlib.pyplot as plt


def _find_zero_strain_folder(base_dir="."):
  for folder in glob(os.path.join(base_dir, "strain_*/")):
      label = os.path.basename(folder.rstrip("/")).split("_", 1)[1]
      try:
          if abs(float(label)) < 1e-12:
              return folder.rstrip("/")
      except ValueError:
          continue
  return None


def _read_last_oszicar_energy(oszicar_path):
  last_energy = None
  last_f_energy = None
  with open(oszicar_path) as handle:
      for line in handle:
          if "E0=" not in line and "F=" not in line:
              continue
          parts = line.split()
          for idx, token in enumerate(parts):
              if token.startswith("E0="):
                  try:
                      last_energy = float(token.split("=", 1)[1])
                  except ValueError:
                      continue
              elif token == "F=" and idx + 1 < len(parts):
                  try:
                      last_f_energy = float(parts[idx + 1])
                  except ValueError:
                      continue
  if last_energy is not None:
      return last_energy
  if last_f_energy is not None:
      return last_f_energy
  raise ValueError("No free-energy line found in {}".format(oszicar_path))


def _read_last_outcar_energy(outcar_path):
  free_energy = None
  entropy_energy = None
  with open(outcar_path) as handle:
      for line in handle:
          if "free  energy   TOTEN" in line:
              parts = line.split()
              if len(parts) >= 5:
                  try:
                      free_energy = float(parts[4])
                  except ValueError:
                      continue
          elif "energy  without entropy=" in line and "energy(sigma->0)" in line:
              parts = line.split()
              if len(parts) >= 7:
                  try:
                      entropy_energy = float(parts[6])
                  except ValueError:
                      continue
  if free_energy is not None:
      return free_energy
  if entropy_energy is not None:
      return entropy_energy
  raise ValueError("No energy line found in {}".format(outcar_path))


def _read_energy(folder):
  outcar_path = os.path.join(folder, "OUTCAR")
  if os.path.exists(outcar_path):
      try:
          return _read_last_outcar_energy(outcar_path)
      except ValueError:
          pass
  oszicar_path = os.path.join(folder, "OSZICAR")
  if os.path.exists(oszicar_path):
      return _read_last_oszicar_energy(oszicar_path)
  raise FileNotFoundError("Neither OUTCAR nor OSZICAR found in {}".format(folder))


def _read_volume_from_outcar(outcar_path):
  volume = None
  with open(outcar_path) as handle:
      for line in handle:
          if "volume of cell" not in line:
              continue
          parts = line.split()
          if not parts:
              continue
          try:
              volume = float(parts[-1])
          except ValueError:
              continue
  if volume is None:
      raise ValueError("No volume line found in {}".format(outcar_path))
  return volume


def _read_volume_from_poscar(poscar_path):
  with open(poscar_path) as handle:
      lines = handle.readlines()
  if len(lines) < 5:
      raise ValueError("POSCAR is too short in {}".format(poscar_path))
  scale = float(lines[1].split()[0])
  lattice = np.array([[float(k) for k in lines[i].split()] for i in range(2, 5)])
  if scale < 0:
      return abs(scale)
  return abs(np.linalg.det(lattice)) * scale ** 3


def _read_volume(base_dir, zero_folder):
  if not os.path.isabs(zero_folder):
      zero_folder = os.path.join(base_dir, zero_folder)
  outcar_path = os.path.join(zero_folder, "OUTCAR")
  if os.path.exists(outcar_path):
      try:
          return _read_volume_from_outcar(outcar_path)
      except ValueError:
          pass
  poscar_path = os.path.join(zero_folder, "POSCAR")
  if os.path.exists(poscar_path):
      return _read_volume_from_poscar(poscar_path)
  raise FileNotFoundError("Neither OUTCAR nor POSCAR found in {}".format(zero_folder))


def coeff(base_dir=".", zero_folder=None):
  d=sorted(glob(os.path.join(base_dir, "*/")))
  strains = []
  energies = []
  for folder in d:
      strain_label = os.path.basename(folder.rstrip("/")).split("_", 1)[1]
      strains.append(strain_label)
      energies.append(_read_energy(folder))
  data = np.array(energies) - np.amin(energies)

  if zero_folder is None:
      zero_folder = _find_zero_strain_folder(base_dir) or os.path.join(base_dir, "strain_0.00")
  Vol = _read_volume(base_dir, zero_folder)

  pp=[float(k) for k in strains]
  out=np.stack((np.array(pp), data), axis=-1)
  out2 = out[out[:,0].argsort()]
  np.savetxt(os.path.join(base_dir, 'Result.txt'), out2,fmt='%5.7f')
  x=out2[:,0]
  y=out2[:,1]/Vol
  z = np.polyfit(x,y,2)
  constant=z[0]*160.2176621
  poly=np.poly1d(z)
  new_x = np.linspace(x[0], x[-1])
  new_y = poly(new_x)
  valeur_T=x
  valeur_min=y
  yhat = poly(valeur_T)
  ybar = sum(valeur_min)/len(valeur_min)
  SST = sum((valeur_min - ybar)**2)
  SSreg = sum((yhat - ybar)**2)
  R2 = SSreg/SST
  image = os.path.join(base_dir, 'plot.jpg')
  lege = 'R^2='+str(R2)
 #=================================
  plt.plot(x, y, "o", new_x, new_y)
  plt.legend([lege])
  plt.savefig(image)
  plt.close()
  return constant
