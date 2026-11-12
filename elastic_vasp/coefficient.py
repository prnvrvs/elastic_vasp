import contextlib
import os
from collections import namedtuple
import re
import warnings
import numpy as np
from glob import glob
import matplotlib.pyplot as plt

from elastic_vasp import poscar


def _find_zero_strain_folder(base_dir="."):
  for folder in glob(os.path.join(base_dir, "strain_*/")):
      label = os.path.basename(folder.rstrip("/")).split("_", 1)[1]
      try:
          if abs(float(label)) < 1e-12:
              return folder.rstrip("/")
      except ValueError:
          continue
  return None


# Energy used for the E(s) fit: "energy(sigma->0)" of the final ionic step in
# OUTCAR (the energy extrapolated to zero smearing). OUTCAR is the only source;
# if it is missing, has no such line, or the run did not finish, post-processing
# stops and says where.
_OUTCAR_SIGMA0 = re.compile(r"energy\(sigma->0\)\s*=\s*(\S+)")
# VASP writes this footer only when a run terminates normally
_OUTCAR_FINISHED = "General timing and accounting"


class MissingEnergyError(RuntimeError):
  """Raised when one or more strain folders have no usable OUTCAR energy."""


# OUTCARs can be hundreds of MB, so during one post-processing run each file is
# scanned once (energy, volume and finished flag together) and the result reused.
_outcar_cache = None


@contextlib.contextmanager
def outcar_cache():
  """Reuse OUTCAR scans inside this block (nesting is fine)."""
  global _outcar_cache
  outermost = _outcar_cache is None
  if outermost:
      _outcar_cache = {}
  try:
      yield
  finally:
      if outermost:
          _outcar_cache = None


OutcarInfo = namedtuple("OutcarInfo", "energy volume finished nsw ibrion relaxed")

_OUTCAR_NSW = re.compile(r"\bNSW\s*=\s*(-?\d+)")
_OUTCAR_IBRION = re.compile(r"\bIBRION\s*=\s*(-?\d+)")
# printed when an ionic relaxation (IBRION = 1, 2, 3) meets EDIFFG
_OUTCAR_RELAXED = "reached required accuracy"
_RELAXATION_IBRION = (1, 2, 3)


def _scan_outcar(outcar_path):
  """Everything needed from OUTCAR, collected in one pass."""
  key = os.path.abspath(outcar_path)
  if _outcar_cache is not None and key in _outcar_cache:
      return _outcar_cache[key]
  energy = volume = nsw = ibrion = None
  finished = relaxed = False
  with open(outcar_path) as handle:
      for line in handle:
          if "energy(sigma->0)" in line:
              match = _OUTCAR_SIGMA0.search(line)
              if match:
                  try:
                      energy = float(match.group(1))
                  except ValueError:
                      pass
          elif "volume of cell" in line:
              try:
                  volume = float(line.split()[-1])
              except (ValueError, IndexError):
                  pass
          elif _OUTCAR_FINISHED in line:
              finished = True
          elif _OUTCAR_RELAXED in line:
              relaxed = True
          elif nsw is None and "NSW" in line:
              match = _OUTCAR_NSW.search(line)
              if match:
                  nsw = int(match.group(1))
          elif ibrion is None and "IBRION" in line:
              match = _OUTCAR_IBRION.search(line)
              if match:
                  ibrion = int(match.group(1))
  result = OutcarInfo(energy, volume, finished, nsw, ibrion, relaxed)
  if _outcar_cache is not None:
      _outcar_cache[key] = result
  return result


def _read_last_outcar_energy(outcar_path):
  """energy(sigma->0) of the final ionic step of a finished (and, for
  relaxations, converged) run."""
  info = _scan_outcar(outcar_path)
  if info.energy is None:
      raise MissingEnergyError(
          "No 'energy(sigma->0)' line found in {} (unfinished or failed run?)".format(outcar_path))
  if not info.finished:
      raise MissingEnergyError(
          "Run did not finish in {}: no '{}' footer (killed, out of walltime or still "
          "running?)".format(outcar_path, _OUTCAR_FINISHED))
  if (info.nsw or 0) > 0 and info.ibrion in _RELAXATION_IBRION and not info.relaxed:
      raise MissingEnergyError(
          "Ionic relaxation did not converge in {}: NSW = {} steps ran out before "
          "'{}' (increase NSW and restart from CONTCAR)".format(outcar_path, info.nsw, _OUTCAR_RELAXED))
  return info.energy


def _read_energy(folder):
  outcar_path = os.path.join(folder, "OUTCAR")
  if not os.path.exists(outcar_path):
      raise MissingEnergyError("OUTCAR not found in {}".format(folder))
  return _read_last_outcar_energy(outcar_path)


def check_energies(folders):
  """Raise MissingEnergyError listing every folder without a usable OUTCAR energy."""
  problems = []
  for folder in folders:
      try:
          _read_energy(folder)
      except MissingEnergyError as exc:
          problems.append(str(exc))
  if problems:
      raise MissingEnergyError(
          "Cannot post-process: {} strain folder(s) have no OUTCAR energy:\n  ".format(len(problems))
          + "\n  ".join(problems))


def _read_volume_from_outcar(outcar_path):
  volume = _scan_outcar(outcar_path).volume
  if volume is None:
      raise ValueError("No volume line found in {}".format(outcar_path))
  return volume


def _read_volume_from_poscar(poscar_path):
  return poscar.volume(poscar_path)


def _read_volume(folder):
  """Cell volume from a folder's OUTCAR, else its POSCAR."""
  outcar_path = os.path.join(folder, "OUTCAR")
  if os.path.exists(outcar_path):
      try:
          return _read_volume_from_outcar(outcar_path)
      except ValueError:
          pass
  poscar_path = os.path.join(folder, "POSCAR")
  if os.path.exists(poscar_path):
      return _read_volume_from_poscar(poscar_path)
  raise FileNotFoundError("Neither OUTCAR nor POSCAR found in {}".format(folder))


def _reference_volume(base_dir, zero_folder):
  """Unstrained volume V0.

  zero_folder (relative to base_dir, or absolute) if given; otherwise the
  strain_0.00-type folder; otherwise the unstrained POSCAR one level up (the
  file pre-processing deformed).
  """
  if zero_folder is not None:
      if not os.path.isabs(zero_folder):
          zero_folder = os.path.join(base_dir, zero_folder)
      return _read_volume(zero_folder)
  found = _find_zero_strain_folder(base_dir)
  if found is not None:
      return _read_volume(found)
  parent_poscar = os.path.join(base_dir, os.pardir, "POSCAR")
  if os.path.exists(parent_poscar):
      return _read_volume_from_poscar(parent_poscar)
  raise FileNotFoundError(
      "No zero-strain folder in {} and no unstrained POSCAR at {}; cannot get the "
      "reference volume (pass zero_folder=...)".format(base_dir, os.path.normpath(parent_poscar)))


def coeff(base_dir=".", zero_folder=None):
  with outcar_cache():
      return _coeff(base_dir, zero_folder)


def _coeff(base_dir, zero_folder):
  # only strain_<value> folders; ignore anything else the user keeps here
  d=sorted(glob(os.path.join(base_dir, "strain_*/")))
  if not d:
      raise FileNotFoundError("No strain_* folders found in {}".format(base_dir))
  if len(d) < 3:
      raise ValueError("{} has only {} strain folder(s); the quadratic E(s) fit needs at "
                       "least 3".format(base_dir, len(d)))
  check_energies(d)
  strains = []
  energies = []
  for folder in d:
      strain_label = os.path.basename(folder.rstrip("/")).split("_", 1)[1]
      strains.append(strain_label)
      energies.append(_read_energy(folder))
  data = np.array(energies) - np.amin(energies)

  Vol = _reference_volume(base_dir, zero_folder)

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
  if SST <= 1e-30:
      # identical energies: the deformation had no effect (e.g. VASP rescaled the
      # cell, or every run returned the same unconverged energy)
      warnings.warn("All energies in {} are identical, so E(s) is flat and the fitted "
                    "coefficient is meaningless. Check these calculations.".format(base_dir),
                    UserWarning, stacklevel=2)
      R2 = float("nan")
  else:
      R2 = SSreg/SST
  image = os.path.join(base_dir, 'plot.jpg')
  lege = 'R^2='+str(R2)
 #=================================
  plt.plot(x, y, "o", new_x, new_y)
  plt.legend([lege])
  plt.savefig(image)
  plt.close()
  return constant
