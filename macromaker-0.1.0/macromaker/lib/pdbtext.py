# encoding: utf-8

__doc__ = """

Text manipulation of PDB protein structure files.

This is a utility function providing common operations
applied to PDB files that can easily be done with text
processing.
"""

solvent_res_types = [
    'HOH', 'WAT', 'TIP', 'SOL',
    'CLA', 'SOD', 'NA', 'CL', 
    'NA+', 'CL-', 'Na', 'Cl',
    'Na+', 'Cl-']


def strip_lines(pdb_txt, tag_func):
  new_lines = []
  for line in pdb_txt.splitlines():
    if tag_func(line):
      continue
    new_lines.append(line)
  return '\n'.join(new_lines)


def strip_hydrogens(pdb_txt):

  def strip_space_and_digits(s):
    result = ""
    for c in s:
      if not (c.isdigit() or c is " "):
        result += c
    return result

  new_lines = []
  for line in pdb_txt.splitlines():
    if line.startswith("ATOM"):
      raw_atom_type = line[12:16]
      element = strip_space_and_digits(raw_atom_type)[0]
      if element == "H":
        continue
    new_lines.append(line)
  return '\n'.join(new_lines)


def strip_solvent(pdb_txt):
  new_lines = []
  for line in pdb_txt.splitlines():
    res_type = line[17:20].strip().upper()
    if not res_type in solvent_res_types:
      new_lines.append(line)
  return '\n'.join(new_lines)


def strip_alternative_atoms(pdb_txt):
  new_lines = []
  for line in pdb_txt.splitlines():
    new_line = line
    if line.startswith('ATOM'):
      alt_loc = line[16]
      if not alt_loc in [' ']:
        if alt_loc in ['A', 'a']:
          new_line = line[:16] + ' ' + line[17:]
        else:
          continue
    new_lines.append(new_line)
  return '\n'.join(new_lines)


def clean_pdb(in_pdb, out_pdb):
  fh_r = open(in_pdb, 'r')
  txt = fh_r.read()
  txt = strip_lines(txt, lambda l: l.startswith('HETATM'))
  txt = strip_solvent(txt)
  txt = strip_alternative_atoms(txt)
  txt = strip_hydrogens(txt)
  fh_r.close()

  fh_w = open(out_pdb, 'w')
  fh_w.write(txt)
  fh_w.close()