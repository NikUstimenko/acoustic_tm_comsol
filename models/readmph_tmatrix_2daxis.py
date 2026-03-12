'''
Read a solved .mph file and produce a .hdf5 file containing the T-matrix and related data.
This script uses the unofficial 'mph' library.

Some of this is adapted Dominik Beutel's original template_comsol.py.

Author: Markus Nyman, 2024

Adapted for acoustic T-matrices of axisymmetric objects by Nikita Ustimenko, 2024.
'''

import numpy as np
import mph
import sys
import h5py
import tmatrix_tools
import argparse

dsname = 'Study 1//Parametric Solutions 1'

parser = argparse.ArgumentParser()
parser.add_argument('modelfile')
parser.add_argument('tmatrixfile')
args = parser.parse_args()

print('Starting client ...')
client = mph.start()
print('Loading model file ...')
model = client.load(args.modelfile)

print('Getting values ...')
# T-matrix coefficients and indexing of the multipoles
m_in = np.int64(model.evaluate('m_in', dataset=dsname))
freq = model.evaluate('freq', dataset=dsname)
uniquefreq, index_uniquefreq = np.unique(freq, return_index=True)


ap = model.evaluate('ap', dataset=dsname)
#am = model.evaluate('am', dataset=dsname)

l_max = np.max(m_in)
n_dim = (l_max+1)**2 #Dimension of T-matrix
it = 0
l_out_eff = np.zeros(n_dim**2 * len(uniquefreq), np.int64)
l_in_eff = np.zeros(n_dim**2 * len(uniquefreq), np.int64)
m_out_eff = np.zeros(n_dim**2 * len(uniquefreq), np.int64)
m_in_eff = np.zeros(n_dim**2 * len(uniquefreq), np.int64)
freq_eff = np.zeros(n_dim**2 * len(uniquefreq), np.float64)


for i in range(l_max + 1):
    for j in range(l_max + 1):
        for l in range(len(uniquefreq)):
            for k in range(-i, i + 1):
                for s in range(-j, j + 1):
                    l_out_eff[it] = i
                    l_in_eff[it] = j
                    m_out_eff[it] = k
                    m_in_eff[it] = s
                    freq_eff[it] = uniquefreq[l]
                    it += 1

dm = 0
with open('tmatrix_coeffs.txt', 'w') as f:
    f.write('% l_in   l_out   m_in  m_out  freq   ap (1)\n')
    for i in range(len(l_in_eff)):
        if m_in_eff[i] == m_out_eff[i]:
            f.write(f'{l_in_eff[i]}   {l_out_eff[i]}   {m_in_eff[i]}  {m_out_eff[i]}  {freq_eff[i]}   {ap[dm]}\n')
            dm += 1
        else:
            f.write(f'{l_in_eff[i]}   {l_out_eff[i]}   {m_in_eff[i]}  {m_out_eff[i]}  {freq_eff[i]}   {0.+0.j}\n')

(
    tmats,
    _,
    ls,
    ms,
    ls_inc,
    ms_inc,
    params,
) = tmatrix_tools.extract_tmatrix_comsol("tmatrix_coeffs.txt")
if np.all(ls == ls_inc) and np.all(ms == ms_inc):
    modes_inc = None
else:
    modes_inc = ls_inc, ms_inc



"""
# Saving the parameters in the hdf5 might be useful to some
#params = model.parameters()
#descriptions = model.descriptions()

# Materials
materials = (model / 'materials').children()
material_tags = [m.tag() for m in materials]
material_names = [m.name() for m in materials]
material_descriptions = ['']*len(materials)
material_keywords = ['']*len(materials)

# The output is for every parameter case, so we just pick one from each
# different frequency.
rho = [model.evaluate(f'{mt}.def.rho', dataset=dsname)[index_uniquefreq] for mt in material_tags]
c = [model.evaluate(f'{mt}.def.c', dataset=dsname)[index_uniquefreq] for mt in material_tags]
#kappa = [model.evaluate(f'{mt}.pg1.kappa', dataset=dsname)[index_uniquefreq] for mt in material_tags]

# Mesh export is defined in the model but running it automatically seems
# to be problematic. So we do it here.
# Dominik's script doesn't work though, so what can we do.
#model.property('exports/Mesh 1', 'filename', 'mesh1.mphtxt')
#model.export('Mesh 1')
"""

print('Building the hdf5 file ...')
# Now put it all in the hdf5 file.
with h5py.File(args.tmatrixfile, 'w') as f:
    tmatrix_tools.base_data(
            f,
            tmatrices=tmats,
            name="Sphere",
            description="Spherical object for testing.",
            keywords="passive",
            freqs=freq,
            ftype="frequency",
            funit='Hz',
            modes=(ls, ms),
            modes_inc=modes_inc,
            format_version="v0.0.1-5-g7a6c03f",
        )
"""
    f.create_group('materials')
    embedding = True
    for name, description, keywords, rho, c in zip(
        material_names,
        material_descriptions,
        material_keywords,
        rho,
        c,
    ):
        f.create_group(f"materials/{name.lower()}")
        tmatrix_tools.isotropic_material(
            f[f"materials/{name.lower()}"],
            name=name,
            description=description,
            keywords=keywords,
            rho=rho,
            c=c,
            embedding=embedding,
        )
        embedding = False
    f.create_group("computation/files")
    tmatrix_tools.computation_data(
        f["computation"],
        name="COMSOL",
        description="",
        keywords="FEM",
        #program_version="comsol=6.1.0.522",
        program_version=model.version(),
        meshfile=None,
    )
"""

print('Done.')
