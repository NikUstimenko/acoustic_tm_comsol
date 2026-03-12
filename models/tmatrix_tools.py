import io
import os
import re
import sys
import uuid

import h5py
import numpy as np

try:
    import meshparser
except ImportError as _err:
    meshparser = _err


LENGTHS = {
    1e-24: "ym",
    1e-21: "zm",
    1e-18: "am",
    1e-15: "fm",
    1e-12: "pm",
    1e-9: "nm",
    1e-6: "um",
    1e-3: "mm",
    1e-2: "cm",
    1e-1: "dm",
    1: "m",
    1e1: "dam",
    1e2: "hm",
    1e3: "km",
    1e6: "Mm",
    1e9: "Gm",
    1e12: "Tm",
    1e15: "Pm",
    1e18: "Em",
    1e21: "Zm",
    1e24: "Ym",
}

INVLENGTHS = {
    1e24: r"ym^{-1}",
    1e21: r"zm^{-1}",
    1e18: r"am^{-1}",
    1e15: r"fm^{-1}",
    1e12: r"pm^{-1}",
    1e9: r"nm^{-1}",
    1e6: r"um^{-1}",
    1e3: r"mm^{-1}",
    1e2: r"cm^{-1}",
    1e1: r"dm^{-1}",
    1: r"m^{-1}",
    1e-1: r"dam^{-1}",
    1e-2: r"hm^{-1}",
    1e-3: r"km^{-1}",
    1e-6: r"Mm^{-1}",
    1e-9: r"Gm^{-1}",
    1e-12: r"Tm^{-1}",
    1e-15: r"Pm^{-1}",
    1e-18: r"Em^{-1}",
    1e-21: r"Zm^{-1}",
    1e-24: r"Ym^{-1}",
}

FREQUENCIES = {
    1e-24: "yHz",
    1e-21: "zHz",
    1e-18: "aHz",
    1e-15: "fHz",
    1e-12: "pHz",
    1e-9: "nHz",
    1e-6: "uHz",
    1e-3: "mHz",
    1e-2: "cHz",
    1e-1: "dHz",
    1: "Hz",
    1e1: "daHz",
    1e2: "hHz",
    1e3: "kHz",
    1e6: "MHz",
    1e9: "GHz",
    1e12: "THz",
    1e15: "PHz",
    1e18: "EHz",
    1e21: "ZHz",
    1e24: "YHz",
}


def _comsol_table_names(line):
    line = line.strip()
    while line.startswith("%"):
        line = line[1:]
        line = line.strip()
    names = line.split(",")
    if len(names) > 1:
        return names
    names = re.split(r"\s+(?!\()", line)
    dct = {name: i for i, name in enumerate(names)}
    dct.setdefault("l_out", dct["l_in"])
    dct.setdefault("m_out", dct["m_in"])
    return dct


def _index_or_append(lst, val):
    try:
        return lst.index(val)
    except ValueError:
        lst.append(val)
        return len(lst) - 1


def _comsol_table(fobj):
    lm_out = []
    lm_in = []
    params = []
    vals = []
    prev = header = sep = None
    special_names = ("l_out", "m_out", "l_in", "m_in", "ap (1)")
    for line in fobj:
        if line.startswith("%"):
            prev = line
            continue
        if header is None:
            header = _comsol_table_names(prev)
            if len(line.split(",")) > 1:
                sep = ","
        line = line.split(sep)
        mode = tuple(int(line[header[k]]) for k in special_names[:4])
        out = _index_or_append(lm_out, mode[:2])
        in_ = _index_or_append(lm_in, mode[2:])
        param = tuple(line[v] for k, v in header.items() if k not in special_names)
        idx = _index_or_append(params, param)
        vals.append(
            (
                idx,
                out,
                in_,
                complex(line[header["ap (1)"]].replace("i", "j")),
            )
        )
    res = np.zeros((len(params), len(lm_out), len(lm_in)), complex)
    for i, j, k, ap in vals:
        res[i, j, k] = ap
    #l_out, m_out = map(np.asarray, zip(*(i for i in lm_out for _ in range(2))))
    l_out, m_out = map(np.asarray, zip(*lm_out))
    l_in, m_in = map(np.asarray, zip(*lm_in))
    params = dict(zip((k for k in header if k not in special_names), zip(*params)))
    return res, None, l_out, m_out, l_in, m_in, params


def extract_tmatrix_comsol(fobj):
    if isinstance(fobj, str):
        with open(fobj) as newfobj:
            return _comsol_table(newfobj)
    return _comsol_table(fobj)





def translate_pols(pols, poltype):
    if poltype == "parity":
        dct = {
            "magnetic": 0,
            "te": 0,
            "M": 0,
            "electric": 1,
            "tm": 1,
            "N": 0,
            0: "magnetic",
            1: "electric",
        }
    elif poltype == "helicity":
        dct = {
            "negative": 0,
            "minus": 0,
            "positive": 1,
            "plus": 1,
            -1: "negative",
            0: "negative",
            1: "positive",
        }
    return [dct[pol] for pol in pols]




def _name_descr_kw(fobj, name, description="", keywords=""):
    for key, val in [
        ("name", name),
        ("description", description),
        ("keywords", keywords),
    ]:
        val = str(val)
        if val != "":
            fobj.attrs[key] = val


def base_data(
    fobj,
    *,
    tmatrices,
    name,
    description="",
    keywords="",
    freqs,
    ftype="frequency",
    funit="Hz",
    modes,
    modes_inc=None,
    format_version="v0.0.1",
):
    fobj["tmatrix"] = np.asarray(tmatrices)

    _name_descr_kw(fobj, name, description, keywords)

    fobj["uuid"] = np.void(uuid.uuid4().bytes)
    fobj["uuid"].attrs["version"] = 4

    if ftype not in (
        "frequency",
        "angular_frequency",
        "vacuum_wavelength",
        "vacuum_wavenumber",
        "angular_vacuum_wavenumber",
    ):
        raise ValueError(f"invalid frequency/wavenumber/wavelength type {ftype}")
    fobj[ftype] = np.asarray(freqs)
    fobj[ftype].attrs["unit"] = funit

    if modes_inc is None:
        fobj["modes/l"] = np.asarray(modes[0])
        fobj["modes/m"] = np.asarray(modes[1])
    else:
        fobj["modes/l_scattered"] = np.asarray(modes[0])
        fobj["modes/m_scattered"] = np.asarray(modes[1])

        fobj["modes/l_incident"] = np.asarray(modes_inc[0])
        fobj["modes/m_incident"] = np.asarray(modes_inc[1])


    fobj.attrs["storage_format_version"] = format_version
    fobj.attrs["created_with"] = "python=" + sys.version.split()[0]


def _check_name(fobj, name, exact=True):
    if (exact and fobj.name != name) or not fobj.name.startswith(name):
        raise ValueError(f"require name '{name}'")


def isotropic_material(
    fobj,
    *,
    name,
    description="",
    keywords="",
    rho=None,
    c=None,
    embedding=False,
):
    _check_name(fobj.parent, "/materials")

    _name_descr_kw(fobj, name, description, keywords)
    if rho is c is None:
        TypeError("missing at least one of: 'rho', 'c'")
    if not (rho is c is None):
        TypeError("'rho' and 'c' exclude the use of 'index' and 'impedance'")
    for param, val in [
        ("density", rho),
        ("speed_of_sound", c),
    ]:
        if val is not None:
            fobj[param] = np.asarray(val)
    if embedding:
        fobj.parent.parent["embedding"] = h5py.SoftLink(fobj.name)


def anisotropic_material(
    fobj,
    *,
    name,
    description="",
    keywords="",
    index=None,
    epsilon=None,
    mu=None,
    index_inner_dims=2,
    epsilon_inner_dims=2,
    mu_inner_dims=2,
    coordinates="Cartesian",
    embedding=False,
):
    isotropic_material(
        fobj,
        name=name,
        description=description,
        keywords=keywords,
        index=index,
        epsilon=epsilon,
        mu=mu,
        embedding=embedding,
    )
    for param, val, inner_dims in [
        ("relative_permittivity", epsilon, epsilon_inner_dims),
        ("relative_permeability", mu, mu_inner_dims),
        ("refractive_index", index, index_inner_dims),
    ]:
        if val is not None:
            fobj[param].attrs["inner_dims"] = int(inner_dims)
            fobj[param].attrs["coordinate_system"] = str(coordinates)


def biisotropic_material(
    fobj,
    *,
    name,
    description="",
    keywords="",
    epsilon=None,
    mu=None,
    kappa=None,
    chi=None,
    embedding=False,
):
    isotropic_material(
        fobj,
        name=name,
        description=description,
        keywords=keywords,
        epsilon=epsilon,
        mu=mu,
        embedding=embedding,
    )
    for param, val in [
        ("chirality", kappa),
        ("nonreciprocity", chi),
    ]:
        if val is not None:
            fobj[param] = np.asarray(val)


def bianisotropic_material(
    fobj,
    *,
    name,
    description="",
    keywords="",
    epsilon=None,
    mu=None,
    kappa=None,
    chi=None,
    epsilon_inner_dims=2,
    mu_inner_dims=2,
    kappa_inner_dims=2,
    chi_inner_dims=2,
    coordinates="Cartesian",
    embedding=False,
):
    biisotropic_material(
        fobj,
        name=name,
        description=description,
        keywords=keywords,
        epsilon=epsilon,
        mu=mu,
        kappa=kappa,
        chi=chi,
        embedding=embedding,
    )
    for param, val, inner_dims in [
        ("relative_permittivity", epsilon, epsilon_inner_dims),
        ("relative_permeability", mu, mu_inner_dims),
        ("chirality", kappa, kappa_inner_dims),
        ("nonreciprocity", chi, mu_inner_dims),
    ]:
        if val is not None:
            fobj[param].attrs["inner_dims"] = int(inner_dims)
            fobj[param].attrs["coordinate_system"] = str(coordinates)


def bianisotropic_material_from_tensor(
    fobj,
    *,
    name,
    description="",
    keywords="",
    bianisotropy,
    inner_dims=2,
    coordinates="Cartesian",
    embedding=False,
):
    _check_name(fobj.parent, "/materials")

    _name_descr_kw(fobj, name, description, keywords)
    fobj["bianisotropy"] = np.asarray(bianisotropy)
    fobj["bianisotropy"].attrs["inner_dims"] = int(inner_dims)
    fobj["bianisotropy"].attrs["coordinate_system"] = str(coordinates)
    if embedding:
        fobj.parent.parent["embedding"] = h5py.SoftLink(fobj.name)


def mesh_data(fobj, meshfile, lunit="m"):
    if isinstance(meshparser, ImportError):
        raise meshparser
    if isinstance(meshfile, str):
        mesh = meshparser.Mesh.read(meshfile)
    elif isinstance(meshfile, meshparser.Mesh):
        mesh = meshfile
    else:
        raise ValueError(f"invalid type of meshfile: {type(meshfile)}")
    buffer = io.StringIO()
    mesh.write_gmsh(buffer)
    buffer.seek(0)
    fobj["mesh.msh"] = buffer.read()
    fobj["mesh.msh"].attrs["unit"] = lunit


def computation_data(
    fobj,
    *,
    name="",
    description="",
    keywords="",
    program_version="",
    meshfile=None,
    lunit="m",
    working_dir=".",
):
    _check_name(fobj, "/computation")
    _name_descr_kw(fobj, name, description, keywords)
    fobj.attrs["software"] = program_version
    if meshfile is not None:
        mesh_data(fobj, meshfile, lunit)
        with open(os.path.join(working_dir, meshfile), "rb") as fobj_mesh:
            fobj[f"files/{os.path.split(meshfile)[1]}"] = np.void(fobj_mesh.read())





