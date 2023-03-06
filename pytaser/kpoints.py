import numpy as np
from pymatgen.core import Structure
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine


def get_kpoint_weights(bandstructure, time_reversal=True, symprec=0.1):
    """
    Function to calculate the kpoint_weights for non-magnetic materials (non-metals).

    Args:
        bandstructure: PMG bandstructure object
        time_reversal:
        symprec: Symmetry precision in Angstrom.(Lower value is more precise, but
            computationally more expensive)
    Returns:
        k-point_weights
    """
    if isinstance(bandstructure, BandStructureSymmLine):
        raise TypeError(
            "This bandstructure object uses a line-mode kpoint mesh instead of a uniform mesh. Kpoint "
            "weighting is not possible"
        )

    kpoints = get_kpoints_from_bandstructure(bandstructure)
    _, _, _, _, _, kp_mapping = expand_kpoints(
        bandstructure.structure,
        kpoints,
        symprec=symprec,
        time_reversal=time_reversal,
        return_mapping=True,
    )
    weights = np.unique(kp_mapping, return_counts=True)[1].astype(float)
    weights /= np.sum(weights)
    return weights


def get_kpoints_from_bandstructure(bandstructure, cartesian=False):
    if cartesian:
        kpoints = np.array([k.cart_coords for k in bandstructure.kpoints])
    else:
        kpoints = np.array([k.frac_coords for k in bandstructure.kpoints])

    return kpoints


def expand_kpoints(
    structure,
    kpoints,
    symprec=0.01,
    return_mapping=False,
    time_reversal=True,
):
    kpoints = np.array(kpoints).round(8)

    # due to limited input precision of the k-points, the mesh is returned as a float
    mesh, is_shifted = get_mesh_from_kpoint_diff(kpoints)

    if is_shifted:
        shift = np.array([1, 1, 1])
    else:
        shift = np.array([0, 0, 0])

    # to avoid issues to limited input precision, recalculate the input k-points
    # so that the mesh is integer and the k-points are not truncated
    # to a small precision
    addresses = np.rint((kpoints + shift / (mesh * 2)) * mesh)
    mesh = np.rint(mesh)
    kpoints = addresses / mesh - shift / (mesh * 2)

    rotations, translations, is_tr = get_reciprocal_point_group_operations(
        structure, symprec=symprec, time_reversal=time_reversal
    )
    len(rotations)
    # rotate all-kpoints
    all_rotated_kpoints = []
    for r in rotations:
        all_rotated_kpoints.append(np.dot(r, kpoints.T).T)
    all_rotated_kpoints = np.concatenate(all_rotated_kpoints)

    # map to first BZ
    all_rotated_kpoints -= np.rint(all_rotated_kpoints)
    all_rotated_kpoints = all_rotated_kpoints.round(8)

    # zone boundary consistent with VASP not with spglib
    all_rotated_kpoints[all_rotated_kpoints == -0.5] = 0.5

    # Find unique points
    unique_rotated_kpoints, unique_idxs = np.unique(
        all_rotated_kpoints, return_index=True, axis=0
    )

    # find integer addresses
    unique_addresses = (unique_rotated_kpoints + shift / (mesh * 2)) * mesh
    unique_addresses -= np.rint(unique_addresses)
    in_uniform_mesh = (np.abs(unique_addresses) < 1e-5).all(axis=1)

    n_mapped = int(np.sum(in_uniform_mesh))
    n_expected = int(np.product(mesh))
    if n_mapped != n_expected:
        raise ValueError(f"Expected {n_expected} points but found {n_mapped}")

    full_kpoints = unique_rotated_kpoints[in_uniform_mesh]
    full_idxs = unique_idxs[in_uniform_mesh]

    if not return_mapping:
        return full_kpoints

    op_mapping = np.floor(full_idxs / len(kpoints)).astype(int)
    kp_mapping = (full_idxs % len(kpoints)).astype(int)

    return full_kpoints, rotations, translations, is_tr, op_mapping, kp_mapping


def get_mesh_from_kpoint_diff(kpoints, ktol=1e-5):
    kpoints = np.array(kpoints)

    # whether the k-point mesh is shifted or Gamma centered mesh
    is_shifted = np.min(np.linalg.norm(kpoints, axis=1)) > 1e-6

    unique_a = np.unique(kpoints[:, 0])
    unique_b = np.unique(kpoints[:, 1])
    unique_c = np.unique(kpoints[:, 2])

    if len(unique_a) == 1:
        na = 1
    else:
        # filter very small changes, with a tol of 5e-4 this means k-point meshes
        # denser than 2000x2000x2000 will be treated as numerical noise. Meshes
        # this dense are extremely unlikely
        diff = np.diff(unique_a)
        diff = diff[diff > ktol]
        na = 1 / np.min(diff[diff > ktol])

    if len(unique_b) == 1:
        nb = 1
    else:
        diff = np.diff(unique_b)
        nb = 1 / np.min(diff[diff > ktol])

    if len(unique_c) == 1:
        nc = 1
    else:
        diff = np.diff(unique_c)
        nc = 1 / np.min(diff[diff > ktol])

    # due to limited precision of the input k-points, the mesh is returned as a float
    return np.array([na, nb, nc]), is_shifted


def get_reciprocal_point_group_operations(
    structure: Structure,
    symprec: float = 0.01,
    time_reversal: bool = True,
):
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

    sga = SpacegroupAnalyzer(structure, symprec=symprec)
    if sga.get_symmetry_dataset() is None:
        # sometimes default angle tolerance doesn't work as expected
        sga = SpacegroupAnalyzer(
            structure, symprec=symprec, angle_tolerance=-1
        )

    rotations = sga.get_symmetry_dataset()["rotations"].transpose((0, 2, 1))
    translations = sga.get_symmetry_dataset()["translations"]
    is_tr = np.full(len(rotations), False, dtype=bool)

    if time_reversal:
        rotations = np.concatenate([rotations, -rotations])
        translations = np.concatenate([translations, -translations])
        is_tr = np.concatenate([is_tr, ~is_tr])

        rotations, unique_ops = np.unique(rotations, axis=0, return_index=True)
        translations = translations[unique_ops]
        is_tr = is_tr[unique_ops]

    # put identity first and time-reversal last
    sort_idx = np.argsort(
        np.abs(rotations - np.eye(3)).sum(axis=(1, 2)) + is_tr * 10
    )

    return rotations[sort_idx], translations[sort_idx], is_tr[sort_idx]
