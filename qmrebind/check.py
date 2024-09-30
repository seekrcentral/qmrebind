"""
Run a series of checks to ensure that outputs are correct.
"""

import os
import sys

import numpy as np
import pandas as pd
import parmed
import simtk

import qmrebind.qmrebind_base as base
import qmrebind.defaults as defaults

RECURSION_LIMIT = 100000

def _organize_energy_decomp_list(prmtop_energy_decomp_non_params, category):
    """
    A commonly-repeated pattern for organizing energy decomposition
    """
    result = base.list_to_dict(
        [
            item for sublist in [
                list(elem) for elem in prmtop_energy_decomp_non_params
            ]
            for item in sublist
        ]
    ).get(category)
    return result

def get_energy_diff(
        forcefield_file, forcefield_file_before_qm, 
        pdb_file_before_qm):

    """
    Use the PARMED module to compare the total
    energy of the system before and after parameterization.

    Parameters
    ----------
    forcefield_file: str
        User-defined topology file (prmtop/parm7 file).

    input_pdb: str
        User-defined PDB file.

    """
    print("Loading files:", forcefield_file_before_qm, "and", pdb_file_before_qm)
    parm_non_params = parmed.load_file(
        forcefield_file_before_qm, pdb_file_before_qm)
    prmtop_energy_decomposition_non_params \
        = parmed.openmm.energy_decomposition_system(
            parm_non_params, parm_non_params.createSystem())
    prmtop_energy_decomp_non_params_value = [
        _organize_energy_decomp_list(
            prmtop_energy_decomposition_non_params, "HarmonicBondForce"),
        _organize_energy_decomp_list(
            prmtop_energy_decomposition_non_params, "HarmonicAngleForce"),
        _organize_energy_decomp_list(
            prmtop_energy_decomposition_non_params, "PeriodicTorsionForce"),
        _organize_energy_decomp_list(
            prmtop_energy_decomposition_non_params, "NonbondedForce"),
    ]
    prmtop_energy_decomp_non_params_list = [
        "HarmonicBondForce",
        "HarmonicAngleForce",
        "PeriodicTorsionForce",
        "NonbondedForce",
    ]
    df_energy_non_params = pd.DataFrame(
        list(
            zip(
                prmtop_energy_decomp_non_params_list,
                prmtop_energy_decomp_non_params_value,
            )
        ),
        columns=["Energy_term", "Energy_non_params"],
    )
    df_energy_non_params = df_energy_non_params.set_index("Energy_term")
    parm_params = parmed.load_file(forcefield_file, pdb_file_before_qm)
    prmtop_energy_decomposition_params \
        = parmed.openmm.energy_decomposition_system(
            parm_params, parm_params.createSystem())
    prmtop_energy_decomposition_params_value = [
        _organize_energy_decomp_list(
            prmtop_energy_decomposition_params, "HarmonicBondForce"),
        _organize_energy_decomp_list(
            prmtop_energy_decomposition_params, "HarmonicAngleForce"),
        _organize_energy_decomp_list(
            prmtop_energy_decomposition_params, "PeriodicTorsionForce"),
        _organize_energy_decomp_list(
            prmtop_energy_decomposition_params, "NonbondedForce"),
    ]
    prmtop_energy_decomposition_params_list = [
        "HarmonicBondForce",
        "HarmonicAngleForce",
        "PeriodicTorsionForce",
        "NonbondedForce",
    ]
    df_energy_params = pd.DataFrame(
        list(
            zip(
                prmtop_energy_decomposition_params_list,
                prmtop_energy_decomposition_params_value,
            )
        ),
        columns=["Energy_term", "Energy_params"],
    )
    df_energy_params = df_energy_params.set_index("Energy_term")
    df_compare = pd.concat([df_energy_non_params, df_energy_params], axis=1)
    df_compare["Energy_difference"] = df_compare["Energy_params"].sub(
        df_compare["Energy_non_params"], axis=0
    )
    print(df_compare)
    return

def run_openmm_sim(input_pdb, forcefield_file, sim_steps, T):

    """
    Run OpenMM simulation to test the validity of the
    reparameterized topology files with the existing
    PDB file.

    Parameters
    ----------
    input_pdb: str
        User-defined PDB file, which is a copy of the PDB file
        for each anchor.

    forcefield_file: str
        User-defined topology file (prmtop/parm7 file), which is
        a copy of the topology file for each anchor.

    sim_steps: int
        Number of simulation steps for the OpenMM MD simulation.

    T: int
        OpenMM simulation temperature for NVT simulation.

    """
    input_pdb_base = os.path.splitext(input_pdb)[0]
    prmtop = simtk.openmm.app.AmberPrmtopFile(forcefield_file)
    pdb = simtk.openmm.app.PDBFile(input_pdb)
    system = prmtop.createSystem(
        nonbondedCutoff=1 * simtk.unit.nanometer, 
        constraints=simtk.openmm.app.HBonds
    )
    integrator = simtk.openmm.LangevinIntegrator(
        T * simtk.unit.kelvin, 1 / simtk.unit.picosecond, 
        0.002 * simtk.unit.picoseconds
    )
    simulation = simtk.openmm.app.Simulation(
        prmtop.topology, system, integrator)
    simulation.context.setPositions(pdb.positions)
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
    simulation.minimizeEnergy(maxIterations=10000)
    print(simulation.context.getState(getEnergy=True).getPotentialEnergy())
    sim_output = f"{input_pdb_base}_openmm_sim.pdb"
    simulation.reporters.append(
        simtk.openmm.app.PDBReporter(sim_output, int(sim_steps / 10))
    )
    simulation.reporters.append(
        simtk.openmm.app.StateDataReporter(
            sys.stdout,
            int(sim_steps / 10),
            step=True,
            potentialEnergy=True,
            temperature=True,
        )
    )
    simulation.step(sim_steps)
    return

def get_charge_diff_file(forcefield_file, ligand_pdb, ligand_charge_diff_file):

    """
    The function creates a new file where the old charges, new charges,
    and the charge difference for the ligand molecule are saved.

    Parameters
    ----------
    forcefield_file: str
        User-defined topology file (prmtop/parm7 file), which is
        a copy of the topology file for each anchor.

    ligand_pdb: str
        User-defined ligand PDB file.

    ligand_charge_diff_file: str
        User-defined charge difference file.

    """
    ff_base = os.path.splitext(forcefield_file)[0]
    ff_ext = os.path.splitext(forcefield_file)[1]
    non_qm_forcefield_file = f"{ff_base}_before_qmmm{ff_ext}"
    
    command = f"diff {non_qm_forcefield_file} {forcefield_file} > diff.txt"
    print("Running command:", command)
    os.system(command)
    with open("diff.txt", "r") as f:
        lines = f.readlines()
        
    charges_list = []
    for line in lines:
        if len(base.extract_charges_str(line)) == 5:
            charges_list.append(base.extract_charges_str(line))
    list_1 = [
        item
        for sublist in charges_list[0 : int(len(charges_list) / 2)]
        for item in sublist
    ]
    list_2 = [
        item
        for sublist in charges_list[int(len(charges_list) / 2) :]
        for item in sublist
    ]
    subtracted = list()
    for item1, item2 in zip(list_1, list_2):
        item = item1 - item2
        subtracted.append(item)
    list1 = []
    list2 = []
    diff = []
    for i in range(len(subtracted)):
        if subtracted[i] != 0:
            list1.append(list_1[i])
            list2.append(list_2[i])
            diff.append(subtracted[i])
    list1 = [i / defaults.CHARGE_CONVERSION for i in list1]
    list2 = [i / defaults.CHARGE_CONVERSION for i in list2]
    diff = [i / defaults.CHARGE_CONVERSION for i in diff]
    print("The sum of charges before parameterization: ", sum(list1))
    print("The sum of charges after parameterization: ", sum(list2))
    #df1 = pd.read_csv(ligand_pdb, header=None, delimiter=r"\s+")
    #df1 = df1[df1.columns[-1]]
    #elements = df1.values.tolist()
    ligand_parmed = parmed.load_file(ligand_pdb)
    elements = []
    for atom in ligand_parmed.atoms:
        elements.append(atom.name)
    df = pd.DataFrame([elements, list1, list2, diff])
    df = df.transpose()
    df.columns = ["Element", "Before QMMM", "After QMMM", "Charge Diff."]
    print("Writing ligand charge difference file:", ligand_charge_diff_file)
    df.to_csv(ligand_charge_diff_file, header=True, index=None, sep=",", 
              mode="w")
    # TODO: check these?
    return

def recurse_atoms(atom, _visited_indices=set()):
    """
    Recursively visit all atoms within a molecule for the purposes
    of determining the molecules (all sets of atoms connected by bonds)
    in the system.
    """
    
    _visited_indices.add(atom.idx)
    for bonded_atom in atom.bond_partners:
        if not bonded_atom.idx in _visited_indices:
            branch_indices = recurse_atoms(bonded_atom, _visited_indices)
            _visited_indices.update(branch_indices)
    return _visited_indices

def check_ligand_same_molecule(pdb_filename, prmtop_filename, ligand_indices):
    """
    The user might accidentally define atom selections that span
    multiple molecules. Check this possibility by finding all molecules
    in the system and ensure that atom selections only exist on one
    molecule.
    """
    sys.setrecursionlimit(RECURSION_LIMIT)
    warnstr1 = """CHECK FAILURE: the atom selection for ligand 
    is split over multiple molecules. Atom index {}, 
    which has the name {} and serial id {}, was the first atom to 
    be detected in a different molecule than the others. Please
    check the structure to ensure that atom indexing is correct, and that
    there aren't multiple ligands."""
    warnstr2 = """CHECK FAILURE: the atom index {} does not exist in the 
    structure."""
    
    structure = parmed.load_file(prmtop_filename, pdb_filename)
    molecules = []
    traversed_atoms = set()
    for index, atom in enumerate(structure.atoms):
        if index in traversed_atoms:
            continue
        molecule = recurse_atoms(atom, set())
        traversed_atoms.update(molecule)
        molecules.append(molecule)
        
    in_molecule = None
    for atom_index in ligand_indices:
        index_found = False
        for mol_id, molecule in enumerate(molecules):
            if atom_index in molecule :
                if in_molecule is None:
                    in_molecule = mol_id
                    index_found = True
                elif in_molecule == mol_id:
                    index_found = True
                else:
                    print(warnstr1.format(
                        atom_index, structure.atoms[atom_index].name, 
                        structure.atoms[atom_index].number))
                    return False
        if not index_found:
            print(warnstr2.format(atom_index))
            return False
        
    return True

def check_region_charge(prmtop_filename, atom_indices, expected_charge, 
                        region="UNKNOWN"):
    """
    Load a prmtop file and see what the charge of a particular region is. 
    If the regional charge doesn't equal the expected charge, throw an error.
    """
    
    total_charge = base.get_region_charge(prmtop_filename, atom_indices)
    if np.isclose(total_charge, expected_charge, atol=0.01):
        return True
    else:
        warnstr = f"""CHECK FAILURE: the total charge of {region} region was
        expected to be {expected_charge}, but was found to be {total_charge}
        in the original forcefield file."""
        print(warnstr)
        return False