"""
All functions related to ORCA runs: inputs and running
"""

import os

from biopandas.pdb import PandasPdb

import qmrebind.qmrebind_base as base

# TODO: return to this function to determine what optimizations, if any, may
# be done to it.
def prepare_orca_pdb(
    input_pdb, ligand_pdb, orca_pdb, ligand_indices, receptor_pdb, 
    cut_off_distance):

    """
    ORCA QM/QM2/MM input file reads the PDB file. The occupancy
    column of the PDB file needs to be readjusted according to
    the definition of the QM/QM2/MM layer. An occupancy value of
    1.00 is defined for the QM region, 2.00 for the QM2 region,
    and 0.00 for the MM region. The function prepares the PDB file
    needed for the QM/QM2/MM calculations by saving a PDB file with
    the above considerations.

    Parameters
    ----------
    input_pdb: str
        User-defined PDB file.

    ligand_pdb: str
        User-defined ligand PDB file.

    orca_pdb: str
        User-defined ORCA PDB file.

    ligand_indices: str
        A list of atom indices representing the ligand.

    recepor_pdb: str
        User-defined receptor PDB coordinate file.

    cut_off_distance: float
        Cut-off distance for the QM2 region within the
        vicinity of the QM region.

    """
    ppdb = PandasPdb().read_pdb(input_pdb)
    # Use occupancy value of 1.00 for the QM region and 2.00 for the QM2 region
    ppdb.df["ATOM"].loc[0:, "occupancy"] = 0.00
    for ligand_index in ligand_indices:
        ppdb.df["ATOM"].loc[ligand_index, "occupancy"] = 1.00
        if ppdb.df["ATOM"].loc[ligand_index, "atom_name"].startswith("Cl"):
            ppdb.df["ATOM"].loc[ligand_index, "element_symbol"] = "Cl"
        
        if ppdb.df["ATOM"].loc[ligand_index, "atom_name"].startswith("CL"):
            ppdb.df["ATOM"].loc[ligand_index, "element_symbol"] = "Cl"
        
    receptor_residues, receptor_indices = base.get_indices_qm2_region(
        ligand_pdb=ligand_pdb, input_pdb=input_pdb, 
        cut_off_distance=cut_off_distance, ligand_indices=ligand_indices
    )
    for receptor_index in receptor_indices:
        ppdb.df["ATOM"].loc[receptor_index, "occupancy"] = 2.00
    
    ppdb.to_pdb(path=orca_pdb, records=None, gz=False, append_newline=True)
    return

def get_amber_to_orca_prms(forcefield_file):

    """
    Use the ORCA built-in command to convert
    the AMBER topology file, i.e., parm7/prmtop files to the
    ORCA forcefield file.

    Parameters
    ----------
    forcefield_file: str
        User-defined topology file (prmtop/parm7 file)

    """
    command = f"orca_mm -convff -AMBER {forcefield_file}"
    print(f"Running command: {command}")
    os.system(command)
    return


def get_orca_input(
    nprocs,
    maxiter,
    qm_method,
    qm_basis_set,
    qm2_method,
    optimization,
    frequency_calculation,
    qm_charge_scheme,
    qm2_charge_scheme,
    qm2_charge,
    qm2_mult,
    forcefield_file,
    input_pdb,
    ligand_indices,
    orca_pdb,
    orca_input_file,
    ligand_pdb,
    receptor_pdb,
    cut_off_distance,
    qm_charge,
    qm_mult,
):

    """
    Create an input file for the ORCA QM/QM2/MM
    calculation. It takes into account all the keywords required
    for performing the multiscale calculation.

    Parameters
    ----------
    nprocs: int
        Number of processors to be used for the MPI
        enabled calculation.

    maxiter: int
        Maximum number of iterations needed for the self
        consistent field to converge.

    qm_method: str
        QM theory to be implemented for the high layer,
        i.e., the QM layer. The options for the QM method
        include the perturbation method, reference method,
        local gradient corrected method, hybrid method,
        and meta GGA/hybrid meta GGA method.

    qm_basis_set: str
        Basis set for the QM level of theory. The options
        include pople basis set, pople polarized basis set,
        pople polarized diffused basis set, correlation
        consistent basis set, DEF2 basis set, DEF2 diffused
        basis set, auxiliary coulomb fitted basis_set, auxiliary
        coulomb fitted and exchange basis set, and auxiliary
        correlation consistent basis set.

    qm2_method: str
        QM theory to be implemented for the lower layer,
        i.e., the QM2 layer. The options for the QM2 method
        include the semiempirical method, the tight binding DFT
        method, and the composite method.

    optimization: str
        The string is set to "True" if the user chooses to
        optimize the geometry of the QM region. If not,
        the string is set to "False".

    frequency_calculation: str
        The string is set to "True" if the user chooses
        a frequency calculation for the QM region. If not,
        the string is set to "False".

    qm_charge_scheme: str
        Charge scheme for the calculation of QM charges for
        the QM region. The options include Hirshfeld,
        CHELPG, Mulliken, and the Loewdin charge calculation
        methods.

    qm2_charge_scheme: str
        Charge scheme for the calculation of QM charges for
        the QM2 region. The options include Hirshfeld,
        CHELPG, Mulliken, and the Loewdin charge calculation
        methods.

    qm2_charge: int
        Charge of the QM2 (low) region. The user does not have
        to provide the charge of the QM region since it is
        calculated by ORCA itself. The user also does not
        have to provide the charge of the MM region since it is
        determined by the ORCA force field file.

    qm2_mult: int
        The multiplicity of the QM2 (low) region. The user does
        not have to provide the multiplicity of the QM region
        since ORCA itself calculates it. The user also
        does not have to provide the multiplicity of the MM
        region since the ORCA forcefield file determines it.

    forcefield_file: str
        User-defined topology file (prmtop/parm7 file).

    input_pdb: str
        User-defined PDB file.

    ligand_resname: str
        Three-letter name for the ligand residue.

    orca_pdb: str
        User-defined ORCA PDB file.

    orca_input_file: str
        User-defined ORCA input file.

    ligand_pdb: str
        User-defined ligand PDB file.

    receptor_pdb: str
        User-defined receptor PDB coordinate file.

    cut_off_distance: int
        Cut-off distance for the QM2 region within the
        vicinity of the QM region.

    qm_charge: int
        Charge of the QM region.

    qm_mult: int
        Multiplicity of the QM region.

    implicit_solv: str
        The Universal Solvation Model (SMD) method, an implicit
        solvation model an improvement over the Conductor-like
        Polarizable Continuum Model (CPCM), since it uses the
        full solute electron density to compute the cavity-dispersion
        contribution instead of the area only.

    """
    ff_base = os.path.splitext(forcefield_file)[0]
    ff_ext = os.path.splitext(forcefield_file)[1]
    if ff_ext == ".prmtop":
        orcaff_file = f"{ff_base}.ORCAFF.prms"
    else:
        orcaff_file = f"{forcefield_file}.ORCAFF.prms"
        
    line_0 = f"%PAL NPROCS {nprocs} END"
    line_1 = f"%scf MaxIter {maxiter} END"
    line_2 = f"!QM/{qm2_method}/MM {qm_method} {qm_basis_set}"
    if optimization:
        line_2 += " OPT"
        
    if frequency_calculation:
        line_2 += " NUMFREQ"
        
    line_3 = f"!{qm_charge_scheme}"
    line_4 = "%QMMM"
    line_5 = "PRINTLEVEL 5"
    line_6 = f"ORCAFFFilename \"{orcaff_file}\""
    receptor_residues, receptor_indices = base.get_indices_qm2_region(
        ligand_pdb=ligand_pdb, input_pdb=input_pdb, 
        cut_off_distance=cut_off_distance, ligand_indices=ligand_indices
    )
    ligand_input_indices = base.make_string_range(ligand_indices)
    receptor_input_indices = base.make_string_range(receptor_indices)
    
    line_7 = f"QMATOMS {{{ligand_input_indices}}} END"
    line_8 = f"QM2ATOMS {{{receptor_input_indices}}} END"
    line_9 = f"CHARGE_MEDIUM {qm2_charge}"
    line_10 = f"MULT_MEDIUM {qm2_mult}"
    line_11 = f"CHARGE_METHOD {qm2_charge_scheme}"
    line_12 = "Use_QM_InfoFromPDB TRUE"
    line_13 = "Use_Active_InfoFromPDB TRUE END"
    line_14 = f"*PDBFILE {qm_charge} {qm_mult} {orca_pdb} "
    # TODO: do this more efficiently
    if nprocs == 0:
        cmd_list = []
    else:
        cmd_list = [line_0]
    
    cmd_list += [line_1, line_2, line_3, line_4, line_5, line_6, line_7,
                 line_8, line_9, line_10, line_11, line_12, line_13, line_14]
    commands = "\n".join(cmd_list)
    print("Writing ORCA input file:", orca_input_file)
    with open(orca_input_file, "w") as f:
        f.write(commands)
    return

def add_xtb_inputs(
    etemp,
    solvation,
    solvent,
    accuracy,
    xtb_memory,
    xtb_nprocs,
    orca_input_file,
    XTB_add_inputs,
):

    """
    Add additional XTB commands in the ORCA input file
    for QM2 calculations.

    Parameters
    ----------
    etemp: int
        Electronic temperature of the system.

    solvation: str
        Implicit solvation of the system.

    solvent: str
        Solvent for the implicit solvation of the
        system.

    accuracy: int
        Accuracy value for XTB optimization, default
        is ORCA’s accuracy x 1.e6.

    xtb_memory: int
        Memory in MB reserved for XTB calculation.

    xtb_nprocs: int
        Number of processors used for running XTB
        calculation.

    orca_input_file: str
        User-defined ORCA input file.

    XTB_add_inputs: bool
        if True, the function will add additional XTB 
        commands in the ORCA input file, and if False,
        the ORCA input file will remain unchanged.

    """
    line_1 = f"%XTB\n"
    line_2 = f"ETEMP        {etemp}\n"
    line_3 = f"DOALPB       {solvation}\n"
    line_4 = f"ALPBSOLVENT  \"{solvent}\"\n"
    line_5 = f"ACCURACY     {accuracy}\n"
    line_6 = f"MAXCORE      {xtb_memory}\n"
    line_7 = f"NPROCS       {xtb_nprocs}\n"
    line_8 = f"END\n"
    with open(orca_input_file, "r") as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        if "*PDBFILE" in line:
            to_begin = int(i)
    last_line_index = len(lines)-1
    lines_I = lines[: to_begin]
    lines_III = lines[to_begin : last_line_index + 1]
    lines_list = [line_1, line_2, line_3, line_4, line_5, line_6, line_7, 
                  line_8]
    all_lines = "".join(lines_list)
    if XTB_add_inputs:
        with open(orca_input_file, "w") as f:
            for line in lines_I:
                f.write(line)
            f.write(all_lines)
            for line in lines_III:
                f.write(line)
    else:
        with open(orca_input_file, "w") as f:
            for line in lines_I:
                f.write(line)
            for line in lines_III:
                f.write(line)
    
    return

def run_orca_qmmm(orca_dir_pwd, orca_input_file, orca_out_file,
                  skip_running_orca=False):

    """
    Run the ORCA QM/QM2/MM calculations
    for the system.

    Parameters
    ----------
    orca_dir_pwd: str
        PWD of the directory where ORCA is installed.

    orca_input_file: str
        User-defined ORCA input file.

    orca_out_file: str
        User-defined ORCA output file.

    """
    orca_cmd = os.path.join(orca_dir_pwd, "orca")
    command = f"{orca_cmd} {orca_input_file} > {orca_out_file}"
    print("Running command:", command)
    if not skip_running_orca:
        os.system(command)
    assert os.path.exists(orca_out_file), \
        "ORCA output file not found: ORCA likely encountered an error when "\
        "running."
    return