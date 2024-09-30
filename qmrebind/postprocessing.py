"""
Extract information from ORCA outputs and change the necessary files.
"""

import os
import shutil

import numpy as np
import pandas as pd

import qmrebind.qmrebind_base as base
import qmrebind.defaults as defaults

def get_qm_charges(
    orca_out_file, qm_charge_file, input_pdb, ligand_indices, qm_charge_scheme
):

    """
    Extract the charges of the QM region of the system
    by parsing the output file of the ORCA QM/QM2/MMM calculations.

    Parameters
    ----------
    orca_out_file: str
        User-defined ORCA output file.

    qm_charge_file: str
        User-defined text file where the charges for the
        atoms of the QM region are saved.

    input_pdb: str
        User-defined PDB file.

    ligand_indices: list
        list of ligand residue indices.

    qm_charge_scheme: str
        Charge scheme for the calculation of QM charges for
        the QM region. The options include Hirshfeld,
        CHELPG, Mulliken, and the Loewdin charge calculation
        methods.

    """
    to_begin = None
    to_end = None
    print("Writing qm charge file:", qm_charge_file)
    with open(orca_out_file, "r") as f:
        lines = f.readlines()
    if qm_charge_scheme == "HIRSHFELD":
        for i, line in enumerate(lines):
            if "HIRSHFELD ANALYSIS" in line:
                to_begin = int(i)
                to_end = len(ligand_indices)
        assert to_begin is not None, \
            "Charges not found in ORCA output: An error likely "\
            "occurred in ORCA. See earlier output."
        charges = lines[to_begin + 7 : to_begin + 7 + to_end]
        with open("temp.txt", "w") as f:
            for charge in charges:
                f.write(charge)
        df_temp = pd.read_csv("temp.txt", delimiter=r"\s+", header=None)
        df_temp.columns = ["Index", "Atom", "Charge", "Spin"]
        df_temp["Charge"].to_csv(qm_charge_file, index=False, header=False, 
                                 sep=" ")
        base.delete_files(["temp.txt"])
    
    elif qm_charge_scheme == "CHELPG":
        for i, line in enumerate(lines):
            if "CHELPG Charges       " in line:
                to_begin = int(i)
            if "CHELPG charges calculated" in line:
                to_end = int(i)
        assert to_begin is not None, \
            "Charges not found in ORCA output: An error likely "\
            "occurred in ORCA. See earlier output."
        
        charges = lines[to_begin + 2 : to_end - 4]
        charge_list = []
        for charge in charges:
            charge_list.append(charge.strip().split())
        charge_list_value = []
        for charge in charge_list:
            charge_list_value.append(charge[3])
        data = list(charge_list_value)
        df_charge = pd.DataFrame(data, columns=["Charge"])
        df_charge.to_csv(qm_charge_file, index=False, header=False, sep=" ")
        
    elif qm_charge_scheme == "MULLIKEN":
        for i, line in enumerate(lines):
            if "MULLIKEN ATOMIC CHARGES" in line:
                to_begin = int(i)
                to_end = len(ligand_indices)
        assert to_begin is not None, \
            "Charges not found in ORCA output: An error likely "\
            "occurred in ORCA. See earlier output."
        charges = lines[to_begin + 2 : to_begin + 2 + to_end]
        with open("temp.txt", "w") as f:
            for charge in charges:
                f.write(charge)
        df_temp = pd.read_csv("temp.txt", sep=":", header=None)
        df_temp.columns = ["Index_Atom", "Charge"]
        df_temp["Charge"].to_csv(qm_charge_file, index=False, header=False, 
                                 sep=" ")
        base.delete_files(["temp.txt"])
        
    elif qm_charge_scheme == "LOEWDIN":
        for i, line in enumerate(lines):
            if "LOEWDIN ATOMIC CHARGES" in line:
                to_begin = int(i)
                to_end = len(ligand_indices)
        assert to_begin is not None, \
            "Charges not found in ORCA output: An error likely "\
            "occurred in ORCA. See earlier output."
        charges = lines[to_begin + 2 : to_begin + 2 + to_end]
        with open("temp.txt", "w") as f:
            for charge in charges:
                f.write(charge)
        df_temp = pd.read_csv("temp.txt", sep=":", header=None)
        df_temp.columns = ["Index_Atom", "Charge"]
        df_temp["Charge"].to_csv(qm_charge_file, index=False, header=False, 
                                 sep=" ")
        base.delete_files(["temp.txt"])
    else:
        raise Exception(f"Invalid qm charge scheme: {qm_charge_scheme}")
    
    return

def get_ff_charges(forcefield_file, ff_charges_file, input_pdb):

    """
    Extract the charges of the QM region of the system
    by parsing the output file of the original forcefield file.

    Parameters
    ----------
    forcefield_file: str
        User-defined topology file (prmtop/parm7 file).

    ff_charges_file: str
        User-defined file where the charges defined in the
        topology file (prmtop/parm7 file) is stored in
        a text file.

    input_pdb: str
        User-defined PDB file.

    """
    to_begin = None
    with open(forcefield_file, "r") as f:
        lines = f.readlines()
    num_atoms = base.get_number_pdb_atoms(input_pdb=input_pdb)
    if num_atoms % 5 == 0:
        lines_to_select = int(num_atoms / 5)
    else:
        lines_to_select = int(num_atoms // 5 + 1)
    for i, line in enumerate(lines):
        if "%FLAG CHARGE" in line:
            to_begin = int(i)
            to_end = int(i) + lines_to_select
    assert to_begin is not None, "CHARGE flag not found: invalid parm7 file?"
    print("to_begin:", to_begin, "to_end:", to_end)
    charges = lines[to_begin + 2 : to_end + 2]
    charge_list = []
    for charge in charges:
        charge_list.append(charge.strip().split())
    charge_list = [item for sublist in charge_list for item in sublist]
    df_charge = pd.DataFrame(charge_list, columns=["Charge"])
    print("Writing ff charge file:", ff_charges_file)
    df_charge.to_csv(ff_charges_file, index=False, header=False, sep=" ")
    return

def get_ff_qm_charges(
    qm_charge_file,
    ff_charges_file,
    ff_charges_qm_fmt_file,
    input_pdb,
    ligand_indices,
):

    """
    Use the previous file that stored the charges of the
    the system obtained from the AMBER force field file, converts
    it into AMBER format, i.e. %FORMAT(5E16.8), and then replaces
    the charges of the QM region with the charges of the QM region
    obtained from the ORCA QM/QM2/MM calculations.

    Parameters
    ----------
    qm_charge_file: str
        User-defined text file where the charges for the
        atoms of the QM region are saved.

    ff_charges_file: str
        User-defined file where the charges defined in the
        topology file (prmtop/parm7 file) is stored in
        a text file.

    ff_charges_qm_fmt_file: str
        User-defined file where the replaced charges in the
        AMBER format is stored.

    input_pdb: str
        User-defined PDB file.

    ligand_resname: str
        Three-letter name for the ligand residue.

    """
    df_qm_charges = pd.read_csv(qm_charge_file, header=None, delimiter=r"\s+")
    df_qm_charges.columns = ["Charge"]
    df_qm_charges["Index"] = ligand_indices
    df_qm_charges = df_qm_charges[["Index", "Charge"]]
    qm_charge_list = df_qm_charges["Charge"].values.tolist()
    # charge conversion factor=18.2223 Units???
    qm_charge_list = [
        i * defaults.CHARGE_CONVERSION for i in qm_charge_list
    ]  
    qm_charge_index_list = df_qm_charges["Index"].values.tolist()
    dict_qm = {
        qm_charge_index_list[i]: qm_charge_list[i]
        for i in range(len(qm_charge_index_list))
    }
    df_ff_charges = pd.read_csv(ff_charges_file, header=None, delimiter=r"\s+")
    df_ff_charges.columns = ["Charge"]
    ff_charges = df_ff_charges["Charge"].values.tolist()
    for i in qm_charge_index_list:
        ff_charges[i] = dict_qm[i]
    ff_charges_fmt = []
    for group in base.group_str(ff_charges, 5):
        ff_charges_fmt.append(
            " ".join("{: .8E}".format(x) for x in group if x is not None) + "\n"
        )
    print("Writing qm charges in ff format:", ff_charges_qm_fmt_file)
    with open(ff_charges_qm_fmt_file, "w") as f:
        for i in ff_charges_fmt:
            f.write(" " + i)
    return

def get_qmrebind_parm(forcefield_file, input_pdb, ff_charges_qm_fmt_file):

    """
    Replace the defined charges in the topology file
    (parm7/prmtop) file with the new set of charges obtained from
    the ORCA QM/QM2/MM calculations.

    Parameters
    ----------
    forcefield_file: str
        User-defined topology file (prmtop/parm7 file).

    input_pdb: str
        User-defined PDB file.

    ff_charges_qm_fmt_file: str
        User-defined file where the replaced charges in the
        AMBER format is stored.

    """
    ff_base = os.path.splitext(forcefield_file)[0]
    ff_ext = os.path.splitext(forcefield_file)[1]
    ff_chg_replacement = f"{ff_base}_before_charge_replacement{ff_ext}"
    shutil.copyfile(forcefield_file, ff_chg_replacement)
    with open(forcefield_file, "r") as f1:
        ff_lines = f1.readlines()
    for i, line in enumerate(ff_lines):
        if "%FLAG CHARGE" in line:
            to_begin = 0
            to_end = int(i)
    parm_lines_a = ff_lines[0 : to_end + 2]
    num_atoms = base.get_number_pdb_atoms(input_pdb=input_pdb)
    if num_atoms % 5 == 0:
        ff_lines_to_select = int(num_atoms / 5)
    else:
        ff_lines_to_select = int(num_atoms // 5 + 1)
    parm_lines_b = ff_lines[to_end + ff_lines_to_select + 2 :]
    with open(ff_charges_qm_fmt_file, "r") as f2:
        qm_ff_lines = f2.readlines()
    print("Replacing charges in file:", forcefield_file)
    with open(forcefield_file, "w") as f:
        for line in parm_lines_a:
            f.write(line)
        for line in qm_ff_lines:
            f.write(line)
        for line in parm_lines_b:
            f.write(line)
            
    return

def get_qmrebind_parm_solvent(input_pdb, forcefield_file, ff_charges_file):

    """
    Replace the defined charges in the original topology
    file (parm7/prmtop) file with the new set of charges obtained from
    the ORCA QM/QM2/MM calculations as defined in the topology file
    with no solvent.

    Parameters
    ----------

    input_pdb: str
        User-defined PDB file.

    forcefield_file: str
        User-defined topology file (prmtop/parm7 file).

    ff_charges_file: str
        User-defined file where the charges defined in the
        topology file (prmtop/parm7 file) is stored in
        a text file.

    """
    input_pdb_base = os.path.splitext(input_pdb)[0]
    ff_base = os.path.splitext(forcefield_file)[0]
    ff_ext = os.path.splitext(forcefield_file)[1]
    ff_chg_base = os.path.splitext(ff_charges_file)[0]
    receptorligand_before_qmmm_pdb = f"{input_pdb_base}_before_qmmm.pdb"
    ff_charges_file_before_qmmm = f"{ff_chg_base}_before_qmmm.txt"
    ff_charges_file_after_qmmm = f"{ff_chg_base}_after_qmmm.txt"
    forcefield_file_before_qmmm = f"{ff_base}_before_qmmm{ff_ext}"
    len_pdb_before_qmmm = base.get_number_pdb_atoms(
        input_pdb=receptorligand_before_qmmm_pdb)
    len_pdb_after_qmmm = base.get_number_pdb_atoms(input_pdb=input_pdb)
    get_ff_charges(
        forcefield_file=forcefield_file_before_qmmm,
        ff_charges_file=ff_charges_file_before_qmmm,
        input_pdb=receptorligand_before_qmmm_pdb,
    )
    get_ff_charges(
        forcefield_file=forcefield_file,
        ff_charges_file=ff_charges_file_after_qmmm,
        input_pdb=input_pdb,
    )
    list_ff_charges_file_before_qmmm \
        = list(np.loadtxt(ff_charges_file_before_qmmm))
    print("ff_charges_file_after_qmmm:", ff_charges_file_after_qmmm)
    list_ff_charges_file_after_qmmm \
        = list(np.loadtxt(ff_charges_file_after_qmmm))
    for i, file_after in enumerate(list_ff_charges_file_after_qmmm):
        list_ff_charges_file_before_qmmm[i] = file_after
    list_ff_charges_file_before_qmmm_8e_format = []
    for charge in list_ff_charges_file_before_qmmm:
        list_ff_charges_file_before_qmmm_8e_format.append(
            "{:.8E}".format(charge))
    df_charge = pd.DataFrame(
        list_ff_charges_file_before_qmmm_8e_format, columns=["Charge"])
    ff_charges_file_after_qmmm_8e_format = (
        f"{ff_chg_base}_after_qmmm_8e_format.txt")
    df_charge.to_csv(
        ff_charges_file_after_qmmm_8e_format, index=False, header=False, 
        sep=" ")
    df_ff_charges = pd.read_csv(
        ff_charges_file_after_qmmm_8e_format, header=None, delimiter=r"\s+")
    df_ff_charges.columns = ["Charge"]
    ff_charges = df_ff_charges["Charge"].values.tolist()
    ff_charges_fmt = []
    for group in base.group_str(ff_charges, 5):
        ff_charges_fmt.append(
            " ".join("{: .8E}".format(x) for x in group if x is not None) + "\n"
        )
    ff_charges_file_after_qmmm_5_8e_format = (
        f"{ff_chg_base}_after_qmmm_5_8e_format.txt")
    with open(ff_charges_file_after_qmmm_5_8e_format, "w") as f:
        for charge in ff_charges_fmt:
            f.write(f" {charge}")
            
    no_solvent_ff_filename = f"{ff_base}_no_solvent{ff_ext}"
    shutil.copyfile(forcefield_file, no_solvent_ff_filename)
    with open(forcefield_file_before_qmmm, "r") as f1:
        ff_lines = f1.readlines()
    for i, line in enumerate(ff_lines):
        if "%FLAG CHARGE" in line:
            to_begin = 0
            to_end = int(i)
    parm_lines_a = ff_lines[0 : to_end + 2]
    num_atoms = base.get_number_pdb_atoms(input_pdb=receptorligand_before_qmmm_pdb)
    if num_atoms % 5 == 0:
        ff_lines_to_select = int(num_atoms / 5)
    else:
        ff_lines_to_select = int(num_atoms // 5 + 1)
    parm_lines_b = ff_lines[to_end + ff_lines_to_select + 2 :]
    with open(ff_charges_file_after_qmmm_5_8e_format, "r") as f2:
        qm_ff_lines = f2.readlines()
    print("Writing new ff file with solvent:", forcefield_file)
    with open(forcefield_file, "w") as f:
        for line in parm_lines_a:
            f.write(line)
        for line in qm_ff_lines:
            f.write(line)
        for line in parm_lines_b:
            f.write(line)
    
    return
