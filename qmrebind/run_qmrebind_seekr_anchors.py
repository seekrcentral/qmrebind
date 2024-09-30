"""
run_qmrebind_seekr_anchors.py

Run QMrebind for each of a set of specified anchors within a SEEKR
filetree.
"""
import os
import glob
import argparse
from shutil import copyfile

import seekr2.modules.common_base as seekr_base

import qmrebind.qmrebind_base as base
import qmrebind.run_qmrebind_amber as run_qmrebind_amber

QMREBIND_MODEL_GLOB = "model_pre_qmrebind_*.xml"
QMREBIND_MODEL_BASE = "model_pre_qmrebind_{}.xml"

def parse_anchor_list(anchor_str):
    if anchor_str == "":
        anchor_list = []
    else:
        anchor_list = anchor_str.split(",")
        for anchor_index in anchor_list:
            anchor_index_int = int(anchor_index)
            assert anchor_index_int >= 0, "anchor indices must only "\
                "include integers greater than or equal to zero."
    
    return list(map(int, anchor_list))

def get_anchor_prmtop_filename(anchor):
    if anchor.amber_params is not None:
        return anchor.amber_params.prmtop_filename
    
    if anchor.forcefield_params is not None:
        raise Exception("Forcefield input not yet implemented in QMrebind.")
        
    if anchor.charmm_params is not None:
        raise Exception("CHARMM input not yet implemented in QMrebind.")

def set_anchor_prmtop_filename(anchor, filename):
    if anchor.amber_params is not None:
        anchor.amber_params.prmtop_filename = filename
    
    return

def save_new_model(model, save_old_model=True):
    """
    At the end of a QMrebind calculation, generate a new model file. The
    old model file(s) will be renamed with a numerical index.
    
    
    Parameters
    ----------
    model : Model()
        The unfilled Seekr2 Model object.
        
    """
    model_path = os.path.join(model.anchor_rootdir, "model.xml")
    if os.path.exists(model_path) and save_old_model:
        # This is expected, because this old model was loaded
        hidr_model_glob = os.path.join(model.anchor_rootdir, 
                                       QMREBIND_MODEL_GLOB)
        num_globs = len(glob.glob(hidr_model_glob))
        new_pre_hidr_model_filename = QMREBIND_MODEL_BASE.format(num_globs)
        new_pre_hidr_model_path = os.path.join(model.anchor_rootdir, 
                                               new_pre_hidr_model_filename)
        print("Renaming model.xml to {}".format(new_pre_hidr_model_filename))
        copyfile(model_path, new_pre_hidr_model_path)
        
    print("Saving new model.xml")
    old_rootdir = model.anchor_rootdir
    model.anchor_rootdir = "."
    seekr_base.save_model(model, model_path)
    model.anchor_rootdir = old_rootdir
    return

def run_qmrebind_seekr_anchors(
        anchor_list, model, ligand_indices=None, ligand_resname="",
        cut_off_distance=3.0, nprocs=1, maxiter=2000, qm_method="B3LYP", 
        qm_basis_set="6-311G", qm_charge_scheme="CHELPG", qm_charge=None, 
        qm_mult=1, qm2_method="XTB", qm2_charge_scheme="CHELPG", 
        qm2_charge=None, qm2_mult=1, orca_dir_pwd=None, skip_checks=False):
    # parse anchor_list
    if anchor_list.lower() != "all":
        anchor_list_int = parse_anchor_list(anchor_list)
    else:
        anchor_list_int = []
    
    for anchor in model.anchors:
        if not anchor.bulkstate:
            if anchor_list.lower() == "all":
                anchor_list_int.append(anchor.index)
    
    # for each anchor in model, define working directory
    starting_dir = os.path.abspath(".")
    for anchor_index in anchor_list_int:
        assert anchor_index in range(model.num_anchors), \
            f"Anchor indices may only be between 0 and {model.num_anchors}."
        anchor = model.anchors[anchor_index]
        building_dir = os.path.join(
            model.anchor_rootdir, anchor.directory, anchor.building_directory)
        os.chdir(building_dir)
        
        forcefield_file = get_anchor_prmtop_filename(anchor)
        input_pdb = seekr_base.get_anchor_pdb_filename(anchor)
        assert input_pdb is not None, \
            f"Anchor {anchor_index} has empty PDB file."
        assert forcefield_file is not None, \
            f"Anchor {anchor_index} has empty Prmtop file."
        
        new_forcefield_filename = run_qmrebind_amber.run_qmrebind_amber(
            input_pdb, forcefield_file, ligand_indices, ligand_resname,
            cut_off_distance=cut_off_distance, nprocs=nprocs, 
            maxiter=max_iterations, qm_method=qm_method, 
            qm_basis_set=qm_basis_set, qm_charge_scheme=qm_charge_scheme, 
            qm_charge=qm_charge, qm_mult=qm_multiplicity, qm2_method=qm2_method, 
            qm2_charge_scheme=qm2_charge_scheme, qm2_charge=qm2_charge, 
            qm2_mult=qm2_mult, orca_dir_pwd=orca_path, work_dir=None,
            skip_checks=skip_checks)
        
        set_anchor_prmtop_filename(anchor, new_forcefield_filename)
        
    os.chdir(starting_dir)
    return

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument(
        "anchor_list", metavar="ANCHOR_LIST", type=str,
        help="The indices of which anchors to run QMrebind for. If the word "\
        "'all' is entered here, then all anchors with existing PDB structures "\
        "will be included")
    argparser.add_argument(
        "model_file", metavar="MODEL_FILE", type=str, 
        help="The name of model XML file for a SEEKR2 calculation. "\
        "The anchors of this SEEKR model file will be reparametrized with "\
        "qmrebind.")
    argparser.add_argument(
        "-l", "--ligand_indices", dest="ligand_indices", 
        metavar="LIGAND_INDICES", type=str, default="",
        help="A comma-separated list of integers defining ligand within the "\
        "input_pdb structure. Ex: -l '1,2,0'. Either the '-l' or '-L' "\
        "arguments must be included.")
    argparser.add_argument(
        "-L", "--ligand_resname", dest="ligand_resname", 
        metavar="LIGAND_RESNAME", type=str, default="",
        help="The residue name of the ligand molecule for automatic index "\
        "selection. Either the '-l' or '-L' arguments must be included.")
    argparser.add_argument(
        "-c", "--cut_off_distance", dest="cut_off_distance", default=3.0,
        help="The cut-off distance (in Angstroms) used to define the QM2 "\
        "region of the ONIOM calculation. Default: 3.0.", type=float)
    argparser.add_argument(
        "-n", "--nprocs", dest="nprocs", default=1,
        help="The number of processors to use for ORCA calculations. "\
        "Default: 1.", type=int)
    argparser.add_argument(
        "-i", "--max_iterations", dest="max_iterations", default=2000,
        help="The maximum number of iterations for ORCA convergence. "\
        "Default: 2000.", type=int)
    argparser.add_argument(
        "-m", "--qm_method", dest="qm_method", default="B3LYP",
        help="The method to use for the QM region of the ONIOM calculation. "\
        "Please see the file orca_methods_basis_sets.pdf for all possible "\
        "options. Default: B3LYP.", type=str)
    argparser.add_argument(
        "-b", "--qm_basis_set", dest="qm_basis_set", default="6-311G",
        help="The basis set to use for the QM region of the ONIOM "\
        "calculation. Please see the file orca_methods_basis_sets.pdf for all "\
        "possible options. Default: 6-311G", type=str)
    argparser.add_argument(
        "-s", "--qm_charge_scheme", dest="qm_charge_scheme", default="CHELPG",
        help="The charge scheme to use for the QM region of the ONIOM "
        "calculation. Please see the file orca_methods_basis_sets.pdf for all "\
        "possible options. Default: CHELPG.", type=str)
    argparser.add_argument(
        "-q", "--qm_charge", dest="qm_charge", default=None,
        help="The total charge of the QM region of the ONIOM calculation. "\
        "If set to None, the quantity will be automatically computed. "\
        "Default: 0.", type=int)
    argparser.add_argument(
        "-u", "--qm_multiplicity", dest="qm_multiplicity", default=1,
        help="The multiplicity of the QM region of the ONIOM calculation. "\
        "Default: 1.", type=int)
    argparser.add_argument(
        "-M", "--qm2_method", dest="qm2_method", default="XTB",
        help="The method to use for the QM2 region of the ONIOM calculation. "\
        "Please see the file orca_methods_basis_sets.pdf for all possible "\
        "options. Default: XTB.", type=str)
    argparser.add_argument(
        "-S", "--qm2_charge_scheme", dest="qm2_charge_scheme", default="CHELPG",
        help="The charge scheme to use for the QM2 region of the ONIOM "
        "calculation. Please see the file orca_methods_basis_sets.pdf for all "\
        "possible options. Default: CHELPG.", type=str)
    argparser.add_argument(
        "-Q", "--qm2_charge", dest="qm2_charge", default=None,
        help="The total charge of the QM2 region of the ONIOM calculation. "\
        "If set to None, the quantity will be automatically computed. "\
        "Default: None.", type=int)
    argparser.add_argument(
        "-U", "--qm2_multiplicity", dest="qm2_multiplicity", default=1,
        help="The multiplicity of the QM2 region of the ONIOM calculation. "\
        "Default: 1.", type=int)
    argparser.add_argument(
        "-O", "--orca_path", dest="orca_path", default=None,
        help="An absolute path to the ORCA program. If not specified, ORCA "\
        "will be found from the shutil.which() command. Default: None.", 
        type=str)
    argparser.add_argument(
        "-x", "--skip_checks", dest="skip_checks", default=False, 
        help="By default, checks will be run at various stages of Qmrebind, "\
        "and if the checks fail, the calculation will not proceed. This "\
        "argument bypasses those checks and allows the calculation to "\
        "proceed anyways. Default: False.", action="store_true")
    
    args = argparser.parse_args()
    args = vars(args)
    anchor_list = args["anchor_list"]
    model_file = args["model_file"]
    ligand_indices = args["ligand_indices"]
    if ligand_indices != "":
        ligand_indices = base.initialize_indices(ligand_indices)
    else:
        ligand_indices = None
    ligand_resname = args["ligand_resname"]
    cut_off_distance = args["cut_off_distance"]
    nprocs = args["nprocs"]
    max_iterations = args["max_iterations"]
    qm_method = args["qm_method"]
    qm_basis_set = args["qm_basis_set"]
    qm_charge_scheme = args["qm_charge_scheme"]
    qm_charge = args["qm_charge"]
    qm_multiplicity = args["qm_multiplicity"]
    qm2_method = args["qm2_method"]
    qm2_charge_scheme = args["qm2_charge_scheme"]
    qm2_charge = args["qm2_charge"]
    qm2_mult = args["qm2_multiplicity"]
    orca_path = args["orca_path"]
    skip_checks = args["skip_checks"]
    
    model = seekr_base.load_model(model_file)
    
    run_qmrebind_seekr_anchors(
        anchor_list, model, ligand_indices, ligand_resname,
        cut_off_distance=cut_off_distance, nprocs=nprocs, 
        maxiter=max_iterations, qm_method=qm_method, qm_basis_set=qm_basis_set,
        qm_charge_scheme=qm_charge_scheme, qm_charge=qm_charge, 
        qm_mult=qm_multiplicity, qm2_method=qm2_method, 
        qm2_charge_scheme=qm2_charge_scheme, qm2_charge=qm2_charge, 
        qm2_mult=qm2_mult, orca_dir_pwd=orca_path, skip_checks=skip_checks)
    
    save_new_model(model)