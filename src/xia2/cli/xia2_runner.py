from __future__ import annotations

import argparse
import json
import os
import stat
import sys
from collections import defaultdict
from pathlib import Path


def generate_scripts(
    data_dirs,
    results_dir,
    options_file,
    conditions_file,
    n_jobs_to_generate,
    dials_source,
):

    protein_options_file = options_file
    protein_conditions_file = conditions_file

    # First, inspect the data and results directories to see what has not been processed
    ready = []
    processing_started = []
    processing_failed = []
    processing_finished = []
    name_to_data_dir = {}

    for parent in data_dirs:
        for dir in parent.iterdir():
            name = dir.name
            name_to_data_dir[name] = dir.parent
            # FIXME - update glob of filenames to h5 and test for if all complete
            n_files = len(list(dir.glob("*.h5")))
            if n_files >= 1:
                # print(f"{name} data fully transferred")
                if Path.is_dir(results_dir / name):
                    if Path.is_file(results_dir / name / "DataFiles" / "merged.mtz"):
                        processing_finished.append(str(name))
                    # FIXME - better test for processing to check if still running rather
                    # than failed part way through - query condor_q?
                    elif Path.is_dir(results_dir / name / "import"):
                        processing_started.append(str(name))
                    elif Path.is_dir(results_dir / name / "DataFiles"):
                        processing_failed.append(str(name))
                    else:
                        ready.append(str(name))
                else:
                    ready.append(str(name))
            elif n_files > 0:
                print(f"{name} data partially transferred")
            else:
                print(f"No data transferred yet for {name}")

    ready = sorted(ready)
    print(f"Ready for processing:\n{', '.join(ready)}")
    print(f"Started processing:\n{', '.join(processing_started)}")
    print(f"Finished processing:\n{', '.join(processing_finished)}")
    if processing_failed:
        print(
            f"Processing failed for some jobs: {', '.join(processing_failed)}"
            + f"\nPlease check results in {results_dir} for cause of failure, and"
            + f"\nthe options specified in {protein_options_file}"
            + "\nThen fix input options, delete the relevant processing directory for that job and try again"
            + "\n"
        )

    # for those ready, make a processing subdirectory and generate a script
    # based on information in the proteinoptions.json and proteinconditions.json

    # First extract the metadata defined in the protein files
    run_to_protein_condition = {}
    with open(protein_conditions_file) as f:
        conditions_dict = json.load(f)
    for prot, conditions in conditions_dict.items():
        for condition, runs in conditions.items():
            for run in runs:
                run_to_protein_condition[run] = [prot, condition]
    with open(protein_options_file) as f:
        options = json.load(f)

    # Now generate processing scripts for each data collection
    directory_to_script = {}
    n_generated = 0
    not_in_conditions = []
    for to_process in ready:

        if to_process not in run_to_protein_condition:
            not_in_conditions.append(to_process)
            continue
        if not Path.is_dir(results_dir / to_process):
            # print(f"Making directory {results_dir / to_process}")
            Path.mkdir(results_dir / to_process)
        image = os.fspath(
            next((name_to_data_dir[to_process] / to_process).glob("*.h5"))
        )
        prot, cond = run_to_protein_condition[to_process]
        script_options = " ".join(o for o in options[prot][cond])
        # Now write the script
        script = f"""#!/bin/bash
source {dials_source}
xia2.ssx image={image} {script_options}
"""
        script_name = "run_xia2.sh"
        script_file = results_dir / to_process / script_name
        with open(script_file, "w") as f:
            print(f"Writing xia2 script to {script_file}")
            f.write(script)
        os.chmod(os.fspath(script_file), stat.S_IRWXU)
        condor_script = """request_cpus = 20
request_memory = 20000M
executable = run_xia2.sh
log = process.log
error   = logs.err
output  = logs.out
queue

"""
        condor_script_name = "condor_submit.sh"
        condor_script_file = results_dir / to_process / condor_script_name
        with open(condor_script_file, "w") as f:
            print(f"Writing condor script to {condor_script_file}")
            f.write(condor_script)
        os.chmod(os.fspath(condor_script_file), stat.S_IRWXU)
        directory_to_script[results_dir / to_process] = condor_script_name
        n_generated += 1
        if n_generated >= n_jobs_to_generate:
            break

    if not_in_conditions:
        not_found = ", ".join(f for f in not_in_conditions)
        print(
            f"Warning, some image files not found in definitions in {protein_conditions_file}:\n"
            + f"{not_found}"
            + "\nDid you remember to add these to the definitions, or is there a typo?"
        )

    # This writes a submit script that can be used to submit %n_jobs_to_generate processing jobs for 'ready' datasets
    # FIXME - update submit.sh for PAL cluster
    if directory_to_script:
        submit_sh = ""
        for dir, script in directory_to_script.items():
            submit_sh += f"""cd {dir}\ncondor_submit {script}"""
        submit_sh += "\n"
        submit_script_name = "submit.sh"
        submit_script_file = results_dir / submit_script_name
        with open(submit_script_file, "w") as f:
            f.write(submit_sh)
        os.chmod(os.fspath(submit_script_file), stat.S_IRWXU)
        print(
            f"\nWritten submit script to process {n_generated} 'ready for processing' jobs to {submit_script_file}"
        )
        print(f"To run the script, just run {submit_script_file}\n")

    # check if submitted same script - perhaps make a separate resubmit_script.sh to redo those with
    # up to date options e.g.

    # This next section generates scripts to merge multiple integrated runs from a given protein/condition.
    prot_condition_to_processed = defaultdict(list)
    for dir in processing_finished:
        prot, cond = run_to_protein_condition[dir]
        prot_condition_to_processed[prot + "_" + cond].append(dir)
    if prot_condition_to_processed:
        for name, merge_group in prot_condition_to_processed.items():
            if len(merge_group) > 1:
                merge_dir_name = "merge_" + name
                if not Path.is_dir(results_dir / merge_dir_name):
                    print(f"Making directory {results_dir / merge_dir_name}")
                    Path.mkdir(results_dir / merge_dir_name)
                # need to get the right reference
                prot, cond = name.split("_")
                reference = None
                for val in options[prot][cond]:
                    if "reference=" in val or "model=" in val:
                        reference = val
                        break
                if not reference:
                    print(
                        f"Can't find reference pdb model in {protein_options_file} for {name}"
                    )
                    continue
                submit_sh = f"""#!/bin/bash
source {dials_source}
xia2.ssx_reduce {reference}"""
                for dataset in merge_group:
                    submit_sh += f" {results_dir / dataset / 'batch_*' / 'integrated*'}"
                submit_sh += "\n"
                submit_script_name = results_dir / merge_dir_name / "merge_all.sh"
                with open(submit_script_name, "w") as f:
                    f.write(submit_sh)
                os.chmod(os.fspath(submit_script_name), stat.S_IRWXU)

                condor_script = """request_cpus = 20
request_memory = 20000M
executable = merge_all.sh
log = process.log
error   = logs.err
output  = logs.out
queue

"""
                condor_script_name = "condor_submit.sh"
                condor_script_file = results_dir / merge_dir_name / condor_script_name
                with open(condor_script_file, "w") as f:
                    f.write(condor_script)
                os.chmod(os.fspath(condor_script_file), stat.S_IRWXU)
                print(f"\nWritten merging script for {name} to {condor_script_file}")
                print(
                    "Please run the script yourself in the directory to merge the data for this protein and condition\n"
                )


def _parse_image_range(name_string: str):
    print(f"interpreting image-name {name_string} as a range separated by ':'")
    try:
        first, last = name_string.split(":")
    except ValueError:
        print("image name option must be in form first:last")
        sys.exit(0)
    else:
        # now strip leading zeroes
        ndigits = len(first)
        assert (
            len(last) == ndigits
        ), "First and last image number don't have same number of leading zeroes"
        last = last.lstrip("0")
        first = first.lstrip("0")
        try:
            last_int = int(last)
            first_int = int(first)
        except ValueError:
            print("First and last image names must be integers")
            sys.exit(0)
        else:
            images = [
                f"{str(i).zfill(ndigits)}" for i in range(first_int, last_int + 1)
            ]
    return images


def run():

    from xia2.cli.runner_options import options as config

    # These options define where files are
    # FIXME - Update data and results directories
    data_dirs = [Path(i) for i in config["data_dirs"]]
    results_dir = Path(config["results_dir"])
    dials_source = config["dials_source"]
    protein_options_file = results_dir / "xia2_options.json"
    protein_conditions_file = results_dir / "conditions.json"

    parser = argparse.ArgumentParser(
        description="Add new runs to options files, generate runner scripts"
    )
    parser.add_argument(
        "--show", action="store_true", help="Print a table of the data stored"
    )
    parser.add_argument(
        "--add",
        type=str,
        nargs=3,
        help="Add a run, or series of runs. e.g --add cld apo 00001 or --add cld apo 00001:00010",
        metavar=("protein", "condition", "image-name"),
    )
    parser.add_argument(
        "--remove",
        type=str,
        nargs=3,
        help="Remove a run, or series of runs.",
        metavar=("protein", "condition", "image-name"),
    )
    parser.add_argument(
        "--generate-scripts",
        metavar=("njobs"),
        nargs="?",
        type=int,
        const=5,
        default=0,
        help="Inspect raw data files, make processing directories \nand generate scripts for $njobs jobs.",
    )
    parser.add_argument(
        "--remove-condition",
        type=str,
        nargs=2,
        help="Remove all data for a protein condition",
        metavar=("protein", "condition"),
    )
    parser.add_argument(
        "--remove-protein",
        type=str,
        help="Remove all data for a protein",
        metavar="protein",
    )

    args = parser.parse_args()

    if args.add:
        with open(protein_conditions_file, "r") as f:
            conditions_dict = json.load(f)
        prot, cond, name = args.add
        if ":" in name:
            images = _parse_image_range(name)
        else:
            images = [name]

        if prot not in conditions_dict:
            conditions_dict[prot] = {}
        if cond not in conditions_dict[prot]:
            conditions_dict[prot][cond] = []
        for img in images:
            conditions_dict[prot][cond].append(img)
        conditions_dict[prot][cond] = list(set(conditions_dict[prot][cond]))
        with open(protein_conditions_file, "w") as f:
            json.dump(conditions_dict, f)
        names = ", ".join(f"{img}" for img in images)
        print(f"Added {names} to protein={prot}, condition={cond}")
        print("Run xia2_runner --show to see full conditions dictionary")

    if args.remove:
        with open(protein_conditions_file, "r") as f:
            conditions_dict = json.load(f)
        prot, cond, name = args.remove
        if prot not in conditions_dict:
            print(f"Unable to find protein {prot} in existing dictionary")
            return
        if cond not in conditions_dict[prot]:
            print(
                f"Unable to find condition {cond} for protein {prot} in existing dictionary"
            )
            return
        if ":" in name:
            images = _parse_image_range(name)
        else:
            images = [name]

        removed = []
        for img in images:
            try:
                index = conditions_dict[prot][cond].index(img)
            except ValueError:
                print(
                    f"Unable to find image {img} in dictionary for protein={prot}, condition={cond}"
                )
            else:
                del conditions_dict[prot][cond][index]
                removed.append(img)
        removed = ", ".join(f"{img}" for img in removed)

        with open(protein_conditions_file, "w") as f:
            json.dump(conditions_dict, f)
        print(f"Removed {removed} from dictionary for protein={prot}, condition={cond}")
        print("Run xia2_runner --show to see full conditions dictionary")

    if args.show:
        with open(protein_conditions_file, "r") as f:
            conditions_dict = json.load(f)
        from dials.util import tabulate

        headers = ["protein", "condition", "collections"]
        rows = []
        for prot, v in conditions_dict.items():
            for condition, namelist in v.items():
                names = ""
                length = 0
                for name in namelist:
                    if length > 50:
                        names += f"\n{name}, "
                        length = len(name)
                    else:
                        names += f"{name}, "
                        length += len(name)
                row = [f"{prot}", f"{condition}", names]
                rows.append(row)
        print(tabulate(rows, headers))

    if args.remove_condition:
        with open(protein_conditions_file, "r") as f:
            conditions_dict = json.load(f)
        prot, cond = args.remove_condition
        if prot not in conditions_dict:
            print(f"Unable to find protein {prot} in existing dictionary")
            return
        if cond not in conditions_dict[prot]:
            print(
                f"Unable to find condition {cond} for protein {prot} in existing dictionary"
            )
            return
        del conditions_dict[prot][cond]
        with open(protein_conditions_file, "w") as f:
            json.dump(conditions_dict, f)
        print(f"Removed {cond} from dictionary for protein={prot}")
        print("Run xia2_runner --show to see full conditions dictionary")

    if args.remove_protein:
        with open(protein_conditions_file, "r") as f:
            conditions_dict = json.load(f)
        prot = args.remove_protein
        if prot not in conditions_dict:
            print(f"Unable to find protein {prot} in existing dictionary")
            return
        del conditions_dict[prot]
        with open(protein_conditions_file, "w") as f:
            json.dump(conditions_dict, f)
        print(f"Removed protein={prot} from dictionary")
        print("Run xia2_runner --show to see full conditions dictionary")

    if args.generate_scripts:
        generate_scripts(
            data_dirs,
            results_dir,
            protein_options_file,
            protein_conditions_file,
            n_jobs_to_generate=args.generate_scripts,
            dials_source=dials_source,
        )


if __name__ == "__main__":
    run()
