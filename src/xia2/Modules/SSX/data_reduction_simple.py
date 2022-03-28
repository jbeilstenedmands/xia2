from __future__ import annotations

import concurrent.futures
import copy
import json
import logging
import os
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional, Tuple

import iotbx.phil
from cctbx import crystal, sgtbx, uctbx
from dials.algorithms.scaling.scaling_library import determine_best_unit_cell
from dials.array_family import flex
from dxtbx.model import ExperimentList
from dxtbx.serialize import load

from xia2.Driver.timing import record_step
from xia2.Handlers.Streams import banner
from xia2.Modules.SSX.data_reduction_base import BaseDataReduction
from xia2.Modules.SSX.data_reduction_programs import (
    FilePair,
    cluster_all_unit_cells,
    inspect_directories,
    merge,
    reference_reindex,
    run_uc_cluster,
    scale,
    scale_against_model,
    scale_cosym,
    select_crystals_close_to,
    split_filtered_data,
)
from xia2.Modules.SSX.reporting import statistics_output_from_scaled_files

xia2_logger = logging.getLogger(__name__)

FilesDict = Dict[
    int, FilePair
]  # A Dict where the keys are an index, corresponding to a filepair


@dataclass
class SimpleReductionParams:
    space_group: sgtbx.space_group
    batch_size: int = 1000
    nproc: int = 1
    d_min: Optional[float] = None
    anomalous: bool = False
    cluster_threshold: float = 1000.0
    absolute_angle_tolerance: float = 0.5
    absolute_length_tolerance: float = 0.2
    model: Optional[Path] = None

    @classmethod
    def from_phil(cls, params: iotbx.phil.scope_extract):
        """Construct from xia2.cli.ssx phil_scope."""
        model = None
        if params.scaling.model:
            model = Path(params.scaling.model).resolve()
        return cls(
            params.space_group,
            params.batch_size,
            params.nproc,
            params.d_min,
            params.scaling.anomalous,
            params.clustering.threshold,
            params.clustering.absolute_angle_tolerance,
            params.clustering.absolute_length_tolerance,
            model,
        )


def assess_for_indexing_ambiguities(
    space_group: sgtbx.space_group_info, unit_cell: uctbx.unit_cell
):
    # if lattice symmetry higher than space group symmetry, then need to
    # assess for indexing ambiguity.
    cs = crystal.symmetry(unit_cell=unit_cell, space_group=sgtbx.space_group())
    # Get cell reduction operator
    cb_op_inp_minimum = cs.change_of_basis_op_to_minimum_cell()
    # New symmetry object with changed basis
    minimum_symmetry = cs.change_basis(cb_op_inp_minimum)

    # Get highest symmetry compatible with lattice
    lattice_group = sgtbx.lattice_symmetry_group(
        minimum_symmetry.unit_cell(),
        max_delta=5,
        enforce_max_delta_for_generated_two_folds=True,
    )
    if lattice_group.order_z() > space_group.group().order_z():
        return True
    else:
        return False


def load_crystal_data_from_new_expts(new_data):
    data = {}
    for file_pair in new_data.values():
        new_expts = load.experiment_list(file_pair.expt, check_format=False)
        data[str(file_pair.expt)] = {
            "identifiers": copy.deepcopy(new_expts.identifiers()),
            "crystals": copy.deepcopy(new_expts.crystals()),
            "keep_all_original": True,
        }
    return data


def determine_best_unit_cell_crystals(crystals_data):
    """Set the median unit cell as the best cell, for consistent d-values across
    experiments."""
    uc_params = [flex.double() for i in range(6)]
    for v in crystals_data.values():
        for c in v["crystals"]:
            unit_cell = c.get_recalculated_unit_cell() or c.get_unit_cell()
            for i, p in enumerate(unit_cell.parameters()):
                uc_params[i].append(p)
    best_unit_cell = uctbx.unit_cell(parameters=[flex.median(p) for p in uc_params])
    return best_unit_cell


class SimpleDataReduction(BaseDataReduction):
    def run(self, reduction_params: SimpleReductionParams) -> None:
        """
        A simple workflow for data reduction. First filter the input data, either
        by clustering on unit cells, or comparing against a previous cell. Then
        reindex data in batches, using cosym, followed by scaling and merging.
        """

        data_already_reindexed: FilesDict = {}
        data_to_reindex: FilesDict = {}

        data_reduction_wd = self._main_directory / "data_reduction"
        filter_wd = data_reduction_wd / "prefilter"
        reindex_wd = data_reduction_wd / "reindex"
        scale_wd = data_reduction_wd / "scale"

        # see if we have any data that has been marked as successfully reindexed.
        reidx_results = reindex_wd / "reindexing_results.json"
        if reidx_results.is_file():
            previous = json.load(reidx_results.open())
            data_already_reindexed = {
                int(i): FilePair(Path(v["expt"]), Path(v["refl"]))
                for i, v in previous["reindexed_files"].items()
            }
            for file_pair in data_already_reindexed.values():
                assert file_pair.check()

        new_data = {}
        if self._new_directories_to_process:
            new_data = inspect_directories(self._new_directories_to_process)

        best_unit_cell = None
        filter_results = filter_wd / "filter_results.json"
        if filter_results.is_file():
            result = json.load(filter_results.open())
            best_unit_cell = uctbx.unit_cell(result["best_unit_cell"])

        new_best_unit_cell = None
        if new_data:
            xia2_logger.notice(banner("Filtering"))  # type: ignore
            crystals_data = load_crystal_data_from_new_expts(new_data)
            if not best_unit_cell:  # i.e. no previous data reduction
                new_best_unit_cell = None
                if reduction_params.cluster_threshold:
                    good_crystals_data = run_uc_cluster(
                        filter_wd,
                        crystals_data,
                        reduction_params.cluster_threshold,
                    )
                    n_cryst = sum(
                        len(v["identifiers"]) for v in good_crystals_data.values()
                    )
                    """cluster_expts, cluster_refls = self.filter_on_unit_cell_clustering(
                        filter_wd,
                        new_data,
                        reduction_params.cluster_threshold,
                    )"""
                elif (
                    reduction_params.absolute_angle_tolerance
                    and reduction_params.absolute_length_tolerance
                ):
                    # calculate the median unit cell
                    """cluster_expts = ExperimentList([])
                    for file_pair in new_data.values():
                        cluster_expts.extend(
                            load.experiment_list(file_pair.expt, check_format=False)
                        )"""
                    new_best_unit_cell = determine_best_unit_cell_crystals(
                        crystals_data
                    )
                    good_crystals_data = select_crystals_close_to(
                        crystals_data,
                        new_best_unit_cell,
                        reduction_params.absolute_angle_tolerance,
                        reduction_params.absolute_length_tolerance,
                    )
                    n_cryst = sum(
                        len(v["identifiers"]) for v in good_crystals_data.values()
                    )
                else:  # join all data for splitting
                    good_crystals_data = crystals_data
                    n_cryst = sum(len(v["identifiers"]) for v in crystals_data.values())
                    xia2_logger.info(
                        f"No unit cell filtering applied, {n_cryst} crystals loaded."
                    )
                if not new_best_unit_cell:
                    new_best_unit_cell = determine_best_unit_cell_crystals(
                        crystals_data
                    )
                data_to_reindex = split_filtered_data(
                    filter_wd, new_data, good_crystals_data, reduction_params.batch_size
                )
                # now write out the best cell
                all_ucs = []
                for v in good_crystals_data.values():
                    all_ucs.extend(
                        [c.get_unit_cell().parameters() for c in v["crystals"]]
                    )
                result = {
                    "best_unit_cell": [
                        round(i, 4) for i in new_best_unit_cell.parameters()
                    ],
                    "n_cryst": n_cryst,
                    "unit_cells": all_ucs,
                }
                with open(filter_wd / "filter_results.json", "w") as fp:
                    json.dump(result, fp)

            else:  # if previous data reduction
                xia2_logger.info(
                    f"Using cell parameters from previous clustering analysis"
                )
                if (
                    reduction_params.absolute_angle_tolerance
                    and reduction_params.absolute_length_tolerance
                ):

                    good_crystals_data = select_crystals_close_to(
                        crystals_data,
                        best_unit_cell,
                        reduction_params.absolute_angle_tolerance,
                        reduction_params.absolute_length_tolerance,
                    )
                    n_cryst = sum(
                        len(v["identifiers"]) for v in good_crystals_data.values()
                    )
                else:
                    good_crystals_data = crystals_data
                    n_cryst = sum(len(v["identifiers"]) for v in crystals_data.values())
                    xia2_logger.info(
                        f"No unit cell filtering applied, {n_cryst} crystals loaded."
                    )

                while n_cryst < reduction_params.batch_size:
                    # we want all jobs to use at least batch_size crystals, so join on
                    # to last reindexed batch until this has been satisfied
                    already_reidx_keys = list(data_already_reindexed.keys())
                    if not already_reidx_keys:
                        break
                    last_batch = data_already_reindexed.pop(already_reidx_keys)[-1]
                    last_expts = load.experiment_list(
                        last_batch.expt, check_format=False
                    )
                    good_crystals_data[str(last_batch.expt)] = {
                        "identifiers": copy.deepcopy(last_expts.identifiers()),
                        "crystals": copy.deepcopy(last_expts.crystals()),
                        "keep_all_original": True,
                    }
                    new_data_last_key = list(new_data.keys)[-1]
                    new_data[new_data_last_key + 1] = last_batch
                    n_cryst += len(last_expts)

                n_already_reindexed_files = len(list(data_already_reindexed.keys()))
                data_to_reindex = split_filtered_data(
                    filter_wd,
                    new_data,
                    good_crystals_data,
                    reduction_params.batch_size,
                    offset=n_already_reindexed_files,
                )

                all_crystals_data = good_crystals_data
                for file_pair in data_already_reindexed.values():
                    expts = load.experiment_list(file_pair.expt, check_format=False)
                    all_crystals_data[file_pair.expt] = {
                        "crystals": copy.deepcopy(expts.crystals())
                    }
                    n_cryst += len(expts)

                new_best_unit_cell = determine_best_unit_cell_crystals(
                    all_crystals_data
                )
                all_ucs = []
                for v in all_crystals_data.values():
                    all_ucs.extend(
                        [c.get_unit_cell().parameters() for c in v["crystals"]]
                    )
                result = {
                    "best_unit_cell": [
                        round(i, 4) for i in new_best_unit_cell.parameters()
                    ],
                    "n_cryst": n_cryst,
                    "unit_cells": all_ucs,
                }
                with open(filter_wd / "filter_results.json", "w") as fp:
                    json.dump(result, fp)
        else:
            new_best_unit_cell = best_unit_cell
            if not new_best_unit_cell:
                raise ValueError(
                    "Error in intepreting new and previous processing results"
                )

        # the assumption in the simple data reducer is that we know the space group
        # and have a good value for the unit cell.
        sym_requires_reindex = assess_for_indexing_ambiguities(
            reduction_params.space_group, new_best_unit_cell
        )
        if sym_requires_reindex:
            files_to_scale = self.reindex(
                reindex_wd,
                data_to_reindex,
                data_already_reindexed,
                space_group=reduction_params.space_group,
                nproc=reduction_params.nproc,
                d_min=reduction_params.d_min,
            )
        else:
            # save new data as already reindexed
            files_to_scale = {**data_already_reindexed, **data_to_reindex}
            output_files_to_scale = {
                k: {"expt": str(v.expt), "refl": str(v.refl)}
                for k, v in files_to_scale.items()
            }
            if not Path.is_dir(reindex_wd):
                Path.mkdir(reindex_wd)
            reidx_results = reindex_wd / "reindexing_results.json"
            with open(reidx_results, "w") as f:
                data = {"reindexed_files": output_files_to_scale}
                json.dump(data, f)

        # if we get here, we have successfully prepared the new data for scaling.
        # So save this to allow reloading in future for iterative workflows.
        data_reduction_progress = {
            "directories_processed": [
                str(i)
                for i in (
                    self._directories_previously_processed
                    + self._new_directories_to_process
                )
            ]
        }
        with open(data_reduction_wd / "data_reduction.json", "w") as fp:
            json.dump(data_reduction_progress, fp)

        if reduction_params.model:
            self.scale_and_merge_using_model(
                scale_wd,
                files_to_scale,
                anomalous=reduction_params.anomalous,
                d_min=reduction_params.d_min,
                model=reduction_params.model,
                nproc=reduction_params.nproc,
                best_unit_cell=new_best_unit_cell,
            )
        else:
            self.scale_and_merge(
                scale_wd,
                files_to_scale,
                anomalous=reduction_params.anomalous,
                d_min=reduction_params.d_min,
                best_unit_cell=new_best_unit_cell,
            )

    @staticmethod
    def filter_on_unit_cell_clustering(
        working_directory: Path,
        new_data: FilesDict,
        threshold: float,
    ) -> Tuple[ExperimentList, flex.reflection_table]:

        """
        Filter the integrated data using dials.cluster_unit_cell. Takes the
        largest cluster found and splits it into batches for reindexing.
        """
        xia2_logger.info("Clustering with dials.cluster_unit_cell")  # type: ignore
        main_cluster_files = cluster_all_unit_cells(
            working_directory,
            new_data,
            threshold,
        )

        # save the results to a json
        cluster_expts = load.experiment_list(
            main_cluster_files.expt,
            check_format=False,
        )
        cluster_refls = flex.reflection_table.from_file(main_cluster_files.refl)
        n_in_cluster = len(cluster_expts)
        uc = determine_best_unit_cell(cluster_expts)
        uc_string = ",".join(f"{i:.2f}" for i in uc.parameters())
        xia2_logger.info(
            f"""Largest cluster contains {n_in_cluster} crystals with cell parameters:
  {uc_string}"""
        )
        result = {
            "best_unit_cell": [round(i, 4) for i in uc.parameters()],
            "n_in_cluster": n_in_cluster,
            "unit_cells": [
                e.crystal.get_unit_cell().parameters() for e in cluster_expts
            ],
        }
        with open(working_directory / "cluster_results.json", "w") as fp:
            json.dump(result, fp)

        return cluster_expts, cluster_refls

    @staticmethod
    def reindex(
        working_directory: Path,
        data_to_reindex: FilesDict,
        data_already_reindexed: FilesDict,
        space_group: sgtbx.space_group,
        nproc: int = 1,
        d_min: float = None,
    ) -> FilesDict:
        """
        Runs dials.scale + dials.cosym on each batch to resolve indexing
        ambiguities. If there is more than one batch, the first batch is used
        as a reference, and dials.reindex is used to reindex the other batches
        against this reference.
        """
        sys.stdout = open(os.devnull, "w")  # block printing from cosym

        if not Path.is_dir(working_directory):
            Path.mkdir(working_directory)

        if 0 in data_already_reindexed:
            reference_files = data_already_reindexed[0]
        else:
            assert 0 in data_to_reindex  # make sure we have something that will
            # become the reference file.
        files_for_reference_reindex: FilesDict = {}
        reindexed_results: FilesDict = {}
        xia2_logger.notice(banner("Reindexing"))  # type: ignore
        with record_step(
            "dials.scale/dials.cosym (parallel)"
        ), concurrent.futures.ProcessPoolExecutor(max_workers=nproc) as pool:
            cosym_futures: Dict[Any, int] = {
                pool.submit(
                    scale_cosym,
                    working_directory,
                    files,
                    index,
                    space_group,
                    d_min,
                ): index
                for index, files in data_to_reindex.items()
            }
            for future in concurrent.futures.as_completed(cosym_futures):
                try:
                    result = future.result()
                except Exception as e:
                    raise ValueError(
                        f"Unsuccessful scaling and symmetry analysis of the new data. Error:\n{e}"
                    )
                else:
                    if list(result.keys()) == [0]:
                        reference_files = result[0]
                        reindexed_results[0] = result[0]
                    else:
                        files_for_reference_reindex.update(result)

        # now do reference reindexing
        with record_step(
            "dials.reindex (parallel)"
        ), concurrent.futures.ProcessPoolExecutor(max_workers=nproc) as pool:
            reidx_futures: Dict[Any, int] = {
                pool.submit(
                    reference_reindex, working_directory, reference_files, files
                ): index
                for index, files in files_for_reference_reindex.items()
            }
            for future in concurrent.futures.as_completed(reidx_futures):
                try:
                    result = future.result()
                    i = reidx_futures[future]
                except Exception as e:
                    raise ValueError(
                        f"Unsuccessful reindexing of the new data. Error:\n{e}"
                    )
                else:
                    reindexed_results[i] = result
                    xia2_logger.info(
                        f"Reindexed batch {i+1} using batch 1 as reference"
                    )

        files_to_scale = {**data_already_reindexed, **reindexed_results}
        output_files_to_scale = {
            k: {"expt": str(v.expt), "refl": str(v.refl)}
            for k, v in files_to_scale.items()
        }

        reidx_results = working_directory / "reindexing_results.json"
        with open(reidx_results, "w") as f:
            data = {"reindexed_files": output_files_to_scale}
            json.dump(data, f)
        sys.stdout = sys.__stdout__  # restore printing
        return files_to_scale

    @staticmethod
    def scale_and_merge(
        working_directory: Path,
        files_to_scale: FilesDict,
        anomalous: bool = True,
        d_min: float = None,
        best_unit_cell: Optional[uctbx.unit_cell] = None,
    ) -> None:
        """Run scaling and merging"""

        if not Path.is_dir(working_directory):
            Path.mkdir(working_directory)
        xia2_logger.notice(banner("Scaling & Merging"))  # type: ignore
        scaled_expts, scaled_table = scale(
            working_directory,
            files_to_scale,
            anomalous,
            d_min,
            best_unit_cell,
        )
        merge(working_directory, scaled_expts, scaled_table, d_min)

    @staticmethod
    def scale_and_merge_using_model(
        working_directory: Path,
        files_to_scale: FilesDict,
        anomalous: bool = True,
        d_min: float = None,
        model: Path = None,
        nproc: int = 1,
        best_unit_cell: Optional[uctbx.unit_cell] = None,
    ) -> None:
        """Run scaling and merging"""

        if not Path.is_dir(working_directory):
            Path.mkdir(working_directory)
        xia2_logger.notice(banner("Scaling using model"))  # type: ignore

        scaled_results: FilesDict = {}
        with record_step(
            "dials.scale (parallel)"
        ), concurrent.futures.ProcessPoolExecutor(max_workers=nproc) as pool:
            scale_futures: Dict[Any, int] = {
                pool.submit(
                    scale_against_model,
                    working_directory,
                    files,
                    index,
                    model,
                    anomalous,
                    d_min,
                    best_unit_cell,
                ): index
                for index, files in files_to_scale.items()
            }
            for future in concurrent.futures.as_completed(scale_futures):
                try:
                    result = future.result()
                    i = scale_futures[future]
                except Exception as e:
                    xia2_logger.warning(f"Unsuccessful scaling of group. Error:\n{e}")
                else:
                    xia2_logger.info(f"Completed scaling of group {i+1}")
                    scaled_results.update(result)

        xia2_logger.notice(banner("Merging"))  # type: ignore
        scaled_expts = ExperimentList([])
        scaled_tables = []
        for file_pair in scaled_results.values():
            scaled_expts.extend(
                load.experiment_list(file_pair.expt, check_format=False)
            )
            scaled_tables.append(flex.reflection_table.from_file(file_pair.refl))
        scaled_table = flex.reflection_table.concat(scaled_tables)

        n_final = len(scaled_expts)
        uc = determine_best_unit_cell(scaled_expts)
        uc_str = ", ".join(str(round(i, 3)) for i in uc.parameters())
        xia2_logger.info(
            f"{n_final} crystals scaled in space group {scaled_expts[0].crystal.get_space_group().info()}\nMedian cell: {uc_str}"
        )
        xia2_logger.info(
            statistics_output_from_scaled_files(scaled_expts, scaled_table, uc)
        )

        merge(working_directory, scaled_expts, scaled_table, d_min)
