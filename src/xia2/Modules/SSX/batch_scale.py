from __future__ import annotations

import copy
import logging
import os
from pathlib import Path

import numpy as np

from cctbx import sgtbx
from dials.algorithms.scaling.algorithm import ScalingAlgorithm
from dials.algorithms.scaling.scaling_library import determine_best_unit_cell
from dials.algorithms.symmetry.reindex_to_reference import (
    determine_reindex_operator_against_reference,
)
from dials.array_family import flex
from dials.util.reference import intensities_from_reference_file
from dxtbx.model import ExperimentList

logger = logging.getLogger("dials")

from typing import List

from dials.util.options import ArgumentParser
from dxtbx.serialize import load
from libtbx import phil

from xia2.Modules.SSX.data_reduction_definitions import FilePair


class BatchScale(object):
    def __init__(
        self, batches: List[FilePair], reference: Path, space_group: sgtbx.space_group
    ):
        phil_scope = phil.parse(
            """
    include scope dials.command_line.scale.phil_scope
""",
            process_includes=True,
        )
        parser = ArgumentParser(phil=phil_scope, check_format=False)
        params, _ = parser.parse_args(args=[], quick_parse=True)
        params.model = "KB"
        params.scaling_options.full_matrix = False
        params.weighting.error_model.error_model = None
        params.scaling_options.outlier_rejection = "simple"
        params.reflection_selection.min_partiality = 0.25
        params.reflection_selection.method = "intensity_ranges"
        params.reflection_selection.Isigma_range = [2, 0]
        params.cut_data.partiality_cutoff = 0.25
        self.params = params
        self.input_sg = space_group
        self.reference = reference
        self.input_batches = batches
        self._experiments = ExperimentList([])
        self._reflections = []
        self._output_expt_files: List[str] = []
        self._output_refl_files: List[str] = []

        # create a single table and expt per batch
        all_expts = ExperimentList([])

        class SANoStats(ScalingAlgorithm):
            def calculate_merging_stats(self):
                pass

        # we modify the input, but these don't get saved
        for i, fp in enumerate(self.input_batches):
            expts = load.experiment_list(fp.expt, check_format=False)
            table = flex.reflection_table.from_file(fp.refl)
            sel = (
                table["intensity.sum.value"] / (table["intensity.sum.variance"] ** 0.5)
            ) > 2.0
            n_sel = sel.count(True)
            if (n_sel / sel.size()) < 0.05 and n_sel > 10:
                params.reflection_selection.Isigma_range = [0, 0]
            else:
                table = table.select(
                    (
                        table["intensity.sum.value"]
                        / (table["intensity.sum.variance"] ** 0.5)
                    )
                    > 2.0
                )
            all_expts.extend(expts)
            expt = copy.deepcopy(expts[0])
            expt.scaling_model = None
            expt.identifier = str(i)
            self._experiments.append(expt)
            if "inverse_scale_factor" in table:
                table["intensity.sum.value"] /= table["inverse_scale_factor"]
                table["intensity.sum.variance"] /= (table["inverse_scale_factor"]) ** 2
                del table["inverse_scale_factor"]
                del table["inverse_scale_factor_variance"]
            for k in list(table.experiment_identifiers().keys()):
                del table.experiment_identifiers()[k]
            table["id"] = flex.int(table.size(), i)
            table.experiment_identifiers()[i] = str(i)
            self._reflections.append(table)

        wavelength = np.mean([expt.beam.get_wavelength() for expt in expts])
        best_unit_cell = determine_best_unit_cell(all_expts)
        for expt in self._experiments:
            expt.crystal.set_unit_cell(best_unit_cell)
            expt.beam.set_wavelength(wavelength)

        self.algorithm = SANoStats(params, self._experiments, self._reflections)

    def run(self):
        self.algorithm.run()
        del self._reflections
        if self.reference:
            wavelength = np.mean(
                [expt.beam.get_wavelength() for expt in self._experiments]
            )
            reference_miller_set = intensities_from_reference_file(
                os.fspath(self.reference), wavelength=wavelength
            )
            test_miller_set = self.algorithm.scaled_miller_array
            change_of_basis_op = determine_reindex_operator_against_reference(
                test_miller_set, reference_miller_set
            )
            for i, fp in enumerate(self.input_batches):
                expts = load.experiment_list(fp.expt, check_format=False)
                refl = flex.reflection_table.from_file(fp.refl)
                for expt in expts:
                    expt.crystal = expt.crystal.change_basis(change_of_basis_op)
                    expt.crystal.set_space_group(self.input_sg)
                expts.as_file(f"processed_{i}.expt")

                refl["miller_index"] = change_of_basis_op.apply(refl["miller_index"])
                refl.as_file(f"processed_{i}.refl")
                self._output_expt_files.append(f"processed_{i}.expt")
                self._output_refl_files.append(f"processed_{i}.refl")
