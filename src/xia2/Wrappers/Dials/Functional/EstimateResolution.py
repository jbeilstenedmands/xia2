from __future__ import annotations

from pathlib import Path

import libtbx.phil
from dials.array_family import flex
from dials.util.resolution_analysis import Resolutionizer, phil_defaults
from dxtbx.model import ExperimentList

from xia2.Driver.timing import record_step
from xia2.lib.bits import _get_number
from xia2.Modules.SSX.util import log_to_file, run_in_directory
from xia2.Wrappers.Dials.Functional import diff_phil_from_params_and_scope


class EstimateResolution:
    def __init__(self, working_directory: Path | None = None):
        ## Working directory is where any output (logfiles, datafiles) will be saved
        if working_directory:
            self._working_directory = working_directory
        else:
            self._working_directory = Path.cwd()

        self.params: libtbx.phil.scope_extract = phil_defaults.extract()
        ## Set any defaults
        self.params.resolution.nbins = 100

        ## Define outputs which are not part of params
        self._resolution_rmerge: float | None = None
        self._resolution_completeness: float | None = None
        self._resolution_cc_half_significance: float | None = None
        self._resolution_cc_half: float | None = None
        self._resolution_isigma: float | None = None
        self._resolution_misigma: float | None = None

    def set_nbins(self, nbins: int) -> None:
        self.params.resolution.nbins = nbins

    def set_limit_rmerge(self, rmerge: float) -> None:
        self.params.resolution.rmerge = rmerge

    def set_limit_completeness(self, completeness: float) -> None:
        self.params.resolution.completeness = completeness

    def set_limit_cc_half(self, cc_half: float) -> None:
        self.params.resolution.cc_half = cc_half

    def set_cc_half_fit(self, cc_half_fit: str) -> None:
        self.params.resolution.cc_half_fit = cc_half_fit

    def set_cc_half_significance_level(self, cc_half_significance_level: float) -> None:
        self.params.resolution.cc_half_significance_level = cc_half_significance_level

    def set_limit_isigma(self, isigma: float) -> None:
        self.params.resolution.isigma = isigma

    def set_limit_misigma(self, misigma: float) -> None:
        self.params.resolution.misigma = misigma

    def set_labels(self, labels: list[str]) -> None:
        self.params.resolution.labels = labels

    def set_batch_range(self, start: int, end: int) -> None:
        self.params.resolution.batch_range = (start, end)

    def _record_results(self, logfile: str):
        with open(logfile, "r") as log:
            for record in log:
                if "Resolution rmerge" in record:
                    self._resolution_rmerge = float(record.split()[-1])
                if "Resolution completeness" in record:
                    self._resolution_completeness = float(record.split()[-1])
                if "Resolution cc_half_significance_level" in record:
                    self._resolution_cc_half_significance = float(record.split()[-1])
                elif "Resolution cc_half" in record:
                    self._resolution_cc_half = float(record.split()[-1])
                if "Resolution I/sig" in record:
                    self._resolution_isigma = float(record.split()[-1])
                if "Resolution Mn(I/sig)" in record:
                    self._resolution_misigma = float(record.split()[-1])

    def run_on_mtz(self, scaled_unmerged_mtz: Path) -> None:
        xpid = _get_number()
        logfile = f"{xpid}_dials.estimate_resolution.log"
        with (
            run_in_directory(self._working_directory),
            record_step("dials.estimate_resolution"),
        ):
            with log_to_file(logfile) as dials_logger:
                diff_phil = diff_phil_from_params_and_scope(self.params, phil_defaults)
                dials_logger.info(diff_phil)
                m = Resolutionizer.from_unmerged_mtz(
                    scaled_unmerged_mtz, self.params.resolution
                )
                m.resolution_auto()
            self._record_results(logfile)

    def run(
        self, experiments: ExperimentList, reflections: flex.reflection_table
    ) -> None:
        xpid = _get_number()
        logfile = f"{xpid}_dials.estimate_resolution.log"
        with (
            run_in_directory(self._working_directory),
            record_step("dials.estimate_resolution"),
        ):
            with log_to_file(logfile) as dials_logger:
                diff_phil = diff_phil_from_params_and_scope(self.params, phil_defaults)
                dials_logger.info(diff_phil)
                m = Resolutionizer.from_reflections_and_experiments(
                    [reflections], experiments, self.params.resolution
                )
                m.resolution_auto()
            self._record_results(logfile)

    @property
    def resolution_completeness(self) -> float | None:
        return self._resolution_completeness

    @property
    def resolution_cc_half(self) -> float | None:
        return self._resolution_cc_half

    @property
    def resolution_rmerge(self) -> float | None:
        return self._resolution_rmerge

    @property
    def resolution_isigma(self) -> float | None:
        return self._resolution_isigma

    @property
    def resolution_misigma(self) -> float | None:
        return self._resolution_misigma
