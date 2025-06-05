from __future__ import annotations

import logging
import os
from pathlib import Path

import iotbx.phil
import libtbx.phil
from dials.array_family import flex
from dials.util.resolution_analysis import Resolutionizer, phil_defaults
from dxtbx.model import ExperimentList

from xia2.Driver.DriverFactory import DriverFactory
from xia2.Driver.timing import record_step
from xia2.lib.bits import _get_number
from xia2.Modules.SSX.util import log_to_file, run_in_directory

logger = logging.getLogger("xia2.Wrappers.XIA.EstimateResolution")


def diff_phil_from_params_and_scope(
    params: libtbx.phil.scope_extract, phil_scope: libtbx.phil.scope
) -> str:
    original = phil_scope.extract()
    diff_phil = ""

    def compare_params(new, original, diff_phil, parent=""):
        for k in [k for k in vars(original).keys() if k[0] != "_"]:
            v2 = getattr(original, k, None)
            v1 = getattr(new, k, None)
            if isinstance(v1, libtbx.phil.scope_extract):
                if parent:
                    diff_phil = compare_params(v1, v2, diff_phil, parent + "." + k)
                else:
                    diff_phil = compare_params(v1, v2, diff_phil, k)
            else:
                if v1 and v1 != v2:
                    diff_phil += f"{parent}.{k} = {v1}\n"
        return diff_phil

    diff_phil = compare_params(params, original, diff_phil)
    pretty = phil_scope.fetch_diff(
        source=phil_scope.fetch(sources=[iotbx.phil.parse(diff_phil)])
    ).as_str()
    pretty = "The following parameters have been modified:\n\n" + pretty
    return pretty


class NewEstimateResolutionWrapper:
    def __init__(self, working_directory: Path):
        self.params: libtbx.phil.scope_extract = phil_defaults.extract()
        self.params.resolution.nbins = 100
        ## Where any output (logfiles, datafiles) will be saved
        self._working_directory: Path = working_directory
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
            with open(logfile, "r") as log:
                for record in log:
                    if "Resolution rmerge" in record:
                        self._resolution_rmerge = float(record.split()[-1])
                    if "Resolution completeness" in record:
                        self._resolution_completeness = float(record.split()[-1])
                    if "Resolution cc_half_significance_level" in record:
                        self._resolution_cc_half_significance = float(
                            record.split()[-1]
                        )
                    elif "Resolution cc_half" in record:
                        self._resolution_cc_half = float(record.split()[-1])
                    if "Resolution I/sig" in record:
                        self._resolution_isigma = float(record.split()[-1])
                    if "Resolution Mn(I/sig)" in record:
                        self._resolution_misigma = float(record.split()[-1])

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


def EstimateResolution(DriverType=None):
    """A factory for EstimateResolutionWrapper classes."""

    DriverInstance = DriverFactory.Driver(DriverType)

    class EstimateResolutionWrapper(DriverInstance.__class__):
        def __init__(self):
            DriverInstance.__class__.__init__(self)
            self.set_executable("dials.estimate_resolution")

            # inputs
            self._hklin = None
            self._reflections = None
            self._experiments = None
            self._limit_rmerge = None
            self._limit_completeness = None
            self._limit_cc_half = None
            self._cc_half_fit = None
            self._cc_half_significance_level = None
            self._limit_isigma = None
            self._limit_misigma = None
            self._nbins = 100
            self._batch_range = None
            self._labels = None

            # outputs
            self._resolution_rmerge = None
            self._resolution_completeness = None
            self._resolution_cc_half = None
            self._resolution_cc_half_significance = None
            self._resolution_isigma = None
            self._resolution_misigma = None
            self._html = None
            self._json = None

        def set_reflections(self, filename):
            self._reflections = filename

        def set_experiments(self, filename):
            self._experiments = filename

        def set_hklin(self, hklin):
            self._hklin = hklin

        def set_nbins(self, nbins):
            self._nbins = nbins

        def set_limit_rmerge(self, limit_rmerge):
            self._limit_rmerge = limit_rmerge

        def set_limit_completeness(self, limit_completeness):
            self._limit_completeness = limit_completeness

        def set_limit_cc_half(self, limit_cc_half):
            self._limit_cc_half = limit_cc_half

        def set_cc_half_fit(self, cc_half_fit):
            self._cc_half_fit = cc_half_fit

        def set_cc_half_significance_level(self, cc_half_significance_level):
            self._cc_half_significance_level = cc_half_significance_level

        def set_limit_isigma(self, limit_isigma):
            self._limit_isigma = limit_isigma

        def set_limit_misigma(self, limit_misigma):
            self._limit_misigma = limit_misigma

        def set_batch_range(self, start, end):
            self._batch_range = (start, end)

        def set_labels(self, labels):
            self._labels = labels

        def get_resolution_rmerge(self):
            return self._resolution_rmerge

        def get_resolution_completeness(self):
            return self._resolution_completeness

        def get_resolution_cc_half(self):
            return self._resolution_cc_half

        def get_resolution_isigma(self):
            return self._resolution_isigma

        def get_resolution_misigma(self):
            return self._resolution_misigma

        def get_resolution_cc_half_significance(self):
            return self._resolution_cc_half_significance

        def get_html(self):
            return self._html

        def get_json(self):
            return self._json

        def run(self):
            assert self._hklin or (self._experiments and self._reflections)
            if self._hklin:
                cl = [self._hklin]
            else:
                cl = [self._experiments, self._reflections]
            cl.append("nbins=%s" % self._nbins)
            cl.append("rmerge=%s" % self._limit_rmerge)
            cl.append("completeness=%s" % self._limit_completeness)
            cl.append("cc_half=%s" % self._limit_cc_half)
            if self._cc_half_fit is not None:
                cl.append("cc_half_fit=%s" % self._cc_half_fit)
            cl.append(
                "cc_half_significance_level=%s" % self._cc_half_significance_level
            )
            cl.append("isigma=%s" % self._limit_isigma)
            cl.append("misigma=%s" % self._limit_misigma)
            if self._batch_range is not None:
                cl.append("batch_range=%i,%i" % self._batch_range)
            if self._labels is not None:
                cl.append("labels=%s" % self._labels)
            for c in cl:
                self.add_command_line(c)
            logger.debug("Resolution analysis: %s", " ".join(cl))

            self._html = os.path.join(
                self.get_working_directory(),
                "%d_dials.estimate_resolution.html" % self.get_xpid(),
            )
            self.add_command_line("output.html=%s" % self._html)

            self._json = os.path.join(
                self.get_working_directory(),
                "%d_dials.estimate_resolution.json" % self.get_xpid(),
            )
            self.add_command_line("output.json=%s" % self._json)

            self.start()
            self.close_wait()
            for record in self.get_all_output():
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

    return EstimateResolutionWrapper()
