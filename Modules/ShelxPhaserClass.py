#!/usr/bin/python
# 
# DEVELOPMENT CODE DO NOT USE
# 
# ShelxPhaserClass.py
# Maintained by G.Winter

import os
import sys

if not os.environ['XIA2_ROOT'] in sys.path:
    sys.path.append(os.environ['XIA2_ROOT'])

from Wrappers.Shelx.Shelxc import Shelxc
from Wrappers.Shelx.Shelxd import Shelxd
from Wrappers.Shelx.Shelxe import Shelxe
from Wrappers.CCP4.F2mtz import F2mtz
from Wrappers.CCP4.Cad import Cad

# interface that this will implement
from HAPhaserClass import HAPhaserClass

class ShelxPhaserClass(HAPhaserClass):

    def __init__(self):
        HAPhaserClass.__init__(self)

        # check we have the things we need
        shelxc = Shelxc()
        shelxd = Shelxd()
        shelxe = Shelxe()

        return

    def phase(self):

        # prepare the data

        shelxc = Shelxc()
        shelxc.write_log_file('shelxc.log')
        shelxc.set_cell(self._cell)
        shelxc.set_symmetry(self._spacegroup)
        shelxc.set_n_sites(self._input_dict['n_sites'])
        
        for wavelength in self._scalepack_files.keys():
            if wavelength.upper() == 'PEAK':
                shelxc.set_peak(self._scalepack_files[wavelength])
            if wavelength.upper() == 'INFL':
                shelxc.set_infl(self._scalepack_files[wavelength])
            if wavelength.upper() == 'LREM':
                shelxc.set_lrem(self._scalepack_files[wavelength])
            if wavelength.upper() == 'HREM':
                shelxc.set_hrem(self._scalepack_files[wavelength])
            if wavelength.upper() == 'SAD':
                shelxc.set_sad(self._scalepack_files[wavelength])

        shelxc.set_name('strawman')
        shelxc.prepare()

        # next look for the sites

        shelxd = Shelxd()
        shelxd.set_name('strawman')
        shelxd.write_log_file('shelxd.log')
        shelxd.set_spacegroup(self._spacegroup)
        sites = shelxd.get_sites()

        # do nothing with these

        shelxe = Shelxe()
        shelxe.write_log_file('shelxe.log')
        shelxe.set_name('strawman')
        shelxe.set_solvent(self._input_dict['solvent'])
        shelxe.phase()

        f = F2mtz()
        f.write_log_file('f2mtz.log')        
        f.set_hklin('strawman.phs')
        f.set_hklout('tmp_strawman.mtz')
        f.set_cell(self._cell)
        f.set_symmetry(self._spacegroup)
        f.f2mtz()

        c = Cad()
        c.write_log_file('cad.log')
        c.add_hklin('tmp_strawman.mtz')
        c.set_hklout('strawman.mtz')
        c.update()

        # phase in enantiomorph spacegroup if enantiomorphic possible

        if self.compute_spacegroup_enantiomorph(self._spacegroup) == \
           self._spacegroup:

            shelxe = Shelxe()
            shelxe.write_log_file('shelxe_oh.log')
            shelxe.set_name('strawman')
            shelxe.set_solvent(self._input_dict['solvent'])
            shelxe.set_enantiomorph()
            shelxe.phase()

            f = F2mtz()
            f.set_hklin('strawman_i.phs')
            f.write_log_file('f2mtz_i.log')        
            f.set_hklout('tmp_strawman_i.mtz')
            f.set_cell(self._cell)
            f.set_symmetry(self._spacegroup)
            f.f2mtz()

            c = Cad()
            c.write_log_file('cad_i.log')
            c.add_hklin('tmp_strawman_i.mtz')
            c.set_hklout('strawman_i.mtz')
            c.update()
            
        return

        
if __name__ == '__main__':
    spc = ShelxPhaserClass()
    spc.phase()

    
        
