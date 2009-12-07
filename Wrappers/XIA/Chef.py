#!/usr/bin/env python
# Chef.py
#   Copyright (C) 2008 CCLRC, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is 
#   included in the root directory of this package.
#
# 5th February 2008
#
# A wrapper for the new program "chef". This has been developed for xia2
# to analyse the bulk properties of intensity measurements, particularly
# looking at how well they agree. The idea is that reflection files with
# DOSE columns added in by DOSER may be inspected to determine the 
# dose / resolution envelope optimal for given analysis processes, viz:
# 
# - substructure determination
# - phase calculation
# - density modification & refinement
# 
# This should give "proper" resolution limits...

import os
import sys
import math

if not os.environ.has_key('XIA2CORE_ROOT'):
    raise RuntimeError, 'XIA2CORE_ROOT not defined'

if not os.environ.has_key('XIA2_ROOT'):
    raise RuntimeError, 'XIA2_ROOT not defined'

if not os.path.join(os.environ['XIA2CORE_ROOT'], 'Python') in sys.path:
    sys.path.append(os.path.join(os.environ['XIA2CORE_ROOT'], 'Python'))

if not os.environ['XIA2_ROOT'] in sys.path:
    sys.path.append(os.environ['XIA2_ROOT'])

from Driver.DriverFactory import DriverFactory
from Decorators.DecoratorFactory import DecoratorFactory
from lib.Guff import transpose_loggraph, mean_sd
from Wrappers.CCP4.Mtzdump import Mtzdump

def Chef(DriverType = None):
    '''A factory for wrappers for the chef.'''

    DriverInstance = DriverFactory.Driver(DriverType)
    CCP4DriverInstance = DecoratorFactory.Decorate(DriverInstance, 'ccp4')

    class ChefWrapper(CCP4DriverInstance.__class__):
        '''Provide access to the functionality in chef.'''

        def __init__(self):

            CCP4DriverInstance.__class__.__init__(self)

            self.set_executable('pychef')

            self._hklin_list = []
            self._anomalous = False
            self._b_width = 0.0
            self._b_max = 0.0
            self._b_labin = None
            self._resolution = 0.0

            self._p_crd = True

            self._completeness = { }
            
            self._title = None

            return

        def add_hklin(self, hklin):
            self._hklin_list.append(hklin)
            return

        def set_anomalous(self, anomalous):
            self._anomalous = anomalous
            return

        def set_resolution(self, resolution):
            self._resolution = resolution
            return

        def set_width(self, width):
            self._b_width = width
            return

        def set_max(self, max):
            self._b_max = max
            return

        def set_labin(self, labin):
            self._b_labin = labin
            return

        def set_title(self, title):
            self._title = title
            return

        def get_completeness(self, wavelength):
            return self._completeness[wavelength]

        def get_completeness_datasets(self):
            return self._completeness.keys()

        def run(self):
            '''Actually run chef...'''

            if len(self._hklin_list) == 0:
                raise RuntimeError, 'HKLIN not defined'

            for j in range(len(self._hklin_list)):
                self.add_command_line('HKLIN%d' % (j + 1))
                self.add_command_line(self._hklin_list[j])

            self.start()

            # this is not needed for pychef
            # self.input('print scp comp')

            if self._anomalous:
                self.input('anomalous on')
            if self._b_width > 0.0:
                self.input('range width %f' % self._b_width)
            if self._b_max > 0.0:
                self.input('range max %f' % self._b_max)

            if self._resolution > 0.0:
                self.input('resolution %.2f' % self._resolution)

            self.input('labin BASE=%s' % self._b_labin)

            if self._title:
                self.input('title %s' % self._title)

            self.close_wait()

            # FIXME should check the status here...

            # read out the completeness curves...

            output = self.get_all_output()

            all_doses = []

            for j in range(len(output)):
                record = output[j]
                if 'Completeness vs. BASELINE' in record:
                    dataset = record.split()[-1]
                    completeness = []
                    k = j + 2
                    record = output[k]
                    while not 'Expected' in record and not '$TABLE' in record:
                        completeness.append((float(record.split()[0]),
                                             float(record.split()[-1])))
                        dose = float(record.split()[0])

                        if not dose in all_doses:
                            all_doses.append(dose)

                        k += 1
                        record = output[k]

                        
                    self._completeness[dataset] = completeness

            # now jimmy these..

            for dataset in self._completeness.keys():
                completeness = self._completeness[dataset]
                cmax = completeness[-1][1]
                cnew = []

                # hash this up
                ctable = { }
                for c in completeness:
                    ctable[c[0]] = c[1]
                
                for dose in all_doses:
                    if dose in ctable:
                        cnew.append((dose, ctable[dose]))
                    else:
                        cnew.append((dose, cmax))

                self._completeness[dataset] = cnew

                # at some point need to figure out how to analyse these
                # results...

            return

        def digest_rd(self, values):
            '''Digest the results of an Rd calculation, working on the
            assumptions that (i) the corresponding dose values are
            meaningless and (ii) we are trying to decide if there is a
            significant gradient there. N.B. does however assume that the
            dose increments are UNIFORM.'''

            sx = 0.0
            sy = 0.0

            n = 0

            for j, v in enumerate(values):

                if not v:
                    continue

                sx += j
                sy += v

                n += 1

            mx = sx / n
            my = sy / n

            sxx = 0.0
            sxy = 0.0

            for j, v in enumerate(values):

                if not v:
                    continue

                sxx += (j - mx) * (j - mx)
                sxy += (j - mx) * (v - my)

            m = sxy / sxx
            c = my - m * mx

            # now calculate residual about this line            

            ss = 0.0

            for j, v in enumerate(values):

                if not v:
                    continue
                
                _v = m * j + c

                ss += (v - _v) * (v - _v)

            sd = math.sqrt(ss / (n - 2))

            # then compute the standard deviation of the population

            var = 0.0

            for j, v in enumerate(values):
                
                if not v:
                    continue

                var += (v - my) * (v - my)

            sigma = math.sqrt(var / (n - 1))

            # return sigma / sd

            return (var / (sd * sd)) / n
        
        def parse(self):
            '''Parse the output of the chef run.'''

            results = self.parse_ccp4_loggraph()

            rd_keys = []
            comp_keys = []

            scp_data = None
            comp_data = { }
            rd_data = { }

            for key in results:
                if 'Completeness vs. ' in key:
                    comp_keys.append(key)
                    comp_data[key.split()[-1]] = transpose_loggraph(
                        results[key])

                elif 'R vs. ' in key:
                    rd_keys.append(key)
                    rd_data[key.split()[-1]] = transpose_loggraph(
                        results[key])                    

                    values = map(float, rd_data[key.split()[-1]]['2_Rd'])
                    print '%s %.2f' % (key.split()[-1], self.digest_rd(values))
                    
                              
                elif 'Cumulative radiation' in key:
                    scp_data = transpose_loggraph(results[key])

            # right, so first work through these to define the limits from
            # where the first set is 50% complete to 90% complete, which
            # will establish the benchmark, then calculate a kinda
            # Z-score for the subsequent Scp values

            lowest_50 = None
            lowest_90 = None

            for dataset in comp_data:
                i_col = '2_I'
                if '5_dI' in comp_data[dataset]:
                    i_col = '4_I'

                completeness = comp_data[dataset][i_col]

                local_50 = None
                local_90 = None

                for j, dose in enumerate(comp_data[dataset]['1_DOSE']):
                    
                    comp = float(completeness[j])

                    if comp > 0.50 and not local_50:
                        local_50 = float(dose)

                    if comp > 0.9 and not local_90:
                        local_90 = float(dose)

            if not lowest_50:
                lowest_50 = local_50
            if local_50 < lowest_50:
                lowest_50 = local_50

            if not lowest_90:
                lowest_90 = local_90
            if local_90 < lowest_90:
                lowest_90 = local_90

            # now build up the reference population

            scp_reference = []

            scp_key = None

            for k in scp_data:
                if 'Scp(d)' in k:
                    scp_key = k
                    
            for j, d in enumerate(scp_data['1_DOSE']):
                dose = float(d)
                if dose >= lowest_50 and dose <= lowest_90:
                    scp_reference.append(float(scp_data[scp_key][j]))

            m, s = mean_sd(scp_reference)

            dose = scp_data['1_DOSE'][0]

            scp_max = 0.0

            for j, d in enumerate(scp_data['1_DOSE']):

                dose = float(d)
                scp = float(scp_data[scp_key][j])
                z = (scp - m) / s

                if dose < lowest_90:
                    scp_max = max(scp, scp_max)
                    continue

                if math.fabs(z) < 1:
                    # add to the population
                    scp_reference.append(scp)
                    m, s = mean_sd(scp_reference)

                if z > 3 and scp > scp_max:
                    break

                scp_max = max(scp, scp_max)

            if dose == float(scp_data['1_DOSE'][-1]):
                print 'Dose limit: use all data'
            else:
                print 'Dose limit: %.1f' % dose
                
    return ChefWrapper()
        
if __name__ == '__main_exec__':
    # then run a test...

    source = os.path.join(os.environ['X2TD_ROOT'], 'Test', 'Chef',
                          'TestData')

    # first find the maximum dose... and minimum resolution range

    dmax = 0.0
    dmin = 0.0

    for hklin in ['TS03_12287_doser_INFL.mtz',
                  'TS03_12287_doser_LREM.mtz',
                  'TS03_12287_doser_PEAK.mtz']:

        md = Mtzdump()
        md.set_hklin(os.path.join(source, hklin))
        md.dump()
        dmax = max(dmax, max(md.get_column_range('DOSE')[:2]))
        dmin = max(dmin, md.get_column_range('DOSE')[2])

    chef = Chef()
    chef.write_log_file('chef.log')

    for hklin in ['TS03_12287_doser_INFL.mtz',
                  'TS03_12287_doser_LREM.mtz',
                  'TS03_12287_doser_PEAK.mtz']:
        chef.add_hklin(os.path.join(source, hklin))

    chef.set_anomalous(True)
    chef.set_width(5.0)
    chef.set_max(1400)
    chef.set_resolution(dmin)
    chef.set_labin('DOSE')

    chef.run()

    # now gather up the completeness curves to figure out the starting
    # point for the analysis

    wavelengths = chef.get_completeness_datasets()

    start_min = 1440

    for w in wavelengths:
        completeness = chef.get_completeness(w)
        for comp in completeness:
            d, c = comp

            if c > 0.5:
                if start_min > d:
                    start_min = d
                break

    # now code this up for MAD - want two relatively complete wavelengths...
    # slightly hacky - add the completeness values together then stop when
    # this gets to say 1.9... assumes we have an inscribed circle on the
    # detector... bugger, this is slightly complicated by the fact that the
    # completeness is not written out beyond the maximum dose. Fudge this.

    total_comp = { }

    for w in wavelengths:
        completeness = chef.get_completeness(w)
        for comp in completeness:
            d, c = comp

            if not d in total_comp:
                total_comp[d] = 0.0

            total_comp[d] += c

    doses = sorted(total_comp.keys())

    end_min = 0.0

    for d in doses:
        if total_comp[d] > 1.8:
            end_min = d
            break

    print 'Establish the baseline from %f to %f' % (start_min, end_min)
        
if __name__ == '__main__':

    chef = Chef()

    # hack for testing - copy input from elsewhere?!

    for record in open(sys.argv[1]):
        chef._standard_output_records.append(record)

    chef.parse()

    
