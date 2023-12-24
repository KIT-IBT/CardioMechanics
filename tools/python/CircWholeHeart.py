import numpy as np
import pandas as pd
import json
import csv
import time

class CircWholeHeart:
    """
    Closed loop circulation model

    References
    ----------
    Gerach, T.; Schuler, S.; Fröhlich, J.; Lindner, L.; Kovacheva, E.; Moss, R.; Wülfers, E.M.;
    Seemann, G.; Wieners, C.; Loewe, A. 
    Electro-Mechanical Whole-Heart Digital Twins: A Fully Coupled Multi-Physics Approach. 
    Mathematics 2021, 9, 1247. https://doi.org/10.3390/math9111247
    """

    def __init__(self, options = dict()):
        
        if isinstance(options, str):
            with open(options, mode='r', newline='') as inputfile:
                options = json.loads(inputfile.read())

        ############ Heartbeat
        self.BPM = float(options.get('BPM', 75)) # [1 / min]
        self.THB = 60. / self.BPM # [s], Heartbeat period
        
        ############ Chambers
        # LA
        options_curr = options.get('LA', dict())
        EA_LA = float(options_curr.get('EA', 0.07)) # [mmHg / ml]
        EB_LA = float(options_curr.get('EB', 0.09)) # [mmHg / ml]
        TC_LA = float(options_curr.get('TC', 0.17)) * self.THB  # [s]
        TR_LA = float(options_curr.get('TR', 0.17)) * self.THB  # [s]
        tC_LA = float(options_curr.get('tC', 0.80)) * self.THB  # [s]
        self.V0_LA = float(options_curr.get('V0', 4.0)) # [ml]
        self.E_LA = self.time_varying_elastance(EA_LA, EB_LA, tC_LA, TC_LA, TR_LA)

        # LV
        options_curr = options.get('LV', dict())
        EA_LV = float(options_curr.get('EA', 2.75)) # [mmHg / ml]
        EB_LV = float(options_curr.get('EB', 0.08)) # [mmHg / ml]
        TC_LV = float(options_curr.get('TC', 0.34)) * self.THB  # [s]
        TR_LV = float(options_curr.get('TR', 0.17)) * self.THB  # [s]
        tC_LV = float(options_curr.get('tC', 0.00)) * self.THB  # [s]
        self.V0_LV = float(options_curr.get('V0', 5.0)) # [ml]
        self.E_LV = self.time_varying_elastance(EA_LV, EB_LV, tC_LV, TC_LV, TR_LV)

        # RA
        options_curr = options.get('RA', dict())
        EA_RA = float(options_curr.get('EA', 0.06)) # [mmHg / ml]
        EB_RA = float(options_curr.get('EB', 0.07)) # [mmHg / ml]
        TC_RA = float(options_curr.get('TC', 0.17)) * self.THB  # [s]
        TR_RA = float(options_curr.get('TR', 0.17)) * self.THB  # [s]
        tC_RA = float(options_curr.get('tC', 0.80)) * self.THB  # [s]
        self.V0_RA = float(options_curr.get('V0', 4.0)) # [ml]
        self.E_RA = self.time_varying_elastance(EA_RA, EB_RA, tC_RA, TC_RA, TR_RA)

        # RV
        options_curr = options.get('RV', dict())
        EA_RV = float(options_curr.get('EA', 0.55)) # [mmHg / ml]
        EB_RV = float(options_curr.get('EB', 0.05)) # [mmHg / ml]
        TC_RV = float(options_curr.get('TC', 0.34)) * self.THB  # [s]
        TR_RV = float(options_curr.get('TR', 0.17)) * self.THB  # [s]
        tC_RV = float(options_curr.get('tC', 0.00)) * self.THB  # [s]
        self.V0_RV = float(options_curr.get('V0', 10.0)) # [ml]
        self.E_RV = self.time_varying_elastance(EA_RV, EB_RV, tC_RV, TC_RV, TR_RV)

        ############ Valves
        options_curr = options.get('valves', dict())
        self.R_MV = float(options_curr.get('R_MV', 0.003)) # [mmHg s / ml]
        self.R_AV = float(options_curr.get('R_AV', 0.006)) # [mmHg s / ml]
        self.R_TV = float(options_curr.get('R_TV', 0.003)) # [mmHg s / ml]
        self.R_PV = float(options_curr.get('R_PV', 0.003)) # [mmHg s / ml]

        ############ Systemic circulation
        options_curr = options.get('SYS', dict())
        self.R_AR_SYS   = float(options_curr.get('R_AR'  , 0.05  )) # [mmHg s /ml]
        self.C_AR_SYS   = float(options_curr.get('C_AR'  , 1.5   )) # [ml / mmHg]
        self.V0_AR_SYS  = float(options_curr.get('V0_AR' , 800.0 )) # [ml]
        self.R_PER_SYS  = float(options_curr.get('R_PER' , 0.9   )) # [mmHg s /ml]
        self.R_VEN_SYS  = float(options_curr.get('R_VEN' , 0.015 )) # [mmHg s /ml]
        self.C_VEN_SYS  = float(options_curr.get('C_VEN' , 90.   )) # [ml / mmHg]
        self.V0_VEN_SYS = float(options_curr.get('V0_VEN', 2830.0)) # [ml]

        ############ Pulmonary circulation
        options_curr = options.get('PUL', dict())
        self.R_AR_PUL   = float(options_curr.get('R_AR'  , 0.015)) # [mmHg s /ml]
        self.C_AR_PUL   = float(options_curr.get('C_AR'  , 6.5  )) # [ml / mmHg]
        self.V0_AR_PUL  = float(options_curr.get('V0_AR' , 120.0)) # [ml]
        self.R_PER_PUL  = float(options_curr.get('R_PER' , 0.07 )) # [mmHg s /ml]
        self.R_VEN_PUL  = float(options_curr.get('R_VEN' , 0.015)) # [mmHg s /ml]
        self.C_VEN_PUL  = float(options_curr.get('C_VEN' , 16.  )) # [ml / mmHg]
        self.V0_VEN_PUL = float(options_curr.get('V0_VEN', 150.0)) # [ml]

        ############ PV relationships
        self.p_LA_func = lambda V, t: self.E_LA(t) * ( V - self.V0_LA )
        self.p_LV_func = lambda V, t: self.E_LV(t) * ( V - self.V0_LV )
        self.p_RA_func = lambda V, t: self.E_RA(t) * ( V - self.V0_RA )
        self.p_RV_func = lambda V, t: self.E_RV(t) * ( V - self.V0_RV )

    def flux_through_valve(self, p1, p2, R):
        return max( (( p1 - p2 ) / R), 0)

    def flux(self, p1, p2, R):
        return ( p1 - p2 ) / R

    def time_varying_elastance(self, EA, EB, time_C, duration_C, duration_R):
        time_R = time_C + duration_C
        e = lambda t: 0.5 * ( 1 - np.cos( np.pi / duration_C * ( np.mod( t - time_C, self.THB ) ) ) ) * ( 0 <= np.mod( t - time_C, self.THB ) ) * ( np.mod( t - time_C, self.THB ) < duration_C ) + \
                      0.5 * ( 1 + np.cos( np.pi / duration_R * ( np.mod( t - time_R, self.THB ) ) ) ) * ( 0 <= np.mod( t - time_R, self.THB ) ) * ( np.mod( t - time_R, self.THB ) < duration_R )
        return lambda t: EA * np.clip(e(t), 0.0, 1.0) + EB

    def initialize(self, initial_state = dict()):
        
        if isinstance(initial_state, str):
            with open(initial_state, mode='r', newline='') as inputfile:
                initial_state = json.loads(inputfile.read())

        self.V_LA      = float(initial_state.get('V_LA'     ,   70.)) #  [ml]
        self.V_LV      = float(initial_state.get('V_LV'     ,  200.)) #  [ml]
        self.V_RA      = float(initial_state.get('V_RA'     ,   98.)) #  [ml]
        self.V_RV      = float(initial_state.get('V_RV'     ,  188.)) #  [ml]

        self.V_AR_SYS  = float(initial_state.get('V_AR_SYS' , 1050.)) #  [ml]
        self.V_VEN_SYS = float(initial_state.get('V_VEN_SYS', 3510.)) #  [ml]
        self.V_AR_PUL  = float(initial_state.get('V_AR_PUL' ,  200.)) #  [ml]
        self.V_VEN_PUL = float(initial_state.get('V_VEN_PUL',  320.)) #  [ml]

        self.update_static_variables(0.)

    def update_static_variables(self, t):
        self.p_LA       = self.p_LA_func(self.V_LA, t)
        self.p_LV       = self.p_LV_func(self.V_LV, t)
        self.p_RA       = self.p_RA_func(self.V_RA, t)
        self.p_RV       = self.p_RV_func(self.V_RV, t)

        self.p_AR_SYS_C = (self.V_AR_SYS - self.V0_AR_SYS) / self.C_AR_SYS
        self.Q_AV       = self.flux_through_valve( self.p_LV, self.p_AR_SYS_C, self.R_AV + self.R_AR_SYS )
        self.p_AR_SYS   = self.p_LV - self.R_AV * self.Q_AV

        self.p_VEN_SYS  = (self.V_VEN_SYS - self.V0_VEN_SYS) / self.C_VEN_SYS
        self.Q_PER_SYS  = self.flux(self.p_AR_SYS_C, self.p_VEN_SYS, self.R_PER_SYS)
        self.Q_VEN_SYS  = self.flux(self.p_VEN_SYS, self.p_RA, self.R_VEN_SYS)
        self.Q_TV       = self.flux_through_valve( self.p_RA, self.p_RV, self.R_TV )

        self.p_AR_PUL_C = (self.V_AR_PUL - self.V0_AR_PUL) / self.C_AR_PUL
        self.Q_PV       = self.flux_through_valve( self.p_RV, self.p_AR_PUL_C, self.R_PV + self.R_AR_PUL )
        self.p_AR_PUL   = self.p_RV - self.R_PV * self.Q_PV
        self.p_VEN_PUL  = (self.V_VEN_PUL - self.V0_VEN_PUL) / self.C_VEN_PUL

        self.Q_PER_PUL  = self.flux(self.p_AR_PUL_C, self.p_VEN_PUL, self.R_PER_PUL)
        self.Q_VEN_PUL  = self.flux(self.p_VEN_PUL, self.p_LA, self.R_VEN_PUL)
        self.Q_MV       = self.flux_through_valve( self.p_LA, self.p_LV    , self.R_MV )
        
    def solve_step_FE(self, t, dt):
        self.update_static_variables(t)

        self.V_LA      += dt * ( self.Q_VEN_PUL - self.Q_MV )
        self.V_LV      += dt * ( self.Q_MV      - self.Q_AV )
        self.V_RA      += dt * ( self.Q_VEN_SYS - self.Q_TV )
        self.V_RV      += dt * ( self.Q_TV      - self.Q_PV )
        self.V_AR_SYS  += dt * ( self.Q_AV      - self.Q_PER_SYS )
        self.V_VEN_SYS += dt * ( self.Q_PER_SYS - self.Q_VEN_SYS )
        self.V_AR_PUL  += dt * ( self.Q_PV      - self.Q_PER_PUL )
        self.V_VEN_PUL += dt * ( self.Q_PER_PUL - self.Q_VEN_PUL )
        
    def solve(self, T = None, num_cycles = None,
              initial_state = None,
              dt = 1e-3,
              dt_eval = None):

        print('Circulation model - running simulation...')
        if (T is None and num_cycles is None) or (T is not None and num_cycles is not None):
            raise Exception('Exactly one among T and num_cycles should be not None.')

        if num_cycles is not None:
            T = self.THB * num_cycles
        if dt_eval is None:
            output_every_n_steps = 1
        else:
            output_every_n_steps = np.round(dt_eval / dt)
        times = np.arange(0, T, dt)

        self.initialize(initial_state = initial_state)
        self.initialize_output()
        self.dump_output(0.0)

        time_start = time.time()

        for iT in range(1, times.shape[0]):
            self.solve_step_FE(times[iT], dt)
            if iT % output_every_n_steps == 0:
                self.dump_output(times[iT])

        duration = time.time() - time_start

        print('Circulation model - elapsed time %1.4f s' % duration)
        return pd.DataFrame(self.results)

    def initialize_output(self):
        self.results = dict()
        self.results['time']    = list()
        self.results['VLA']     = list()
        self.results['VLV']     = list()
        self.results['VRA']     = list()
        self.results['VRV']     = list()
        self.results['VARSYS']  = list()
        self.results['VVENSYS'] = list()
        self.results['VARPUL']  = list()
        self.results['VVENPUL'] = list()
        self.results['pARSYS']  = list()
        self.results['pVENSYS'] = list()
        self.results['pARPUL']  = list()
        self.results['pVENPUL'] = list()
        self.results['QPERSYS'] = list()
        self.results['QVENSYS'] = list()
        self.results['QPERPUL'] = list()
        self.results['QVENPUL'] = list()
        self.results['pLA']     = list()
        self.results['pLV']     = list()
        self.results['pRA']     = list()
        self.results['pRV']     = list()
        self.results['ELA']     = list()
        self.results['ELV']     = list()
        self.results['ERA']     = list()
        self.results['ERV']     = list()
        self.results['QMV']     = list()
        self.results['QAV']     = list()
        self.results['QTV']     = list()
        self.results['QPV']     = list()

    def dump_output(self, t):
        self.results['time'   ].append(t)
        self.results['VLA'    ].append(self.V_LA)
        self.results['VLV'    ].append(self.V_LV)
        self.results['VRA'    ].append(self.V_RA)
        self.results['VRV'    ].append(self.V_RV)
        self.results['VARSYS' ].append(self.V_AR_SYS)
        self.results['VVENSYS'].append(self.V_VEN_SYS)
        self.results['VARPUL' ].append(self.V_AR_PUL)
        self.results['VVENPUL'].append(self.V_VEN_PUL)
        self.results['pARSYS' ].append(self.p_AR_SYS)
        self.results['pVENSYS'].append(self.p_VEN_SYS)
        self.results['pARPUL' ].append(self.p_AR_PUL)
        self.results['pVENPUL'].append(self.p_VEN_PUL)
        self.results['QPERSYS'].append(self.Q_PER_SYS)
        self.results['QVENSYS'].append(self.Q_VEN_SYS)
        self.results['QPERPUL'].append(self.Q_PER_PUL)
        self.results['QVENPUL'].append(self.Q_VEN_PUL)
        self.results['pLA'    ].append(self.p_LA)
        self.results['pLV'    ].append(self.p_LV)
        self.results['pRA'    ].append(self.p_RA)
        self.results['pRV'    ].append(self.p_RV)
        self.results['ELA'    ].append(self.E_LA(t))
        self.results['ELV'    ].append(self.E_LV(t))
        self.results['ERA'    ].append(self.E_RA(t))
        self.results['ERV'    ].append(self.E_RV(t))
        self.results['QMV'    ].append(self.Q_MV)
        self.results['QAV'    ].append(self.Q_AV)
        self.results['QTV'    ].append(self.Q_TV)
        self.results['QPV'    ].append(self.Q_PV)

    def save_state(self, filename):
        
        with open(filename, mode='w', newline='') as outfile:
            state = dict()
            state['V_LA']      = float(self.V_LA)
            state['V_LV']      = float(self.V_LV)
            state['V_RA']      = float(self.V_RA)
            state['V_RV']      = float(self.V_RV)
            state['p_LA']      = float(self.p_LA)
            state['p_LV']      = float(self.p_LV)
            state['p_RA']      = float(self.p_RA)
            state['p_RV']      = float(self.p_RV)
            state['p_AR_SYS']  = float(self.p_AR_SYS)
            state['p_VEN_SYS'] = float(self.p_VEN_SYS)
            state['p_AR_PUL']  = float(self.p_AR_PUL)
            state['p_VEN_PUL'] = float(self.p_VEN_PUL)
            state['Q_PER_SYS'] = float(self.Q_PER_SYS)
            state['Q_VEN_SYS'] = float(self.Q_VEN_SYS)
            state['Q_PER_PUL'] = float(self.Q_PER_PUL)
            state['Q_VEN_PUL'] = float(self.Q_VEN_PUL)
            json.dump(state, outfile, indent=2)

    def print_info(self):
        
        print('V_LA      = %4.2f mL' % self.V_LA)
        print('V_LV      = %4.2f mL' % self.V_LV)
        print('V_RA      = %4.2f mL' % self.V_RA)
        print('V_RV      = %4.2f mL' % self.V_RV)
        print('V_AR_SYS  = %4.2f mL' % self.V_AR_SYS)
        print('V_VEN_SYS = %4.2f mL' % self.V_VEN_SYS)
        print('V_AR_PUL  = %4.2f mL' % self.V_AR_PUL)
        print('V_VEN_PUL = %4.2f mL' % self.V_VEN_PUL)

        V_tot_HEART = self.V_LA + self.V_RA + self.V_LV + self.V_RV
        V_tot_SYS   = self.V_AR_SYS + self.V_VEN_SYS
        V_tot_PUL   = self.V_AR_PUL + self.V_VEN_PUL
        V_tot       = V_tot_HEART + V_tot_SYS + V_tot_PUL

        print('======================')
        print('V (HEART) = %4.2f mL' % V_tot_HEART)
        print('V (SYS)   = %4.2f mL' % V_tot_SYS)
        print('V (PUL)   = %4.2f mL' % V_tot_PUL)
        print('======================')
        print('V         = %4.2f mL' % V_tot)

