from dependencies import split, remove_elem
import pandas as pd
import numpy as np
import math
import constants
import preferences

#Import fission and nubase data. Define prefix and unit conversions.
fiss_data = pd.read_csv('./databases/fyu235thermal.txt', sep="	", header=0)
fiss_data.columns = ["Z", "A", "Level", "YI", "YI uncert"]
fiss_data = fiss_data.sort_values('A')
Z, A, Level, Yeild, Yeild_uncert = fiss_data['Z'].tolist(), fiss_data['A'].tolist(), fiss_data['Level'].tolist(), fiss_data['YI'].tolist(), fiss_data['YI uncert'].tolist()
nubase = pd.read_excel('./databases/Parsed_nubase2016.xlsx', header = 0)
nuA, nuZ= nubase['A'].tolist(), nubase['Z'].tolist()
prefix = ['y','z','a','f','p','n','u','m','c','d','da','h','k','M','G','T','P','E','Z','Y']
prefix_conv = [10**(-24),10**(-21),10**(-18),10**(-15),10**(-12),10**(-9),10**(-6),10**(-3),10**(-2),10**(-1),10**(0),10**(1),10**(2),10**(3),10**(6),10**(9),10**(12),10**(15),10**(18),10**(21),10**(24)]
unit = ['s','m','h','d','y']
unit_conv = [1,60,3600,24*3600,365*24*3600]
energy_units = 'keV'
#'is' in decay mode is natural abundance
#? indicates decay mode is energetically allowed
#=? indicates

class isomer():
    MISSING_NUBASE_ISO = []
    def __init__(self,Z,A,level,yeild,yeild_uncert):
        '''
        Initialise class element with fission yield data
        
        @param: self
        @param: Z, float - Proton number
        @param: A, float - Atomic number
        @param: level, float - isomer level
        @param: yeild, float - spelt wrong to avoid collision with python yield function. Fission yield
        @param: yeild_uncert, float - Fission yield uncertainty 
        '''
        self.Z = Z
        self.A = A
        self.N = self.A-self.Z
        self.level = level
        self.yeild = yeild
        self.yeild_uncert = yeild_uncert
        self.NU = int(self.Z*10000000+self.A*10000+level)
    
    def find_nubase(self):
        '''
        Search through nubase data to match atomic number and proton number. Once found, grab all corresponding
        Z values for all branch modes, including excited states.

        @param: self

        generates:
            self.index - int, index of nubase dataframe where A and Z first meet input conditions
            self.nubase_data - dataframe, nubase data for only the corresponding A and Z values, including excitation energy information.
        '''
        self.index = nuA.index(self.A)
        self.index = nuZ[self.index:].index(self.Z)+self.index
        for i,val in enumerate(nuZ[self.index:]):
            if val != self.Z:
                end_index = i
                break
        self.nubase_data = nubase.iloc[self.index:self.index+end_index]
        return self

    def excit_states(self):
        '''
        Labels nubase data with excitation levels, then matches the excitation level in nubase with input level.

        @param: self

        generates:
            self.states - list, energy levels from nubase data
            self.excit_index - int, generates the excitation index of the input level. If level = 0 then this name is misleading as it is ground state.
            self.isomer_states - dataframe, isomer rows
        '''
        rows = len(self.nubase_data)
        self.states = np.zeros(rows).tolist()
        excit_col = self.nubase_data['Excitation energy']
        for i in range(len(excit_col)):
            if i>0:
                if math.isnan(float(excit_col[self.index+i])) == False or float(excit_col[self.index+i])>float(excit_col[self.index+i-1]):
                    self.states[i] = self.states[i-1]+1
                if float(excit_col[self.index+i]) == float(excit_col[self.index+i-1]):
                    self.states[i] = self.states[i-1]
        level_index = self.states.index(self.level)
        level_count = self.states.count(self.level)
        self.ground_level_count = self.states.count(0)
        self.excit_index = self.index+level_index
        self.isomer_states = self.nubase_data.loc[self.excit_index:self.excit_index+level_count-1]
        return self
    
    def half_life(self):
        '''
        Goes through half lives, converts them to units of seconds, as well as erasing any unwanted symbols.
        If the value is a ?, output is NaN. Yet to impliment flags for fixing later.

        @param: self

        generates:
            self.HL - float, Half life in seconds
            self.HL_uncert - float, Half life uncertainty in seconds
            self.HL_unit - str, Half life units in seconds 's'
            self.excit_HL - float, isomer half life, returns nan if level is not 0, in seconds
            self.excit_HL_uncert - float, isomer half life uncertainty, returns nan if level is not 0, in seconds
            self.excit_HL_unit - str, isomer half life unit, returns nan if level is not 0, in seconds
        '''
        if len(self.isomer_states['Decay mode'].tolist())>1 and self.isomer_states['Decay mode'].loc[self.excit_index] == 'is':
            self.excit_index = self.excit_index + 1
        self.HL = self.isomer_states['Half life'].loc[self.excit_index]
        self.HL_uncert = self.isomer_states['Half life uncert'].loc[self.excit_index]
        self.HL_unit = self.isomer_states['unit'].loc[self.excit_index]
        if self.level>0:
            self.excit_HL = self.isomer_states['Half life 2'].loc[self.excit_index]
            self.excit_HL_uncert = self.isomer_states['Half life uncert 2'].loc[self.excit_index]
            self.excit_HL_unit = self.isomer_states['unit 2'].loc[self.excit_index]
        else:
            self.excit_HL = np.nan
            self.excit_HL_uncert = np.nan
            self.excit_HL_unit = np.nan

        def fix_val(val,uni):
            '''
            Fixes value by removing unwanted symbols and converting to seconds.

            @param: val, str
            @param: uni, str

            @return: val_float, float
            '''
            unwanted_sym = ['?','>','<','~',' ']
            split_val = split(val)
            for i,sym in enumerate(unwanted_sym):
                if sym in split_val:
                    split_val.remove(sym)
            for i,value in enumerate(split_val):
                if value != '.':
                    try:
                        float(value)
                        if i == len(split_val)-1:
                            val_float = ''
                            val_float = float(val_float.join(split_val))
                            unit_split = split(uni)
                            if len(unit_split)>1:
                                val_float = val_float*prefix_conv[prefix.index(unit_split[0])]*unit_conv[unit.index(unit_split[1])]
                            else:
                                val_float = val_float*unit_conv[unit.index(unit_split[0])]
                        elif split_val[-1] == '.':
                            split_val.remove('.')
                            val_float = ''
                            val_float = float(val_float.join(split_val))
                            unit_split = split(uni)
                            if len(unit_split)>1:
                                val_float = val_float*prefix_conv[prefix.index(unit_split[0])]*unit_conv[unit.index(unit_split[1])]
                            else:
                                val_float = val_float*unit_conv[unit.index(unit_split[0])]
                    except ValueError:
                        from_end = i-len(split_val)
                        val_float = ''
                        val_float = float(val_float.join(split_val[:from_end]))
                        unit_list = split_val[from_end:]
                        if len(unit_list)>1:
                            val_float = val_float*prefix_conv[prefix.index(unit_list[0])]*unit_conv[unit.index(unit_list[1])]
                        elif len(unit_list) == 1:
                            val_float = val_float*unit_conv[unit.index(unit_list[0])]
                        else:
                            unit_split = split(uni)
                            if len(unit_split)>1:
                                val_float = val_float*prefix_conv[prefix.index(unit_split[0])]*unit_conv[unit.index(unit_split[1])]
                            else:
                                val_float = val_float*unit_conv[unit.index(unit_split[0])]
                        break

            return val_float

        def check_nan(val):
            '''
            Checks if value is nan, flags as True if it is.

            @param: val, str
            
            @returns: flag
            '''
            try:
                value = float(val)
                if math.isnan(value) == True:
                    flag = True
                else:
                    flag = False
            except ValueError:
                flag = False
            return flag
        
        def convert(val,uni):
            '''
            Uses fix_val function after checking value isn't nan, stable or ?, else it returns nan

            @param: val, str
            @param: uni, str

            @return: val, float
            @return: fix_unit, str
            '''
            if val != 'stable' and val != '?' and check_nan(val) == False:
                val = fix_val(val,uni)
                fix_unit = 's'
            elif val =='?':
                #Should change later. Assuming goes to stable
                val = np.inf
                fix_unit = np.nan
            else:
                val = np.inf
                fix_unit = np.nan
            return val, fix_unit
        
        self.HL,_ = convert(self.HL,self.HL_unit)
        self.HL_uncert,self.HL_unit = convert(self.HL_uncert,self.HL_unit)
        self.excit_HL,_ = convert(self.excit_HL,self.excit_HL_unit)
        self.excit_HL_uncert,self.excit_HL_unit = convert(self.excit_HL_uncert,self.excit_HL_unit)

        if self.level == 0:
            self.halfLife = self.HL
        else:
            self.halfLife = self.excit_HL
        return self
        
    
    def missing_isomer(self):
        '''
        If there is a missing isomer in nubase, it's atomic and proton numbers will be saved as a tuple in MISSING_NUBASE_ISO list

        @param: self

        updates:
            self.MISSING_NUBASE_ISO, list
            self.nubase_data, empty list []
        '''
        if self.NU not in self.MISSING_NUBASE_ISO:
            self.MISSING_NUBASE_ISO.append(self.NU)
            self.nubase_data = []
        return self
    
    def decay_info(self):
        '''
        Get decay information, including decay modes and branch ratios. If the isomer is in an
        excited state, it will take the ground state data incase excited states transition to
        ground state. Filters out # and flags as an estimate. If decay mode is 'is', then it is
        stable and the branch ratio is the natural abundance. This is accounted for.

        @param: self  

        generates:
            self.ground_dacay_modes, list - ground state decay modes
            self.decay_modes, list - excited decay modes, returns empty list if isomer input is ground.
            self.ground_branch_ratio, list - ground state branch ratios
            self.ground_branch_ratio_uncert, list - corresponding uncertainty
            self.branch_ratio, list - excited states decay branch ratios, returns empty list if isomer input is ground.
            self.branch_ratio_uncert, list - corresponding uncertainty.
            self.ground_states, list - ground state information.
            self.excit_energy_allowed, list - exited state energy allows decay mode
            self.excit_estimate_flags, list - excited state branch ratio is an estimate
            self.ground_energy_allowed, list - ground state energy allows decay mode
            self.ground_estimate_flags, list - ground state estimate flag
        '''
        if self.level == 0:
            self.ground_decay_modes = [val for val in self.isomer_states['Decay mode']]
            self.decay_modes = []
            self.ground_branch_ratio = [val for val in self.isomer_states['Branch ratio/Isoptope abund']]
            self.branch_ratio = []
            self.ground_branch_ratio_uncert = [val for val in self.isomer_states['Branch ratio/Isoptope abund uncert']]
            self.branch_ratio_uncert = []
            if 'is' in self.ground_decay_modes:
                isIndex = self.ground_decay_modes.index('is')
                self.natural_abundance = self.ground_branch_ratio[isIndex]
                self.natural_abundance_uncert = self.ground_branch_ratio_uncert[isIndex]
                del self.ground_decay_modes[isIndex]
                del self.ground_branch_ratio[isIndex]
                del self.ground_branch_ratio_uncert[isIndex]
            else:
                self.natural_abundance = []
                self.natural_abundance_uncert = []
        else:
            ground_data = self.nubase_data.loc[self.index]
            self.ground_decay_modes = [ground_data['Decay mode']]
            self.decay_modes = [val for val in self.isomer_states['Decay mode 2']]
            self.ground_states = self.nubase_data.loc[self.index:self.index+self.ground_level_count]
            self.ground_branch_ratio = [val for val in self.ground_states['Branch ratio/Isoptope abund']]
            self.branch_ratio = [val for val in self.isomer_states['Branch ratio 2']]
            self.ground_branch_ratio_uncert = [val for val in self.ground_states['Branch ratio/Isoptope abund uncert']]
            self.branch_ratio_uncert = [val for val in self.isomer_states['Branch ratio uncert 2']]
            if 'is' in self.ground_decay_modes:
                isIndex = self.ground_decay_modes.index('is')
                self.natural_abundance = self.ground_branch_ratio[isIndex]
                self.natural_abundance_uncert = self.ground_branch_ratio_uncert[isIndex]
                del self.ground_decay_modes[isIndex]
                del self.ground_branch_ratio[isIndex]
                del self.ground_branch_ratio_uncert[isIndex]
            else:
                self.natural_abundance = []
                self.natural_abundance_uncert = []
        if self.level == 0:
            self.decayModes = self.ground_decay_modes
        else:
            self.decayModes = self.decay_modes

        def fix_br(val):
            '''
            Fixes Branch ratio value by first checking if it is a ?, if so then set value to nan
            and flag as energetically allowed. If it isn't a ?, then check if it contains a #. If so
            then remove that symbol and flag as an estimate based on surrounding N and Z numbers.

            @param: val, str

            @returns:   val, float
                        E_flag, True/False
                        estimate_flag, True/False
            '''
            if val == '?':
                val = np.nan
                E_flag = True
                estimate_flag = True
            else:
                E_flag = False
                val_split = split(val)
                if '#' in val_split:
                    val_split.remove('#')
                    estimate_flag = True
                    val = ''
                    val = float(val.join(val_split))
                else:
                    val = float(val)
                    estimate_flag = False
                    E_flag = False
            return val, E_flag, estimate_flag

        self.excit_energy_allowed, self.excit_estimate_flags = [np.zeros(len(self.branch_ratio)).tolist() for i in range(2)]
        self.ground_energy_allowed, self.ground_estimate_flags = [np.zeros(len(self.ground_branch_ratio)).tolist() for i in range(2)]
        for i,val in enumerate(self.branch_ratio):
            self.branch_ratio[i], self.excit_energy_allowed[i], self.excit_estimate_flags[i] = fix_br(str(val))
        for i,val in enumerate(self.ground_branch_ratio):
            self.ground_branch_ratio[i],self.ground_energy_allowed[i],self.ground_estimate_flags[i] = fix_br(str(val))
        if self.level == 0:
            self.branchRatio = self.ground_branch_ratio
        else:
            self.branchRatio = self.branch_ratio
        return self

    def decay_constant(self):
        '''
        Simply calculates the decay constant from the half life.

        @param: self

        generates:  self.decay_const - decay constant
                    self.decay_const_uncert - uncertainty in decay constant

        '''
        ln2 = np.log(2)
        if self.level > 0:
            self.decay_const = ln2/self.excit_HL
            self.decay_const_uncert = self.decay_const*(self.excit_HL_uncert/self.excit_HL)
        else:
            if self.HL != 'stable' and math.isinf(self.HL) == False:
                self.decay_const = ln2/self.HL
                self.decay_const_uncert = self.decay_const*(self.HL/self.HL_uncert)
            else:
                self.decay_const = 0
                self.decay_const_uncert = 0
        return self

    def daughter(self):
        '''
        Finds daughter isomers from decay modes. Daughter isotopes are in ground state. Not true for real events.

        @param: self
        
        generates: self.daughters, list - saves atomic data in NU format
        '''
        if self.decay_const != 0:
            if self.level == 0:
                decayModes = self.ground_decay_modes
            else:
                decayModes = self.decay_modes
            self.daughters = np.zeros(len(decayModes)).tolist()
            for i,decay in enumerate(decayModes):
                if decay == 'b-':
                    #Beta minus decay. Emits antineutrino
                    self.daughters[i] = self.NU+10000000
                elif decay == '2b-':
                    #2 Beta minus decay. Emits 2 antineutrinos
                    self.daughters[i] = self.NU+20000000
                elif decay == 'it':
                    #isomeric transition. Emits gamma radiation
                    self.daughters[i] = self.NU-self.NU%10000
                elif decay == 'b+':
                    #Beta plus decay. Emits neturino
                    self.daughters[i] = self.NU-10000000
                elif decay == '2b+':
                    #Beta plus decay. Emits 2 neutrinos
                    self.daughters[i] = self.NU-20000000
                elif decay == 'ec':
                    #Electron capture. Emits neutrino, daughter in exited state
                    self.daughters[i] = self.NU-10000000
                elif decay == 'a':
                    #alpha decay
                    self.daughters[i] = self.NU-20040000
                elif decay =='b-n':
                    #beta and neutron emission
                    self.daughters[i] = self.NU+10000000-10000
                elif decay == 'p':
                    #proton emission
                    self.daughters[i] = self.NU-10010000
                elif decay =='2p':
                    #2 proton emission
                    self.daughters[i] = self.NU-20020000
                elif decay == 'n':
                    #neutron emission
                    self.daughters[i] = self.NU-10000
                elif decay == '2n':
                    #2 neutron emission
                    self.daughters[i] = self.NU-20000
                elif decay == 'b-2n':
                    #beta and 2 neutron emission
                    self.daughters[i] = self.NU+10000000-20000
                elif decay == 'b+n':
                    self.daughters[i] = self.NU-10000000-10000
                elif decay == 'b+2n':
                    self.daughters[i] = self.NU-10000000-20000
                if preferences.set_ground_state == True:
                    self.daughters[i] = self.daughters[i]-self.daughters[i]%10000
        else:
            self.daughters = []
        return self