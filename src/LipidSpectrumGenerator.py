import tkinter as tk
from itertools import combinations_with_replacement as cwr
from itertools import product
from itertools import compress
import math
import time

Charge_Options = {'POS':                     False,     
                  'NEG':                     True}

Lipid_Options = {'Plasmenyl':                True,
                'Odd_Chains':                True,
                'Variable_SB_Length':        False} # SB = Sphingoid base

lipid_definitions = {'MAG':                 [False, [1, 'Glycerol', 'MAG'],          ["[M+H]+", "[M+H-H2O]+", "[M+Na]+", "[M+NH4]+"]], # [Number of lipid tails, Lipid body, Headgroup/Type]
                     'DAG':                 [False, [2, 'Glycerol', 'DAG'],          ["[M+H]+", "[M+H-H2O]+", "[M+Na]+", "[M+NH4]+"]],           
                     'TAG':                 [False, [3, 'Glycerol', 'TAG'],          ["[M+Na]+", "[M+NH4]+"]],

                     'PA':                  [False, [2, 'Glycerol', 'PA'],           ["[M-H]-", "[M+H]+", "[M+Na]+"]],
                     'Lyso PA':             [False, [1, 'Glycerol', 'PA'],           ["[M-H]-", "[M+H]+", "[M+Na]+"]],
                     'PC':                  [False, [2, 'Glycerol', 'PC'],           ["[M+H]+", "[M+Na]+"]],
                     'Lyso PC':             [False, [1, 'Glycerol', 'PC'],           ["[M+H]+", "[M+Na]+"]],
                     'PE':                  [False, [2, 'Glycerol', 'PE'],           ["[M-H]-", "[M+H]+", "[M+Na]+"]],
                     'Lyso PE':             [False, [1, 'Glycerol', 'PE'],           ["[M-H]-", "[M+H]+", "[M+Na]+"]],
                     'PG':                  [False, [2, 'Glycerol', 'PG'],           ["[M-H]-", "[M+H]+", "[M+Na]+"]],
                     'Lyso PG':             [False, [1, 'Glycerol', 'PG'],           ["[M-H]-", "[M+H]+", "[M+Na]+"]],
                     'PI':                  [True,  [2, 'Glycerol', 'PI'],           ["[M-H]-"]],
                     'Lyso PI':             [True,  [1, 'Glycerol', 'PI'],           ["[M-H]-"]],
                     'PIP':                 [False, [2, 'Glycerol', 'PIP'],          []],
                     'PIP2':                [False, [2, 'Glycerol', 'PIP2'],         []],
                     'PS':                  [False, [2, 'Glycerol', 'PS'],           ["[M-H]-", "[M+H]+", "[M+Na]+"]],
                     'Lyso PS':             [False, [1, 'Glycerol', 'PS'],           ["[M-H]-", "[M+H]+", "[M+Na]+"]],

                     'MGDG':                [False, [2, 'Glycerol', 'MGDG'],         ["[M+Na]+", "[M+NH4]+"]],
                     'DGDG':                [False, [2, 'Glycerol', 'DGDG'],         ["[M+Na]+", "[M+NH4]+"]],
                     
                     'Sphingosine':         [False, [1, 'SB', 'Sphingosine'],        ["[M-H]-", "[M+Cl]-", "[M+H]+", "[M+Na]+"]], # LipidBlast seemed to have only "[M+Na]+" spectra. Seen "[M+NH4]+" and "[M-H]-" in some articles.
                     'N-acylsphingosine':   [False, [2, 'SB', 'N-acylsphingosine'],  ["[M-H]-", "[M+Cl]-", "[M+H]+", "[M+Na]+"]], # Lipidblast also only has 1 fragment for these, and there aren't many characteristic studies.
                     
                     'Ceramide PA':         [False, [2, 'SB', 'PA'],                 ["[M-H]-", "[M+Cl]-", "[M+H]+", "[M+Na]+"]],
                     'Ceramide PC':         [False, [2, 'SB', 'PC'],                 ["[M-H]-", "[M+Cl]-", "[M+H]+", "[M+Na]+"]],
                     'Ceramide PE':         [False, [2, 'SB', 'PE'],                 ["[M-H]-", "[M+Cl]-", "[M+H]+", "[M+Na]+"]],
                     'Ceramide PG':         [False, [2, 'SB', 'PG'],                 ["[M-H]-", "[M+Cl]-", "[M+H]+", "[M+Na]+"]],
                     'Ceramide PI':         [False, [2, 'SB', 'PI'],                 ["[M-H]-", "[M+Cl]-", "[M+H]+", "[M+Na]+"]],
                     'Ceramide PS':         [False, [2, 'SB', 'PS'],                 ["[M-H]-", "[M+Cl]-", "[M+H]+", "[M+Na]+"]],
                     'Ceramide MG':         [False, [2, 'SB', 'MGDG'],               ["[M-H]-", "[M+Cl]-", "[M+H]+", "[M+Na]+"]],
                     'Ceramide DG':         [False, [2, 'SB', 'DGDG'],               ["[M-H]-", "[M+Cl]-", "[M+H]+", "[M+Na]+"]]}     

Adduct_Masses = {"[M-H]-":                  [-1.007276,'N'],
                 "[M+Cl]-":                 [34.969401, 'N'],
                 "[M+H]+":                  [1.007276,'P'],
                 "[M+H-H2O]+":              [-17.003289, 'P'],
                 "[M+Na]+":                 [22.989221, 'P'],
                 "[M+NH4]+":                [18.033826, 'P']}

def calculate_peaks(species, tails, Adduct, lipids):

    def new_peak_list(Name, PorN, PrecursorMZ, Adduct, Peak_Count, Peaks):

        a = {"Name:":Name,
            "Spectrum_type:":"MS2",
            "Ion_mode:":PorN, # Positive or Negative
            "PrecursorMZ:":f"{PrecursorMZ:.6f}",
            "Precursor_type:":Adduct,
            "Num Peaks:":Peak_Count,
            "peaks":Peaks} # Fragmentation spectra here as [[mz, intensity],[...]]
        return a
     
    name =  []
    Peaks = [] # This is where the fragments will be stored
    body = lipid_definitions[species][1][1]
    head = lipid_definitions[species][1][2]

    Masses = {'H':            1.007276,
              'H2O':          18.010565,
              'Glycerol':     92.047344}

    headgroup_mass = {'PA':   97.976895,  # Masses are headgroup + phosphate
                      'PC':   183.066044, # PC has -H to maintain neutral charge
                      'PE':   141.019094,
                      'PG':   172.013674,
                      'PI':   260.029718,
                      'PS':   185.008923,
                      'MGDG': 180.063388, # Masses are just of the galactose ring
                      'DGDG': 342.116212}

    if body == 'Glycerol':
        Mass = Masses['Glycerol']
        for tail in tails:
            Mass += (tail[1] - Masses['H2O'])
            name.append(tail[0])
        if len(tails) < 2:
            name.append('0:0')  # Adds 0:0 to the name for Lyso GPLs 
        if species in ('MAG', 'DAG'):
            name.append('0:0')  # Adds another 0:0 to the name for MAGs and DAGs
        try:
            Mass += (headgroup_mass[head] - Masses['H2O'])
        except: pass

    if body == 'SB':
        Mass = tails[0][1]
        name.append(tails[0][0])
        try: 
            Mass += tails[1][1] - Masses['H2O']
            name.append(tails[1][0])
        except: pass
        try:
            Mass += (headgroup_mass[head] - Masses['H2O'])
        except: pass

    nl_Headgroup = {'PA':97.976895, # These are separate because I'm too tired to think of a clever way to derive them from the headgroup masses
                    'PI':180.063388 - Masses['H2O'],
                    'PE':0,
                    'PS':89.047678 - 2*Masses['H'],
                    'PG':92.047344 - Masses['H2O']} # Neutral loss masses

    gpl_characteristic = {'PA':{'N':[[78.959053, 5], [96.969618, 5],[152.995833, 10]],                                                                                            'P':[[Mass - headgroup_mass['PA'] + Masses['H'], 100]]},
                          'PC':{'N':[[78.959053, 5], [96.969618, 5],[152.995833, 10]],                                                                                            'P':[[headgroup_mass['PC'] + Masses['H'], 100]]},
                          'PE':{'N':[[78.959053, 5], [96.969618, 5],[152.995833, 10], [140.011817, 3], [196.038032, 5]],                                                          'P':[[Mass - headgroup_mass['PE'] + Masses['H'], 100]]},
                          'PG':{'N':[[78.959053, 5], [96.969618, 5],[152.995833, 10], [209.022048, 5], [227.032612, 5], [245.043177, 5]],                                         'P':[[Mass - headgroup_mass['PG'] + Masses['H'], 100]]},
                          'PI':{'N':[[78.959053, 5], [96.969618, 5],[152.995833, 10], [223.001312, 5], [241.011877, 25], [259.022442, 5], [297.038092, 5], [315.048656, 5]],      'P':[]},
                          'PS':{'N':[[78.959053, 5], [96.969618, 5],[152.995833, 10]],                                                                                            'P':[[Mass - headgroup_mass['PS'] + Masses['H'], 100]]},
                          'MGDG':{'N':[],                                                                                                                                         'P':[[Mass - headgroup_mass['MGDG'] + Masses['H2O'] + Masses['H'], 50]]},
                          'DGDG':{'N':[],                                                                                                                                         'P':[[Mass - headgroup_mass['DGDG'] + Masses['H2O'] + Masses['H'], 50]]}}

    spl_characteristic = {'PA':{'N':[],                                                                                          'P':[]},
                          'PC':{'N':[],                                                                                          'P':[[headgroup_mass['PC'] + Masses['H'], 100], [Mass - 59.073499 + Masses['H'], 25], [Mass - 59.073499 - Masses['H2O'] + Masses['H'], 25], [Mass - 183.066044 + Masses['H'], 25]]},
                          'PE':{'N':[[140.011817, 3], [196.038032, 5]],                                                          'P':[]},
                          'PG':{'N':[[209.022048, 5], [227.032612, 5], [245.043177, 5]],                                         'P':[]},
                          'PI':{'N':[[223.001312, 5], [241.011877, 25], [259.022442, 5], [297.038092, 5], [315.048656, 5]],      'P':[]},
                          'PS':{'N':[],                                                                                          'P':[]},
                          'MGDG':{'N':[],                                                                                        'P':[[Mass - headgroup_mass['MGDG'] + Masses['H2O'] + Masses['H'], 50]]},
                          'DGDG':{'N':[],                                                                                        'P':[[Mass - headgroup_mass['DGDG'] + Masses['H2O'] + Masses['H'], 50]]}}

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    if (body == ('Glycerol')) and (head in ('MAG', 'DAG', 'TAG')): ##### J Am Soc Mass Spectrom. 2010 April ; 21(4): 657–669. https://doi.org/10.1016/j.jasms.2010.01.007
                                                                   ##### Metabolites 2016, 6(3), 25; https://doi.org/10.3390/metabo6030025
        i = {'MAG':[100, 50,  50, 50, 0,  10, 0,  100],
             'DAG':[10,  100, 50, 1,  25,  5, 10, 25],
             'TAG':[0,   0,   0,  10, 100, 5, 10, 25]}

        if Adduct in ("[M+NH4]+", "[M+H]+"): Peaks.extend([[Mass + Adduct_Masses["[M+H]+"][0], i[head][0]], [Mass + Adduct_Masses["[M+H-H2O]+"][0], i[head][1]]]) # Lipidblast includes "[M+H]+", "[M+H-H2O]+" parents in the "[M+NH4]" spectra.
        if Adduct in ("[M+Na]+"): Peaks.append([Mass - Masses['H2O'] + Adduct_Masses[Adduct][0], i[head][2]])
        for tail in tails: 
            Peaks.extend([[(tail[1] - Masses['H2O'] + Masses['H']), i[head][3]], [(Mass - tail[1] + Masses['H']), i[head][4]]])
            if Adduct == "[M+Na]+":
                Peaks.extend([[(tail[1] - Masses['H2O'] + Adduct_Masses[Adduct][0]), i[head][5]], [(Mass - tail[1] + Adduct_Masses[Adduct][0]), i[head][6]]])        
        Peaks.append([Mass + Adduct_Masses[Adduct][0], i[head][7]]) # Add [M+Adduct]+ to the spectra

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    if (body == ('Glycerol')) and (head in ('PA', 'PC', 'PE', 'PG', 'PI', 'PIP', 'PIP2', 'PS', 'MGDG', 'DGDG')): ##### J Chromatogr B Analyt Technol Biomed Life Sci. 2009 September 15; 877(26): 2673–2695. doi:10.1016/j.jchromb.2009.02.033.

        i = {'PA':  [100, 15, 5, 15, 5, 25,    8, 0,  0,   0, 0, 0, 0, 0, 4, 25],
             'PC':  [100, 15, 5, 15, 5, 25,    0, 0,  0,   0, 0, 8, 4, 8, 4, 25],
             'PE':  [100, 15, 5, 15, 5, 25,    8, 0,  0,   0, 0, 0, 0, 0, 4, 25],
             'PG':  [100, 15, 5, 15, 5, 25,    8, 0,  0,   0, 0, 0, 0, 0, 4, 25],
             'PI':  [100, 15, 5, 15, 5, 25,    0, 0,  0,   0, 0, 0, 0, 0, 4, 25],
             'PIP': [100, 15, 5, 15, 5, 25,    0, 0,  0,   0, 0, 0, 0, 0, 4, 25],
             'PIP2':[100, 15, 5, 15, 5, 25,    0, 0,  0,   0, 0, 0, 0, 0, 4, 25],
             'PS':  [100, 15, 5, 15, 5, 25,    8, 0,  0,   0, 0, 0, 0, 0, 4, 25],
             'MGDG':[100, 15, 5, 15, 5, 25,    0, 50, 100, 4, 8, 0, 0, 0, 4, 25],
             'DGDG':[100, 15, 5, 15, 5, 25,    0, 50, 100, 4, 8, 0, 0, 0, 4, 25]}

        if Adduct_Masses[Adduct][1] == 'N' and Charge_Options['NEG'] == True: # Negative spectra
            for fragment in gpl_characteristic[head]['N']: Peaks.append(fragment)
            for tail in tails: 
                if 'O-' not in tail[0]: Peaks.extend([[tail[1] - Masses['H'], i[head][0]], [Mass - tail[1] - Masses['H'], i[head][1]], [Mass - tail[1] + Masses['H2O'] - Masses['H'], i[head][2]], [Mass - nl_Headgroup[head] - tail[1] - Masses['H'], i[head][3]], [Mass - nl_Headgroup[head] - tail[1] + Masses['H2O'] - Masses['H'], i[head][4]], [Mass + Adduct_Masses[Adduct][0], i[head][5]]])
        elif Adduct_Masses[Adduct][1] == 'P' and Charge_Options['POS'] == True: # Positive Spectra
            if Adduct == "[M+Na]+": counter_ion = Adduct_Masses[Adduct][0]
            else: counter_ion = Masses['H']
            for fragment in gpl_characteristic[head]['P']: Peaks.append(fragment)  
            for tail in tails:
                Peaks.extend([[gpl_characteristic[head]['P'][0][0] - tail[1] + (Masses['H2O'] - Masses['H']) + counter_ion, i[head][6]], [gpl_characteristic[head]['P'][0][0] - Masses['H2O'], i[head][7]], [gpl_characteristic[head]['P'][0][0] - tail[1] - Masses['H'] + counter_ion, i[head][8]], [Mass - tail[1] + Masses['H2O'] + counter_ion, i[head][9]], [Mass - tail[1] + counter_ion, i[head][10]]])
                if Adduct == "[M+Na]+": Peaks.append([Mass - 59.073499 + counter_ion, i[head][11]]) #(Loss of N(CH3)3) Only happens for Alkaline salts?
                Peaks.extend([[Mass - tail[1] + Masses['H2O'] + counter_ion, i[head][12]], [Mass - tail[1] + counter_ion, i[head][13]], [tail[1] - Masses['H2O'] + counter_ion, i[head][14]]]) # For every tail, add a fragment corresponding to [RC=O]+ to the spectra
            Peaks.append([Mass + Adduct_Masses[Adduct][0], i[head][15]]) 

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    if (body == ('SB')) and (head in ('Sphingosine', 'N-acylsphingosine')):
        if Adduct_Masses[Adduct][1] == 'N' and Charge_Options['NEG'] == True: # Negative spectra
            Peaks.extend([[Mass - 30.010565 - Masses['H'], 25], [Mass - 32.026215 - Masses['H'], 12], [Mass - 30.010565 - Masses['H2O'] - Masses['H'], 12]])
            try: Peaks.extend([[tails[1][1] - Masses['H'], 12], [Mass - (tails[0][1] - 61.052764) - Masses['H2O'] - Masses['H'], 100], [Mass - (tails[0][1] - 61.052764) - 3*Masses['H'], 12]])
            except: pass
        elif Adduct_Masses[Adduct][1] == 'P' and Charge_Options['POS'] == True: # Positive Spectra
            Peaks.extend([[Mass - Masses['H2O'] + Adduct_Masses[Adduct][0], 60], [Mass - 2*Masses['H2O'] + Adduct_Masses[Adduct][0], 25]])
            try: Peaks.extend([[Mass - tails[1][1] + Adduct_Masses[Adduct][0], 25], [Mass - tails[1][1] - Masses['H2O'] + Adduct_Masses[Adduct][0], 25]])
            except: pass
        Peaks.append([Mass + Adduct_Masses[Adduct][0], 25])

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    if (body == ('SB')) and (head in ('PA', 'PC', 'PE', 'PG', 'PI', 'PIP', 'PIP2', 'PS', 'MGDG', 'DGDG')):
        if Adduct_Masses[Adduct][1] == 'N' and Charge_Options['NEG'] == True: # Negative spectra
            for fragment in spl_characteristic[head]['N']: Peaks.append(fragment)
            Peaks.extend([[Mass - 30.010565 - Masses['H'], 25], [Mass - 32.026215 - Masses['H'], 12], [Mass - 30.010565 - Masses['H2O'] - Masses['H'], 12]])
            try: Peaks.extend([[tails[1][1] - Masses['H'], 12], [Mass - (tails[0][1] - 61.052764) - Masses['H2O'] - Masses['H'], 100], [Mass - (tails[0][1] - 61.052764) - 3*Masses['H'], 12]])   
            except: pass
        if Adduct == "[M+Na]+": # This is to account for the extra mass when the +'ve charge comes from Na+, whereas NH4, H, H-H2O all give H+
            counter_ion = Adduct_Masses[Adduct][0]
        else: counter_ion = Masses['H']
        if Adduct_Masses[Adduct][1] == 'P' and Charge_Options['POS'] == True: # Positive Spectra ##### 10.1194/jlr.M067199
            for fragment in spl_characteristic[head]['P']: Peaks.append(fragment)
            Peaks.extend([[Mass - Masses['H2O'] + counter_ion, 60], [Mass - 2*Masses['H2O'] + counter_ion, 10], [Mass - headgroup_mass[head] + counter_ion, 10]])
            try: Peaks.extend([[Mass - tails[1][1] - Masses['H'] + counter_ion, 20], [Mass - headgroup_mass[head] - tails[1][1] + counter_ion, 20]])
            except: pass
        Peaks.append([Mass + Adduct_Masses[Adduct][0], 25])

    # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

    Peaks.sort()
    new_Peaks = [] # Work around to remove duplicates from fragment list  
    max_intensity = 100 #max([peak[-1] for peak in Peaks])
    for peak in Peaks:
        if (f"{peak[0]:.6f}" not in (frag[0] for frag in new_Peaks)) and (peak[1] >= 1): # Warning, do not look at this.
            new_Peaks.append([f"{peak[0]:.6f}", 100*(peak[1]/max_intensity)])
    Peaks = new_Peaks

    new_list = new_peak_list(f"{species} {'_'.join(str(nam) for nam in name)} {Adduct}", Adduct_Masses[Adduct][1], Mass + Adduct_Masses[Adduct][0], Adduct, len(Peaks), Peaks)
    lipids[new_list["Name:"]] = new_list
    return

def new_spectrum_library():

    def calculate_glycero(tail):

        if tail[0] >= 2:
            Name = f"{tail[0]}:{tail[1]}"
            Mass = (31.98982926 + tail[0]*14.01565007 - tail[1]*2.01565007)
        else: pass
        return [Name, Mass]

    def calculate_glycero_plasmenyl(tail):

        if tail[0] >= 2:
            Name = f"O-{tail[0]}:{tail[1]}"
            Mass = (15.994915 + tail[0]*14.01565007 + 2.01565007 - tail[1]*2.01565007)
        else: pass
        return [Name, Mass]

    def calculate_sphingoid(tail):

        tails = []
        if Lipid_Options['Variable_SB_Length'] == True:
            c_d = tail[0]
        else: c_d = [18, 1]

        if c_d[0] >= 3:
            Name = f"Sp_{c_d[0]}:{c_d[1]}"
            Mass = (91.063329 + (c_d[0] - 3)*14.01565007 - c_d[1]*2.01565007)
            if [Name, Mass] not in tails: tails.append([Name, Mass])
            else: pass
        else: pass

        try: 
            c_d = tail[1]
            if c_d[0] >= 2:
                Name = f"N_{c_d[0]}:{c_d[1]}"
                Mass = (31.98982926 + c_d[0]*14.01565007 - c_d[1]*2.01565007)
                tails.append([Name, Mass])
            else: pass
        except: pass

        return tails

    Cmin = 12   # Min and Max number of Carbons
    Cmax = 20 + 1
    Dmin = 0    # Min and Max number of Desaturated bonds
    Dmax = 6 + 1

    glycero_tails = []
    sb_tails = []
    lipids = {}

    for C in range(Cmin, Cmax): # Limits all fatty acids within the ranges set
        for D in range(Dmin, Dmax):
            if D <= ((C-1)/2): # Limits the amount of desaturation a fatty acid can have
                if (Lipid_Options['Odd_Chains'] == False) and (C % 2 != 0):
                    continue
                else:
                    glycero_tails.append(calculate_glycero([C, D]))
                    sb_tails.append([C, D])
                    if Lipid_Options['Plasmenyl'] == True: glycero_tails.append(calculate_glycero_plasmenyl([C, D]))
                    else: pass
            else: continue

    for species in list(compress([lipid for lipid in lipid_definitions],[int(lipid_definitions[lipid][0]) for lipid in lipid_definitions])):
        body = lipid_definitions[species][1][1]
        for Adduct in lipid_definitions[species][2]:
            if body == 'Glycerol':
                for comb in cwr(glycero_tails, lipid_definitions[species][1][0]):
                    calculate_peaks(species, comb, Adduct, lipids)
            elif body == 'SB':
                for comb in product(sb_tails, repeat = lipid_definitions[species][1][0]):
                    calculate_peaks(species, calculate_sphingoid(comb), Adduct, lipids)

    return lipids

def buttonpush():

    def write_to_file(spectrum_library):

        if Charge_Options['POS'] == True: Poslibrary = open("Positive_Library.msp", "a")
        if Charge_Options['NEG'] == True: Neglibrary = open("Negative_Library.msp", "a")

        count = 0
        for lipid in spectrum_library: # For each lipid
            if spectrum_library[lipid]['Ion_mode:'] == 'P' and Charge_Options['POS'] == True:
                for line in spectrum_library[lipid]: # Write the lines
                    if line == "peaks": Poslibrary.writelines([f"{peak[0]} {peak[1]}\n" for peak in spectrum_library[lipid][line]]) # if spectra, whole thing can be written at once
                    else: Poslibrary.write(f"{line} {spectrum_library[lipid][line]}\n") # Information is written one line at a time to look for where the spectra is  
                Poslibrary.write("\n") # After it finishes a spectrum, new line ready for the next
                count += 1

            if spectrum_library[lipid]['Ion_mode:'] == 'N' and Charge_Options['NEG'] == True:
                for line in spectrum_library[lipid]: # Write the lines
                    if line == "peaks": Neglibrary.writelines([f"{peak[0]} {peak[1]}\n" for peak in spectrum_library[lipid][line]]) # if spectra, whole thing can be written at once
                    else: Neglibrary.write(f"{line} {spectrum_library[lipid][line]}\n") # Information is written one line at a time to look for where the spectra is      
                Neglibrary.write("\n") # After it finishes a spectrum, new line ready for the next
                count += 1
        
        try: Poslibrary.close()
        except: pass
        try: Neglibrary.close()
        except: pass

        return count

    t0 = time.time()
    spectrum_library = new_spectrum_library()
    print("Writing...")
    #for lipid in spectrum_library:
    count = write_to_file(spectrum_library)      
    t1 = time.time()
    print(f"Produced {count} spectra in {t1 - t0:.4f} seconds")   

root = tk.Tk()
button = tk.Button(root, text = "Push me", command = buttonpush)
button.pack()
tk.mainloop()