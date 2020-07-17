import tkinter as tk
from itertools import combinations_with_replacement as cwr
from itertools import compress
import time

Charge_Options = {'POS':True,
                  'NEG':False}
Glycerolipids = {'MAG':False, 
                 'DAG':False,
                 'TAG':False}
Glycerophospholipids = {'PA':False,
                        'PC':False,
                        'PE':False,
                        'PG':False,
                        'PI':False,
                        'PIP':False,
                        'PIP2':False,
                        'PS':False}
Glycerophospholipids_Options = {'Lyso':False,
                                'Diacyl':False}
FA_Options = {'Plasmenyl':False,
              'Odd_Chains':False}
Sphingolipids = {'Sphingosine':False,
                 'N-acylsphingosine':True}
Phosphosphingolipids = {'PA':False,
                        'PC':False,
                        'PE':False,
                        'PG':False,
                        'PI':False,
                        'PIP':False,
                        'PIP2':False,
                        'PS':False}
Sphingo_Options = {'Variable_Sphingosine_Length':False}
Chain_Number = {'MAG':1,
                'DAG':2,
                'TAG':3,
                'Sphingosine':1,
                'N-acylsphingosine':2}
Adduct_Definitions = {'MAG':["[M+H]+", "[M+H-H2O]+", "[M+Na]+", "[M+NH4]+"],
                      'DAG':["[M+H]+", "[M+H-H2O]+", "[M+Na]+", "[M+NH4]+"],
                      'TAG':["[M+Na]+", "[M+NH4]+"],
                      'PA':["[M-H]-", "[M+H]+", "[M+Na]+"],
                      'PC':["[M+H]+", "[M+Na]+"],
                      'PE':["[M-H]-", "[M+H]+", "[M+Na]+"],
                      'PG':["[M-H]-", "[M+H]+", "[M+Na]+"],
                      'PI':["[M-H]-"],
                      'PIP':[],
                      'PIP2':[],
                      'PS':["[M-H]-", "[M+H]+", "[M+Na]+"],
                      'Sphingosine':["[M-H]-", "[M+Cl]-", "[M+H]+", "[M+Na]+"],
                      'N-acylsphingosine':["[M-H]-", "[M+Cl]-", "[M+H]+", "[M+Na]+"]}
Adduct_Masses = {"[M-H]-":[-1.007276,'N'],
                 "[M+Cl]-":[34.969401, 'N'],
                 "[M+H]+":[1.007276,'P'],
                 "[M+H-H2O]+":[-17.003289, 'P'],
                 "[M+Na]+":[22.989221, 'P'],
                 "[M+NH4]+":[18.033826, 'P']}

def new_peak_list(Name, PorN, PrecursorMZ, Adduct, Peak_Count, Peaks):

    a = {"Name:":Name,
        "Spectrum_type:":"MS2",
        "Ion_mode:":PorN, # Positive or Negative
        "PrecursorMZ:":f"{PrecursorMZ:.6f}",
        "Precursor_type:":Adduct,
        "Num Peaks:":Peak_Count,
        "peaks":Peaks} # Fragmentation spectra here as [[mz, intensity],[...]]
    return a

def calculate_peaks(path, species, tails, Adduct, lipids):

    name = []
    Peaks = [] # This is where the fragments will be stored
    Masses = {'H':1.007276,
              'H2O':18.010565,
              'Glycerol':92.047344,
              'Sphingosine':299.282429} # Spicy numbers
    headgroup_mass = {'PA':97.976895, # Masses are headgroup + phosphate
                      'PC':183.066044, # PC has -H to maintain neutral charge
                      'PE':141.019094,
                      'PG':172.013674,
                      'PI':260.029718,
                      'PS':185.008923}
    nl_Headgroup = {'PA':97.976895, # These are separate because I'm too tired to think of a clever way to derive them from the headgroup masses
                    'PI':180.063388 - Masses['H2O'],
                    'PE':0,
                    'PS':89.047678 - 2*Masses['H'],
                    'PG':92.047344 - Masses['H2O']} # Neutral loss masses

    def characteristic(Mass):
        characteristic = {'PA':{'N':[],                                                                                            'P':[[Mass - headgroup_mass['PA'] + Masses['H'], 100]]},
                          'PC':{'N':[],                                                                                            'P':[[headgroup_mass['PC'] + Masses['H'], 100]]},
                          'PE':{'N':[[140.011817, 3], [196.038032, 5]],                                                            'P':[[Mass - headgroup_mass['PE'] + Masses['H'], 100]]},
                          'PG':{'N':[[209.022048, 5], [227.032612, 5], [245.043177, 5]],                                           'P':[[Mass - headgroup_mass['PG'] + Masses['H'], 100]]},
                          'PI':{'N':[[223.001312, 5], [241.011877, 25], [259.022442, 5], [297.038092, 5], [315.048656, 5]],        'P':[]},
                          'PS':{'N':[],                                                                                            'P':[[Mass - headgroup_mass['PS'] + Masses['H'], 100]]}}
        return characteristic
        
    if path == 1 and Charge_Options['POS'] == True: # For MAG, DAG or TAG. Only positive spectra generated, so no negative path. ##### J Am Soc Mass Spectrom. 2010 April ; 21(4): 657–669. doi:10.1016/j.jasms.2010.01.007
        Mass = Masses['Glycerol'] # Begin with a glycerol backbone
        for tail in tails: # For every tail being added to the backbone, mass is increased by the mass of the tail and a water is removed to account for the bonding
            Mass += (tail[1] - Masses['H2O'])
            name.append(tail[0])
        if Adduct == "[M+Na]+": # Na+ spectra are slightly different, they have peaks for both [DAG+Na+] and [DAG+H+]
            for tail in tails: 
                Peaks.append([(tail[1] - Masses['H2O'] + Adduct_Masses[Adduct][0]), 5]) # For every tail, add a fragment corresponding to [RC=O]+ to the spectra
                Peaks.append([(Mass - tail[1] + Masses['H']), 40])
                Peaks.append([(Mass - tail[1] + Adduct_Masses[Adduct][0]), 100]) # For every tail, add a fragment corresponding to the neutral loss of RCOOH to the spectra
        else:
            for tail in tails: 
                Peaks.append([(tail[1] - Masses['H2O'] + Masses['H']), 5]) # For every tail, add a fragment corresponding to [RC=O]+ to the spectra
                Peaks.append([(Mass - tail[1] + Masses['H']), 100]) # For every tail, add a fragment corresponding to the neutral loss of RCOOH to the spectra
        Peaks.append([Mass + Adduct_Masses[Adduct][0], 25]) # Add [M+Adduct]+ to the spectra

    if path == 2: # For PA, PC, PE, PG, PI & PS GPLs ##### J Chromatogr B Analyt Technol Biomed Life Sci. 2009 September 15; 877(26): 2673–2695. doi:10.1016/j.jchromb.2009.02.033.
        Mass = Masses['Glycerol']
        for tail in tails: # For every tail being added to the backbone, mass is increased by the mass of the tail and a water is removed to account for the bonding
            Mass += (tail[1] - Masses['H2O'])
            name.append(tail[0])
        if len(tails) < 2:
            name.append('0:0')  # Adds 0:0 to the name for Lyso GPLs 
        Mass += (headgroup_mass[species] - Masses['H2O'])
        if Adduct_Masses[Adduct][1] == 'N' and Charge_Options['NEG'] == True: # Negative spectra
            Peaks.append([78.959053, 5])    # PO3-
            Peaks.append([96.969618, 5])    # H2PO4-
            Peaks.append([152.995833, 10])  # Glycerol-3-phosphate, -H2O
            for fragment in characteristic(Mass)[species]['N']: Peaks.append(fragment)
            for tail in tails:
                Peaks.append([tail[1] - Masses['H'], 100])
                Peaks.append([Mass - tail[1] - Masses['H'], 15])
                Peaks.append([Mass - tail[1] + Masses['H2O'] - Masses['H'], 5])
                Peaks.append([Mass - nl_Headgroup[species] - tail[1] - Masses['H'], 15])
                Peaks.append([Mass - nl_Headgroup[species] - tail[1] + Masses['H2O'] - Masses['H'], 5])
            Peaks.append([Mass + Adduct_Masses[Adduct][0], 25]) 
        elif Adduct_Masses[Adduct][1] == 'P' and Charge_Options['POS'] == True: # Positive Spectra
            if Adduct == "[M+Na]+": # This is to account for the extra mass when the +'ve charge comes from Na+, whereas NH4, H, H-H2O all give H+
                counter_ion = Adduct_Masses[Adduct][0]
            else: counter_ion = Masses['H']
            for fragment in characteristic(Mass)[species]['P']: Peaks.append(fragment)  
            for tail in tails:
                if species != 'PC': 
                    Peaks.append([characteristic(Mass)[species]['P'][0][0] - tail[1] + (Masses['H2O'] - Masses['H']) + counter_ion, 8])
                if species == 'PC': 
                    Peaks.append([Mass - tail[1] + Masses['H2O'] + counter_ion, 4])
                    Peaks.append([Mass - tail[1] + counter_ion, 8])
                Peaks.append([tail[1] - Masses['H2O'] + counter_ion, 4]) # For every tail, add a fragment corresponding to [RC=O]+ to the spectra
            Peaks.append([Mass + Adduct_Masses[Adduct][0], 25]) 

    if path == 3:
        Mass = tails[0][1]
        try: Mass += tails[1][1] - Masses['H2O']
        except: pass
        for tail in tails:
            name.append(tail[0])
        print(tails)


    new_Peaks = [] # Work around to remove duplicates from fragment list   
    for peak in Peaks:
        if f"{peak[0]:.6f}" not in (frag[0] for frag in new_Peaks): # Warning, do not look at this.
            new_Peaks.append([f"{peak[0]:.6f}", peak[1]])
    Peaks = new_Peaks

    new_list = new_peak_list(f"{species} {'_'.join(str(nam) for nam in name)} {Adduct}", Adduct_Masses[Adduct][1], Mass + Adduct_Masses[Adduct][0], Adduct, len(Peaks), Peaks)
    lipids[new_list["Name:"]] = new_list
    return

def calculate_acyl_tail(tail):

    tails = []
    for c_d in tail:
        if c_d[0] >= 2:
            Name = f"{c_d[0]}:{c_d[1]}"
            Mass = (31.98982926 + c_d[0]*14.01565007 - c_d[1]*2.01565007)
            tails.append([Name, Mass])
        else: pass
    return tails

def calculate_sphingo_tail(tail):

    tails = []

    c_d = tail[0]
    if c_d[0] >= 3:
        Name = f"S_{c_d[0]}:{c_d[1]}"
        Mass = (91.063329 + (c_d[0] - 3)*14.01565007 - c_d[1]*2.01565007)
        tails.append([Name, Mass])
    else: pass

    try: 
        c_d = tail[1]
        if c_d[0] >= 2:
            Name = f"{c_d[0]}:{c_d[1]}"
            Mass = (31.98982926 + c_d[0]*14.01565007 - c_d[1]*2.01565007)
            tails.append([Name, Mass])
        else: pass
    except: pass

    return tails

def new_spectrum_library():

    Cmin = 16   # Min and Max number of Carbons
    Cmax = 18 + 1
    Dmin = 0    # Min and Max number of Desaturated bonds
    Dmax = 2 + 1

    tails = []
    lipids = {}

    for C in range(Cmin, Cmax): # Limits all fatty acids within the ranges set
        for D in range(Dmin, Dmax):
            if D <= ((C-1)/2): # Limits the amount of desaturation a fatty acid can have
                if (FA_Options['Odd_Chains'] == False) and (C % 2 != 0):
                    continue
                else:
                    tails.append([C, D])
            else: continue

    for lipid in list(compress([lipid for lipid in Glycerolipids],[int(Glycerolipids[lipid]) for lipid in Glycerolipids])): # For MAG, DAG, TAG
        for Adduct in Adduct_Definitions[lipid]:
            for comb in cwr(tails, Chain_Number[lipid]):
                calculate_peaks(1, lipid, calculate_acyl_tail(comb), Adduct, lipids)

    for lipid in list(compress([lipid for lipid in Glycerophospholipids],[int(Glycerophospholipids[lipid]) for lipid in Glycerophospholipids])): # For Glycerophospholipids
        if Glycerophospholipids_Options['Lyso'] == True:
            for Adduct in Adduct_Definitions[lipid]:
                for comb in cwr(tails, 1):
                    calculate_peaks(2, lipid, calculate_acyl_tail(comb), Adduct, lipids)
        if Glycerophospholipids_Options['Diacyl'] == True:
            for Adduct in Adduct_Definitions[lipid]:
                for comb in cwr(tails, 2):
                    calculate_peaks(2, lipid, calculate_acyl_tail(comb), Adduct, lipids)

    for lipid in list(compress([lipid for lipid in Sphingolipids],[int(Sphingolipids[lipid]) for lipid in Sphingolipids])): # For Sphingolipids and N-acyl-Sphingolipids Currently Doesnt Work.
        for Adduct in Adduct_Definitions[lipid]:
            for comb in cwr(tails, Chain_Number[lipid]):
                calculate_peaks(3, lipid, calculate_sphingo_tail(comb), Adduct, lipids)
        
    return lipids

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

def buttonpush():

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

