Charge_Opts = {

    'POS':  True,  
    'NEG':  True}

Masses = {
    "[M-H]-":     
    [-1.007276, 'N', 1],

    "[M-2H]2-":   
    [-2.014552, 'N', 2],

    "[M-3H]3-":   
    [-3.021828, 'N', 3],

    "[M+Cl]-":    
    [34.969401, 'N', 1],

    "[M+H]+":     
    [1.007276,  'P', 1],

    "[M+H-H2O]+": 
    [-17.003289,'P', 1],

    "[M+Na]+":    
    [22.989221, 'P', 1],

    "[M+NH4]+":   
    [18.033826, 'P', 1],

    "H":           
    1.007276,

    "H2O":         
    18.010565}

Chain_parameters = [8, 20, 0, 4]

classes = { # Structural and spectra information defined here.
    
    'MAG':                                          # Class
     [True,                                        # [Generate/Don't generate,
     [1, 'Glycerol', None],                         # [# of tails, Body type ('Glycerol' for glycerolipids or 'SB' for sphingolipids), Headgroup mass (None if no headgroup)]]
     {"[M+H]+":     [['MA', 'MA_sub_H2O', 'FFAk'],  # { Adduct:   [[Types of fragments to generate],
                     [100, 50, 50]],                #              [Intensities of fragment]] }                        
      "[M+H-H2O]+": [['MA', 'FFAk'], 
                     [100, 50]],
      "[M+Na]+":    [['MA', 'MA_sub_H2O','FFAk', 'FFAkA'],
                     [100, 50, 50, 10]],
      "[M+NH4]+":   [['MA', 'MH', 'MH_sub_H2O', 'FFAk'],
                     [100, 100, 50, 50]]}], 
                                                    
    'DAG':                                          ##### Metabolites 2016, 6(3), 25; https://doi.org/10.3390/metabo6030025
     [True,                                        ##### J Am Soc Mass Spectrom. 2010 April ; 21(4): 657–669. https://doi.org/10.1016/j.jasms.2010.01.007
     [2, 'Glycerol', None],           
     {"[M+H]+":     [['MA', 'MH_sub_H2O', 'FFAk', 'MH_sub_FA'],
                     [25, 100, 1, 25]],
      "[M+H-H2O]+": [['MA', 'FFAk', 'MH_sub_FA'],
                     [25, 1, 25]],
      "[M+Na]+":    [['MA', 'MA_sub_H2O', 'FFAk', 'MH_sub_FA', 'FFAkA', 'MA_sub_FA'],
                     [25, 50, 1, 25, 5, 10]],
      "[M+NH4]+":   [['MA', 'MH', 'MH_sub_H2O', 'FFAk', 'MH_sub_FA'],
                     [25, 10, 100, 1, 25]]}],

    'TAG':                 
     [True, [3, 'Glycerol', None],           
     {"[M+Na]+":   [['MA', 'MH_sub_FA', 'MA_sub_FA', 'FFAk', 'FFAkA'], 
                    [25, 100, 10, 10, 5]],
     "[M+NH4]+":   [['MA', 'MH_sub_FA', 'FFAk'], 
                    [25, 100, 10]]}],

    'PA':                  
     [True, [2, 'Glycerol', 97.976895],           
     {"[M-H]-":    [['MA', 'gPA', 'FFA', 'MH_sub_FA', 'MH_sub_FAk'],
                    [25, 0, 100, 15, 5]],
      "[M+H]+":    [['MA', 'gPA', 'FFAk'],
                    [25, 0, 4]],
      "[M+Na]+":   [['MA', 'gPA', 'FFAkA'], 
                    [25, 0, 4]]}],

    'Lyso PA':      
     [True, [1, 'Glycerol', 97.976895],           
     {"[M-H]-":    [['MA', 'gPA', 'FFA', 'MH_sub_FA', 'MH_sub_FAk'],
                    [25, 0, 100, 15, 5]],
      "[M+H]+":    [['MA', 'gPA', 'FFAk'],
                    [25, 0, 4]],
      "[M+Na]+":   [['MA', 'gPA', 'FFAkA'],  # Should Lyso GPLs have [M+H-H2O]+ ?
                    [25, 0, 4]]}],

    'PC':                  
     [True, [2, 'Glycerol', 183.066044], # headgroup mass has -H to maintain neutral charge              
     {"[M+H]+":    [['MA', 'gPC', 'MH_sub_FA', 'MH_sub_FAk', 'FFAk'],
                    [25, 0, 4, 8, 4]],
      "[M+Na]+":   [['MA', 'gPC', 'MA_sub_FA', 'MA_sub_FAk', 'FFAkA'],
                    [25, 0, 4, 8, 4]]}],

    'Lyso PC':             
     [True, [1, 'Glycerol', 183.066044], # headgroup mass has -H to maintain neutral charge              
     {"[M+H]+":    [['MA', 'gPC', 'MH_sub_FA', 'MH_sub_FAk', 'FFAk'],
                    [25, 0, 4, 8, 4]],
      "[M+Na]+":   [['MA', 'gPC', 'MA_sub_FA', 'MA_sub_FAk', 'FFAkA'],
                    [25, 0, 4, 8, 4]]}],

    'PE':                  
     [True, [2, 'Glycerol', 141.019094],           
     {"[M-H]-":    [['MA', 'gPE', 'FFA', 'MH_sub_FA', 'MH_sub_FAk'],
                   [25, 0, 100, 15, 5]],
      "[M+H]+":    [['MA', 'gPE', 'FFAk'],
                   [25, 0, 4]],
      "[M+Na]+":   [['MA', 'gPE', 'FFAkA'],
                   [25, 0, 4]]}],

    'Lyso PE':             
     [True, [1, 'Glycerol', 141.019094],           
     {"[M-H]-":    [['MA', 'gPE', 'FFA', 'MH_sub_FA', 'MH_sub_FAk'],
                   [25, 0, 100, 15, 5]],
      "[M+H]+":    [['MA', 'gPE', 'FFAk'],
                   [25, 0, 4]],
      "[M+Na]+":   [['MA', 'gPE', 'FFAkA'],
                   [25, 0, 4]]}],

    'PG':                  
     [True, [2, 'Glycerol', 172.013674],           
     {"[M-H]-":   [['MA', 'gPG', 'FFA', 'MH_sub_FA', 'MH_sub_FAk'],
                   [25, 0, 100, 15, 5]],
      "[M+H]+":   [['MA', 'gPG', 'FFAk'],
                   [25, 0, 4]],
      "[M+Na]+":  [['MA', 'gPG', 'FFAkA'],
                   [25, 0, 4]]}],

    'Lyso PG':             
     [True, [1, 'Glycerol', 172.013674],           
     {"[M-H]-":   [['MA', 'gPG', 'FFA', 'MH_sub_FA', 'MH_sub_FAk'],
                   [25, 0, 100, 15, 5]],
      "[M+H]+":   [['MA', 'gPG', 'FFAk'],
                   [25, 0, 4]],
      "[M+Na]+":  [['MA', 'gPG', 'FFAkA'],
                   [25, 0, 4]]}],                                   # Perhaps I should add PGly ?
                                                    
    'PI':                  
     [False, [2, 'Glycerol', 260.029718],           
     {"[M-H]-":   [[],
                   []]}],

    'Lyso PI':             
     [False, [1, 'Glycerol', 260.029718],           
     {"[M-H]-":   [[],
                   []]}],

    'PIP':                 
     [False, [2, 'Glycerol', 339.996048],          
     {"[M-H]-":[],
      "[M-2H]2-":[]}],

    'PIP2':                
     [False, [2, 'Glycerol', 419.962378],
     {"[M-H]-":[],
      "[M-2H]2-":[],
       "[M-3H]3-":[]}],

    'PS':                  
     [False, [2, 'Glycerol', 185.008923],           
     {"[M-H]-":[],
      "[M+H]+":[],
       "[M+Na]+":[]}],

    'Lyso PS':             
     [False, [1, 'Glycerol', 185.008923],           
     {"[M-H]-":[],
      "[M+H]+":[],
       "[M+Na]+":[]}],

    'MGDG':                
     [False, [2, 'Glycerol', 180.063388], # headgroup mass is just of the galactose ring        
     {"[M+Na]+":[],
      "[M+NH4]+":[]}],

    'DGDG':                
     [False, [2, 'Glycerol', 342.116212], # headgroup mass is just of the galactose rings       
     {"[M+Na]+":[],
      "[M+NH4]+":[]}],

    'Sphingosine': # LipidBlast seemed to have only "[M+Na]+" spectra. Seen "[M+NH4]+" and "[M-H]-" in some articles.         
     [False, [0, 'SB', None],        
     {"[M-H]-":[],
      "[M+Cl]-":[],
       "[M+H]+":[],
        "[M+Na]+":[]}],
    
    'N-acylsphingosine': # Lipidblast also only has 1 fragment for these, and there aren't many characteristic studies.   
     [False, [1, 'SB', None],  
     {"[M-H]-":[],
      "[M+Cl]-":[],
       "[M+H]+":[],
        "[M+Na]+":[]}],             
    
    'Ceramide PA':         
     [False, [1, 'SB', 97.976895],                 
     {"[M-H]-":[],
      "[M+Cl]-":[],
       "[M+H]+":[],
        "[M+Na]+":[]}],
    
    'Ceramide PC':         
     [False, [1, 'SB', 183.066044], # headgroup mass has -H to maintain neutral charge                 
     {"[M-H]-":[],
      "[M+Cl]-":[],
       "[M+H]+":[],
        "[M+Na]+":[]}],
    
    'Ceramide PE':         
     [False, [1, 'SB', 141.019094],                 
     {"[M-H]-":[],
      "[M+Cl]-":[],
       "[M+H]+":[],
        "[M+Na]+":[]}],
    
    'Ceramide PG':         
     [False, [1, 'SB', 172.013674],                 
     {"[M-H]-":[],
      "[M+Cl]-":[],
       "[M+H]+":[],
        "[M+Na]+":[]}],
    
    'Ceramide PI':         
     [False, [1, 'SB', 260.029718],                 
     {"[M-H]-":[],
      "[M+Cl]-":[],
       "[M+H]+":[],
        "[M+Na]+":[]}], # Perhaps PIP and PIP2 should be included.
    
    'Ceramide PS':         
     [False, [1, 'SB', 185.008923],                
     {"[M-H]-":[],
      "[M+Cl]-":[],
       "[M+H]+":[],
        "[M+Na]+":[]}],
    
    'Ceramide MG':         
     [False, [1, 'SB', 180.063388], # headgroup mass is just of the galactose ring               
     {"[M-H]-":[],
      "[M+Cl]-":[],
       "[M+H]+":[],
        "[M+Na]+":[]}],
    
    'Ceramide DG':         
     [False, [1, 'SB', 342.116212], # headgroup mass is just of the galactose rings            
     {"[M-H]-":[],
      "[M+Cl]-":[],
       "[M+H]+":[],
        "[M+Na]+":[]}]}



### Below are the functions used to generate fragments for the above lipid species ###



def generate_peaks(key, mass, lipid, adduct): # Was meant to be simple...
    
    peaks = [] # Fragmentation spectra stored here as [[mz, intensity], [...], ... ] 
   
    def MA(intensity): # A == Adduct
        peaks.append([round((mass + Masses[adduct][0])
        /Masses[adduct][2], 6), intensity])

    def MA_sub_H2O(intensity):
        peaks.append([round((mass + Masses[adduct][0] - Masses['H2O'])
        /Masses[adduct][2], 6), intensity])

    def MA_sub_2H2O(intensity):
        peaks.append([round((mass + Masses[adduct][0] - 2*Masses['H2O'])
        /Masses[adduct][2], 6), intensity])

    def MA_sub_FA(intensity):
        peaks.extend([[round(mass + Masses[adduct][0] - tail[0], 6), intensity] 
         for tail in lipid])

    def MA_sub_FAk(intensity):
        peaks.extend([[round(mass + Masses[adduct][0] + Masses['H2O'] - tail[0], 6), intensity] 
         for tail in lipid])  


    def MH(intensity): # H == Hydrogen. Sometimes the adduct is lost during fragmentation, replaced w/ H+.
        if Masses[adduct][1] == 'P':
            peaks.append([round(mass + Masses['H'], 6), intensity])
        if Masses[adduct][1] == 'N': 
            peaks.append([round(mass - Masses['H'], 6), intensity])  
    
    def MH_sub_H2O(intensity):
        if Masses[adduct][1] == 'P':
            peaks.append([round(mass + Masses['H'] - Masses['H2O'], 6), intensity])
        if Masses[adduct][1] == 'N':
            peaks.append([round(mass - Masses['H'] - Masses['H2O'], 6), intensity])
    
    def MH_sub_2H2O(intensity):
        if Masses[adduct][1] == 'P':
            peaks.append([round(mass + Masses['H'] - 2*Masses['H2O'], 6), intensity])
        if Masses[adduct][1] == 'N':
            peaks.append([round(mass - Masses['H'] - 2*Masses['H2O'], 6), intensity])
    
    def MH_sub_FA(intensity):
        if Masses[adduct][1] == 'P':
            peaks.extend([[round(mass + Masses['H'] - tail[0], 6), intensity] 
             for tail in lipid])
        if Masses[adduct][1] == 'N':
            peaks.extend([[round(mass - Masses['H'] - tail[0], 6), intensity] 
             for tail in lipid])
    
    def MH_sub_FAk(intensity):
        if Masses[adduct][1] == 'P':
            peaks.extend([[round(mass + Masses['H'] + Masses['H2O'] - tail[0], 6), intensity] 
             for tail in lipid])
        if Masses[adduct][1] == 'N':
            peaks.extend([[round(mass - Masses['H'] + Masses['H2O'] - tail[0], 6), intensity] 
             for tail in lipid])


    
    def FA_sub_H(intensity): 
        peaks.extend([[round(tail[0] - Masses['H'], 6), intensity]
         for tail in lipid if 'E' not in tail[1]])
    
    def FAk_H(intensity):
        if Masses[adduct][1] == 'P':
            peaks.extend([[round(tail[0] - Masses['H2O'] + Masses['H'], 6), intensity]
             for tail in lipid if 'E' not in tail[1]])
        if Masses[adduct][1] == 'N':
            peaks.extend([[round(tail[0] - Masses['H2O'] - Masses['H'], 6), intensity]
             for tail in lipid if 'E' not in tail[1]])
    
    def FAk_A(intensity):
        peaks.extend([[round(tail[0] - Masses['H2O'] + Masses[adduct][0], 6), intensity]
         for tail in lipid if 'E' not in tail[1]])
    
    def HG_FA_NL(HG): # Fragments common to some GPLs. Loss of headgroup and tail.
        if Masses[adduct][1] == 'P':
            peaks.extend([[round(mass - HG - tail[0] + Masses['H2O'] + Masses[adduct][0], 6), 8]
            for tail in lipid])
        if Masses[adduct][1] == 'N':
            peaks.extend([[round(mass - HG - tail[0] - Masses['H'], 6), 15] # Not sure if these should be + Masses[adduct][0] rather than - Masses['H']. If the prior, could be done for both polarities?
            for tail in lipid])
            peaks.extend([[round(mass - HG - tail[0] + Masses['H2O'] - Masses['H'], 6), 5]
            for tail in lipid])


    def PA_gpl_Character(Place_Holder):
        HG_FA_NL(97.976895)
        if Masses[adduct][1] == 'P':
            peaks.extend([[round(mass - 97.976895 + Masses['H'], 6), 100]])
        if Masses[adduct][1] == 'N':
            peaks.extend([[78.959053,   5],
                          [96.969618,   5],
                          [152.995833, 10]])

    def PC_gpl_Character(Place_Holder):
        #HG_FA_NL(183.066044) # ?
        if Masses[adduct][1] == 'P':
            peaks.extend([[round(183.066044 + Masses['H'], 6), 100]])
            if adduct == "[M+Na]+":
                peaks.extend([[round(mass - 59.073499  + Masses[adduct][0], 6), 25], #(Loss of N(CH3)3) Only happens for Alkaline salts?
                              [round(mass - 183.066044 + Masses[adduct][0], 6),  4],
                              [round(mass - 183.066044 + Masses['H'],       6),  4]]) # Should this be for all positive adducts?
        if Masses[adduct][1] == 'N':
            peaks.extend([[78.959053,  5],
                         [96.969618,   5],
                         [152.995833, 10]])

    def PE_gpl_Character(Place_Holder):
        HG_FA_NL(141.019094)
        if Masses[adduct][1] == 'P':
            peaks.extend([[round(mass - 141.019094 + Masses['H'], 6), 100]])
        if Masses[adduct][1] == 'N':
            peaks.extend([[78.959053,   5],
                          [96.969618,   5],
                          [152.995833, 10],
                          [140.011817,  3],
                          [196.038032,  5]])

    def PG_gpl_Character(Place_Holder):
        
        if Masses[adduct][1] == 'P':
            HG_FA_NL(172.013674) # Neutral loss is headgroup mass
            peaks.extend([[round(mass - 172.013674 + Masses['H'], 6), 100]])
        if Masses[adduct][1] == 'N':
            HG_FA_NL(74.036779) # Neutral loss is different to headgroup mass
            peaks.extend([[78.959053,   5],
                          [96.969618,   5],
                          [152.995833, 10],
                          [209.022048,  5],
                          [227.032612,  5],
                          [245.043177,  5]])

    def PI_gpl_Character(Place_Holder):
        if Masses[adduct][1] == 'P':
            pass
        if Masses[adduct][1] == 'N':
            peaks.extend([[78.959053,   5],
                          [96.969618,   5],
                          [152.995833, 10],
                          [223.001312,  5],
                          [241.011877, 25],
                          [259.022442,  5],
                          [297.038092,  5],
                          [315.048656,  5]])

    def PIP_gpl_Character(Place_Holder):
        if Masses[adduct][1] == 'P':
            pass
        if Masses[adduct][1] == 'N':
            peaks.extend([[78.959053,   5],
                          [96.969618,   5],
                          [152.995833, 10],
                          [223.001312,  5],
                          [241.011877, 25],
                          [259.022442,  5],
                          [297.038092,  5],
                          [315.048656,  5],
                          [302.957642,  5],
                          [320.978207,  5],
                          [395.014986,  5],
                          [159.985465,  5],
                          [168.9907475, 5],
                          [mass - 80.973606, 5],
                          [mass - 80.973606 - Masses['H2O'], 5],
                          [mass - Masses['H'] - Masses['H2O'], 5]])
    
    def PIP2_gpl_Character(Place_Holder):
        if Masses[adduct][1] == 'P': pass
        if Masses[adduct][1] == 'N':
            peaks.extend([[78.959053,   5],
                          [96.969618,   5],
                          [152.995833, 10],
                          [223.001312,  5],
                          [241.011877, 25],
                          [259.022442,  5],
                          [297.038092,  5],
                          [315.048656,  5],
                          [mass - 80.973606, 5],
                          [mass - 80.973606 - Masses['H2O'], 5],
                          [mass - Masses['H'] - Masses['H2O'], 5]])
    
    def PS_gpl_Character(Place_Holder):
        if Masses[adduct][1] == 'P':
            peaks.extend([[round(mass - 185.008923 + Masses['H'], 6), 100]])
        if Masses[adduct][1] == 'N':
            peaks.extend([[78.959053,   5],
                          [96.969618,   5],
                          [152.995833, 10]])
    
    def MG_gpl_Character():
        if Masses[adduct][1] == 'P':
            peaks.extend([[round(mass - 180.063388 + Masses['H2O'] + Masses['H']), 50]])
        if Masses[adduct][1] == 'N':
            pass
    
    def DG_gpl_Character():
        if Masses[adduct][1] == 'P':
            peaks.extend([[round(mass - 342.116212 + Masses['H2O'] + Masses['H']), 50]])
        if Masses[adduct][1] == 'N':
            pass

    f = {
        'MA':           MA,
        'MA_sub_H2O':   MA_sub_H2O,
        'MA_sub_2H2O':  MA_sub_2H2O,
        'MA_sub_FA':    MA_sub_FA,
        'MA_sub_FAk':   MA_sub_FAk,

        'MH':           MH,
        'MH_sub_H2O':   MH_sub_H2O,
        'MH_sub_2H2O':  MH_sub_2H2O,
        'MH_sub_FA':    MH_sub_FA,
        'MH_sub_FAk':   MH_sub_FAk,
        
        'FFA':          FA_sub_H,
        'FFAk':         FAk_H,
        'FFAkA':        FAk_A,

        'gPA':           PA_gpl_Character,
        'gPC':           PC_gpl_Character,
        'gPE':           PE_gpl_Character,
        'gPG':           PG_gpl_Character,
        'gPI':           PI_gpl_Character,
        'gPIP':          PIP_gpl_Character,
        'gPIP2':         PIP2_gpl_Character,
        'gPS':           PS_gpl_Character,
        'gMG':           MG_gpl_Character,
        'gDG':           DG_gpl_Character
    }

    spectra = classes[key][2][adduct]
    for func, intens in zip(spectra[0], spectra[1]):
        f[func](intens)

    keep = [] # Remove duplicates
    seen = set([])
    for peak in peaks:
        if peak[0] not in seen:
            keep.append(peak)
            seen.add(peak[0])
    keep.sort() # Sort fragments based on mass
    peaks = keep

    return peaks