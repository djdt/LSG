import Definitions as D
import tkinter as tk
from itertools import combinations_with_replacement as cwr
from itertools import chain
from itertools import product
from itertools import compress
import math
import time

Charge_Opts = D.Charge_Opts
Masses = D.Masses
classes = D.classes

def generate_lipids(n): # requires limits for saturation (d0, d1) and limits for carbon length (c0, c1), stored in n

    def Sp():
        a = list(compress([classes[lipid][1][1] for lipid in classes],
         [int(classes[lipid][0]) for lipid in classes]))
        if 'SB' in a: return True
        else: return False

    def generate_tails():
        Tail_Type =   {'Ac':   True,
                       'Pl':   False,
                       'Sp':   Sp()}
        Tail_Opts =   {'Odd':  True}
        loc = {'Ac':0, 'Pl':0, 'Sp':1}
        tails = [[],[]]
        types = list(compress([t for t in Tail_Type], # Compresses keys for types based on true/false into a list
         [int(Tail_Type[t]) for t in Tail_Type]))

        def mass_name(c, d):
            definitions = {'Ac':[31.98982926, 0, 0,  ''], # Acyl chains
                           'Pl':[15.994915, 0, 1,  'E'], # Plasmenyl chains
                           'Sp':[91.063329, 3, 0, 'Sp_']} # Sphingoid bases
            i = definitions[type]
            return [i[0] + (c-i[1])*14.01565007 - (d-i[2])*2.01565007, f"{i[3]}{c}:{d}"]
        for type in types: 
            tails[loc[type]].extend([mass_name(c, d) # Separates glycero and sphingo lipids
                for c in range(n[0], n[1] + 1) 
                    if (Tail_Opts['Odd'] == False and (c % 2 == 0))
                    or (Tail_Opts['Odd'] == True)
                        for d in range(n[2], n[3] + 1) if d <= ((c-1)/2)]) # Limits desaturation based on C number
        return tails
    
    def package_lipid(key):
        body = classes[key][1][1]
        N_Tails = classes[key][1][0]
        if body == 'Glycerol': lst = [list(comb) 
         for comb in cwr(tails[0], N_Tails)]
        elif body == 'SB': lst = [list(filter(None, comb))
         for comb in product(tails[1], [tail * N_Tails
          for tail in tails[0]])] # Hacky, only works for 0 and 1.
        return lst

    tails = generate_tails() # Tails created here
    species = {key:package_lipid(key) # Lipids packaged here
     for key in compress([lipid for lipid in classes], [int(classes[lipid][0]) 
      for lipid in classes])}

    return species

def new_library(lipids):

    def new_spectra(key, lipid, adducts):

        def new_list(): # .mgf format
            lst = {
            "Name:":f"{key} {'_'.join(str(nam) for nam in name)} {adduct}",
            "Spectrum_type:":"MS2",
            "Ion_mode:":Masses[adduct][1],
            "PrecursorMZ:":f"{mass + Masses[adduct][0]:.6f}",
            "Precursor_type:":adduct,
            "Num Peaks:":len(peaks),
            "peaks":peaks} 
            return lst       

        Mass = { # If statements results into tuples, whereas dictionary gives floats.
         'Glycerol': 92.047344 + sum((tail[0] - Masses['H2O']) 
           for tail in lipid), # Sums tails for glycerolipids
         'SB': Masses['H2O'] + sum(tail[0] - Masses['H2O'] 
           for tail in lipid)} # Sums tails for sphingolipids
        mass = Mass[classes[key][1][1]] # classes[key][1][1] is backbone
        if classes[key][1][2]: mass += (classes[key][1][2] - Masses['H2O']) # Adds headgroup mass

        name = [tail[1] for tail in lipid] # Put together the name from tail composition and class
        if len(name) < 2: name.append('0:0')  # Adds 0:0 to the name for Lyso GPLs / Sphingosines
        if key in ('MAG', 'DAG'): name.append('0:0')  # Adds another 0:0 to the name for MAGs and DAGs to complement TAGs
        
        lst = []
        for adduct in adducts:
            if (D.Masses[adduct][1] == 'P' and D.Charge_Opts['POS']) or (D.Masses[adduct][1] == 'N' and D.Charge_Opts['NEG']):
                peaks = D.generate_peaks(key, mass, lipid, adduct)
                lst.append(new_list())
            else: pass
        return lst

    library = {} # This is where the fragments will be stored
    for key in list(lipids):
        library[key] = []
        for lipid in lipids[key]:
            library[key].append(new_spectra(key, lipid, classes[key][2]))

    return(library)

def buttonpush():

    def write_to_file(spectrum_library):

        if D.Charge_Opts['POS'] == True: Poslibrary = open("Positive_Library.msp", "a")
        if D.Charge_Opts['NEG'] == True: Neglibrary = open("Negative_Library.msp", "a")

        count = 0

        for classes in spectrum_library: # For each lipid
            for lipid in spectrum_library[classes]:
                for adduct in lipid:

                    if adduct['Ion_mode:'] == 'P' and D.Charge_Opts['POS'] == True:
                        for line in adduct: # Write the lines
                            if line == "peaks": Poslibrary.writelines([f"{peak[0]} {peak[1]}\n" for peak in adduct[line]]) # if spectra, whole thing can be written at once
                            else: Poslibrary.write(f"{line} {adduct[line]}\n") # Information is written one line at a time to look for where the spectra is  
                        Poslibrary.write("\n") # After it finishes a spectrum, new line ready for the next
                        count += 1

                    if adduct['Ion_mode:'] == 'N' and D.Charge_Opts['NEG'] == True:
                        for line in adduct: # Write the lines
                            if line == "peaks": Neglibrary.writelines([f"{peak[0]} {peak[1]}\n" for peak in adduct[line]]) # if spectra, whole thing can be written at once
                            else: Neglibrary.write(f"{line} {adduct[line]}\n") # Information is written one line at a time to look for where the spectra is      
                        Neglibrary.write("\n") # After it finishes a spectrum, new line ready for the next
                        count += 1
        
        try: Poslibrary.close()
        except: pass
        try: Neglibrary.close()
        except: pass

        return count

    t0 = time.time()

    print('Generating...')
    lib = new_library(generate_lipids(D.Chain_parameters))
    print("Writing...")
    count = write_to_file(lib)

    t1 = time.time()
    
    print(f"Produced {count} spectra in {t1 - t0:.4f} seconds")   


root = tk.Tk()

button = tk.Button(root, text = "Push me", command = buttonpush)
button.pack()
tk.mainloop()