#!/usr/bin/env python
# coding: utf-8

"""

Autor: Jose Cleydson F Silva
Data: 05/12/2023
Versão: 1.1 delta

"""

from Bio.PDB import *
import nglview as nv


parser = PDBParser()
structure = parser.get_structure("NLR", "6j5t.pdb")
atoms = structure.get_atoms()
residues = structure.get_residues()
chains = structure.get_chains()

hydrogen_B = 0
lines_total = 0
lines_total2 = 0
refined = 0
HC_total2 = 0.0
VDW_total = 0.0
repulsive = 0
london = 0.0
RT = 0.0
surface_tension = 0.0
hydrophobicity = 0.0


# In[33]:


def split(s, c, v):
    i = 0
    j = s.find(c)
    while j != -1:
        a = s[i:j]
        if a != "":
            v.append(a)
        i = j + 1
        j = s.find(c, i)
        if j == -1:
            a = s[i:]
            if a != "" and a != "\n":
                v.append(a)


# In[34]:


def create_map_tables(param_file_name):
    param_map = {}
    with open(param_file_name) as param_file:
        for param_line in param_file:
            param_v = []
            split(param_line, ':', param_v)
            param_map[param_v[0]] = float(param_v[1])
            # print(param_v[0], ": ", param_v[1])
    return param_map


# In[35]:


def surface_tension_hydrophobicity_calculator(hydro_map, tension_map, infile_name):
    with open(infile_name) as infile:
        aa_numbers = []
        counter = 0
        for p_dist_line in infile:
            v = []
            split(p_dist_line, ' ', v)
            exists = v[4] in aa_numbers
            aa_numbers.append(v[4])
            counter += 1
            if not exists:
                global hydrophobicity
                hydrophobicity += hydro_map[v[2]]
                global surface_tension
                surface_tension += tension_map[v[2]]


# In[36]:


def fncomp(lhs, rhs):
    return lhs < rhs

class classcomp:
    def __call__(self, lhs, rhs):
        return lhs < rhs


# In[37]:


def replace_all(str, from_str, to_str):
    while True:
        start_pos = str.find(from_str)
        if start_pos == -1:
            break
        end_pos = start_pos + len(from_str)
        str = str[:start_pos] + to_str + str[end_pos:]
    return str


# In[39]:


protein_name = ""
ligand_name = ""
ligand_pdb = ""
ligand_type = ""

def summation(a, b):
    return sum(range(int(a), int(b) + 1))


# In[40]:


def open_ligand():
    achou_mol2 = 0
    try:
        inlig = open(ligand_name, "r")
    except FileNotFoundError:
        print("Error opening input file...\n\nPlease, be sure that the *.PDB file you want to open is placed in the same\ndirectory of this program, that you wrote the complete and correct name and\n.pdb extention, and try again.\n")
        exit(1)
    pos = 0
    outlig = open("ligand.tmp", "w")
    line = ''
    while True:
        line = inlig.readline()
        if not line:
            break
        if line == "@<TRIPOS>ATOM":
            pos = 0
            line = ' ' + inlig.readline()
            while line[0] == ' ':
                line = inlig.readline()
                if line[0] != '@':
                    # ... get 52 characters (example)
                    # 1 N          15.0880   10.7980   23.5470 N.4       1 MOL         0.2328
                    outlig.write(line[:52] + "\n")
            achou_mol2 = 1

    outlig.write("\n") #remove this once I fix the bug that makes it necessary to have an extra new line when reading the ligand input.
    if achou_mol2 == 0:
        print(f"Input error: {ligand_name} doesn't look to be a *.mol2 file. Please check input files and order of input data in comand line.")
        exit(1)
    inlig.close()
    outlig.close()


# In[43]:


def distancia():
    
    import os

    pdistmp = open("P_DIST.tmp", "w")
    ldistmp = open("L_DIST.tmp", "w")
    ligtmp = open("ligand.tmp", "r")

    #-----------------------CALCULOS----------------------------
    n = 0 #conta posicao do caractere
    i = 0

    #-----------------------------------------DADOS DO LIGANTE----------------------
    atomo = ""
    nome_atomo = ""
    x = ""
    y = ""
    z = ""
    enter = ""
    num_atomo = 0
    num_atomo_copia = 0
    num_x = 0.0
    num_y = 0.0
    num_z = 0.0

    #distancia permitida
    permitido = 0.0
    #RADII 1 E RADII 2
    r1 = 0.0
    r2 = 0.0

    #------------------------------------------DADOS DA PROTEINA----------------------------------
    p_atomo = ""
    p_nome_atomo = ""
    p_aminoacido = ""
    p_num_aminoacido = ""
    p_x = ""
    p_y = ""
    p_z = ""
    chain = ""
    num_p_atomo = 0
    num_p_atomo_copia = 0
    num_px = 0.0
    num_py = 0.0
    num_pz = 0.0
    counter = 0
    protmp = open("protein.tmp", "r") #vai abrir o arquivo temporario da proteina para poder encontrar a raiz do doador ou aceitador e doador ou aceitador
    RT_found = open("RT_found.tmp", "w")
    HC_lig_contact = open("HC_lig_contact.tmp", "w")
    HC_lig_contact_lim = open("HC_lig_contact_lim.tmp", "w")

    #------------------------------------------------------------------------------------------
    line = "*"
    hibr = ""
    v = []
    found_hibr = 0
    found_p = 0
    
    while not ligtmp.eof():
        found_hibr = 0
        line = ligtmp.readline()
        v = []
        line = line.replace("-", " -")
        split_line = line.split(' ')
        if line.startswith('') or line.startswith('\n') or line.startswith('\0'):
            break
        
        p_atomo = v[0];
        num_p_atomo = atoi(v[0].c_str()); # Preciso ver essa funcao
        p_nome_atomo = v[1];
        p_aminoacido = v[2];
        chain = v[3];
        p_num_aminoacido = v[4];
        num_px = atof(v[5].c_str());
        p_x = v[5];
        (/cout, <<, num_px, <<, endl;)
        num_py = atof(v[6].c_str());
        p_y = v[6];
        (/cout, <<, num_py, <<, endl;)
        num_pz = atof(v[7].c_str());
        p_z = v[7];
        
        DX = (num_x - num_px);
        DY = (num_y - num_py);
        DZ = (num_z - num_pz);
        EDX = (DX*DX);
        EDY = (DY*DY);
        EDZ = (DZ*DZ);
        SD = (EDX+EDY+EDZ);
        DIS = sqrt(SD);

        if hibr.startswith("N.3"):
            r2 = 1.87
            found_hibr = 1
        elif hibr.startswith("N.1") or hibr.startswith("N.2") or hibr.startswith("N.ar") or hibr.startswith("N.pl3"):
            r2 = 1.86
            found_hibr = 1
        elif hibr.startswith("N.am"):
            r2 = 1.83
            found_hibr = 1
        elif hibr.startswith("C.3"):
            r2 = 1.94
            found_hibr = 1
        elif hibr.startswith("C.1") or hibr.startswith("C.2"):
            r2 = 1.90
            found_hibr = 1
        elif hibr.startswith("C.ar"):
            r2 = 1.85
            found_hibr = 1
        elif hibr.startswith("O.3"):
            r2 = 1.74
            found_hibr = 1
        elif hibr.startswith("O.2") or hibr.startswith("O.co2"):
            r2 = 1.66
            found_hibr = 1
        elif hibr.startswith("F") or hibr.startswith("O.w"):
            r2 = 1.77
            found_hibr = 1
        elif hibr.startswith("S"):
            r2 = 2.01
            found_hibr = 1
        elif hibr.startswith("S.3"):
            r2 = 2.09
            found_hibr = 1
        elif hibr.startswith("Cl"):
            r2 = 2.00
            found_hibr = 1
        elif hibr.startswith("Br"):
            r2 = 2.22
            found_hibr = 1
        elif hibr.startswith("I"):
            r2 = 2.42
            found_hibr = 1
        elif hibr.startswith("P"):
            r2 = 2.03
            found_hibr = 1
        DIS2 = 0
        linep = ""

        #--------LIGANTE-------------RADII-Raio-de-VDW-PERMITIDOS
        #cout << "*" << atomo << "*" << espaco_1 << "*" << nome_atomo << "*" << espaco_2 << "*" << espaco_3 << "*" << x << "*" << y << "*" << z << "*" << espaco_5 << "*" << hibr << "*" << enter;
        #cout << endl;
        #cout << "number: " << num_x << "name: " << v[2] << endl;
        DIS2 = 0
        linep = ""
        while protmp:
            found_p = 0
            linep = protmp.readline()
            #print(linep)
            v = []
            line = line.replace("-", " -")
            v = linep.split()
            if v[0] == "" or v[0] == "\0" or v[0] == "\n":
                break
            #print("*" + p_atomo + "*" + p_nome_atomo + "*" + p_aminoacido + "*" + p_num_aminoacido + "*" + px + "*" + py + "*" + pz + "*" + p_enter)
            if linep[:1] == "" or linep[:1] == "\n" or linep[:1] == "\0":
                break
            #print("")
            #print("convertendo caracteres da proteina\n")
            p_atomo = v[0]
            num_p_atomo = int(v[0])
            p_nome_atomo = v[1]
            p_aminoacido = v[2]
            chain = v[3]
            p_num_aminoacido = v[4]
            num_px = float(v[5])
            p_x = v[5]
            #print(num_px)
            num_py = float(v[6])
            p_y = v[6]
            #print(num_py)
            num_pz = float(v[7])
            p_z = v[7]
            #double DDA = sqrt(((num_x-num_px)*(num_x-num_px))+((num_y-num_py)*(num_y-num_py))+((num_z-num_pz)+(num_z-num_pz))); 
            #distancia entre o atomo do doador e aceitador  
            #angulos
            #cos-1=[(AB)2 + (AC)2 - (BC)2 / 2(AB)*(AC)]
            DX = (num_x - num_px)
            DY = (num_y - num_py)
            DZ = (num_z - num_pz)
            EDX = (DX*DX)
            EDY = (DY*DY)
            EDZ = (DZ*DZ)
            SD = (EDX+EDY+EDZ)
            DIS = sqrt(SD)
            #print(DIS)
            DIS2 = DIS
            #print("")
            #CONTINUACAO ACHA DISTANCIA PERMITIDA
            #print("\n*" + p_nome_atomo + "*\n")
            #print("\n*" + nome_atomo + "*\n")
            if p_nome_atomo[:1] == "N":
                r1 = 1.65; found_p = 1
                #print(" found : " + p_nome_atomo + "\n")
            if p_nome_atomo[:1] == "O":
                r1 = 1.40; found_p = 1
            if p_nome_atomo[:2] == "NZ":
                r1 = 1.50; found_p = 1
            if p_nome_atomo[:1] == "S":
                
            if p_nome_atomo.startswith('C'):
                r1 = 1.76
                found_p = 2

            if p_nome_atomo.startswith(('CA', 'CB', 'CD', 'CE', 'CG', 'CH', 'CZ')):
                r1 = 1.87
                found_p = 2

            '''
            1.65    N, NE, NH1, NH2, ND2, NE2, ND1
            1.40    O, OD1, OD2, OE1, OE2, OH, 
            1.85    SG
            1.50    NZ
            '''

            permitido = r1 + r2
            permitido1 = permitido - 0.7
            extended = 0.7

            # RT
            if DIS2 <= permitido + extended and not nome_atomo.startswith('H') and found_hibr == 1 and found_p == 1:
                for i in range(6 - len(atomo)):
                    RT_found.write(" ")
                RT_found.write(atomo + '\n')

            # VDW
            # [JC] van der Waals interactions (VDWs)
            if not hibr.startswith('H') and not atomo.startswith('H') and found_hibr == 1 and found_p != 0:
                # function for VDW
                VDW_radii = r1 + r2
                var_0 = VDW_radii / DIS2
                exp1 = 8.0
                exp2 = 4.0
                var_1 = var_0 ** exp1
                var_2 = var_0 ** exp2
                VDW_int = var_1 - (2 * var_2)

                if VDW_int >= 100:
                    VDW_int = 100
                if VDW_int <= 100:
                    VDW_total += VDW_int
                VDW_int = 0

            # HC
            # [JC] Hydrophobic contacts (HCs)
            HC = 0
            HC2 = 0
            HC_exp = 0

            if p_nome_atomo.startswith('C') and nome_atomo.startswith('C') and found_hibr == 1 and found_p != 0:
                HC_VDW = permitido
                HC_allowed_1 = HC_VDW + 0.5
                HC_allowed_2 = HC_VDW + 2.0

            # ---------------------HC2--------------------------------
            if DIS2 <= HC_allowed_1:
                HC2 = 1

            if DIS2 > HC_allowed_1 and DIS2 <= HC_allowed_2:
                HC2 = (1 / 1.5) * ((HC_VDW + 2.0) ** 2 - DIS2 ** 2)

            if DIS2 > HC_allowed_2:
                HC2 = 0

            HC_total2 += HC2
            if DIS2 <= permitido and nome_atomo[:1] != "H" and found_hibr == 1:
                repulsive += 1
            if london + 1 / DIS2 != float("inf"):
                london += 1 / DIS2

            # ---------------------HC2--------------------------------
            # [JC] Hydrophobic contacts (HCs) 2
            if DIS2 >= permitido1 and DIS2 <= permitido and found_hibr == 1 and found_p == 1:
                if (p_nome_atomo[:1] == "N" or p_nome_atomo[:1] == "O" or p_nome_atomo[:1] == "S") and (nome_atomo[:1] == "N" or nome_atomo[:1] == "O" or nome_atomo[:1] == "S"):
                    n = 0
                    if len(atomo) < 7:
                        for n in range(7 - len(atomo)):
                            ldistmp.write(" ")
                    ldistmp.write(atomo + "  " + nome_atomo)
                    if len(nome_atomo) < 7:
                        for n in range(7 - len(nome_atomo)):
                            ldistmp.write(" ")
                    for n in range(10 - len(x)):
                        ldistmp.write(" ")
                    ldistmp.write(x)
                    for n in range(10 - len(y)):
                        ldistmp.write(" ")
                    ldistmp.write(y)
                    for n in range(10 - len(z)):
                        ldistmp.write(" ")
                    ldistmp.write(z + " " + hibr)
                    if len(hibr) < 5:
                        for n in range(5 - len(hibr)):
                            ldistmp.write(" ")
                    ldistmp.write("\n")
                    n = 0
                    if len(p_atomo) < 6:
                        for n in range(6 - len(p_atomo)):
                            pdistmp.write(" ")
                    pdistmp.write(p_atomo + " " + p_nome_atomo)
                    if len(p_nome_atomo) < 4:
                        for n in range(4 - len(p_nome_atomo)):
                            pdistmp.write(" ")
                    if len(p_aminoacido) < 4:
                        for n in range(4 - len(p_aminoacido)):
                            pdistmp.write(" ")
                    pdistmp.write(p_aminoacido)
                    if len(chain) < 4:
                        for n in range(4 - len(chain)):
                            pdistmp.write(" ")
                    pdistmp.write(chain)
                    if len(p_num_aminoacido) < 4:
                        for n in range(4 - len(p_num_aminoacido)):
                            pdistmp.write(" ")
                    if len(p_num_aminoacido) >= 4 and len(p_num_aminoacido) < 9:
                        for n in range(9 - len(p_num_aminoacido)):
                            pdistmp.write(" ")
                    pdistmp.write(p_num_aminoacido)
                    if len(p_x) < 9:
                        for n in range(9 - len(p_x)):
                            pdistmp.write(" ")
                    pdistmp.write(p_x)
                    if len(p_y) <= 8:
                        for n in range(8 - len(p_y)):
                            pdistmp.write(" ")
                    pdistmp.write(p_y)
                    if len(p_z) <= 8:
                        for n in range(8 - len(p_z)):
                            pdistmp.write(" ")
                    pdistmp.write(p_z + " \n")
                    n = 0

        protmp.clear();              # forget we hit the end of file
        protmp.seekg(0, ios::beg);   # move to the start of the file

	protmp.close(); 
	pdistmp.close();
	ldistmp.close();  
	RT_found.close();
	HC_lig_contact.close();
	HC_lig_contact_lim.close();


# In[ ]:


def raiz_proteina():
    root = open("P_ROOT_0.tmp", "w")
    AD = open("P_DIST.tmp", "r")
    p_atomo = ""
    p_nome_atomo = ""
    p_aminoacido = ""
    p_num_aminoacido = ""
    px = ""
    py = ""
    pz = ""
    p_enter = ""
    num_p_num_aminoacido = 0
    num_pp_num_aminoacido = 0
    num_p_atomo = 0
    num_p_atomo_copia = 0
    num_px = 0.0
    num_py = 0.0
    num_pz = 0.0
    n = 0
    pp_atomo = ""
    pp_nome_atomo = ""
    pp_aminoacido = ""
    pp_num_aminoacido = ""
    ppx = ""
    ppy = ""
    ppz = ""
    pp_enter = ""
    num_pp_atomo = 0
    num_pp_atomo_copia = 0
    num_ppx = 0.0
    num_ppy = 0.0
    num_ppz = 0.0
    # tentar abrir arquivo da proteina
    while not AD.eof():
        pp_atomo = AD.read(8)
        pp_nome_atomo = AD.read(5)
        pp_aminoacido = AD.read(4)
        pp_num_aminoacido = AD.read(8)
        ppx = AD.read(9)
        ppy = AD.read(9)
        ppz = AD.read(9)
        pp_enter = AD.read(3)
        #print("1", pp_atomo, "2", pp_nome_atomo, "3", pp_aminoacido, "4", pp_num_aminoacido, "5", ppx, "6", ppy, "7", ppz, "8", pp_enter)
        #print(pp_nome_atomo[0], pp_nome_atomo[1], pp_nome_atomo[2])
        AD_root = open("protein.tmp")
        aa = 0
        
        while True:
            aa += 1
            p_atomo = AD_root.read(7)
            p_nome_atomo = AD_root.read(4)
            p_aminoacido = AD_root.read(3)
            p_num_aminoacido = AD_root.read(7)
            px = AD_root.read(8)
            py = AD_root.read(8)
            pz = AD_root.read(8)
            p_enter = AD_root.read(2)
            # print("1", p_atomo, "2", p_nome_atomo, "3", p_aminoacido, "4", p_num_aminoacido, "5", px, "6", py, "7", pz, "8", p_enter)
            # -------------------------------N-PROLINA
            # [JC]
            
            if (pp_nome_atomo[0]=='N' && pp_aminoacido[0]=='P' && pp_aminoacido[1]=='R' && pp_aminoacido[2]=='O'):
                if aa!=1:
                    num_p_num_aminoacido = int(p_num_aminoacido)
                    num_pp_num_aminoacido = int(pp_num_aminoacido)
                    if p_nome_atomo[0]=='C' and p_nome_atomo[1]==' ' and num_p_num_aminoacido==(num_pp_num_aminoacido -1 ):
                        # print("R", p_atomo, "*", p_nome_atomo, "*", p_aminoacido, "*", p_num_aminoacido,  "*", px,  "*", py, "*", pz, "*", p_enter) 
                        # cria arquivo temp com as raizes
                        root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)
                    if p_nome_atomo[0]=='C' and p_nome_atomo[1]=='A' and p_nome_atomo[3]==' ' and pp_aminoacido[0]==p_aminoacido[0] and pp_aminoacido[1]==p_aminoacido[1] and pp_aminoacido[2]==p_aminoacido[2] and pp_num_aminoacido[0]==p_num_aminoacido[0] and pp_num_aminoacido[1]==p_num_aminoacido[1] and pp_num_aminoacido[2]==p_num_aminoacido[2] and pp_num_aminoacido[3]==p_num_aminoacido[3] and pp_num_aminoacido[4]==p_num_aminoacido[4] and pp_num_aminoacido[5]==p_num_aminoacido[5] and pp_num_aminoacido[6]==p_num_aminoacido[6] and pp_num_aminoacido[7]==p_num_aminoacido[7]:
                        # print("R", p_atomo, "*", p_nome_atomo, "*", p_aminoacido, "*", p_num_aminoacido,  "*", px,  "*", py, "*", pz, "*", p_enter)   
                        root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)
                    if p_nome_atomo[0] == 'C' and p_nome_atomo[1] == 'D' and p_nome_atomo[3] == ' ' and pp_aminoacido[0] == p_aminoacido[0] and pp_aminoacido[1] == p_aminoacido[1] and pp_aminoacido[2] == p_aminoacido[2] and pp_num_aminoacido[0:8] == p_num_aminoacido[0:8]:
                        # print("R" + p_atomo + "*" + p_nome_atomo + "*" + p_aminoacido + "*" + p_num_aminoacido + "*" + px + "*" + py + "*" + pz + "*" + p_enter) 
                        root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)
            
                if aa == 1:  # se não for o primeiro aminoácido
                    if p_nome_atomo[0] == 'C' and p_nome_atomo[1] == 'A' and p_nome_atomo[3] == ' ' and pp_aminoacido[0] == p_aminoacido[0] and pp_aminoacido[1] == p_aminoacido[1] and pp_aminoacido[2] == p_aminoacido[2] and pp_num_aminoacido[0] == p_num_aminoacido[0] and pp_num_aminoacido[1] == p_num_aminoacido[1] and pp_num_aminoacido[2] == p_num_aminoacido[2] and pp_num_aminoacido[3] == p_num_aminoacido[3] and pp_num_aminoacido[4] == p_num_aminoacido[4] and pp_num_aminoacido[5] == p_num_aminoacido[5] and pp_num_aminoacido[6] == p_num_aminoacido[6] and pp_num_aminoacido[7] == p_num_aminoacido[7]:
                        # print("R", p_atomo, "*", p_nome_atomo, "*", p_aminoacido, "*", p_num_aminoacido, "*", px, "*", py, "*", pz, "*", p_enter)
                        root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)
                        # print("1", p_atomo, "2", p_nome_atomo, "3", p_aminoacido, "4", p_num_aminoacido, "5", px, "6", py, "7", pz, "8", p_enter)
                        # print("1", pp_atomo, "2", pp_nome_atomo, "3", pp_aminoacido, "4", pp_num_aminoacido, "5", ppx, "6", ppy, "7", ppz, "8", pp_enter)

            # N-TODAS-PROTEINAS-EXCETO-PRO
            if pp_nome_atomo[0] == 'N' and (pp_aminoacido[0] != 'P' or pp_aminoacido[1] != 'R' or pp_aminoacido[2] != 'O'):
                # N
                if pp_nome_atomo[1] == ' ':
                    if aa != 1:  # se nao for o prmeiro aminoacido
                        num_p_num_aminoacido = int(p_num_aminoacido)
                        num_pp_num_aminoacido = int(pp_num_aminoacido)
                        num_pp_num_aminoacido_2 = num_pp_num_aminoacido - 1
                        if p_nome_atomo[0] == 'C' and p_nome_atomo[1] == ' ' and num_p_num_aminoacido == num_pp_num_aminoacido_2:
                            # cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter << endl;
                            root += f"R{p_atomo}{p_nome_atomo}{p_aminoacido}{p_num_aminoacido}{px}{py}{pz}{p_enter}"
                        # num_p_num_aminoacido == (num_pp_num_aminoacido - 1))
                        if (p_nome_atomo[0] == 'C' and p_nome_atomo[1] == 'A' and p_nome_atomo[2] == ' ' and 
                            pp_aminoacido[0] == p_aminoacido[0] and pp_aminoacido[1] == p_aminoacido[1] and 
                            pp_aminoacido[2] == p_aminoacido[2] and pp_num_aminoacido[0] == p_num_aminoacido[0] and 
                            pp_num_aminoacido[1] == p_num_aminoacido[1] and pp_num_aminoacido[2] == p_num_aminoacido[2] and 
                            pp_num_aminoacido[3] == p_num_aminoacido[3] and pp_num_aminoacido[4] == p_num_aminoacido[4] and 
                            pp_num_aminoacido[5] == p_num_aminoacido[5] and pp_num_aminoacido[6] == p_num_aminoacido[6] and 
                            pp_num_aminoacido[7] == p_num_aminoacido[7]):
                            # cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;
                            root += f"R{p_atomo}{p_nome_atomo}{p_aminoacido}{p_num_aminoacido}{px}{py}{pz}{p_enter}"
                    # if aa!=1
                    if aa == 1:  # se for o primeiro aminoáido
                        if (p_nome_atomo[0] == 'C' and p_nome_atomo[1] == 'A' and p_nome_atomo[2]

                # ----------------------------------NZ-LYS
                if pp_nome_atomo[1] == 'Z' and pp_nome_atomo[2] == ' ' and pp_aminoacido[0] == 'L' and pp_aminoacido[1] == 'Y' and pp_aminoacido[2] == 'S':
                    if p_nome_atomo[0] == 'C' and p_nome_atomo[1] == 'E' and p_nome_atomo[2] == ' ' and \
                        pp_aminoacido[0] == p_aminoacido[0] and pp_aminoacido[1] == p_aminoacido[1] and pp_aminoacido[2] == p_aminoacido[2] and \
                        pp_num_aminoacido[0] == p_num_aminoacido[0] and pp_num_aminoacido[1] == p_num_aminoacido[1] and pp_num_aminoacido[2] == p_num_aminoacido[2] and \
                        pp_num_aminoacido[3] == p_num_aminoacido[3] and pp_num_aminoacido[4] == p_num_aminoacido[4] and pp_num_aminoacido[5] == p_num_aminoacido[5] and \
                        pp_num_aminoacido[6] == p_num_aminoacido[6] and pp_num_aminoacido[7] == p_num_aminoacido[7]:
                            # print("R" + str(p_atomo) + "*" + p_nome_atomo + "*" + p_aminoacido + "*" + p_num_aminoacido + "*" + px + "*" + py + "*" + pz + "*" + p_enter)
                            root.write("R" + str(p_atomo) + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)
        
                #----------------------------------NE-ARG
                if pp_nome_atomo[1]=='E' and pp_nome_atomo[2]==' ' and pp_aminoacido[0]=='A' and pp_aminoacido[1]=='R' and pp_aminoacido[2]=='G':
                    if p_nome_atomo[0]=='C' and p_nome_atomo[1]=='D' and p_nome_atomo[2]==' ' and pp_aminoacido[0]==p_aminoacido[0] and pp_aminoacido[1]==p_aminoacido[1] and pp_aminoacido[2]==p_aminoacido[2] and pp_num_aminoacido[0]==p_num_aminoacido[0] and pp_num_aminoacido[1]==p_num_aminoacido[1] and pp_num_aminoacido[2]==p_num_aminoacido[2] and pp_num_aminoacido[3]==p_num_aminoacido[3] and pp_num_aminoacido[4]==p_num_aminoacido[4] and pp_num_aminoacido[5]==p_num_aminoacido[5] and pp_num_aminoacido[6]==p_num_aminoacido[6] and pp_num_aminoacido[7]==p_num_aminoacido[7]:
                        #print("R" + p_atomo + "*" + p_nome_atomo + "*" + p_aminoacido + "*" + p_num_aminoacido +  "*" + str(px) +  "*" + str(py) + "*" + str(pz) + "*" + p_enter)   
                        root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + str(px) + str(py) + str(pz) + p_enter)
                    if p_nome_atomo[0]=='C' and p_nome_atomo[1]=='Z' and p_nome_atomo[2]==' ' and pp_aminoacido[0]==p_aminoacido[0] and pp_aminoacido[1]==p_aminoacido[1] and pp_aminoacido[2]==p_aminoacido[2] and pp_num_aminoacido[0]==p_num_aminoacido[0] and pp_num_aminoacido[1]==p_num_aminoacido[1] and pp_num_aminoacido[2]==p_num_aminoacido[2] and pp_num_aminoacido[3]==p_num_aminoacido[3] and pp_num_aminoacido[4]==p_num_aminoacido[4] and pp_num_aminoacido[5]==p_num_aminoacido[5] and pp_num_aminoacido[6]==p_num_aminoacido[6] and pp_num_aminoacido[7]==p_num_aminoacido[7]:
                        #print("R" + p_atomo + "*" + p_nome_atomo + "*" + p_aminoacido + "*" + p_num_aminoacido +  "*" + str(px) +  "*" + str(py) + "*" + str(pz) + "*" + p_enter)   
                                        root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + str(px) + str(py) + str(pz) + p_enter)
                
                # ----------------------------------NH1/2-ARG
                if (pp_nome_atomo[1] == 'H' and (pp_nome_atomo[2] == '1' or pp_nome_atomo[2] == '2') and pp_aminoacido[0] == 'A' and pp_aminoacido[1] == 'R' and pp_aminoacido[2] == 'G'):
                    if (p_nome_atomo[0] == 'C' and p_nome_atomo[1] == 'Z' and p_nome_atomo[2] == ' ' and pp_aminoacido[0] == p_aminoacido[0] and pp_aminoacido[1] == p_aminoacido[1] and pp_aminoacido[2] == p_aminoacido[2] and pp_num_aminoacido[0] == p_num_aminoacido[0] and pp_num_aminoacido[1] == p_num_aminoacido[1] and pp_num_aminoacido[2] == p_num_aminoacido[2] and pp_num_aminoacido[3] == p_num_aminoacido[3] and pp_num_aminoacido[4] == p_num_aminoacido[4] and pp_num_aminoacido[5] == p_num_aminoacido[5] and pp_num_aminoacido[6] == p_num_aminoacido[6] and pp_num_aminoacido[7] == p_num_aminoacido[7]):
                        # cria arquivo temp com as raizes
                        root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)

                #----------------------------------ND1-HIS
                if (pp_nome_atomo[1]=='D' and pp_nome_atomo[2]=='1' and pp_aminoacido[0]=='H' and pp_aminoacido[1]=='I' and pp_aminoacido[2]=='S'):
                    if (p_nome_atomo[0]=='C' and p_nome_atomo[1]=='G' and p_nome_atomo[2]==' ' and pp_aminoacido[0]==p_aminoacido[0] and pp_aminoacido[1]==p_aminoacido[1] and pp_aminoacido[2]==p_aminoacido[2] and pp_num_aminoacido[0]==p_num_aminoacido[0] and pp_num_aminoacido[1]==p_num_aminoacido[1] and pp_num_aminoacido[2]==p_num_aminoacido[2] and pp_num_aminoacido[3]==p_num_aminoacido[3] and pp_num_aminoacido[4]==p_num_aminoacido[4] and pp_num_aminoacido[5]==p_num_aminoacido[5] and pp_num_aminoacido[6]==p_num_aminoacido[6] and pp_num_aminoacido[7]==p_num_aminoacido[7]):
                        #print("R" + p_atomo + "*" + p_nome_atomo + "*" + p_aminoacido + "*" + p_num_aminoacido + "*" + px + "*" + py + "*" + pz + "*" + p_enter)
                        root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)

                    if (p_nome_atomo[0]=='C' and p_nome_atomo[1]=='E' and p_nome_atomo[2]=='1' and pp_aminoacido[0]==p_aminoacido[0] and pp_aminoacido[1]==p_aminoacido[1] and pp_aminoacido[2]==p_aminoacido[2] and pp_num_aminoacido[0]==p_num_aminoacido[0] and pp_num_aminoacido[1]==p_num_aminoacido[1] and pp_num_aminoacido[2]==p_num_aminoacido[2] and pp_num_aminoacido[3]==p_num_aminoacido[3] and pp_num_aminoacido[4]==p_num_aminoacido[4] and pp_num_aminoacido[5]==p_num_aminoacido[5] and pp_num_aminoacido[6]==p_num_aminoacido[6] and pp_num_aminoacido[7]==p_num_aminoacido[7]):
                        #print("R" + p_atomo + "*" + p_nome_atomo + "*" + p_aminoacido + "*" + p_num_aminoacido + "*" + px + "*" + py + "*" + pz + "*" + p_enter)
                        root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)

                # ----------------------------------NE2-HIS
                if pp_nome_atomo[1] == 'E' and pp_nome_atomo[2] == '2' and pp_aminoacido[0] == 'H' and pp_aminoacido[1] == 'I' and pp_aminoacido[2] == 'S':
                    if p_nome_atomo[0] == 'C' and p_nome_atomo[1] == 'D' and p_nome_atomo[2] == '2' and pp_aminoacido[0] == p_aminoacido[0] and pp_aminoacido[1] == p_aminoacido[1] and pp_aminoacido[2] == p_aminoacido[2] and pp_num_aminoacido[0] == p_num_aminoacido[0] and pp_num_aminoacido[1] == p_num_aminoacido[1] and pp_num_aminoacido[2] == p_num_aminoacido[2] and pp_num_aminoacido[3] == p_num_aminoacido[3] and pp_num_aminoacido[4] == p_num_aminoacido[4] and pp_num_aminoacido[5] == p_num_aminoacido[5] and pp_num_aminoacido[6] == p_num_aminoacido[6] and pp_num_aminoacido[7] == p_num_aminoacido[7]:
                        # cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;
                        root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)
                    if p_nome_atomo[0] == 'C' and p_nome_atomo[1] == 'E' and p_nome_atomo[2] == '1' and pp_aminoacido[0] == p_aminoacido[0] and pp_aminoacido[1] == p_aminoacido[1] and pp_aminoacido[2] == p_aminoacido[2] and pp_num_aminoacido[0] == p_num_aminoacido[0] and pp_num_aminoacido[1] == p_num_aminoacido[1] and pp_num_aminoacido[2] == p_num_aminoacido[2] and pp_num_aminoacido[3] == p_num_aminoacido[3] and pp_num_aminoacido[4] == p_num_aminoacido[4] and pp_num_aminoacido[5] == p_num_aminoacido[5] and pp_num_aminoacido[6] == p_num_aminoacido[6] and pp_num_aminoacido[7] == p_num_aminoacido[7]:
                        # cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido <<  "*" << px <<  "*" << py << "*" << pz << "*" << p_enter;
                        root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)

                # ---------------------------ND2--ASN
                if (pp_nome_atomo[1]=='D' and pp_nome_atomo[2]=='2' and pp_aminoacido[0]=='A' and pp_aminoacido[1]=='S' and pp_aminoacido[2]=='N'):
                    if (p_nome_atomo[0]=='C' and p_nome_atomo[1]=='G' and p_nome_atomo[2]==' ' and pp_aminoacido[0]==p_aminoacido[0] and pp_aminoacido[1]==p_aminoacido[1] and pp_aminoacido[2]==p_aminoacido[2] and pp_num_aminoacido[0]==p_num_aminoacido[0] and pp_num_aminoacido[1]==p_num_aminoacido[1] and pp_num_aminoacido[2]==p_num_aminoacido[2] and pp_num_aminoacido[3]==p_num_aminoacido[3] and pp_num_aminoacido[4]==p_num_aminoacido[4] and pp_num_aminoacido[5]==p_num_aminoacido[5] and pp_num_aminoacido[6]==p_num_aminoacido[6] and pp_num_aminoacido[7]==p_num_aminoacido[7]):
                        root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)
                
                # ---------------------------NE1--GLN
                if (pp_nome_atomo[1]=='E' and pp_nome_atomo[2]=='2' and pp_aminoacido[0]=='G' and pp_aminoacido[1]=='L' and pp_aminoacido[2]=='N'):
                    if (p_nome_atomo[0]=='C' and p_nome_atomo[1]=='D' and p_nome_atomo[3]==' ' and pp_aminoacido[0]==p_aminoacido[0] and pp_aminoacido[1]==p_aminoacido[1] and pp_aminoacido[2]==p_aminoacido[2] and pp_num_aminoacido[0]==p_num_aminoacido[0] and pp_num_aminoacido[1]==p_num_aminoacido[1] and pp_num_aminoacido[2]==p_num_aminoacido[2] and pp_num_aminoacido[3]==p_num_aminoacido[3] and pp_num_aminoacido[4]==p_num_aminoacido[4] and pp_num_aminoacido[5]==p_num_aminoacido[5] and pp_num_aminoacido[6]==p_num_aminoacido[6] and pp_num_aminoacido[7]==p_num_aminoacido[7]):
                        root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)
                
                # ---------------------------NE1--TRP                   
                if pp_nome_atomo[1] == 'E' and pp_nome_atomo[2] == '1' and pp_aminoacido[0] == 'T' and pp_aminoacido[1] == 'R' and pp_aminoacido[2] == 'P':
                    if p_nome_atomo[0] == 'C' and p_nome_atomo[1] == 'D' and p_nome_atomo[2] == '1' and pp_aminoacido[0] == p_aminoacido[0] and pp_aminoacido[1] == p_aminoacido[1] and pp_aminoacido[2] == p_aminoacido[2] and pp_num_aminoacido[0] == p_num_aminoacido[0] and pp_num_aminoacido[1] == p_num_aminoacido[1] and pp_num_aminoacido[2] == p_num_aminoacido[2] and pp_num_aminoacido[3] == p_num_aminoacido[3] and pp_num_aminoacido[4] == p_num_aminoacido[4] and pp_num_aminoacido[5] == p_num_aminoacido[5] and pp_num_aminoacido[6] == p_num_aminoacido[6] and pp_num_aminoacido[7] == p_num_aminoacido[7]:
                        # print("R" + p_atomo + "*" + p_nome_atomo + "*" + p_aminoacido + "*" + p_num_aminoacido +  "*" + px +  "*" + py + "*" + pz + "*" + p_enter)
                        root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)
                    if p_nome_atomo[0] == 'C' and p_nome_atomo[1] == 'E' and p_nome_atomo[3] == '2' and pp_aminoacido[0] == p_aminoacido[0] and pp_aminoacido[1] == p_aminoacido[1] and pp_aminoacido[2] == p_aminoacido[2] and pp_num_aminoacido[0] == p_num_aminoacido[0] and pp_num_aminoacido[1] == p_num_aminoacido[1] and pp_num_aminoacido[2] == p_num_aminoacido[2] and pp_num_aminoacido[3] == p_num_aminoacido[3] and pp_num_aminoacido[4] == p_num_aminoacido[4] and pp_num_aminoacido[5] == p_num_aminoacido[5] and pp_num_aminoacido[6] == p_num_aminoacido[6] and pp_num_aminoacido[7] == p_num_aminoacido[7]:
                        # print("R" + p_atomo + "*" + p_nome_atomo + "*" + p_aminoacido + "*" + p_num_aminoacido +  "*" + px +  "*" + py + "*" + pz + "*" + p_enter)
                        root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)

                if pp_nome_atomo[0] == 'O':
                # ----------------------------------O
                    if pp_nome_atomo[1] == ' ':
                        if p_nome_atomo[0] == 'C' and p_nome_atomo[1] == ' ' and p_nome_atomo[3] == ' ' and pp_aminoacido[0] == p_aminoacido[0] and pp_aminoacido[1] == p_aminoacido[1] and pp_aminoacido[2] == p_aminoacido[2] and pp_num_aminoacido[0] == p_num_aminoacido[0] and pp_num_aminoacido[1] == p_num_aminoacido[1] and pp_num_aminoacido[2] == p_num_aminoacido[2] and pp_num_aminoacido[3] == p_num_aminoacido[3] and pp_num_aminoacido[4] == p_num_aminoacido[4] and pp_num_aminoacido[5] == p_num_aminoacido[5] and pp_num_aminoacido[6] == p_num_aminoacido[6] and pp_num_aminoacido[7] == p_num_aminoacido[7]:
                            # cout << "R" << p_atomo << "*" << p_nome_atomo << "*" << p_aminoacido << "*" << p_num_aminoacido << "*" << px << "*" << py << "*" << pz << "*" << p_enter;
                            root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)
                    # ---------------------------OD1/2--ASP
                    if pp_nome_atomo[1] == 'D' and (pp_nome_atomo[2] == '1' or pp_nome_atomo[2] == '2') and pp_aminoacido[0] == 'A' and pp_aminoacido[1] == 'S' and pp_aminoacido[2] == 'P':
                        if p_nome_atomo[0] == 'C' and p_nome_atomo[1] == 'G' and p_nome_atomo[3] == ' ' and pp_aminoacido[0] == p_aminoacido[0] and pp_aminoacido[1] == p_aminoacido[1] and pp_aminoacido[2] == p_aminoacido[2] and pp_num_aminoacido[0] == p_num_aminoacido[0] and pp_num_aminoacido[1] == p_num_aminoacido[1] and pp_num_aminoacido[2] == p_num_aminoacido[2] and pp_num_aminoacido[3] == p_num_aminoacido[3] and pp_num_aminoacido[4] == p_num_aminoacido[4] and pp_num_aminoacido[5] == p_num_aminoacido[5] and pp_num_aminoacido[6] == p_num_aminoacido[6] and pp_num_aminoacido[7] == p_num_aminoacido[7]:
                            root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)
                    if (pp_nome_atomo[1]=='E' and (pp_nome_atomo[2]=='1' or pp_nome_atomo[2]=='2') and pp_aminoacido[0]=='G' and pp_aminoacido[1]=='L' and pp_aminoacido[2]=='U'):
                            root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)
                    
                    if (p_nome_atomo[0]=='C' and p_nome_atomo[1]=='G' and p_nome_atomo[3]==' ' and pp_aminoacido[0]==p_aminoacido[0] and pp_aminoacido[1]==p_aminoacido[1] and pp_aminoacido[2]==p_aminoacido[2] and pp_num_aminoacido[0]==p_num_aminoacido[0] and pp_num_aminoacido[1]==p_num_aminoacido[1] and pp_num_aminoacido[2]==p_num_aminoacido[2] and pp_num_aminoacido[3]==p_num_aminoacido[3] and pp_num_aminoacido[4]==p_num_aminoacido[4] and pp_num_aminoacido[5]==p_num_aminoacido[5] and pp_num_aminoacido[6]==p_num_aminoacido[6] and pp_num_aminoacido[7]==p_num_aminoacido[7]):
                            root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)

                    # ---------------------------OD1--ASN    
                    if (pp_nome_atomo[1]=='D' and pp_nome_atomo[2]=='1' and pp_aminoacido[0]=='A' and pp_aminoacido[1]=='S' and pp_aminoacido[2]=='N'):
                        if (p_nome_atomo[0]=='C' and p_nome_atomo[1]=='G' and p_nome_atomo[2]==' ' and pp_aminoacido[0]==p_aminoacido[0] and pp_aminoacido[1]==p_aminoacido[1] and pp_aminoacido[2]==p_aminoacido[2] and pp_num_aminoacido[0]==p_num_aminoacido[0] and pp_num_aminoacido[1]==p_num_aminoacido[1] and pp_num_aminoacido[2]==p_num_aminoacido[2] and pp_num_aminoacido[3]==p_num_aminoacido[3] and pp_num_aminoacido[4]==p_num_aminoacido[4] and pp_num_aminoacido[5]==p_num_aminoacido[5] and pp_num_aminoacido[6]==p_num_aminoacido[6] and pp_num_aminoacido[7]==p_num_aminoacido[7]):
                            root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)

                    # ---------------------------OG--SER    
                    if (pp_nome_atomo[1]=='G' and pp_nome_atomo[2]==' ' and pp_aminoacido[0]=='S' and pp_aminoacido[1]=='E' and pp_aminoacido[2]=='
                            root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)
                    
                    #---------------------------OG1--THR	
                    if (pp_nome_atomo[1]=='G' && pp_nome_atomo[2]==' ' && pp_aminoacido[0]=='S' && pp_aminoacido[1]=='E' && pp_aminoacido[2]=='R'):
                        if (p_nome_atomo[0]=='C' and p_nome_atomo[1]=='B' and p_nome_atomo[3]==' ' and pp_aminoacido[0]==p_aminoacido[0] and pp_aminoacido[1]==p_aminoacido[1] and pp_aminoacido[2]==p_aminoacido[2] and pp_num_aminoacido[0]==p_num_aminoacido[0] and pp_num_aminoacido[1]==p_num_aminoacido[1] and pp_num_aminoacido[2]==p_num_aminoacido[2] and pp_num_aminoacido[3]==p_num_aminoacido[3] and pp_num_aminoacido[4]==p_num_aminoacido[4] and pp_num_aminoacido[5]==p_num_aminoacido[5] and pp_num_aminoacido[6]==p_num_aminoacido[6] and pp_num_aminoacido[7]==p_num_aminoacido[7]):
                            root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)
                    
                    # ---------------------------OH--TYR
                    if pp_nome_atomo[1] == 'G' and pp_nome_atomo[2] == '1' and pp_aminoacido[0] == 'T' and pp_aminoacido[1] == 'H' and pp_aminoacido[2] == 'R':
                        if p_nome_atomo[0] == 'C' and p_nome_atomo[1] == 'B' and p_nome_atomo[3] == ' ' and pp_aminoacido[0] == p_aminoacido[0] and pp_aminoacido[1] == p_aminoacido[1] and pp_aminoacido[2] == p_aminoacido[2] and pp_num_aminoacido[0] == p_num_aminoacido[0] and pp_num_aminoacido[1] == p_num_aminoacido[1] and pp_num_aminoacido[2] == p_num_aminoacido[2] and pp_num_aminoacido[3] == p_num_aminoacido[3] and pp_num_aminoacido[4] == p_num_aminoacido[4] and pp_num_aminoacido[5] == p_num_aminoacido[5] and pp_num_aminoacido[6] == p_num_aminoacido[6] and pp_num_aminoacido[7] == p_num_aminoacido[7]:
                            root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)
                    
                    # ---------------------------OE1--GLN
                    if pp_nome_atomo[1] == 'E' and pp_nome_atomo[2] == '1' and pp_aminoacido[0] == 'G' and pp_aminoacido[1] == 'L' and pp_aminoacido[2] == 'N':
                        if p_nome_atomo[0] == 'C' and p_nome_atomo[1] == 'D' and p_nome_atomo[3] == ' ' and pp_aminoacido[0] == p_aminoacido[0] and pp_aminoacido[1] == p_aminoacido[1] and pp_aminoacido[2] == p_aminoacido[2] and pp_num_aminoacido[0] == p_num_aminoacido[0] and pp_num_aminoacido[1] == p_num_aminoacido[1] and pp_num_aminoacido[2] == p_num_aminoacido[2] and pp_num_aminoacido[3] == p_num_aminoacido[3] and pp_num_aminoacido[4] == p_num_aminoacido[4] and pp_num_aminoacido[5] == p_num_aminoacido[5] and pp_num_aminoacido[6] == p_num_aminoacido[6] and pp_num_aminoacido[7] == p_num_aminoacido[7]:
                            root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)

                    # ---------------------------OXT--TODOS		
                    if pp_nome_atomo[1] == 'X' and pp_nome_atomo[2] == 'T':
                        if p_nome_atomo[0] == 'C' and p_nome_atomo[1] == ' ' and p_nome_atomo[3] == ' ' and                            pp_aminoacido[0] == p_aminoacido[0] and pp_aminoacido[1] == p_aminoacido[1] and                            pp_aminoacido[2] == p_aminoacido[2] and pp_num_aminoacido[0] == p_num_aminoacido[0] and                            pp_num_aminoacido[1] == p_num_aminoacido[1] and pp_num_aminoacido[2] == p_num_aminoacido[2] and                            pp_num_aminoacido[3] == p_num_aminoacido[3] and pp_num_aminoacido[4] == p_num_aminoacido[4] and                            pp_num_aminoacido[5] == p_num_aminoacido[5] and pp_num_aminoacido[6] == p_num_aminoacido[6] and                            pp_num_aminoacido[7] == p_num_aminoacido[7]:
                            root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)
                # -------------------------------------S
                # --------------------------------SD-MET
                if (pp_nome_atomo[0]=='S' and pp_nome_atomo[1]=='D' and pp_nome_atomo[2]==' ' and pp_aminoacido[0]=='M' and pp_aminoacido[1]=='E' and pp_aminoacido[2]=='T'):
                    if (p_nome_atomo[0]=='C' and p_nome_atomo[1]=='E' and p_nome_atomo[3]==' ' and pp_aminoacido[0]==p_aminoacido[0] and pp_aminoacido[1]==p_aminoacido[1] and pp_aminoacido[2]==p_aminoacido[2] and pp_num_aminoacido[0]==p_num_aminoacido[0] and pp_num_aminoacido[1]==p_num_aminoacido[1] and pp_num_aminoacido[2]==p_num_aminoacido[2] and pp_num_aminoacido[3]==p_num_aminoacido[3] and pp_num_aminoacido[4]==p_num_aminoacido[4] and pp_num_aminoacido[5]==p_num_aminoacido[5] and pp_num_aminoacido[6]==p_num_aminoacido[6] and pp_num_aminoacido[7]==p_num_aminoacido[7]):
                        root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)
                    if (p_nome_atomo[0]=='C' and p_nome_atomo[1]=='G' and p_nome_atomo[3]==' ' and pp_aminoacido[0]==p_aminoacido[0] and pp_aminoacido[1]==p_aminoacido[1] and pp_aminoacido[2]==p_aminoacido[2] and pp_num_aminoacido[0]==p_num_aminoacido[0] and pp_num_aminoacido[1]==p_num_aminoacido[1] and pp_num_aminoacido[2]==p_num_aminoacido[2] and pp_num_aminoacido[3]==p_num_aminoacido[3] and pp_num_aminoacido[4]==p_num_aminoacido[4] and pp_num_aminoacido[5]==p_num_aminoacido[5] and pp_num_aminoacido[6]==p_num_aminoacido[6] and pp_num_aminoacido[7]==p_num_aminoacido[7]):
                        root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)
        
                if pp_nome_atomo[0]=='S' and (pp_nome_atomo[1]==' ' or pp_nome_atomo[1]=='G') and pp_nome_atomo[2]==' ' and pp_aminoacido[0]=='C' and pp_aminoacido[1]=='Y' and pp_aminoacido[2]=='S':
                    if p_nome_atomo[0]=='C' and p_nome_atomo[1]=='B' and p_nome_atomo[3]==' ' and pp_aminoacido[0]==p_aminoacido[0] and pp_aminoacido[1]==p_aminoacido[1] and pp_aminoacido[2]==p_aminoacido[2] and pp_num_aminoacido[0]==p_num_aminoacido[0] and pp_num_aminoacido[1]==p_num_aminoacido[1] and pp_num_aminoacido[2]==p_num_aminoacido[2] and pp_num_aminoacido[3]==p_num_aminoacido[3] and pp_num_aminoacido[4]==p_num_aminoacido[4] and pp_num_aminoacido[5]==p_num_aminoacido[5] and pp_num_aminoacido[6]==p_num_aminoacido[6] and pp_num_aminoacido[7]==p_num_aminoacido[7]:
                        root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)
        root.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)        
        AD_root.close();
	AD.close();
	root.close();           


# In[ ]:


def P_centro_geom():
    # r1 = raiz1
    pos = 0
    conta_r = 0
    tipo = [2]
    p_atomo = [8]
    p_nome_atomo = [5]
    p_aminoacido = [4]
    p_num_aminoacido = [8]
    px = [9]
    py = [9]
    pz = [9]
    p_enter = [3]
    # raiz 1 
    r1_atomo = [9]
    r1_nome_atomo = [5]
    r1_aminoacido = [4]
    r1_num_aminoacido = [8]
    r1x = [9]
    r1y = [9]
    r1z = [9]
    r1_enter = [3]
    # raiz 2 
    r2_atomo = [9]
    r2_nome_atomo = [5]
    r2_aminoacido = [4]
    r2_num_aminoacido = [8]
    r2x = [9]
    r2y = [9]
    r2z = [9]
    r2_enter = [3]
    #--------converte caracteres para 2 raizes
    num_r1x = 0
    num_r1y = 0
    num_r1z = 0
    num_r2x = 0
    num_r2y = 0
    num_r2z = 0
    #--------converte caracteres para 3 raizes
    num_r3x = 0
    num_r3y = 0
    num_r3z = 0
    # raiz 3 
    r3_atomo = [9]
    r3_nome_atomo = [5]
    r3_aminoacido = [4]
    r3_num_aminoacido = [8]
    r3x = [9]
    r3y = [9]
    r3z = [9]
    r3_enter = [3]
    """
    # raiz 4 
    r4_atomo = [9]
    r4_nome_atomo = [5]
    r4_aminoacido = [4]
    r4_num_aminoacido = [8]
    r4x = [9]
    r4y = [9]
    r4z = [9]
    r4_enter = [3]
    """
    #--------valores para calcular o centro geometrico (geometric center)
    gcx2 = 0
    gcy2 = 0
    gcz2 = 0
    gcx2t = 0
    gcy2t = 0
    gcz2t = 0
    gcx3 = 0
    gcy3 = 0
    gcz3 = 0
    gcx3t = 0
    gcy3t = 0
    gcz3t = 0
    p_root_f = open('p_root_f.txt', 'w')  # final
    
    p_root_f = open("P_ROOT.tmp", "w")
    
    with open("P_ROOT_0.tmp", "r") as p_root:
        if not p_root:
            print("Could not create temp file!")

        while not p_root.eof():
            # pega uma linha
            get:
            tipo = p_root.read(sizeof(tipo)).decode('utf-8')
            p_atomo = p_root.read(sizeof(p_atomo)).decode('utf-8')
            p_nome_atomo = p_root.read(sizeof(p_nome_atomo)).decode('utf-8')
            p_aminoacido = p_root.read(sizeof(p_aminoacido)).decode('utf-8')
            p_num_aminoacido = p_root.read(sizeof(p_num_aminoacido)).decode('utf-8')
            px = p_root.read(sizeof(px)).decode('utf-8')
            py = p_root.read(sizeof(py)).decode('utf-8')
            pz = p_root.read(sizeof(pz)).decode('utf-8')
            p_enter = p_root.read(sizeof(p_enter)).decode('utf-8')
            # se nao for raiz e tiver so uma raiz
            # ... o resto do código aqui ...
            if tipo[0] == 'B' and conta_r == 1:
                conta_r = 0
                p_root_f.write("R" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + rx + ry + rz + p_enter)
                p_root_f.write("B" + p_atomo + p_nome_atomo + p_aminoacido + p_num_aminoacido + px + py + pz + p_enter)
                goto get # Ver isso
            
            # se a primeira linha for uma raiz
            if tipo[0] == 'R' and conta_r == 0:
                conta_r += 1
                #------copia do numero do atomo
                pos = 0
                while pos != 9:
                    r1_atomo[pos] = p_atomo[pos]
                    pos += 1
                #------copia do nome do atomo
                pos = 0
                while pos != 6:
                    r1_nome_atomo[pos] = p_nome_atomo[pos]
                    pos += 1
                #------copia do nome do aa
                pos = 0
                while pos != 5:
                    r1_aminoacido[pos] = p_aminoacido[pos]
                    pos += 1
                #------copia do numero do atomo
                pos = 0
                while pos != 9:
                    r1_num_aminoacido[pos] = p_num_aminoacido[pos]
                    pos += 1
                #------copia posicao x
                pos = 0
                while pos != 10:
                    r1x[pos] = px[pos]
                    pos += 1
                #------copia posicao y
                pos = 0
                while pos != 10:
                    r1y[pos] = py[pos]
                    pos += 1
                #------copia posicao z
                pos = 0
                while pos != 10:
                    r1z[pos] = pz[pos]
                    pos += 1
                #------copia final
                pos = 0
                while pos != 5:
                    r1_enter[pos] = p_enter[pos]
                    pos += 1
                pos = 0
                #print("*" + r1_atomo + "*" + r1_nome_atomo + "*" + r1_aminoacido + "*" + r1_num_aminoacido + "*" + r1x + "*" + r1y + "*" + r1z + "*" + r1_enter)
                #print(conta_r)
                get() # Ver isso

            # se a primeira linha for uma raiz  
            if tipo[0] == 'R' and conta_r == 1:
                conta_r += 1
                # ------copia do numero do atomo
                pos = 0
                while pos != 9:
                    r2_atomo[pos] = p_atomo[pos]
                    pos += 1
                # ------copia do nome do atomo
                pos = 0
                while pos != 6:
                    r2_nome_atomo[pos] = p_nome_atomo[pos]
                    pos += 1
                # ------copia do nome do aa
                pos = 0
                while pos != 5:
                    r2_aminoacido[pos] = p_aminoacido[pos]
                    pos += 1
                # ------copia do numero do atomo
                pos = 0
                while pos != 9:
                    r2_num_aminoacido[pos] = p_num_aminoacido[pos]
                    pos += 1
                # ------copia posicao x
                pos = 0
                while pos != 10:
                    r2x[pos] = px[pos]
                    pos += 1
                # ------copia posicao y
                pos = 0
                while pos != 10:
                    r2y[pos] = py[pos]
                    pos += 1
                # ------copia posicao z
                pos = 0
                while pos != 10:
                    r2z[pos] = pz[pos]
                    pos += 1
                # ------copia final
                pos = 0
                while pos != 5:
                    r2_enter[pos] = p_enter[pos]
                    pos += 1
                pos = 0
                # print("*" + r2_atomo + "*" + r2_nome_atomo + "*" + r2_aminoacido + "*" + r2_num_aminoacido + "*" + r2x + "*" + r2y + "*" + r2z + "*" + r2_enter)
                # print(conta_r)
                goto.get() # Ver isso
                
            if tipo[0] == 'B' and (conta_r == 2 or conta_r == 3):
                # elimina erros dos vetores da primeira raiz
                if r1x[0] == ' ':
                    r1x[0], r1x[1], r1x[2], r1x[3], r1x[4], r1x[5], r1x[6], r1x[7], r1x[8], r1x[9] = r1x[1], r1x[2], r1x[3], r1x[4], r1x[5], r1x[6], r1x[7], r1x[8], r1x[9], r1x[10]
                if r1y[0] == ' ':
                    r1y[0], r1y[1], r1y[2], r1y[3], r1y[4], r1y[5], r1y[6], r1y[7], r1y[8], r1y[9] = r1y[1], r1y[2], r1y[3], r1y[4], r1y[5], r1y[6], r1y[7], r1y[8], r1y[9], r1y[10]
                if r1z[0] == ' ':
                    r1z[0], r1z[1], r1z[2], r1z[3], r1z[4], r1z[5], r1z[6], r1z[7], r1z[8], r1z[9] = r1z[1], r1z[2], r1z[3], r1z[4], r1z[5], r1z[6], r1z[7], r1z[8], r1z[9], r1z[10]

                # elimina erros dos vetores da segunda raiz
                if r2x[0] == ' ':
                    r2x[0], r2x[1], r2x[2], r2x[3], r2x[4], r2x[5], r2x[6], r2x[7], r2x[8], r2x[9] = r2x[1], r2x[2], r2x[3], r2x[4], r2x[5], r2x[6], r2x[7], r2x[8], r2x[9], r2x[10]
                if r2y[0] == ' ':
                    r2y[0], r2y[1], r2y[2], r2y[3], r2y[4], r2y[5], r2y[6], r2y[7], r2y[8], r2y[9] = r2y[1], r2y[2], r2y[3], r2y[4], r2y[5], r2y[6], r2y[7], r2y[8], r2y[9], r2y[10]
                if r2z[0] == ' ':
                    r2z[0], r2z[1], r2z[2], r2z[3], r2z[4], r2z[5], r2z[6], r2z[7], r2z[8], r2z[9] = r2z[1], r2z[2], r2z[3], r2z[4], r2z[5], r2z[6], r2z[7], r2z[8], r2z[9], r2z[10]

                num_r1x = float(r1x)
                num_r1y = float(r1y)
                num_r1z = float(r1z)
                
                # print(num_r1x, num_r1y, num_r1z)
                num_r2x = float(r2x)
                num_r2y = float(r2y)
                num_r2z = float(r2z)
                
                # print(num_r2x, num_r2y, num_r2z)
                if conta_r == 2:
                    gcx2t = num_r1x + num_r2x
                    gcx2 = gcx2t / 2
                    gcy2t = num_r1y + num_r2y
                    gcy2 = gcy2t / 2
                    gcz2t = num_r1z + num_r2z
                    gcz2 = gcz2t / 2
                    # print("resultado(2): ", gcx3, "/", gcy3, "/", gcz3)
                    p_root_f.write("R   geom. center(2)   ")
                    letrasx = [' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']
                    letrasy = [' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']
                    letrasz = [' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ', ' ']
                    if -99.9999 <= gcx2 <= 999.9999:
                        letrasx = f"{gcx2:.4f}"
                    if -999.9999 <= gcx2 <= -100.0000:
                        letrasx = f"{gcx2:.3f}"
                    if -99.9999 <= gcy2 <= 999.9999:
                        letrasy = f"{gcy2:.4f}"
                    if -999.9999 <= gcy2 <= -100.0000:
                        letrasy = f"{gcy2:.3f}"
                    if -99.9999 <= gcz2 <= 999.9999:
                        letrasz = f"{gcz2:.4f}"
                    if -999.9999 <= gcz2 <= -100.0000:
                        letrasz = f"{gcz2:.3f}"
                    # ---------------X
                    if letrasx[8] == ' ':
                        letrasx[0:8] = letrasx[1:9]
                        letrasx[9] = ' '
                    if letrasx[8] == ' ':
                        letrasx[0:8] = letrasx[1:9]
                        letrasx[9] = ' '
                    if letrasx[8] == ' ':
                        letrasx[0:8] = letrasx[1:9]
                        letrasx[9] = ' '
                    # ---------------Y
                    if letrasy[8] == ' ':
                        letrasy[0:8] = letrasy[1:9]
                        letrasy[9] = ' '
                    if letrasy[8] == ' ':
                        letrasy[0:8] = letrasy[1:9]
                        letrasy[9] = ' '
                    if letrasy[8] == ' ':
                        letrasy[0:8] = letrasy[1:9]
                        letrasy[9] = ' '
                    # ---------------Z
                    if letrasz[8] == ' ':
                        letrasz[0:8]
                        
                    if letrasz[8] == ' ':
                        letrasz[8] = letrasz[7]
                        letrasz[7] = letrasz[6]
                        letrasz[6] = letrasz[5]
                        letrasz[5] = letrasz[4]
                        letrasz[4] = letrasz[3]
                        letrasz[3] = letrasz[2]
                        letrasz[2] = letrasz[1]
                        letrasz[1] = letrasz[0]
                        letrasz[0] = ' '
                    if letrasz[8] == ' ':
                        letrasz[8] = letrasz[7]
                        letrasz[7] = letrasz[6]
                        letrasz[6] = letrasz[5]
                        letrasz[5] = letrasz[4]
                        letrasz[4] = letrasz[3]
                        letrasz[3] = letrasz[2]
                        letrasz[2] = letrasz[1]
                        letrasz[1] = letrasz[0]
                        letrasz[0] = ' '
                    # cout << letrasx << endl;
                    # cout << letrasx[0] << letrasx[1] << letrasx[2] << letrasx[3] << letrasx[4] << letrasx[5] << letrasx[6] << letrasx[7];
                    # p_root_f << letrasx[0] << letrasx[1] << letrasx[2] << letrasx[3] << letrasx[4] << letrasx[5] << letrasx[6] << letrasx[7];
                    p_root_f.write(letrasx + letrasy + letrasz + ' ')
                    p_root_f.write('\n')
                    conta_r = 0
                    
            if conta_r == 3:
            # elimina erros dos vetores da segunda raiz
            if r3x[0] == ' ':
                r3x[0], r3x[1], r3x[2], r3x[3], r3x[4], r3x[5], r3x[6], r3x[7], r3x[8], r3x[9] =                     r3x[1], r3x[2], r3x[3], r3x[4], r3x[5], r3x[6], r3x[7], r3x[8], r3x[9], r3x[10]

            if r3y[0] == ' ':
                r3y[0], r3y[1], r3y[2], r3y[3], r3y[4], r3y[5], r3y[6], r3y[7], r3y[8], r3y[9] =                     r3y[1], r3y[2], r3y[3], r3y[4], r3y[5], r3y[6], r3y[7], r3y[8], r3y[9], r3y[10]

            if r3z[0] == ' ':
                r3z[0], r3z[1], r3z[2], r3z[3], r3z[4], r3z[5], r3z[6], r3z[7], r3z[8], r3z[9] =                     r3z[1], r3z[2], r3z[3], r3z[4], r3z[5], r3z[6], r3z[7], r3z[8], r3z[9], r3z[10]




# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:


def calcula_RT():
    with open("limit_type.tmp") as limit_type:
        line_count = 0
        for line in limit_type:
            line_count += 1
        
        atom_0 = line_count[8]
        atom_1 = line_count[8]
        total_atom_l = line_count
        
        limit_type.seek(0)
        current_line = 0
        for line in limit_type:
            atom_1[current_line] = line.strip()
            current_line += 1

    with open("limit_type.tmp") as limit_type:
        line_count = 0
        for line in limit_type:
            line_count += 1
        atom_0 = [""] * line_count
        atom_1 = [""] * line_count
        limit_type.seek(0)
        current_line = 0
        for line in limit_type:
            atom_1[current_line] = line.strip()
            current_line += 1

    with open("RT_found.tmp") as RT_found:
        rt_string = []
        for line in RT_found:
            rt_string.append(line.strip())

    aux = rt_string.copy()
    for i in range(lig_line):
        for j in range(lig_line):
            if rt_string[i] == aux[j] and i != j:
                aux[i] = "REP"

    lines_total2 = -1
    lines_total2 -= 1

    limit_type.close()
    lig_line = 0
    RT_found = open("RT_found.tmp")
    rt_line = bytearray(8)
    while True:
        line = RT_found.read(8)
        if not line:
            break
        rt_line[:len(line)] = line
        lig_line += 1
    total_lig_line = lig_line
    rt_string = bytearray(total_lig_line*8)
    lig_line = 0
    RT_found.seek(0)
    while True:
        line = RT_found.read(8)
        if not line:
            break
        rt_string[lig_line*8:(lig_line+1)*8] = line
        lig_line += 1
    i = 0
    j = 0
    aux = rt_string.copy()
    for i in range(lig_line):
        for j in range(lig_line):
            if i != j and rt_string[i*8:(i+1)*8] == aux[j*8:(j+1)*8]:
                aux[i*8:(i+1)*8] = b'REP'
    lines_total2 = 0
    lines_total2 -= 1
    j = 0
    rt_string_copy = bytearray(total_lig_line*8)
    for i in range(lig_line):
        if aux[i*8:(i+1)*8] != b'REP':
            rt_string_copy[j*8:(j+1)*8] = aux[i*8:(i+1)*8]
            j += 1
            lines_total2 += 1
    current_line_0 = 0
    current_line_1 = 0
    lig_line = 0
    ligand = open("ligand.tmp")
    ligand_line = bytearray(55)
    atom_type = bytearray(line_count//2*12)
    while True:
        line = ligand.read(55)
        if not line:
            break
        ligand_line[:len(line)] = line
        atom_type[lig_line*12:(lig_line+1)*12] = (ligand_line[1:8] + ligand_line[47:52])
        lig_line += 1
    lig_line = 0
    ligand.close()
    mol2 = open("ligand.tmp")
    mol2_line = mol2_line[54]
    while True:
        line = mol2.read(54)
        if not line:
            break
        mol2_line[:len(line)] = line
        type_l += 1
    mol2.clear()
    mol2.seek(0)
    mol2_type = mol2_type[type_l][10]
    type_l = 0
    
    mol2 = open("file.mol2", "r")
    mol2_lines = mol2.readlines()
    mol2.close()

    mol2_type = [[''] * 10 for i in range(len(mol2_lines))]
    type_l = 0
    for mol2_line in mol2_lines:
        mol2_type[type_l][0] = mol2_line[0]
        mol2_type[type_l][1] = mol2_line[1]
        mol2_type[type_l][2] = mol2_line[2]
        mol2_type[type_l][3] = mol2_line[3]
        mol2_type[type_l][4] = mol2_line[4]
        mol2_type[type_l][5] = mol2_line[5]
        mol2_type[type_l][6] = mol2_line[6]
        mol2_type[type_l][7] = mol2_line[47]
        mol2_type[type_l][8] = mol2_line[48]
        mol2_type[type_l][9] = mol2_line[49]
        type_l += 1

    sp2_count = 0
    atom_count = 0
    marker = 0
    pos = 0
    pair_count1 = 0
    pair_count2 = 1
    num_atom1 = 0
    num_atom2 = 0
    atom1 = ''
    atom2 = ''
    number_count = 0
    pair = 0
    next = 0
    while current_line_0 < lines_total2:
        if marker < line_count and current_line_0 < lines_total2:
            while marker < line_count:
                number_count += 1
                next = marker + 1
                if rt_string[current_line_0] == atom_1[marker][:6]:
                    atom_count += 1
                    if atom_1[next][5] == '1':
                        atom1 = atom_1[marker][:6]
                        num_atom1 = int(atom1)
                        if mol2_type[num_atom1][9] == '2':
                            sp2_count += 1
                        if number_count % 2 != 0:
                            pair = 2
                        if number_count % 2 == 0:
                            pair = -2
                        atom2 = atom_1[marker + pair][:6]
                        num_atom2 = int(atom2)
                        if mol2_type[num_atom2][9] == '2':
                            sp2_count += 1
                        if sp2_count < 2 and mol2_type[num_atom2][7] != 'H' and mol2_type[num_atom1][7] != 'H':
                            s_bond += 1
                marker += 2
        marker = 0

        if atom_count > 1:
            if s_bond == 0:
                RT += 0
            elif s_bond == 1 or s_bond >= 3:
                RT += 0.5
            elif s_bond == 2:
                RT += 1

        s_bond = 0
        marker = 0
        sp2_count = 0
        atom_count = 0


# In[42]:


def open_protein():
    atoms = structure.get_atoms()
    residues = structure.get_residues()
    chains = structure.get_chains()
    
    return atoms,residues,chains


# In[44]:


def result_score_calc():
    import os
    import subprocess

    HB = 0  # HB = hydrogen bonds = C
    HB = refined
    result_file_name = ""
    v = []
    v.clear()
    v = ligand_name.split(".")
    result_file_name = v[0] + "_result.txt"
    result_score = open(result_file_name, "a")

    result_score.write("Hydrophobic contacts:        " + str(HC_total2) + "\n")
    result_score.write("Van der waals:               " + str(VDW_total) + "\n")
    result_score.write("Deformation effect:          " + str(RT) + "\n")
    result_score.write("Hydrogen bonds (HB):         " + str(hydrogen_B) + "\n")
    result_score.write("Repulsive VDW score:         " + str(repulsive) + "\n")
    result_score.write("London dispersion force:     " + str(london) + "\n")

    test = 1
    test2 = 5
    result = summation(test, test2)

    ASA1 = 0
    ASA2 = 0
    sendme = ""
    result1 = ""
    result2 = ""
    result3 = ""
    str1 = "python asa/asa.py " + protein_name
    str2 = "python asa/asa.py " + ligand_pdb
    fp1 = subprocess.Popen(str1, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result1 = fp1.communicate()[0]
    fp1.stdout.close()
    fp1.stderr.close()
    fp2 = subprocess.Popen(str2, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    result2 = fp2.communicate()[0]
    fp2.stdout.close()
    fp2.stderr.close()

    ASA1 = float(result1)
    ASA2 = float(result2)
    str3 = ligand_name + ".interaction_terms.txt"
    type = ""
    if ligand_type == "DNA" or ligand_type == "RNA" or ligand_type == "DNA/RNA" or ligand_type == "nucleotide" or ligand_type == "protein-DNA" or ligand_type == "protein-RNA" or ligand_type == "protein-DNA/RNA" or ligand_type == "protein-nucleotide" or ligand_type == "nucl":
        str4 = "Rscript R/DNA/script.R "
        type = "DNA"
        # print("found DNA")
    else:
        if ligand_type == "small" or ligand_type == "small molecule" or ligand_type == "drug" or ligand_type == "prtoein-small" or ligand_type == "protein-small molecule" or ligand_type == "protein-drug":
            str4 = "Rscript R/small_molecule/script.R "
            # print("found small molecule")
            type = "small"
        else:
            if ligand_type == "protein" or ligand_type == "protein-protein" or ligand_type == "peptide" or ligand_type == "polypetide" or ligand_type == "prot":
                str4 = "Rscript R/protein/script.R "
                # print("found protein")
                type = "protein"
            else:
                str4 = "Rscript R/DNA/script.R "
                if ligand_type == "unknown":
                    print("No ligand type provided, using default scoring function (protein-DNA/RNA).\n")
                    type = "DNA"
                else:
                    print("Found unknown ligand type parameter " + ligand_type + ", using default (DNA).\n")
                    type = "DNA"

    # start: calculates surface tension and hydrophobic
    # start: calculates surface tension and hydrophobicity
    hydro_file = "hydrophobicity.param"
    tension_file = "tension.param"
    hydro_map = {}
    tension_map = {}
    with open(hydro_file, 'r') as f1, open(tension_file, 'r') as f2:
        for line in f1:
            k, v = line.split()
            hydro_map[k] = float(v)
        for line in f2:
            k, v = line.split()
            tension_map[k] = float(v)

    protein_file_name = "protein.tmp"
    p_dist_file_name = "P_DIST.tmp"
    total_hydrophobicity = 0
    total_surface_tension = 0
    contact_hydrophobicity = 0
    contact_surface_tension = 0
    hydrophobicity = 0
    surface_tension = 0

    surface_tension_hydrophobicity_calculator(hydro_map, tension_map, protein_file_name)
    total_hydrophobicity = hydrophobicity
    total_surface_tension = surface_tension

    hydrophobicity = 0
    surface_tension = 0
    surface_tension_hydrophobicity_calculator(hydro_map, tension_map, p_dist_file_name)
    contact_hydrophobicity = hydrophobicity
    contact_surface_tension = surface_tension

    # end: calculates surface tension and hydrophobicity
    with open(str3, 'w') as interaction_terms:
        if type == "DNA":
            # DNA
            interaction_terms.write("V2\tV3\tV4\tV5\tV6\tV7\tV18\tV19\tV20\tV21\tV22\tV23\n")
            interaction_terms.write(f"{HC_total2}\t{VDW_total}\t{RT}\t{hydrogen_B}\t{ASA1}\t{ASA2}\t{repulsive}\t{london}\t{contact_hydrophobicity}\t{total_hydrophobicity}\t{contact_surface_tension}\t{total_surface_tension}\n")
        elif type == "small":
            # small
            interaction_terms.write("V2\tV3\tV4\tV5\tV6\tV7\tV17\tV18\tV20\tV21\tV22\tV23\n")
            interaction_terms.write(f"{HC_total2}\t{VDW_total}\t{RT}\t{hydrogen_B}\t{ASA1}\t{ASA2}\t{repulsive}\t{london}\t{contact_hydrophobicity}\t{total_hydrophobicity}\t{contact_surface_tension}\t{total_surface_tension}\n")
        elif type == "protein":
            # protein
            interaction_terms.write("V2\tV3\tV4\tV5\tV6\tV7\tV16\tV17\tV20\tV21\tV22\tV23\n")
            interaction_terms.write(f"{HC_total2}\t{VDW_total}\t{RT}\t{hydrogen_B}\t{ASA1}\t{ASA2}\t{repulsive}\t{london}\t{contact_hydrophobicity}\t{total_hydrophobicity}\t{contact_surface_tension}\t{total_surface_tension}\n")

    # execute shell command and get output
    os.system(f"{str4} {str3} > result.txt")
    with open('result.txt', 'r') as f:
        result3 = f.read()
    print(float(result3))
    
    


# In[ ]:


import sys

def main(argv):
    if len(argv) < 4:
        print("usage: GLM-Score <protein_file PDB> <ligand_file PDB> <ligand_file MOL2> <Ligand Type> \n\n<Ligand Type> is the type of ligand molecule, which supports 3 alternative options:\n\nDNA (or RNA), small-molecule (or drug), and protein.\nIf no <Ligand Type> is provided, protein-DNA/RNA scoring function will be used as default.\n\n")
        sys.exit(1)
    
    protein_name = argv[1]
    ligand_pdb = argv[2]
    ligand_name = argv[3]
    
    if len(argv) == 4:
        ligand_type = "unknown"
    elif len(argv) == 5:
        ligand_type = argv[4]

    # 1) processa entreadas
    # [JC] Open the file ligand.mol2
    open_ligand()
    # [JC] Open PDF file
    open_protein()

    # 2) Calculos
    # 2.1 [JC] Hydrophobic contacts (HCs) / van der Waals interactions (VDWs) / Repulsive interactions (RIs)	
    distancia()

    # 2.2 [JC]
    raiz_proteina()
    raiz_ligante()
    P_centro_geom()
    L_centro_geom()
    angulos()

    # 3) Saida
    salva_ligante()
    saida_PDB()
    salva_proteina()

    # 4) refinamento do HB
    if hydrogen_B != 0:
        limit()
    calcula_RT()
    result_score_calc()
    deleta_temp()
    # system("pause")

if __name__ == "__main__":
    main(sys.argv[1:])

