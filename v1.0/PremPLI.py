#!/usr/bin/python
# coding=utf-8
import os, sys, os.path, re, getopt, time
from string import ascii_uppercase
from string import ascii_lowercase
from collections import defaultdict,Counter
import pandas as pd
import rpy2.robjects as robjects
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import Selection
from sklearn.neighbors import KDTree
r = robjects.r
r('''library(randomForest)''')
ascii_cases = ascii_uppercase + ascii_lowercase

path = Your working directory
pathvmd = path for running VMD software  # /usr/local/bin/vmd
pathmkdssp = path for running DSSP software  # /usr/local/bin/mkdssp
pathpsiblast = path for running PSI-BLAST software  # /usr/local/bin/blast/psiblast
pathblastdb = path for blastdb  # /usr/local/bin/blastdb/nr
pathprovean = path for PROVEAN software  # /usr/bin/provean.sh
patharpeggio = path for Arpeggio software  # /usr/local/bin/arpeggio/arpeggio.py
pathxlogp3 = path for running XLOGP3 software # /usr/local/bin/
pathinput = path+ 'inputfiles'

# set up input and output file name for each job
jobid = ''
myopts, args = getopt.getopt(sys.argv[1:], "i:")
for o, a in myopts:
    if o == '-i':
        jobid = a
    else:
        print ("Usage: %s -i jobid" % sys.argv[0])

jobpath = path + jobid
pathoutput = path + jobid + "out"  # can change
os.system('mkdir %s'%(pathoutput))
in_file = jobpath + '/' + jobid + '.input'
out_file = jobpath + '/' + jobid + '.sunddg'
os.chdir(path)

# map residue name three letters to one
map_three_one = {"GLY": "G", "ALA": "A", "SER": "S", "THR": "T", "CYS": "C",
                 "VAL": "V", "LEU": "L", "ILE": "I", "MET": "M", "PRO": "P",
                 "PHE": "F", "TYR": "Y", "TRP": "W", "ASP": "D", "GLU": "E",
                 "ASN": "N", "GLN": "Q", "HIS": "H", "LYS": "K", "ARG": "R",
                 "ASX": "X", "GLX": "X", "CSO": "X", "HIP": "X", "MSE": "X",
                 "UNK": "X", "SEC": "X", "PYL": "X", "SEP": "X", "TPO": "X",
                 "PTR": "X", "XLE": "X", "XAA": "X", "HSD": "H", "HID": "H",
                 "HSE": "H"}

normal_format_pro = ['CYS','GLN','ILE','SER','VAL','MET','ASN','PRO','LYS','THR','PHE','ALA','HIS','GLY','ASP','LEU','ARG','TRP','GLU','TYR']

# standard state amino acid suface areas.REF: Hydrophobicity of amino acid residues in globular proteins
map_surface = {'A':118.1,'R':256.0,'N':165.5,'D':158.7,'C':146.1,'Q':193.2,'E':186.2,'G':88.1,'H':202.5,'I':181.0,'L':193.1,'K':225.8,'M':203.4,'F':222.8,'P':146.8,'S':129.8,'T':152.5,'W':266.3,'Y':236.8,'V':164.5,'X':88.1}


def ProPDB1():
    pdball = []
    f = open(in_file, 'r')
    f.next()
    for line in f:
        ff = line.split("\t")
        pdb = ff[0].split('.')[0]
        suffix = ff[0].split('.')[1].strip()
        ligand = ff[6]
        lig_pos = ff[7].split("\n")[0].split('.')
        if pdb not in pdball:
            pdball.append(pdb)
            ffpdb = open(pathoutput + '/' + jobid + '_p.pdb', 'w')
            try:
                fpdb = open(jobpath + '/' + pdb.upper() + '.' + suffix, 'r')
            except:
                fpdb = open(jobpath + '/' + pdb.lower() + '.' + suffix, 'r')
            if suffix != 'pdb':
                ST_play = False
                for linepdb in fpdb:
                    if linepdb[0:5] == "MODEL":
                        CountModel = linepdb.split()[1]
                        ST_play = True
                        continue
                    if ST_play:
                        if (linepdb[:4] == "ATOM") and (linepdb[17:20].strip() in normal_format_pro):
                            ffpdb.write("%s                 %s\n" % (linepdb[0:55].strip('\r\n'), str(linepdb[21:22]) + '_' + str(CountModel)))
                        if (linepdb[0:6] == "HETATM") and (ligand in linepdb[17:20]):
                            flag = str(linepdb[21:22]) + '_' + str(CountModel)+'__'+str(linepdb[22:27].strip())
                            if flag in lig_pos:
                                ffpdb.write("HETATM%s%s                 %s\n" % (linepdb[6:20], linepdb[20:55].strip('\r\n'), str(linepdb[21:22]) + '_' + str(CountModel)))
            else:
                ST_play = True 
                for linepdb in fpdb:
                    line_list = re.split(r'\s+', linepdb)
                    if (line_list[0] == 'MODEL') and ST_play:
                        countmodel = line_list[1]
                        ST_play = False
                    if (line_list[0] == 'MODEL') and (line_list[1] != countmodel):
                        break
                    if (linepdb[:4] == "ATOM") and (linepdb[17:20].strip() in normal_format_pro):
                        ffpdb.write("%s                 %s\n" % (
                        linepdb[0:55].strip('\r\n'), str(linepdb[21:22]) + '_' + str(1)))
                    if (linepdb[0:6] == "HETATM") and (ligand in linepdb[17:20]):
                        flag = str(linepdb[21:22]) + '_' + str(1) + '__' +str(linepdb[22:27].strip())
                        if flag in lig_pos:
                            ffpdb.write("HETATM%s%s                 %s\n" % (linepdb[6:20], linepdb[20:55].strip('\r\n'), str(linepdb[21:22]) + '_' + str(1)))
            fpdb.close()
            ffpdb.close()
        else:
            continue
    f.close()


def del_unknown_incomplete():
    residues = [k for k, v in map_three_one.iteritems() if v != 'X']
    f = open(in_file, 'r')
    f.next()
    for line in f:
        ff = line.split("\t")
        pdbid = ff[0].split('.')[0].lower()
        ligand = ff[6]
        pdbfile = open(pathoutput + '/' + jobid + '_p.pdb').readlines()
        # delete unknow resdues
        pdbfile_del_unknown = [i for i in pdbfile if i[17:20] not in ['UNK']]
        # delete incomplete residues
        final_row = pdbfile_del_unknown[-1]
        last = ''
        above = []
        allresidues = []
        for row in pdbfile_del_unknown:
            if row[17:26] == last and row == final_row:  # when read final row，append it if equal to last,MET A   1
                above.append(row)
                atoms = [i[13:16].strip() for i in above]  # C27
                if set(['C', 'N', 'O', 'CA']).issubset(set(atoms)) and ('ATOM' in [i[0:6].strip() for i in above]):  # (\s)ALA
                    allresidues.append(above)
                elif (ligand in [i[17:20].strip() for i in above]) and ('HETATM' in [i[0:6].strip() for i in above]):
                    allresidues.append(above)
            elif row[17:26] == last and row != final_row:  # when read same residue, but not last row
                above.append(row)
            else:  # when read different residue
                if len(above) >= 4:
                    atoms = [i[13:16].strip() for i in above]
                    if (set(['C', 'N', 'O', 'CA']).issubset(set(atoms))) and ('ATOM' in [i[0:6].strip() for i in above]):
                        allresidues.append(above)
                    elif (ligand in [i[17:20].strip() for i in above]) and ('HETATM' in [i[0:6].strip() for i in above]):
                        allresidues.append(above)
                above = [row]
            last = row[17:26]
        # write out
        with open(pathoutput + '/' + jobid + '_p_test.pdb', 'w') as fw:
            fw.write(''.join([y for x in allresidues for y in x]))
        break
    os.system('mv %s/%s_p_test.pdb %s/%s_p.pdb' % (pathoutput, jobid, pathoutput, jobid))
    f.close()


#split pdbfile by chains
def splitchain():
    pdball = []
    with open(in_file, 'r') as f:
        f.next()
        for line in f:
            ff = line.split('\t')
            pdbfile = ff[0]  # 1A43.pdb
            pdb = pdbfile.split(".")[0].upper()
            partner1 = ff[1].split(".")
            partner2 = ff[2].split(".")
            mutachain = ff[3].split("_")
            if pdbfile not in pdball:
                for chain in list(partner1 + partner2):
                    f1 = open(pathoutput + '/' + jobid + '_p.pdb', 'r')
                    fw = open(pathoutput + '/' + pdb + '_' + chain + '.pdb', 'w')
                    for line1 in f1:
                        if chain == line1[72:].strip():
                            fw.write(line1)
                    f1.close()
                    fw.close()
                pdball.append(pdbfile)
            else:
                continue
    f.close()


#split ligand by chains
def splitligand():
    f = open(in_file, 'r')
    _unused = f.next()
    for line in f:
        ff = line.split("\t")
        pdbfile = ff[0]
        pdb = pdbfile.split(".")[0].upper()
        partner1 = ff[1].split(".")
        partner2 = ff[2].split(".")
        ligand = ff[6]
        for chain in list(partner1):
            fpdb = open(pathoutput + '/' + pdb + '_' + chain + '.pdb', 'r')
            fw = open(pathoutput + '/' + pdb + '_PROTEIN_' + chain + '.pdb', 'w')
            for line1 in fpdb:
                if line1[0:6] != 'HETATM':
                    fw.write(line1)
                else:
                    continue
            fpdb.close()
            fw.close()
        for chain_lig in list(partner2):
            fpdb = open(pathoutput + '/' + pdb + '_' + chain_lig + '.pdb', 'r')
            fp = open(pathoutput + '/' + pdb + '_' + ligand + '_' + chain_lig + '.pdb', 'w')
            for line2 in fpdb:
                if line2[0:6] == 'HETATM':
                    fp.write(line2)
                else:
                    continue
            fpdb.close()
            fp.close()

    f.close()


#renumber chain and position
def CleanPdb():
    need = [k for k, v in map_three_one.items() if v != 'X']
    first_line = open(in_file, 'r').readlines()[0][:-1]
    fw = open(in_file + ".cleaned", "w")
    fw.write("%s\t%s\t%s\t%s\t%s\n" % (first_line, "PDBid", "NewPartner1", "NewPartner2", "Mutation_cleaned"))
    second_line = open(in_file, 'r').readlines()[1]  
    ff = second_line.split("\t")
    pdb = ff[0].split(".")[0]  # 1AAY
    partner1 = ff[1].split(".")  # [A_1]
    partner2 = ff[2].split(".")  # [B_1,C_1]
    ligand = ff[6]
    mapchainarray = []
    counti = 0
    for chains in list(partner1 + partner2):
        if chains not in dict(iter(mapchainarray)).keys():
            cc = (chains, ascii_cases[counti])
            mapchainarray.append(cc)
            counti += 1
            mapchaindict = dict(iter(mapchainarray))
        else:
            continue
    same_chain = []
    choice = [i for i in ascii_uppercase if i not in mapchaindict.values()]
    for chain_new in list(partner2):
        if chain_new in list(partner1):
            random1 = choice.pop()
            mapchainarray.append((random1 + '_1', random1))
            mapchaindict = dict(iter(mapchainarray))
            same_chain.append((chain_new, random1 + '_1'))
            same_chain_dict = dict(iter(same_chain))
        else:
            continue

    newpartner1 = ''
    for chains in list(partner1):
        newpartner1 += mapchaindict[chains]
    newpartner2 = ''
    for chain_new in list(partner2):
        if chain_new not in list(partner1):
            newpartner2 += mapchaindict[chain_new]
        else:
            newpartner2 += mapchaindict[same_chain_dict[chain_new]]
    countchain = 1
    for chains in list(partner1):
        fwpdb = open(pathoutput + "/" + pdb + "_PROTEIN_CH" + str(countchain) + ".pdb", "w")
        fvar = open(pathoutput + "/" + pdb + "_" + chains + ".var", "w")
        countchain += 1
        count = 1
        fpdb = open(pathoutput + "/" + pdb + "_PROTEIN_" + chains + ".pdb", "r")
        resname = fpdb.readlines()[0][17:20].strip()
        fpdb = open(pathoutput + "/" + pdb + "_PROTEIN_" + chains + ".pdb", "r")
        resnum = fpdb.readlines()[0][22:27].strip()
        line1 = open(pathoutput + "/" + pdb + "_PROTEIN_" + chains + ".pdb", "r").readlines()[0]
        fpdb = open(pathoutput + "/" + pdb + "_PROTEIN_" + chains + ".pdb", "r")
        atomname = fpdb.readlines()[0][13:16].strip()
        fwpdb.write("%s %s%s   %s %s" % (line1[0:16], line1[17:21], mapchaindict[chains], str(count), line1[27:]))

        mutchainall = []
        f = open(in_file, 'r')
        _unused = f.next()
        for line in f:
            ff = line.split("\t")
            partner2 = ff[2].split('.')
            mut = ff[4].upper()
            mutchain = ff[3]
            if str(chains) == str(mutchain) and mutchain not in mutchainall:
                if resname in need:
                    fseq = open(pathoutput + "/" + pdb + "_" + mutchain + ".seq", "w")
                    fseq.write("%s %s\n" % (">", pdb + mutchain))
                    fseq.write("%s" % (map_three_one[resname]))
                else:
                    continue
            mutchainall.append(mutchain)
            if (str(chains) == str(mutchain)):
                if (str(resnum) == str(mut[1:-1])) and (str(map_three_one[resname]) == str(mut[0:1])):
                    fw.write("%s\t%s\t%s\t%s\t%s\n" % (line[:-1], pdb, newpartner1, newpartner2, str(map_three_one[resname]) + mapchaindict[chains] + str(count) + mut[-1:]))
                    fvar.write("%s\n" % (str(map_three_one[resname]) + str(count) + mut[-1:]))
        f.close()

        fpdb = open(pathoutput + "/" + pdb + "_PROTEIN_" + chains + ".pdb", "r")
        for linepdb in fpdb.readlines()[1:]:
            if linepdb[16:17] == " " or linepdb[16:17] == "A":
                resnamepdb = linepdb[17:20].strip()
                resnumpdb = linepdb[22:27].strip()
                atomnamepdb = linepdb[13:16].strip()
                if resnamepdb == resname and resnumpdb == resnum:
                    if atomnamepdb != atomname:
                        if count < 10:
                            fwpdb.write("%s %s%s   %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                        if count >= 10 and count < 100:
                            fwpdb.write("%s %s%s  %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                        if count >= 100 and count < 1000:
                            fwpdb.write("%s %s%s %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                        if count >= 1000 and count < 10000:
                            fwpdb.write("%s %s%s%s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                    else:
                        continue
                else:
                    if atomnamepdb != atomname:
                        count += 1
                        if count < 10:
                            fwpdb.write("%s %s%s   %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                        if count >= 10 and count < 100:
                            fwpdb.write("%s %s%s  %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                        if count >= 100 and count < 1000:
                            fwpdb.write("%s %s%s %s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))
                        if count >= 1000 and count < 10000:
                            fwpdb.write("%s %s%s%s %s" % (linepdb[0:16], linepdb[17:21], mapchaindict[chains], str(count), linepdb[27:]))

                        mutchainall = []
                        with open(in_file, 'r') as f:
                            _unused = f.next()
                            for line in f:
                                ff = line.split("\t")
                                mut = ff[4]
                                mutchain = ff[3]
                                if( str(chains) == str(mutchain)) and (mutchain not in mutchainall):
                                    if resnamepdb in need:
                                        fseq.write("%s" % (map_three_one[resnamepdb]))
                                    else:
                                        continue
                                mutchainall.append(mutchain)
                                if str(resnumpdb) == str(mut[1:-1]) and str(map_three_one[resnamepdb]) == str(mut[0:1]) and str(chains) == str(mutchain):
                                    fw.write("%s\t%s\t%s\t%s\t%s\n" % (
                                    line[:-1], pdb, newpartner1, newpartner2, str(map_three_one[resnamepdb]) + mapchaindict[chains] + str(count) + mut[-1:]))
                                    fvar.write("%s\n" % (str(map_three_one[resnamepdb]) + str(count) + mut[-1:]))
                    else:
                        continue
                resname = linepdb[17:20].strip()
                resnum = linepdb[22:27].strip()
                atomname = linepdb[13:16].strip()
        fpdb.close()
        fwpdb.close()
        fvar.close()
    fseq.close()
    fw.close()

    for chain_new in list(partner2):
        fwlig = open(pathoutput + '/' + pdb + '_' + ligand + '_' + chain_new + '.pdb', "r")
        fclean = open(pathoutput + '/' + pdb + '_' + ligand + '_TEST_CH' + str(countchain) + '.pdb', "w")
        count = 1
        for line in fwlig.readlines():
            if line[16:17] == " " or line[16:17] == "A":
                if chain_new not in list(partner1):
                    if count < 10:
                        fclean.write("HETATM    %s%s %s%s%s" % (str(count), line[11:16], line[17:21], mapchaindict[chain_new], line[22:]))
                    if count >= 10 and count < 100:
                        fclean.write("HETATM   %s%s %s%s%s" % (str(count), line[11:16], line[17:21], mapchaindict[chain_new], line[22:]))
                    if count >= 100 and count < 1000:
                        fclean.write("HETATM  %s%s %s%s%s" % (str(count), line[11:16], line[17:21], mapchaindict[chain_new], line[22:]))
                    if count >= 1000 and count < 10000:
                        fclean.write("HETATM %s%s %s%s%s" % (str(count), line[11:16], line[17:21], mapchaindict[chain_new], line[22:]))
                    if count >= 10000 and count < 100000:
                        fclean.write("HETATM%s%s %s%s%s" % (str(count), line[11:16], line[17:21], mapchaindict[chain_new], line[22:]))
                else:
                    if count < 10:
                        fclean.write("HETATM    %s%s %s%s%s" % (str(count), line[11:16], line[17:21], mapchaindict[same_chain_dict[chain_new]], line[22:]))
                    if count >= 10 and count < 100:
                        fclean.write("HETATM   %s%s %s%s%s" % (str(count), line[11:16], line[17:21], mapchaindict[same_chain_dict[chain_new]], line[22:]))
                    if count >= 100 and count < 1000:
                        fclean.write("HETATM  %s%s %s%s%s" % (str(count), line[11:16], line[17:21], mapchaindict[same_chain_dict[chain_new]], line[22:]))
                    if count >= 1000 and count < 10000:
                        fclean.write("HETATM %s%s %s%s%s" % (str(count), line[11:16], line[17:21], mapchaindict[same_chain_dict[chain_new]], line[22:]))
                    if count >= 10000 and count < 100000:
                        fclean.write("HETATM%s%s %s%s%s" % (str(count), line[11:16], line[17:21], mapchaindict[chain_new], line[22:]))
                count = count + 1
        fwlig.close()
        fclean.close()
        county = 1
        fclean1 = open(pathoutput + '/' + pdb + '_' + ligand + '_TEST_CH' + str(countchain) + '.pdb', "r")
        print(pathoutput + '/' + pdb + '_' + ligand + '_TEST_CH' + str(countchain) + '.pdb')
        resname = fclean1.readlines()[0][17:20].strip()
        fclean1 = open(pathoutput + '/' + pdb + '_' + ligand + '_TEST_CH' + str(countchain) + '.pdb', "r")
        resnum = fclean1.readlines()[0][22:27].strip()
        fclean1 = open(pathoutput + '/' + pdb + '_' + ligand + '_TEST_CH' + str(countchain) + '.pdb', "r")
        atomname = fclean1.readlines()[0][13:16].strip()
        line1 = open(pathoutput + '/' + pdb + '_' + ligand + '_TEST_CH' + str(countchain) + '.pdb', "r").readlines()[0]
        fclean_new = open(pathoutput + '/' + pdb + '_' + ligand + '_CH' + str(countchain) + '.pdb', "w")
        fclean_new.write("%s %s   %s %s" % (line1[0:16], line1[17:22], str(county), line1[27:]))

        fclean1 = open(pathoutput + '/' + pdb + '_' + ligand + '_TEST_CH' + str(countchain) + '.pdb', "r")
        for linepdb in fclean1.readlines()[1:]:
            if linepdb[16:17] == " " or linepdb[16:17] == "A":
                resnamepdb = linepdb[17:20].strip()
                resnumpdb = linepdb[22:27].strip()
                atomnamepdb = linepdb[13:16].strip()
                if resnamepdb == resname and resnumpdb == resnum:
                    if atomnamepdb != atomname:
                        if county < 10:
                            fclean_new.write( "%s %s   %s %s" % (linepdb[0:16], linepdb[17:22], str(county), linepdb[27:]))
                        if county >= 10 and county < 100:
                            fclean_new.write("%s %s  %s %s" % (linepdb[0:16], linepdb[17:22], str(county), linepdb[27:]))
                        if county >= 100 and county < 1000:
                            fclean_new.write("%s %s %s %s" % (linepdb[0:16], linepdb[17:22], str(county), linepdb[27:]))
                        if county >= 1000 and county < 10000:
                            fclean_new.write("%s %s%s %s" % (linepdb[0:16], linepdb[17:22], str(county), linepdb[27:]))
                    else:
                        continue
                else:
                    if atomnamepdb != atomname:
                        county += 1
                        if county < 10:
                            fclean_new.write("%s %s   %s %s" % (linepdb[0:16], linepdb[17:22], str(county), linepdb[27:]))
                        if county >= 10 and county < 100:
                            fclean_new.write("%s %s  %s %s" % (linepdb[0:16], linepdb[17:22], str(county), linepdb[27:]))
                        if county >= 100 and county < 1000:
                            fclean_new.write("%s %s %s %s" % (linepdb[0:16], linepdb[17:22], str(county), linepdb[27:]))
                        if county >= 1000 and county < 10000:
                            fclean_new.write("%s %s%s %s" % (linepdb[0:16], linepdb[17:22], str(county), linepdb[27:]))
                    else:
                        continue
                resname = linepdb[17:20].strip()
                resnum = linepdb[22:27].strip()
                atomname = linepdb[13:16].strip()
        fclean1.close()
        fclean_new.close()
        # chimera addh
        os.chdir(pathoutput)
        os.system('chimera --nogui --script "%s/runchimera.py %s %s %s"' % (pathinput, pdb, ligand, str(countchain)))
        countchain = countchain + 1


# produce renumbered and processed 1A43.pdb wild-type pdb structure
#only protein ,no ligand
# prepare for foldx,generate muta structure
def wtpdb():
    pdball = []
    with open(in_file + ".cleaned", 'r') as f:
        f.next()
        for line in f:
            ff = line.split("\t")
            pdb = ff[8]  # 1A43
            partner2 = ff[2].split(".")
            mutachain = ff[3].split("_")[0]
            if pdb not in pdball:
                os.system('cat %s/%s_PROTEIN_CH*.pdb > %s/%s.pdb' % (pathoutput, pdb, pathoutput, pdb))
                pdball.append(pdb)


# produce psf and pdb files of wild-type with vmd.
# produce 1A43_vmd.psf and 1A43_vmd.pdb using 1A43_CH1.pdb and 1A43_CH2.pdb.
def vmd_wt():
    pdball = []
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    template = open(pathinput + '/vmd.pgn').read()
    for line in f:
        ff = line.split("\t")
        pdb = ff[8]
        partner1 = ff[9]
        NumChain = int(len(partner1))
        if pdb not in pdball:
            vmd_pdb = template.replace('protname', pdb+'_PROTEIN').replace('NumChain', str(NumChain)).replace('pathinput', pathinput).replace('pathoutput', pathoutput)
            with open(pathoutput + '/vmd_protein_' + pdb + '.pgn', 'w') as fw:
                fw.write(vmd_pdb)
            os.system('%s -dispdev text -e %s/vmd_protein_%s.pgn' % (pathvmd,pathoutput, pdb))
            pdball.append(pdb)
        else:
            continue
    f.close()


def generate_wt():
    pdball = []
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    for line in f:
        ff = line.split("\t")
        pdb = ff[8].lower()
        partner1 = ff[9]
        partner2 = ff[10]
        ligand = ff[6].lower()
        mut = ff[11][:-1].lower()
        protname = pdb + '_' + mut
        countp = 1
        if pdb not in pdball:
            for chains in (list(partner1)):
                # split file through diff chain, and generate pro_ch* file
                os.system('grep "^.\{21\}%s" %s/%s_vmd.pdb > %s/%s_%s.pdb' % (
                chains, pathoutput, pdb.upper()+'_PROTEIN', pathoutput, pdb+'_protein', 'ch' + str(countp)))
                ffpdb = open(pathoutput + '/' + pdb+'_protein' + '_ch' + str(countp) + '.pdb', 'a')
                ffpdb.write("%s\n" % ('END'))
                ffpdb.close()
                countp = countp + 1
            # generate ligand_ch file
            for chain2 in partner2:
                frlig = open(pathoutput + '/' + pdb.upper() + '_' + ligand.upper() + '_ADDH_CH' + str(countp) + '.pdb', 'r')
                fwlig = open(pathoutput + '/' + pdb + '_' + ligand + '_ch' + str(countp) + '.pdb', 'w')
                lastres = 1
                for line in frlig:
                    if (line[0:6] == 'HETATM') and (countp < 10):
                        fwlig.write('ATOM  %s      CH%s   \n' % (line[6:66], str(countp)))
                    if (line[0:6] == 'ATOM  ') and (countp < 10):
                        fwlig.write('ATOM  %s      CH%s   \n' % (line[6:66], str(countp)))
                    if (line[0:6] == 'HETATM') and (countp >= 10) and (countp < 100):
                        fwlig.write('ATOM  %s     CH%s   \n' % (line[6:66], str(countp)))
                    if (line[0:6] == 'ATOM  ') and (countp >= 10) and (countp < 100):
                        fwlig.write('ATOM  %s     CH%s   \n' % (line[6:66], str(countp)))
                fwlig.write("%s\n" % ('END'))
                fwlig.close()
                frlig.close()
                countp = countp + 1
        else:
            continue
    f.close()


# make inputfile for foldx4 with mutchain for calculating folding free energy
def inputfoldx():
    os.chdir(path)
    f = open(in_file + ".cleaned", 'r')
    f.next()
    for line in f:
        ff = line.split("\t")
        mut = ff[11][:-1]  # GA9A
        pdb = ff[8]  # 1A43
        with open('individual_list_' + jobid + '_' + mut + '.txt', 'w') as fic:
            fic.write('%s' % (mut + ';'))
        with open('foldx_buildmodel_' + jobid + '_' + mut + '.txt', 'w') as fsc:
            fsc.write('command=BuildModel\npdb=%s\nmutant-file=%s' % (jobid + '_' + mut + '.pdb', 'individual_list_' + jobid + '_' + mut + '.txt'))
        os.system("cp %s/%s.pdb %s_%s.pdb" % (pathoutput, pdb, jobid, mut))
    f.close()


# build model with foldx, produce Dif_1A43_A_Repair_GA9A.fxout
def runfoldx_mut():
    f = open(in_file + ".cleaned", 'r')
    f.next()
    for line in f:
        ff = line.split("\t")
        mut = ff[11][:-1]  # GA9A
        pdb = ff[8]  # 1A43
        os.system('./foldx -f foldx_buildmodel_%s_%s.txt' % (jobid, mut))
        os.system("rm WT_%s_%s_1.pdb" % (jobid, mut))
        os.system("rm individual_list_%s_%s.txt" % (jobid, mut))
        os.system("rm foldx_buildmodel_%s_%s.txt" % (jobid, mut))
        os.system("rm Average_%s_%s.fxout" % (jobid, mut))
        os.system("rm Raw_%s_%s.fxout" % (jobid, mut))
        os.system("rm PdbList_%s_%s.fxout" % (jobid, mut))
        os.system("rm %s_%s.fxout" % (jobid, mut))
        os.system("mv Dif_%s_%s.fxout %s" % (jobid, mut, pathoutput))
        os.system("rm %s_%s.pdb" % (jobid, mut))
        os.system("mv %s_%s_1.pdb %s/%s_%s.pdb" % (jobid, mut, pathoutput, pdb, mut))
    f.close()


def splitchain_mut():
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    for line in f:
        ff = line.split("\t")
        mut = ff[11][:-1]  # GA9A
        pdb = ff[8]  # 1A43
        partner1 = ff[9]
        pdb = pdb + '_' + mut
        countchain = 1
        for chains in list(partner1):
            os.system('grep "^.\{21\}%s" %s/%s.pdb > %s/%s_%s.pdb' % (chains, pathoutput, pdb, pathoutput, pdb, 'CH' + str(countchain)))
            countchain += 1
    f.close()


# produce psf and pdb files of mutant with vmd 
def vmd_mut():
    f = open(in_file + ".cleaned", 'r')
    _unused = f.next()
    template = open(pathinput + '/vmd.pgn').read()
    for line in f:
        ff = line.split("\t")
        mut = ff[11][:-1]  # GA9A
        pdb = ff[8]  # 1A43
        partner1 = ff[9]
        protname = pdb + '_' + mut
        NumChain = int(len(partner1))
        vmd_pdb = template.replace('protname', protname).replace('NumChain', str(NumChain)).replace('pathinput', pathinput).replace('pathoutput', pathoutput)
        with open(pathoutput + '/vmd_' + protname + '.pgn', 'w') as fw:
            fw.write(vmd_pdb)
        os.system('%s -dispdev text -e %s/vmd_%s.pgn' % (pathvmd,pathoutput, protname))
    f.close()


def wtpdb_all_addh():
    pdball = []
    with open(in_file + ".cleaned", 'r') as f:
        f.next()
        for line in f:
            ff = line.split("\t")
            pdb = ff[8].lower()  # 1A43
            ligand = ff[6].lower()
            if pdb not in pdball:
                os.system('cat %s/%s_protein_ch*.pdb %s/%s_%s_ch*.pdb > %s/%s_wt_test.pdb' % (pathoutput, pdb, pathoutput, pdb, ligand, pathoutput, pdb))
                fnew = open(pathoutput + '/' + pdb +'_wt.pdb','w')
                f_old = open(pathoutput + '/' + pdb +'_wt_test.pdb','r')
                count = 1
                for line1 in f_old:
                    if line1[0:6].strip() == 'ATOM':
                        if count < 10:
                            fnew.write('ATOM      %s%s' % (str(count),line1[11:]))
                        if count >= 10 and count < 100:
                            fnew.write('ATOM     %s%s' % (str(count),line1[11:]))
                        if count >= 100 and count < 1000:
                            fnew.write('ATOM    %s%s' % (str(count),line1[11:]))
                        if count >= 1000 and count < 10000:
                            fnew.write('ATOM   %s%s' % (str(count),line1[11:]))
                        if count >= 10000 and count < 100000:
                            fnew.write('ATOM  %s%s' % (str(count),line1[11:]))
                        count = count + 1
                fnew.write("%s\n" % ('END'))
                fnew.close()
                f_old.close()

                pdball.append(pdb)
            else:
                continue
            os.system('rm %s/%s_wt_test.pdb' % (pathoutput, pdb))


def mtpdb_all_addh():
    pdball = []
    with open(in_file + ".cleaned", 'r') as f:
        f.next()
        for line in f:
            ff = line.split("\t")
            pdb = ff[8].lower()  # 1A43
            ligand = ff[6].lower()
            mut = ff[11][:-1].lower()
            pdb = pdb + '_' + mut
            if pdb not in pdball:
                os.system('cat %s/%s_vmd.pdb %s/%s_%s_ch*.pdb > %s/%s_mt_test.pdb' % (pathoutput, pdb.upper(), pathoutput, pdb.split('_')[0], ligand, pathoutput, pdb))
                fnew = open(pathoutput + '/' + pdb +'_mt.pdb','w')
                f_old = open(pathoutput + '/' + pdb +'_mt_test.pdb','r')
                count = 1
                for line1 in f_old:
                    if line1[0:6].strip() == 'ATOM':
                        if count < 10:
                            fnew.write('ATOM      %s%s' % (str(count),line1[11:]))
                        if count >= 10 and count < 100:
                            fnew.write('ATOM     %s%s' % (str(count),line1[11:]))
                        if count >= 100 and count < 1000:
                            fnew.write('ATOM    %s%s' % (str(count),line1[11:]))
                        if count >= 1000 and count < 10000:
                            fnew.write('ATOM   %s%s' % (str(count),line1[11:]))
                        if count >= 10000 and count < 100000:
                            fnew.write('ATOM  %s%s' % (str(count),line1[11:]))
                        count = count + 1
                fnew.write("%s\n" % ('END'))
                fnew.close()
                f_old.close()
                pdball.append(pdb)
            else:
                continue
            os.system('rm %s/%s_mt_test.pdb' % (pathoutput, pdb))


# Arpeggio
def arpeggio():
    pdball = []
    with open(in_file + ".cleaned", 'r') as f:
        f.next()
        for line in f:
            ff = line.split("\t")
            pdbid = ff[8].lower()  # 1A43
            ligand = ff[6].lower()
            mut = ff[11][:-1].lower()
            # wild
            if pdbid not in pdball:
                pdball.append(pdbid)
                vmd_temp = open('{}/{}_wt.pdb'.format(pathoutput, pdbid)).read().replace('HSD','HIS')
                with open('{}/{}_wt_temp.pdb'.format(pathoutput, pdbid),'w') as fw:
                    fw.write(vmd_temp)
                # run arpeggio
                os.system('python {} -he {}/{}_wt_temp.pdb'.format(patharpeggio,pathoutput, pdbid))
            # mut
            print(pdbid+'_'+mut)
            vmd_temp_mut =  open('{}/{}_mt.pdb'.format(pathoutput, pdbid+'_'+mut)).read().replace('HSD','HIS')
            with open('{}/{}_mt_temp.pdb'.format(pathoutput, pdbid+'_'+mut),'w') as fw:
                fw.write(vmd_temp_mut)
            # run arpeggio
            os.system('python {} -he {}/{}_mt_temp.pdb'.format(patharpeggio,pathoutput, pdbid+'_'+mut))
            # delete unused files 
            if os.path.exists('{}/{}_wt_temp.contacts'.format(pathoutput,pdbid)) and os.path.exists('{}/{}_mt_temp.contacts'.format(pathoutput,pdbid+'_'+mut)):
                os.system('rm {}/*_temp.pdb'.format(pathoutput))
                os.system('rm {}/*.bs_contacts'.format(pathoutput))
                os.system('rm {}/*.specific.sift'.format(pathoutput))
                os.system('rm {}/*.sift'.format(pathoutput))
                os.system('rm {}/*.specific.siftmatch'.format(pathoutput))
                os.system('rm {}/*.siftmatch'.format(pathoutput))
                os.system('rm {}/*.specific.polarmatch'.format(pathoutput))
                os.system('rm {}/*.polarmatch'.format(pathoutput))
                os.system('rm {}/*.ri'.format(pathoutput))
                os.system('rm {}/*.rings'.format(pathoutput))
                os.system('rm {}/*.ari'.format(pathoutput))
                os.system('rm {}/*.amri'.format(pathoutput))
                os.system('rm {}/*.amam'.format(pathoutput))
                os.system('rm {}/*.residue_sifts'.format(pathoutput))
            else:
                print('arpeggio run error :{}'.format(jobid))


def arpeggio_extract():
    first_line = open(in_file + ".cleaned", 'r').readlines()[0][:-1]
    fclean = pd.read_csv(in_file + ".cleaned",sep='\t',dtype=str)
    namelist = ['Atom1','Atom2','Clash','Covalent','VdWClash','VdW','Proximal','HydrogenBond','WeakHydrogenBond','HalogenBond','Ionic','MetalComplex','Aromatic','Hydrophobic','Carbonyl','Polar','WeakPolar','InteractingEntities']
    allnamelist = ['Arpeggio_'+i+'_sitep2_byAtom_wt' for i in namelist[2:-1]]+['Arpeggio_'+i+'_sitep2_byAtom_mut' for i in namelist[2:-1]] 
    with open(in_file + ".cleaned", 'r') as f, open(pathoutput+'/{}.arpeggio'.format(jobid),'w') as fw:
        fw.write('{}\t{}\t{}\n'.format('jobid',first_line,'\t'.join(allnamelist)))
        f.next()
        for line in f:
            ff = line.split("\t")
            pdbfile = ff[0]
            pdbid = ff[8].lower()  # 1A43
            mut = ff[11][:-1]
            partner1 = ff[9]
            partner2 = ff[10]
            p1_chain = [i for i in partner1]
            p2_chain = [j for j in partner2]
            mutchain = ff[3]
            col_dict = {i:namelist[i] for i in range(18)}
            # wild
            result = pd.read_csv('{}/{}_wt_temp.contacts'.format(pathoutput,pdbid), sep='\t', header=None)
            result = result.rename(columns = col_dict)
            result = result[result['Atom1']!='atom_bgn']
            for m in namelist[2:-1]:
                result[m] = result[m].astype(int)
            contact1 = result['Atom1'].str.split('/', expand=True).rename(columns = {0:'chain1',1:'loc1',2:'atom1'})
            contact2 = result['Atom2'].str.split('/', expand=True).rename(columns = {0:'chain2',1:'loc2',2:'atom2'})
            allcontacts = pd.concat([contact1, contact2,result], axis=1)
            allcontacts['chain_loc1'] = allcontacts['chain1']+allcontacts['loc1']
            allcontacts['chain_loc2'] = allcontacts['chain2']+allcontacts['loc2']
            sitep2 = allcontacts[((allcontacts['chain1'].isin([mut[1]])) & (allcontacts['loc1'].isin([mut[2:-1]])) & (allcontacts['chain2'].isin(p2_chain))) | ((allcontacts['chain2'].isin([mut[1]])) & (allcontacts['loc2'].isin([mut[2:-1]])) & (allcontacts['chain1'].isin(p2_chain)))]
            ## sitep2
            atom_sitep2 = [sum(sitep2[name]) for name in namelist[2:-1]]
            # mut
            resultmut = pd.read_csv('{}/{}_mt_temp.contacts'.format(pathoutput,pdbid+'_'+mut.lower()), sep='\t', header=None)
            resultmut = resultmut.rename(columns = col_dict)
            resultmut = resultmut[resultmut['Atom1']!='atom_bgn']
            for n in namelist[2:-1]:
                resultmut[n] = resultmut[n].astype(int)
            contact1mut = resultmut['Atom1'].str.split('/', expand=True).rename(columns = {0:'chain1',1:'loc1',2:'atom1'})
            contact2mut = resultmut['Atom2'].str.split('/', expand=True).rename(columns = {0:'chain2',1:'loc2',2:'atom2'})
            allcontactsmut = pd.concat([contact1mut, contact2mut,resultmut], axis=1)
            allcontactsmut['chain_loc1'] = allcontactsmut['chain1']+allcontactsmut['loc1']
            allcontactsmut['chain_loc2'] = allcontactsmut['chain2']+allcontactsmut['loc2']
            sitep2mut = allcontactsmut[((allcontactsmut['chain1'].isin([mut[1]])) & (allcontactsmut['loc1'].isin([mut[2:-1]])) & (allcontactsmut['chain2'].isin(p2_chain))) | ((allcontactsmut['chain2'].isin([mut[1]])) & (allcontactsmut['loc2'].isin([mut[2:-1]])) & (allcontactsmut['chain1'].isin(p2_chain)))]
            ## sitep2
            atom_sitep2mut = [sum(sitep2mut[name]) for name in namelist[2:-1]]
            fw.write('{}\t{}\t{}\n'.format(jobid, line[0:-1], '\t'.join([str(i) for i in atom_sitep2+atom_sitep2mut])))


class pdb_network:
    def __init__(self,pdb,chain=[],uid='0'):
        self.pdbid = pdb
        try:
            file_name = pathoutput+'/'+self.pdbid + ".pdb"
            print(file_name)
            f = open(file_name,'r')
            self.pdb_content  = f.read()
        except:
            print('no such pdb file')
        p = PDBParser(PERMISSIVE=1)
        s = p.get_structure(self.pdbid, file_name)[0]
        if chain != []:
            self.chain = chain
        else:
            chain_list = Selection.unfold_entities(s, 'C')  #C represents chain
            real_list = []
            for i in chain_list:
                real_list.append(i.id)
            self.chain = real_list
        self.res_map = {}
        for i in self.chain:
             self.res_map[i] = {}
        self.get_atom_info()
        self.build_side_chain_network()
        self.buid_space_tree()

    def get_atom_info(self):
        pdb_text = self.pdb_content
        chain = self.chain
        lines = pdb_text.split('\n')
        all_atom_dict = {}
        model_count = 0
        chain_ass = ''
        for l in lines:
            if l[0:5] == 'MODEL': #check model
                model_count += 1 
                if model_count ==1:
                    chain_ass = ''
                else:
                    chain_ass = '_' + str(model_count-1)
            if l[0:4] == 'ATOM':
                serial = int(l[7:11].replace(' ',''))
                #name = l[13:16].replace(' ','')
                resName = l[17:20].replace(' ','')
                chain_id = l[20:22].replace(' ','') + chain_ass
                resSeq = str(l[22:26].replace(' ',''))
                insider_code = str(l[26]).replace(' ','')
                if insider_code:
                    resSeq += insider_code
                x = float(l[31:38].replace(' ',''))
                y = float(l[39:46].replace(' ',''))
                z = float(l[47:54].replace(' ',''))
                element = l[11:17].replace(' ','')
                element_type = l[11:17].replace(' ','')
                if chain_id in self.chain:
                    if resSeq not in self.res_map[chain_id]:
                        self.res_map[chain_id][resSeq] = resName
                resSeq = str(resSeq) + '_' + resName
                if (chain_id in self.chain) and (element[0] != 'H'):
                    if chain_id in all_atom_dict.keys():
                        if resSeq in all_atom_dict[chain_id]:
                            all_atom_dict[chain_id][resSeq].append([x,y,z,element_type])
                        else:
                            all_atom_dict[chain_id][resSeq] = []
                            all_atom_dict[chain_id][resSeq].append([x,y,z,element_type])
                    else:
                        all_atom_dict[chain_id] = {}
                        all_atom_dict[chain_id][resSeq] = []
                        all_atom_dict[chain_id][resSeq].append([x,y,z,element_type])
        self.atom_dict = all_atom_dict

    def build_side_chain_network(self):
        main_chain = ['C', 'N', 'O', 'OT1', 'OT2', 'CA', 'HA', 'HN', 'HT1', 'HT2', 'HT3']
        self.residue_network ={}
        self.have_no_side_chain = []
        for chain_id in self.atom_dict:
            self.residue_network[chain_id] = {}
            for res_id  in self.atom_dict[chain_id]:
                res_atom_crood_list = self.atom_dict[chain_id][res_id]
                x_t = 0
                y_t = 0
                z_t = 0
                x_main = 0
                y_main = 0
                z_main = 0
                atom_count = 0
                main_count = 0
                for atom_crood in res_atom_crood_list:
                    if atom_crood[3] not in main_chain:
                        x_t += atom_crood[0]
                        y_t += atom_crood[1]
                        z_t += atom_crood[2]
                        atom_count += 1
                    else:
                        x_main += atom_crood[0]
                        y_main += atom_crood[1]
                        z_main += atom_crood[2]
                        main_count += 1
                if atom_count != 0:
                    self.residue_network[chain_id][res_id] = [round(x_t/atom_count,4),round(y_t/atom_count,4),round(z_t/atom_count,4)]
                elif main_chain != 0:
                    self.residue_network[chain_id][res_id] = [round(x_main/main_count,4),round(y_main/main_count,4),round(z_main/main_count,4)]
                    self.have_no_side_chain.append(res_id)

    def buid_space_tree(self):
        self.side_map = {}
        side_X = []
        side_uid = 0
        for chain in self.residue_network:
            for res in self.residue_network[chain]:
                side_X.append(self.residue_network[chain][res])
                self.side_map[side_uid] = chain + ':' + res
                side_uid += 1
        self.atom_map = {}
        atom_X = []
        atom_uid = 0
        for chain in self.atom_dict:
            for res in self.atom_dict[chain]:
                for atom in self.atom_dict[chain][res]:
                    if atom[3][0] != 'H':
                        atom_X.append([atom[0],atom[1],atom[2]])
                        self.atom_map[atom_uid] = chain + ':' + res
                        atom_uid += 1
        self.side_tree = KDTree(side_X)
        self.atom_tree = KDTree(atom_X)
        #print(self.side_map)

    def get_nearby(self,chain,pos,distance,n_type='atom'):
        pos = str(pos)
        result_list = []
        if n_type == 'atom':
            pos_complex = pos + '_' + self.res_map[chain][pos]
            atom_list = self.atom_dict[chain][pos_complex]
            for atom_crood in atom_list:
                if atom_crood[3][0] != 'H':
                    indices = self.atom_tree.query_radius([[atom_crood[0],atom_crood[1],atom_crood[2]]], r=distance)
                    for pos in indices[0]:
                        near_res = self.atom_map[int(pos)]
                        if near_res not in result_list:
                            result_list.append(near_res)
        return result_list


def mutsite_interface_residue(): # Interface/Non-interface
    first_line = open(in_file + ".cleaned", 'r').readlines()[0][:-1]
    fnum = open(pathoutput+'/'+jobid+'_mutsite_in_dis5A.txt','w')
    fnum.write('%s\t%s\t%s\n'%('jobid',first_line,'Interface'))
    #wt
    fw =open(in_file + ".cleaned", 'r')
    for needline in fw.readlines()[1:]:
        pdb_low = needline.split('\t')[8].lower()
        ligand = needline.split('\t')[6]
        lig_chain = needline.split('\t')[10]
        mut = needline.split('\t')[11].strip('\n').lower()
        lig_pos_list = []
        fx = open(pathoutput+'/'+pdb_low+ '_wt.pdb','r')
        for line in fx:
            if (line[0:3]!='END') and (line[17:20].strip()==ligand) and (line[20:22].strip() in lig_chain):
                lig_pos = line[20:22].strip()+';'+line[22:26].strip()
                if lig_pos not in lig_pos_list:
                    lig_pos_list.append(lig_pos)
        #print lig_pos_list
        fx.close()

        all_atom_info = pdb_network(pdb_low+'_wt')
        result_list_all = []
        for x in lig_pos_list:
            ligand_chain = x.split(';')[0]
            ligand_pos = x.split(';')[1]
            result_list = all_atom_info.get_nearby(ligand_chain,ligand_pos,5,n_type='atom')
            #print result_list
            result_list_all.append(result_list)
        result_list_final = list(set([j for i in result_list_all for j in i]))
        judge = [1 for i in result_list_final if (mut.upper()[1] == i.split(':')[0]) and (mut.upper()[2:-1] == i.split(':')[1].split('_')[0])]
        if judge == []:
            judge_new = 0
        elif judge == [1]:
            judge_new = 1
        fnum.write('%s\t%s\t%s\n'%(jobid,needline[0:-1],judge_new))
    fw.close()
    fnum.close()


def lig_property():
    ff = pd.read_csv(in_file + ".cleaned",sep='\t',dtype=str)
    lig = list(set(ff['LIG_ID']))[0].upper()
    os.system('cp %s/%s_ideal.sdf %s/' % (jobpath, lig, pathoutput))
    os.system('%sxlogp3.lnx.x86_64 -p %s/%s_ideal.sdf %s/%s_lig_xlogp3.out %sDEFAULT.TTDB' % (pathxlogp3,pathoutput, lig, pathoutput, lig, pathxlogp3))
    fx = pd.read_csv(pathoutput+'/'+lig+'_lig_xlogp3.out',names=['tag'])
    fy = open(pathoutput+'/'+jobid+'_xlogp3_result.txt','w')
    fy.write('LIG_ID\tMWt\tN_NO\n')
    mweight = ';'.join([line.split(':')[1].strip() for line in fx['tag'] if 'Molecular Weight' in line])
    noatom = ';'.join([line.split(':')[1].strip() for line in fx['tag'] if 'Nitrogen and Oxygen Atoms' in line])
    fy.write('%s\t%s\t%s\n'%(lig.upper(),mweight,noatom))


# calculate secondary structure  with DSSP using wild type crystal structure of mutchain, 1A43_A_1.pdb, produce 1A43.dssp
def dssp():
    pdball = []
    with open(in_file + ".cleaned", 'r') as f:
        f.next()
        for line in f:
            ff = line.strip().split('\t')
            pdb = ff[8]  # 1A43
            mut = ff[11]  # GA9A
            # ./mkdssp -i 1crn.pdb -o 1crn.dssp
            if pdb not in pdball:
                os.system('%s -i %s/%s.pdb -o %s/%s.dssp' % (pathmkdssp,pathoutput, pdb, pathoutput, pdb))


# hydro infomation 
def hydro():
    rotamers = pd.read_csv(pathinput + '/Rotamers_hydro_Margo', header=0, index_col=0, dtype=str, sep='\t')
    with open('{}/{}_hydro.txt'.format(pathoutput, jobid), 'w') as fw, open(in_file + '.cleaned') as f:
        f.next()
        fw.write('{}\n'.format('\t'.join(['PDB_ID', 'Mutation_PDB', 'Hydro'])))
        for line in f:
            ff = line.strip().split('\t')
            pdb = ff[8]  # 1A43
            Mutation_PDB = ff[4]  # G9A
            wild_residue = Mutation_PDB[0]  # R
            Hwt = rotamers.loc['H', wild_residue]
            fw.write('{}\n'.format('\t'.join([pdb, Mutation_PDB, Hwt])))


def run_provean():
    mutchainall = []
    with open(in_file + ".cleaned", 'r') as f1:
        f1.next()
        for line in f1:
            ff = line.strip().split('\t')
            mutchain = ff[3]
            pdb = ff[8]
            mutation_cleaned = ff[-1]
            if mutchain not in mutchainall:
                mutchainall.append(mutchain)
                os.system('%s -q %s/%s_%s.seq -v %s/%s_%s.var > %s/provean_%s_%s.out  --num_threads 30' % (pathprovean,pathoutput,pdb,mutchain,pathoutput,pdb,mutchain,pathoutput,pdb,mutchain))


# P_Q
def solart():
    pdball = []
    allclass = ['buried_Q']
    with open(in_file+'.cleaned') as f, open(pathoutput+'/'+jobid+'.solart','w') as fw:
        f.next()
        fw.write('%s\t%s\n' %('\t'.join(['PDBid','Mutation_cleaned']),'\t'.join([i+'_wt' for i in allclass])))
        for line in f:
            linelist = line.strip().split('\t')
            pdb = linelist[8]
            chains = linelist[9]
            mut = linelist[11]
            # Denominator: len_protein
            count_pro = set()
            with open(pathoutput + '/' + jobid + '_p.pdb', 'r') as f:
                for row in f:
                    count_pro.add((row[72:-1], row[22:27].strip()))
                pro_len = len(count_pro)
            # wt
            if pdb not in pdball:
                pdball.append(pdb)
                dssc = defaultdict(list)
                daac_wt = defaultdict(list)
                for chain in chains:
                    with open(pathoutput+'/'+pdb+'.dssp') as fdssp:
                        for ldssp in fdssp:
                            if re.match(r'^\s+\d+\s+\d+ {}'.format(chain), ldssp):
                                res = ldssp[13]
                                resnum = ldssp[5:10].strip()
                                reschain = ldssp[11]
                                ss = ldssp[16]
                                if ss not in ['H','E']:
                                    ss = 'C'
                                resacc = ldssp[35:39].strip()
                                resacc_tri = float(resacc)/map_surface[res]
                                if resacc_tri <= 0.2:
                                    resloc = 'buried'
                                elif resacc_tri < 0.5:
                                    resloc = 'mod_buried'
                                else:
                                    resloc = 'exposed'
                                dssc[resloc+'_'+ss].append(reschain+'_'+resnum)
                                daac_wt[resloc+'_'+res].append(reschain+'_'+resnum)
            # Amino acid composition
            aacresult = []
            # wt
            # 20中aa组分
            aacdic_wt = defaultdict(float)
            for aac in allclass:
                if aac in daac_wt.keys():
                    aacdic_wt[aac] = float(len(daac_wt[aac]))/pro_len
                else:
                    aacdic_wt[aac] = 0
                aacresult.append(aacdic_wt[aac])
            # write
            aacresult = [str(i) for i in aacresult]
            fw.write(pdb+'\t'+mut+'\t'+'\t'.join(aacresult)+'\n')


# PSSM 
def run_PSSM():
    with open(in_file+'.cleaned') as f:
        f.next()
        pdball = []
        for line in f:
            linelist = line.strip().split('\t')
            pdb = linelist[8]
            mutchain = linelist[3]
            inputseq_wt = pathoutput+'/'+pdb+'_'+mutchain+'.seq' # !!!
            output_wt = pathoutput+'/'+pdb+'_'+mutchain
            # run PSSM
            if pdb+mutchain not in pdball:
                pdball.append(pdb+mutchain)
                os.system('{} -query {} -db {} -num_iterations 3 -out_ascii_pssm {}.pssm'.format(pathpsiblast,inputseq_wt,pathblastdb,output_wt))


# get result form PSSM
def get_PSSM():
    dindex = {'A':0,'R':1,'N':2,'D':3,'C':4,'Q':5,'E':6,'G':7,'H':8,'I':9,'L':10,'K':11,'M':12,'F':13,'P':14,'S':15,'T':16,'W':17,'Y':18,'V':19}
    with open(in_file+'.cleaned') as f,open(pathoutput+'/'+jobid+'.pssm','w') as fw:
        title = f.next()
        fw.write(title.strip()+'\tPSSM\n')
        for line in f:
            linelist = line.strip().split('\t')
            pdb = linelist[8]
            mutchain = linelist[3]
            mut = linelist[11]
            # wild
            with open(pathoutput+'/'+pdb+'_'+mutchain+'.pssm') as fpssm:
                for ffp in fpssm:
                    if re.match(r'\s+'+mut[2:-1]+'\s+'+mut[0], ffp):  # !!!
                        position_score_wtwt = ffp.strip().split()[dindex[mut[0]]+2]
            # write
            fw.write(line.strip()+'\t'+position_score_wtwt+'\n')


# P_RKDE  
def AAcompo():
    residues = ['K', 'D', 'R', 'E']
    namelist_wt = ['AAcompo_RKDE_wt']
    namelist = namelist_wt
    with open(in_file+'.cleaned') as f, open(pathoutput+'/amino_acid_composition_mutchain.txt','w') as fw:
        fw.write('PDBid\tMutation_cleaned\t'+'\t'.join(namelist)+'\n')
        f.next()
        pdball = []
        for line in f:
            linelist = line.strip().split('\t')
            pdb = linelist[8]
            mutchain = linelist[3]
            mut = linelist[11]
            location = int(mut[2:-1])
            # wild
            if pdb+mutchain not in pdball:
                pdball.append(pdb+mutchain)
                seqwt = list(open('{}/{}_{}.seq'.format(pathoutput,pdb,mutchain)).readlines()[1])
                length = len(seqwt)
                composition_wt = dict(Counter(seqwt))
                dwt = {'AAcompo_'+i+'_wt':((composition_wt[i]+0.00)/length) if i in composition_wt.keys() else 0.00 for i in residues}
                dwt.update({'AAcompo_RKDE_wt':dwt['AAcompo_K_wt']+dwt['AAcompo_R_wt']+dwt['AAcompo_D_wt']+dwt['AAcompo_E_wt']})
                print(dwt)
            fw.write('{}\t{}\t{}\n'.format(pdb, mut,'\t'.join([str(dwt[i]) for i in namelist_wt])))


# get all features to 1a43.input.cleanedpro.outdata
def getenergy():
    # make the header for aaindex
    aaindex_2 = pd.read_csv(pathinput + '/aaindex2.txt', sep='\t', index_col=None)
    aaindex_list = [i for i in aaindex_2.columns][1:]
    aaindex_header = '\t'.join(['M_AA1','M_AA2'])
    fw = open(in_file + ".cleaned.outdata", 'w')
    first_line = open(in_file + ".cleaned", 'r').readlines()[0][:-1]
    fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (first_line,'DCS','Hydro','PSSM',aaindex_header,'P_RKDE','P_Q'))
    f = open(in_file + ".cleaned", 'r')
    f.next()
    for line in f:
        ff = line.strip().split("\t")
        pdb = ff[8].lower()  # 1a43
        partner1 = ff[9]
        Mutation_PDB = ff[4].strip()  # G156A
        mutchain_cleaned = ff[3]  # A_1
        mut = ff[11].lower()  # ga9a
        mutchain = mut[1:2]  # a
        resnum = mut[2:-1]  # 9
        protname = pdb + '_' + mut  # 1a43_ga9a
        aaindex = ''
        pdball = []
        count = 1
        mapchains = defaultdict(str)
        for i in list(partner1):
            mapchains[i] = 'CH'+str(count)
            count+=1

        # add aaindex
        sub_type = Mutation_PDB[0] + Mutation_PDB[-1]
        with open(pathinput + '/aaindex2.txt') as faai:
            faai.next()
            for laai in faai:
                laai_list = laai.strip().split('\t')
                if laai_list[0] == sub_type or laai_list[0] == sub_type[::-1]:
                    laai_list = laai_list[1:]
                    if '-' in laai_list:
                        laai_list[laai_list.index('-')] = 'NONE'
                        aaindex = '\t'.join(laai_list)
                    else:
                        aaindex = '\t'.join(laai_list)

        # add provean score
        ST_play = False
        with open(pathoutput+'/provean_'+pdb.upper()+'_'+mutchain_cleaned+'.out','r') as fprovean:
            for ffprovean in fprovean:
                if ffprovean[2:11] == "VARIATION":
                    ST_play = True
                    continue
                if ST_play:
                    ffp = ffprovean.strip().split("\t")
                    if ffp[0] == (mut[0]+mut[2:]).upper():
                        provean_score = ffp[1]

        # hydro
        with open(pathoutput + '/{}_hydro.txt'.format(jobid)) as frm:
            frm.next()
            for lrm in frm:
                lrm_list = lrm.strip().split('\t')
                if Mutation_PDB == lrm_list[1]:
                    hydro_wt = '\t'.join(lrm_list[2:])

        # PSSM
        with open(pathoutput + '/{}.pssm'.format(jobid)) as fpssm:
            fpssm.next()
            for lpssm in fpssm:
                lpssm_list = lpssm.strip().split('\t')
                if Mutation_PDB == lpssm_list[4]:
                    pssm_out = '\t'.join(lpssm_list[12:])

        # solart    
        with open(pathoutput + '/{}.solart'.format(jobid)) as fsol:
            fsol.next()
            for lsol in fsol:
                lsol_list = lsol.strip().split('\t')
                if mut.upper() == lsol_list[1]:
                    sol_ratio = '\t'.join(lsol_list[2:])

        aacomp = pd.read_csv(pathoutput+'/amino_acid_composition_mutchain.txt', header=0, index_col=None, sep='\t',dtype=str)
        aacomp_cols = ['AAcompo_RKDE_wt']
        aacomp_new = aacomp[aacomp['Mutation_cleaned']==mut.upper()].loc[:,aacomp_cols]
        aacomp_features = '\t'.join(str(aacomp_new[i].tolist()[0].strip()) for i in aacomp_cols)
        # write final file
        fw.write(
            "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (line.strip(),provean_score,hydro_wt,pssm_out,aaindex,aacomp_features,sol_ratio))
    f.close()
    fw.close()


def get_features():
    fclean = pd.read_csv(in_file + ".cleaned",sep='\t',dtype=str)
    f_arpeggio = pd.read_csv(pathoutput + '/' + jobid + '.arpeggio',sep='\t',dtype=str)
    f_mutsite_interface = pd.read_csv(pathoutput+'/'+jobid+'_mutsite_in_dis5A.txt',sep='\t',dtype=str)
    f_output1 = pd.merge(f_arpeggio,f_mutsite_interface,on=['jobid', 'PDBfile', 'Partner1', 'Partner2', 'LIG_ID', 'Mutation_cleaned','NewPartner1','NewPartner2','MutChain','Mutation_PDB','Result_Id','LIG_POS','PDBid'],how='left')
    lig = list(set(f_output1['LIG_ID']))[0].upper()
    f_ligand = pd.read_csv(pathoutput+'/'+jobid+'_xlogp3_result.txt',sep='\t',dtype=str)
    feature_new = pd.merge(f_output1,f_ligand,on=['LIG_ID'],how='left')
    f_output2 = pd.read_csv(in_file + ".cleaned.outdata",sep='\t',dtype=str)
    f_output3 = pd.merge(feature_new,f_output2,on = ['NewPartner1', 'PDBfile', 'Partner2', 'PDBid', 'NewPartner2', 'Result_Id', 'LIG_POS', 'Mutation_cleaned', 'Mutation_PDB', 'LIG_ID', 'MutChain', 'Partner1'],how='left')
    f_output3.rename(columns={'Arpeggio_Proximal_sitep2_byAtom_wt':'Prox_wt','Arpeggio_Proximal_sitep2_byAtom_mut':'Prox_mut'},inplace=True)
    fmodel = f_output3.loc[:,['PDBfile','Partner1','Partner2','MutChain','Mutation_PDB','Result_Id','LIG_ID','LIG_POS','PDBid','NewPartner1','NewPartner2','Mutation_cleaned','Interface','Prox_wt','Prox_mut','DCS','PSSM','P_RKDE','MWt','N_NO','Hydro','P_Q','M_AA1','M_AA2']]
    fmodel.to_csv(in_file + ".cleaned.output.model", sep='\t', index=False)


def Prediction():
    outdata = in_file + '.cleaned.output.model'
    robjects.globalenv["outdata"] = outdata
    robjects.globalenv["workdir"] = path
    r('''test = read.table(outdata,header=T,sep="\t")''')
    r('''filename = paste(workdir, 'inputfiles/PremPLI.RData',sep = '')''')
    r('''load(file = filename)''')
    PredR = r('''predict(prempli.rf,test)''')
    cutoff = 1
    robjects.globalenv["PredR"] = PredR
    predDDG = r('''PredR''')
    first_line = open(in_file, 'r').readlines()[0][:-1]
    fw = open(out_file, "w")
    fw.write("%s\t%s\t%s\n" % (first_line, "PremPLI", "Interface"))
    f = open(outdata, 'r')
    _unused = f.next()
    count = 0
    for line in f:
        ff = line.split("\t")
        if str(ff[12]).strip() == '1':
            a = 'Yes'
        elif str(ff[12]).strip() == '0':
            a = 'No'
        fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%3.2f\t%s\n" % (ff[0], ff[1], ff[2], ff[3], ff[4], ff[5], ff[6], ff[7], predDDG[count], a))
        count += 1
    f.close()
    fw.close()
    os.system("rm %s/%s.input.cleaned.*" % (jobpath,jobid))


def main():
    ProPDB1()
    del_unknown_incomplete()
    splitchain()
    splitligand()
    CleanPdb()
    wtpdb()
    vmd_wt()
    generate_wt()
    inputfoldx()
    runfoldx_mut()
    splitchain_mut()
    vmd_mut()
    wtpdb_all_addh()
    mtpdb_all_addh()
    arpeggio()
    arpeggio_extract()
    mutsite_interface_residue()
    lig_property()
    dssp()
    hydro()
    run_provean()
    solart()
    run_PSSM()
    get_PSSM()
    AAcompo()
    getenergy()
    get_features()
    Prediction()

if __name__ == '__main__':
    main()
