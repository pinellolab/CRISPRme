#! /usr/bin/env python
#  -*- coding: utf-8 -*-
import json
import sys
try:
    import Tkinter as tk
except ImportError:
    import tkinter as tk

try:
    import ttk
    py3 = False
except ImportError:
    import tkinter.ttk as ttk
    py3 = True

from tkinter import filedialog, END, IntVar, messagebox
import gzip
import os
from shutil import copy, rmtree, move
# from change_dict import updateDictionary
# import change_dict


class WaitingWindow(tk.Toplevel):
    def __init__(self, top=None):

        tk.Toplevel.__init__(self, top)
        self.root = top

        self.geometry("600x344+650+150")
        self.minsize(1, 1)
        self.maxsize(1825, 970)
        self.resizable(1, 1)
        self.title("Progress")

        self.FrameDict = tk.Frame(self)
        self.FrameDict.place(relx=0.015, rely=0.029,
                             relheight=0.131, relwidth=0.975)
        self.FrameDict.configure(relief='groove')
        self.FrameDict.configure(borderwidth="2")
        self.FrameDict.configure(relief="groove")

        self.LabelDict = tk.Label(self.FrameDict)
        self.LabelDict.place(relx=0.017, rely=0.222, height=25, width=155)
        self.LabelDict.configure(text='''Updating dictionaries''')

        self.LabelDictStatus = tk.Label(self.FrameDict)
        self.LabelDictStatus.place(
            relx=0.632, rely=0.222, height=25, width=118)

        self.LabelDictStatus.configure(text='''Pending...''')

        self.FrameVariants = tk.Frame(self)
        self.FrameVariants.place(
            relx=0.017, rely=0.174, relheight=0.131, relwidth=0.975)
        self.FrameVariants.configure(relief='groove')
        self.FrameVariants.configure(borderwidth="2")
        self.FrameVariants.configure(relief="groove")

        self.LabelVariants = tk.Label(self.FrameVariants)
        self.LabelVariants.place(relx=0.017, rely=0.222, height=25, width=109)
        self.LabelVariants.configure(text='''Adding Variants''')

        self.LabelVariantsStatus = tk.Label(self.FrameVariants)
        self.LabelVariantsStatus.place(
            relx=0.632, rely=0.222, height=25, width=118)
        self.LabelVariantsStatus.configure(text='''Pending...''')

        self.FrameIdx = tk.Frame(self)
        self.FrameIdx.place(relx=0.017, rely=0.32,
                            relheight=0.131, relwidth=0.975)

        self.FrameIdx.configure(relief='groove')
        self.FrameIdx.configure(borderwidth="2")
        self.FrameIdx.configure(relief="groove")

        self.LabelIdx = tk.Label(self.FrameIdx)
        self.LabelIdx.place(relx=0.015, rely=0.222, height=25, width=199)
        self.LabelIdx.configure(text='''Re-indexing enriched genome''')

        self.LabelIdxStatus = tk.Label(self.FrameIdx)
        self.LabelIdxStatus.place(relx=0.632, rely=0.222, height=25, width=118)
        self.LabelIdxStatus.configure(text='''Pending...''')

        self.FrameSample = tk.Frame(self)
        self.FrameSample.place(relx=0.017, rely=0.465,
                               relheight=0.131, relwidth=0.975)
        self.FrameSample.configure(relief='groove')
        self.FrameSample.configure(borderwidth="2")
        self.FrameSample.configure(relief="groove")

        self.LabelSample = tk.Label(self.FrameSample)
        self.LabelSample.place(relx=0.017, rely=0.222, height=25, width=99)
        self.LabelSample.configure(text='''Adding samples''')

        self.LabelSampleStatus = tk.Label(self.FrameSample)
        self.LabelSampleStatus.place(
            relx=0.632, rely=0.222, height=25, width=118)
        self.LabelSampleStatus.configure(text='''Pending...''')

        self.ButtonDismiss = tk.Button(self)
        #self.ButtonDismiss.place(relx=0.433, rely=0.814, height=35, width=80)
        self.ButtonDismiss.configure(text='''Dismiss''')

    def doneLabel(self, label):
        if label == "dict":
            label = self.LabelDictStatus
        elif label == "var":
            label = self.LabelVariantsStatus
        elif label == "idx":
            label = self.LabelIdxStatus
        elif label == "sample":
            label = self.LabelSampleStatus
        label.configure(text="DONE!")
        self.update()

    def placeDismiss(self):
        self.ButtonDismiss.place(relx=0.433, rely=0.814, height=35, width=80)
        self.ButtonDismiss.configure(command=self.closeAll)
        self.update()

    def closeAll(self):
        self.destroy()
        self.quit()
        self.root.deiconify()
        self.root.destroy()
        self.root.quit()


def updateDictionary(oldDictionaryFile, newVCFFile):
    '''
    oldDictionaryFile : dictionary .json that will be updated
    newVCFFile: .vcf.gz file containing the new entries to be added to the old .json
    '''
    isSorted = True
    verbose = False
    oldDict = None
    oldDict = json.load(open(oldDictionaryFile))

    chr_dict = dict()
    start_time = time.time()
    oldEntry = 0
    newEntry = 0
    with gzip.open(newVCFFile, 'rb') as targets:
        if verbose:
            logFile = open("logUpdateDictionaries/log__"+os.path.basename(
                oldDictionaryFile)+"__"+os.path.basename(newVCFFile)+".log", 'w')
        # Skip vcf header
        for line in targets:
            line = line.decode('ascii')
            if ('#CHROM') in line:
                # Save this header for retrieving sample id
                column_vcf = line.strip().split('\t')
                break
        # Save CHROM [0], POS[1], REF [3], ALT [4], List of Samples [9:]
        for line in targets:
            line = line.decode('ascii').strip().split('\t')
            list_samples = []
            list_chars = []
            # if sample has 1|1 0|1 or 1|0, #NOTE may change for different vcf
            for pos, i in enumerate(line[9:]):
                if ('1' in i):
                    list_samples.append(column_vcf[pos + 9])
            if "chr" not in line[0]:
                line[0] = "chr"+str(line[0])
            chr_pos_string = line[0] + ',' + line[1]
            # Add in last two position the ref and alt nucleotide, eg: chrX,100 -> sample1,sample5,sample10;A,T
            # If no sample was found, the dict is chrX,100 -> ;A,T
            list_chars.append(line[3])
            list_chars.append(line[4])

            if chr_pos_string in oldDict:  # entry already present
                oldEntry += 1
                variations = findVars(oldDict[chr_pos_string])
                varToFind = list_chars[0]+","+list_chars[1]
                if varToFind in variations:  # variant already present
                    entry, start, end = getSamples(
                        oldDict[chr_pos_string], variations, varToFind)
                    entry = entry.split(",")
                    goodNewSamples = 0
                    for sample in list_samples:
                        if isSorted:
                            entry, inserted = putIntoEntrySorted(entry, sample)
                        else:
                            entry, inserted = putIntoEntry(entry, sample)
                        if inserted:
                            goodNewSamples += 1
                    if goodNewSamples > 0:
                        if verbose:
                            print("New Sample(s) in "+chr_pos_string +
                                  ": #"+str(goodNewSamples)+" new", file=logFile)
                        oldDict[chr_pos_string] = oldDict[chr_pos_string][:start] + \
                            ','.join(entry) + oldDict[chr_pos_string][end:]
                else:  # variant not present
                    if verbose:
                        print("New Variation in "+chr_pos_string+": "+varToFind +
                              " with #"+str(len(list_samples))+" samples", file=logFile)
                    oldDict[chr_pos_string] = oldDict[chr_pos_string] + \
                        "/" + ','.join(list_samples) + ';' + \
                        ','.join(list_chars)
            else:  # entry not present
                newEntry += 1
                if verbose:
                    print("New Entry "+chr_pos_string+" with " +
                          str(len(list_samples))+" samples", file=logFile)
                try:
                    oldDict[chr_pos_string] = ','.join(
                        list_samples) + ';' + ','.join(list_chars)
                except:
                    oldDict[chr_pos_string] = ';' + ','.join(list_chars)
        if verbose:
            logFile.close()

    with open(os.path.dirname(oldDictionaryFile) + "/my_dict_" + str(line[0]) + ".json", 'w') as f:
        json.dump(oldDict, f)
    print('Updated ' + oldDictionaryFile + ' in', time.time() - start_time)
    print("Old entries updated: "+str(oldEntry),
          "New entries inserted: "+str(newEntry))


def startUpdateDict(pathDir):
    root = tk.Tk()
    root.withdraw()
    top = TopUpdateDict(root, pathDir)
    #top = WaitingWindow(root)
    root.mainloop()
    root.quit()


def getAssociatedPAM(genomeName, pathDir):
    pams, bMaxs = [], []
    for d in os.listdir(pathDir+"/genome_library/"):
        if genomeName in d:
            parts = d.split("_")
            pams.append(parts[0])
            bMaxs.append(parts[1])
    return pams, bMaxs


class TopUpdateDict(tk.Toplevel):
    def __init__(self, top=None, pathDir=None):

        tk.Toplevel.__init__(self, top)
        self.protocol("WM_DELETE_WINDOW", self.closeAll)
        self.pathDir = pathDir
        self.root = top
        self.geometry("600x250+640+212")
        self.minsize(1, 1)
        #top.maxsize(1825, 970)
        self.resizable(1, 1)
        self.title("Change Dictionary")
        self.configure(highlightcolor="black")

        self.FrameStatus = tk.Frame(self)
        self.FrameStatus.place(relx=0.017, rely=0.003,
                               relheight=0.099, relwidth=0.975)

        self.LabelStatus = tk.Label(self.FrameStatus)
        self.LabelStatus.place(relx=0.014, rely=0.179, height=22, width=569)
        self.LabelStatus.configure(activebackground="#f9f9f9")
        self.LabelStatus.configure(
            text='''Please choose the following directories and file''')

        self.FrameDict = tk.Frame(self)
        self.FrameDict.place(relx=0.033, rely=0.086,
                             relheight=0.28, relwidth=0.942)

        """
        self.LabelDict = tk.Label(self.FrameDict)
        self.LabelDict.place(relx=0.018, rely=0.13, height=19, width=176)
        self.LabelDict.configure(activebackground="#f9f9f9")
        self.LabelDict.configure(cursor="fleur")
        self.LabelDict.configure(text='Select dictionary folder')
        """

        self.ButtonDict = tk.Button(self.FrameDict)
        self.ButtonDict.place(relx=0.018, rely=0.307, height=35, width=130)
        self.ButtonDict.configure(activebackground="#f9f9f9")
        self.ButtonDict.configure(text='''Choose Dictionary''')
        self.ButtonDict.configure(command=self.pickDictFile)

        self.TextDict = tk.Text(self.FrameDict)
        self.TextDict.place(relx=0.257, rely=0.307,
                            relheight=0.362, relwidth=0.772)
        self.TextDict.configure(background="white")
        self.TextDict.configure(font="TkTextFont")
        self.TextDict.configure(selectbackground="#c4c4c4")
        self.TextDict.configure(wrap="word")
        self.TextDict.configure(state="disabled")

        self.FrameVCF = tk.Frame(self)
        self.FrameVCF.place(relx=0.033, rely=0.326,
                            relheight=0.28, relwidth=0.942)

        """
        self.LabelVCF = tk.Label(self.FrameVCF)
        self.LabelVCF.place(relx=0.007, rely=0.119, height=17, width=139)
        self.LabelVCF.configure(activebackground="#f9f9f9")
        self.LabelVCF.configure(text='Select VCF folder')
        """

        self.ButtonVCF = tk.Button(self.FrameVCF)
        self.ButtonVCF.place(relx=0.018, rely=0.307, height=35, width=130)
        self.ButtonVCF.configure(activebackground="#f9f9f9")
        self.ButtonVCF.configure(text='Choose VCFs')
        self.ButtonVCF.configure(command=self.pickVCFDir)

        self.TextVCF = tk.Text(self.FrameVCF)
        self.TextVCF.place(relx=0.257, rely=0.307,
                           relheight=0.424, relwidth=0.772)
        self.TextVCF.configure(background="white")
        self.TextVCF.configure(font="TkTextFont")
        self.TextVCF.configure(selectbackground="#c4c4c4")
        self.TextVCF.configure(wrap="word")
        self.TextVCF.configure(state="disabled")

        self.FrameSample = tk.Frame(self)
        self.FrameSample.place(relx=0.033, rely=0.566,
                               relheight=0.28, relwidth=0.942)

        """
        self.LabelSample = tk.Label(self.FrameSample)
        self.LabelSample.place(relx=0.007, rely=0.119, height=17, width=171)
        self.LabelSample.configure(activebackground="#f9f9f9")
        self.LabelSample.configure(text='Select sampleID to add')
        """

        self.ButtonSample = tk.Button(self.FrameSample)
        self.ButtonSample.place(relx=0.018, rely=0.307, height=35, width=130)
        self.ButtonSample.configure(activebackground="#f9f9f9")
        self.ButtonSample.configure(text='Choose SamplesID')
        self.ButtonSample.configure(command=self.pickSampleFile)

        self.TextSample = tk.Text(self.FrameSample)
        self.TextSample.place(relx=0.257, rely=0.307,
                              relheight=0.424, relwidth=0.772)
        self.TextSample.configure(background="white")
        self.TextSample.configure(font="TkTextFont")
        self.TextSample.configure(selectbackground="#c4c4c4")
        self.TextSample.configure(wrap="word")
        self.TextSample.configure(state="disabled")

        """
        self.FrameChoice = tk.Frame(top)
        self.FrameChoice.place(relx=0.033, rely=0.559, relheight=0.307
                , relwidth=0.942)
        self.FrameChoice.configure(relief='groove')
        self.FrameChoice.configure(borderwidth="2")
        self.FrameChoice.configure(relief="groove")

        self.LabelOverwrite = tk.Label(self.FrameChoice)
        self.LabelOverwrite.place(relx=0.018, rely=0.042, height=21, width=229)
        self.LabelOverwrite.configure(text='Overwrite existing dictionary?')

        self.CheckbuttonOverwrite = tk.Checkbutton(self.FrameChoice)
        self.CheckbuttonOverwrite.place(relx=0.416, rely=0.083, relheight=0.135
                , relwidth=0.028)
        self.CheckbuttonOverwrite.configure(activebackground="#f9f9f9")
        self.CheckbuttonOverwrite.configure(justify='left')
        self.varCheck = IntVar()
        self.CheckbuttonOverwrite.configure(var=self.varCheck, command=self.disableFrameName)

        self.FrameName = tk.Frame(self.FrameChoice)
        self.FrameName.place(relx=0.016, rely=0.302, relheight=0.677
                , relwidth=0.965)

        self.LabelName = tk.Label(self.FrameName)
        self.LabelName.place(relx=0.006, rely=0.323, height=18, width=179)
        self.LabelName.configure(activebackground="#f9f9f9")
        self.LabelName.configure(text='Updated dictionary name:')

        self.EntryName = tk.Entry(self.FrameName)
        self.EntryName.place(relx=0.33, rely=0.323,height=17, relwidth=0.653)
        self.EntryName.configure(background="white")
        self.EntryName.configure(font="TkFixedFont")
        self.EntryName.configure(selectbackground="#c4c4c4")
        """
        self.ButtonNext = tk.Button(self)
        self.ButtonNext.place(relx=0.433, rely=0.850, height=25, width=70)
        self.ButtonNext.configure(activebackground="#f9f9f9")
        self.ButtonNext.configure(text='Next')
        self.ButtonNext.configure(command=self.confirm)

        self.oldDicts = ""
        self.VCFDir = ""
        self.sampleFile = ""

    def pickDictFile(self):
        # filedialog.askopenfilename(filetypes=(("json files","*.json"),("all files","*.*")))
        self.oldDicts = filedialog.askdirectory()
        self.TextDict.configure(state="normal")
        self.TextDict.delete(1.0, END)
        self.TextDict.insert(END, self.oldDicts)
        self.TextDict.configure(state="disabled")

    def pickVCFDir(self):
        self.VCFDir = filedialog.askdirectory()  # filedialog.askopenfilename()
        self.TextVCF.configure(state="normal")
        self.TextVCF.delete(1.0, END)
        self.TextVCF.insert(END, self.VCFDir)
        self.TextVCF.configure(state="disabled")

    def pickSampleFile(self):
        self.sampleFile = filedialog.askopenfilename()
        self.TextSample.configure(state="normal")
        self.TextSample.delete(1.0, END)
        self.TextSample.insert(END, self.sampleFile)
        self.TextSample.configure(state="disabled")

    def disableFrameName(self):
        if self.varCheck.get() == 0:
            for ele in self.FrameName.winfo_children():
                ele.configure(state='normal')
        elif self.varCheck.get() == 1:
            for ele in self.FrameName.winfo_children():
                ele.configure(state='disabled')

    # TODO to prevent pam with no guide len in the name, put 25 Ns , check if positions are correct
    def createTempPAM(self, pam, pamFile):
        parts = pamFile.split("-")
        if parts[0] != pam:
            with open(self.pathDir+"/pam/tempPAM.txt", "w") as tempPAM:
                fullPAM = "".join(["N"]*int(parts[0][0:2])) + \
                    pam + " " + str(len(pam))
                tempPAM.write(fullPAM)
        else:
            with open(self.pathDir+"/pam/tempPAM.txt", "w") as tempPAM:
                fullPAM = pam + \
                    "".join(["N"]*int(parts[1][0:2])) + " -" + str(len(pam))
                tempPAM.write(fullPAM)

    def confirm(self):
        if not self.oldDicts:
            self.LabelStatus.configure(fg='Red')
            self.LabelStatus.configure(text="Please select a dictionary")
            return
        elif not self.VCFDir:
            self.LabelStatus.configure(fg='Red')
            self.LabelStatus.configure(text="Please select a VCF")
            return
        elif not self.sampleFile:
            self.LabelStatus.configure(fg='Red')
            self.LabelStatus.configure(text="Please select a sampleID file")
            return
        """
        elif self.varCheck.get() == 0 and self.EntryName.get() == "":
            self.LabelStatus.configure(fg='Red')
            self.LabelStatus.configure(text="Please write a name for the new dictionary or overwrite the original")
            return
        if self.varCheck.get() == 1:
            print("#####################################")
            print("Overwriting old dictionary")
            print("#####################################")
            os.system("python file_per_crispritz/change_dict.py "+self.oldDicts+" "+self.VCFDir+" "+self.oldDicts)
        else:
            print("#####################################")
            print("Writing updated dictionary")
            print("#####################################")
            os.system("python change_dict.py "+self.oldDicts+" "+self.VCFDir+" ../dictionaries/"
                      +self.EntryName.get())
        """
        ww = WaitingWindow(self.root)
        ww.lift()
        ww.update()

        print("#####################################")
        print("Updating dictionaries")
        print("#####################################")
        dict_vcf = []
        for i, VCFFile in enumerate(os.listdir(self.VCFDir)):
            with gzip.open(self.VCFDir+"/"+VCFFile, 'rb') as targets:
                # Skip vcf header
                for line in targets:
                    line = line.decode('ascii')
                    if ('#CHROM') in line:
                        # Save this header for retrieving sample id
                        column_vcf = line.strip().split('\t')
                        break

                # Save CHROM [0], POS[1], REF [3], ALT [4], List of Samples [9:]
                line = targets.readline()
                line = line.decode('ascii').strip().split('\t')
                chrom = line[0]
                # print(chrom)
                if "chr" not in chrom:
                    chrom = "chr"+chrom
                for oldDict in os.listdir(self.oldDicts):
                    if oldDict[8:len(oldDict)-5] == chrom:
                        dict_vcf.append((oldDict, VCFFile))
                        break
            if len(dict_vcf) != i+1:
                print("WARNING: No dictionaries associated to "+VCFFile)
                messagebox.showerror(
                    "WARNING!", "No dictionaries associated to "+VCFFile)
                ww.destroy()
                return
        # print(dict_vcf)

        for oldDict, VCFFile in dict_vcf:
            print("Updating "+oldDict)
            updateDictionary(self.oldDicts+"/"+oldDict,
                             self.VCFDir+"/"+VCFFile)
            # os.system("python change_dict.py "+self.oldDicts+"/"+oldDict+" "+self.VCFDir+"/"+VCFFile)

        ww.doneLabel("dict")

        os.chdir(self.pathDir)
        print("#####################################")
        print("Updating enriched genome with new variants")
        print("#####################################")
        genomeEnr = os.path.basename(self.oldDicts).split('dictionary_')[-1]
        for item in os.listdir(self.pathDir+"/Genomes/"+genomeEnr+"/"):
            chrom = item.split('.')[0]
            os.rename(self.pathDir+"/Genomes/"+genomeEnr+"/"+item,
                      self.pathDir+"/Genomes/"+genomeEnr+"/"+chrom+".fa")
        os.system("crispritz.py add-variants "+self.VCFDir +
                  "/ "+self.pathDir+"/Genomes/"+genomeEnr+"/")
        rmtree(self.pathDir+"/Genomes/"+genomeEnr+"/")
        # if not os.path.exists("Genomes/"+genomeEnr):
        #     os.mkdir(self.pathDir+"/Genomes/"+genomeEnr)
        # for item in os.listdir(self.pathDir+"/variants_genome/SNPs_genome/"):
        #     copy(self.pathDir+"/variants_genome/SNPs_genome/"+item, self.pathDir+"/Genomes/"+genomeEnr+"/"+item)
        move(self.pathDir+"/variants_genome/SNPs_genome/" +
             genomeEnr + '_enriched', self.pathDir+"/Genomes/"+genomeEnr)
        # rmtree("variants_genome/")
        ww.doneLabel("var")

        print("#####################################")
        print("Indexing enriched genome with new variants")
        print("#####################################")
        pams, bMaxs = getAssociatedPAM(genomeEnr, self.pathDir)
        for pam, bMax in zip(pams, bMaxs):
            fullPAM = ""
            for d in os.listdir(self.pathDir+"/pam/"):
                if pam in d:
                    fullPAM = d
                    break
            if fullPAM == "":
                print("PAM NOT FOUND!")
                raise("PAM not found")
                return
            # print(fullPAM)
            self.createTempPAM(pam, fullPAM)
            os.system("crispritz.py index-genome " +
                      genomeEnr+" " +
                      self.pathDir+"/Genomes/"+genomeEnr+"/ " +
                      self.pathDir+"/pam/tempPAM.txt"+" -bMax "+str(bMax))
            os.remove(self.pathDir+"/pam/tempPAM.txt")
        ww.doneLabel("idx")

        print("#####################################")
        print("Adding samplesID")
        print("#####################################")
        with open(self.pathDir+'/samplesID/samples_'+genomeEnr+".txt", 'a') as oldSampleFile:
            with open(self.sampleFile, 'r') as newSampleFile:
                newSampleFile.readline()  # header
                for line in newSampleFile:
                    oldSampleFile.write("\n"+line.strip())
        ww.doneLabel("sample")

        print("#####################################")
        print("Procedure Finished")
        print("#####################################")
        ww.placeDismiss()

    def closeAll(self):
        self.destroy()
        self.quit()
        self.root.deiconify()
        self.root.destroy()
        self.root.quit()


if __name__ == '__main__':
    startUpdateDict(sys.argv[1])
