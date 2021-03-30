#! /usr/bin/env python
#  -*- coding: utf-8 -*-

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

from tkinter import filedialog, END
from shutil import copy, rmtree, copytree, move
import os
import sys
from os import listdir                      #for getting directories
from os.path import isfile, isdir,join      #for getting directories

class WaitingWindow(tk.Toplevel):
    def __init__(self, top=None, alertVCF=None):
        
        tk.Toplevel.__init__(self,top)
        self.root = top
        
        self.geometry("600x344+650+150")
        self.minsize(1, 1)
        self.maxsize(1825, 970)
        self.resizable(1, 1)
        self.title("Progress")

        self.FrameCopy = tk.Frame(self)
        self.FrameCopy.place(relx=0.015, rely=0.029, relheight=0.131
                , relwidth=0.975)
        self.FrameCopy.configure(relief='groove')
        self.FrameCopy.configure(borderwidth="2")
        self.FrameCopy.configure(relief="groove")

        self.LabelFiles = tk.Label(self.FrameCopy)
        self.LabelFiles.place(relx=0.017, rely=0.222, height=25, width=109)
        self.LabelFiles.configure(text='''Copying files''')

        self.LabelFilesStatus = tk.Label(self.FrameCopy)
        self.LabelFilesStatus.place(relx=0.632, rely=0.222, height=25, width=118)

        self.LabelFilesStatus.configure(text='''Pending...''')

        self.FrameVariants = tk.Frame(self)
        self.FrameVariants.place(relx=0.017, rely=0.174, relheight=0.131
                , relwidth=0.975)
        self.FrameVariants.configure(relief='groove')
        self.FrameVariants.configure(borderwidth="2")
        self.FrameVariants.configure(relief="groove")

        self.LabelVariants = tk.Label(self.FrameVariants)
        self.LabelVariants.place(relx=0.017, rely=0.222, height=25, width=119)
        self.LabelVariants.configure(text='''Adding Variants''')

        self.LabelVariantsStatus = tk.Label(self.FrameVariants)
        self.LabelVariantsStatus.place(relx=0.632, rely=0.222, height=25
                , width=118)
        if alertVCF:
            self.LabelVariantsStatus.configure(text='''SKIPPED''')
        else:
            self.LabelVariantsStatus.configure(text='''Pending...''')

        self.FrameRefGen = tk.Frame(self)
        self.FrameRefGen.place(relx=0.017, rely=0.32, relheight=0.131, relwidth=0.975)

        self.FrameRefGen.configure(relief='groove')
        self.FrameRefGen.configure(borderwidth="2")
        self.FrameRefGen.configure(relief="groove")

        self.LabelRefGen = tk.Label(self.FrameRefGen)
        self.LabelRefGen.place(relx=0.015, rely=0.222, height=25, width=189)
        self.LabelRefGen.configure(text='''Indexing reference genome''')

        self.LabelRefGenStatus = tk.Label(self.FrameRefGen)
        self.LabelRefGenStatus.place(relx=0.632, rely=0.222, height=25
                , width=118)
        self.LabelRefGenStatus.configure(text='''Pending...''')

        self.FrameEnrGen = tk.Frame(self)
        self.FrameEnrGen.place(relx=0.017, rely=0.465, relheight=0.131
                , relwidth=0.975)
        self.FrameEnrGen.configure(relief='groove')
        self.FrameEnrGen.configure(borderwidth="2")
        self.FrameEnrGen.configure(relief="groove")

        self.LabelEnrGen = tk.Label(self.FrameEnrGen)
        self.LabelEnrGen.place(relx=0.017, rely=0.222, height=25, width=179)
        self.LabelEnrGen.configure(text='''Indexing enriched genome''')

        self.LabelEnrGenStatus = tk.Label(self.FrameEnrGen)
        self.LabelEnrGenStatus.place(relx=0.632, rely=0.222, height=25
                , width=118)
        if alertVCF:
            self.LabelEnrGenStatus.configure(text='''SKIPPED''')
        else:
            self.LabelEnrGenStatus.configure(text='''Pending...''')

        self.FrameDict = tk.Frame(self)
        self.FrameDict.place(relx=0.017, rely=0.61, relheight=0.131, relwidth=0.975)

        self.FrameDict.configure(relief='groove')
        self.FrameDict.configure(borderwidth="2")
        self.FrameDict.configure(relief="groove")

        self.LabelDict = tk.Label(self.FrameDict)
        self.LabelDict.place(relx=0.017, rely=0.222, height=25, width=160)
        self.LabelDict.configure(text='''Writing dictionaries''')

        self.LabelDictStatus = tk.Label(self.FrameDict)
        self.LabelDictStatus.place(relx=0.632, rely=0.222, height=25, width=118)
        if alertVCF:
            self.LabelDictStatus.configure(text='''SKIPPED''')
        else:
            self.LabelDictStatus.configure(text='''Pending...''')

        self.ButtonDismiss = tk.Button(self)
        #self.ButtonDismiss.place(relx=0.433, rely=0.814, height=35, width=80)
        self.ButtonDismiss.configure(text='''Dismiss''')

    def doneLabel(self, label):
        if label == "copy":
            label = self.LabelFilesStatus
        elif label == "var":
            label = self.LabelVariantsStatus
        elif label == "ref":
            label = self.LabelRefGenStatus
        elif label == "enr":
            label = self.LabelEnrGenStatus
        elif label == "dict":
            label = self.LabelDictStatus
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

def startChooseFiles(datiFolder, initialDir):
    root = tk.Tk()
    root.withdraw()
    #top1 = WaitingWindow(root)
    top = TopChooseFiles(root, datiFolder, initialDir)#os.path.dirname(os.path.abspath(__file__))
    root.mainloop()
    root.quit()
            
    
class TopChooseFiles(tk.Toplevel):
    def __init__(self, top=None, pathDir=None, initialDir=None):
        
        self.initialDir = initialDir
        self.thisTop = tk.Toplevel.__init__(self, top)
        
        self.geometry("600x480+624+239")
        self.minsize(1, 1)
        self.maxsize(1825, 970)
        self.resizable(1, 1)
        self.title("Add Genome")

        self.protocol("WM_DELETE_WINDOW", self.closeAll)
        self.pathDir = pathDir
        self.root = top
        
        os.chdir(self.pathDir)
        
        self.FrameStatus = tk.Frame(self)
        self.FrameStatus.place(relx=0.033, rely=0.011, relheight=0.036
                , relwidth=0.942)

        self.LabelStatus = tk.Label(self.FrameStatus)
        self.LabelStatus.place(relx=0.012, rely=0.24, height=15, width=549)
        self.LabelStatus.configure(text="Please select the following directories and files")

        self.FrameGenome = tk.Frame(self)
        self.FrameGenome.place(relx=0.033, rely=0.057, relheight=0.224
                , relwidth=0.942)
        
        """
        self.LabelGenome = tk.Label(self.FrameGenome)
        self.LabelGenome.place(relx=0.012, rely=0.003, height=18, width=249)
        self.LabelGenome.configure(text="Select reference genome directory")
        """
        
        self.ButtonGenome = tk.Button(self.FrameGenome)
        self.ButtonGenome.place(relx=0.023, rely=0.242, height=35, width=125)
        self.ButtonGenome.configure(text="Choose Reference\nGenome")
        self.ButtonGenome.configure(command=self.pickGenomeDir)

        self.TextGenome = tk.Text(self.FrameGenome)
        self.TextGenome.place(relx=0.255, rely=0.242, relheight=0.3
                , relwidth=0.772)
        self.TextGenome.configure(background="white")
        self.TextGenome.configure(font="TkTextFont")
        self.TextGenome.configure(selectbackground="#c4c4c4")
        self.TextGenome.configure(state='disabled')
        self.TextGenome.configure(wrap="word")
        
        self.comboGenRef = ttk.Combobox(self.FrameGenome)
        self.comboGenRef.place(relx=0.023, rely=0.58, height=27)
        gen_refs = [""]
        for sub in os.listdir(self.pathDir+"/Genomes"):
            if '+' not in sub and sub[0] != ".":
                gen_refs.append(sub)
        self.comboGenRef.configure(values=gen_refs)
        self.comboGenRef.bind("<<ComboboxSelected>>", self.pickComboGenome)
        self.comboGenRef.configure(state='readonly')

        self.FrameVCF = tk.Frame(self)
        self.FrameVCF.place(relx=0.033, rely=0.24, relheight=0.124
                , relwidth=0.942)
        
        """
        self.LabelVCF = tk.Label(self.FrameVCF)
        self.LabelVCF.place(relx=0.012, rely=0.123, height=18, width=155)
        self.LabelVCF.configure(text="Select VCF directory")
        """
        
        self.ButtonVCF = tk.Button(self.FrameVCF)
        self.ButtonVCF.place(relx=0.023, rely=0.402, height=35, width=125)
        self.ButtonVCF.configure(text="Choose VCFs")
        self.ButtonVCF.configure(command=self.pickVCFDir)
        
        self.TextVCF = tk.Text(self.FrameVCF)
        self.TextVCF.place(relx=0.255, rely=0.4, relheight=0.538
                , relwidth=0.772)
        self.TextVCF.configure(background="white")
        self.TextVCF.configure(font="TkTextFont")
        self.TextVCF.configure(selectbackground="#c4c4c4")
        self.TextVCF.configure(state='disabled')
        self.TextVCF.configure(wrap="word")

        self.FramePAM = tk.Frame(self)
        self.FramePAM.place(relx=0.033, rely=0.370, relheight=0.124
                , relwidth=0.942)

        """
        self.LabelPAM = tk.Label(self.FramePAM)
        self.LabelPAM.place(relx=0.012, rely=0.123, height=18, width=125)
        self.LabelPAM.configure(text="Select PAM file")
        """
    
        self.ButtonPAM = tk.Button(self.FramePAM)
        self.ButtonPAM.place(relx=0.023, rely=0.402, height=35, width=125)
        self.ButtonPAM.configure(text="Choose PAM")
        self.ButtonPAM.configure(command=self.pickPAMFile)

        self.TextPAM = tk.Text(self.FramePAM)
        self.TextPAM.place(relx=0.255, rely=0.4, relheight=0.538
                , relwidth=0.772)
        self.TextPAM.configure(background="white")
        self.TextPAM.configure(font="TkTextFont")
        self.TextPAM.configure(selectbackground="#c4c4c4")
        self.TextPAM.configure(state='disabled')
        self.TextPAM.configure(wrap="word")

        self.FrameAnnotation = tk.Frame(self)
        self.FrameAnnotation.place(relx=0.033, rely=0.495, relheight=0.124
                , relwidth=0.942)
        
        """
        self.LabelAnnotation = tk.Label(self.FrameAnnotation)
        self.LabelAnnotation.place(relx=0.012, rely=0.123, height=18, width=170)
        self.LabelAnnotation.configure(text="Select annotation file")
        """

        self.ButtonAnnotation = tk.Button(self.FrameAnnotation)
        self.ButtonAnnotation.place(relx=0.023, rely=0.402, height=35, width=125)
        self.ButtonAnnotation.configure(text="Choose Annotation")
        self.ButtonAnnotation.configure(command=self.pickAnnotationFile)

        self.TextAnnotation = tk.Text(self.FrameAnnotation)
        self.TextAnnotation.place(relx=0.255, rely=0.4, relheight=0.538
                , relwidth=0.772)
        self.TextAnnotation.configure(background="white")
        self.TextAnnotation.configure(font="TkTextFont")
        self.TextAnnotation.configure(selectbackground="#c4c4c4")
        self.TextAnnotation.configure(state='disabled')
        self.TextAnnotation.configure(wrap="word")

        self.FrameSampleList = tk.Frame(self)
        self.FrameSampleList.place(relx=0.033, rely=0.620, relheight=0.124
                , relwidth=0.942)
        
        """
        self.LabelSampleList = tk.Label(self.FrameSampleList)
        self.LabelSampleList.place(relx=0.012, rely=0.123, height=18, width=180)
        self.LabelSampleList.configure(text="Select sample list file")
        """
        
        self.ButtonSampleList = tk.Button(self.FrameSampleList)
        self.ButtonSampleList.place(relx=0.023, rely=0.402, height=35, width=125)
        self.ButtonSampleList.configure(text="Choose SamplesID")
        self.ButtonSampleList.configure(command=self.pickSampleListFile)
        
        self.TextSampleList = tk.Text(self.FrameSampleList)
        self.TextSampleList.place(relx=0.255, rely=0.4, relheight=0.538
                , relwidth=0.772)
        self.TextSampleList.configure(background="white")
        self.TextSampleList.configure(font="TkTextFont")
        self.TextSampleList.configure(selectbackground="#c4c4c4")
        self.TextSampleList.configure(state='disabled')
        self.TextSampleList.configure(wrap="word")
        
        self.disableFrame(self.FrameSampleList)
        
        self.FrameBMax = tk.Frame(self)
        self.FrameBMax.place(relx=0.182, rely=0.745, relheight=0.1
                , relwidth=0.658)
        
        self.LabelBMax = tk.Label(self.FrameBMax)
        self.LabelBMax.place(relx=0.003, rely=0.222, height=25, width=159)
        self.LabelBMax.configure(text="Choose # of bulges:")
        
        self.comboBMax = ttk.Combobox(self.FrameBMax)
        self.comboBMax.place(relx=0.403, rely=0.222, height=27)
        self.comboBMax.configure(values=[0,1,2])
        self.comboBMax.configure(state='readonly')

        self.FrameNameGenome = tk.Frame(self)
        self.FrameNameGenome.place(relx=0.182, rely=0.835, relheight=0.1
                , relwidth=0.658)

        self.LabelGenomeName = tk.Label(self.FrameNameGenome)
        self.LabelGenomeName.place(relx=0.003, rely=0.222, height=25, width=159)
        self.LabelGenomeName.configure(text="Enriched genome name:")

        self.EntryGenomeName = tk.Entry(self.FrameNameGenome)
        self.EntryGenomeName.place(relx=0.403, rely=0.222, height=27
                , relwidth=0.572)
        self.EntryGenomeName.configure(background="white")
        self.EntryGenomeName.configure(font="TkFixedFont")
        
        self.disableFrame(self.FrameNameGenome)

        self.ButtonConfirm = tk.Button(self)
        self.ButtonConfirm.place(relx=0.443, rely=0.928, height=25, width=56)
        self.ButtonConfirm.configure(text="Next")
        self.ButtonConfirm.configure(command=self.confirm)
        
        self.genomeDir = ""
        self.VCFDir = ""
        self.PAMFile = ""
        self.annotationFile = ""
        self.sampleFile = ""
        

    def disableFrame(self, frame):
        for ele in frame.winfo_children():
                ele.configure(state='disabled')
                
    def enableFrame(self, frame):
        for ele in frame.winfo_children():
                ele.configure(state='normal')
        
    def pickGenomeDir(self):
        self.genomeDir = filedialog.askdirectory()
        self.TextGenome.configure(state="normal")
        self.TextGenome.delete(1.0, END)
        self.TextGenome.insert(END, self.genomeDir)
        self.TextGenome.configure(state="disabled")
        if self.genomeDir != "":
            self.comboGenRef.configure(state='disabled')
        else:
            self.comboGenRef.configure(state='readonly')
            
    def pickComboGenome(self, event):
        self.genomeDir = self.pathDir+"/Genomes/"+self.comboGenRef.get()
        if self.comboGenRef.get() != "":
            self.TextGenome.configure(state="disabled")
            self.ButtonGenome.configure(state="disabled")
        else:
            self.genomeDir = ""
            #self.TextGenome.configure(state="normal")
            self.ButtonGenome.configure(state="normal")
        
    def pickVCFDir(self):
        self.VCFDir = filedialog.askdirectory()
        self.TextVCF.configure(state="normal")
        self.TextVCF.delete(1.0, END)
        self.TextVCF.insert(END, self.VCFDir)
        self.TextVCF.configure(state="disabled")
        if not self.VCFDir:
            self.disableFrame(self.FrameSampleList)
            self.disableFrame(self.FrameNameGenome)
        else:
            self.enableFrame(self.FrameSampleList)
            self.enableFrame(self.FrameNameGenome)
        
    def pickPAMFile(self):
        self.PAMFile = filedialog.askopenfilename()
        self.TextPAM.configure(state="normal")
        self.TextPAM.delete(1.0, END)
        self.TextPAM.insert(END, self.PAMFile)
        self.TextPAM.configure(state="disabled")
    
    def pickAnnotationFile(self):
        self.annotationFile = filedialog.askopenfilename(initialdir = self.pathDir + 'annotations/')
        self.TextAnnotation.configure(state="normal")
        self.TextAnnotation.delete(1.0, END)
        self.TextAnnotation.insert(END, self.annotationFile)
        self.TextAnnotation.configure(state="disabled")
    
    def pickSampleListFile(self):
        self.sampleFile = filedialog.askopenfilename()
        self.TextSampleList.configure(state="normal")
        self.TextSampleList.delete(1.0, END)
        self.TextSampleList.insert(END, self.sampleFile)
        self.TextSampleList.configure(state="disabled")
        
    def confirm(self):
        alertVCF = False
        if not self.VCFDir:
            alertVCF = True
        if not self.genomeDir:
            self.LabelStatus.configure(fg='Red')
            self.LabelStatus.configure(text="Please select a genome directory")
            return
        elif not self.PAMFile:
            self.LabelStatus.configure(fg='Red')
            self.LabelStatus.configure(text="Please select a PAM file")
            return
        elif not self.annotationFile:
            self.LabelStatus.configure(fg='Red')
            self.LabelStatus.configure(text="Please select an annotation file")
            return
        elif not alertVCF and not self.sampleFile:
            self.LabelStatus.configure(fg='Red')
            self.LabelStatus.configure(text="Please select a sample ID file")
            return
        elif not self.comboBMax.get():
            self.LabelStatus.configure(fg='Red')
            self.LabelStatus.configure(text="Please select a # of bulges")
            return
        elif not alertVCF and not self.EntryGenomeName.get():
            self.LabelStatus.configure(fg='Red')
            self.LabelStatus.configure(text="Please write a name for the enriched genome")
            return
        
        
        bMax = self.comboBMax.get()
        if self.comboGenRef.get() != "":
            cleanName = self.comboGenRef.get()[0:len(self.comboGenRef.get())-4]
        else:
            cleanName = ''.join([char for char in os.path.basename(os.path.normpath(self.genomeDir)) if char != "+"])
        if not alertVCF:
            for f in os.listdir(self.pathDir+"/Genomes/"):
                if "+" in f:
                    genEnr = f.split("+")[1]
                    if cleanName+"_"+self.EntryGenomeName.get() == genEnr:
                        self.LabelStatus.configure(fg='Red')
                        self.LabelStatus.configure(text="The name for the enriched genome is already used, please change it")
                        return
        """
        elif self.comboBMax.get() not in ("0","1","2"):
            self.LabelStatus.configure(fg='Red')
            self.LabelStatus.configure(text="Please select a # of bulges between 0, 1 and 2")
            return
        """
        
        if alertVCF:
            choice = tk.messagebox.askquestion("Add genome", "Do you want to add only the reference genome without variants?", icon='warning')
            if choice != "yes":
                self.LabelStatus.configure(fg='Red')
                self.LabelStatus.configure(text="Please select a VCF directory")
                return
            
        ww = WaitingWindow(self.root, alertVCF)
        ww.lift()
        ww.update()
        
        os.chdir(self.pathDir)
        if not os.path.exists(self.pathDir+"/Genomes/"+cleanName+"_ref"):
            os.mkdir(self.pathDir+"/Genomes/"+cleanName+"_ref")
        for item in os.listdir(self.genomeDir+"/"):
            if not os.path.exists(self.pathDir+"/Genomes/"+cleanName+
                 "_ref/"+item):
                copy(self.genomeDir+"/"+item, self.pathDir+"/Genomes/"+cleanName+
                     "_ref/"+item)
        # if not os.path.exists(self.pathDir+"/pam/"+os.path.basename(self.PAMFile)):
        #     copy(self.PAMFile, self.pathDir+"/pam/"+os.path.basename(self.PAMFile))
        with open (self.PAMFile) as pam:
            line = pam.read().strip()
            pam = line.split(' ')[0]
            len_pam = int(line.split(' ')[1])
            guide_len = len(pam) - len_pam
            pos_beg = 0
            pos_end = None
            pam_begin = 0
            pam_end = len_pam * (-1)
            if len_pam < 0:
                guide_len = len(pam) + len_pam
                pam = pam[: (len_pam * (-1))]
                len_pam = len_pam * (-1)
                pos_beg = len_pam
                pos_end = None
                pam_begin = 0
                pam_end = len_pam
                pam_at_beginning = True
            else:
                pam = pam[(len_pam * (-1)):]
                pos_beg = 0
                pos_end = len_pam * (-1)
                pam_begin = len_pam * (-1)
                pam_end = None
        pam_available = [f for f in listdir(self.pathDir + '/pam') if isfile(join(self.pathDir + '/pam', f))]
        if next ((s for s in pam_available if pam in s), None):  #Return None if the selected pam is not in the pam directory
            pass
        else:       #Add the pam file to the pam directory
            with open (self.pathDir+"/pam/"+os.path.basename(self.PAMFile), "w") as pam_file:
                pam_file.write(pam+" "+str(len_pam))
        if not os.path.exists(self.pathDir+"/annotations/"+cleanName+"_ref.annotations.bed"):
            copy(self.annotationFile, self.pathDir+"/annotations/"+cleanName+"_ref.annotations.bed")
        if not alertVCF and not os.path.exists(self.pathDir+"/samplesID/samples_"+cleanName+"_ref+"+cleanName+"_"+self.EntryGenomeName.get()+".txt"):
            copy(self.sampleFile, self.pathDir+"/samplesID/samples_"+cleanName+"_ref+"+cleanName+"_"+self.EntryGenomeName.get()+".txt")
        ww.doneLabel("copy")
        if not alertVCF:
            print("#####################################")
            print("Creating enriched genome")
            print("#####################################")
            os.system("crispritz.py add-variants "+self.VCFDir+"/ "+self.genomeDir+"/")
            # if not os.path.exists(self.pathDir + "Genomes/"+cleanName+"_ref+"+cleanName+"_"+self.EntryGenomeName.get()):
            #     os.mkdir(self.pathDir+"/Genomes/"+cleanName+"_ref+"+cleanName+"_"+self.EntryGenomeName.get())
            # for item in os.listdir(self.pathDir+"/variants_genome/SNPs_genome/"):
            #     copytree(self.pathDir+"/variants_genome/SNPs_genome/"+item, self.pathDir+"/Genomes/"+
            #          cleanName+"_ref+"+cleanName+"_"+self.EntryGenomeName.get()+"/"+item)
            move(self.pathDir+"/variants_genome/SNPs_genome/" + cleanName + '_enriched', self.pathDir+"/Genomes/"+
                     cleanName+"_ref+"+cleanName+"_"+self.EntryGenomeName.get())
            rmtree("variants_genome/")
            ww.doneLabel("var")
        print("#####################################")
        print("Creating indexes for reference genome")
        print("#####################################")
        os.system("crispritz.py index-genome "+
                  cleanName+"_ref "+
                  self.pathDir+"/Genomes/"+cleanName+"_ref/ "+
                  self.PAMFile+" -bMax "+str(bMax))
        ww.doneLabel("ref")
        if not alertVCF:
            print("#####################################")
            print("Creating indexes for enriched genome")
            print("#####################################")
            os.system("crispritz.py index-genome "+cleanName+"_ref+"+cleanName+"_"+self.EntryGenomeName.get()+
                      " "+self.pathDir+"/Genomes/"+cleanName+"_ref+"+cleanName+"_"+self.EntryGenomeName.get()+
                      "/ "+self.PAMFile+" -bMax "+str(bMax))
            ww.doneLabel("enr")
            print("#####################################")
            print("Creating dictionaries")
            print("#####################################")
            os.chdir(self.initialDir)
            if not os.path.exists(self.pathDir+"/dictionaries/dictionary_"+cleanName+"_ref+"+cleanName+"_"+self.EntryGenomeName.get()):
                os.mkdir(self.pathDir+"/dictionaries/dictionary_"+cleanName+"_ref+"+cleanName+"_"+self.EntryGenomeName.get())
            for item in os.listdir(self.VCFDir+"/"):
                print("Creating dictionary for file "+item)
                os.system("python creazione_dizionari.py "+self.VCFDir+"/"+item+" "+
                          self.pathDir+"/dictionaries/dictionary_"+cleanName+"_ref+"+cleanName+
                          "_"+self.EntryGenomeName.get()+"/"+item)
            ww.doneLabel("dict")
        # with open(self.pathDir+"/logGenomes/"+os.path.basename(os.path.normpath(self.genomeDir))+"+"+
        #                self.EntryGenomeName.get()+".log","w") as logFile:
        #     logFile.writelines("Reference Genome\t"+self.genomeDir+"\n")
        #     logFile.writelines("VCF\t"+self.VCFDir+"\n")
        #     logFile.writelines("PAM\t"+self.PAMFile+"\n")
        #     logFile.writelines("Annotation\t"+self.annotationFile+"\n")
        #     logFile.writelines("Sample List\t"+self.sampleFile+"\n")
        #     logFile.writelines("# Bulges\t"+str(bMax)+"\n")
        #     logFile.writelines("Name Enriched\t"+self.EntryGenomeName.get()+"\n")
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
    initialDir = os.path.dirname(os.path.abspath(__file__))
    startChooseFiles(sys.argv[1], initialDir)





