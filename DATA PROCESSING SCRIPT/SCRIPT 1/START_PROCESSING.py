"""
KineMed Automated Data Extraction Script (Script 1 of 2 in Workflow)
 
This script will generate the following files/folders in the KineMed data 
processing pipeline:
    
    1) Mass Hunter Method (Created in the 'C:\Methods' directory)
    2) MIDA database file (Created in the 'C:\Peaklists\DB files' directory)
    3) A new directory for subsequent compound reports with:
        a) A copy of the original peptide list
        b) A copy of the data processing script
        c) A copy of the MIDA database
        d) A 'Parameters Template' file, with all sample and parameter details
           that is needed for the data processing script
        e) A note of the script version used
        f) A note if the retention time has been modified from the
           value (0.8 min)
           
The script requires 3 specific files within this directory in order to be run:
    
    1) A copy of the Spectrum Mill peptide list
    2) A copy of the sample submission form (version 2 or higher)
    3) The 'SETTINGS' file (All analysis parameters can be modified here)
    
"""

# Import all modules
from __future__ import division
import pandas as pd
import os
from collections import OrderedDict
from subprocess import call,CREATE_NEW_CONSOLE
import sys
import xml.etree.ElementTree as ET
import shutil
import xlrd
from datetime import date

# Get the type of method to be used, depending on the given setting
def getMethod(drive,setting):
    return {
        '5K': os.path.join(drive,'Methods','FBF011816_5K.m'), #'FBF031714_5K.m'
        '10K': os.path.join(drive,'Methods','FBF011816_10K.m'), #'FBF031714_5K.m'
        '15K': os.path.join(drive,'Methods','FBF011816_15K.m'), #'FBF031714_5K.m'
        '30K': os.path.join(drive,'Methods','FBF011816_30K.m') #'FBF031714_5K.m'
    }.get(setting, '') 

# Get specific version and settings to be used for MIDA script, depending on the
# given setting
def getSettings(setting):
    return {
        'Default': [6.22,50,1,False,0.1,6,10,'generate_emp.py'],
        'D20 0-7': [6.22,12,7,False,0.1,8,15,'generate_emp_BW07.py'],
        'Low': [6.22,12,6,False,0.1,6,12,'generate_emp.py'],
        'Lumican': [6.22,12,7,False,0.1,5,5,'generate_emp.py'],
        'In Vitro': [6.22,12,7,False,0.1,8,12,'generate_emp.py'],
        'D3 Leucine': [6.22,12,7,False,0.1,8,0,'generate_emp_d3.py'],
        'CKMM': [6.22,12,6,False,0.1,7,7,'generate_emp.py'],
        'Acetyl': [6.22,12,7,True,0.1,8,0,'generate_emp_acetyl.py'],
        '15N': [6.22,12,7,False,0.1,8,10,'generate_emp_15N.py']
    }.get(setting, [6.22,12,7,False,0.1,8,15,'generate_emp.py']) 
    
# Create the MIDA database
def createDatabase(mergedList):
    # Assign values from peptide list to variable names
    ppmDelta = mergedList['delta_parent_mass_ppm'].values
    checkScore = mergedList['score'].values
    fwdRevScore = mergedList['deltaForwardReverseScore'].values/mergedList['score'].values
    checkN = mergedList['n'].values
    checkMod = mergedList['modifications'].fillna('none').values
    numAA = mergedList['sequence'].map(len)
    
    # Filter peptide list 
    mergedList['PASS'] = ['OK' if 
        ppmDelta[i]>-maxPPMError and 
        ppmDelta[i]<maxPPMError and 
        checkScore[i]>=minPeptideScore and 
        fwdRevScore[i]>=fwdRevScoreThreshold and
        checkN[i]>=minN and
        numAA[i] >= minAA and 
        (useAcetyl==True or checkMod[i].find('Acetyl')==-1)
        else '' for i in xrange(len(mergedList))]

    mergedList = mergedList[mergedList['PASS']=='OK'].reset_index()
    
    # Assign column values
    formula = mergedList.groupby(['accession_number','sequence'],sort=False)['composition'].first()
    rt = mergedList.groupby(['accession_number','sequence'],sort=False)['retentionTimeMin'].mean().round(2).values
    matchedMass = mergedList.groupby(['accession_number','sequence'],sort=False)['matched_parent_mass'].mean().values
    cpd = mergedList.groupby(['accession_number','sequence'],sort=False)['accession_number'].first() + "+" + mergedList.groupby(['accession_number','sequence'],sort=False)['entry_name'].first()
    sequence = mergedList.groupby(['accession_number','sequence'],sort=False)['sequence'].first()
    newSequence = [sequence[i].replace('p', 'xp').replace('q', 'xq').replace('m', 'xm').replace('k', 'xk') for i in xrange(len(sequence))]
    mass = [matchedMass[i]-1.0073+(newSequence[i].count('k')*42.011)+(newSequence[i].count('m')*15.995)+(newSequence[i].count('q')*-17.027)+(newSequence[i].count('p')*15.995) for i in xrange(len(newSequence))]
    maxAA = mergedList.groupby(['accession_number','sequence'],sort=False)['StartAA'].max().values
    modifications = mergedList.groupby(['accession_number','sequence'],sort=False)['modifications'].first().fillna('none').map(str)
    notes = [str(maxAA[i])+'+'+newSequence[i]+'-'+modifications[i] for i in xrange(len(newSequence))]
    score = mergedList.groupby(['accession_number','sequence'],sort=False)['score'].max().values
    composition = formula
    n = mergedList.groupby(['accession_number','sequence'],sort=False)['n'].max().values
    m0 = mergedList.groupby(['accession_number','sequence'],sort=False)['M0'].max().values
    m1 = mergedList.groupby(['accession_number','sequence'],sort=False)['M1'].max().values
    m2 = mergedList.groupby(['accession_number','sequence'],sort=False)['M2'].max().values
    m3 = mergedList.groupby(['accession_number','sequence'],sort=False)['M3'].max().values
    m4 = mergedList.groupby(['accession_number','sequence'],sort=False)['M4'].max().values
    em0cc3 = mergedList.groupby(['accession_number','sequence'],sort=False)['EM0 cubic coeff 3'].max().values
    em0cc2 = mergedList.groupby(['accession_number','sequence'],sort=False)['EM0 cubic coeff 2'].max().values
    em0cc1 = mergedList.groupby(['accession_number','sequence'],sort=False)['EM0 cubic coeff 1'].max().values
    em1cc3 = mergedList.groupby(['accession_number','sequence'],sort=False)['EM1 cubic coeff 3'].max().values
    em1cc2 = mergedList.groupby(['accession_number','sequence'],sort=False)['EM1 cubic coeff 2'].max().values
    em1cc1 = mergedList.groupby(['accession_number','sequence'],sort=False)['EM1 cubic coeff 1'].max().values
    em2cc3 = mergedList.groupby(['accession_number','sequence'],sort=False)['EM2 cubic coeff 3'].max().values
    em2cc2 = mergedList.groupby(['accession_number','sequence'],sort=False)['EM2 cubic coeff 2'].max().values
    em2cc1 = mergedList.groupby(['accession_number','sequence'],sort=False)['EM2 cubic coeff 1'].max().values
    em3cc3 = mergedList.groupby(['accession_number','sequence'],sort=False)['EM3 cubic coeff 3'].max().values
    em3cc2 = mergedList.groupby(['accession_number','sequence'],sort=False)['EM3 cubic coeff 2'].max().values
    em3cc1 = mergedList.groupby(['accession_number','sequence'],sort=False)['EM3 cubic coeff 1'].max().values
    em4cc3 = mergedList.groupby(['accession_number','sequence'],sort=False)['EM4 cubic coeff 3'].max().values
    em4cc2 = mergedList.groupby(['accession_number','sequence'],sort=False)['EM4 cubic coeff 2'].max().values
    em4cc1 = mergedList.groupby(['accession_number','sequence'],sort=False)['EM4 cubic coeff 1'].max().values
    allSwissProt = mergedList.groupby(['accession_number','sequence'],sort=False)['accession_numbers'].first()
    multipleIDs = [1 if allSwissProt[i].find('|')!=-1 else 0 for i in xrange(len(allSwissProt))]

    # Define dictionary to construct dataframe to be output
    midaColumns = OrderedDict([
        ('Formula', formula),
        (' RT', rt),
        (' Mass', mass),
        (' Cpd', cpd),
        ('Notes', notes),
        ('sequence', sequence),
        ('modifications', modifications),
        ('score', score),
        ('composition', composition),
        ('n', n),
        ('M0', m0),
        ('M1', m1),
        ('M2', m2),
        ('M3', m3),
        ('M4', m4),
        ('EM0 cubic coeff 3', em0cc3),
        ('EM0 cubic coeff 2', em0cc2),
        ('EM0 cubic coeff 1', em0cc1),
        ('EM1 cubic coeff 3', em1cc3),
        ('EM1 cubic coeff 2', em1cc2),
        ('EM1 cubic coeff 1', em1cc1),
        ('EM2 cubic coeff 3', em2cc3),
        ('EM2 cubic coeff 2', em2cc2),
        ('EM2 cubic coeff 1', em2cc1),
        ('EM3 cubic coeff 3', em3cc3),
        ('EM3 cubic coeff 2', em3cc2),
        ('EM3 cubic coeff 1', em3cc1),
        ('EM4 cubic coeff 3', em4cc3),
        ('EM4 cubic coeff 2', em4cc2),
        ('EM4 cubic coeff 1', em4cc1),
        ('All Swissprot IDs', allSwissProt),
        ('multiple IDs (1=yes)', multipleIDs)])
        
    midaDatabase = pd.DataFrame(midaColumns)
    return midaDatabase

# Check peptide sequence for any invalud characters
def CheckSequence(sequence):
    if any(c in set('BZJX') for c in sequence):
        return True
    else:
        return False  

# Construct Mass Hunter method and output method and MIDA database
def OutputDatabaseAndMethod(drive,DBFile,Filename,MethodTemplate,RetentionTime,CRDirectory):
    DBFile.to_csv(os.path.join(drive,'Compound Reports',CRDirectory,Filename+'_DB.csv'),index = False)
    DBFile.to_csv(os.path.join(drive,'Peaklists\DB files',Filename+'_DB.csv'),index = False)
    
    # Determine new name for method
    NewMethodName = MethodTemplate.split('\\')[-1]
    NewMethodFilename = os.path.join(drive,'Methods',Filename+'_'+NewMethodName)
    if os.path.exists(NewMethodFilename) == False:
        shutil.copytree(MethodTemplate,NewMethodFilename)
    
    # Copy method template and modify settings
    xmlTree = ET.parse(os.path.join(NewMethodFilename,'DaMethod\Qual\qualitative.xml'))
    root = xmlTree.getroot()
    for parameter in root.iter('Parameter'):
        idVal = parameter.get('id')
        if idVal == 'RetentionTimeTolerance':
            parameter.find('Value').text = str(RetentionTime)
        if idVal == 'RTExpansionWidth':
            parameter.find('Value').text = str(RetentionTime*1.5)
        if idVal == 'DatabasePath':
            parameter.find('Value').text = str(os.path.join(drive,'Peaklists\DB files',Filename+'_DB.csv'))
        if idVal == 'PeakHeightAbsThreshold':
            parameter.find('Value').text = str(int(MethodTemplate.split('\\')[-1].split('_')[1].split('K')[0])*1000)
        if idVal == 'MethodName':
            parameter.find('Value').text = NewMethodName
        xmlTree.write(os.path.join(NewMethodFilename,'DaMethod\Qual\qualitative.xml'))
        # Create, save and remove dummy file so that modified date is updated
        DummyNote = open(os.path.join(NewMethodFilename,'dummy.txt'),'w')
        DummyNote.close()
        os.remove(os.path.join(NewMethodFilename,'dummy.txt'))
    return

# Create Parameters Template file based on information from the sample
# submission form    
def makeParametersFile(FileName,OutputFolder,ParameterTable):
    # Read data from sample submission form
    ExcelFile = xlrd.open_workbook(FileName, on_demand=True)
    SheetNames = [i for i in ExcelFile.sheet_names() if i.find('Page')!=-1]    
    ProjectDetails = ExcelFile.sheet_by_name(SheetNames[0])    
    instrument = [str(i) for i in ProjectDetails.row_values(0) if i.startswith('Q')][0]
    projectLeader = str(ProjectDetails.row_values(1)[2])
    processedBy = str(ProjectDetails.row_values(2)[2])
    projectCode = str(ProjectDetails.row_values(2)[-6])
    submitDate = str(date(*xlrd.xldate_as_tuple(ProjectDetails.row_values(3)[2],ExcelFile.datemode)[:3]))
    notebookCode = str(ProjectDetails.row_values(3)[-6])
    tissueFluid = str(ProjectDetails.row_values(4)[-6])
    prep = str(ProjectDetails.row_values(5)[-6])
    species = str(ProjectDetails.row_values(4)[2])[0]    
    ExcelFile.unload_sheet(SheetNames[0])
    #import pdb; pdb.set_trace()
    
    AllPages = []
    for i in xrange(len(SheetNames)):
        Page = pd.read_excel(FileName,sheetname=SheetNames[i],skiprows=17)
        Page.columns = [i.strip() for i in Page.columns]
        Page = Page[pd.notnull(Page['4 char Subject Code'])]
        AllPages.append(Page)
    
    SampleData = pd.concat(AllPages).reset_index()
    SampleData.columns = [i.strip() for i in SampleData.columns]  
    DefaultParameters = [[ParameterTable.index[i],ParameterTable.ix[i,'VALUE']] for i in xrange(len(ParameterTable))]
    DefaultParameters[5] = DefaultParameters[5] + [round(i,2) for i in ParameterTable.ix['M Difference Criteria',:].tolist()][1:]
   
    # Construct output formatting for Parameters Template
    SampleDetails = [
        ['Instrument:',instrument,'','','','',''],
        ['Project Leader:',projectLeader,'','','','',''],
        ['Processed By:',processedBy,'','','','',''],
        ['Submit Date:',submitDate,'','','','',''],
        ['Project Code & Title:',projectCode,'','','','',''],
        ['Notebook Code:',notebookCode,'','','','',''],
        ['Tissue/Fluid:',tissueFluid,'','','','',''],
        ['Prep',prep,'','','','',''],
        ['Species:',species,'','','','',''],
        ['','','','','','',''],
        ['#','8 char Sample Code','9 char Notebook Code','4 char Subject Code','Tissue & Fraction','Timepoint, treatment','%BW (exact %)']]

    for i in xrange(len(SampleData)):
        SampleDetails.append(['','','',str(SampleData.ix[i,4]),'',SampleData.ix[i,6],float(SampleData.ix[i,7])])
    
    # Save new Parameters Template file
    # Page 1 - Project/Sample information
    # Page 2 - Parameter Names/Values     
    writer = pd.ExcelWriter(os.path.join(OutputFolder,notebookCode +' Parameters Template.xlsx'))
    sheet1 = pd.DataFrame(SampleDetails,columns = None)
    sheet2 = pd.DataFrame(DefaultParameters, columns = None)
    sheet1.to_excel(writer,'Sample',index=False,header=False)
    sheet2.to_excel(writer,'Parameters',index=False,header=False)
    writer.save()
        
src = os.path.dirname(os.path.abspath(__file__))
# Intiialize computer drive
# Workstations 113/097 --> C:\
# Workstation 103 --> E:\ (External hard drive)
Drive = 'D:\\'

print "IMPORTING PEPTIDE LIST AND SAMPLE SUBMISSION FORM..."
# Import peptide list, sample submission form, settings file and data processing script
PeptideFileName = os.path.join(src,[i for i in os.listdir(src) if i[1]!='.' and (i.endswith('.txt') or i.endswith('.ssv'))][0])
SampleSubmissionFileName = os.path.join(src,[i for i in os.listdir(src) if i.endswith('.xls') and i.find('SETTINGS')==-1][0])
SettingsFileName = os.path.join(src,[i for i in os.listdir(src) if i.endswith('.xlsx') and i.find('SETTINGS')!=-1][0])
ScriptFileName = os.path.join(src,[i for i in os.listdir(src) if i.endswith('.py') and i.find('DataProcessingScript')!=-1][0])

# Identify current version of data processing script, to be used to create note
Current_Script_Version = ScriptFileName.split('\\')[-1].split('.py')[0].split('_')[1]
Update_Date = '09-08-16'

print "INITIALIZING MIDA DATABASE..."
# Create 'INPUT.csv' file to be passed to MIDA database script
peptideListOriginal = pd.read_table(PeptideFileName,delimiter=";")
peptideList = peptideListOriginal[peptideListOriginal['sequence'].map(CheckSequence)==False].reset_index()
midaInput = pd.DataFrame([peptideList.sequence[i] for i in xrange(len(peptideList.sequence))],columns=['sequence'])
midaInput.to_csv('INPUT.csv',index=False)

print "IMPORTING SETTINGS..."
# Read in settings from 'SETTING.xlsx' file
SettingsTable = pd.read_excel(SettingsFileName, sheetname='Settings',index_col=0)
MassHunterMethod = getMethod(Drive,SettingsTable.ix['Method','VALUE'])
MIDASettings = getSettings(SettingsTable.ix['MIDA','VALUE'])
RTWindow = float(SettingsTable.ix['RT Window','VALUE'])

rtMin = MIDASettings[0]
maxPPMError = MIDASettings[1]
minPeptideScore = MIDASettings[2]
useAcetyl = MIDASettings[3]
fwdRevScoreThreshold = MIDASettings[4]
minAA = MIDASettings[5]
minN = MIDASettings[6]
MidaScriptType = MIDASettings[7]

print "EXECUTING MIDA SCRIPT..."
# Call MIDA database script, passing 'INPUT'csv' and writing output to 'OUTPUT.csv'
call([sys.executable, MidaScriptType, 'INPUT.csv', 'OUTPUT.csv'], creationflags=CREATE_NEW_CONSOLE)

print "CREATING MIDA DATABASE..."
# Create and filter MIDA database
MidaOutputFileName = os.path.join(src,[i for i in os.listdir(src) if i.endswith('.csv') and i.upper().find('OUTPUT')!=-1][0])
midaList = pd.read_csv(MidaOutputFileName)
del peptideList['sequence']
mergedList = pd.concat([peptideList, midaList], axis=1)
midaDatabase = createDatabase(mergedList)
midaDatabase = midaDatabase.sort(columns=[' Mass','sequence',' Cpd'],ascending=[True,False,True]).reset_index(drop=True)
checkMass = ['']

for i in range(len(midaDatabase))[1:len(midaDatabase)-1]:
    if (midaDatabase[' Mass'][i] == midaDatabase[' Mass'][i-1] and (midaDatabase[' Cpd'][i][0:6] == midaDatabase[' Cpd'][i-1][0:6] or midaDatabase['sequence'][i] == midaDatabase['sequence'][i-1])):
        checkMass.append('')
    elif (midaDatabase[' Mass'][i] == midaDatabase[' Mass'][i-1] and abs(midaDatabase[' RT'][i-1] - midaDatabase[' RT'][i])<0.3):
        checkMass.append('BAD') 
    elif (midaDatabase[' Mass'][i] == midaDatabase[' Mass'][i+1] and (midaDatabase[' Cpd'][i][0:6] == midaDatabase[' Cpd'][i+1][0:6] or midaDatabase['sequence'][i] == midaDatabase['sequence'][i+1])):
        checkMass.append('')
    elif (midaDatabase[' Mass'][i] == midaDatabase[' Mass'][i+1] and abs(midaDatabase[' RT'][i] - midaDatabase[' RT'][i+1])<0.3):
        checkMass.append('BAD')
    else:
        checkMass.append('')

checkMass.append('')
midaDatabase['MASS CHECK'] = checkMass
midaDatabase = midaDatabase[midaDatabase['MASS CHECK']!='BAD'].reset_index(drop=True)
del midaDatabase['MASS CHECK']

# Determine project details, to be used for folder and file naming
NotebookCode = PeptideFileName.split('\\')[-1].split('_')[1].replace('-','')
Fraction = PeptideFileName.split('\\')[-1].split('_')[2]
Month = PeptideFileName.split('\\')[-1].split('_')[0].split('-')[1]
Day = PeptideFileName.split('\\')[-1].split('_')[0].split('-')[2]
Year = PeptideFileName.split('\\')[-1].split('_')[0].split('-')[0]
Instrument = PeptideFileName.split('\\')[-1].split('_')[3][0] + str(PeptideFileName.split('\\')[-1].split('_')[3][-1])
NewFolderName = '_'.join([NotebookCode,Fraction,'-'.join([Month,Day,Year[-2:]])])+Instrument

if os.path.exists(os.path.join(Drive,'Compound Reports',NewFolderName)):
    FolderList = [i[0] for i in os.walk(os.path.join(Drive,'Compound Reports')) if i[0].find(NewFolderName)!=-1]
    FolderList.sort()
    if FolderList[-1].find('Version')==-1:
        NewFolderName = NewFolderName+'_Version1'
    else:
        NewFolderName = NewFolderName+'_Version'+str(int(FolderList[-1].split('_')[-1].split('Version')[-1])+1)
    os.makedirs(os.path.join(Drive,'Compound Reports',NewFolderName))# ,mode=0777)
else:
    os.makedirs(os.path.join(Drive,'Compound Reports',NewFolderName))# ,mode=0777)

FullMethodName = '_'.join(['-'.join([Year[-4:],Month,Day]),NotebookCode,Fraction,Instrument])

print "CREATING MASS HUNTER METHOD..."
# Create and save Mass Hunter method file
OutputDatabaseAndMethod(Drive,midaDatabase,FullMethodName,MassHunterMethod,RTWindow,NewFolderName)

print "CREATING COMPOUND REPORT FOLDER..."
# Create new folder in 'Compound Reports' directory and add notes
peptideListOriginal.to_csv(os.path.join(Drive,'Compound Reports',NewFolderName,PeptideFileName.split('\\')[-1].split('.')[0]+'.txt'), sep=';', mode='a')
if RTWindow != 0.8:
    RTNote = open(os.path.join(Drive,'Compound Reports',NewFolderName,'NOTE - Increased RT window to ' + str(RTWindow) + ' min.txt'),'w')
    RTNote.close()
ScriptNote = open(os.path.join(Drive,'Compound Reports',NewFolderName,'NOTE - Generated results with script ' + Current_Script_Version + ' updated ' + Update_Date +'.txt'),'w')
ScriptNote.close()
shutil.copy(ScriptFileName,os.path.join(Drive,'Compound Reports',NewFolderName,ScriptFileName.split('\\')[-1]))

print "CREATING PARAMETERS TEMPLATE..."
# Create parameters template file
makeParametersFile(SampleSubmissionFileName,os.path.join(Drive,'Compound Reports',NewFolderName),SettingsTable)

# Remove temporary MIDA script files, peptide list and sample submission form
os.remove('INPUT.csv')
os.remove('OUTPUT.csv')
os.remove(PeptideFileName)
#os.remove(SampleSubmissionFileName)

raw_input("SCRIPT COMPLETE!..PRESS ENTER TO EXIT")
