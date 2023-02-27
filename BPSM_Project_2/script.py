#!/usr/bin/python3
import os, subprocess, re
import pandas as pd

#Because there will a lot of prompts in this code I wanted to highlight to the user when the program is waiting for his input with a certain colour:
black   = "\033[0;30m"
red     = "\033[0;31m" #This will highlight errors
green   = "\033[0;32m" #This will highlight successfully exectued tasks
yellow  = "\033[0;33m" #This will highlight prompts (asking for input)
cyan = "\033[0;36m" #This will highlight instructions (advice on how to use the script)
blue = "\033[0;34m"
white   = "\033[0;37m"


print("")
print("           Welcome to the script!")
print("Coded by a student from the University of Edinburgh")
print("___________________ BB226659 ______________________" )
print("")

#I created different functions for input and esearch in order to have a good error trap control later and be able to redo searches and/or inputs
def input_function():
  #First I need to ask the user for 2 variables (protein family and the taxonomic group)
  p_family=input(yellow+"Please input the protein family: "+white)# or "pyruvate dehydrogenase" #Ask for the input that will be saved in a variable
  while p_family == "": #While the variable is empty
    print(red+"The input cannot be empty!")
    p_family=input(yellow+"Please input the protein family: "+white) #Ask again for the input and put in in the variable
      
  t_group=input(yellow+"Please input the taxonomic group: "+white)# or "ascomycete fungi"
  while t_group == "":
    print(red+"The input cannot be empty!"+white)
    t_group=input(yellow+"Please input the taxonomic group: "+white) 
  
  
  
  #Create a linux command which will be run trhough os module
  cmd="esearch -db protein -query "+p_family+" NOT PARTIAL | efilter -query "+t_group+""
  return cmd




def search_function(command): #get count function
  #assert command==""  
  #search_o=subprocess.check_output(cmd).decode("utf-8")
  search_o=os.popen(command).read() #Save the output of the search into a variable

  #Import a xml module "cElementTree", to easier extract a count result
  # The following code was modified by me, but inspired by the following: https://stackoverflow.com/questions/7691514/extracting-text-from-xml-using-python
  from xml.etree import cElementTree as readXML
  root = readXML.fromstring(search_o)
  count=root.find('Count').text
  print("Your input resulted in", count, "search results")
  return count


#assert search_function()

#I also decided to create a function, that will nest all the steps of chossing the correct set of sequences desired by the user (input_function and search_function).
#This allows me to easily repeat the search process if the user is not satisfied with the outcome of his/her search.
def step1_search():
  

  #Get the input and prepare command which will be run on linux
  cmd=input_function() 


  #Run the search and check the outcome, if the number of results is acceptable
  #esearch can either fail - SOLVED (in except AttributeError:)
  #Or it can return less than two results - SOLVED
  #Or it can return too many hits- SOLVED
  while True: #Retry in case of exception
    try: #Try the search
      count = search_function(cmd)
      if int(count)<1000 and int(count)>2: #Test the number of counts
        break #Exit the loop is succesfull and no. of hits acceptable
      elif int(count)<2:
        print(red+"Your input resulted in an empty search. Change your input!"+white)
        cmd=input_function()
      else: 
        print(red+"Your input resulted in a search with too many hits. Narrow down your search"+white)
        cmd=input_function()
    except AttributeError: #In case if the search fails, ask the user for a different input
      print(red+"Your input resulted in no hits :("+white)
      print("Change the input and run another search :)")
      cmd=input_function()

  
  print("Downloading the sequences, please wait :)")
  cmd_fetch="{} | efetch -format fasta > sequences.fa".format(cmd) #prepare a command which will be run in a shell and save output in a file
  os.popen(cmd_fetch).read()

  if os.path.exists("sequences.fa") == False or os.path.getsize('sequences.fa') == 0: #Check if the file was created and if it has any contents, if not ask for the input again.
    print(red+"There seems to be an error, please try again :("+white)
    step1_search()  
  



  #The next step is to create a list with all the species that came up in the search in order to allow the user to change the input in order to allow a more/
  #useful processing
  with open('sequences.fa') as sequences : #read the file
    species = []
    for line in sequences: #read lines of the file
      
      if line[0]==">" and "[" in line: #If the line starts with a ">" and has a "[" character
        #print(line)
        start, end = "[", "]"
        species_name=line[line.find(start)+len(start):line.find(end)] #Get the name of the species, which written in a fasta file "[species]'
        
        species.append(species_name) #Append the species names into a list

  #uniq_species=set(species)#Testing the outputs
  #print(uniq_species)
  no_uniq_species=len(set(species)) # This line counts how many unique species the dataset contains

  print("--------------------------------------------------------------------------------")
  print(green, count,"sequences downloaded from",no_uniq_species,"different species"+white)
  print("--------------------------------------------------------------------------------")

  def ask_2():
    #assert len(ask)==0 Checking the function
    ask=input(yellow+"Would you like to:\n1. Continue with the current dataset?\n2. Change the input and acquire a different dataset?\nChoose 1 or 2: "+white)
    if ask=="1":
      return
    elif ask=="2":
      print("Removing the old dataset!")
      os.remove("sequences.fa")
      step1_search()  
    else:
      print(red+"Please choose one of the options (1 or 2)"+white)
      ask_2()

  
  ask_2()

  
    
      










#The script code starts executing here:
#-------------------------------------------------------
#Step 1
start = os.getcwd()
step1_search()







#--------------------------------------------------------
#STEP 2
#CLUSTALO
print("Aligning the sequences")
os.system('clustalo --force --verbose --threads 32 -i sequences.fa -o aligned_seqs.fa') #Maybe add -maxnumseq 1000 (set a limit to how many sequences should be aligned)
# Added --verbose, because it can take a long time, so it is nice to show what is happening behind the scenes
# --force to overwrite the old fasta file if script was run more than once
# --threads to multithread the procces, which hopefully speeds it up :D




#I can also create a distance matrix for which I prepared the code, but I decided to not use it, as I am struggling to understand it's value in this case
# and it just takes too long if the number of sequences is high
'''
os.system('clustalo --force --verbose --threads 32 -i sequences.fa -o aligned_seqs.fa --percent-id --distmat-out=distance_matrix.txt --full')
'''






#EMBOSS PLOTCON -I have put it in a function, to allow user to change the setting on which the plot will be calculated and repeat the process
def plotting_plotcon(): 

  win_size=input(yellow+"Choose the window size "+cyan+"(any +ve integer number)"+yellow+", which will be used to calculate the similarity: "+white)# or "5"
  
  #assert int(win_size)>0 #The number has to larger than 0
  if win_size.isnumeric() and int(win_size)>0: #Checking if the conditions are met (positive integer)
    os.system('plotcon -sequence aligned_seqs.fa -winsize ' +win_size+ ' -graph pdf') 
  elif win_size=="":
    print(red+"The input has to be an integer number! It can't be empty"+white)
    plotting_plotcon()
    return
  else:
    print(red+"The input has to be an integer number!"+white) #If the user types in text or boolean, the default value of 4 will be used
    plotting_plotcon()
    return
  

  #os.wait()
  # Because I have already limited the number of sequences that the user can process, there is no need to limit it here
  # I chosed pdf as my desired output as the ps was returning the graph vertically, while it is nicer to read in a horizontal way

  #I want to open the plotcon output graph for a user to see
  """ 
  I tried different ways to open the images, but all of them were very slow...

  1.
  #os.system('gio open plotcon.1.png') #Takes ages

  2.
  from matplotlib import pyplot as plt
  from matplotlib import image as mpimg
  image = mpimg.imread("plotcon.1.png")
  plt.imshow(image)
  plt.show()
  """
  #The best method, but it returned a lot of text into the terminal, which looked unprofessional and confusing, it also returned an error
  # if a window containing the image was closed before code continued. Therefore I had to make the errors quiet.
  print(cyan+'Close the plot window and then press enter in order to continue!'+white)
  subprocess.run('gs -q plotcon.pdf', shell=True, stderr=subprocess.DEVNULL)

  def ask_wSize():
    ask2=input(yellow+"1. Are you happy with the plot and would like to continue?\n2. Would you like to change the Window Size and create a new plot?\nChoose (1 or 2): "+white)
    
    if ask2 == "1":
      return
    elif ask2 == "2":
      plotting_plotcon()

    #else:
    #  print(red+"Please choose one of the options (1 or 2)"+white)
    #  ask_wSize() #Ask for input untill either option 1 or 2 is given
    
  ask_wSize()


print("")
print("_________________________________________________________________________________________________")
print(cyan+"We are now going to plot the level of conservation between the sequences from out set"+white)
input(yellow+"Press Enter to continue..."+white)
print("")


plotting_plotcon()


#Showing off with my error traping abilities (various ways of doing it) :D
#try:
#  plotting_plotcon()
#except AssertionError: 
#  print("The value of the window size has to be larger than 0")
#  plotting_plotcon()





#--------------------------------------------------------
#STEP 3
# I need to split my multi-fasta file into separate fasta files, I will need that to run each sequence against prosite databse searching for motifs
print("")
print("___________________________________________________________________________________________________________")
print(cyan+"We will now determine whether any known motifs (domains) are associated with our subset of sequences"+white)
input(yellow+"Press Enter to continue..."+white)
print("")


print("Splitting the sequences into separate files...")
try:
  os.mkdir("Individual") #Create a new folder for the individual protein fasta sequences
except FileExistsError: #If the script was run for the second time, ingore the error as the directory will already be created
  pass




os.system('rm -f ./Individual/*.fa') #Delete the .fa files if the script is being run for the second time on a different (smaller dataset)

#pullseq or seqretsplit Emboss programe could be used instead, but I wanted to represent my python abilities

c=0 #Counter will be used in naming files
with open('sequences.fa') as sequences : 
    for line in sequences:
        if line[0]==">":
            c += 1
            
            f = open("./Individual/indseq"+str(c)+".fa", "w")
            f.write(line)
        else:
            f.write(line)

os.system('rm -f ./Individual/*.table') #Delete the .table files if the script is being run for the second time


print("Individual sequences are now stored in a new direcotry 'Individual'")


import warnings 
warnings.simplefilter(action='ignore', category=FutureWarning) #New versions of pandas do not like the append command, but I tried using concat, but it was not doing what I wanted


print("The sequences are now being scanned for the known motifs. Please wait :)")
#Run patmatmotifs on every single individual file
for file in os.listdir("Individual"):
    if file.endswith(".fa"):
      #os.system("patmatmotifs -sequence './Individual/"+file+"' -outfile './Individual/" +file+".table' -rformat excel -raccshow2 Y -verbose N") #Does not silence the output, so i did subprocess in order to be able to avoid the function printing unneccesary stuff
      subprocess.run("patmatmotifs -sequence './Individual/"+file+"' -outfile './Individual/" +file+".table' -rformat excel -raccshow2 Y -verbose N", shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
#Prepare a patmatmotif output table for the user
co=0
for f in os.listdir("Individual"):# loop through all the files in the folder...
    if f.endswith(".table"): #But only process the files that end with .table
        co += 1 #If it is a first file, create a dataframe
        if co==1:
            prosite=pd.read_csv('./Individual/'+f+'', sep='\t', header=0)
            
        else: #if the dataframe already exists, just append a new file to an already existing dataframe
            new_row= pd.read_csv('./Individual/'+f+'', sep='\t')
            #pd.concat(prosite, new_row)
            prosite= prosite.append(new_row)

prosite.to_csv('prosite.tsv', index=False, sep ='\t')
print(prosite) #The code returns a table, with the name of a sequence and various information about the motif (start, end, score, stand, motif)


print(green+"This table has been saved to your current directory as prosite.tsv"+white)
display1=input(yellow+" Would you like to open the table in the Vim editor for a better view? \n1. Yes, display the full table \n2. No, just continue  \nChoose (1 or 2): "+white)
if display1=="y" or display1=="1" or display1=="Y" or display1=="Yes" or display1=="yes":
  input(cyan+"If you want to exit the table and return to the program, use "+white+"[Shift+Z+Z]"+cyan+" keybord combination\nNow press Enter to open the table"+white)
  subprocess.run("vim -R prosite.tsv", shell=True)  #This will allow the user to see all the results and scroll if they do not fit on the page Read-Only









#--------------------------------------------------------
#STEP 4 - Wildcard option
# I will try to predict the secondary structure of all the sequences fetched by the user. Because protein structure is more conserved than the sequence,
# secondary structure can be used to improve sequence alignment quality

#Garnier Osguthorpe Robson algorithm
#I will use the GOT (garnier) method, despite it not being the most accurare, but it can be used on most of everyday computers
print("")
print("")
print("_________________________________________________________________________________________________________________________")
print(cyan+"We are now going to predict the secondary structure of the proteins using Garnier Osguthorpe Robson algorithm"+white)
input(yellow+"Press Enter to continue..."+white)
print("")

print("The secondary structure prediciton is being calculated right now, please wait :)")
def secondary_structure():
  cols=['File_name', 'Sequence_name', 'Total_res_perc']
  sname=[]
  bline=[]
  filenames=[]
  
  for f in os.listdir("Individual"):
    if f.endswith(".fa"):
      subprocess.run("garnier -sequence './Individual/"+f+"' -outfile './Individual/"+f+"_o' ", shell=True, stderr=subprocess.DEVNULL, stdout=subprocess.DEVNULL)
      try:
        with open(start+"/Individual/"+f+"_o", "r+") as fp:
            for line in fp:
                if line.startswith("# Sequence:"):  
                    sname.append(line[11:24])
                    filenames.append(f)
                elif line.startswith("#         percent:"):
                    bline.append(line[19:50])
      except FileNotFoundError:
        continue        
      #datafr=pd.DataFrame([filenames, sname, bline], columns=['File_name', 'Sequence_name', 'Total_res_perc'])

  secondary_str = pd.DataFrame(
      {'File_name': filenames,
      'Sequence_name': sname,
      'Total_res_perc': bline
      })


  secondary_str.sort_values(by=['File_name'], inplace=True, ascending=True)
  print(secondary_str)
  return(secondary_str)
  

  
table=secondary_structure()

table.to_csv('secondary_structure.tsv', index=False, sep ='\t')


print("")
print(green+"This table has been saved to your current directory as secondary_structure.tsv"+white)
display2=input(yellow+" Would you like to open the table now in the Vim editor for a better view? \n1. Yes, display the full table \n2. No, just continue  \nChoose (1 or 2): "+white)
if display2=="1" or display2=="y" or display2=="Yes" or display2=="yes" or display2=="Y":
  input(cyan+"If you want to exit the table and return to the program, use "+white+"[Shift+Z+Z]"+cyan+" keybord combination\nNow press Enter to open the table"+white)
  subprocess.run("vim -R secondary_structure.tsv", shell=True) #Read only




















print("")
print("\033[0;35m")
print("\033[3m")
print("Thank you for using this program <3"+white) #There is no point asking the user if they want to run it again, as it would be just easier for them to call the script again



