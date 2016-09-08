#!/usr/bin/env python

#Useful utilities for automate.py

import os
import time
import datetime
import xml.etree.ElementTree as ET 
import hgapi
import re

from sumatra.projects import Project
from sumatra.projects import load_project
from sumatra.parameters import build_parameters


def check_file_exists(input_file, output_file):
    '''Checks if file exists'''
    if not os.path.exists(output_file):
        return True, "Missing file %s for %s" % (output_file, input_file)
    else:
        return False, "File %s exists for %s" % (output_file, input_file) 
   
def spawn_job(jobname, SAMPLE, LOG_PATH, RESULTS_EMAIL, SCHEDULER, walltime, queue, nodes, ppn, memory, script, args_list):
    '''Spawns a job on a cluster (HPC or AWS Star Cluster) using DRMAA'''
    import drmaa
    s = drmaa.Session()
    s.initialize()
    print 'Creating job template for ' + jobname
    jt = s.createJobTemplate()
    print 'Job template created'
    jt.jobName = jobname + "_" + str(SAMPLE)
    print "job name is " + jt.jobName
    jt.outputPath = LOG_PATH
    #print "Error path is" + jt.outputPath
    jt.errorPath = LOG_PATH
    jt.email = RESULTS_EMAIL
    #print "email is" + str(jt.email)
    if SCHEDULER == "PBS":
        jt.hardWallclockTimeLimit = walltime
        jt.softWallClockTimeLimit = walltime
        jt.hardRunDurationLimit = walltime
        jt.softRunDurationLimit = walltime
        jt.nativeSpecification = "-q " + queue + " -l " + "nodes=" + str(nodes) + ":ppn=" + str(ppn) + " -l mem=" + memory  
        print "native specification is" + jt.nativeSpecification
    elif SCHEDULER == "LSF":
        print "LSF currently unsupported. Please contact us to request this feature."
    elif SCHEDULER == "Slurm":
        print "Slurm currently unsupported. Please contact us to request this feature."
    elif SCHEDULER == "SGE":
        jt.nativeSpecification = "-V -cwd -pe local " + str(ppn) + " -l h_vmem=" + re.sub("gb","G",memory)
        print "native specification is " + jt.nativeSpecification
    else:
        print "Scheduler unsupported. Please make sure you have a parameter in your parameter file SCHEDULER with the options PBS, SGE, LSF, StarCluster or Slurm"

    jt.remoteCommand = os.getcwd() + script 
    #print "remote command is" + jt.remoteCommand   #THIS PRINTS THEN HANGS
    jt.args = args_list 
    print "ArgsList is " + str(jt.args)
    jt.joinFiles = True
    jobid = s.runJob(jt)
    now = datetime.datetime.now()
    date = now.strftime("%Y-%m-%d %H:%M")
    print "Date/Time: " + date 
    print "Job has been submitted with id" + jobid + " at Date/Time: " + date
    retval = s.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)
    now = datetime.datetime.now()
    date = now.strftime("%Y-%m-%d %H:%M")
    print "Job: " + str(retval.jobId) + ' finished with status: ' + str(retval.hasExited) + ' and exit status: ' + str(retval.exitStatus) + " at Date/Time: " + date
    print "Date/Time: " + date 
    print 'Cleaning up'
    s.deleteJobTemplate(jt)
    s.exit()
    return 

def job_status(jobname, resultspath, SAMPLE, outputfilename, FLAG_PATH):
    '''Checks to see if a job has successfully finished by checking if the specified output file exists.'''
    stage = jobname
    try:
        size = os.path.getsize(resultspath + "/" + outputfilename) 
    except OSError:
        print("Looking for file: " + resultspath + "/" + outputfilename)
        print("%s failed to produce any output files" % stage) 
    if 'size' in locals():
        if size == 0:
            print "Job Failed!"
            print('%s produced an empty output file' % stage)
        else:
            print("%s Finished and Successfully produced an output file of size %s" % (stage,size))
            flag_file = "%s/%s_%s_completed.flag" % (FLAG_PATH, stage, SAMPLE)
            open(flag_file, 'w').close()
    else:
        print("%s failed to produce any output files" % stage)   
    return

def job_status_nfsc(jobname, resultspath, SAMPLE, outputfilename, FLAG_PATH):
    '''Checks to see if a job as succesfully finished by checking if the specified output file exits.'''
    stage = jobname
    try:
        size = os.path.getsize(resultspath + "/" + outputfilename)
    except OSError:
        print("Looking for file: " + resultspath + "/" + outputfilename)
        print("%s failed to produce any output files" % stage)
    if 'size' in locals():
        print("%s Finished and Successfully produced an output file of size %s" % (stage,size))
        flag_file = "%s/%s_%s_completed.flag" % (FLAG_PATH, stage, SAMPLE)
        open(flag_file, 'w').close()
    else:
        print("%s failed to produce any output files" % stage)
    return
    


def decompress(file_directory, compression_type):
    '''Unzips (bzip or gzip) files'''
    if compression_type == "bzip":
        os.chdir(file_directory)
        os.system("bunzip2 *.bz2")
        print "Bzip files successfully uncompressed"
    elif compression_type == "gzip":
        os.chdir(file_directory)
        os.system("gunzip *.gz")
        print "Gzip files successfully uncompressed"
    else:
        print "Working with uncompressed files"
    return

def compress(file_directory, compression_type):
    '''Compresses (bzip or gzip) files'''
    if compression_type == "bzip":
        os.chdir(file_directory)
        os.system("bzip2 *.fastq")
        print "Bzip files successfully compressed"
    elif compression_type == "gzip":
        os.chdir(file_directory)
        os.system("gzip *.fastq")
        print "Gzip files successfully compressed"
    else:
        print "Working with uncompressed files" 
    return


class Bunch(object):
    '''Bunches parameters into a dictionary'''
    def __init__(self, adict):
        self.__dict__.update(adict)
        
        
def check_create_dir(directory):
    '''Creates a directory if it does not exist'''
    if not os.path.exists(directory):
        os.makedirs(directory)
    return


def parse_xml(file):
    '''Parses XML file to extract sample names from TCGA XML manifest'''
    tree = ET.parse(file)
    root = tree.getroot()
    for child in root:
        tree = ET.ElementTree(child)
        for id in tree.findall('analysis_id'):
            print id.text    
            name = "/gpfs/home/kfisch/" + id.text + ".xml"
    tree.write(name)     
    return
  

def get_TCGA_ID(file):
    '''Gets TCGA ID from TCGA XML manifest and creates a sample list'''
    tree = ET.parse(file)
    root = tree.getroot()
    sample_list = []
    for id in root.findall('Result'):
        analysis_id = id.find('analysis_id').text
        analysis_id_string = "TCGA_" + str(analysis_id)
        sample_list.append(analysis_id_string)
    print sample_list    
    return sample_list   


def make_params(step, sample_list, flag_path):
    '''Creates parameter input lists for each step in pipeline'''
    vars()['inputList_' + step] = []
    for sample in sample_list:
        vars()['inputList_' + step].append([sample, "%s/%s_%s_completed.flag" % (flag_path, step, sample)])
        #print vars()['inputList_' + step] 
    return vars()['inputList_' + step] 


def get_samples_from_txt_file(file):
    '''Creates sample list from text file of samples'''
    global SAMPLE_LIST
    sample_file = open(file, 'r')
    reader = csv.reader(sample_file)
    sample_list = [row for row in reader]
    return SAMPLE_LIST

def sumatra_start(repository, sumatra_db_path, results_path, working_dir, hg_username, sumatra_run_name, parameters):
    '''Clones the Omics Pipe repository from Bitbucket, creates a Sumatra project, and creates a Sumatra record for the current run'''
    print "sumatra_db_path is " + sumatra_db_path
    print type(sumatra_db_path)
    check_create_dir(sumatra_db_path)
    os.chdir(sumatra_db_path)
    repo1 = hgapi.Repo(repository)
    repo_path = sumatra_db_path +"/omics_pipe"
    repo= {"url":repo_path, 
           "type":"sumatra.versioncontrol._mercurial.MercurialRepository",
           "upstream":repository}
    executable= {"path":"",
                 "version": "",
                 "type":"sumatra.programs.PythonExecutable",
                 "options":"",
                 "name": "Python"}
    sumatra_launch_mode = {"working_directory": working_dir, "type": "sumatra.launch.SerialLaunchMode"}
    data_store1 = {"root":results_path, "type": "sumatra.datastore.filesystem.FileSystemDataStore"}
    database_path = sumatra_db_path + "/records/recordstore.db"
    record_store1 = {"db_file": database_path, "type": "sumatra.recordstore.django_store.DjangoRecordStore"}
    input_datastore1 = {"root": results_path, "type": "sumatra.datastore.filesystem.FileSystemDataStore"}
    while True:
        try:
            repo1.hg_clone(url = repository, path=repo_path)
            with open(repo_path + "/.hg/hgrc", "a") as myfile:
                myfile.write("[ui]\nusername= " + hg_username)         
            print "Omics pipe repository cloned to : " + repo_path
            break
        except hgapi.hgapi.HgException:
            print "Omics pipe repository already exists."
            break
    while True:
        try:
            Project(sumatra_run_name, default_repository=repo, default_executable=executable, 
                    default_launch_mode = sumatra_launch_mode, on_changed='store-diff',
                    data_store=data_store1, record_store=record_store1, input_datastore=input_datastore1)            
            print "Sumatra project created: " + sumatra_run_name + " in directory: " + sumatra_db_path
            break
        except Exception:
            print "Sumatra project already exists, loading project: " + sumatra_run_name
            break
    project = load_project(path=sumatra_db_path)
    print project
    sumatra_params = build_parameters(parameters)
    print sumatra_params
    os.chdir(repo_path)
    repo_main = "omics_pipe/main.py"
    record = project.new_record(parameters=sumatra_params, main_file=repo_main)
    print record
    return record,project

def sumatra_end(start_time, record, project):
    '''Saves the Sumatra project to the database'''
    #file1 = open("/gpfs/home/kfisch/test/test.txt", "w")
    #file1.write("test")
    #file1.close()
    record.duration = time.time() - start_time
    print record.duration
    record.output_data = record.datastore.find_new_data(record.timestamp)
    print record.output_data
    project.add_record(record)
    print project
    project.save()
    return
