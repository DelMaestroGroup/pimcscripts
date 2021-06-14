# ==============================================================================
# Script for uploading files to mediaWiki site.
#
# Author:  Bohdan Kulchytskyy
# Last Modified: 21-AUG-2013 by Max Graves
# ==============================================================================

import wikitools
#import parser
import argparse,os,shutil,glob

# ==============================================================================
def parseCMD():
    parser = argparse.ArgumentParser(description='runs saveMCbatch')
    parser.add_argument('fileNames', help='Files to be uploaded',
            nargs='+')
    parser.add_argument('-F', '--Force', action='store_true', dest='Force',
            default=False,
            help='Do we want to write over a file with the same name?')
       
    return parser.parse_args()
 
# ==============================================================================
def upload_file(file_name, force=True):
    '''
    Upload a file from a HD to a mediawiki.
    Set force to False to avoid overriding an existing file.
    '''
 
    wiki_url      = 'http://group.delmaestro.org/api.php'
    wiki_username = ""
    wiki_password = ""
 
    if not(file_name.split('.')[-1] in  ['png', 'gif', 'jpg', 'jpeg', 'pdf', 
        'tiff', 'bmp', 'ps', 'odt', 'ods', 'odp', 'odg', 'txt', 'zip',
        'pbs', 'mpg', 'mp4', 'avi', 'svg']):
       print 'Unsupported file type'
       return 1

    # Connect and login
    try:     wiki = wikitools.wiki.Wiki(wiki_url)
    except:  print "Cannot connect with wiki. Check the URL"
    try:     wiki.login(username=wiki_username,password=wiki_password)
    except:  print "Invalid Username or Password"
 
    # File Upload       
    file_object = open(file_name,"r")
    fwiki = wikitools.wikifile.File(wiki=wiki, title=file_name)
    if force : 
        rmessage = fwiki.upload(fileobj=file_object,comment="Auto upload",
                ignorewarnings=True)
    else:      
        rmessage = fwiki.upload(fileobj=file_object,comment="Auto upload", 
                ignorewarnings=False)
 
    # Checking    
    result = rmessage['upload']['result']
    if result == 'Success':
        print 'Upload: ' + result
    else: 
        print rmessage['upload']['warnings']
        return 1
    return 0

# =============================================================================
def main():

    args = parseCMD()
    fileNames = args.fileNames
    Force = args.Force

    print Force
    #tempDirName = 'kitties'
    #if not os.path.exists(tempDirName):
    #    os.makedirs(tempDirName)

    #tempDirPath = os.getcwd()+'/'+tempDirName

    #for f in fileNames:
    #    shutil.copy2(f, tempDirPath)

    #os.chdir(tempDirPath)
    #print os.getcwd()

    #files = glob.glob('*')
    #print files

    scriptDir = os.getcwd()

    for f in fileNames:
        os.chdir(os.path.split(f)[0])
        upload_file(os.path.split(f)[1], Force)
        os.chdir(scriptDir)

    #shutil.rmtree(tempDirName)

# =============================================================================
if __name__=='__main__':
    main()
