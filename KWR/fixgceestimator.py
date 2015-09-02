import argparse
import os

def fixFile(File):    
    inFile = open(File, 'r')
    lines = inFile.readlines()
    outname = File.rstrip('.dat')+'-fixed.dat'
    outFile = open(outname, 'w')
    outFile.write(lines[0])
    line2 = lines[1].rstrip('\n')
    outFile.write(line2+'        diagonal\n')
    data = lines[2]
    data = data.split() 
    low = 0
    high = 10
    while high < len(data):
        rowList = data[low:high]
        rowstring = '%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s%16s' % (rowList[0],rowList[1],rowList[2],rowList[3],
        rowList[4],rowList[5],rowList[6],rowList[7],rowList[8],rowList[9],'0')
        rowstring = rowstring+'\n' 
        outFile.write(rowstring)
        low += 10
        high += 10 
    print 'Fixed file %s \n' % File

def main():
    parser = argparse.ArgumentParser(description='Fixes gce estimator files')
    parser.add_argument('-p', metavar='Input Path',help='Path to estimator files')
    parser.add_argument('-o', metavar='Output Path',help='Path to desired output directory')
    args = parser.parse_args()
    
    #Move to desired directory
    os.chdir(args.p)
    
    # Get all estimator files
    fileList = []
    for fn in os.listdir(os.getcwd()):
        if os.path.isfile(fn):
            if ('estimator' and '.dat' and 'gce') in fn:
                print fn 
                fileList.append(fn)
                
    for File in fileList:  
        fixFile(File) 
    
if __name__=='__main__':
    main()