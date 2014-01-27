# THIS ONLY WORKS FOR SINGLE ESTIMATORS RIGHT NOW!!

import pylab as pl
import MTG_jkTools as jk

def main():

    args = jk.parseCMD()

    files = args.fileNames

    for f in files:

        headers = jk.getHeadersFromFile(f)

        AVG = pl.array([])
        STD = pl.array([])

        for n in range(int(len(headers))):

            avgs,stds,bins = pl.genfromtxt(f, 
                    usecols=(0+3*n,1+3*n,2+3*n), 
                    unpack=True, delimiter=',')

            # get rid of any items which are non numbers
            bins = bins[pl.logical_not(pl.isnan(bins))]
            stds = stds[pl.logical_not(pl.isnan(stds))]
            avgs = avgs[pl.logical_not(pl.isnan(avgs))]
            
            weights = bins/pl.sum(bins)

            avgs *= weights
            stds *= weights

            avg = pl.sum(avgs)
            stdErr = pl.sum(stds)

            print avg,' +/- ',stdErr

            AVG = pl.append(AVG, avg)
            STD = pl.append(STD, stdErr)

        pl.errorbar(headers,AVG,STD,fmt='o',color='Pink')
        
        pl.show()
    


if __name__=='__main__':
    main()
