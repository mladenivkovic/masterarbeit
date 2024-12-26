#!/usr/bin/env python3


#====================================================
# print z in every output
# execute in parent directory of output_XXXXX dirs
#====================================================


import numpy as np




#====================
def main():
#====================

    from os import getcwd
    from os import listdir

    workdir = getcwd()
    filelist = listdir(workdir)

    outputlist = []
    outdirtemp = 'output_'
    for filename in filelist:
        if filename.startswith(outdirtemp):
            outputlist.append(filename)

    if len(outputlist)<1:
        print("I didn't find any output_XXXXX directories in current working directory.")
        print("Are you in the correct workdir?")
        quit()

    outputlist.sort()
    noutput = len(outputlist)

    outputnrs = np.array(range(1, noutput+1))
    a_exp = np.zeros(noutput)

    for i,out in enumerate(outputlist):
        fn = out+'/info_'+out[-5:]+'.txt'
        print("Reading", fn)
        try:
            infofile = open(fn)
            for j in range(9):
                infofile.readline() # skip first 9 lines

            # get expansion factor
            aline = infofile.readline()
            astring, equal, aval = aline.partition("=")
            afloat = float(aval)
            a_exp[i] = afloat

        except IOError: # If file doesn't exist
            print("Didn't find any info data in ", out)
            break

    redshift = 1.0/a_exp - 1


    print(a_exp)

    print("Output", "\t", "redshift")
    for i in range(noutput):
        print(outputnrs[i], "\t", redshift[i])



    return


#==================================
if __name__ == "__main__":
#==================================

    main()
