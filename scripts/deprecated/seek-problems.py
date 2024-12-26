#!/usr/bin/python2



#=============================================================
def search_for_problems(descendants, progenitors, params):
#=============================================================
    """
    Seek for descendants with main progenitor = 0
    and other progenitors that are marked as
    merged into said descendant
    """

    import numpy as np


    for out in range(len(descendants)):
        
        print "checking in output", params.noutput-out

        descs = descendants[out]
        if (len(descs)>0):
            is_main = np.ones(len(descs), dtype=bool)
            is_main[:] = True
            
            for i in range(len(descs)):
                if (descs[i]<0):
                    descs[i] *= -1
                    is_main[i] = False
                
            # sort lists
            sort_ind=list(range(len(descs)))
            descs, sort_ind = zip(*sorted(zip(descs, sort_ind)))

            descs = np.array(descs)
            sort_ind = np.array(sort_ind)
            progs = np.array(progenitors[out])
            progs = progs[sort_ind]
            is_main = is_main[sort_ind]

            uniqs, counts = np.unique(progs, return_counts=True)

            for i in range(counts.shape[0]):
                if counts[i] > 1 and uniqs[i]!=0:
                    print "problem 3: multiple instances of progenitor", uniqs[i], "found"

            for i in range(len(progs)):
                if progs[i]!=0 and descs[i]==0:
                    print "problem 4: prog with descendant 0", "D:", descs[i], "P:", progs[i]


            start = 0
            while start < len(descs)-1:
                # find start and stop for any descendant
                if descs[start]==descs[start+1]:
                    end = start+1
                    while end < len(descs)-1:
                        if descs[start] == descs[end+1]:
                            end += 1
                        else:
                            break


                    thisprogs = progs[start:end+1]
                    thismains = is_main[start:end+1]

                    counter = 0
                    for i in range(len(thisprogs)):
                        if is_main[i] and thisprogs[i]==0:
                            print "found problem 1: main_prog=0 with mergers in output", params.noutput-out, "D:", descs[start], "Progs:", thisprogs
                            counter += 1

                    if counter > 1:
                        print "found problem 2: multiple mains in output", params.noutput-out, "D:", descs[start], "Progs:", thisprogs

                    
                   

                start += 1


    return





#==============
def main():
#==============
    """
    Main function. Calls all the rest.
    """

    import treeplot_module as tm


    #-----------------------
    # Set up
    #-----------------------

    params = tm.global_params()
    #  params.read_cmdlineargs()
    params.get_output_info()
    params.set_outputfilename()


    descendants, progenitors, progenitor_outputnrs, outputnrs, t = tm.read_mergertree_data(params)


    search_for_problems(descendants, progenitors, params)



#=================================
if __name__ == "__main__":
#=================================
    main()
