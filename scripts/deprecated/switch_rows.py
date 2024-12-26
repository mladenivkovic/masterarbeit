#!/usr/bin/python3



fn = 'ntracers_global_stats.tex'

f = open(fn, 'r')
lines = f.readlines()
f.close()


f = open(fn, 'w')
for l in lines:
    if "&" in l:
        split = l.split("&")
        
        line = split[0]+"&"+split[1]+"&"+split[3]+"&"+split[2]+"&"+split[4]+"\n"

    else:
        line = l

    f.write(line)
