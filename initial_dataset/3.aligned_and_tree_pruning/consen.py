import sys
filename=sys.argv[1]
outname=sys.argv[2]
f = open(filename, 'r')
data=f.read()
d=data.replace('\n\n', '@').replace('\n', '.')
ar1 = []
splat = d.split("@")
finalconsen=''
finalnuc=''
##most_frequent function learned and imported from https://www.geeksforgeeks.org/python-find-most-frequent-element-in-a-list/
def most_frequent(List):
    counter = 0
    num = List[0]
    l=len(List)
    for i in List:
        curr_frequency = List.count(i)
        if(curr_frequency> counter):
            counter = curr_frequency
            num = i
    if counter/l > 0.25:
        return(num)
    else:
        return('-')
for number, paragraph in enumerate(splat, 1):
        ar1 = [paragraph]
        print(number)
        print(ar1)
        newar=ar1[0].split('.')
        if '' in newar:
                newar.remove('')
        newar=[ x for x in newar if ' ' not in x ]
        leng=len(newar[0])
        numele=len(newar)
        for pos in range(0,leng):
                consen=[]
                for ele in range(0,numele):
                        obj=newar[ele][pos]
                        consen.append(obj)
                consen=[ x for x in consen if '#' not in x ]
                print(consen)
                print(finalnuc)
                if not consen:
                        consen.append('-')
                finalnuc=most_frequent(consen)
                finalconsen = finalconsen + finalnuc
print(finalconsen)
print(finalconsen, file=open(str(outname)+"_consensus.txt", "w"))
