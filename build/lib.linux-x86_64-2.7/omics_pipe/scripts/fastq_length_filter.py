#!/usr/bin/python

import sys

def usageMessage():
    print("Usage: fastq_length_filter file.fastq ")
    print("Note: ")
    print("Filters out the reads > 35bp.")
    return()

# Main
    
if len(sys.argv) > 2:
    print("Too many arguments.")
    usageMessage()
    sys.exit()

if len(sys.argv) < 2:
    print("Too few arguments.")
    usageMessage()
    sys.exit()

sequencesFile = sys.argv[1]


f = open( sequencesFile, 'r' )

thisFname = sequencesFile.split('_trimmed.fastq')[0]
thisFname = thisFname + '.fastq'
g = open( thisFname, 'w' )

   
#  barcodeSize = len( barcodes.keys()[0] )
#@HISEQ:89:D0DU5ABXX:7:1101:2384:1935 1:N:0:CGATGT
#12345678901234567890123456789012345678901234567890
#TAGCTTATCAGACTGATGTTGACNNNNNNNNNNNNNNNNNNNNNNNNNNN
#+
#@@@DDDDDB?DFFDHB>BHE93<hhhhhhhhhhhhhhhhhhhhhhhhhhh


done = False

while not(done):
    line = f.readline()
    if line == '':
        done = True
        break
    if line[0] == '@':
        currentHeader = line
        currentSequence = f.readline()
        currentQHeader = f.readline()
        currentQScores = f.readline()
        sequenceSize = len(currentSequence)
#        print 'length read ',sequenceSize
#  final bp is position sequenceSize-2
        finalIndex = sequenceSize-2
#        print ' final=  ',  str(currentSequence[49])
        pos = finalIndex
        found = True
        skip = False 
        while found:
             if pos >= 35 or pos <= 10:
                 skip = True
                 break
             if currentSequence[pos] == 'N':
                   pos = pos -1
             else:
                   found = False 
        if skip:
             continue
        else:
             pos += 1
             currentSequence = currentSequence[0:pos] + '\n'
             currentQScores = currentQScores[0:pos] + '\n'
             g.write( currentHeader )
             g.write( currentSequence )
             g.write( currentQHeader )
             g.write( currentQScores )

f.close()
g.close()
