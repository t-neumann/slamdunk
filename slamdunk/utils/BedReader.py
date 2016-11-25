from intervaltree import IntervalTree

def bedToIntervallTree(bed):
    utrs = {}
        
    for utr in BedIterator(bed):

        if (not utrs.has_key(utr.chromosome)) :
            utrs[utr.chromosome] = IntervalTree()
        
        utrs[utr.chromosome][utr.start:(utr.stop + 1)] = utr.name
        
    return utrs


class BedEntry:

    def __init__(self):
        self.chromosome = ""
        self.start = 0
        self.stop = 0
        self.name = ""
        self.score = "."
        self.strand = "."
        
    def __repr__(self):
        return (self.chromosome + "\t" + str(self.start) + "\t" + str(self.stop) + "\t" + self.name)
    
    def getLength(self):
        return self.stop - self.start
    
    def hasStrand(self):
        return self.strand == "+" or self.strand == "-" 
    
    def hasNonEmptyName(self):
        return self.name != ""

class BedIterator:
    
    def __init__(self, filename):
        self._bedFile = open(filename, "r")
        
    def __iter__(self):
        return self
    
    def _toBED(self, line):
        cols = line.rstrip().split("\t")
        bedEntry = BedEntry()
        bedEntry.chromosome = cols[0]
        bedEntry.start = int(cols[1]) 
        bedEntry.stop = int(cols[2]) 
        bedEntry.name = cols[3]
        
        if (len(cols) > 4) :
            bedEntry.score = cols[4]
        # Add strand info if available
        if (len(cols) > 5) :
            bedEntry.strand = cols[5]
                
        return bedEntry
    
    def next(self):
        try:
            return self._toBED(self._bedFile.next())
        except StopIteration:
            self._bedFile.close()
            raise StopIteration
