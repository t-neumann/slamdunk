class BedEntry:
    
    chromosome = ""
    start = 0
    stop = 0
    name = ""
    strand = "."
    
    def __repr__(self):
        return (self.chromosome + "\t" + str(self.start) + "\t" + str(self.stop) + "\t" + self.name)
    
    def getLength(self):
        return self.stop - self.start
    
    def hasStrand(self):
        return self.strand != "."
    

class BedIterator:
    
    _bedFile = None
    
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
        
        # Add strand info if available
        if (len(cols) > 4) :
            bedEntry.strand = cols[5]
                
        return bedEntry
    
    def next(self):
        try:
            return self._toBED(self._bedFile.next())
        except StopIteration:
            self._bedFile.close()
            raise StopIteration
