# Copyright (c) 2015 Tobias Neumann, Philipp Rescheneder.
#
# This file is part of Slamdunk.
#
# Slamdunk is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as
# published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# Slamdunk is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from intervaltree import IntervalTree

def bedToIntervallTree(bed):
    utrs = {}

    for utr in BedIterator(bed):

        if (not utr.chromosome in utrs) :
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

    def __next__(self):
        try:
            return self._toBED(self._bedFile.__next__())
        except StopIteration:
            self._bedFile.close()
            raise StopIteration
