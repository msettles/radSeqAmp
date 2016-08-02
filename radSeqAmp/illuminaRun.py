# illuminaRun.py
#
# illuminaRun.py handles input and output of illumina reads from/to files.
# Input reads are in 4, 2, 1 read sets and output is in 2 and 1 reads setss
import os
import glob
import gzip
import sys
from radSeqAmp import TwoSequenceReadSet
from radSeqAmp import OneSequenceReadSet
from radSeqAmp import misc


class TwoReadIlluminaRun:
    """
    Class to open/close and read a two read illumina sequencing run. Data is expected to be in
    fastq format (possibly gzipped) first processed with radSeqAmp preprocess subroutine
    """
    def __init__(self, read1, read2):
        """
        Initialize a TwoReadIlluminaRun object with expandible paths (with glob) to the two
        sequencing read files. A vector of multiple files per read is allowed.
        """
        self.isOpen = False
        self.mcount = 0
        self.fread1 = []
        self.fread2 = []
        try:
            for fread in read1:
                self.fread1.extend(glob.glob(fread))
                if len(self.fread1) == 0 or not all(os.path.isfile(f) for f in self.fread1):
                    sys.stderr.write('ERROR:[TwoReadIlluminaRun] read1 file(s) not found\n')
                    raise
            if read2 is None:
                for fread in self.fread1:
                    self.fread2.append(misc.infer_read_file_name(fread, "2"))
            else:
                for fread in read2:
                    self.fread2.extend(glob.glob(fread))
                    if len(self.fread2) == 0 or not all(os.path.isfile(f) for f in self.fread2):
                        sys.stderr.write('ERROR:[TwoReadIlluminaRun] read2 file not found\n')
                        raise
            if len(self.fread1) != len(self.fread2):
                sys.stderr.write('ERROR:[TwoReadIlluminaRun] Inconsistent number of files for each read\n')
                raise
        except:
            raise
        # record the number of files per read
        self.numberoffiles = len(self.fread1)

    def open(self):
        """
        Open a OneReadIlluminaRun file set, if file ends in .gz, open will use gzip
        """
        if self.isOpen:
            self.close()
        if self.numberoffiles > 0:
            try:
                read1 = self.fread1.pop()
                if read1.split(".")[-1] == "gz":
                    self.R1 = misc.sp_gzip_read(read1)
                else:
                    self.R1 = open(read1, 'r')
                read2 = self.fread2.pop()
                if read2.split(".")[-1] == "gz":
                    self.R2 = misc.sp_gzip_read(read2)
                else:
                    self.R2 = open(read2, 'r')
            except:
                sys.stderr.write('ERROR:[TwoReadIlluminaRun] cannot open input files\n')
                raise
            self.isOpen = True
            self.numberoffiles -= 1
            return 0
        else:
            return 1

    def close(self):
        """
        Close a TwoReadIlluminaRun file set
        """
        self.R1.close()
        self.R2.close()
        self.isOpen = False

    def count(self):
        """
        Provide the current count of reads read in
        """
        return self.mcount

    def next(self, ncount=1):
        """
        Extract and store the next [count] reads into a TwoSequenceReadSet object.
        If the file object is not open, or if 'next' reaches the end of a file, it will
        attempt to open the file in the list, or gracefully exit
        """
        if not self.isOpen:
            try:
                if self.open() == 1:
                    sys.stderr.write('ERROR:[TwoReadIlluminaRun] ERROR Opening files for reading\n')
                    raise
            except:
                raise
        reads = []
        i = 0
        while i < ncount:
            try:
                # process read_1 (Read1)
                name_1 = self.R1.next().rstrip()  # name
                read_1 = self.R1.next().rstrip()  # read
                self.R1.next()  # '+'
                qual_1 = self.R1.next().rstrip()  # qual
                # process read _2 (Read2)
                name_2 = self.R2.next().rstrip()  # name
                read_2 = self.R2.next().rstrip()  # read
                self.R2.next()  # '+'
                qual_2 = self.R2.next().rstrip()  # qual
                # add it to the stack
                if name_1.split(" ")[0] != name_2.split(" ")[0]:  # check name
                    sys.stderr.write('ERROR:[TwoReadIlluminaRun] Read names do not match each other\n')
                    raise
                reads.append(TwoSequenceReadSet(name_1=name_1, read_1=read_1, qual_1=qual_1, name_2=name_2, read_2=read_2, qual_2=qual_2))
                self.mcount += 1
            except StopIteration:
                if self.numberoffiles > 0:
                    try:
                        if self.open() == 1:
                            sys.stderr.write('ERROR:[TwoReadIlluminaRun] ERROR Opening files for reading\n')
                            raise
                    except:
                        raise
                    continue
                break
            except:
                sys.stderr.write('ERROR:[TwoReadIlluminaRun] Error reading next read\n')
                raise
            i += 1
        return reads


class OneReadIlluminaRun:
    """
    Class to open/close and read a one read illumina sequencing run. Data is expected to be in
    fastq format (possibly gzipped), first processed with radSeqAmp preprocess subroutine and
    then joined using some method like flash.
    """
    def __init__(self, read1):
        self.isOpen = False
        self.mcount = 0
        self.fread1 = []
        try:
            for fread in read1:
                self.fread1.extend(glob.glob(fread))
                if len(self.fread1) == 0 or not all(os.path.isfile(f) for f in self.fread1):
                    sys.stderr.write('ERROR:[OneReadIlluminaRun] read1 file(s) not found\n')
                    raise
        except:
            raise
        # record the number of files per read
        self.numberoffiles = len(self.fread1)

    def open(self):
        """
        Open a OneReadIlluminaRun file set, if file ends in .gz, open will use gzip
        """
        if self.isOpen:
            self.close()
        if self.numberoffiles > 0:
            try:
                read1 = self.fread1.pop()
                if read1.split(".")[-1] == "gz":
                    self.R1 = misc.sp_gzip_read(read1)
                else:
                    self.R1 = open(read1, 'r')
            except:
                sys.stderr.write('ERROR:[OneReadIlluminaRun] cannot open input files\n')
                raise
            self.isOpen = True
            self.numberoffiles -= 1
            return 0
        else:
            return 1

    def close(self):
        """
        Close a OneReadIlluminaRun file set
        """
        self.R1.close()
        self.isOpen = False

    def count(self):
        """
        Provide the current count of reads read in
        """
        return self.mcount

    def next(self, ncount=1):
        if not self.isOpen:
            try:
                if self.open() == 1:
                    sys.stderr.write('ERROR:[OneReadIlluminaRun] ERROR Opening files for reading\n')
                    raise
            except:
                raise
        reads = []
        i = 0
        while i < ncount:
            try:
                # process read_1 (Read1)
                name_1 = self.R1.next().rstrip()  # name
                read_1 = self.R1.next().rstrip()  # read
                self.R1.next()  # '+'
                qual_1 = self.R1.next().rstrip()  # qual
                # add it to the stack
                reads.append(OneSequenceReadSet(name_1=name_1, read_1=read_1, qual_1=qual_1))
                self.mcount += 1
            except StopIteration:
                if self.numberoffiles > 0:
                    try:
                        if self.open() == 1:
                            sys.stderr.write('ERROR:[OneReadIlluminaRun] ERROR Opening files for reading')
                            raise
                    except:
                        raise
                    continue
                break
            except:
                sys.stderr.write('ERROR:[OneReadIlluminaRun] Error reading next read')
                raise
            i += 1
        return reads


class IlluminaTwoReadOutput:
    """
    Given Paired-end reads, output them to a paired files (possibly gzipped)
    """
    def __init__(self, output_prefix, uncompressed):
        """
        Initialize an IlluminaTwoReadOutput object with output_prefix and whether or not
        output should be compressed with gzip [uncompressed True/False]
        """
        self.isOpen = False
        self.output_prefix = output_prefix
        self.uncompressed = uncompressed
        self.R1 = []
        self.R2 = []
        self.mcount = 0
        if self.uncompressed is True:
            if os.path.isfile(self.output_prefix + '_R1.fastq'):
                sys.stderr.write('WARNING:[IlluminaTwoReadOutput] File with prefix: %s exists, DELETING\n' % self.output_prefix)
                try:
                    os.remove(self.output_prefix + '_R1.fastq')
                    os.remove(self.output_prefix + '_R2.fastq')
                except:
                    sys.stderr.write('WARNING:[IlluminaTwoReadOutput] Cannot delete file with prefix: %s\n' % self.output_prefix)
                    raise
        else:
            if os.path.isfile(self.output_prefix + '_R1.fastq.gz'):
                sys.stderr.write('WARNING:[IlluminaTwoReadOutput] File with prefix: %s exists, DELETING\n' % self.output_prefix)
                try:
                    os.remove(self.output_prefix + '_R1.fastq.gz')
                    os.remove(self.output_prefix + '_R2.fastq.gz')
                except:
                    sys.stderr.write('WARNING:[IlluminaTwoReadOutput] Cannot delete file with prefix: %s\n' % self.output_prefix)
                    raise

    def open(self):
        """
        Open the two read files for writing, appending _R1.fastq and _R2.fastq to the output_prefix.
        Create directories as needed.
        """
        if self.isOpen:
            self.close()
        try:
            misc.make_sure_path_exists(os.path.dirname(self.output_prefix))
            if self.uncompressed is True:
                self.R1f = open(self.output_prefix + '_R1.fastq', 'a')
                self.R2f = open(self.output_prefix + '_R2.fastq', 'a')
            else:
                self.R1f = gzip.open(self.output_prefix + '_R1.fastq.gz', 'ab')
                self.R2f = gzip.open(self.output_prefix + '_R2.fastq.gz', 'ab')
        except:
            sys.stderr.write('ERROR:[IlluminaTwoReadOutput] Cannot write reads to file with prefix: %s\n' % self.output_prefix)
            raise
        self.isOpen = True
        return 0

    def close(self):
        """
        Close an IlluminaTwoReadOutput file set
        """
        try:
            self.R1f.close()
            self.R2f.close()
        except:
            raise
        self.isOpen = False

    def count(self):
        """
        Provide the current read count for the file output
        """
        return self.mcount

    def addRead(self, read):
        """
        Add a pair of reads to the output queue
        """
        self.R1.append(read[0])
        self.R2.append(read[1])
        self.mcount += 1

    def writeReads(self):
        """
        Write the paired reads in the queue to the output files
        """
        if (len(self.R1) == 0):
            pass
        else:
            if not self.isOpen:
                try:
                    if self.open() == 1:
                        sys.stderr.write('ERROR:[IlluminaTwoReadOutput] ERROR Opening files for writing\n')
                        raise
                except:
                    raise
            try:
                self.R1f.write('\n'.join(self.R1) + '\n')
                self.R2f.write('\n'.join(self.R2) + '\n')
            except:
                sys.stderr.write('ERROR:[IlluminaTwoReadOutput] Cannot write reads to file with prefix: %s\n' % self.output_prefix)
                raise
            self.R1 = []
            self.R2 = []
            self.close()


class IlluminaOneReadOutput:
    """
    Given single reads, output them to a file (possibly gzipped)
    """
    def __init__(self, output_prefix, uncompressed):
        """
        Initialize an IlluminaOneReadOutput object with output_prefix and whether or not
        output should be compressed with gzip [uncompressed True/False]
        """
        self.isOpen = False
        self.output_prefix = output_prefix
        self.uncompressed = uncompressed
        self.mcount = 0
        self.R1 = []
        if self.uncompressed is True:
            if os.path.isfile(self.output_prefix + '_R1.fastq'):
                sys.stderr.write('WARNING:[IlluminaOneReadOutput] File with prefix: %s exists, DELETING\n' % self.output_prefix)
                try:
                    os.remove(self.output_prefix + '_R1.fastq')
                except:
                    sys.stderr.write('WARNING:[IlluminaOneReadOutput] Cannot delete file with prefix: %s\n' % self.output_prefix)
                    raise
        else:
            if os.path.isfile(self.output_prefix + '_R1.fastq.gz'):
                sys.stderr.write('WARNING:[IlluminaOneReadOutput] File with prefix: %s exists, DELETING\n' % self.output_prefix)
                try:
                    os.remove(self.output_prefix + '_R1.fastq.gz')
                except:
                    sys.stderr.write('WARNING:[IlluminaOneReadOutput] Cannot delete file with prefix: %s\n' % self.output_prefix)
                    raise

    def open(self):
        """
        Open the read file for writing, appending .fastq to the output_prefix.
        Create directories as needed.
        """
        if self.isOpen:
            self.close()
        try:
            misc.make_sure_path_exists(os.path.dirname(self.output_prefix))
            if self.uncompressed is True:
                self.R1f = open(self.output_prefix + '.fastq', 'a')
            else:
                self.R1f = gzip.open(self.output_prefix + '.fastq.gz', 'ab')
        except:
            sys.stderr.write('ERROR:[IlluminaOneReadOutput] Cannot write reads to file with prefix: %s\n' % self.output_prefix)
            raise
        self.isOpen = True
        return 0

    def close(self):
        """
        Close an IlluminaOneReadOutput file set
        """
        try:
            self.R1f.close()
        except:
            raise
        self.isOpen = False

    def count(self):
        """
        Provide the current read count for the file output
        """
        return self.mcount

    def addRead(self, read):
        """
        Add a pair of reads to the output queue
        """
        self.R1.append(read[0])
        self.mcount += 1

    def writeReads(self):
        """
        Write the reads in the queue to the output files
        """
        if (len(self.R1) == 0):
            pass
        else:
            if not self.isOpen:
                try:
                    if self.open() == 1:
                        sys.stderr.write('ERROR:[IlluminaOneReadOutput] ERROR Opening file for writing\n')
                        raise
                except:
                    raise
            try:
                self.R1f.write('\n'.join(self.R1) + '\n')
            except:
                sys.stderr.write('ERROR:[IlluminaOneReadOutput] Cannot write read to file with prefix: %s\n' % self.output_prefix)
                raise
            self.R1 = []
            self.close()


class IlluminaFastaOutput:
    """
    Given single reads, output them to a fasta file (possibly gzipped)
    """
    def __init__(self, output_prefix):
        """
        Initialize an IlluminaFastaOutput object with output_prefix
        """
        self.isOpen = False
        self.output_prefix = output_prefix + '.fasta'
        self.mcount = 0
        self.R1 = []

    def open(self):
        """
        Open the read files for writing, appending .fastq to the output_prefix.
        Create directories as needed.
        """
        if self.isOpen:
            self.close()
        try:
            misc.make_sure_path_exists(os.path.dirname(self.output_prefix))
            self.R1f = open(self.output_prefix, 'a')
        except:
            sys.stderr.write('ERROR:[IlluminaFastaOutput] Cannot write reads to file with prefix: %s\n' % self.output_prefix)
            raise
        self.isOpen = True
        return 0

    def close(self):
        """
        Close an IlluminaFastaOutput file set
        """
        try:
            self.R1f.close()
        except:
            raise
        self.isOpen = False

    def count(self):
        """
        Provide the current read count for the file output
        """
        return self.mcount

    def addRead(self, read):
        """
        Add a pair of reads to the output queue
        """
        self.R1.append(read[0])
        self.mcount += 1

    def writeReads(self):
        """
        Write the reads in the queue to the output files
        """
        if (len(self.R1) == 0):
            pass
        else:
            if not self.isOpen:
                try:
                    if self.open() == 1:
                        sys.stderr.write('ERROR:[IlluminaFastaOutput] ERROR Opening file for writing\n')
                        raise
                except:
                    raise
            try:
                self.R1f.write('\n'.join(self.R1) + '\n')
            except:
                sys.stderr.write('ERROR:[IlluminaFastaOutput] Cannot write read to file with prefix: %s\n' % self.output_prefix)
                raise
            self.R1 = []
            self.close()
