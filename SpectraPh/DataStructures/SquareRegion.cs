using System;
using System.Linq;
using System.Text.RegularExpressions;
using SpectraPh.Actions;

namespace SpectraPh.DataStructures
{
    public class SquareRegion
    {
        public uint leftStart, leftEnd, rightStart, rightEnd;
        public int leftSNPidx = -1, rightSNPidx = -1; //This is the nearest bin, by SNPs, which the end falls in. Not necessarily overlapping!
        public byte leftVariant, rightVariant; //0 - new. 1/2 - known variants. 3 - no variant. 4 - unknown variant
        public short lchr, rchr;
        public object locker = new object();
        private bool initialized = false;

        public double LeftMidpoint()
        {
            return leftStart + ((leftEnd - leftStart) / 2);
        }
        public double RightMidpoint()
        {
            return rightStart + ((rightEnd - rightStart) / 2);
        }

        public void SetOrderOfVariants(short chrl, short chrr, uint leftstart, uint rightstart)
        {
            if (initialized)
            {
                if (lchr != chrl || rchr != chrr || leftStart != leftstart || rightStart != rightstart)
                    throw new Exception("bug in my code");
            }
            lchr = chrl;
            rchr = chrr;
            leftStart = leftstart;
            rightStart = rightstart;
            initialized = true;

        }

        public void SetLeftVariant(string seq, string cigar, string phred=null)
        {
            leftEnd = leftStart + ParseCIGARlen(cigar);
            if (!IsSpanningSnp(lchr, leftStart, leftEnd, out leftSNPidx))
            {
                leftVariant = 3;
                return;
            }
            var SNPpos = Program.SNPmap[lchr - 1].ElementAt(leftSNPidx);
            if (!string.IsNullOrEmpty(phred))
            {
                var qualityOfVariant = CharToPhred(ParseCIGAR(phred, (int)(SNPpos.Key - leftStart), cigar));
                if (qualityOfVariant < (int)PSAMBIN2CONTIGS.Thresholds.MinSnpPhred)
                {
                    PSAMBIN2CONTIGS.filteredCounts.AddOrUpdate("Snp Phred quality < " + PSAMBIN2CONTIGS.Thresholds.MinSnpPhred, t => 1, (a, b) => b + 1);
                    rightVariant = 4;
                    return;
                }
            }
            var Variant = ParseCIGAR(seq, (int)(SNPpos.Key - leftStart), cigar);
            leftVariant = SNPpos.Value.TupleId(Variant);
        }

        public void SetRightVariant(string seq, string cigar, string phred=null)
        {
            rightEnd = rightStart + ParseCIGARlen(cigar);
            if (!IsSpanningSnp(rchr, rightStart, rightEnd, out rightSNPidx))
            {
                rightVariant = 3;
                return;
            }
            var SNPpos = Program.SNPmap[rchr - 1].ElementAt(rightSNPidx);
            if (string.IsNullOrEmpty(phred))
            {
                var qualityOfVariant = CharToPhred(ParseCIGAR(phred, (int) (SNPpos.Key - rightStart), cigar));
                if (qualityOfVariant < (int)PSAMBIN2CONTIGS.Thresholds.MinSnpPhred)
                {
                    PSAMBIN2CONTIGS.filteredCounts.AddOrUpdate("Snp Phred quality < " + PSAMBIN2CONTIGS.Thresholds.MinSnpPhred, t => 1, (a, b) => b + 1);
                    rightVariant = 4;
                    return;
                }
            }
            var Variant = ParseCIGAR(seq, (int)(SNPpos.Key - rightStart), cigar);
            rightVariant = SNPpos.Value.TupleId(Variant);
        }

        public static int CharToPhred(char chr)
        {
            return ((int) chr) - ((int) '!');
        }

        private bool IsSpanningSnp(short chr, long start, long end, out int nearestSNPidx)
        {
            //Determine if contains SNP:
            nearestSNPidx = Program.SNPmap[chr - 1].BinarySearch(start);
            if (nearestSNPidx >= Program.SNPmap[ chr - 1].Count) return false;
            long locPostSNP = Program.SNPmap[ chr - 1].ElementAt(nearestSNPidx).Key;
            return locPostSNP >= start && locPostSNP < end;
        }

        

        private uint ParseCIGARlen(string cigar)
        {
            uint genomicLoc = 0;
            var chunks = Helpers.CigarRegex.Matches(cigar);
            foreach (Match chunk in chunks)
            {
                var chunkLen = Convert.ToByte(chunk.Value.Substring(0, chunk.Value.Length - 1));
                var chunkType = chunk.Value.Last();
                switch (chunkType)
                {
                    case 'M':
                        //alignment match (can be a sequence match or mismatch)
                        genomicLoc += chunkLen;
                        break;
                    case 'I':
                        //insertion to the reference
                        genomicLoc += chunkLen; //Requested genomic position is shifted in this read
                        break;
                    case 'D':
                        //deletion from the reference
                        genomicLoc -= chunkLen; //Requested genomic position is shifted in this read
                        break;
                    case 'N':
                        //skipped region from the reference
                        break;
                    case 'S':
                        // soft clipping (clipped sequences present in SEQ)
                        genomicLoc += chunkLen; //Requested genomic position is shifted in this read
                        break;
                    case 'H':
                        //hard clipping (clipped sequences NOT present in SEQ)
                        break;
                    case 'P':
                        //padding (silent deletion from padded reference)
                        break;
                    case '=':
                        //Sequence match
                        break;
                    case 'X':
                        //Sequence mismatch
                        break;
                    default:
                        Helpers.Write((int)PrintDetail.Debug, String.Format("Unsupported character ({0}) in CIGAR string!", chunkType));
                        break;
                }
             }
            return genomicLoc;
        }

        
        public static char ParseCIGAR(string seq, int pos, string cigar)
        {
            var shift = 0;
            var genomicLoc = 0;
            var chunks = Helpers.CigarRegex.Matches(cigar);
            foreach (Match chunk in chunks)
            {
                var chunkLen = Convert.ToByte(chunk.Value.Substring(0, chunk.Value.Length - 1));
                var chunkType = chunk.Value.Last();
                switch (chunkType)
                {
                    case 'M':
                        //alignment match (can be a sequence match or mismatch)
                        genomicLoc += chunkLen;
                        break;
                    case 'I':
                        //insertion to the reference
                        genomicLoc += chunkLen; //Requested genomic position is shifted in this read
                        shift += chunkLen;
                        break;
                    case 'D':
                        //deletion from the reference
                        genomicLoc -= chunkLen; //Requested genomic position is shifted in this read
                        shift -= chunkLen;
                        break;
                    case 'N':
                        //skipped region from the reference
                        break;
                    case 'S':
                        // soft clipping (clipped sequences present in SEQ)
                        genomicLoc += chunkLen; //Requested genomic position is shifted in this read
                        break;
                    case 'H':
                        //hard clipping (clipped sequences NOT present in SEQ)
                        break;
                    case 'P':
                        //padding (silent deletion from padded reference)
                        break;
                    case '=':
                        //Sequence match
                        break;
                    case 'X':
                        //Sequence mismatch
                        break;
                    default:
                        Helpers.Write((int)PrintDetail.Debug, String.Format("Unsupported character ({0}) in CIGAR string!", chunkType));
                        return ' ';
                }
                if (genomicLoc < pos) continue;
                return seq.Length > pos + shift && pos + shift >= 0 ? seq[pos + shift] : ' ';
            }
            return ' ';
        }


        public bool MissingLeft()
        {
            return this.leftVariant == 0;
        }

        public bool MissingRight()
        {
            return this.rightVariant == 0;
        }

        public bool HasUnknownVariant()
        {
            return (this.leftVariant == 4 || this.rightVariant == 4);
        }

    }
}
