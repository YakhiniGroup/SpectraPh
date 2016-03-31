using System;
using System.Linq;
using System.Text.RegularExpressions;

namespace SpectraPh.DataStructures
{
    public class MateRead
    {
        public int Length;
        public short Chromosome;
        public long Position;
        public char Variant = '?';
        public byte VariantId = 0;
        private byte shift = 0;

        public long EndPosition
        {
            get
            {
                return Position + Length + shift;
            }
            private set { }
        }

        public MateRead(string seq)
        {
            
        }

        public MateRead(short chr, long pos, string seq, string CIGAR)
        {
            //Determine if contains SNP:
            int nearestSNPidx = Program.SNPmap[(short)chr - 1].BinarySearch(pos);
            long locPostSNP = Program.SNPmap[(short)chr - 1].ElementAt(nearestSNPidx).Key;
            if (locPostSNP >= pos && locPostSNP < (pos + seq.Length))
            {
                //We have a read spanning a SNP! determine variant
                var SNPpos = Program.SNPmap[(short)chr - 1].ElementAt(nearestSNPidx).Key;
                ParseCIGAR(seq, (int)(SNPpos - pos), CIGAR);
                var tVariantId = Program.SNPmap[(short)chr - 1].ElementAt(nearestSNPidx).Value.TupleId(Variant);
                VariantId = tVariantId < 3 ? tVariantId : (byte) 0;
                if (VariantId == 0)
                {
                    Program.unmappedEnd++;
                    Helpers.Write((int)PrintDetail.Debug, string.Format("{0}-{1}-{2}\t{3}\t{4}",
                        seq.Substring(Math.Max(0, (int)(SNPpos - pos - 5)), 4), Variant,
                        seq.Substring(Math.Min(seq.Length, (int)(SNPpos - pos + 1)),
                            Math.Max(0, Math.Min((int)(seq.Length - (SNPpos - pos + 5)), 4))),
                        Program.SNPmap[(short)chr - 1].ElementAt(nearestSNPidx).Value, CIGAR));

                    Helpers.Write((int)PrintDetail.Debug, string.Format(
                            "Found variant not considered in haplotype. Dropping it from list ({0}/{1}).", Program.mappedEnd,
                            Program.mappedEnd + Program.unmappedEnd));
                    VariantId = 3;
                }
                else
                {
                    Program.mappedEnd++;
                }
            }
        }

        public void ParseCIGAR(string seq, int pos, string cigar)
        {
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
                        Variant = ' ';
                        break;
                }
                if (genomicLoc >= pos) //we have reached the chunk containing pos.
                {
                    if (seq.Length > pos + shift)
                        Variant = seq[pos + shift];
                    else
                        Variant = ' ';
                    break;
                }
            }
            Variant = ' ';
        }
    }

    public class PairedRead
    {
        //public string readName;//, seqA, seqB;
        public MateRead Mate1, Mate2;
        private int _endCounter = 0;
        

        public PairedRead(string readname)
        {
            //readName = readname;
        }

        public bool IsComplete()
        {
            return (Mate1 != null && Mate2!=null);
        }

        public MateRead GetLastEnd()
        {
            switch (_endCounter)
            {
                case 1:
                    return Mate1;
                case 2:
                    return Mate2;
                default:
                    return null;
            }
        }

        public MateRead SetEnd(int chr, long pos, string seq, string CIGAR)
        {
            _endCounter++;

            
            switch (_endCounter)
            {
                case 1:
                    Mate1 = new MateRead((short)chr, pos, seq, CIGAR);
                    return Mate1;
                    if (Mate1.VariantId == 3) _endCounter = 3;
                    break;
                case 2:
                    Mate1 = new MateRead((short)chr, pos, seq, CIGAR);
                    if (Mate2.VariantId == 3) _endCounter = 3;
                    return Mate2;
                    break;
                default:
                    //throw new Exception("Read has more than 2 ends.");
                    break;
            }
            return null;
        }

        


        internal int GetRelevantEndLength(int chr, long pos)
        {
            return (chr == Mate1.Chromosome && pos == Mate1.Position)
                ? Mate1.Length
                : (chr == Mate2.Chromosome && pos == Mate2.Position) ? Mate2.Length : -1;
        }

        public byte GetEndVariantId(int chr, long pos)
        {
            return (chr == Mate1.Chromosome && pos == Mate1.Position)
                ? Mate1.VariantId
                : (chr == Mate2.Chromosome && pos == Mate2.Position) ? Mate2.VariantId : (byte)0;
        }

        public MateRead GetRelevantEnd(int chr, long pos)
        {
            if (_endCounter > 2)
                return null;
            return (chr == Mate1.Chromosome && pos == Mate1.Position)
                ? Mate1
                : (chr == Mate2.Chromosome && pos == Mate2.Position) ? Mate2 : null;
        }
    }
}
