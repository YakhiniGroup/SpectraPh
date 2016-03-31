using System;
using System.Collections.Concurrent;
using System.Collections.Generic;

namespace SpectraPh.DataStructures
{
    public class BlockMatrix
    {
        //Our goal is to have a list of locis, bin-centers, which are either a SNP, or a SNPless region covered by a number of fixed reads such that near bins are not identical in count.
        private ConcurrentDictionary<Tuple<long, long>, HashSet<string>>[,] PosReadsConcurrent;
        public SortedList<long, SortedList <long,HashSet<string>>>[,] PosReads; //Genomic position, Set of read names which start here
        public ConcurrentDictionary<string, PairedRead> NamedReads; //Read name, position, SNP variants (- if none), 
        private int CHROMOSOMES;
        
        public BlockMatrix(int chr)
        {
            CHROMOSOMES = chr;
            PosReads = new SortedList<long, SortedList<long, HashSet<string>>>[CHROMOSOMES, CHROMOSOMES];
            PosReadsConcurrent = new ConcurrentDictionary<Tuple<long, long>, HashSet<string>>[CHROMOSOMES, CHROMOSOMES];
            for (var i = 0; i < CHROMOSOMES; i++)
                for (var j = 0; j < CHROMOSOMES; j++)
                {
                    PosReads[i, j] = new SortedList<long, SortedList<long, HashSet<string>>>();
                    PosReadsConcurrent[i, j] = new ConcurrentDictionary<Tuple<long, long>, HashSet<string>>();
                }
            NamedReads = new ConcurrentDictionary<string, PairedRead>();
        }

        public void AddPair(string pairName, long posA, long posB, int chrA, int chrB, string readContents, string CIGAR)
        {
            PosReadsConcurrent[chrA - 1, chrB - 1].AddOrUpdate(new Tuple<long, long>(posA, posB),
                t => new HashSet<string>() { pairName },
                (a, b) =>
                {
                    b.Add(pairName);
                    return b;
                });

            PosReadsConcurrent[chrB - 1, chrA - 1].AddOrUpdate(new Tuple<long, long>(posB, posA),
                t => new HashSet<string>() { pairName },
                (a, b) =>
                {
                    b.Add(pairName);
                    return b;
                });

            NamedReads.AddOrUpdate(pairName, t =>
            {
                var pr = new PairedRead(pairName);
                pr.SetEnd(chrA, posA, readContents, CIGAR);
                return pr;
            }, (a, b) =>
            {
                b.SetEnd(chrA, posA, readContents, CIGAR);
                return b;
            });
        }

        public void FinalizeAdded()
        {
            for(var i=0; i<CHROMOSOMES;i++)
                for (var j = 0; j < CHROMOSOMES; j++)
                {
                    PosReads[i,j]= new SortedList<long, SortedList<long, HashSet<string>>>();
                    foreach (var kvp in PosReadsConcurrent[i, j])
                    {
                        if (!PosReads[i,j].ContainsKey(kvp.Key.Item1))
                            PosReads[i,j].Add(kvp.Key.Item1,new SortedList<long, HashSet<string>>());
                        PosReads[i, j][kvp.Key.Item1].Add(kvp.Key.Item2, kvp.Value);
                    }
                }
            PosReadsConcurrent = null;
        }
        
    }
}
