using System;
using System.Collections;
using System.Collections.Generic;

namespace SpectraPh.DataStructures
{
    public class SNPMap
    {
        private SortedList snpMap;
        private List<long> snpList;
        private Dictionary<long, bool> phasedSNPlist;
        private bool issorted = false;

        public SNPMap()
        {
            var comp = new snpMapCompare();
            snpMap = new SortedList(comp);
            snpList = new List<long>();
            phasedSNPlist=new Dictionary<long, bool>();
        }

        public void Add(long loc, Tuple<char, char> variants)
        {
            snpMap.Add(loc, variants);
            snpList.Add(loc);
        }
        public void SetPhase(long pos, Tuple<char, char> ord)
        {
            snpMap.SetByIndex(snpMap.IndexOfKey(pos), ord);
        }

        public void FinalizeSort()
        {
            snpList.Sort();
            snpMap = SortedList.Synchronized(snpMap);
            issorted = true;
        }

        public int BinarySearch(long posA)
        {
            if (!issorted)
                throw new Exception("Did not finalize SNPmap!");
            //Determine if contains SNP:
            var nearestSNPidx = snpList.BinarySearch(posA);
            nearestSNPidx = nearestSNPidx > 0 ? nearestSNPidx : ~nearestSNPidx;
            //find next snp index
            if (nearestSNPidx < snpList.Count)
                nearestSNPidx = ElementAt(nearestSNPidx).Key < posA ? nearestSNPidx + 1 : nearestSNPidx;
            return nearestSNPidx;
        }

        public KeyValuePair<long, Tuple<char, char>> ElementAt(int i)
        {
            return snpMap.Count <= i
                ? new KeyValuePair<long, Tuple<char, char>>(long.MaxValue, new Tuple<char, char>('-', '-'))
                : new KeyValuePair<long, Tuple<char, char>>((long) snpMap.GetKey(i),
                    (Tuple<char, char>) snpMap.GetByIndex(i));
        }

        public int Count { get { return snpList.Count; } }
        public List<long> GetList()
        {
            return snpList;
        }

        public void RemoveSnp(long pos)
        {
            if (snpList.Contains(pos))
                snpList.Remove(pos);
            if (snpMap.ContainsKey(pos))
                snpMap.Remove(pos);
        }
    }
    public class snpMapCompare : IComparer
    {
        public int Compare(object x, object y)
        {
            var ox = (long)x;
            var oy = (long)y;
            return ox.CompareTo(oy);
        }

    }


}
