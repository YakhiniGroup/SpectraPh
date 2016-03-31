using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.Configuration;
using System.IO;
using System.Linq;
using System.Threading;
using System.Threading.Tasks;
using SpectraPh.DataStructures;

namespace SpectraPh.Actions
{
    public class PSAMBIN2CONTIGS : IHaploParserAction
    {
        public class Thresholds
        {
            public static long ContigLengths = 2000000;
            public static double ConvolutionSigmaDensity = 0.2;
            public static long REMaxDistance = 10000000000;
            public static int MinAvgPhred = 10;
            public static int MinSnpPhred = 10;
            public static int ReadMinDistances = 1000;
            public static int ReadMinBinDistances = 1;
            public static bool PhaseNonSnps = false;
        }

        private static Dictionary<int,Dictionary<long, Tuple<int, char, char>>> KnownPhasing; //chr->pos->[contigid,ref,alt]
        private static Dictionary<int, List<KeyValuePair<long, int>>> contigPos;
        private static Dictionary<int, int> contigIds;
        private static ConcurrentDictionary<int, long>[] contigLengths;

        private static kvpCompare kvpcomp = new kvpCompare();
        private static long Sigma;
        private static Dictionary<int, List<uint>> REsites;
        public static ConcurrentDictionary<string, int> filteredCounts;

        public string Usage()
        {
            string message = "<PSAM InputFolder> <HAP InputFolder> <Partial Phasing> <OutputFolder>";
            message += Environment.NewLine + "Iterates on PSAM files and HAPlotype maps, given a , generating Sparse Hi-C read maps.";
            return message;
        }

        public string Command()
        {
            return "PSAMBIN2CONTIGS";
        }

        public void Go(string[] args)
        {
            Program.CHROMOSOMES = Convert.ToInt32(ConfigurationManager.AppSettings["PARAM:Chromosomes"]);
            Program.PLOIDITY = Convert.ToInt32(ConfigurationManager.AppSettings["PARAM:Ploidity"]);
            Program.PRINTDETAIL = (int)Enum.Parse(typeof(PrintDetail), ConfigurationManager.AppSettings["PARAM:PrintDetail"]);

            Thresholds.ConvolutionSigmaDensity = Convert.ToDouble(ConfigurationManager.AppSettings["PARAM:ConvolutionSigmaDensity"]);
            Thresholds.MinAvgPhred = Convert.ToInt32(ConfigurationManager.AppSettings["PARAM:MinAvgPhred"]);
            Thresholds.MinSnpPhred = Convert.ToInt32(ConfigurationManager.AppSettings["PARAM:MinSnpPhred"]);
            Thresholds.ReadMinDistances = Convert.ToInt32(ConfigurationManager.AppSettings["PARAM:ReadMinDistances"]);
            Thresholds.ReadMinBinDistances = Convert.ToInt32(ConfigurationManager.AppSettings["PARAM:ReadMinBinDistances"]);
            Thresholds.PhaseNonSnps = Convert.ToBoolean(ConfigurationManager.AppSettings["PARAM:PhaseNonSnps"]);
        

            if (Helpers.GetFilesFromPath(args[0], "*.psam").Any())
            {
                
                Helpers.Write((int)PrintDetail.Important, "Initializing data structures.");
                InitializeMaps();
                
                Helpers.Write((int)PrintDetail.Important, "Initializing data structures.");
                LoadREsites(args);

                Helpers.Write((int)PrintDetail.Important, "Loading haplotype SNPs.");
                LoadSNPs(args);

                Helpers.Write((int)PrintDetail.Important, "Loading partial phasing.");
                //LoadPphasing(args);
                LoadHTBlockPositions(args);
                for (var i = 0; i < Program.CHROMOSOMES; i++)
                    Program.SNPmap[i].FinalizeSort();
                Sigma = (long)(Thresholds.ConvolutionSigmaDensity * contigLengths.SelectMany(t => t.Values).Average());

                Helpers.Write((int)PrintDetail.Important, "Binning PSAM (Paired-SAM) file by contigs.");
                SerializeScatteredMapFromPSAMs(args);


                Helpers.Write((int)PrintDetail.Important, "Done binning PSAM files.");
            }
            else
                Helpers.Write((int)PrintDetail.Important, "No files found.");
        }


        private static void LoadHTBlockPositions(string[] args)
        {
            contigPos = new Dictionary<int, List<KeyValuePair<long, int>>>();
            contigIds = new Dictionary<int, int>();
            contigLengths = new ConcurrentDictionary<int, long>[Program.CHROMOSOMES];
            using (var fileStream = new StreamReader(args[2]))
                while (!fileStream.EndOfStream)
                {
                    string[] splitline = fileStream.ReadLine().Split(',');
                    int chr = Helpers.ChromNum(splitline[0]);
                    int HTblock = int.Parse(splitline[1]);
                    long pos = long.Parse(splitline[2]);
                    long posEnd = long.Parse(splitline[3]);
                    var alleleOrder = new Tuple<char, char>(splitline[3][0], splitline[4][0]);
                    if (!contigPos.ContainsKey(chr))
                    {
                        contigPos.Add(chr, new List<KeyValuePair<long, int>>());
                        contigIds.Add(chr, 0);
                        contigLengths[chr-1] = new ConcurrentDictionary<int, long>();
                    }
                    contigPos[chr].Add(new KeyValuePair<long, int>());
                    contigIds[chr] = Math.Max(HTblock, contigIds[chr]);
                    contigLengths[chr-1].GetOrAdd(contigIds[chr], t => posEnd-pos);
                }
            
        }

        private static void LoadPphasing(string[] args)
        {
            KnownPhasing = new Dictionary<int, Dictionary<long, Tuple<int, char, char>>>();
            var GroundTruth = new Dictionary<int, Dictionary<long, Tuple<char, char>>>();
            Helpers.Write((int)PrintDetail.Detailed, "Loading ground truth file:" + args[2]);
            using (var fileStream = new StreamReader(args[2]))
                while (!fileStream.EndOfStream)
                {
                    string[] splitline = fileStream.ReadLine().Split(',');
                    if (splitline[0].StartsWith("#")) continue;
                    int chr = Helpers.ChromNum(splitline[0]);
                    long pos = long.Parse(splitline[1]);
                    if (!GroundTruth.ContainsKey(chr)) GroundTruth.Add(chr, new Dictionary<long, Tuple<char, char>>());
                    GroundTruth[chr].Add(pos, new Tuple<char, char>(splitline[2][0], splitline[3][0]));
                }

            Helpers.Write((int)PrintDetail.Detailed, "Building contigs.");
            var loccontigPos = new Dictionary<int, SortedList<long, int>>();
            contigPos = new Dictionary<int, List<KeyValuePair<long, int>>>();
            contigIds=new Dictionary<int, int>();
            contigLengths = new ConcurrentDictionary<int, long>[GroundTruth.Keys.Count];
            foreach (var dct in GroundTruth)
            {
                var nonPhased = Program.SNPmap[dct.Key - 1].GetList().Except(dct.Value.Keys).ToList();
                foreach (var snpPos in nonPhased)
                    Program.SNPmap[dct.Key - 1].RemoveSnp(snpPos);
                contigIds.Add(dct.Key,0);
                bool incontig = false;
                KnownPhasing.Add(dct.Key, new Dictionary<long, Tuple<int, char, char>>());
                loccontigPos[dct.Key] = new SortedList<long, int>();
                contigLengths[dct.Key-1] = new ConcurrentDictionary<int, long>();;
                long maxcontigpos = 0;
                for (var i = 0; i < Program.SNPmap[dct.Key - 1].Count; i++)
                {
                    var pos = Program.SNPmap[dct.Key - 1].ElementAt(i).Key;
                    //var isphased = dct.Value.ContainsKey(pos);
                    //if (isphased)
                    //{
                    Program.SNPmap[dct.Key - 1].SetPhase(pos, dct.Value[pos]); //Set order by known phasing
                    if (!incontig)
                    {
                        incontig = true;
                        contigIds[dct.Key]++;
                        loccontigPos[dct.Key].Add(pos, contigIds[dct.Key]);
                    }
                    maxcontigpos = Math.Max(maxcontigpos, pos);
                    KnownPhasing[dct.Key].Add(pos, new Tuple<int, char, char>(contigIds[dct.Key], dct.Value[pos].Item1, dct.Value[pos].Item2));
                    //}
                    if (loccontigPos[dct.Key].Last().Key + Thresholds.ContigLengths > pos) //Grace.
                    {
                        maxcontigpos = Math.Max(maxcontigpos, pos);
                    }
                    else
                    {
                        incontig = false;
                        if (maxcontigpos >= loccontigPos[dct.Key].Last().Key)
                            contigLengths[dct.Key-1].GetOrAdd(contigIds[dct.Key],
                                t => maxcontigpos - loccontigPos[dct.Key].Last().Key);
                        if (!loccontigPos[dct.Key].ContainsKey(maxcontigpos))
                            loccontigPos[dct.Key].Add(maxcontigpos, contigIds[dct.Key]);
                    }
                }
                contigPos.Add(dct.Key, new List<KeyValuePair<long, int>>());
                contigPos[dct.Key] =
                    new List<KeyValuePair<long, int>>(
                        loccontigPos[dct.Key].Select(t => new KeyValuePair<long, int>(t.Key, t.Value)));
                contigPos[dct.Key].Sort((x, y) => x.Key.CompareTo(y.Key));
                var chromcontigs = KnownPhasing[dct.Key].GroupBy(t => t.Value.Item1)
                    .Select(t => new {Contig = t.Key, First = t.Min(u => u.Key), Last = t.Max(u => u.Key), Size=t.Count()})
                    .OrderBy(t => t.Contig).ToList();
                if (chromcontigs.Any())
                {
                    Task.Run(() =>
                    {
                        Helpers.Write((int)PrintDetail.Important, string.Format("Serializing contig map for {0}. phased loci={1} contigs={2}", dct.Key, dct.Value.Count, contigIds[dct.Key]));
                        using (var fileout = new StreamWriter(Path.Combine(Path.GetFullPath(args[3]), string.Format("contigMap.{0}.{1}.csv", dct.Key, Thresholds.ContigLengths))))
                            foreach (var contig in chromcontigs)
                                fileout.WriteLine("{0},{1},{2}", contig.Contig, contig.First, contig.Last);
                    });
                }
            }
        }
                
        private static void LoadSNPs(string[] args)
        {
            //Parallel.ForEach(Helpers.GetFilesFromPath(args[1], "*.hap"), filename =>
            foreach (string filename in Helpers.GetFilesFromPath(args[1], "*.hap"))
            {
                Helpers.Write((int) PrintDetail.Detailed, "Loading file:" + filename);
                int chr = -1;
                using (var fileStream = new StreamReader(filename))
                    while (!fileStream.EndOfStream)
                    {
                        string[] splitline = fileStream.ReadLine().Split('\t');
                        chr = Helpers.ChromNum(splitline[0]);
                        Program.SNPmap[chr - 1].Add(Convert.ToInt64(splitline[1]),
                            new Tuple<char, char>(splitline[3].First(), splitline[4].First()));
                    }
            }//);

        }

        private static void InitializeMaps()
        {
            Program.SNPmap = new SNPMap[Program.CHROMOSOMES];
            for (int i = 0; i < Program.CHROMOSOMES; i++)
                Program.SNPmap[i] = new SNPMap();
            Program.Matrix = new BlockMatrix(Program.CHROMOSOMES);
        }

        private static void LoadREsites(string[] args)
        {
            var files = Helpers.GetFilesFromPath(args[4], "*.resites.txt").ToList();
            REsites = new Dictionary<int, List<uint>>();
            foreach (var filename in files)
            {
                int chr;
                using (var file = new StreamReader(filename))
                {
                    var header = file.ReadLine().Substring(1);
                    chr = Helpers.ChromNum(header);
                    REsites.Add(chr, new List<uint>());
                    while (!file.EndOfStream)
                    {
                        REsites[chr].Add(uint.Parse(file.ReadLine()));
                    }
                }
                REsites[chr].Sort();
            }
        }

        private static long DistanceFromClosestREsite(int chr, uint pos)
        {
            var idx = REsites[chr].BinarySearch(pos);
            idx = idx >= 0 ? idx : ~idx;
            idx = Math.Min(idx, REsites[chr].Count-1);
            var res = Math.Abs(REsites[chr][idx] - pos);
            return res;
        }

        private static void SerializeScatteredMapFromPSAMs(string[] args)
        {
            var files = Helpers.GetFilesFromPath(args[0], "*.psam").ToList();

            foreach (var file in files)
            {
                #region Read file to data structure
                filteredCounts = new ConcurrentDictionary<string, int>();
                var filenameSplit = Path.GetFileNameWithoutExtension(file).Split('_');

                // Our thread-safe collection used for the handover.
                var lines = new BlockingCollection<string>(20000);
                var chrA = Helpers.ChromNum(filenameSplit[0]);
                var chrB = Helpers.ChromNum(filenameSplit[1]);
                //Thread for parsing file data into region-tracking data structure
                int numReads = 0, badReadCount = 0; //interlocked to convert readName to number (less ram)
                var concurrentNameTracker = new ConcurrentDictionary<string, int>(50, 1000000);
                var regionTracker = new ConcurrentDictionary<int, SquareRegion>(50, 1000000);
                //defines region parameters
                //var sparseMat = new ConcurrentDictionary<Tuple<int, int>, int>[Program.PLOIDITY + 1, Program.PLOIDITY + 1];
                var sparseSnpMap = new ConcurrentDictionary<Tuple<int, int>, ConcurrentDictionary<Tuple<int, int>, int>>[Program.PLOIDITY + 1, Program.PLOIDITY + 1];
                var sparseMatMapper = new ConcurrentDictionary<Tuple<int, int>, LockableDictionary<Tuple<double, double>, int>>[Program.PLOIDITY + 1, Program.PLOIDITY + 1];
                for (var i = 0; i < Program.PLOIDITY + 1; i++)
                    for (var j = 0; j < Program.PLOIDITY + 1; j++)
                    {
                        //sparseMat[i, j] = new ConcurrentDictionary<Tuple<int, int>, int>();
                        sparseSnpMap[i, j] = new ConcurrentDictionary<Tuple<int, int>, ConcurrentDictionary<Tuple<int, int>, int>>();
                        sparseMatMapper[i, j] = new ConcurrentDictionary<Tuple<int, int>, LockableDictionary<Tuple<double, double>, int>>();
                    }
                var ignoredReads = new ConcurrentBag<string>();
                var tasklst = new List<Task>();
                var positions = new ConcurrentDictionary<Tuple<int, int>, byte>();
                Helpers.Write((int)PrintDetail.Important, string.Format("Reading pair {0},{1}.", chrA, chrB));

                //Thread for file reading into blockingCollection
                var stage1 = Task.Run(() =>
                {
                    using (var reader = new StreamReader(file))
                    {
                        string line;
                        while ((line = reader.ReadLine()) != null)
                        {
                            ParseRead(line, concurrentNameTracker, ref numReads, regionTracker, ignoredReads, chrA, chrB,
                                positions, sparseSnpMap, sparseMatMapper, ref badReadCount);
                            //lines.Add(line);
                        }
                    }
                    lines.CompleteAdding();
                });

                var stage2 = Task.Run(() =>
                {
                    foreach (var line in lines.GetConsumingEnumerable())
                    {
                        if (tasklst.Count > 200)
                        {
                            var fid = Task.WaitAny(tasklst.ToArray());
                            tasklst.RemoveAt(fid);
                        }
                        var line1 = line;
                        tasklst.Add(Task.Run(() => ParseRead(line1, concurrentNameTracker, ref numReads, regionTracker, ignoredReads, chrA,
                            chrB, positions, sparseSnpMap, sparseMatMapper, ref badReadCount)));
                    }
                    Task.WaitAll(tasklst.ToArray());
                });
                Console.WriteLine("");
                using (var dbgfile = new StreamWriter("debug.ignoredVariants.csv"))
                    for (var i = 0; i < badReadCount; i++)
                    {
                        string item;
                        if (ignoredReads.TryTake(out item))
                            dbgfile.WriteLine(item);
                    }
                
                Task.WaitAll(stage1, stage2);
                regionTracker.Clear();
                Helpers.Write((int)PrintDetail.Important,
                    string.Format("Done reading pair {0},{1}. Threw out {2} unknown variants, had {3} reads. Furtherst snp bins = {4}", chrA, chrB, badReadCount, numReads, sparseSnpMap[0, 0].Max(t => t.Key.Item2 - t.Key.Item1)));


                #endregion

                #region Convert data structure for output

                var th = new Thresholds();
                using (
                    var fileout =
                        new StreamWriter(Path.Combine(Path.GetFullPath(args[3]),
                            string.Format("filterCounts_{0}.{1}.{2}.csv", Thresholds.ContigLengths, chrA, chrB))))
                {
                    fileout.WriteLine("Params:");
                    foreach (var field in typeof(Thresholds).GetFields())
                        fileout.WriteLine("{0},{1}", field.Name, field.GetValue(th));
                    fileout.WriteLine("\nFilterCounts:");
                    foreach (var kvp in filteredCounts)
                        fileout.WriteLine("{0},{1}",kvp.Key,kvp.Value);
                    fileout.WriteLine("\nThrew out {0} / {1} reads.", badReadCount, numReads);
                }

                using (var fileout = new StreamWriter(Path.Combine(Path.GetFullPath(args[3]), string.Format("snpReadCounts_{0}.{1}.{2}.csv", Thresholds.ContigLengths, chrA, chrB))))
                    for (var i = 0; i < Program.PLOIDITY; i++)
                    {
                        var baseA = i * contigIds[chrA];
                        for (var j = 0; j < Program.PLOIDITY; j++)
                        {
                            var baseB = j * contigIds[chrB];
                            foreach (var row in sparseMatMapper[i, j].OrderBy(t => t.Key.Item1).ThenBy(t => t.Key.Item2))
                                fileout.WriteLine("{0},{1},{2},{3}", baseA + row.Key.Item1, baseB + row.Key.Item2,
                                    sparseSnpMap[i,j][row.Key].Count, row.Value.Sum(t=>t.Value));
                        }

                    }
                var phasedMat = PhaseByGaussianConvolution(sparseMatMapper);

                //Collect relevant columns per row for processing in blockingCollection

                Helpers.Write((int)PrintDetail.Important, string.Format("Serializing pair {0},{1}.", chrA, chrB));
                
                using (var fileout = new StreamWriter(Path.Combine(Path.GetFullPath(args[3]), string.Format("phaseMat_{0}.{1}.{2}.csv", Thresholds.ContigLengths, chrA, chrB))))
                    for (var i = 0; i < Program.PLOIDITY; i++)
                    {
                        var baseA = i * contigIds[chrA];
                        for (var j = 0; j < Program.PLOIDITY; j++)
                        {
                            var baseB = j * contigIds[chrB];
                            foreach (var row in phasedMat[i, j].OrderBy(t => t.Key.Item1).ThenBy(t => t.Key.Item2))
                                fileout.WriteLine("{0},{1},{2:R}", baseA + row.Key.Item1, baseB + row.Key.Item2,
                                    row.Value/getNormalizer(row.Key.Item1, row.Key.Item2, chrA, chrB));
                        }

                    }
                #endregion
            }

            Helpers.Write((int)PrintDetail.Important, string.Format("Done serializing."));
        }

        private static void ParseRead(string line1, ConcurrentDictionary<string, int> concurrentNameTracker, ref int numReads,
            ConcurrentDictionary<int, SquareRegion> regionTracker, ConcurrentBag<string> ignoredReads, int chrA, int chrB, ConcurrentDictionary<Tuple<int, int>, byte> positions,
            ConcurrentDictionary<Tuple<int, int>, ConcurrentDictionary<Tuple<int, int>, int>>[,] sparseSnpMap, ConcurrentDictionary<Tuple<int, int>, LockableDictionary<Tuple<double, double>, int>>[,] sparseMatMapper, ref int badReadCount)
        {
            var splitline = line1.Split('\t');
            //Parse line info
            var rname = splitline[0];
            var readId = concurrentNameTracker.GetOrAdd(rname, Interlocked.Increment(ref numReads));
            //Read name to id
            var mychr = Helpers.ChromNum(splitline[1]);
            var mypos = Convert.ToUInt32(splitline[2]);
            var cigar = splitline[3];
            var otherchr = splitline[4].Equals("=") ? mychr : Helpers.ChromNum(splitline[4]);
            var otherpos = Convert.ToUInt32(splitline[5]);
            var myseq = splitline[6];
            var seqquality = splitline[7];
            var avgPhred = seqquality.Select(SquareRegion.CharToPhred).Average();
            bool badread = false;
            //Add to / update data structure
            var added = regionTracker.GetOrAdd(readId, t => new SquareRegion());
            Tuple<int, int> locus = null;
            lock (added.locker)
            {
                if (mychr > otherchr || (mychr == otherchr && mypos >= otherpos))
                {
                    // this is right side
                    added.SetOrderOfVariants((short) otherchr, (short) mychr, otherpos, mypos);
                    added.SetRightVariant(myseq, cigar, seqquality);
                }
                else
                {
                    added.SetOrderOfVariants((short) mychr, (short) otherchr, mypos, otherpos);
                    added.SetLeftVariant(myseq, cigar, seqquality);
                }

                if (added.MissingLeft() || added.MissingRight()) //not finished parsing both mates 
                {
                    filteredCounts.AddOrUpdate("Missing End", t => 1, (a, b) => b + 1);
                    return;
                }
                filteredCounts.AddOrUpdate("Missing End", t => 0, (a, b) => b - 1);
                if (added.HasUnknownVariant())
                {
                    filteredCounts.AddOrUpdate("Unknown Variant", t => 1, (a, b) => b + 1);
                    badread = true;
                }
                if (mychr == otherchr && (added.rightStart-added.leftEnd) < Thresholds.ReadMinDistances)// && added.leftVariant == added.rightVariant && added.leftVariant<3)
                {
                    if (added.rightStart < added.leftEnd)
                        throw new Exception("wrong end order!");
                    filteredCounts.AddOrUpdate("Read distance < " + Thresholds.ReadMinDistances, t => 1, (a, b) => b + 1);
                    badread = true;
                }
                
                if (avgPhred < (int)Thresholds.MinAvgPhred)
                {
                    filteredCounts.AddOrUpdate("Average Phred < " + Thresholds.MinAvgPhred, t => 1, (a, b) => b + 1);
                    badread = true;
                }
                locus = new Tuple<int, int>(GetContigId(added.lchr, added.leftStart),
                    GetContigId(added.rchr, added.rightStart));
                if (mychr == otherchr && locus.Item1 >= 0 && locus.Item2 >= 0 && (locus.Item2 - locus.Item1) < Thresholds.ReadMinBinDistances)
                {
                    if (locus.Item2 < locus.Item1)
                        throw new Exception("wrong end order!");
                    filteredCounts.AddOrUpdate("Read bin distance < " + Thresholds.ReadMinBinDistances, t => 1, (a, b) => b + 1);
                    badread = true;
                }
                if (locus.Item1 == -1 || locus.Item2 == -1) //not in contig
                {
                    filteredCounts.AddOrUpdate("Read not mapped to contig", t => 1, (a, b) => b + 1);
                    badread = true;
                }
                if (DistanceFromClosestREsite(added.lchr, added.leftStart) > (long)Thresholds.REMaxDistance ||
                    DistanceFromClosestREsite(added.rchr, added.rightStart) > (long)Thresholds.REMaxDistance)
                {
                    filteredCounts.AddOrUpdate("Read distance from RE site > " + Thresholds.REMaxDistance, t => 1, (a, b) => b + 1);
                    //badread = true;
                }
                if (badread) //Bad read
                {
                    SquareRegion junkregion;
                    int junkid;
                    regionTracker.TryRemove(readId, out junkregion);
                    concurrentNameTracker.TryRemove(rname, out junkid);
                    ignoredReads.Add(rname);
                    Interlocked.Increment(ref badReadCount);
                    return;
                }
            }
            //add to count matrix
            positions.GetOrAdd(locus, 0);
            sparseSnpMap[added.leftVariant - 1, added.rightVariant - 1].GetOrAdd(locus,
                t => new ConcurrentDictionary<Tuple<int, int>, int>())
                .AddOrUpdate(new Tuple<int, int>(added.leftSNPidx, added.rightSNPidx),
                    t => 1, (a, b) => b + 1);
            sparseMatMapper[added.leftVariant - 1, added.rightVariant - 1].GetOrAdd(locus,
                t => new LockableDictionary<Tuple<double, double>, int>())
                .AddOrUpdate(new Tuple<double, double>(added.LeftMidpoint(), added.RightMidpoint()),
                    t => 1, (a, b) => b + 1);
        }

        private static double getNormalizer(int left, int right, int chrA, int chrB)
        {
            var normalizer = (contigLengths[chrA-1][left]*contigLengths[chrB-1][right])/1000000000.0;
            return normalizer;
        }


        private static double ComputeGaussianDensity(double x, double y, double x0, double y0)
        {
            var res = Math.Exp(-((x - x0)*(x - x0) + (y - y0)*(y - y0))/(2*Sigma*Sigma));

            if (double.IsNaN(res))
            {
                throw new Exception("problem");
            }
            return res;
        }

        private static ConcurrentDictionary<Tuple<int, int>, double>[,] PhaseByGaussianConvolution(ConcurrentDictionary<Tuple<int, int>, LockableDictionary<Tuple<double, double>, int>>[,] sparseMatMapper)
        {
            var ratioMat = new ConcurrentDictionary<Tuple<int, int>, ConcurrentDictionary<Tuple<double, double>, double[,]>>();

            var phasedMat = new ConcurrentDictionary<Tuple<int, int>, double>[Program.PLOIDITY, Program.PLOIDITY];
            var bins = //observed contig pairs
                new HashSet<Tuple<int, int>>(sparseMatMapper[0, 0].Keys.Union(sparseMatMapper[0, 1].Keys)
                    .Union(sparseMatMapper[1, 1].Keys)
                    .Union(sparseMatMapper[1, 0].Keys));
            Parallel.ForEach(bins, new ParallelOptions(){},
                bin =>
                {
                    var binp1 = sparseMatMapper[0, 0].ContainsKey(bin) ? new HashSet<Tuple<double, double>>(sparseMatMapper[0, 0][bin].Keys) : new HashSet<Tuple<double, double>>();
                    var binp2 = sparseMatMapper[0, 1].ContainsKey(bin) ? new HashSet<Tuple<double, double>>(sparseMatMapper[0, 1][bin].Keys) : new HashSet<Tuple<double, double>>();
                    var binp3 = sparseMatMapper[1, 0].ContainsKey(bin) ? new HashSet<Tuple<double, double>>(sparseMatMapper[1, 0][bin].Keys) : new HashSet<Tuple<double, double>>();
                    var binp4 = sparseMatMapper[1, 1].ContainsKey(bin) ? new HashSet<Tuple<double, double>>(sparseMatMapper[1, 1][bin].Keys) : new HashSet<Tuple<double, double>>();
                    var binposes = new HashSet<Tuple<double, double>>(binp1.Union(binp2).Union(binp3).Union(binp4));
                    Parallel.ForEach(binposes, new ParallelOptions() {  }
                        , binpos =>
                            {
                                var int1 = (sparseMatMapper[0, 0].ContainsKey(bin) &&
                                            sparseMatMapper[0, 0][bin].ContainsKey(binpos))
                                    ? sparseMatMapper[0, 0][bin][binpos]
                                    : 0.0;
                                var int2 = (sparseMatMapper[0, 1].ContainsKey(bin) &&
                                            sparseMatMapper[0, 1][bin].ContainsKey(binpos))
                                    ? sparseMatMapper[0, 1][bin][binpos]
                                    : 0.0;
                                var int3 = (sparseMatMapper[1, 0].ContainsKey(bin) &&
                                            sparseMatMapper[1, 0][bin].ContainsKey(binpos))
                                    ? sparseMatMapper[1, 0][bin][binpos]
                                    : 0.0;
                                var int4 = (sparseMatMapper[1, 1].ContainsKey(bin) &&
                                            sparseMatMapper[1, 1][bin].ContainsKey(binpos))
                                    ? sparseMatMapper[1, 1][bin][binpos]
                                    : 0.0;
                                var bgsum = (int1 + int2 + int3 + int4);
                                ratioMat.GetOrAdd(bin, new ConcurrentDictionary<Tuple<double, double>, double[,]>())
                                    .GetOrAdd(binpos, new[,] {{int1/bgsum, int2/bgsum}, {int3/bgsum, int4/bgsum}});
                            });
                });

            //Go over phased data summing it into the bin
            for (var i = 0; i < Program.PLOIDITY; i++)
                for (var j = 0; j < Program.PLOIDITY; j++)
                {
                        var i1 = i;
                        var j1 = j;
                        phasedMat[i1, j1] = new ConcurrentDictionary<Tuple<int, int>, double>();
                        Parallel.ForEach(sparseMatMapper[i1, j1].Keys, new ParallelOptions() {  },
                            binpos =>
                                phasedMat[i1, j1].AddOrUpdate(binpos,
                                    t => sparseMatMapper[i1, j1][binpos].Sum(u => u.Value),
                                    (a, b) => sparseMatMapper[i1, j1][binpos].Sum(u => u.Value)));
                }

            for (var i = 0; i < Program.PLOIDITY+1; i++)
                for (var j = 0; j < Program.PLOIDITY+1; j++)
                if (i < Program.PLOIDITY && j < Program.PLOIDITY)
                    continue;
                else if (i < Program.PLOIDITY) //Reads where left end is phased
                {
                    int ti = i;
                    int tj = j;
                    Parallel.ForEach(sparseMatMapper[i, j].Keys, new ParallelOptions() {  },
                        binpos => Parallel.ForEach(sparseMatMapper[ti, tj][binpos], new ParallelOptions() {  }, reads =>
                        {
                            var contrib = new double[2];
                            for (var j1 = 0; j1 < Program.PLOIDITY; j1++)
                            {
                                //Ratio of reads in variant at locus * #reads overlapping current position in bin * density of gaussian kernel at variant locus
                                var readContribution = ratioMat.ContainsKey(binpos)
                                    ? ratioMat[binpos].Sum(
                                        locus =>
                                            locus.Value[ti, j1]*
                                            ComputeGaussianDensity(locus.Key.Item1, locus.Key.Item2, reads.Key.Item1,
                                                reads.Key.Item2))
                                                : Thresholds.PhaseNonSnps ? 0.5 : 0.0;
                                contrib[j1] = readContribution;
                            }
                            var contribSum = contrib.Sum();
                            if (contribSum>0)
                            for (var j1 = 0; j1 < Program.PLOIDITY; j1++)
                            {
                                var res = reads.Value * contrib[j1] / contribSum;
                                if (double.IsNaN(res))
                                    throw new Exception("problem");
                                if (res > 0)
                                phasedMat[ti, j1].AddOrUpdate(binpos, t => res,
                                    (a, b) => b + res);
                            }

                        }));
                }
                else if (j < Program.PLOIDITY) //Reads where right end is phased
                {
                    int ti = i;
                    int tj = j;
                    Parallel.ForEach(sparseMatMapper[i, j].Keys, new ParallelOptions() {  },
                        binpos => Parallel.ForEach(sparseMatMapper[ti, tj][binpos], new ParallelOptions() {  }, reads =>
                        {
                            var contrib = new double[2];
                            for (var i1 = 0; i1 < Program.PLOIDITY; i1++)
                            {
                                //Ratio of reads in variant at locus * #reads overlapping current position in bin * density of gaussian kernel at variant locus
                                var readContribution =
                                    ratioMat.ContainsKey(binpos)
                                        ? ratioMat[binpos].Sum(
                                            locus =>
                                                locus.Value[i1, tj]*
                                                ComputeGaussianDensity(locus.Key.Item1, locus.Key.Item2, reads.Key.Item1,
                                                    reads.Key.Item2))
                                        : Thresholds.PhaseNonSnps ? 0.5 : 0.0;
                                contrib[i1] = readContribution;
                            }
                            var contribSum = contrib.Sum();
                            if (contribSum > 0)
                            for (var i1 = 0; i1 < Program.PLOIDITY; i1++)
                            {
                                var res = reads.Value * contrib[i1] / contribSum;
                                if (double.IsNaN(res))
                                    throw new Exception("problem");
                                if (res > 0)
                                phasedMat[i1, tj].AddOrUpdate(binpos, t => res,
                                    (a, b) => reads.Value * (b + res));
                            }

                        }));
                }
                else //Reads where neither end is phased
                {
                    int ti = i;
                    int tj = j;
                    Parallel.ForEach(sparseMatMapper[ti, tj].Keys, new ParallelOptions() {  },
                        binpos => Parallel.ForEach(sparseMatMapper[ti, tj][binpos], new ParallelOptions() {  }, reads =>
                        {
                            var contrib = new double[2, 2];
                            for (var i1 = 0; i1 < Program.PLOIDITY; i1++)
                                for (var j1 = 0; j1 < Program.PLOIDITY; j1++)
                                {
                                    //Ratio of reads in variant at locus * #reads overlapping current position in bin * density of gaussian kernel at variant locus
                                    var readContribution = ratioMat.ContainsKey(binpos)
                                        ? ratioMat[binpos].Sum(
                                            locus =>
                                                locus.Value[i1, j1]*
                                                ComputeGaussianDensity(locus.Key.Item1, locus.Key.Item2, reads.Key.Item1,
                                                    reads.Key.Item2))
                                        : Thresholds.PhaseNonSnps ? 0.25 : 0.0;
                                    contrib[i1, j1] = readContribution;
                                }
                            var contribSum = contrib.Cast<double>().Sum();
                            if (contribSum > 0)
                            for (var i1 = 0; i1 < Program.PLOIDITY; i1++)
                                for (var j1 = 0; j1 < Program.PLOIDITY; j1++)
                                {
                                    var res = reads.Value*contrib[i1, j1]/contribSum;
                                    if(double.IsNaN(res)) 
                                        throw new Exception("problem");
                                    if (res>0)
                                    phasedMat[i1, j1].AddOrUpdate(binpos, t => res,
                                        (a, b) => b + res);
                                }

                        }));
                }
            return phasedMat;
        }


        private static ConcurrentDictionary<Tuple<int, int>, double>[,] PhaseByBinCounts(ConcurrentDictionary<Tuple<int, int>, int>[,] sparseMat, ConcurrentDictionary<Tuple<int, int>, byte> positions)
        {
            var ratios = positions.AsParallel().Select(kvp =>
                {
                    var t = kvp.Key;
                    var int1 = sparseMat[0, 0].ContainsKey(t) ? sparseMat[0, 0][t] : 0;
                    var int2 = sparseMat[0, 1].ContainsKey(t) ? sparseMat[0, 1][t] : 0;
                    var int3 = sparseMat[1, 1].ContainsKey(t) ? sparseMat[1, 1][t] : 0;
                    var int4 = sparseMat[1, 0].ContainsKey(t) ? sparseMat[1, 0][t] : 0;
                    if ((new[] {int1, int2, int3, int4}).Any(v => v != 0))
                        return new
                        {
                            Position = t,
                            IntraAB = int1 / (double)(int1 + int2 + int3 + int4),
                            InterAbeta = int2 / (double)(int1 + int2 + int3 + int4),
                            IntraAlphabeta = int3 / (double)(int1 + int2 + int3 + int4),
                            InterAlphaB = int4 / (double)(int1 + int2 + int3 + int4)
                        };
                    return new
                    {
                        Position = t,
                        IntraAB = 0.25,
                        InterAbeta = 0.25,
                        IntraAlphabeta = 0.25,
                        InterAlphaB = 0.25
                    };
                })
                    .ToDictionary(t => t.Position,
                        t => new[,] {{t.IntraAB, t.InterAbeta}, {t.InterAlphaB, t.IntraAlphabeta}});
            Console.WriteLine(ratios.Count);
            //reassign non-snp reads by ratios
            var phasedMat = new ConcurrentDictionary<Tuple<int, int>, double>[Program.PLOIDITY, Program.PLOIDITY];
            for (var i = 0; i < Program.PLOIDITY + 1; i++)
                for (var j = 0; j < Program.PLOIDITY + 1; j++)
                {
                    var i1 = i;
                    var j1 = j;
                    if (i < Program.PLOIDITY && j < Program.PLOIDITY) //Reads where both ends are phased
                    {
                        phasedMat[i, j] = new ConcurrentDictionary<Tuple<int, int>, double>();
                        Parallel.ForEach(sparseMat[i, j],
                            t => phasedMat[i1, j1].AddOrUpdate(t.Key, u => t.Value, (a, b) => b + t.Value));
                    }
                    else if (i < Program.PLOIDITY) //Reads where left end is phased
                    {
                        Parallel.ForEach(sparseMat[i, j],
                            t =>
                            {
                                var varA = ratios[t.Key][0, 0] + ratios[t.Key][0, 1];
                                var varB = ratios[t.Key][1, 0] + ratios[t.Key][1, 1];
                                var tratios = new[] {varA/(varA + varB), varB/(varA + varB)};
                                for (var a = 0; a < Program.PLOIDITY; a++)
                                    phasedMat[i1, a].AddOrUpdate(t.Key, u => tratios[a]*t.Value,
                                        (x, y) => y + (tratios[a]*t.Value));
                            });
                    }
                    else if (j < Program.PLOIDITY) //Reads where right end is phased
                    {
                        Parallel.ForEach(sparseMat[i, j],
                            t =>
                            {
                                var varA = ratios[t.Key][0, 0] + ratios[t.Key][1, 0];
                                var varB = ratios[t.Key][0, 1] + ratios[t.Key][1, 1];
                                var tratios = new[] {varA/(varA + varB), varB/(varA + varB)};
                                for (var a = 0; a < Program.PLOIDITY; a++)
                                    phasedMat[a, j1].AddOrUpdate(t.Key, u => tratios[a]*t.Value,
                                        (x, y) => y + (tratios[a]*t.Value));
                            });
                    }
                    else //Reads where neither end is phased
                        Parallel.ForEach(sparseMat[i, j],
                            t =>
                            {
                                for (var a = 0; a < Program.PLOIDITY; a++)
                                    for (var b = 0; b < Program.PLOIDITY; b++)
                                        phasedMat[a, b].AddOrUpdate(t.Key, u => ratios[t.Key][a, b]*t.Value,
                                            (x, y) => y + (ratios[t.Key][a, b]*t.Value));
                            });
                }
            return phasedMat;
        }

        private static int GetContigId(int chr, uint pos)
        {
            //Find contigid
            var idx = contigPos[chr].BinarySearch(new KeyValuePair<long, int>(pos, 0), kvpcomp);
            if (idx<0)
            {
                idx = ~idx;
                if (idx==0 || idx==contigPos[chr].Count) return - 1; //before /after first contig
                if (contigPos[chr].ElementAt(idx - 1).Value == contigPos[chr].ElementAt(idx).Value)
                    return contigPos[chr].ElementAt(idx - 1).Value;
                return -1; //between contigs
            }
            return contigPos[chr].ElementAt(idx).Value; //on top of contig
        }

        public class kvpCompare : IComparer<KeyValuePair<long, int>>
        {
            public int Compare(KeyValuePair<long, int> x, KeyValuePair<long, int> y)
            {
                return x.Key.CompareTo(y.Key);
            }

        }
    }
}
