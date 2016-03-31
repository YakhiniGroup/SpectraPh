using System;
using System.Collections.Concurrent;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using System.Threading.Tasks;

namespace SpectraPh.Actions
{
    public class SAM2PSAM : IHaploParserAction
    {
        public string Usage()
        {
            string message = "<InputFolder> <OutputFolder>";
            message += Environment.NewLine + "Iterates on SAM files, generating paired SAMs.";
            return message;
        }

        public string Command()
        {
            return "SAM2PSAM";
        }

        public void Go(string[] args)
        {
            Program.CHROMOSOMES = 23;
            if (!Helpers.GetFilesFromPath(args[0], "*.psam").Any())
            {
                Helpers.Write((int)PrintDetail.Important, "Translating SAM files to PSAM (Paired-SAM) file.");
                TranslateSamToReadPairs(args);
                Helpers.Write((int)PrintDetail.Important, "Done translating PSAM (Paired-SAM) file.");
            }
            else
                Helpers.Write((int)PrintDetail.Important, "Found precomputed PSAM files.");
        }

        private static void TranslateSamToReadPairs(string[] args)
        {
            var files = Helpers.GetFilesFromPath(args[0], "*.sam").ToList();
            var outfilearray = new TextWriter[files.Count, files.Count];
            for (var fidA = 0; fidA < files.Count; fidA++)
                for (var fidB = fidA; fidB < files.Count; fidB++)
                {
                    outfilearray[fidA, fidB] = TextWriter.Synchronized(new StreamWriter(args[1] + string.Format("chr{0}_chr{1}.psam", fidA + 1, fidB + 1), false));
                }
            var linecollection = new BlockingCollection<string>(1000000);
            var producer = Task.Run(() =>
            {
                Parallel.ForEach(files, filename =>
                {
                    using (var fileStream = new StreamReader(filename))
                    {
                        Helpers.Write((int)PrintDetail.Detailed, string.Format("Starting to parse {0}", filename));
                        while (!fileStream.EndOfStream)
                        {
                            linecollection.Add(fileStream.ReadLine());
                        }
                    }
                    Helpers.Write((int)PrintDetail.Detailed, string.Format("Done parsing {0}", filename));
                });
                linecollection.CompleteAdding();
            });
            var tsklst = new List<Task>();
            foreach (var line in linecollection.GetConsumingEnumerable())
            {
                if (tsklst.Count > 100)
                {
                    var tid = Task.WaitAny(tsklst.ToArray());
                    tsklst.RemoveAt(tid);
                }
                tsklst.Add(Task.Run(() =>
                {
                    var splitline = line.Split('\t');
                    string chrA = splitline[2];
                    string chrB = splitline[6];
                    chrB = chrB.Equals("=", StringComparison.InvariantCulture) ? chrA : chrB;
                    int chrAid = Helpers.ChromNum(chrA);
                    int chrBid = Helpers.ChromNum(chrB);
                    var chr1st = Math.Min(chrAid, chrBid);
                    var chr2nd = Math.Max(chrAid, chrBid);
                    var mapQ = int.Parse(splitline[4]);
                    //var PhredLikelihood = int.Parse(splitline[17].Split(':')[2]);
                    if (chrAid > 0 && chrBid > 0 && mapQ > 30) //minimum Mapping quality
                    {
                        outfilearray[chr1st - 1, chr2nd - 1].WriteLine(string.Join("\t",
                            splitline.GetSubsetById(new[] {0, 2, 3, 5, 6, 7, 9, 10})));
                    }
                }));
            }
            for (var fidA = 0; fidA < files.Count; fidA++)
                for (var fidB = fidA; fidB < files.Count; fidB++)
                {
                    outfilearray[fidA, fidB].Close();
                }
            Helpers.Write((int)PrintDetail.Detailed, string.Format("Done parsing all"));
        }

    }
}
